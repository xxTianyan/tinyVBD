//
// Created by 徐天焱 on 2025/11/6.
//

#include "Mesh.h"
#include <fstream>
#include <queue>

static std::string trim(const std::string& str) {
    const size_t a = str.find_first_not_of("\t\r\n");
    if (a == std::string::npos)
        return "";
    const size_t b = str.find_last_not_of("\t\r\n");
    return str.substr(a, b - a + 1);
}

static bool start_with(const std::string& str, const char* pfx) {
    const size_t n = std::strlen(pfx);
    return str.size() >= n && std::equal(pfx, pfx + n, str.begin());
}

void ParseMSH(const std::string& path, mesh_on_cpu* cpu_mesh) {

    if (!cpu_mesh) throw std::runtime_error("ParseMSH: cpu_mesh is null");

    std::ifstream in(path);
    if (!in) throw std::runtime_error("Could not open file " + path);

    // clear old date
    cpu_mesh->clear_topology();
    size_t num_nodes = 0;

    std::string line;
    bool have_nodes = false;
    bool have_tris  = false;
    bool have_tets  = false;

    auto read_vertex_id_0based = [&](VertexId& out) {
        unsigned long v_l = 0;
        in >> v_l;

        // 基本合法性：Gmsh 节点编号从 1 开始
        if (v_l == 0) throw std::runtime_error("Invalid vertex id (0) in msh element");
        if (num_nodes > 0 && v_l > num_nodes) {
            throw std::runtime_error("Invalid vertex id (out of range) in msh element");
        }

        out = static_cast<VertexId>(v_l - 1); // 1-based -> 0-based
    };

    while (std::getline(in, line)) {
        line = trim(line);
        // -------- Nodes (MSH v2 ASCII) --------
        if (line == "$Nodes") {
            in >> num_nodes;
            if (num_nodes > INVALID_VERTEX_ID) {
                throw std::runtime_error("Too many nodes for a single mesh.");
            }

            cpu_mesh->resize(num_nodes);

            for (size_t i = 0; i < num_nodes; i++) {
                unsigned long node_id = 0;
                float x, y, z;
                in >> node_id >> x >> y >> z;
                cpu_mesh->p[i] = Vec3(x, y, z);
                // here we don't ask node id to be continued, but if you want it strictly continue：
                // if (node_id != i+1) ...
            }
            have_nodes = true;
            continue;
        }
        // -------- Elements (MSH v2 ASCII) --------
        if (line == "$Elements") {
            size_t num_elements = 0;
            in >> num_elements;

            // 粗略 reserve（真实数量未知，先按总数 reserve，避免频繁扩容）
            cpu_mesh->m_edges.reserve(num_elements);
            cpu_mesh->m_tris.reserve(num_elements);
            cpu_mesh->m_tets.reserve(num_elements);

            for (size_t e = 0; e < num_elements; ++e) {
                unsigned long elm_id = 0;
                unsigned long elm_type   = 0;
                unsigned long num_tags   = 0;

                in >> elm_id >> elm_type >> num_tags;

                // 跳过 tags（物理组/几何实体等）
                for (unsigned long t = 0; t < num_tags; ++t) {
                    unsigned long tag_dummy = 0;
                    in >> tag_dummy;
                }

                // 根据 elementType 读取节点
                // Gmsh v2: 1=edge(2), 2=tri(3), 4=tet(4)
                if (elm_type == 1) {
                    VertexId v0, v1;
                    read_vertex_id_0based(v0);
                    read_vertex_id_0based(v1);
                    cpu_mesh->m_edges.emplace_back(v0, v1);
                }
                else if (elm_type == 2) {
                    VertexId v0, v1, v2;
                    read_vertex_id_0based(v0);
                    read_vertex_id_0based(v1);
                    read_vertex_id_0based(v2);
                    cpu_mesh->m_tris.emplace_back(v0, v1, v2);
                    have_tris = true;
                }
                else if (elm_type == 4) {
                    VertexId v0, v1, v2, v3;
                    read_vertex_id_0based(v0);
                    read_vertex_id_0based(v1);
                    read_vertex_id_0based(v2);
                    read_vertex_id_0based(v3);
                    cpu_mesh->m_tets.emplace_back(v0, v1, v2, v3);
                    have_tets = true;
                }
                else {
                    // 其他 elementType：不支持，丢弃该行剩余内容
                    //（在 MSH v2 ASCII 中，每个元素占一行，getline 丢掉剩余 token 即可）
                    std::getline(in, line);
                }
            }
            continue;
        }
    }

    if (!have_nodes) {
        throw std::runtime_error("MSH parse error: no $Nodes section");
    }

    // shrink to show accurate number of elements
    cpu_mesh->m_edges.shrink_to_fit();
    cpu_mesh->m_tris.shrink_to_fit();
    cpu_mesh->m_tets.shrink_to_fit();

    // what is necessary
    if (have_tris)
        cpu_mesh->m_surface_tris = BuildSurfaceTriangles(cpu_mesh->m_tris);
    else if (have_tets)
        cpu_mesh->m_surface_tris = BuildSurfaceTriangles(cpu_mesh->m_tets);
    else
        throw std::runtime_error("NO tris or tets for building surface triangles");

}

struct FaceKey {
    VertexId a, b, c; // sorted ascending
    bool operator==(const FaceKey& o) const { return a==o.a && b==o.b && c==o.c; }
};

struct FaceKeyHash {
    size_t operator()(const FaceKey& k) const {
        return static_cast<size_t>(k.a)*73856093u
             ^ static_cast<size_t>(k.b)*19349663u
             ^ static_cast<size_t>(k.c)*83492791u;
    }
};

struct EdgeKey {
    VertexId a, b; // sorted ascending
    bool operator==(const EdgeKey& o) const { return a==o.a && b==o.b; }
};

struct EdgeKeyHash {
    size_t operator()(const EdgeKey& k) const {
        return static_cast<size_t>(k.a)*73856093u ^ static_cast<size_t>(k.b)*19349663u;
    }
};

static inline FaceKey MakeFaceKey(VertexId i0, VertexId i1, VertexId i2) {
    if (i0 > i1) std::swap(i0, i1);
    if (i1 > i2) std::swap(i1, i2);
    if (i0 > i1) std::swap(i0, i1);
    return FaceKey{i0, i1, i2};
}

static inline EdgeKey MakeEdgeKey(VertexId i0, VertexId i1) {
    if (i0 > i1) std::swap(i0, i1);
    return EdgeKey{i0, i1};
}

IndexBuffer BuildSurfaceTriangles(const std::vector<triangle>& tris) {
    IndexBuffer out;
    out.reserve(tris.size() * 3);

    if (tris.empty()) return out;

    struct HalfEdgeInfo {
        uint32_t tri;     // triangle index
        VertexId from;    // directed edge: from -> to (as appears in triangle order)
        VertexId to;
    };

    // 对每条无向边，记录出现过的半边（通常 1 或 2 个；非流形可能 >2）
    std::unordered_map<EdgeKey, std::vector<HalfEdgeInfo>, EdgeKeyHash> edge_map;
    edge_map.reserve(tris.size() * 3);

    auto add_half_edge = [&](const uint32_t ti, const VertexId u, const VertexId v) {
        const EdgeKey ek = MakeEdgeKey(u, v);
        edge_map[ek].push_back(HalfEdgeInfo{ti, u, v}); // 方向按原三角形顶点顺序记录
    };

    for (uint32_t ti = 0; ti < static_cast<uint32_t>(tris.size()); ++ti) {
        const auto& v = tris[ti].vertices; // v0,v1,v2
        add_half_edge(ti, v[0], v[1]);
        add_half_edge(ti, v[1], v[2]);
        add_half_edge(ti, v[2], v[0]);
    }

    // 构建三角形邻接图：边共享则连边，并记录 same_dir（两三角对该边的方向是否相同）
    struct NeighborRel {
        uint32_t nb;
        bool same_dir; // true: 两个三角在该共享边上方向相同；一致绕序要求共享边方向相反
    };

    std::vector<std::vector<NeighborRel>> adj(tris.size());
    for (auto&[fst, snd] : edge_map) {
        auto& inc = snd;
        if (inc.size() < 2) continue;

        // 将该无向边的所有 incident 三角连接到第一个（处理非流形边时至少不崩）
        const auto& base = inc[0];
        for (size_t j = 1; j < inc.size(); ++j) {
            const auto&[tri, from, to] = inc[j];

            const bool same_dir = (base.from == from && base.to == to);
            adj[base.tri].push_back(NeighborRel{tri, same_dir});
            adj[tri].push_back(NeighborRel{base.tri, same_dir});
        }
    }

    // BFS/DFS 传播 flip：要求共享边方向相反。
    // 推导关系：flipB = flipA XOR same_dir
    std::vector<int8_t> visited(tris.size(), 0);
    std::vector<int8_t> flip(tris.size(), 0);

    std::queue<uint32_t> q;
    for (uint32_t s = 0; s < static_cast<uint32_t>(tris.size()); ++s) {
        if (visited[s]) continue;
        visited[s] = 1;
        flip[s] = 0; // 以该连通块第一片三角的输入绕序为“基准正面”
        q.push(s);

        while (!q.empty()) {
            const uint32_t u = q.front(); q.pop();
            for (const auto&[nb, same_dir] : adj[u]) {
                const uint32_t v = nb;
                const auto desired = static_cast<int8_t>(flip[u] ^ (same_dir ? 1 : 0));
                if (!visited[v]) {
                    visited[v] = 1;
                    flip[v] = desired;
                    q.push(v);
                } else {
                    // 非流形/自相交/输入绕序混乱时，可能出现矛盾；此处不抛异常，只保持已有结果
                    // 如需严格，可在这里检测 flip[v] != desired 并报错
                }
            }
        }
    }

    // 输出 index buffer：flip==1 时交换 v1,v2 反转绕序
    for (uint32_t ti = 0; ti < static_cast<uint32_t>(tris.size()); ++ti) {
        const auto& v = tris[ti].vertices;
        out.push_back(v[0]);
        if (flip[ti] == 0) {
            out.push_back(v[1]);
            out.push_back(v[2]);
        } else {
            out.push_back(v[2]);
            out.push_back(v[1]);
        }
    }

    return out;
}


IndexBuffer BuildSurfaceTriangles(const std::vector<tetrahedron>& tets) {
    IndexBuffer out;
    if (tets.empty()) return out;

    struct FaceRec {
        uint32_t count = 0;
        VertexId oa = 0, ob = 0, oc = 0; // oriented (as we want to emit if boundary)
    };

    std::unordered_map<FaceKey, FaceRec, FaceKeyHash> mp;
    mp.reserve(tets.size() * 4 * 2);

    auto add_face = [&](const VertexId a, const VertexId b, const VertexId c,
                        const VertexId oa, const VertexId ob, const VertexId oc) {
        FaceKey k = MakeFaceKey(a, b, c);
        if (const auto it = mp.find(k); it == mp.end()) {
            FaceRec rec;
            rec.count = 1;
            rec.oa = oa; rec.ob = ob; rec.oc = oc;
            mp.emplace(k, rec);
        } else {
            it->second.count += 1;
        }
    };

    // 约定：tet 顶点顺序 v0,v1,v2,v3 为“正向/右手/正体积”
    // 则外向三角面（与常见 FEM/几何约定一致）可取：
    // face opposite v0: (v1, v2, v3)
    // face opposite v1: (v0, v3, v2)
    // face opposite v2: (v0, v1, v3)
    // face opposite v3: (v0, v2, v1)
    for (const tetrahedron& t : tets) {
        const auto& v = t.vertices;

        // 用 (a,b,c) 做 key（无向），用 (oa,ob,oc) 记录有向输出
        add_face(v[1], v[2], v[3],  v[1], v[2], v[3]);
        add_face(v[0], v[2], v[3],  v[0], v[3], v[2]);
        add_face(v[0], v[1], v[3],  v[0], v[1], v[3]);
        add_face(v[0], v[1], v[2],  v[0], v[2], v[1]);
    }

    out.reserve(mp.size() * 3);

    for (const auto&[fst, snd] : mp) {
        if (const auto&[count, oa, ob, oc] = snd; count == 1) {
            out.push_back(oa);
            out.push_back(ob);
            out.push_back(oc);
        }
    }

    return out;
}

std::vector<float> ComputeNormal(mesh_on_cpu* cpu_mesh) {
    std::vector<float> normals;
    normals.resize((cpu_mesh->size() * 3));
    const size_t T = cpu_mesh->m_surface_tris.size() / 3;

    for (size_t t = 0; t < T; t++) {
        const size_t v1 = cpu_mesh->m_surface_tris[t*3+0];
        const size_t v2 = cpu_mesh->m_surface_tris[t*3+1];
        const size_t v3 = cpu_mesh->m_surface_tris[t*3+2];

        Vec3 edge1 = cpu_mesh->p[v2] - cpu_mesh->p[v1];
        Vec3 edge2 = cpu_mesh->p[v3] - cpu_mesh->p[v1];
        const Vec3 face_n = edge1.cross(edge2);

        auto accum = [&](const size_t tri_idx) {
            normals[3*tri_idx + 0] += face_n.x();
            normals[3*tri_idx + 1] += face_n.y();
            normals[3*tri_idx + 2] += face_n.z();
        };
        accum(v1);accum(v2);accum(v3);
    }

    for (size_t i = 0; i < cpu_mesh->size(); i++) {
        Vec3 v(normals[3*i+0], normals[3*i+1], normals[3*i+2]);
        float norm = v.norm();
        v = v / norm;
        if (norm > 1e-12f) v /= norm;
        else v = Vec3(0, 1, 0); // 默认向上或其他安全值
        normals[3*i + 0] = v.x();
        normals[3*i + 1] = v.y();
        normals[3*i + 2] = v.z();
        cpu_mesh->n[i] = v;
    }

    return normals;
}

// assemble Vertex position (XYZ - 3 components per vertex) (shader-location = 0)
std::vector<float> assemble_vertices(const mesh_on_cpu* cpu_mesh) {
    std::vector<float> vertices;
    const size_t num_nodes = cpu_mesh->size();
    vertices.resize(num_nodes * 3);

    for (size_t i = 0; i < num_nodes; i++) {
        vertices[3*i + 0] = cpu_mesh->p[i].x();
        vertices[3*i + 1] = cpu_mesh->p[i].y();
        vertices[3*i + 2] = cpu_mesh->p[i].z();
    }
    return vertices;
}

void BuildAdjacency(mesh_on_cpu* mesh) {

    if (!mesh)
        return;

    const size_t num_nodes = mesh->size();
    auto& adj = mesh->adjacencyInfo;
    auto assign_offsets = [num_nodes](std::vector<uint32_t>& offsets) {
        offsets.assign(num_nodes+1, 0u);
    };

    auto build_vertex_incident_csr =
    [num_nodes](auto const& elems,
                const uint32_t verts_per_elem,
                auto&& get_v,
                std::vector<uint32_t>& offsets,
                std::vector<uint32_t>& incident) {
        offsets.assign(num_nodes + 1, 0u);
        incident.clear();

        if (elems.empty()) {
            // offsets 全 0，incident 空
            return;
        }

        // 1) count degree
        for (uint32_t ei = 0; ei < static_cast<uint32_t>(elems.size()); ++ei) {
            const auto& e = elems[ei];
            for (uint32_t k = 0; k < verts_per_elem; ++k) {
                const auto v = static_cast<uint32_t>(get_v(e, k));
                assert(v < num_nodes);
                offsets[v + 1] += 1u;
            }
        }

        // 2) prefix sum
        for (size_t i = 1; i < offsets.size(); ++i) {
            offsets[i] += offsets[i - 1];
        }

        // 3) fill
        incident.resize(offsets.back());
        std::vector<uint32_t> cursor = offsets;

        for (uint32_t ei = 0; ei < static_cast<uint32_t>(elems.size()); ++ei) {
            const auto& e = elems[ei];
            for (uint32_t k = 0; k < verts_per_elem; ++k) {
                const auto v = static_cast<uint32_t>(get_v(e, k));
                const uint32_t dst = cursor[v]++;  // 写入位置
                incident[dst] = ei;                // 存 incident 元素编号
            }
        }
    };


    // build edge adjacency (vertex -> incident edges)
    if (!mesh->m_edges.empty()) {
        build_vertex_incident_csr(
            mesh->m_edges,
            /*verts_per_elem=*/2u,
            [](auto const& edge, uint32_t k) -> uint32_t {
                return static_cast<uint32_t>(edge.vertices[k]);
            },
            adj.v_adj_edges_offsets,
            adj.v_adj_edges
        );
    }
    else {
        // 没有 edges：保持 CSR 合法
        assign_offsets(adj.v_adj_edges_offsets);
        adj.v_adj_edges.clear();
    }

    // build triangle adjacency (vertex -> incident faces)
    if (!mesh->m_tris.empty()) {
        build_vertex_incident_csr(
            mesh->m_tris,
            /*verts_per_elem=*/3u,
            [](auto const& tri, uint32_t k) -> uint32_t {
                return static_cast<uint32_t>(tri.vertices[k]);
            },
            adj.v_adj_faces_offsets,
            adj.v_adj_faces
        );
    }
    else {
        assign_offsets(adj.v_adj_faces_offsets);
        adj.v_adj_faces.clear();
    }

    // build tetrahedron adjacency (vertex -> incident tets)
    if (!mesh->m_tets.empty()) {
        build_vertex_incident_csr(
            mesh->m_tets,
            /*verts_per_elem=*/4u,
            [](auto const& tet, uint32_t k) -> uint32_t {
                return static_cast<uint32_t>(tet.vertices[k]);
            },
            adj.v_adj_tet_offsets,
            adj.v_adj_tets
        );
    }
    else {
        assign_offsets(adj.v_adj_tet_offsets);
        adj.v_adj_tets.clear();
    }
}









