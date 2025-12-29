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
                cpu_mesh->pos[i] = Vec3(x, y, z);
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
    VertexId a,b,c; // sorted ascending
    bool operator==(const FaceKey& o) const { return a==o.a && b==o.b && c==o.c; }
};


struct FaceKeyHash {
    size_t operator()(const FaceKey& k) const {
        // simple mix
        return static_cast<size_t>(k.a)*73856093u
             ^ static_cast<size_t>(k.b)*19349663u
             ^ static_cast<size_t>(k.c)*83492791u;
    }
};

IndexBuffer BuildSurfaceTriangles(const std::vector<tetrahedron>& tets) {
    IndexBuffer out;

    auto make_face = [](VertexId i0, VertexId i1, VertexId i2) {
        if (i0>i1) std::swap(i0,i1);
        if (i1>i2) std::swap(i1,i2);
        if (i0>i1) std::swap(i0,i1);
        return FaceKey{i0,i1,i2};
    };

    std::unordered_map<FaceKey, uint8_t, FaceKeyHash> cnt;
    cnt.reserve(tets.size()*4*2);

    auto add = [&](const VertexId a, const VertexId b, const VertexId c){
        FaceKey k = make_face(a,b,c);
        if (const auto it = cnt.find(k); it == cnt.end()) cnt.emplace(k, 1);
        else ++it->second;
    };

    for (const tetrahedron& t : tets) {
        const auto &v = t.vertices;
        add(v[0], v[1], v[2]);
        add(v[0], v[1], v[3]);
        add(v[0], v[2], v[3]);
        add(v[1], v[2], v[3]);
    }

    out.reserve(cnt.size());

    auto emit_if_boundary = [&](const VertexId a, const VertexId b, const VertexId c){
        if (cnt[make_face(a,b,c)] == 1) {
            out.push_back(a);
            out.push_back(b);
            out.push_back(c);
        }
    };

    for (const tetrahedron& t : tets) {
        const auto &v = t.vertices; // 假定为 [v0,v1,v2,v3] 且整体右手/正向
        emit_if_boundary(v[1], v[2], v[3]);
        emit_if_boundary(v[0], v[3], v[2]);
        emit_if_boundary(v[0], v[1], v[3]);
        emit_if_boundary(v[0], v[2], v[1]);
    }

    return out;
}


/*
 * TODO: his functions is extremely incorrect.
 * Maybe use m_tris index as surface_tri is a good way.
 */

IndexBuffer BuildSurfaceTriangles(const std::vector<triangle>& tris) {
    if (tris.empty()) return {};

    const auto num_tris = static_cast<uint32_t>(tris.size());

    // --- 1. 构建边邻接表 ---
    struct HalfEdge {
        VertexId v0, v1; // 原始方向
        uint32_t tri_idx;
        VertexId edge_key[2]; // 排序后的 ID，用于匹配共享边

        bool operator<(const HalfEdge& o) const {
            if (edge_key[0] != o.edge_key[0]) return edge_key[0] < o.edge_key[0];
            return edge_key[1] < o.edge_key[1];
        }
    };

    std::vector<HalfEdge> all_edges;
    all_edges.reserve(num_tris * 3);
    for (uint32_t i = 0; i < num_tris; ++i) {
        for (int j = 0; j < 3; ++j) {
            VertexId u = tris[i].vertices[j];
            VertexId v = tris[i].vertices[(j + 1) % 3];
            const VertexId k0 = std::min(u, v);
            const VertexId k1 = std::max(u, v);
            all_edges.push_back({u, v, i, {k0, k1}});
        }
    }
    std::sort(all_edges.begin(), all_edges.end());

    struct Neighbor {
        uint32_t tri_idx;
        bool same_direction; // 两个三角形在该边上的采样方向是否一致
    };
    std::vector<std::vector<Neighbor>> adj(num_tris);

    for (size_t i = 0; i < all_edges.size(); ) {
        size_t j = i + 1;
        while (j < all_edges.size() &&
               all_edges[i].edge_key[0] == all_edges[j].edge_key[0] &&
               all_edges[i].edge_key[1] == all_edges[j].edge_key[1]) {

            // 发现共享边：连接两个三角形
            const uint32_t t1 = all_edges[i].tri_idx;
            const uint32_t t2 = all_edges[j].tri_idx;
            // 如果两个三角形在共享边上方向相同，说明其中一个需要翻转
            const bool same_dir = (all_edges[i].v0 == all_edges[j].v0);

            adj[t1].push_back({t2, same_dir});
            adj[t2].push_back({t1, same_dir});
            j++;
        }
        i = j;
    }

    // --- 2. BFS 统一绕序 ---
    std::vector<int8_t> flip(num_tris, 0); // 0: 原样, 1: 翻转
    std::vector visited(num_tris, false);
    std::queue<uint32_t> q;

    for (uint32_t i = 0; i < num_tris; ++i) {
        if (visited[i]) continue;
        visited[i] = true;
        q.push(i);

        while (!q.empty()) {
            const uint32_t u = q.front();
            q.pop();

            for (const auto&[tri_idx, same_direction] : adj[u]) {
                if (!visited[tri_idx]) {
                    visited[tri_idx] = true;
                    // 如果 A 和 B 共享边方向相同，且 A 没翻转，则 B 必须翻转
                    // 逻辑：flip[B] = flip[A] XOR same_direction
                    flip[tri_idx] = flip[u] ^ (same_direction ? 1 : 0);
                    q.push(tri_idx);
                }
            }
        }
    }

    // 构建初步 IndexBuffer
    IndexBuffer out;
    out.reserve(num_tris * 3);
    for (uint32_t i = 0; i < num_tris; ++i) {
        const auto& v = tris[i].vertices;
        if (flip[i] == 0) {
            out.push_back(v[0]); out.push_back(v[1]); out.push_back(v[2]);
        } else {
            out.push_back(v[0]); out.push_back(v[2]); out.push_back(v[1]);
        }
    }

    return out;
}

std::vector<float> ComputeNormal(mesh_on_cpu* cpu_mesh) {
    if (!cpu_mesh) return {};

    const size_t nV = cpu_mesh->size();
    std::vector normals(nV * 3, 0.0f);

    const auto& tris = cpu_mesh->m_surface_tris;
    // Debug 下建议 assert
    // assert(tris.size() % 3 == 0);

    const size_t T = tris.size() / 3;
    for (size_t t = 0; t < T; ++t) {
        const size_t v1 = tris[t * 3 + 0];
        const size_t v2 = tris[t * 3 + 1];
        const size_t v3 = tris[t * 3 + 2];

        if (v1 >= nV || v2 >= nV || v3 >= nV) continue;

        Vec3 e1 = cpu_mesh->pos[v2] - cpu_mesh->pos[v1];
        Vec3 e2 = cpu_mesh->pos[v3] - cpu_mesh->pos[v1];
        Vec3 fn = e1.cross(e2);

        // 可选：跳过退化面
        if (fn.squaredNorm() < 1e-24f) continue;

        auto accum = [&](const size_t vi) {
            normals[3 * vi + 0] += fn.x();
            normals[3 * vi + 1] += fn.y();
            normals[3 * vi + 2] += fn.z();
        };
        accum(v1); accum(v2); accum(v3);
    }

    for (size_t i = 0; i < nV; ++i) {
        Vec3 v(normals[3*i+0], normals[3*i+1], normals[3*i+2]);
        if (const float len = v.norm(); len > 1e-12f) v /= len;
        else v = Vec3(0, 1, 0);

        normals[3*i+0] = v.x();
        normals[3*i+1] = v.y();
        normals[3*i+2] = v.z();
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
        vertices[3*i + 0] = cpu_mesh->pos[i].x();
        vertices[3*i + 1] = cpu_mesh->pos[i].y();
        vertices[3*i + 2] = cpu_mesh->pos[i].z();
    }
    return vertices;
}

void BuildAdjacency(mesh_on_cpu& mesh) {

    const size_t num_nodes = mesh.size();
    auto&[vertex_edges, vertex_faces, vertex_tets] = mesh.adjacencyInfo;
    auto assign_offsets = [num_nodes](std::vector<uint32_t>& offsets) {
        offsets.assign(num_nodes+1, 0u);
    };

    auto build_vertex_incident_csr =
    [num_nodes](auto const& elems,
                const uint32_t verts_per_elem,
                auto&& get_v,
                AdjacencyCSR& adj) {
        auto& offsets = adj.offsets;
        auto& incidents = adj.incidents;

        offsets.assign(num_nodes + 1, 0u);
        incidents.clear();

        if (elems.empty())
            // offsets 全 0，incident 空
            return;

        // 本实现使用 pack(ei,k)，k 占低2位 => k 必须 < 4
        if (verts_per_elem > 4u) {
            throw std::runtime_error("verts_per_elem > 4 is not supported by AdjacencyCSR::pack");
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
        incidents.resize(offsets.back());
        std::vector<uint32_t> cursor = offsets;

        for (size_t ei = 0; ei < elems.size(); ++ei) {
            const auto& e = elems[ei];
            for (uint32_t k = 0; k < verts_per_elem; ++k) {
                const auto v = static_cast<uint32_t>(get_v(e, k));
                const uint32_t dst = cursor[v]++;  // 写入位置
                incidents[dst] = AdjacencyCSR::pack(ei, k);                // 存 incident 元素编号
            }
        }
    };


    // build edge adjacency (vertex -> incident edges)
    if (!mesh.m_edges.empty()) {
        build_vertex_incident_csr(
            mesh.m_edges,
            /*verts_per_elem=*/2u,
            [](auto const& edge, uint32_t k) -> uint32_t {
                return static_cast<uint32_t>(edge.vertices[k]);
            },
            vertex_edges
        );
    }
    else {
        // 没有 edges：保持 CSR 合法
        assign_offsets(vertex_edges.offsets);
        vertex_edges.incidents.clear();
    }

    // build triangle adjacency (vertex -> incident faces)
    if (!mesh.m_tris.empty()) {
        build_vertex_incident_csr(
            mesh.m_tris,
            /*verts_per_elem=*/3u,
            [](auto const& tri, uint32_t k) -> uint32_t {
                return static_cast<uint32_t>(tri.vertices[k]);
            },
            vertex_faces
        );
    }
    else {
        assign_offsets(vertex_faces.offsets);
        vertex_faces.incidents.clear();
    }

    // build tetrahedron adjacency (vertex -> incident tets)
    if (!mesh.m_tets.empty()) {
        build_vertex_incident_csr(
            mesh.m_tets,
            /*verts_per_elem=*/4u,
            [](auto const& tet, uint32_t k) -> uint32_t {
                return static_cast<uint32_t>(tet.vertices[k]);
            },
            vertex_tets
        );
    }
    else {
        assign_offsets(vertex_tets.offsets);
        vertex_tets.incidents.clear();
    }
}









