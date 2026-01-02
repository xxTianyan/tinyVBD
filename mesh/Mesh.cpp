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

    auto read_vertex_id_0based = [&](VertexID& out) {
        unsigned long v_l = 0;
        in >> v_l;

        // 基本合法性：Gmsh 节点编号从 1 开始
        if (v_l == 0) throw std::runtime_error("Invalid vertex id (0) in msh element");
        if (num_nodes > 0 && v_l > num_nodes) {
            throw std::runtime_error("Invalid vertex id (out of range) in msh element");
        }

        out = static_cast<VertexID>(v_l - 1); // 1-based -> 0-based
    };

    while (std::getline(in, line)) {
        line = trim(line);
        // -------- Nodes (MSH v2 ASCII) --------
        if (line == "$Nodes") {
            in >> num_nodes;
            if (num_nodes >= INVALID_VERTEX_ID) {
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
                // Gmsh v2: 2=tri(3), 4=tet(4)
                if (elm_type == 2) {
                    VertexID v0, v1, v2;
                    read_vertex_id_0based(v0);
                    read_vertex_id_0based(v1);
                    read_vertex_id_0based(v2);
                    const auto& v0_pos = cpu_mesh->pos[v0];
                    const auto& v1_pos = cpu_mesh->pos[v1];
                    const auto& v2_pos = cpu_mesh->pos[v2];
                    cpu_mesh->m_tris.emplace_back(v0, v1, v2, v0_pos, v1_pos, v2_pos);
                    have_tris = true;
                }
                else if (elm_type == 4) {
                    VertexID v0, v1, v2, v3;
                    read_vertex_id_0based(v0);
                    read_vertex_id_0based(v1);
                    read_vertex_id_0based(v2);
                    read_vertex_id_0based(v3);
                    const auto& v0_pos = cpu_mesh->pos[v0];
                    const auto& v1_pos = cpu_mesh->pos[v1];
                    const auto& v2_pos = cpu_mesh->pos[v2];
                    const auto& v3_pos = cpu_mesh->pos[v3];
                    cpu_mesh->m_tets.emplace_back(v0, v1, v2, v3, v0_pos, v1_pos, v2_pos, v3_pos);
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
    VertexID a,b,c; // sorted ascending
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

    auto make_face = [](VertexID i0, VertexID i1, VertexID i2) {
        if (i0>i1) std::swap(i0,i1);
        if (i1>i2) std::swap(i1,i2);
        if (i0>i1) std::swap(i0,i1);
        return FaceKey{i0,i1,i2};
    };

    std::unordered_map<FaceKey, uint8_t, FaceKeyHash> cnt;
    cnt.reserve(tets.size()*4*2);

    auto add = [&](const VertexID a, const VertexID b, const VertexID c){
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

    auto emit_if_boundary = [&](const VertexID a, const VertexID b, const VertexID c){
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


IndexBuffer BuildSurfaceTriangles(const std::vector<triangle>& tris) {
    IndexBuffer out;
    if (tris.empty()) return out;

    const auto num_tris = static_cast<uint32_t>(tris.size());

    struct HalfEdge {
        VertexID v0, v1;       // 有向边 v0->v1（来自三角形绕序）
        uint32_t tri_idx;

        VertexID k0, k1;       // 无向边 key：min/max，用来匹配共享边

        bool operator<(const HalfEdge& o) const
        {
            if (k0 != o.k0) return k0 < o.k0;
            return k1 < o.k1;
        }
    };

    // 1) 收集所有 half-edges
    std::vector<HalfEdge> edges;
    edges.reserve(static_cast<size_t>(num_tris) * 3);

    for (uint32_t t = 0; t < num_tris; ++t) {
        const auto& v = tris[t].vertices;
        for (int e = 0; e < 3; ++e) {
            VertexID a = v[e];
            VertexID b = v[(e + 1) % 3];
            VertexID k0 = std::min(a, b);
            VertexID k1 = std::max(a, b);
            edges.push_back(HalfEdge{a, b, t, k0, k1});
        }
    }

    std::sort(edges.begin(), edges.end());

    struct Neighbor {
        uint32_t tri_idx;
        bool same_dir; // 两个三角形在共享边上的有向边方向是否相同
    };

    std::vector<std::vector<Neighbor>> adj(num_tris);

    // 2) 按无向边分组，建立邻接与方向约束
    bool has_nonmanifold = false;
    for (size_t i = 0; i < edges.size(); ) {
        size_t j = i + 1;
        while (j < edges.size() && edges[j].k0 == edges[i].k0 && edges[j].k1 == edges[i].k1) ++j;

        const size_t group_size = j - i;

        // group_size == 1: 边界边（只有一个三角形拥有该边），无邻接，OK
        // group_size == 2: 流形共享边，理想情况
        // group_size  > 2: 非流形边，可能不可定向/数据异常，需要标记并尽量加约束
        if (group_size >= 2) {
            if (group_size > 2) has_nonmanifold = true;

            // 组内两两建立约束（k 通常很小，O(k^2)可接受）
            for (size_t a = i; a < j; ++a) {
                for (size_t b = a + 1; b < j; ++b)
                {
                    const auto& e0 = edges[a];
                    const auto& e1 = edges[b];
                    const uint32_t t0 = e0.tri_idx;
                    const uint32_t t1 = e1.tri_idx;

                    // 同一条无向边下，e0 与 e1 要么同向(a->b & a->b)，要么反向(a->b & b->a)
                    const bool same_dir = (e0.v0 == e1.v0 && e0.v1 == e1.v1);

                    adj[t0].push_back(Neighbor{t1, same_dir});
                    adj[t1].push_back(Neighbor{t0, same_dir});
                }
            }
        }

        i = j;
    }

    // 3) BFS/二染色：flip[t] = 0/1
    // 使用 -1 表示未赋值；遇到矛盾时记录但不崩溃
    std::vector<int8_t> flip(num_tris, int8_t(-1));
    std::queue<uint32_t> q;

    bool has_conflict = false;

    for (uint32_t seed = 0; seed < num_tris; ++seed) {
        if (flip[seed] != -1) continue;

        flip[seed] = 0; // 该连通分量的“基准方向”任意选一个
        q.push(seed);

        while (!q.empty()) {
            const uint32_t u = q.front();
            q.pop();

            for (const auto& nb : adj[u]) {
                const uint32_t v = nb.tri_idx;

                // 约束：若共享边同向，则两三角形必须一翻一不翻；若反向，则翻转状态相同
                // 即 flip[v] = flip[u] XOR same_dir
                const auto required = static_cast<int8_t>(flip[u] ^ (nb.same_dir ? 1 : 0));

                if (flip[v] == -1) {
                    flip[v] = required;
                    q.push(v);
                }
                else if (flip[v] != required) {
                    // 数据存在矛盾：不可定向或网格/拓扑有问题
                    // 这里不强行修改已定值，避免振荡；仅记录冲突
                    has_conflict = true;
                }
            }
        }
    }

    // 4) 输出 index buffer（根据 flip 统一绕序）
    out.reserve(static_cast<size_t>(num_tris) * 3);

    for (uint32_t t = 0; t < num_tris; ++t) {
        const auto& v = tris[t].vertices;
        if (flip[t] == 1) {
            // 翻转：交换 (1,2)
            out.push_back(v[0]);
            out.push_back(v[2]);
            out.push_back(v[1]);
        }
        else {
            out.push_back(v[0]);
            out.push_back(v[1]);
            out.push_back(v[2]);
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
std::vector<float> AssembleVertices(const mesh_on_cpu* cpu_mesh) {
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
            return;

        // we use low 2 bit to mask vertex order, thus the number of vertex must less than 4
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

void DistributeMass(mesh_on_cpu& mesh) {
    auto& inv_mass = mesh.inv_mass;
    std::fill(inv_mass.begin(), inv_mass.end(), 1.0f);
}

void InitMesh(mesh_on_cpu& mesh) {
    BuildAdjacency(mesh);
    DistributeMass(mesh);
}





