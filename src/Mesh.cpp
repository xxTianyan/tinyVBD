//
// Created by 徐天焱 on 2025/11/6.
//

#include "Mesh.h"
#include <fstream>

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
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Could not open file " + path);

    std::string line;
    bool have_nodes = false, have_tets = false;

    while (std::getline(in, line)) {
        line = trim(line);
        if (line == "$Nodes") {
            size_t num_nodes = 0;
            in >> num_nodes;
            if (num_nodes > INVALID_VERTEX_ID) throw std::runtime_error("Too Much nodes for a single mesh.");
            cpu_mesh->resize(num_nodes);
            for (size_t i = 0; i < num_nodes; i++) {
                size_t dummy;
                in >> dummy >> cpu_mesh->px[i] >> cpu_mesh->py[i] >> cpu_mesh->pz[i];
            }
            have_nodes = true;
        }
        else if (line == "$Elements") {
            size_t num_tets = 0;
            in >> num_tets;
            cpu_mesh->m_tets_local.reserve(num_tets);
            for (size_t i = 0; i < num_tets; i++) {
                unsigned long idx_l, v1_l, v2_l, v3_l, v4_l, dummy_l;

                in >> idx_l >> dummy_l >> dummy_l >> v1_l >> v2_l >> v3_l >> v4_l;

                if (!isValidVertexId(v1_l) ||
                    !isValidVertexId(v2_l) ||
                    !isValidVertexId(v3_l) ||
                    !isValidVertexId(v4_l)) {
                    throw std::runtime_error("Invalid vertex id");
                    }

                const auto v1 = static_cast<VertexId>(v1_l);
                const auto v2 = static_cast<VertexId>(v2_l);
                const auto v3 = static_cast<VertexId>(v3_l);
                const auto v4 = static_cast<VertexId>(v4_l);

                cpu_mesh->m_tets_local.emplace_back(v1 - 1, v2 - 1, v3 - 1, v4 - 1);
            }
            have_tets = true;
        }

        // only vertex and tetrahedron in msh file
        if (have_nodes && have_tets) break;
    }

    if (!have_nodes) throw std::runtime_error("MSH parse error: no $Nodes section");
    if (!have_tets) throw std::runtime_error("MSH parse error: no $Elements section");

    // what is necessary
    cpu_mesh->m_surface_tris_local = BuildSurfaceTriangles(cpu_mesh->m_tets_local);
}

NodeTetAdj BuildNodeTetAdj(const size_t num_nodes, const std::vector<tetrahedron>& tets) {

    NodeTetAdj adj;
    adj.offsets.assign(num_nodes+1, 0);

    for (const auto & tet : tets) {
        for (size_t k = 0; k < 4; k++) {
            const VertexId v = tet.vertices[k];
            adj.offsets[v+1]++;
        }
    }

    for (size_t i = 1; i < num_nodes+1; i++) {
        adj.offsets[i] += adj.offsets[i-1];
    }

    adj.incidentTets.resize(adj.offsets.back());

    auto cursor = adj.offsets;

    for (size_t i = 0; i < tets.size(); i++) {
        for (size_t k = 0; k < 4; k++) {
            const VertexId v = tets[i].vertices[k];
            const uint32_t c = cursor[v]++;
            adj.incidentTets[c] = static_cast<uint32_t>(i);
        }
    }
    return adj;
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
    IndexBuffer tri_indices;

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

    tri_indices.reserve(cnt.size());

    auto emit_if_boundary = [&](const VertexId a, const VertexId b, const VertexId c){
        if (cnt[make_face(a,b,c)] == 1) {
            tri_indices.push_back(a);
            tri_indices.push_back(b);
            tri_indices.push_back(c);
        }
    };

    for (const tetrahedron& t : tets) {
        const auto &v = t.vertices; // 假定为 [v0,v1,v2,v3] 且整体右手/正向
        emit_if_boundary(v[1], v[2], v[3]);
        emit_if_boundary(v[0], v[3], v[2]);
        emit_if_boundary(v[0], v[1], v[3]);
        emit_if_boundary(v[0], v[2], v[1]);
    }

    return tri_indices;
}

std::vector<float> ComputeNormal(mesh_on_cpu* cpu_mesh) {
    std::vector<float> normals;
    normals.resize((cpu_mesh->size() * 3));
    const size_t T = cpu_mesh->m_surface_tris_local.size() / 3;

    for (size_t t = 0; t < T; t++) {
        const size_t v1 = cpu_mesh->m_surface_tris_local[t*3+0];
        const size_t v2 = cpu_mesh->m_surface_tris_local[t*3+1];
        const size_t v3 = cpu_mesh->m_surface_tris_local[t*3+2];

        Vec3 a(cpu_mesh->px[v1], cpu_mesh->py[v1], cpu_mesh->pz[v1]);
        Vec3 b(cpu_mesh->px[v2], cpu_mesh->py[v2], cpu_mesh->pz[v2]);
        Vec3 c(cpu_mesh->px[v3], cpu_mesh->py[v3], cpu_mesh->pz[v3]);

        Vec3 n = (b - a).cross(c - a);

        if (constexpr float EPS2 = 1e-30f; n.squaredNorm() < EPS2) continue;

        auto accum = [&](const size_t tri_idx) {
            normals[3*tri_idx + 0] += n.x();
            normals[3*tri_idx + 1] += n.y();
            normals[3*tri_idx + 2] += n.z();
        };
        accum(v1);accum(v2);accum(v3);
    }

    for (size_t i = 0; i < cpu_mesh->size(); i++) {
        Vec3 v(normals[3*i+0], normals[3*i+1], normals[3*i+2]);
        float norm = v.norm();
        v = v / norm;
        normals[3*i + 0] = cpu_mesh->nx[i] = v.x();
        normals[3*i + 1] = cpu_mesh->ny[i] = v.y();
        normals[3*i + 2] = cpu_mesh->nz[i] = v.z();
    }

    return normals;
}

// assemble Vertex position (XYZ - 3 components per vertex) (shader-location = 0)
std::vector<float> assemble_vertices(const mesh_on_cpu* cpu_mesh) {
    std::vector<float> vertices;
    const size_t num_nodes = cpu_mesh->size();
    vertices.resize(num_nodes * 3);

    for (size_t i = 0; i < num_nodes; i++) {
        vertices[3*i + 0] = cpu_mesh->px[i];
        vertices[3*i + 1] = cpu_mesh->py[i];
        vertices[3*i + 2] = cpu_mesh->pz[i];
    }
    return vertices;
}

// Mesh.cpp

void mesh_on_cpu::InitializePhysics(const float density) {
    tet_Dm_inv.resize(m_tets_local.size());
    tet_vol.resize(m_tets_local.size());
    std::fill(mass.begin(), mass.end(), 0.0f);

    for (size_t i = 0; i < m_tets_local.size(); ++i) {
        const auto& tet = m_tets_local[i];
        const VertexId id0 = tet.vertices[0];
        const VertexId id1 = tet.vertices[1];
        const VertexId id2 = tet.vertices[2];
        const VertexId id3 = tet.vertices[3];

        Vec3 p0(px[id0], py[id0], pz[id0]);
        Vec3 p1(px[id1], py[id1], pz[id1]);
        Vec3 p2(px[id2], py[id2], pz[id2]);
        Vec3 p3(px[id3], py[id3], pz[id3]);

        // Compute rest shape matrix Dm (using first three edge vectors)
        // deformation gradient F = Ds Dm^{-1}
        Mat3 Dm;
        Dm.col(0) = p1 - p0;
        Dm.col(1) = p2 - p0;
        Dm.col(2) = p3 - p0;

        // Compute inverse of Dm
        if (std::abs(Dm.determinant()) < 1e-6) {
            tet_Dm_inv[i] = Mat3::Identity();
            tet_vol[i] = 0.0f;
        }
        else {
            tet_Dm_inv[i] = Dm.inverse();
            tet_vol[i] = std::abs(Dm.determinant()) / 6.0f;
        }

        // Distribute mass to vertices (simple volume/4 approach)
        const float pm = tet_vol[i] * density / 4.0f;
        mass[id0] += pm; mass[id1] += pm;
        mass[id2] += pm; mass[id3] += pm;
    }

    // Compute inverse mass, set fixed points to infinity (inv_mass=0)
    for(size_t i=0; i<size(); ++i) {
        if(mass[i] > 1e-9) inv_mass[i] = 1.0f / mass[i];
        else inv_mass[i] = 0.0f;
    }

    inited = true;
}