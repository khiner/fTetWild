// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#ifndef FLOATTETWILD_LOCALOPERATIONS_H
#define FLOATTETWILD_LOCALOPERATIONS_H

#include <floattetwild/Mesh.hpp>
#include <floattetwild/AABBWrapper.h>
#include <floattetwild/Predicates.hpp>

#include <unordered_set>

namespace floatTetWild {
    extern bool use_old_energy;

    // The position of corner j of tet t_id.
    inline const Vector3& tet_pos(const Mesh& mesh, int t_id, int j) {
        return mesh.tet_vertices[mesh.tets[t_id][j]].pos;
    }

    int get_opp_t_id(const Mesh& mesh, int t_id, int j);
    inline int get_local_f_id(int t_id, int v1_id, int v2_id, int v3_id, const Mesh &mesh) {
        for (int j = 0; j < 4; j++) {
            if (mesh.tets[t_id][j] != v1_id && mesh.tets[t_id][j] != v2_id && mesh.tets[t_id][j] != v3_id)
                return j;
        }
        assert(false);
        return -1;
    }

    // The face on the other side of face j of tet t_id: the neighbouring tet and the local index
    // the shared face has there. Returns false and leaves both alone when face j is on the
    // boundary, which every caller has to handle.
    inline bool get_opp_face(const Mesh& mesh, int t_id, int j, int& opp_t_id, int& opp_j) {
        const int opp = get_opp_t_id(mesh, t_id, j);
        if (opp < 0)
            return false;
        opp_t_id = opp;
        opp_j    = get_local_f_id(opp,
                                  mesh.tets[t_id][(j + 1) % 4],
                                  mesh.tets[t_id][(j + 2) % 4],
                                  mesh.tets[t_id][(j + 3) % 4],
                                  mesh);
        return true;
    }

    // The normal of face j of tet t_id, pointing away from the corner the face is opposite to.
    // Not normalised: the callers either normalise it or only look at the sign of a dot product.
    inline Vector3 get_face_normal(const Mesh& mesh, int t_id, int j) {
        const Vector3& v1 = tet_pos(mesh, t_id, (j + 1) % 4);
        const Vector3& v2 = tet_pos(mesh, t_id, (j + 2) % 4);
        const Vector3& v3 = tet_pos(mesh, t_id, (j + 3) % 4);
        if (Predicates::orient_3d(v1, v2, v3, tet_pos(mesh, t_id, j)) ==
            Predicates::ORI_POSITIVE)
            return (v2 - v1).cross(v3 - v1);
        return (v3 - v1).cross(v2 - v1);
    }

    // A face a local operation has just created carries none of the old marks.
    inline void clear_face_marks(MeshTet& t, int j) {
        t.is_surface_fs[j] = NOT_SURFACE;
        t.surface_tags[j]  = NO_SURFACE_TAG;
        t.is_bbox_fs[j]    = NOT_BBOX;
    }

    // The two endpoints in ascending order, so an edge has one spelling wherever it is collected.
    inline std::array<int, 2> sorted_edge(int a, int b) {
        return a < b ? std::array<int, 2>{{a, b}} : std::array<int, 2>{{b, a}};
    }

    // The three corners of face j of a tet, in the local order that winds away from corner j.
    template<typename Tet>
    std::array<int, 3> face_corners(const Tet& t, int j) {
        return {{t[(j + 1) % 4], t[(j + 2) % 4], t[(j + 3) % 4]}};
    }

    // The same three, ascending, so a face has one spelling wherever it is collected.
    template<typename Tet>
    std::array<int, 3> sorted_face(const Tet& t, int j) {
        std::array<int, 3> f = face_corners(t, j);
        std::sort(f.begin(), f.end());
        return f;
    }

    // The six edges of a tet, in the order every caller has always pushed them.
    template<typename Tet>
    void push_tet_edges(const Tet& t, std::vector<std::array<int, 2>>& out) {
        for (int j = 0; j < 3; j++) {
            out.push_back(sorted_edge(t[0], t[j + 1]));
            out.push_back(sorted_edge(t[j + 1], t[(j + 1) % 3 + 1]));
        }
    }

    // The three edges of a triangle, likewise.
    template<typename Tri>
    void push_tri_edges(const Tri& f, std::vector<std::array<int, 2>>& out) {
        for (int j = 0; j < 3; j++)
            out.push_back(sorted_edge(f[j], f[(j + 1) % 3]));
    }

    void get_all_edges(const Mesh& mesh, std::vector<std::array<int, 2>>& edges);

    Scalar get_edge_length_2(const Mesh& mesh, int v1_id, int v2_id);

    // What the split and collapse thresholds are scaled by for a given edge: the sizing field
    // averaged over its two ends.
    inline Scalar avg_sizing_scalar(const Mesh& mesh, int v1_id, int v2_id) {
        return (mesh.tet_vertices[v1_id].sizing_scalar + mesh.tet_vertices[v2_id].sizing_scalar) / 2;
    }

    Scalar get_quality(const Mesh& mesh, const MeshTet& t);
    Scalar get_quality(const Mesh& mesh, int t_id);
    // The largest quality among the given tets, 0 when there are none.
    Scalar get_max_quality(const Mesh& mesh, const std::vector<int>& t_ids);
    Scalar get_quality(const MeshVertex& v0, const MeshVertex& v1, const MeshVertex& v2, const MeshVertex& v3);
    Scalar get_quality(const Vector3& v0, const Vector3& v1, const Vector3& v2, const Vector3& v3);

    bool is_inverted(const Mesh& mesh, int t_id);
    bool is_inverted(const Mesh& mesh, const MeshTet& t);
    bool is_inverted(const Mesh& mesh, int t_id, int j, const Vector3& new_p);
    bool is_inverted(const MeshVertex& v0, const MeshVertex& v1, const MeshVertex& v2, const MeshVertex& v3);
    bool is_inverted(const Vector3& v0, const Vector3& v1, const Vector3& v2, const Vector3& v3);

    bool is_out_envelope(Mesh& mesh, int v_id, const Vector3& new_pos, const AABBWrapper& tree);
    bool is_out_boundary_envelope(const Mesh& mesh, int v_id, const Vector3& new_pos, const AABBWrapper& tree);
    void sample_triangle(const std::array<Vector3, 3>& vs, std::vector<geo::vec3>& ps, Scalar sampling_dist);
    bool sample_triangle_and_check_is_out(const std::array<Vector3, 3>& vs, Scalar sampling_dist,
            Scalar eps_2, const AABBWrapper& tree, geo::index_t& prev_facet);

    bool is_bbox_edge(const Mesh& mesh, int v1_id, int v2_id, const std::vector<int>& n12_t_ids);
    bool is_surface_edge(const Mesh& mesh, int v1_id, int v2_id, const std::vector<int>& n12_t_ids);
    bool is_boundary_edge(const Mesh& mesh, int v1_id, int v2_id, const AABBWrapper& tree);
    bool is_valid_edge(const Mesh& mesh, int v1_id, int v2_id);
    bool is_valid_edge(const Mesh& mesh, int v1_id, int v2_id, const std::vector<int>& n12_t_ids);

    bool is_isolate_surface_point(const Mesh& mesh, int v_id);
    bool is_point_out_boundary_envelope(const Mesh& mesh, const Vector3& p, const AABBWrapper& tree);

    void get_new_tet_slots(Mesh& mesh, int n, std::vector<int>& new_conn_tets);
    int get_new_vertex_slot(Mesh& mesh, const MeshVertex& v);
    void split_tets_at(Mesh& mesh, int v_id, int v1_id, int v2_id,
                       const std::vector<int>& old_t_ids, std::vector<int>& new_t_ids);

    inline Scalar get_area(const Vector3& a, const Vector3& b, const Vector3& c) {
        return tri_normal(a, b, c).norm();
    }

    template<typename T>
    void vector_unique(std::vector<T>& v){
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
    }
    template<typename T>
    bool vector_erase(std::vector<T>& v, const T& t){
        auto it = std::find(v.begin(), v.end(), t);
        if(it == v.end())
            return false;
        v.erase(it);
        return true;
    }
    template<typename T>
    bool vector_contains(const std::vector<T>& v, const T& t){
        return std::find(v.begin(), v.end(), t) != v.end();
    }

    // The edges of the given tets, each once, appended to out. This is what an operation re-queues
    // after it rebuilds a patch.
    template<typename Tets>
    void collect_tet_edges(const Tets& ts, std::vector<std::array<int, 2>>& out) {
        out.reserve(out.size() + ts.size() * 6);
        for (const auto& t : ts)
            push_tet_edges(t, out);
        vector_unique(out);
    }

    inline void collect_tet_edges(const Mesh& mesh, const std::vector<int>& t_ids,
                                  std::vector<std::array<int, 2>>& out) {
        out.reserve(out.size() + t_ids.size() * 6);
        for (int t_id : t_ids)
            push_tet_edges(mesh.tets[t_id], out);
        vector_unique(out);
    }

    // The corners of the given tets, each once, ascending.
    inline void collect_tet_vertices(const Mesh& mesh, const std::vector<int>& t_ids,
                                     std::vector<int>& out) {
        out.reserve(out.size() + t_ids.size() * 4);
        for (int t_id : t_ids) {
            for (int j = 0; j < 4; j++)
                out.push_back(mesh.tets[t_id][j]);
        }
        vector_unique(out);
    }

    void set_intersection(const std::unordered_set<int>& s1, const std::unordered_set<int>& s2, std::vector<int>& v);
    void set_intersection(const std::vector<int>& s1, const std::vector<int>& s2, std::vector<int>& v);

    // The tets that share the face with these three corners: two of them, or one on the boundary.
    void get_face_tets(const Mesh& mesh, int v1_id, int v2_id, int v3_id, std::vector<int>& t_ids);

    class ElementInQueue{
    public:
        std::array<int, 2> v_ids;
        Scalar weight;

        ElementInQueue(){}
        ElementInQueue(const std::array<int, 2>& ids, Scalar w): v_ids(ids), weight(w){}
    };
    struct cmp_l {
        bool operator()(const ElementInQueue &e1, const ElementInQueue &e2) {
            if (e1.weight == e2.weight)
                return e1.v_ids > e2.v_ids;
            return e1.weight < e2.weight;
        }
    };
    struct cmp_s {
        bool operator()(const ElementInQueue &e1, const ElementInQueue &e2) {
            if (e1.weight == e2.weight)
                return e1.v_ids < e2.v_ids;
            return e1.weight > e2.weight;
        }
    };

    // The same edge can sit in the queue more than once. After taking one, drop the copies behind
    // it.
    template<typename Queue>
    void pop_duplicates(Queue &queue, const std::array<int, 2> &v_ids) {
        while (!queue.empty() && queue.top().v_ids == v_ids)
            queue.pop();
    }

    Scalar AMIPS_energy(const std::array<Scalar, 12>& T);
    void AMIPS_jacobian(const std::array<Scalar, 12>& T, Vector3& result_0);
    void AMIPS_hessian(const std::array<Scalar, 12>& T, Matrix3& result_0);
}

#endif //FLOATTETWILD_LOCALOPERATIONS_H
