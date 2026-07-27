#include "geo_mesh.h"
#include "geo_permutation.h"

namespace floatTetWild {
namespace geo {

    void MeshVertices::permute_elements(vector<index_t>& permutation) {
        // geogram keeps the coordinates in a vertex attribute, so its permute_elements() moves
        // them along with everything else bound to a vertex. Here they are a plain array, so
        // they are moved explicitly first.
        {
            std::vector<double> reordered(points_.size());
            for(index_t v = 0; v < permutation.size(); ++v) {
                const double* from = point_ptr(permutation[v]);
                std::copy(from, from + DIMENSION, reordered.begin() + v * DIMENSION);
            }
            points_.swap(reordered);
        }
        attributes_.apply_permutation(permutation);

        Permutation::invert(permutation);

        MeshFacetCorners& facet_corners = mesh_.facet_corners;
        for(index_t c = 0; c < facet_corners.nb(); ++c) {
            facet_corners.set_vertex(c, permutation[facet_corners.vertex(c)]);
        }
    }

    index_t MeshFacets::create_facets(
        index_t nb_facets, index_t nb_vertices_per_polygon
    ) {
        if(nb_vertices_per_polygon != 3) {
            is_not_simplicial();
        }

        index_t first_facet = nb();
        index_t co = facet_corners_.nb();
        facet_corners_.create_sub_elements(nb_facets * nb_vertices_per_polygon);

        nb_ += nb_facets;
        attributes_.resize(nb_);
        if(!is_simplicial_) {
            facet_ptr_.resize(nb_ + 1);
            for(index_t f = first_facet; f <= first_facet + nb_facets; ++f) {
                facet_ptr_[f] = co;
                co += nb_vertices_per_polygon;
            }
        }
        return first_facet;
    }

    void MeshFacets::is_not_simplicial() {
        if(is_simplicial_) {
            is_simplicial_ = false;
            facet_ptr_.resize(nb() + 1);
            for(index_t f = 0; f < facet_ptr_.size(); ++f) {
                facet_ptr_[f] = 3 * f;
            }
        }
    }

    void MeshFacets::assign_triangle_mesh(vector<index_t>& triangle_vertex_index) {
        index_t nb_triangles = index_t(triangle_vertex_index.size() / 3);
        is_simplicial_ = true;
        facet_ptr_.clear();
        nb_ = nb_triangles;
        facet_corners_.corner_vertex_.assign(
            triangle_vertex_index.begin(), triangle_vertex_index.end()
        );
        facet_corners_.corner_adjacent_facet_.assign(nb_triangles * 3, NO_FACET);
        // geogram zeroes the attributes here rather than carrying them across the retriangulation,
        // so a facet tag set before mesh_repair() does not survive it.
        attributes_.clear(true);
        attributes_.resize(nb_);
        facet_corners_.attributes_.clear(true);
        facet_corners_.attributes_.resize(nb_triangles * 3);
    }

    void MeshFacets::triangulate() {
        if(is_simplicial_) {
            return;
        }
        index_t nb_triangles = 0;
        for(index_t f = 0; f < nb(); f++) {
            nb_triangles += (nb_vertices(f) - 2);
        }
        vector<index_t> new_corner_vertex_index;
        new_corner_vertex_index.reserve(nb_triangles * 3);
        for(index_t f = 0; f < nb(); f++) {
            index_t v0 = facet_corners_.vertex(corners_begin(f));
            for(index_t c = corners_begin(f) + 1; c + 1 < corners_end(f); ++c) {
                new_corner_vertex_index.push_back(v0);
                new_corner_vertex_index.push_back(facet_corners_.vertex(c));
                new_corner_vertex_index.push_back(facet_corners_.vertex(c + 1));
            }
        }
        assign_triangle_mesh(new_corner_vertex_index);
    }

    namespace {
        index_t next_corner_around_facet(
            const MeshFacets& facets, index_t f, index_t c
        ) {
            return c + 1 == facets.corners_end(f) ? facets.corners_begin(f) : c + 1;
        }

        index_t prev_corner_around_facet(
            const MeshFacets& facets, index_t f, index_t c
        ) {
            return c == facets.corners_begin(f) ? facets.corners_end(f) - 1 : c - 1;
        }
    }

    void MeshFacets::connect() {
        if(nb() == 0) {
            return;
        }

        index_t c_begin = corners_begin(0);
        index_t c_end = corners_end(nb() - 1);
        index_t v_begin = 0;
        index_t v_end = mesh_.vertices.nb();

        // Gives for each corner the facet incident to it
        // (or use c/3 if the surface is triangulated).
        vector<index_t> c2f;
        if(!is_simplicial_) {
            c2f.assign(c_end - c_begin, NO_FACET);
            for(index_t f = 0; f < nb(); ++f) {
                for(index_t c = corners_begin(f); c < corners_end(f); ++c) {
                    c2f[c - c_begin] = f;
                }
            }
        }

        for(index_t c = c_begin; c < c_end; ++c) {
            facet_corners_.set_adjacent_facet(c, NO_FACET);
        }

        // Gives for each vertex a corner incident to it.
        vector<index_t> v2c(v_end - v_begin, NO_CORNER);

        // Chains the corners around each vertex.
        vector<index_t> next_corner_around_vertex(c_end - c_begin, NO_CORNER);

        // Step 1: chain corners around vertices and compute v2c
        for(index_t f = 0; f < nb(); ++f) {
            for(index_t c = corners_begin(f); c < corners_end(f); ++c) {
                index_t v = facet_corners_.vertex(c);
                next_corner_around_vertex[c - c_begin] = v2c[v - v_begin];
                v2c[v - v_begin] = c;
            }
        }

        // Step 2: connect
        for(index_t f1 = 0; f1 < nb(); ++f1) {
            for(index_t c1 = corners_begin(f1); c1 < corners_end(f1); ++c1) {
                if(facet_corners_.adjacent_facet(c1) == NO_FACET) {

                    index_t nb_candidates = 0;
                    index_t c_candidate = NO_CORNER;

                    index_t v1 = facet_corners_.vertex(c1);
                    index_t v2 = facet_corners_.vertex(
                        next_corner_around_facet(*this, f1, c1)
                    );

                    //   Traverse all the corners c2 incident to v1, and
                    // find among them the one(s) that is opposite to c1
                    for(
                        index_t c2 = v2c[v1 - v_begin];
                        c2 != NO_CORNER;
                        c2 = next_corner_around_vertex[c2 - c_begin]
                    ) {
                        if(c2 != c1) {
                            index_t f2 = is_simplicial_ ? c2 / 3 : c2f[c2 - c_begin];
                            index_t c2_prev = prev_corner_around_facet(*this, f2, c2);

                            index_t v3 = facet_corners_.vertex(c2);
                            index_t v4 = facet_corners_.vertex(c2_prev);

                            geo_assert(v1 == v3);

                            if(
                                v4 == v2 && (
                                    facet_corners_.adjacent_facet(c2_prev) == NO_FACET
                                )
                            ) {
                                c_candidate = c2_prev;
                                ++nb_candidates;
                            }
                        }
                    }
                    // If there were more than 1 candidate, do not connect.
                    if(nb_candidates == 1) {
                        index_t c2 = c_candidate;
                        index_t f2 = is_simplicial_ ? (c2 / 3) : c2f[c2 - c_begin];
                        facet_corners_.set_adjacent_facet(c1, f2);
                        facet_corners_.set_adjacent_facet(c2, f1);
                    }
                }
            }
        }
    }

    void MeshFacets::clear(bool keep_attributes, bool keep_memory) {
        nb_ = 0;
        is_simplicial_ = true;
        facet_ptr_.clear();
        if(!keep_memory) {
            facet_ptr_.shrink_to_fit();
        }
        attributes_.clear(keep_attributes);
        facet_corners_.clear(keep_attributes, keep_memory);
    }

    void MeshFacets::permute_elements(vector<index_t>& permutation) {
        attributes_.apply_permutation(permutation);

        std::vector<index_t>& corner_vertex = facet_corners_.corner_vertex_;
        std::vector<index_t>& corner_adjacent_facet = facet_corners_.corner_adjacent_facet_;

        if(facet_corners_.attributes().size() != 0) {
            vector<index_t> facet_corners_permutation;
            facet_corners_permutation.reserve(facet_corners_.nb());

            for(index_t new_f = 0; new_f < nb(); ++new_f) {
                index_t old_f = permutation[new_f];
                for(
                    index_t old_c = corners_begin(old_f);
                    old_c < corners_end(old_f); ++old_c) {
                    facet_corners_permutation.push_back(old_c);
                }
            }

            facet_corners_.attributes().apply_permutation(facet_corners_permutation);
        }

        if(is_simplicial_) {
            // If the surface is triangulated,
            // everything can be done in-place (great !!)

            Permutation::apply(
                corner_vertex.data(), permutation, index_t(sizeof(index_t) * 3)
            );

            Permutation::apply(
                corner_adjacent_facet.data(), permutation, index_t(sizeof(index_t) * 3)
            );

            Permutation::invert(permutation);

            for(index_t c = 0; c < corner_adjacent_facet.size(); ++c) {
                if(corner_adjacent_facet[c] != NO_FACET) {
                    corner_adjacent_facet[c] = permutation[corner_adjacent_facet[c]];
                }
            }

        } else {

            {
                std::vector<index_t> new_corner_vertex;
                new_corner_vertex.reserve(corner_vertex.size());
                std::vector<index_t> new_corner_adjacent_facet;
                new_corner_adjacent_facet.reserve(corner_adjacent_facet.size());
                std::vector<index_t> new_facet_ptr;
                new_facet_ptr.reserve(nb() + 1);

                new_facet_ptr.push_back(0);
                for(index_t new_f = 0; new_f < nb(); ++new_f) {
                    index_t old_f = permutation[new_f];
                    for(
                        index_t old_c = corners_begin(old_f);
                        old_c < corners_end(old_f); ++old_c
                    ) {
                        new_corner_vertex.push_back(facet_corners_.vertex(old_c));
                        new_corner_adjacent_facet.push_back(
                            facet_corners_.adjacent_facet(old_c)
                        );
                    }
                    new_facet_ptr.push_back(
                        new_facet_ptr[new_facet_ptr.size() - 1] + nb_vertices(old_f)
                    );
                }

                corner_vertex.swap(new_corner_vertex);
                corner_adjacent_facet.swap(new_corner_adjacent_facet);
                facet_ptr_.swap(new_facet_ptr);
            }

            Permutation::invert(permutation);

            for(index_t c = 0; c < corner_adjacent_facet.size(); ++c) {
                if(corner_adjacent_facet[c] != NO_FACET) {
                    corner_adjacent_facet[c] = permutation[corner_adjacent_facet[c]];
                }
            }
        }
    }

    double bbox_diagonal(const Mesh& M) {
        double xyzmin[3];
        double xyzmax[3];
        for(index_t c = 0; c < 3; ++c) {
            xyzmin[c] = Numeric::max_float64();
            xyzmax[c] = Numeric::min_float64();
        }
        for(index_t v = 0; v < M.vertices.nb(); ++v) {
            const double* p = M.vertices.point_ptr(v);
            for(index_t c = 0; c < 3; ++c) {
                xyzmin[c] = std::min(xyzmin[c], p[c]);
                xyzmax[c] = std::max(xyzmax[c], p[c]);
            }
        }
        return ::sqrt(
            geo_sqr(xyzmax[0] - xyzmin[0]) +
            geo_sqr(xyzmax[1] - xyzmin[1]) +
            geo_sqr(xyzmax[2] - xyzmin[2])
        );
    }
}
}
