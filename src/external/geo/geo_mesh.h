// A stand-in for the slice of geogram/mesh/mesh.h that fTetWild and the vendored spatial sort,
// loaders and AABB use: vertices, facets, facet corners and their attributes.
//
// Reimplemented rather than copied, unlike its neighbours here. A container addresses elements by
// index and has no comparisons of its own, so there is no tie-breaking to preserve; what it does
// have to match is the element order geogram's operations leave behind, which is why
// permute_elements(), triangulate() and connect() follow geogram's implementations step for step.
//
// geogram's cells are not here: fTetWild's only use of them was Mesh::partition(), which has no
// live caller.

#pragma once

#include "geo_attributes.h"
#include "geo_geometry.h"

namespace floatTetWild {
namespace geo {

    static const index_t NO_VERTEX = index_t(-1);
    static const index_t NO_FACET = index_t(-1);
    static const index_t NO_CORNER = index_t(-1);

    class Mesh;

    /**
     * \brief The vertices of a Mesh, stored as one contiguous array of coordinates.
     */
    class MeshVertices {
    public:
        explicit MeshVertices(Mesh& mesh) : mesh_(mesh) {
        }

        index_t nb() const {
            return index_t(points_.size() / DIMENSION);
        }

        /**
         * \brief geogram supports other dimensions. fTetWild only ever builds 3d meshes.
         */
        index_t dimension() const {
            return DIMENSION;
        }

        double* point_ptr(index_t v) {
            geo_debug_assert(v < nb());
            return points_.data() + v * DIMENSION;
        }

        const double* point_ptr(index_t v) const {
            geo_debug_assert(v < nb());
            return points_.data() + v * DIMENSION;
        }

        vec3& point(index_t v) {
            return *reinterpret_cast<vec3*>(point_ptr(v));
        }

        const vec3& point(index_t v) const {
            return *reinterpret_cast<const vec3*>(point_ptr(v));
        }

        index_t create_vertex() {
            return create_vertices(1);
        }

        index_range::iterator begin() const {
            return index_range(0, nb()).begin();
        }

        index_range::iterator end() const {
            return index_range(0, nb()).end();
        }

        index_t create_vertices(index_t nb_to_create) {
            index_t first = nb();
            points_.resize((first + nb_to_create) * DIMENSION, 0.0);
            attributes_.resize(nb());
            return first;
        }

        void clear(bool keep_attributes = true, bool keep_memory = false) {
            points_.clear();
            if(!keep_memory) {
                points_.shrink_to_fit();
            }
            attributes_.clear(keep_attributes);
        }

        /**
         * \brief Reorders the vertices so that vertex \p k becomes the vertex that was at
         *  \p permutation[k], and remaps everything that refers to a vertex.
         * \details \p permutation is inverted in place, as in geogram.
         */
        void permute_elements(vector<index_t>& permutation);

        AttributesManager& attributes() {
            return attributes_;
        }

    private:
        static const index_t DIMENSION = 3;

        Mesh& mesh_;
        std::vector<double> points_;
        AttributesManager attributes_;
    };

    /**
     * \brief The corners of the facets of a Mesh: for each one, the vertex it points at and the
     *  facet on the other side of the edge that starts there.
     */
    class MeshFacetCorners {
    public:
        MeshFacetCorners() {
        }

        index_t nb() const {
            return index_t(corner_vertex_.size());
        }

        index_t vertex(index_t c) const {
            geo_debug_assert(c < nb());
            return corner_vertex_[c];
        }

        void set_vertex(index_t c, index_t v) {
            geo_debug_assert(c < nb());
            corner_vertex_[c] = v;
        }

        index_t adjacent_facet(index_t c) const {
            geo_debug_assert(c < nb());
            return corner_adjacent_facet_[c];
        }

        void set_adjacent_facet(index_t c, index_t f) {
            geo_debug_assert(c < nb());
            corner_adjacent_facet_[c] = f;
        }

        AttributesManager& attributes() {
            return attributes_;
        }

        index_range::iterator begin() const {
            return index_range(0, nb()).begin();
        }

        index_range::iterator end() const {
            return index_range(0, nb()).end();
        }

    private:
        friend class MeshVertices;
        friend class MeshFacets;

        index_t create_sub_elements(index_t nb_to_create) {
            index_t first = nb();
            corner_vertex_.resize(first + nb_to_create, NO_VERTEX);
            corner_adjacent_facet_.resize(first + nb_to_create, NO_FACET);
            attributes_.resize(nb());
            return first;
        }

        void clear(bool keep_attributes, bool keep_memory) {
            corner_vertex_.clear();
            corner_adjacent_facet_.clear();
            if(!keep_memory) {
                corner_vertex_.shrink_to_fit();
                corner_adjacent_facet_.shrink_to_fit();
            }
            attributes_.clear(keep_attributes);
        }

        std::vector<index_t> corner_vertex_;
        std::vector<index_t> corner_adjacent_facet_;
        AttributesManager attributes_;
    };

    /**
     * \brief The facets of a Mesh.
     * \details As in geogram, a mesh whose facets are all triangles keeps them in a simplicial
     *  representation with the corners of facet f at 3f, 3f+1, 3f+2 and no facet pointer array.
     */
    class MeshFacets {
    public:
        explicit MeshFacets(MeshFacetCorners& facet_corners) :
            facet_corners_(facet_corners) {
        }

        index_t nb() const {
            return nb_;
        }

        bool are_simplices() const {
            return is_simplicial_;
        }

        index_t corners_begin(index_t f) const {
            geo_debug_assert(f < nb());
            return is_simplicial_ ? 3 * f : facet_ptr_[f];
        }

        index_t corners_end(index_t f) const {
            geo_debug_assert(f < nb());
            return is_simplicial_ ? 3 * (f + 1) : facet_ptr_[f + 1];
        }

        index_t nb_vertices(index_t f) const {
            return corners_end(f) - corners_begin(f);
        }

        index_range corners(index_t f) const {
            return index_range(corners_begin(f), corners_end(f));
        }

        index_range::iterator begin() const {
            return index_range(0, nb()).begin();
        }

        index_range::iterator end() const {
            return index_range(0, nb()).end();
        }

        index_t vertex(index_t f, index_t lv) const {
            geo_debug_assert(lv < nb_vertices(f));
            return facet_corners_.vertex(corners_begin(f) + lv);
        }

        void set_vertex(index_t f, index_t lv, index_t v) {
            geo_debug_assert(lv < nb_vertices(f));
            facet_corners_.set_vertex(corners_begin(f) + lv, v);
        }

        index_t adjacent(index_t f, index_t le) const {
            geo_debug_assert(le < nb_vertices(f));
            return facet_corners_.adjacent_facet(corners_begin(f) + le);
        }

        void set_adjacent(index_t f, index_t le, index_t f2) {
            geo_debug_assert(le < nb_vertices(f));
            facet_corners_.set_adjacent_facet(corners_begin(f) + le, f2);
        }

        /**
         * \brief Creates a contiguous chunk of facets with the same number of vertices each.
         * \return the index of the first facet
         */
        index_t create_facets(index_t nb_facets, index_t nb_vertices_per_polygon);

        index_t create_triangles(index_t nb_triangles) {
            return create_facets(nb_triangles, 3);
        }

        /**
         * \brief Creates one facet.
         * \return the index of the created facet
         */
        index_t create_polygon(index_t nb_facet_vertices) {
            return create_facets(1, nb_facet_vertices);
        }



        void clear(bool keep_attributes = true, bool keep_memory = false);

        /**
         * \brief Reorders the facets so that facet \p k becomes the facet that was at
         *  \p permutation[k].
         * \details \p permutation is inverted in place, as in geogram.
         */
        void permute_elements(vector<index_t>& permutation);

        AttributesManager& attributes() {
            return attributes_;
        }

    private:
        friend class Mesh;

        /**
         * \brief Switches to the general representation, materialising facet_ptr_.
         */
        void is_not_simplicial();

        /**
         * \brief Replaces the facets with the triangles given as a flat list of vertex indices.
         */
        void assign_triangle_mesh(vector<index_t>& triangle_vertex_index);

        MeshFacetCorners& facet_corners_;
        index_t nb_ = 0;
        bool is_simplicial_ = true;
        std::vector<index_t> facet_ptr_;
        AttributesManager attributes_;
    };

    /**
     * \brief A surface mesh.
     */
    class Mesh {
    public:
        Mesh() : vertices(*this), facets(facet_corners) {
        }

        Mesh(const Mesh&) = delete;
        Mesh& operator= (const Mesh&) = delete;

        void clear(bool keep_attributes = true, bool keep_memory = false) {
            vertices.clear(keep_attributes, keep_memory);
            facets.clear(keep_attributes, keep_memory);
        }

        MeshVertices vertices;
        MeshFacetCorners facet_corners;
        MeshFacets facets;
    };

    /**
     * \brief Returns the length of the diagonal of the bounding box of \p M.
     */
    double bbox_diagonal(const Mesh& M);
}
}
