#include <floattetwild/AABBWrapper.h>
#include <floattetwild/LocalOperations.h>
#include <floattetwild/TriangleInsertion.h>
#include <floattetwild/geo_mesh_reorder.h>

namespace floatTetWild {
namespace {

// The envelope trees are built over segments, which geogram's AABB stores as degenerate triangles
// with the first corner repeated. An empty edge list still needs a mesh to query, so it gets one
// triangle collapsed to the origin.
void build_segment_mesh(geo::Mesh& mesh, const std::vector<std::array<Vector3, 2>>& segments)
{
    mesh.vertices.clear();
    mesh.facets.clear();

    if (segments.empty()) {
        mesh.vertices.create_vertices(1);
        mesh.vertices.point(0) = geo::vec3(0, 0, 0);
        mesh.facets.create_triangles(1);
        for (int lv = 0; lv < 3; lv++)
            mesh.facets.set_vertex(0, lv, 0);
        return;
    }

    mesh.vertices.create_vertices((int)segments.size() * 2);
    for (int i = 0; i < segments.size(); i++) {
        for (int j = 0; j < 2; j++) {
            geo::vec3& p = mesh.vertices.point(i * 2 + j);
            for (int k = 0; k < 3; k++)
                p[k] = segments[i][j][k];
        }
    }

    mesh.facets.create_triangles((int)segments.size());
    for (int i = 0; i < segments.size(); i++) {
        mesh.facets.set_vertex(i, 0, i * 2);
        mesh.facets.set_vertex(i, 1, i * 2);
        mesh.facets.set_vertex(i, 2, i * 2 + 1);
    }
}

}  // namespace

void AABBWrapper::init_b_mesh_and_tree(const std::vector<Vector3>&  input_vertices,
                                       const std::vector<Vector3i>& input_faces,
                                       Mesh&                        mesh)
{
    b_mesh.clear(false, false);

    std::vector<std::pair<std::array<int, 2>, std::vector<int>>> _;
    std::vector<bool>                                           _1;
    std::vector<std::array<int, 2>>                             b_edges;
    find_boundary_edges(input_vertices,
                        input_faces,
                        std::vector<bool>(input_faces.size(), true),
                        std::vector<bool>(input_faces.size(), true),
                        _,
                        _1,
                        b_edges);

    std::vector<std::array<Vector3, 2>> segments;
    segments.reserve(b_edges.size());
    for (const auto& e : b_edges)
        segments.push_back({{input_vertices[e[0]], input_vertices[e[1]]}});

    build_segment_mesh(b_mesh, segments);
    mesh_reorder(b_mesh, geo::MESH_ORDER_MORTON);
    b_tree = std::make_shared<MeshFacetsAABBWithEps>(b_mesh);

    if (b_edges.empty())
        mesh.is_closed = true;
}

void AABBWrapper::init_tmp_b_mesh_and_tree(const std::vector<Vector3>&            input_vertices,
                                           const std::vector<Vector3i>&           input_faces,
                                           const std::vector<std::array<int, 2>>& b_edges1,
                                           const Mesh&                            mesh,
                                           const std::vector<std::array<int, 2>>& b_edges2)
{
    std::vector<std::array<Vector3, 2>> segments;
    segments.reserve(b_edges1.size() + b_edges2.size());
    for (const auto& e : b_edges1)
        segments.push_back({{input_vertices[e[0]], input_vertices[e[1]]}});
    for (const auto& e : b_edges2)
        segments.push_back({{mesh.tet_vertices[e[0]].pos, mesh.tet_vertices[e[1]].pos}});

    build_segment_mesh(tmp_b_mesh, segments);
    mesh_reorder(tmp_b_mesh, geo::MESH_ORDER_MORTON);
    tmp_b_tree = std::make_shared<MeshFacetsAABBWithEps>(tmp_b_mesh);
}

}  // namespace floatTetWild
