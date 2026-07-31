# The slice of geogram that fTetWild reaches, copied into floatTetWild::geo rather than
# reimplemented: Delaunay's output depends on vertex index order and on the symbolic
# perturbation that breaks ties, so a rewrite that produced the same vertex set in a
# different order would still change the mesh. See LICENSE.geogram.
set(SOURCES
	geo_basic.h
	geo_determinant.h
	geo_vecg.h
	geo_geometry.h
	geo_geometry_nd.h

	geo_algorithm.h
	geo_parallel.h
	geo_permutation.h
	geo_attributes.h

	geo_mesh.h
	geo_mesh.cpp
	geo_mesh_reorder.h
	geo_mesh_reorder.cpp

	geo_multi_precision.h
	geo_multi_precision.cpp

	geo_nn_search.h
	geo_nn_search.cpp
	geo_kd_tree.h
	geo_kd_tree.cpp

	geo_cavity.h
	geo_delaunay.h
	geo_delaunay_3d.h
	geo_delaunay_3d.cpp

	predicates_orient3d.h
	predicates_side4h.h
	geo_predicates.h
	geo_predicates.cpp
)

prepend_current_path(SOURCES)
float_tetwild_copy_headers(${SOURCES})
float_tetwild_set_source_group(${SOURCES})
target_sources(FloatTetwild PRIVATE ${SOURCES})
