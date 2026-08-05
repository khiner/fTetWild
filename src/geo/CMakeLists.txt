# The slice of geogram that fTetWild reaches, copied into floatTetWild::geo rather than
# reimplemented: Delaunay's output depends on vertex index order and on the symbolic
# perturbation that breaks ties, so a rewrite that produced the same vertex set in a
# different order would still change the mesh. See LICENSE.geogram.
set(SOURCES
	geo_basic.h
	geo_parallel.h
	geo_geometry.h

	geo_mesh.h
	geo_mesh_reorder.h
	geo_mesh_reorder.cpp

	geo_multi_precision.h
	geo_multi_precision.cpp

	geo_kd_tree.h
	geo_kd_tree.cpp

	geo_delaunay_3d.h
	geo_delaunay_3d.cpp

	geo_predicates.h
	geo_predicates.cpp
)

prepend_current_path(SOURCES)
float_tetwild_copy_headers(${SOURCES})
float_tetwild_set_source_group(${SOURCES})
target_sources(FloatTetwild PRIVATE ${SOURCES})
