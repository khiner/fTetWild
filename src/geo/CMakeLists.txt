# The slice of geogram that fTetWild reaches, copied into floatTetWild::geo rather than
# reimplemented: Delaunay's output depends on vertex index order and on the symbolic
# perturbation that breaks ties, so a rewrite that produced the same vertex set in a
# different order would still change the mesh. See LICENSE.geogram.
target_sources(FloatTetwild PRIVATE
	geo_mesh_reorder.cpp
	geo_multi_precision.cpp
	geo_kd_tree.cpp
	geo_delaunay_3d.cpp
	geo_predicates.cpp
)
