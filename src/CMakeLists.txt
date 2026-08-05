set(SOURCES
	# The mesher, from the top down.
	FloatTetwild.h
	FloatTetwild.cpp
	Simplification.h
	Simplification.cpp
	FloatTetDelaunay.h
	FloatTetDelaunay.cpp
	TriangleInsertion.h
	TriangleInsertion.cpp
	CutMesh.h
	CutMesh.cpp
	MeshImprovement.h
	MeshImprovement.cpp

	# The local operations the optimization loop runs.
	LocalOperations.h
	LocalOperations.cpp
	EdgeSplitting.h
	EdgeSplitting.cpp
	EdgeCollapsing.h
	EdgeCollapsing.cpp
	EdgeSwapping.h
	EdgeSwapping.cpp
	VertexSmoothing.h
	VertexSmoothing.cpp

	# Supporting pieces.
	AABBWrapper.h
	AABBWrapper.cpp
	auto_table.hpp
	auto_table.cpp
	CSGTreeParser.hpp
	CSGTreeParser.cpp
	intersections.h
	intersections.cpp
	Logger.hpp
	Logger.cpp
	ParallelFor.hpp
	ParallelFor.cpp
	Matrix.hpp
	Mesh.hpp
	Parameters.h
	Statistics.h
	Timer.h
	Types.hpp
)

# Code from other projects, kept under its own licence, listed in THIRD_PARTY.md and repeated in
# each file's header. No subdirectory because the include path is flat either way:
# float_tetwild_copy_headers puts every header in <floattetwild/> wherever it came from.
set(VENDORED_SOURCES
	# Shewchuk's exact predicates, public domain, and fTetWild's wrapper over them.
	predicates.c
	Predicates.hpp
	Predicates.cpp

	# geogram, BSD 3-clause. See geo/LICENSE.geogram.
	mesh_AABB.h
	mesh_AABB.cpp

	# Guigue-Devillers triangle-triangle overlap, from the JGT paper. No licence stated upstream.
	triangle_triangle_intersection.cpp

	# Peak RSS, CC BY 3.0. Read by the stats csv summary row in main.cpp.
	getRSS.c

	# libigl, MPL 2.0, the same licence as fTetWild.
	MeshCleanup.hpp
	MeshCleanup.cpp

	# The fast winding number, Side Effects Software, MIT. Runs on ParallelFor.hpp.
	FastWindingNumber.h
	FastWindingNumber.cpp
)

list(APPEND SOURCES ${VENDORED_SOURCES})

prepend_current_path(SOURCES)
float_tetwild_copy_headers(${SOURCES})
float_tetwild_set_source_group(${SOURCES})
target_sources(FloatTetwild PRIVATE ${SOURCES})

# File I/O is not part of the core library: the binary owns it. writeOBJ.h is vendored libigl,
# MPL 2.0, reached only from main.cpp.
set(IO_SOURCES
	MeshIO.hpp
	MeshIO.cpp
	writeOBJ.h
)
prepend_current_path(IO_SOURCES)
float_tetwild_copy_headers(${IO_SOURCES})
float_tetwild_set_source_group(${IO_SOURCES})
if(TARGET FloatTetwild_bin)
	target_sources(FloatTetwild_bin PRIVATE ${IO_SOURCES})
endif()

# Third-party code we do not edit, so its warnings are not ours to fix. TARGET_DIRECTORY because
# both targets are defined one level up.
if(NOT MSVC)
	set_source_files_properties(
		${CMAKE_CURRENT_SOURCE_DIR}/predicates.c
		${CMAKE_CURRENT_SOURCE_DIR}/getRSS.c
		${CMAKE_CURRENT_SOURCE_DIR}/mesh_AABB.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/geo/geo_mesh_reorder.cpp
		TARGET_DIRECTORY FloatTetwild PROPERTIES COMPILE_OPTIONS
		"-Wno-cast-qual;-Wno-deprecated-non-prototype;-Wno-strict-prototypes;-Wno-switch-default;-Wno-unused-function")
endif()

################################################################################
# Subfolders
################################################################################

add_subdirectory(geo)
