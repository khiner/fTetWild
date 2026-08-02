set(SOURCES
	AABBWrapper.h
	AABBWrapper.cpp


	auto_table.hpp
	auto_table.cpp

	CSGTreeParser.hpp
	CSGTreeParser.cpp

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

	MeshImprovement.h
	MeshImprovement.cpp


	FloatTetwild.h
	FloatTetwild.cpp

	FloatTetDelaunay.h
	FloatTetDelaunay.cpp
		TriangleInsertion.cpp
		TriangleInsertion.h

		CutMesh.cpp
		CutMesh.h

	intersections.h
	intersections.cpp

	Simplification.h
	Simplification.cpp

	Mesh.cpp
	Mesh.hpp


	Logger.hpp
	Logger.cpp

	ParallelFor.hpp
	ParallelFor.cpp

	Parameters.h

	Matrix.hpp
	Matrix.ipp

	Types.hpp

	Statistics.h

	Timer.h
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

	# Peak RSS, CC BY 3.0. Read through get_peak_mem in Statistics.h.
	getRSS.c

	# libigl, MPL 2.0, the same licence as fTetWild.
	MeshCleanup.hpp
	MeshCleanup.cpp
	writeOBJ.h

	# The fast winding number: libigl's wrapper, MPL 2.0, over Side Effects' implementation, MIT.
	# parallel_for.h keeps igl::parallel_for's signature, so those call sites need no edit, but
	# runs on ParallelFor.hpp.
	parallel_for.h
	FastWindingNumberForSoups.h
	fast_winding_number.h
)

list(APPEND SOURCES ${VENDORED_SOURCES})

prepend_current_path(SOURCES)
float_tetwild_copy_headers(${SOURCES})
float_tetwild_set_source_group(${SOURCES})
target_sources(FloatTetwild PRIVATE ${SOURCES})

# File I/O is not part of the core library: the binary owns it. The headers are still copied so
# that main.cpp can include them, only the translation units move. MshSaver is PyMesh's .msh
# writer, reached only from MeshIO, so it belongs on this side too.
set(IO_SOURCES
	MeshIO.hpp
	MeshIO.cpp
	SurfaceMeshLoad.hpp
	SurfaceMeshLoad.cpp

	MshSaver.h
	MshSaver.cpp
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
		${CMAKE_CURRENT_SOURCE_DIR}/MshSaver.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/geo/geo_mesh_reorder.cpp
		TARGET_DIRECTORY FloatTetwild PROPERTIES COMPILE_OPTIONS
		"-Wno-cast-qual;-Wno-deprecated-non-prototype;-Wno-strict-prototypes;-Wno-switch-default;-Wno-unused-function")
endif()

################################################################################
# Subfolders
################################################################################

add_subdirectory(geo)
