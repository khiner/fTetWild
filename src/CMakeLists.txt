set(SOURCES
	# The mesher, from the top down.
	FloatTetwild.cpp
	Simplification.cpp
	TriangleInsertion.cpp
	MeshImprovement.cpp

	# The local operations the optimization loop runs.
	LocalOperations.cpp
	EdgeSplitting.cpp
	EdgeCollapsing.cpp
	EdgeSwapping.cpp
	VertexSmoothing.cpp

	# Supporting pieces.
	AABBWrapper.cpp
	auto_table.cpp
	CSGTreeParser.cpp
	intersections.cpp
	Logger.cpp
	ParallelFor.cpp
)

# Code from other projects, kept under its own licence, listed in THIRD_PARTY.md and repeated in
# each file's header. The Guigue-Devillers triangle-triangle overlap routine, no licence stated
# upstream, sits at the end of intersections.cpp above.
set(VENDORED_SOURCES
	# Shewchuk's exact predicates, public domain, and fTetWild's wrapper over them.
	predicates.c
	Predicates.cpp

	# libigl, MPL 2.0, the same licence as fTetWild.
	MeshCleanup.cpp

	# The fast winding number, Side Effects Software, MIT. Runs on ParallelFor.hpp.
	FastWindingNumber.cpp
)

target_sources(FloatTetwild PRIVATE ${SOURCES} ${VENDORED_SOURCES})

# File I/O is not part of the core library: the binary owns it.
target_sources(FloatTetwild_bin PRIVATE MeshIO.cpp)

# Third-party code we do not edit, so its warnings are not ours to fix. TARGET_DIRECTORY because
# both targets are defined one level up.
if(NOT MSVC)
	set_source_files_properties(
		${CMAKE_CURRENT_SOURCE_DIR}/predicates.c
		${CMAKE_CURRENT_SOURCE_DIR}/geo/geo_mesh_reorder.cpp
		TARGET_DIRECTORY FloatTetwild PROPERTIES COMPILE_OPTIONS
		"-Wno-cast-qual;-Wno-deprecated-non-prototype;-Wno-strict-prototypes;-Wno-switch-default;-Wno-unused-function")
endif()

add_subdirectory(geo)
