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
)

prepend_current_path(SOURCES)
float_tetwild_copy_headers(${SOURCES})
float_tetwild_set_source_group(${SOURCES})
target_sources(FloatTetwild PRIVATE ${SOURCES})

# File I/O is not part of the core library: the binary owns it. The headers are still copied so
# that main.cpp can include them, only the translation units move.
set(IO_SOURCES
	MeshIO.hpp
	MeshIO.cpp
	SurfaceMeshLoad.hpp
	SurfaceMeshLoad.cpp
)
prepend_current_path(IO_SOURCES)
float_tetwild_copy_headers(${IO_SOURCES})
float_tetwild_set_source_group(${IO_SOURCES})
if(TARGET FloatTetwild_bin)
	target_sources(FloatTetwild_bin PRIVATE ${IO_SOURCES})
endif()

################################################################################
# Subfolders
################################################################################

add_subdirectory(external)
