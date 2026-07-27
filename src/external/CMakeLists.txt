if(FLOAT_TETWILD_WITH_EXACT_ENVELOPE )
	set(SOURCES
	predicates.c
	predicates.h
		Predicates.cpp
		Predicates.hpp

		mesh_AABB.h
		mesh_AABB.cpp



		triangle_triangle_intersection.cpp

		getRSS.c
		get_mem.h
		get_mem.cpp

		bfs_orient.h
		bfs_orient.cpp

		Timer.h
		round.h
		sortrows.h
		unique_rows.h
		remove_duplicate_vertices.h

		writeOBJ.h

		default_num_threads.h
		parallel_for.h
		FastWindingNumberForSoups.h
		fast_winding_number.h

		sort.h
		vertex_components.h
		orientable_patches.h

#			WindingNumber.h
#		FastWindingNumber.hpp
#		FastWindingNumber.cpp

		Rational.h
	)
else()
	set(SOURCES
			predicates.c
			predicates.h
			Predicates.cpp
			Predicates.hpp

			mesh_AABB.h
			mesh_AABB.cpp



			triangle_triangle_intersection.cpp

			getRSS.c
			get_mem.h
			get_mem.cpp

			bfs_orient.h
			bfs_orient.cpp

			Timer.h
			round.h
			sortrows.h
			unique_rows.h
			remove_duplicate_vertices.h

			writeOBJ.h

			default_num_threads.h
			parallel_for.h
			FastWindingNumberForSoups.h
			fast_winding_number.h

			sort.h
			vertex_components.h
			orientable_patches.h

			Rational.h
			)
endif()

prepend_current_path(SOURCES)
float_tetwild_copy_headers(${SOURCES})
float_tetwild_set_source_group(${SOURCES})
target_sources(FloatTetwild PRIVATE ${SOURCES})

# PyMesh .msh reader and writer, used only by MeshIO and by main.cpp for --bg-mesh, so they
# belong with the binary rather than the core library.
set(IO_SOURCES
	MshSaver.h
	MshSaver.cpp
	MshLoader.h
	MshLoader.cpp
	Exception.h
)
prepend_current_path(IO_SOURCES)
float_tetwild_copy_headers(${IO_SOURCES})
float_tetwild_set_source_group(${IO_SOURCES})
if(TARGET FloatTetwild_bin)
	target_sources(FloatTetwild_bin PRIVATE ${IO_SOURCES})
endif()
