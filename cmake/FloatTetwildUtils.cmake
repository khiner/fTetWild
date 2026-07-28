# Copy float tetwild header files into the build directory
function(float_tetwild_copy_headers)
	foreach(filepath IN ITEMS ${ARGN})
		get_filename_component(filename "${filepath}" NAME)
		if(${filename} MATCHES ".*\.(hpp|h|ipp)$")
			configure_file(${filepath} ${PROJECT_BINARY_DIR}/include/floattetwild/${filename})
		endif()
	endforeach()
endfunction()

# Set source group for IDE like Visual Studio or XCode
function(float_tetwild_set_source_group)
	foreach(filepath IN ITEMS ${ARGN})
		get_filename_component(folderpath "${filepath}" DIRECTORY)
		get_filename_component(foldername "${folderpath}" NAME)
		source_group("${foldername}" FILES "${filepath}")
	endforeach()
endfunction()
