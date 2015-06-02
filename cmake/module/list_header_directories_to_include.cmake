macro(list_header_directories_to_include
  project_name
  headers
  )

################################################################################
#
# Usage: list_header_directories_to_include(project_name, header_path1, header_path2 ...)
# list all the different path found in all given header_path, and add them to
# ${project_name}_INCLUDE_DIRS
#
################################################################################

## #############################################################################
## List header directories to include
## #############################################################################

foreach(headers_root ${ARGV})
    file(GLOB_RECURSE
	HEADERS
        ${headers_root}/*.h
        ${headers_root}/*.hpp
        ${headers_root}/*.hxx
    )

    foreach(headers_path ${HEADERS})
        get_filename_component(dir_path ${headers_path} PATH)
        set(${project_name}_INCLUDE_DIRS
            ${dir_path}
            ${${project_name}_INCLUDE_DIRS}
        )
    endforeach()
endforeach()

list(REMOVE_DUPLICATES ${project_name}_INCLUDE_DIRS)

## #############################################################################
## add it to the parent scope for future usage.
## #############################################################################

set(${project_name}_INCLUDE_DIRS
    ${${project_name}_INCLUDE_DIRS}
#    PARENT_SCOPE
    )

endmacro()
