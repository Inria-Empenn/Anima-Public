
cmake_path(GET CMAKE_CURRENT_SOURCE_DIR STEM DIRNAME  )
set(PROJECT_PREFIX anima${DIRNAME})

################################################################################
# Here go the add_subdirectories, no code should be at the root of the project
################################################################################

add_subdirectory(Arithmetic)
add_subdirectory(CommonTools)
add_subdirectory(NonLocalBase)
add_subdirectory(StatisticalTests)
add_subdirectory(Statistics)
