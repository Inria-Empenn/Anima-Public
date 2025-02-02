cmake_path(GET CMAKE_CURRENT_SOURCE_DIR STEM DIRNAME)
set(PROJECT_PREFIX ${PROJECT_PREFIX}${DIRNAME})

add_subdirectory(ByKMeansAlgorithm)
add_subdirectory(ByNyulAlgorithm)
