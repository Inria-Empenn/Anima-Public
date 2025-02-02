cmake_path(GET CMAKE_CURRENT_SOURCE_DIR STEM DIRNAME)
set(PROJECT_PREFIX anima${DIRNAME})

add_subdirectory(Estimator)
add_subdirectory(Generator)
add_subdirectory(Reductor)
