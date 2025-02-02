cmake_path(GET CMAKE_CURRENT_SOURCE_DIR STEM DIRNAME)
set(PROJECT_PREFIX anima${DIRNAME})

add_subdirectory(BiasFieldCorrection)
add_subdirectory(Standardization)
