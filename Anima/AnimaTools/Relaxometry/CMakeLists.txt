################################################################################
# Here go the add_subdirectories, no code should be at the root of the project
################################################################################

cmake_path(GET CMAKE_CURRENT_SOURCE_DIR STEM DIRNAME  )
set(PROJECT_PREFIX anima${DIRNAME})

add_subdirectory(Estimation)
add_subdirectory(Simulation)
