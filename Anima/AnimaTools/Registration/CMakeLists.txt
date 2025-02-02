cmake_path(GET CMAKE_CURRENT_SOURCE_DIR STEM DIRNAME)
set(PROJECT_PREFIX anima${DIRNAME})

################################################################################
# Here go the add_subdirectories, no code should be at the root of the project
################################################################################

add_subdirectory(AnatomicalCommands)
add_subdirectory(ArithmeticTools)
add_subdirectory(CommonTools)
add_subdirectory(DiffusionCommands)
add_subdirectory(Resamplers)
add_subdirectory(Symmetry)
