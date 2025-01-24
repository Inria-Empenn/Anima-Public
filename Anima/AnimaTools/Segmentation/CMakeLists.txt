cmake_path(GET CMAKE_CURRENT_SOURCE_DIR STEM DIRNAME  )
set(PROJECT_PREFIX anima${DIRNAME})

################################################################################
# Here go the add_subdirectories, no code should be at the root
################################################################################

add_subdirectory(BasicTools)
add_subdirectory(Commands)
add_subdirectory(GraphCut)
add_subdirectory(Isosurface)
add_subdirectory(LesionSimulation)
add_subdirectory(Staple)
add_subdirectory(TissuesEMClassification)
add_subdirectory(ValidationTools)
