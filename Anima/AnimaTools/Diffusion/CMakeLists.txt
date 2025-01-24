################################################################################
# Here go the add_subdirectories, no code should be at the root of the project
################################################################################

cmake_path(GET CMAKE_CURRENT_SOURCE_DIR STEM DIRNAME  )
set(PROJECT_PREFIX anima${DIRNAME})

add_subdirectory(dti)
add_subdirectory(MCM)

if(USE_VTK AND VTK_FOUND)
  add_subdirectory(FibersAnalysis)
endif()

add_subdirectory(odf)
