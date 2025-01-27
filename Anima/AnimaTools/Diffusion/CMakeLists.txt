################################################################################
# Here go the add_subdirectories, no code should be at the root
################################################################################

cmake_path(GET CMAKE_CURRENT_SOURCE_DIR STEM DIRNAME  )
set(PROJECT_PREFIX anima${DIRNAME})

add_subdirectory(dti)
add_subdirectory(MCM)
add_subdirectory(FibersAnalysis)
add_subdirectory(odf)
