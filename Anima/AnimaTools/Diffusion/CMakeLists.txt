################################################################################
# Here go the add_subdirectories, no code should be at the root of the project
################################################################################

add_subdirectory(dti)
add_subdirectory(MCMEstimator)
add_subdirectory(MCMEstimator/LowMemory)
add_subdirectory(MCM_ModelAveraging)
add_subdirectory(MCMTools)

if(USE_VTK AND VTK_FOUND)
  add_subdirectory(FibersAnalysis)
endif()

add_subdirectory(odf)
