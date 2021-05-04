project(ANIMA-DIFFUSION)

################################################################################
# Here go the add_subdirectories, no code should be at the root of the project
################################################################################

add_subdirectory(dti_estimator)
add_subdirectory(dti_tools)
add_subdirectory(mcm)
add_subdirectory(mcm_estimator)
add_subdirectory(mcm_tools)

if(USE_VTK AND VTK_FOUND)
  add_subdirectory(fibers_analysis)
  add_subdirectory(tractography)
endif()

add_subdirectory(odf)
