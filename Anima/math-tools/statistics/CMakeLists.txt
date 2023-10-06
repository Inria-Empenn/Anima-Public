## #############################################################################
## Subdirs exe directories
## #############################################################################

add_subdirectory(build_samples)
add_subdirectory(gamma_estimation)
add_subdirectory(dirichlet_estimation)

if(BUILD_TESTING)
  add_subdirectory(gamma_estimation_test)
  add_subdirectory(dirichlet_estimation_test)
endif()

add_subdirectory(boot_strap_4d_volume)
add_subdirectory(local_patch_covariance_distance)
add_subdirectory(local_patch_mean_distance)
add_subdirectory(low_memory_tools)
add_subdirectory(roi_intensities_stats)
