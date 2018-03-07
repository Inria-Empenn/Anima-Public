project(ANIMA-FILTERING)

################################################################################
# Here go the add_subdirectories, no code should be at the root of the project
################################################################################

add_subdirectory(bias_correction)
add_subdirectory(denoising)
add_subdirectory(dti_tools)

if (USE_GMP AND GMP_FOUND AND USE_MPFR AND MPFR_FOUND AND USE_FLINT AND FLINT_FOUND AND USE_ARB AND ARB_FOUND) 
  add_subdirectory(gaussianizer)
endif()

add_subdirectory(noise_generator)
add_subdirectory(nyul_standardization)
add_subdirectory(regularization)
