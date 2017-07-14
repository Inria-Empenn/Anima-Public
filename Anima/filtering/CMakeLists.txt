project(ANIMA-FILTERING)

################################################################################
# Here go the add_subdirectories, no code should be at the root of the project
################################################################################

add_subdirectory(denoising)
add_subdirectory(intensity_normalization)
add_subdirectory(noise_generation)
add_subdirectory(regularization)
add_subdirectory(dti_tools)
add_subdirectory(bias_correction)
