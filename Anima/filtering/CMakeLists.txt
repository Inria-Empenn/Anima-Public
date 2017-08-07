project(ANIMA-FILTERING)

################################################################################
# Here go the add_subdirectories, no code should be at the root of the project
################################################################################

add_subdirectory(bias_correction)
add_subdirectory(denoising)
add_subdirectory(dti_tools)
add_subdirectory(noise_generator)
add_subdirectory(nyul_standardization)
add_subdirectory(regularization)
