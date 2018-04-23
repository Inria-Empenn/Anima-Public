project(ANIMA-MATHS)

################################################################################
# Here go the add_subdirectories, no code should be at the root of the project
################################################################################

add_subdirectory(arithmetic)
add_subdirectory(common_tools)

if (USE_VTK AND VTK_FOUND)
  add_subdirectory(data_io)
endif()

add_subdirectory(multi_compartment_base)
add_subdirectory(optimizers)
add_subdirectory(special_functions)
add_subdirectory(sphere_operations)
add_subdirectory(spherical_harmonics)
add_subdirectory(statistics)
add_subdirectory(statistical_tests)
