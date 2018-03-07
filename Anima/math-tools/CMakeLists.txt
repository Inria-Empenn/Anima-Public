project(ANIMA-MATHS)

################################################################################
# Here go the add_subdirectories, no code should be at the root of the project
################################################################################

if (USE_GMP AND GMP_FOUND AND USE_MPFR AND MPFR_FOUND AND USE_FLINT AND FLINT_FOUND AND USE_ARB AND ARB_FOUND)
  add_subdirectory(arb_special_functions)
endif()

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

if (BUILD_TESTING)
  add_subdirectory(matrix_operations/qr_test)
endif()