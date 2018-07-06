add_subdirectory(average_images)

if (USE_RPI AND RPI_FOUND)
  add_subdirectory(dense_transform_arithmetic)
endif()

add_subdirectory(image_arithmetic)
add_subdirectory(linear_transform_arithmetic)
