project(AnimaDiffDTIBackgroundNoiseVarianceEstimatorFilter)

## #############################################################################
## List Sources
## #############################################################################

list_source_files(${PROJECT_NAME}
  ${CMAKE_CURRENT_SOURCE_DIR}
)

## #############################################################################
## add executable
## #############################################################################

add_library(${PROJECT_NAME} INTERFACE
  ${${PROJECT_NAME}_CFILES}
)

## #############################################################################
## Link
## #############################################################################

target_link_libraries(${PROJECT_NAME} INTERFACE
  ${ITK_TRANSFORM_LIBRARIES}
)

## ##############################################################################
## set target properties
## ##############################################################################

target_include_directories(${PROJECT_NAME} INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}
)
