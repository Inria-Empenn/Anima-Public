# Should not be here if VTK not found but you're never too careful
if(USE_VTK AND VTK_FOUND)

project(AnimaDiffLogTensorImageFilter)

## #############################################################################
## List Sources
## #############################################################################

list_source_files(${PROJECT_NAME}
  ${CMAKE_CURRENT_SOURCE_DIR}
)

## #############################################################################
## add lib
## #############################################################################

add_library(${PROJECT_NAME} INTERFACE
  ${${PROJECT_NAME}_CFILES}
)

## ##############################################################################
## set target properties
## ##############################################################################

target_include_directories(${PROJECT_NAME} INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}
)

## #############################################################################
## Link
## #############################################################################

target_link_libraries(${PROJECT_NAME} INTERFACE
  ${ITKIO_LIBRARIES}
  AnimaMathsMatrixOperations
)

endif() #USE_VTK AND VTK_FOUND
