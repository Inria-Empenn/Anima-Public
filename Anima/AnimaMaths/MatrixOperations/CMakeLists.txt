project(AnimaMathsMatrixOperations)

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


# ##############################################################################
# set target properties
# ##############################################################################

#set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "${${PROJECT_NAME}_HEADERS} ${${PROJECT_NAME}_TEMPLATES}")
target_include_directories(${PROJECT_NAME}  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})


## #############################################################################
## Link
## #############################################################################

target_link_libraries(${PROJECT_NAME} INTERFACE
  AnimaMathsSpecialFunctions
  AnimaRegistrationTransformations
  )
  