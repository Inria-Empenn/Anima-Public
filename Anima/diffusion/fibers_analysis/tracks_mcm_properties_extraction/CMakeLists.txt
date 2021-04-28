if(BUILD_TOOLS)

project(animaTracksMCMPropertiesExtraction)


## #############################################################################
## List Sources
## #############################################################################

list_source_files(${PROJECT_NAME}
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

## #############################################################################
## add executable
## #############################################################################

add_executable(${PROJECT_NAME}
  ${${PROJECT_NAME}_CFILES}
  )

## #############################################################################
## Link
## #############################################################################

target_link_libraries(${PROJECT_NAME}
  ${ITKIO_LIBRARIES}
  AnimaMCMPrivate
  AnimaMCMPrivateBase
  AnimaMCM
  AnimaMCMBase
  AnimaDataIO
  AnimaTractography
  AnimaOptimizers
  ${TinyXML2_LIBRARY}
  ${VTK_PREFIX}CommonCore
  ${VTK_PREFIX}IOXML
  ${VTK_PREFIX}IOLegacy
  )

## #############################################################################
## install
## #############################################################################

set_exe_install_rules(${PROJECT_NAME})

endif()
