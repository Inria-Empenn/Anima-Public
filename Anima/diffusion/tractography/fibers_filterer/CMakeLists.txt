if(BUILD_TOOLS AND USE_VTK AND VTK_FOUND)

project(animaFibersFilterer)

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
  AnimaTractography
  ${ITKIO_LIBRARIES}
  vtkCommonCore
  vtkFiltersCore
  vtkIOXML
  vtkIOLegacy
  vtksys
  )

## #############################################################################
## install
## #############################################################################

set_exe_install_rules(${PROJECT_NAME})

endif()
