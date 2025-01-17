if(BUILD_TOOLS)

project(animaRelaxoGMMT2Estimation)

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
  ITKOptimizers
  ${NLOPT_LIBRARY}
  AnimaCore
  AnimaMathsOptimizers
  AnimaRelaxoEstimation
)

## #############################################################################
## install
## #############################################################################

set_exe_install_rules(${PROJECT_NAME})

endif()
