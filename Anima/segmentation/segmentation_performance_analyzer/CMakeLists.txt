if(BUILD_TOOLS)

PROJECT (animaSegPerfAnalyzer)

## #############################################################################
## List Sources
## #############################################################################

list_source_files(${PROJECT_NAME}
	${CMAKE_CURRENT_SOURCE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}/Analyzer
	${CMAKE_CURRENT_SOURCE_DIR}/Results
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
  )

endif()
