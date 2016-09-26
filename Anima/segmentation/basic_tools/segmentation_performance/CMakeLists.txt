if(BUILD_TOOLS)

PROJECT (animaSegPerfAnalyzer)
	
## #############################################################################
## List Sources
## #############################################################################


list_source_files(${PROJECT_NAME}
	${CMAKE_CURRENT_SOURCE_DIR}/src
	${CMAKE_CURRENT_SOURCE_DIR}/src/Analyzer
	${CMAKE_CURRENT_SOURCE_DIR}/src/Results)


## ############################################################################
## Creation of groups of sources for IDE
## ############################################################################

source_group("Header Files" FILES)

source_group("Source Files" FILES
    ${PROJECT_SOURCE_DIR}/src/animaSegmentationPerformance.cpp
    ${PROJECT_SOURCE_DIR}/src/SegPerfApp.cpp
    ${PROJECT_SOURCE_DIR}/src/SegPerfApp.h)

source_group("Source Files\\Analyzer" FILES
    ${PROJECT_SOURCE_DIR}/src/Analyzer/Analyzer.cpp
    ${PROJECT_SOURCE_DIR}/src/Analyzer/Analyzer.h
    ${PROJECT_SOURCE_DIR}/src/Analyzer/SegmentationMeasuresImageFilter.hxx
    ${PROJECT_SOURCE_DIR}/src/Analyzer/SegmentationMeasuresImageFilter.h)

source_group("Source Files\\Results" FILES
    ${PROJECT_SOURCE_DIR}/src/Results/Results.cpp
    ${PROJECT_SOURCE_DIR}/src/Results/Results.h)

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
