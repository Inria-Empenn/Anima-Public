project(libAnimaSmoothing)

## #############################################################################
## List Sources
## #############################################################################

list_source_files_recurse(${PROJECT_NAME} ${CMAKE_CURRENT_SOURCE_DIR})
list_include_dirs(${PROJECT_NAME})

## #############################################################################
## add lib
## #############################################################################

add_library(${PROJECT_NAME} INTERFACE ${${PROJECT_NAME}_CFILES})

## ##############################################################################
## set target properties
## ##############################################################################

target_include_directories(${PROJECT_NAME} INTERFACE ${${PROJECT_NAME}_INCLUDE_PATHS})
