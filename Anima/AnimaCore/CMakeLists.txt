project(AnimaCore)

## ##############################################################################
## Here go the add_subdirectories, no code should be at the root of the project
## ##############################################################################

add_subdirectory(BasicImageFilters)

if(USE_VTK AND VTK_FOUND)
  add_subdirectory(DataIO)
endif()


