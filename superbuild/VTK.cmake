set (proj VTK)

set (cmake_args
  ${common_cache_args}
  -DBUILD_EXAMPLES:BOOL=OFF
  -DBUILD_TESTING:BOOL=OFF
  )

set (location "")
if (NOT DEFINED ${proj}_SRC_DIR)
  set(tag v6.3.0)
  set(location GIT_REPOSITORY ${GITHUB_PREFIX}Kitware/VTK.git GIT_TAG ${tag})
endif()

ExternalProject_Add(${proj}
  ${location}
  PREFIX ${CMAKE_BINARY_DIR}/External-Projects/${proj}
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/External-Projects/${proj}
  CMAKE_GENERATOR ${cmake_gen}
  CMAKE_ARGS ${cmake_args}
  BUILD_IN_SOURCE 0
  BINARY_DIR ${CMAKE_BINARY_DIR}/${proj}
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  )

ExternalProject_Get_Property(${proj} binary_dir)
set(${proj}_BUILD_DIR ${binary_dir})
set(${proj}_SRC_DIR ${CMAKE_SOURCE_DIR}/External-Projects/${proj})

set(Anima_DEPS "${Anima_DEPS};${proj}")

# Update custom target
set (GIT_COMMAND ${GIT_BIN} pull --ff-only)
add_custom_target(update-${proj} 
  COMMAND ${GIT_COMMAND}
  WORKING_DIRECTORY ${${proj}_SRC_DIR}
  COMMENT "Updating '${proj}' with '${GIT_COMMAND}'"
  )

set(Update_Repositories "${Update_Repositories};update-${proj}")

# Build custom target
add_custom_target(build-${proj} 
  COMMAND ${CMAKE_COMMAND} --build . --config ${CMAKE_BUILD_TYPE}
  WORKING_DIRECTORY ${${proj}_BUILD_DIR}
  COMMENT "build '${proj}' with '${CMAKE_COMMAND} --build . --config ${CMAKE_BUILD_TYPE}'"
  )
