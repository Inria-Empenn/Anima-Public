set (proj Anima-Private)

set (cmake_args
  ${common_cache_args}
  -DANIMA_DIR:PATH=${Anima_BUILD_DIR}
  -DEXECUTABLE_OUTPUT_PATH=${CMAKE_BINARY_DIR}/bin
  -DLIBRARY_OUTPUT_PATH=${CMAKE_BINARY_DIR}/lib
  -DBUILD_TOOLS:BOOL=${BUILD_ANIMA_TOOLS}
  -DBUILD_TESTING:BOOL=${BUILD_ANIMA_TESTING}
  )

set (location "")
if (NOT DEFINED ${proj}_SRC_DIR)
  set(location GIT_REPOSITORY ${GITHUB_PREFIX}ocommowi/Anima.git GIT_TAG anima-private)
endif()

ExternalProject_Add(${proj}
  ${location}
  DEPENDS Anima
  PREFIX ${CMAKE_BINARY_DIR}/${proj}
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/${proj}
  CMAKE_GENERATOR ${cmake_gen}
  CMAKE_ARGS ${cmake_args}
  BUILD_IN_SOURCE 0
  BINARY_DIR ${CMAKE_BINARY_DIR}/${proj}
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  )

# Update custom target
set (GIT_COMMAND ${GIT_BIN} pull --ff-only)
add_custom_target(update-${proj} 
  COMMAND ${GIT_COMMAND}
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/${proj}
  COMMENT "Updating '${proj}' with '${GIT_COMMAND}'"
  )

set(Update_Repositories "${Update_Repositories};update-${proj}")

# Build custom target
add_custom_target(build-${proj} 
  COMMAND ${CMAKE_GENERATOR} --build . --config ${CMAKE_BUILD_TYPE}
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/${proj}
  COMMENT "build '${proj}' with '${CMAKE_GENERATOR} --build . --config ${CMAKE_BUILD_TYPE}'"
  )

set(Build_Targets "${Build_Targets};build-${proj}")
