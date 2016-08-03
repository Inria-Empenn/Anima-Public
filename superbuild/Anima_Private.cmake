set (proj Anima_Private)

set (cmake_args
  ${common_cache_args}
  -DANIMA_DIR:PATH=${Anima_BUILD_DIR}
  -DEXECUTABLE_OUTPUT_PATH=${CMAKE_BINARY_DIR}/bin
  -DLIBRARY_OUTPUT_PATH=${CMAKE_BINARY_DIR}/lib
  -DBUILD_TOOLS:BOOL=${BUILD_ANIMA_TOOLS}
  -DBUILD_TESTING:BOOL=${BUILD_ANIMA_TESTING}
  -DBUILD_DOCUMENTATION:BOOL=${BUILD_ANIMA_DOCUMENTATION}
  )

set (location "")
if (NOT DEFINED ${proj}_SRC_DIR)
  set(location GIT_REPOSITORY ${GITHUB_PREFIX}Inria-Visages/Anima.git GIT_TAG origin/yogesh_normalization)
endif()

ExternalProject_Add(${proj}
  ${location}
  DEPENDS Anima
  PREFIX ${CMAKE_BINARY_DIR}/${proj}
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/${proj}
  CMAKE_GENERATOR ${cmake_gen}
  CMAKE_ARGS ${cmake_args}
  BUILD_ALWAYS 1
  BINARY_DIR ${CMAKE_BINARY_DIR}/${proj}
  INSTALL_COMMAND ""
  )
