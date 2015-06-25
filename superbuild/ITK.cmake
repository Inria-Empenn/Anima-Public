set (proj ITK)

set (cmake_args
  -DCMAKE_C_FLAGS:STRING=${cmake_flags}
  -DCMAKE_CXX_FLAGS:STRING=${cmake_flags}
  -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
  -DBUILD_EXAMPLES:BOOL=OFF
  -DBUILD_TESTING:BOOL=OFF
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  )

set (location "")
if (NOT DEFINED ${proj}_SRC_DIR)
  set(location GIT_REPOSITORY ${GITHUB_PREFIX}InsightSoftwareConsortium/ITK.git)
endif()

ExternalProject_Add(${proj}
  ${location}
  PREFIX ${CMAKE_BINARY_DIR}/${proj}
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/Projects/${proj}
  CMAKE_ARGS ${cmake_args}
  BUILD_IN_SOURCE 0
  BINARY_DIR ${CMAKE_BINARY_DIR}/${proj}
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  )

ExternalProject_Get_Property(${proj} binary_dir)
set(${proj}_BUILD_DIR ${binary_dir})
set(${proj}_SRC_DIR ${CMAKE_SOURCE_DIR}/Projects/${proj})
