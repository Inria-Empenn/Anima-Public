set (proj RPI)

set (cmake_args
  -DCMAKE_C_FLAGS:STRING=${cmake_flags}
  -DCMAKE_CXX_FLAGS:STRING=${cmake_flags}
  ${MACOSX_RPATH_OPTION}
  -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
  -DRPI_BUILD_EXAMPLES:BOOL=OFF
  -DITK_DIR:PATH=${ITK_BUILD_DIR}
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  )

set (location "")
if (NOT DEFINED ${proj}_SRC_DIR)
  set(location GIT_REPOSITORY ${GITHUB_PREFIX}Inria-Asclepios/RPI.git)
endif()

set(${proj}_DEPS "")
if (NOT USE_SYSTEM_ITK)
  set(${proj}_DEPS "${${proj}_DEPS};ITK")
endif()

ExternalProject_Add(${proj}
  ${location}
  DEPENDS ${${proj}_DEPS}
  PREFIX ${CMAKE_BINARY_DIR}/External-Projects/${proj}
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/External-Projects/${proj}
  CMAKE_ARGS ${cmake_args}
  BUILD_IN_SOURCE 0
  BINARY_DIR ${CMAKE_BINARY_DIR}/${proj}
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  )

ExternalProject_Get_Property(${proj} binary_dir)
set(${proj}_BUILD_DIR ${binary_dir})
set(${proj}_SRC_DIR ${CMAKE_SOURCE_DIR}/External-Projects/${proj})
