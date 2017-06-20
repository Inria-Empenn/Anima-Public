set (proj NLOPT)

set (cmake_args
  ${common_cache_args}
  -DBUILD_GUILE:BOOL=OFF
  -DBUILD_MATLAB:BOOL=OFF
  -DBUILD_OCTAVE:BOOL=OFF
  -DBUILD_PYTHON:BOOL=OFF
  -DUSE_SWIG:BOOL=OFF
  -DWITH_CXX:BOOL=ON
  )

set (location "")
if (NOT DEFINED ${proj}_SRC_DIR)
  set(location GIT_REPOSITORY ${GITHUB_PREFIX}stevengj/nlopt.git GIT_TAG origin/master)
endif()

ExternalProject_Add(${proj}
  ${location}
  PREFIX ${CMAKE_BINARY_DIR}/External-Projects/${proj}
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/External-Projects/${proj}
  CMAKE_GENERATOR ${cmake_gen}
  CMAKE_ARGS ${cmake_args}
  BUILD_ALWAYS 1
  BINARY_DIR ${CMAKE_BINARY_DIR}/${proj}
  INSTALL_COMMAND ""
  )

ExternalProject_Get_Property(${proj} binary_dir)
set(${proj}_BUILD_DIR ${binary_dir})
set(${proj}_SRC_DIR ${CMAKE_SOURCE_DIR}/External-Projects/${proj}/api)
set(${proj}_AUTOGEN_DIR ${binary_dir}/api)

set(Anima_DEPS "${Anima_DEPS};${proj}")
