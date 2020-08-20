set (proj VTK)

set (cmake_args
  ${common_cache_args}
  -DBUILD_EXAMPLES:BOOL=OFF
  -DBUILD_TESTING:BOOL=OFF
  -DVTK_LEGACY_REMOVE:BOOL=ON
  )

set (location "")
if (NOT DEFINED ${proj}_SRC_DIR)
  set(tag origin/release)
  set(location GIT_REPOSITORY ${GITHUB_PREFIX}Kitware/VTK.git GIT_TAG ${tag})
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
set(${proj}_SRC_DIR ${CMAKE_SOURCE_DIR}/External-Projects/${proj})

set(Anima_DEPS "${Anima_DEPS};${proj}")
