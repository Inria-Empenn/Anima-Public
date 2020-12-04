set (proj ITK)

set(VTK_DEP_ARGS "")
set(VTK_PROJ_DEP "")

if (USE_VTK)
  set(VTK_DEP_ARGS
    -DModule_ITKVtkGlue:BOOL=ON
    -DVTK_DIR:PATH=${VTK_BUILD_DIR}
  )

  set(VTK_PROJ_DEP
    VTK
  )
endif()

set(ITK_OPTIM_FLAGS "")

set (cmake_args
  ${common_cache_args}
  -DBUILD_EXAMPLES:BOOL=OFF
  -DBUILD_TESTING:BOOL=OFF
  -DModule_ITKReview:BOOL=ON
  -DITK_C_OPTIMIZATION_FLAGS:STRING=${ITK_OPTIM_FLAGS}
  -DITK_CXX_OPTIMIZATION_FLAGS:STRING=${ITK_OPTIM_FLAGS}
  ${VTK_DEP_ARGS}
  )

set (location "")
if (NOT DEFINED ${proj}_SRC_DIR)
  set(tag origin/release)
  set(location GIT_REPOSITORY ${GITHUB_PREFIX}InsightSoftwareConsortium/ITK.git GIT_TAG ${tag})
endif()

ExternalProject_Add(${proj}
  ${location}
  DEPENDS ${VTK_PROJ_DEP}
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
