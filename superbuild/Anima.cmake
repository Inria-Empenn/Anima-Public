set (proj Anima)

set (cmake_args
  ${common_cache_args}
  -DEXECUTABLE_OUTPUT_PATH=${CMAKE_BINARY_DIR}/bin
  -DLIBRARY_OUTPUT_PATH=${CMAKE_BINARY_DIR}/lib
  -DBUILD_TOOLS:BOOL=${BUILD_ANIMA_TOOLS}
  -DBUILD_TESTING:BOOL=${BUILD_ANIMA_TESTING}
  -DBUILD_ALL_MODULES:BOOL=ON
  -DBUILD_MODULE_MATHS:BOOL=ON
  -DBUILD_MODULE_FILTERING:BOOL=ON
  -DBUILD_MODULE_REGISTRATION:BOOL=ON
  -DBUILD_MODULE_DIFFUSION:BOOL=ON
  -DBUILD_MODULE_SEGMENTATION:BOOL=ON
  -DBUILD_MODULE_QUANTITATIVE_MRI:BOOL=ON
  -DTinyXML2_INCLUDE_DIR:PATH=${TinyXML2_SRC_DIR}
  -DTinyXML2_LIBRARY_DIR:PATH=${TinyXML2_BUILD_DIR}
  -DUSE_NLOPT=${USE_NLOPT}
  -DNLOPT_INCLUDE_DIR:PATH=${NLOPT_SRC_DIR}
  -DNLOPT_LIBRARY_DIR:PATH=${NLOPT_BUILD_DIR}
  -DUSE_RPI:BOOL=${USE_RPI}
  -DUSE_VTK:BOOL=${USE_VTK}
  -DRPI_DIR:PATH=${RPI_BUILD_DIR}
  -DITK_DIR:PATH=${ITK_BUILD_DIR}
  -DVTK_DIR:PATH=${VTK_BUILD_DIR}
  -DBoost_INCLUDE_DIR:PATH=${Boost_SRC_DIR}
  -DTCLAP_INCLUDE_DIR:PATH=${TCLAP_SRC_DIR}/include
  )

ExternalProject_Add(${proj}
  DEPENDS ${${proj}_DEPS}
  PREFIX ${CMAKE_BINARY_DIR}/${proj}
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/${proj}
  CMAKE_GENERATOR ${cmake_gen}
  CMAKE_ARGS ${cmake_args}
  BUILD_IN_SOURCE 0
  BINARY_DIR ${CMAKE_BINARY_DIR}/${proj}
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  )

ExternalProject_Get_Property(${proj} binary_dir)
set(${proj}_BUILD_DIR ${binary_dir})

# Update custom target
set (GIT_COMMAND ${GIT_BIN} pull --ff-only)
add_custom_target(update-${proj} 
  COMMAND ${GIT_COMMAND}
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  COMMENT "Updating '${proj}' with '${GIT_COMMAND}'"
  )

set(Update_Repositories "${Update_Repositories};update-${proj}")

# Build custom target
add_custom_target(build-${proj} 
  COMMAND ${CMAKE_COMMAND} --build . --config ${CMAKE_BUILD_TYPE}
  WORKING_DIRECTORY ${${proj}_BUILD_DIR}
  COMMENT "build '${proj}' with '${CMAKE_COMMAND} --build . --config ${CMAKE_BUILD_TYPE}'"
  )

set(Build_Targets "${Build_Targets};build-${proj}")
