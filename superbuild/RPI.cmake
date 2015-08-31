set (proj RPI)

set (cmake_args
  ${common_cache_args}
  -DRPI_BUILD_EXAMPLES:BOOL=OFF
  -DITK_DIR:PATH=${ITK_BUILD_DIR}
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
foreach (dep ${${proj}_DEPS})
    set(build-${proj}_deps build-${dep} ${build-${proj}_deps})
endforeach()

add_custom_target(build-${proj} 
  COMMAND ${CMAKE_COMMAND} --build . --config ${CMAKE_BUILD_TYPE}
  WORKING_DIRECTORY ${${proj}_BUILD_DIR}
  COMMENT "build '${proj}' with '${CMAKE_COMMAND} --build . --config ${CMAKE_BUILD_TYPE}'"
  DEPENDS ${build-${proj}_deps}
  )
