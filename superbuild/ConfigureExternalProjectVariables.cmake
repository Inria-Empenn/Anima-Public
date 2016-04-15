if(${USE_GITHUB_SSH})
  set(GITHUB_PREFIX git@github.com:)
else()
  set(GITHUB_PREFIX https://github.com/)
endif()

# Handle Mac compiling option if CMake version greater or equal to 3
set(MACOSX_RPATH_OPTION "-DCMAKE_MACOSX_RPATH:BOOL=OFF")

set(common_c_flags 
  "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_INIT} ${ADDITIONAL_C_FLAGS}"
  )

set(common_cxx_flags 
  "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_INIT} ${ADDITIONAL_CXX_FLAGS}"
  )

set(c_flags "${common_c_flags}")
set(cxx_flags "${common_cxx_flags}")
if (UNIX AND NOT BUILD_SHARED_LIBS AND "${CMAKE_SYSTEM_PROCESSOR}" MATCHES 64)
  set (c_flags "${c_flags} -fPIC")
  set (cxx_flags "${cxx_flags} -fPIC")
endif()

set(common_cache_args
  -DCMAKE_C_FLAGS:STRING=${c_flags}
  -DCMAKE_CXX_FLAGS:STRING=${cxx_flags}
  -DCMAKE_CXX_STANDARD:STRING=${CMAKE_CXX_STANDARD}
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
  -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
  -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
  -DCMAKE_MAKE_PROGRAM:FILEPATH=${CMAKE_MAKE_PROGRAM}
  ${MACOSX_RPATH_OPTION}
)

if (CMAKE_EXTRA_GENERATOR)
  set(cmake_gen "${CMAKE_EXTRA_GENERATOR} -G ${CMAKE_GENERATOR}")
else()
  set(cmake_gen "${CMAKE_GENERATOR}")
endif()
