# Must be run with LC_MESSAGES=C set on gpu-visages

set( $ENV{LC_MESSAGES}    "C" )
set(CTEST_SOURCE_DIRECTORY "$ENV{HOME}/anima-src")
set(CTEST_BINARY_DIRECTORY "$ENV{HOME}/anima-build")

set(CTEST_SITE "gpu-visages")
set(CTEST_BUILD_NAME "Linux-GCC")
set(CTEST_PROJECT_NAME "Anima")
set(CTEST_NIGHTLY_START_TIME "00:00:00 EST")

set(CTEST_START_WITH_EMPTY_BINARY_DIRECTORY TRUE)
set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
set(CTEST_BUILD_CONFIGURATION "Release")

set(CTEST_BUILD_OPTIONS " \"-DBUILD_DOCUMENTATION:BOOL=ON\" -DITK_DIR:PATH=$ENV{HOME}/ITK/shared-build \"-DUSE_RPI:BOOL=ON\" \"-DUSE_VTK:BOOL=ON\" \"-DVTK_DIR:PATH=$ENV{HOME}/VTK/shared-build\" \"-DCMAKE_C_FLAGS:STRING=${FLAGS}\" \"-DCMAKE_CXX_FLAGS:STRING=${FLAGS}\" \"-DUSE_RPI:BOOL=ON\" \"-DRPI_DIR:STRING=$ENV{HOME}/asclepiospublic/shared-build\" \"-DUSE_VTK:BOOL=ON\" \"-DNLOPT_INCLUDE_DIR:PATH=$ENV{HOME}/nlopt-install/include\" \"-DNLOPT_LIBRARY_DIR:PATH=$ENV{HOME}/nlopt-install/lib\" \"-DTinyXML2_INCLUDE_DIR:PATH=$ENV{HOME}/tinyxml2/src\" \"-DTinyXML2_LIBRARY_DIR:PATH=$ENV{HOME}/tinyxml2/build\"")

#######################################################################

ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})

find_program(CTEST_GIT_COMMAND NAMES git)

if(NOT EXISTS "${CTEST_SOURCE_DIRECTORY}")
  set(CTEST_CHECKOUT_COMMAND "${CTEST_GIT_COMMAND} clone git@github.com:Inria-Visages/Anima.git ${CTEST_SOURCE_DIRECTORY}")
endif()

set(CTEST_UPDATE_COMMAND "${CTEST_GIT_COMMAND}")

set(CTEST_CONFIGURE_COMMAND "${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE:STRING=${CTEST_BUILD_CONFIGURATION}")
set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} ${CTEST_BUILD_OPTIONS}")
set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} \"-G${CTEST_CMAKE_GENERATOR}\"")
set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} \"${CTEST_SOURCE_DIRECTORY}\"")

ctest_start("Nightly")
ctest_update(RETURN_VALUE count)
ctest_configure()
ctest_build()
ctest_test()
ctest_submit()

message("Number of update changes: ${count}") 
if (${count} GREATER 0)
  execute_process(COMMAND /home/animabot/scriptUpdateDoxygen.sh)
endif()
