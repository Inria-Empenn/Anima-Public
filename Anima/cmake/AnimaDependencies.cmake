# ITK
find_package(ITK REQUIRED)

set(_ITKVersionString "${ITK_VERSION_MAJOR}.${ITK_VERSION_MINOR}.${ITK_VERSION_PATCH}" )
math(EXPR _ITKVersionNum "${ITK_VERSION_MAJOR}*100*100 + ${ITK_VERSION_MINOR}*100 + ${ITK_VERSION_PATCH}")
  
if(_ITKVersionNum LESS 41300)
  message(SEND_ERROR "The ITK version you want to use (${_ITKVersionString}) is not supported by this project. Please use a more recent version of ITK. The minimum required version is 4.13.0")
else()
  include(${ITK_USE_FILE})
endif()

# Boost
find_package(Boost 1.40.0 REQUIRED)

# TCLAP
option(BUILD_TOOLS "Build command line executables" ON)
mark_as_advanced(BUILD_TOOLS)
if(BUILD_TOOLS)
  find_package(TCLAP REQUIRED)
endif()

# NL-OPT
option(USE_NLOPT 
  "Use NL-OPT external library (necessary for some optimizers)" 
  OFF
  )

if(USE_NLOPT)
  find_package(NLOPT REQUIRED)
endif()

# TinyXML2
if (BUILD_MODULE_REGISTRATION OR BUILD_MODULE_DIFFUSION)
	find_package(TinyXML2 REQUIRED)
endif()

# VTK
if (BUILD_MODULE_REGISTRATION OR BUILD_MODULE_DIFFUSION OR BUILD_MODULE_MATHS)
  option(USE_VTK "Use VTK libraries (necessary for some registration tools and tractography)" OFF)

  if(USE_VTK)
    find_package(VTK REQUIRED)

    if (VTK_VERSION_MAJOR LESS 6)
      message(SEND_ERROR "VTK has to be version 6 or higher.")
    endif()

    include(${VTK_USE_FILE})
  endif()
endif()

# Asclepios RPI
if (BUILD_MODULE_REGISTRATION)
  option(USE_RPI "Use Asclepios RPI library (necessary for Block-Matching)" OFF)

  if(USE_RPI)
    find_package(RPI REQUIRED)
    if(RPI_FOUND)
      include(${RPI_USE_FILE})
    endif(RPI_FOUND)
  endif()
endif()
