# Try to find the nlopt library
# Once done this will define
#
# NLOPT_FOUND         - system has nlopt and it can be used
# NLOPT_INCLUDE_DIR   - directory where the header file can be found
# NLOPT_LIBRARY       - Path where nlopt library file can be found
#

SET (NLOPT_FOUND FALSE)

FIND_PATH(NLOPT_INCLUDE_DIR nlopt.h
  /usr/include
  /usr/local/include
)

FIND_PATH(NLOPT_AUTOGEN_INCLUDE_DIR nlopt.hpp
  /usr/include
  /usr/local/include
)

set(NLOPT_LIBRARY_DIR "" CACHE PATH "NLOPT library folder")

FIND_LIBRARY(NLOPT_LIBRARY NAMES nlopt_cxx libnlopt_cxx
  HINTS "${NLOPT_LIBRARY_DIR}" "${NLOPT_LIBRARY_DIR}/Release" "${NLOPT_LIBRARY_DIR}/Debug"
)

mark_as_advanced(NLOPT_LIBRARY)

IF( EXISTS "${NLOPT_INCLUDE_DIR}" AND EXISTS "${NLOPT_AUTOGEN_INCLUDE_DIR}" AND EXISTS "${NLOPT_LIBRARY}" )
  set(NLOPT_FOUND TRUE)
  set(NLOPT_INCLUDE_DIRS "${NLOPT_INCLUDE_DIR}" "${NLOPT_AUTOGEN_INCLUDE_DIR}" CACHE PATH "NLOPT include directories")
  mark_as_advanced(NLOPT_INCLUDE_DIRS)
else()
  message(FATAL_ERROR "NLOPT not found")
endif()

