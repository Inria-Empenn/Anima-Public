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

set(NLOPT_LIBRARY_DIR "" CACHE PATH "NLOPT library folder")

FIND_LIBRARY(NLOPT_LIBRARY NAMES nlopt libnlopt
  HINTS "${NLOPT_LIBRARY_DIR}" "${NLOPT_LIBRARY_DIR}/Release" "${NLOPT_LIBRARY_DIR}/Debug"
)

mark_as_advanced(NLOPT_LIBRARY)

IF( EXISTS "${NLOPT_INCLUDE_DIR}" AND EXISTS "${NLOPT_LIBRARY}" )
  SET(NLOPT_FOUND TRUE)
else()
  message(FATAL_ERROR "NLOPT not found")
endif()

