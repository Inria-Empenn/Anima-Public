# Try to find the mpfr library
# Once done this will define
#
# MPFR_FOUND         - system has mpfr and it can be used
# MPFR_INCLUDE_DIR   - directory where the header file can be found
# MPFR_LIBRARY       - Path where mpfr library file can be found
#

SET (MPFR_FOUND FALSE)

FIND_PATH(MPFR_INCLUDE_DIR mpfr.h
  /usr/include
  /usr/local/include
)

set(MPFR_LIBRARY_DIR "" CACHE PATH "MPFR library folder")

FIND_LIBRARY(MPFR_LIBRARY NAMES mpfr libmpfr
  HINTS "${MPFR_LIBRARY_DIR}"
)

mark_as_advanced(MPFR_LIBRARY)

IF( EXISTS "${MPFR_INCLUDE_DIR}" AND EXISTS "${MPFR_LIBRARY}" )
  set(MPFR_FOUND TRUE)
else()
  message(FATAL_ERROR "MPFR not found")
endif()

