# Try to find the gmp library
# Once done this will define
#
# GMP_FOUND         - system has gmp and it can be used
# GMP_INCLUDE_DIR   - directory where the header file can be found
# GMP_LIBRARY       - Path where gmp library file can be found
#

SET (GMP_FOUND FALSE)

FIND_PATH(GMP_INCLUDE_DIR gmp.h
  /usr/include
  /usr/local/include
)

set(GMP_LIBRARY_DIR "" CACHE PATH "GMP library folder")

FIND_LIBRARY(GMP_LIBRARY NAMES gmp libgmp
  HINTS "${GMP_LIBRARY_DIR}"
)

mark_as_advanced(GMP_LIBRARY)

IF( EXISTS "${GMP_INCLUDE_DIR}" AND EXISTS "${GMP_LIBRARY}" )
  set(GMP_FOUND TRUE)
else()
  message(FATAL_ERROR "GMP not found")
endif()

