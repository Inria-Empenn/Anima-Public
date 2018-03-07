# Try to find the flint library
# Once done this will define
#
# FLINT_FOUND         - system has flint and it can be used
# FLINT_INCLUDE_DIR   - directory where the header file can be found
# FLINT_LIBRARY       - Path where flint library file can be found
#

SET (FLINT_FOUND FALSE)

FIND_PATH(FLINT_INCLUDE_DIR flint/flint.h
  /usr/include
  /usr/local/include
)

set(FLINT_LIBRARY_DIR "" CACHE PATH "FLINT library folder")

FIND_LIBRARY(FLINT_LIBRARY NAMES flint libflint
  HINTS "${FLINT_LIBRARY_DIR}"
)

mark_as_advanced(FLINT_LIBRARY)

IF( EXISTS "${FLINT_INCLUDE_DIR}" AND EXISTS "${FLINT_LIBRARY}" )
  set(FLINT_FOUND TRUE)
else()
  message(FATAL_ERROR "FLINT not found")
endif()

