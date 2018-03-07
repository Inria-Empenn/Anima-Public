# Try to find the arb library
# Once done this will define
#
# ARB_FOUND         - system has arb and it can be used
# ARB_INCLUDE_DIR   - directory where the header file can be found
# ARB_LIBRARY       - Path where arb library file can be found
#

SET (ARB_FOUND FALSE)

FIND_PATH(ARB_INCLUDE_DIR arb.h
  /usr/include
  /usr/local/include
)

set(ARB_LIBRARY_DIR "" CACHE PATH "ARB library folder")

FIND_LIBRARY(ARB_LIBRARY NAMES arb libarb
  HINTS "${ARB_LIBRARY_DIR}"
)

mark_as_advanced(ARB_LIBRARY)

IF( EXISTS "${ARB_INCLUDE_DIR}" AND EXISTS "${ARB_LIBRARY}" )
  set(ARB_FOUND TRUE)
else()
  message(FATAL_ERROR "ARB not found")
endif()

