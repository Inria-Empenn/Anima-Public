# Try to find the TinyXML2 library
# Once done this will define
#
# TinyXML2_FOUND         - system has TinyXML2 and it can be used
# TinyXML2_INCLUDE_DIR   - directory where the header file can be found
# TinyXML2_LIBRARY       - Path where TinyXML2 library file can be found
#

SET (TinyXML2_FOUND FALSE)

FIND_PATH(TinyXML2_INCLUDE_DIR tinyxml2.h
  /usr/include
  /usr/local/include
)

set(TinyXML2_LIBRARY_DIR "" CACHE PATH "TinyXML2 library folder")

FIND_LIBRARY(TinyXML2_LIBRARY NAMES tinyxml2 libtinyxml2
  HINTS "${TinyXML2_LIBRARY_DIR}" "${TinyXML2_LIBRARY_DIR}/Release" "${TinyXML2_LIBRARY_DIR}/Debug"
)

mark_as_advanced(TinyXML2_LIBRARY)

IF( EXISTS "${TinyXML2_INCLUDE_DIR}" AND EXISTS "${TinyXML2_LIBRARY}" )
  SET(TinyXML2_FOUND TRUE)
else()
  message(FATAL_ERROR "TinyXML2 not found")
endif()
