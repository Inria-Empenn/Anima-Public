set (proj Boost)

set (location "")
if (NOT DEFINED ${proj}_SRC_DIR)
  set(location URL http://freefr.dl.sourceforge.net/project/boost/boost/1.58.0/boost_1_58_0.tar.gz)
endif()

ExternalProject_Add(${proj}
  ${location}
  PREFIX ${CMAKE_BINARY_DIR}/${proj}
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/Projects/${proj}
  BUILD_IN_SOURCE 0
  BINARY_DIR ${CMAKE_BINARY_DIR}/${proj}
  CMAKE_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND ""
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  )

set(${proj}_SRC_DIR ${CMAKE_SOURCE_DIR}/Projects/${proj})
