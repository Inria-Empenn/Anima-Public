
## MANDATORY VARIABLES
# SRC_PATERNS
# PATH_DESTINATION

## OPTIONALS VARIABLES
# COPY_ON_ERROR
# SYMBOLIC

######################################################################
## 0. Check if SRC_PATERNS and PATH_DESTINATION are defined
if(NOT DEFINED SRC_PATERNS)
  message(FATAL_ERROR "SRC_PATERNS is not defined")
endif()

if(NOT DEFINED PATH_DESTINATION)
  message(FATAL_ERROR "PATH_DESTINATION is not defined")
endif()

message("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
message("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
message("XXXXXXXX    ${SRC_PATERNS}   XXXXXXXXXXXXXXX")
message("XXXXXXXX    ${PATH_DESTINATION}   XXXXXXXXXXXXXXX")
message("XXXXXXXX    ${COPY_ON_ERROR}   XXXXXXXXXXXXXXX")
message("XXXXXXXX    ${SYMBOLIC}   XXXXXXXXXXXXXXX")
message("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
message("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")

######################################################################
## 1. Check if COPY_ON_ERROR and SYMBOLIC are defined
if(NOT DEFINED COPY_ON_ERROR)
  set(COPY_ON_ERROR "COPY_ON_ERROR")
else()
  if(COPY_ON_ERROR)
    set(COPY_ON_ERROR "COPY_ON_ERROR")
  else()
    set(COPY_ON_ERROR "")
  endif()
endif()

if(NOT DEFINED SYMBOLIC)
  set(SYMBOLIC "")
else()
  if(SYMBOLIC)
    set(SYMBOLIC "SYMBOLIC")
  else()
    set(SYMBOLIC "")
  endif()
endif()


######################################################################
## 2. Check if PATH_DESTINATION is a directory or a file
if(EXISTS "${PATH_DESTINATION}")
  if(IS_DIRECTORY ${PATH_DESTINATION})
	set(DEST_PATH_BASE ${PATH_DESTINATION})
	set(DEST_FILE_NAME "")
  else()
	file(REMOVE ${PATH_DESTINATION})
    cmake_path(GET PATH_DESTINATION PARENT_PATH DEST_PATH_BASE)
    cmake_path(GET PATH_DESTINATION FILENAME DEST_FILE_NAME)
  endif()
else()
  cmake_path(GET PATH_DESTINATION PARENT_PATH DEST_PATH_BASE)
  cmake_path(GET PATH_DESTINATION FILENAME    DEST_FILE_NAME)
endif()


######################################################################
## 3. Create links
foreach(MATCH_PATERN ${SRC_PATERNS})
  file(GLOB MATCH_FILES LIST_DIRECTORIES false ${MATCH_PATERN})

  foreach(FILE_TO_LINK ${MATCH_FILES})  
    cmake_path(GET FILE_TO_LINK FILENAME FILE_NAME)
    cmake_path(GET FILE_TO_LINK PARENT_PATH PATH_FILE)
	
    if("${DEST_FILE_NAME}" STREQUAL "")
        set(DEST_PATH "${DEST_PATH_BASE}/${FILE_NAME}")
    else()
        set(DEST_PATH "${DEST_PATH_BASE}/${DEST_FILE_NAME}")
    endif()

    if(EXISTS ${DEST_PATH})
      file(REMOVE "${DEST_PATH}")
	endif()

	file(CREATE_LINK ${FILE_TO_LINK} ${DEST_PATH} RESULT RES ${COPY_ON_ERROR} ${SYMBOLIC})
  endforeach()

endforeach()

