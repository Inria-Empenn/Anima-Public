

###############################################################################
## Set Options and CACHE variables for user's choices
if(WIN32)
    set(_LINK_TYPE OFF)
else()
    set(_LINK_TYPE ON )
endif (WIN32)

set(OLD_NAME_COMPATIBILITY        ON                                  CACHE BOOL "If this is activate compilation will generate also binaries with old names")
set(OLD_NAME_LINK_TYPE            ${_LINK_TYPE}                       CACHE BOOL "If OLD_NAME_LINK_TYPE is activate old names will be Symlink, else hardlink. By default HardLink on Windows and SymLink on Unix like")
set(OLD_NAME_IN_SPECIFIC_FOLDER   ON                                  CACHE BOOL "If this is activate old names will be generated into specific folder")
set(OLD_NAME_SPECIFIC_FOLDER_PATH "${CMAKE_BINARY_DIR}/compatibility" CACHE PATH "Path where compatibility binaries will be stored.")

if(${OLD_NAME_IN_SPECIFIC_FOLDER})
    set(OLD_NAME_SPECIFIC_FOLDER_PATH "${OLD_NAME_SPECIFIC_FOLDER_PATH}" CACHE PATH     "Path where compatibility binaries will be stored." FORCE)
    set(OLD_NAME_PATH "${OLD_NAME_SPECIFIC_FOLDER_PATH}")
else()
    set(OLD_NAME_SPECIFIC_FOLDER_PATH "${OLD_NAME_SPECIFIC_FOLDER_PATH}" CACHE INTERNAL "Path where compatibility binaries will be stored." FORCE)
    set(OLD_NAME_PATH "${CMAKE_BINARY_DIR}/bin")
endif()




## #############################################################################
## Generate OLD NAME LINK COMMAND for an target. To use with add_executable.
## #############################################################################
function(GenerateOldExeNameCommand target oldName oldPath OutVar)
    #####################################################
    ##  determine sourceFile var value
    get_property(GENERATOR_MULTI_CONFIG GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)

    if(${GENERATOR_MULTI_CONFIG})
        get_target_property(_outputDir_debug   ${target} RUNTIME_OUTPUT_DIRECTORY_DEBUG          )
        get_target_property(_outputDir_release ${target} RUNTIME_OUTPUT_DIRECTORY_RELEASE        )
        get_target_property(_outputDir_minsize ${target} RUNTIME_OUTPUT_DIRECTORY_MINSIZEREL     )
        get_target_property(_outputDir_reldbg  ${target} RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO )
        
        set(sourceFile  $<$<CONFIG:debug>:${_outputDir_debug}/${target}${CMAKE_EXECUTABLE_SUFFIX}>
                        $<$<CONFIG:release>:${_outputDir_release}/${target}${CMAKE_EXECUTABLE_SUFFIX}>
                        $<$<CONFIG:MinSizeRel>:${_outputDir_minsize}/${target}${CMAKE_EXECUTABLE_SUFFIX}>
                        $<$<CONFIG:RelWithDebInfo>:${_outputDir_reldbg}/${target}${CMAKE_EXECUTABLE_SUFFIX}>)
    else()                                                                                 
        get_target_property(sourceFile ${target} RUNTIME_OUTPUT_DIRECTORY )
        set(sourceFile "${sourceFile}/${target}${CMAKE_EXECUTABLE_SUFFIX}")
    endif()    

    #####################################################
    ##  determine linkDest var value
    set(linkDest "${oldPath}/${oldName}${CMAKE_EXECUTABLE_SUFFIX}")

    #####################################################
    ##  determine SYMBOLIC var value
    if(${OLD_NAME_LINK_TYPE})
            set(SYMBOLIC "SYMBOLIC")
    else()
            set(SYMBOLIC "")
    endif ()

    #####################################################
    ##  Build OLD_NAME_LINK_COMMAND var
    set(OLD_NAME_LINK_COMMAND ${CMAKE_COMMAND}
        -DSRC_PATERNS:STRING=${sourceFile}
        -DPATH_DESTINATION:STRING=${linkDest}
        -DSYMBOLIC:STRING="${SYMBOLIC}"
        -P
        ${CMAKE_SOURCE_DIR}/cmake/fsLink.cmake)
        
    #####################################################
    ##  Up OLD_NAME_LINK_COMMAND to parent scop
    set(${OutVar} ${OLD_NAME_LINK_COMMAND} PARENT_SCOPE)
	
    #####################################################
    ##  Write informations files on pairring
	set(pairing_command 
	"file(APPEND \"${_outputDir_release}/newToOld.txt\" \"${target};${oldName}\\r\")\r \
	 file(APPEND \"${oldPath}/oldToNew.txt\"            \"${oldName};${target}\\r\")")
	
	message("--------------------------------- \r${pairing_command}\r---------------------------------")
	file(WRITE toto.cmake "${pairing_command}")
	#####################################################
    ##  Build OLD_NAME_LINK_COMMAND var
    set(OLD_NAME_LINK_COMMAND ${CMAKE_COMMAND}
        -P
        ${CMAKE_SOURCE_DIR}/cmake/fsLink.cmake)

endfunction()



macro(legacy_name target OLD_BIN_NAME)

    if(${OLD_NAME_COMPATIBILITY})
    
        if(NOT EXISTS  ${OLD_NAME_PATH})
            file(MAKE_DIRECTORY ${OLD_NAME_PATH})
        endif()
        
        GenerateOldExeNameCommand( ${target} ${OLD_BIN_NAME} ${OLD_NAME_PATH} _OLD_NAME_LINK_COMMAND)
        
        add_custom_command(TARGET ${target} POST_BUILD
		COMMAND ${_OLD_NAME_LINK_COMMAND})
      
    endif()

endmacro()