
###############################################################################
## Set Options and CACHE variables for user's choices
if(WIN32)
  set(_LINK_TYPE OFF) # Default to hard links on Windows
else()
  set(_LINK_TYPE ON ) # Default to symbolic links on Unix-like systems
endif(WIN32)

# Enable/disable old name compatibility.
set(OLD_NAME_COMPATIBILITY        ON            CACHE BOOL "If this is activated, compilation will generate also binaries with old names")
# Choose between symbolic and hard links for old names.
set(OLD_NAME_LINK_TYPE            ${_LINK_TYPE} CACHE BOOL "If OLD_NAME_LINK_TYPE is activated, old names will be Symlink, else hardlink. By default HardLink on Windows and SymLink on Unix like")
# Enable/disable generating old names in a specific folder.
set(OLD_NAME_IN_SPECIFIC_FOLDER   ON            CACHE BOOL "If this is activated, old names will be generated into specific folder")
# Path to the specific folder for old names.
set(OLD_NAME_SPECIFIC_FOLDER_PATH "${CMAKE_BINARY_DIR}/compatibility" CACHE PATH "Path where compatibility binaries will be stored.")

# Determine the path for old names.
if(${OLD_NAME_IN_SPECIFIC_FOLDER})
  set(OLD_NAME_SPECIFIC_FOLDER_PATH "${OLD_NAME_SPECIFIC_FOLDER_PATH}" CACHE PATH "Path where compatibility binaries will be stored." FORCE)
  set(OLD_NAME_PATH "${OLD_NAME_SPECIFIC_FOLDER_PATH}")
else()
  set(OLD_NAME_SPECIFIC_FOLDER_PATH "${OLD_NAME_SPECIFIC_FOLDER_PATH}" CACHE INTERNAL "Path where compatibility binaries will be stored." FORCE)
  set(OLD_NAME_PATH "${CMAKE_BINARY_DIR}/bin") # Default to the binary directory's bin folder
endif()




## #############################################################################
## Generate OLD NAME LINK COMMAND for a target. To use with add_executable.
## #############################################################################
function(GenerateOldExeNameCommand target oldName oldPath OutVar_LINK_COMMAND)
  # Generates a CMake command to create a compatibility link for an executable with an old name.
  # This function is used to create a link from an old executable name to the current executable target.
  #
  # Parameters:
  #   target              - The CMake target name (the current executable).
  #   oldName             - The old name under which the executable should be accessible.
  #   oldPath             - The directory where the compatibility link will be created.
  #   OutVar_LINK_COMMAND - The output variable that will contain the generated CMake command.
  #
  # Dependencies:
  #   OLD_NAME_LINK_TYPE - A CMake cache variable that determines whether to create symbolic or hard links.
  #   CMAKE_EXECUTABLE_SUFFIX - The executable suffix for the current platform.
  #   CMAKE_SOURCE_DIR/cmake/fsLink.cmake - The path to the fsLink.cmake script that creates the link.
  #
  # Behavior:
  #   1. Determines the source file path of the executable based on the build configuration.
  #   2. Determines the destination path for the compatibility link.
  #   3. Determines the type of link (symbolic or hard) based on OLD_NAME_LINK_TYPE.
  #   4. Builds a CMake command that uses fsLink.cmake to create the compatibility link.
  #   5. Sets the generated command in the output variable OutVar_LINK_COMMAND in the parent scope.
  #   6. Writes information to pairing files using pairing_csv_${target}.cmake.

  #####################################################
  ## Determine sourceFile variable value
  get_property(GENERATOR_MULTI_CONFIG GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)

  if(${GENERATOR_MULTI_CONFIG})
    get_target_property(_outputDir_debug ${target} RUNTIME_OUTPUT_DIRECTORY_DEBUG           )
    get_target_property(_outputDir_release ${target} RUNTIME_OUTPUT_DIRECTORY_RELEASE       )
    get_target_property(_outputDir_minsize ${target} RUNTIME_OUTPUT_DIRECTORY_MINSIZEREL    )
    get_target_property(_outputDir_reldbg ${target} RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO )

    set(sourceFile $<$<CONFIG:debug>:${_outputDir_debug}/${target}${CMAKE_EXECUTABLE_SUFFIX}>
                 $<$<CONFIG:release>:${_outputDir_release}/${target}${CMAKE_EXECUTABLE_SUFFIX}>
                 $<$<CONFIG:MinSizeRel>:${_outputDir_minsize}/${target}${CMAKE_EXECUTABLE_SUFFIX}>
                 $<$<CONFIG:RelWithDebInfo>:${_outputDir_reldbg}/${target}${CMAKE_EXECUTABLE_SUFFIX}>)

    set(sourceDir ${_outputDir_release})
				 
  else()
    get_target_property(sourceFile ${target} RUNTIME_OUTPUT_DIRECTORY )
    set(sourceDir  "${sourceFile}")
    set(sourceFile "${sourceFile}/${target}${CMAKE_EXECUTABLE_SUFFIX}")
  endif()

  #####################################################
  ## Determine linkDest variable value
  set(linkDest "${oldPath}/${oldName}${CMAKE_EXECUTABLE_SUFFIX}")

  #####################################################
  ## Determine SYMBOLIC variable value
  if(${OLD_NAME_LINK_TYPE})
    set(SYMBOLIC "SYMBOLIC") # Use symbolic links
  else()
    set(SYMBOLIC "") # Use hard links
  endif()

  #####################################################
  ## Build OLD_NAME_LINK_COMMAND variable
  set(OLD_NAME_LINK_COMMAND ${CMAKE_COMMAND}
                           -DSRC_PATERNS:STRING=${sourceFile}
                           -DPATH_DESTINATION:STRING=${linkDest}
                           -DSYMBOLIC:STRING="${SYMBOLIC}"
                           -P
                           ${CMAKE_SOURCE_DIR}/cmake/fsLink.cmake)

  #####################################################
  ## Up OLD_NAME_LINK_COMMAND to parent scope
  set(${OutVar_LINK_COMMAND} ${OLD_NAME_LINK_COMMAND} PARENT_SCOPE)

  #####################################################
  ## Write information files for pairing
  set(pairing_command "file(APPEND \"${sourceDir}/newToOld.csv\" \"${target};${oldName}\\n\")\nfile(APPEND \"${oldPath}/oldToNew.csv\" \"${oldName};${target}\\n\")\n")

  file(WRITE ${CMAKE_SOURCE_DIR}/cmake/pairing/pairing_csv_${target}.cmake "${pairing_command}")
endfunction()



macro(legacy_name target OLD_BIN_NAME)
  # This macro creates a link (symbolic or hard) for an executable with an old name.
  # It is used to ensure compatibility with older versions of Anima.
  #
  # target: The CMake target name (the executable) for which to create the link.
  # OLD_BIN_NAME: The old name of the executable.
  #
  # Documentation:
  # Creates a compatibility link for an executable, allowing it to be run under an older name.
  #
  # Parameters:
  #   target      - The CMake target representing the executable.
  #   OLD_BIN_NAME - The old name under which the executable should be accessible.
  #
  # Dependencies:
  #   OLD_NAME_COMPATIBILITY - A CMake cache variable that enables or disables the creation of compatibility links.
  #   OLD_NAME_PATH - The directory where the compatibility links will be created.
  #
  # Behavior:
  #   If OLD_NAME_COMPATIBILITY is enabled, creates a link from OLD_BIN_NAME to the executable target.
  #   The type of link (symbolic or hard) is determined by the OLD_NAME_LINK_TYPE variable.
  #   Also adds a dependency to ensure the compatibility mapping files are generated after the target is built.

  if(${OLD_NAME_COMPATIBILITY})
    # Checks if old name compatibility is enabled.

    if(NOT EXISTS ${OLD_NAME_PATH})
      # If the destination directory for the links does not exist, creates it.
      file(MAKE_DIRECTORY ${OLD_NAME_PATH})
    endif()

    # Generates the command to create the link with the old name.
    GenerateOldExeNameCommand(${target} ${OLD_BIN_NAME} ${OLD_NAME_PATH} _OLD_NAME_LINK_COMMAND)

    # Adds a custom command that will be executed after the target is built.
    # This command executes the link creation command generated previously.
    add_custom_command(TARGET ${target} POST_BUILD
                       COMMAND ${_OLD_NAME_LINK_COMMAND})

    # Adds a dependency to ensure the mapping files are created after the target is built.
    add_dependencies(Z_CREATE_CSV_COMPATIBILITY_MAPPING ${target})
  endif()
endmacro()




if(${OLD_NAME_COMPATIBILITY})
  # This if block is executed only if old name compatibility is enabled.
  # Adds a custom target rebuild mapping files.
  # This target is used as a dependency to ensure the files are created after the executables are built.

  # Defines the command to create the mapping files (newToOld.csv and oldToNew.csv).
  set(PAIRING_CSV_COMMAND ${CMAKE_COMMAND}
                         -DNEW_TO_OLD_PATH:STRING=${CMAKE_BINARY_DIR}/bin
                         -DOLD_TO_NEW_PATH:STRING=${OLD_NAME_PATH}
                         -DPAIRING_DIR:STRING=${CMAKE_SOURCE_DIR}/cmake/pairing/
                         -P
                         ${CMAKE_SOURCE_DIR}/cmake/pairing/pairing.cmake)

  # Adds a custom target to execute the mapping files creation command.
  add_custom_target(Z_CREATE_CSV_COMPATIBILITY_MAPPING ALL COMMAND ${PAIRING_CSV_COMMAND})
endif()
