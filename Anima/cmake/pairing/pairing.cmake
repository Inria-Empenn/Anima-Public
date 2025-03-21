## MANDATORY VARIABLES
# NEW_TO_OLD_PATH: Path to the directory containing new executable files.
# OLD_TO_NEW_PATH: Path to the directory where old-to-new mapping files will be stored.
# PAIRING_DIR: Path to the directory containing pairing_csv_*.cmake files.

# Remove existing mapping files if they exist.
if(EXISTS ${NEW_TO_OLD_PATH}/newToOld.csv)
  file(REMOVE ${NEW_TO_OLD_PATH}/newToOld.csv)
endif()
if(EXISTS ${OLD_TO_NEW_PATH}/oldToNew.csv)
  file(REMOVE ${OLD_TO_NEW_PATH}/oldToNew.csv)
endif()

# Get a list of all executable files in the NEW_TO_OLD_PATH directory.
file(GLOB new_exe_files "${NEW_TO_OLD_PATH}/*")

# Extract the executable names from the file paths.
foreach(exe_file ${new_exe_files})
  cmake_path(GET exe_file STEM exe_name)
  list(APPEND exe_names ${exe_name})
endforeach()

# Remove duplicate executable names from the list (due to debug info files or other stuff like that).
list(REMOVE_DUPLICATES exe_names)

# Process each executable name.
foreach(exe_name ${exe_names})
  # Include the corresponding pairing_csv_*.cmake file if it exists.
  if(EXISTS "${PAIRING_DIR}/pairing_csv_${exe_name}.cmake")
    include("${PAIRING_DIR}/pairing_csv_${exe_name}.cmake")
  endif()
endforeach()