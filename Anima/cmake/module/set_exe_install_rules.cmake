# set_outputs_and_install_rules
macro(set_exe_install_rules
  target
  )

################################################################################
#
# Usage: set_exe_install_rules(target)
# set rules for the executable designed by the target
#
################################################################################
get_property(GENERATOR_MULTI_CONFIG GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)

if(${GENERATOR_MULTI_CONFIG})
  set_target_properties( ${target} PROPERTIES RUNTIME_OUTPUT_DIRECTORY_DEBUG          ${CMAKE_BINARY_DIR}/bin)
  set_target_properties( ${target} PROPERTIES RUNTIME_OUTPUT_DIRECTORY_RELEASE        ${CMAKE_BINARY_DIR}/bin)
  set_target_properties( ${target} PROPERTIES RUNTIME_OUTPUT_DIRECTORY_MINSIZEREL     ${CMAKE_BINARY_DIR}/bin)
  set_target_properties( ${target} PROPERTIES RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO ${CMAKE_BINARY_DIR}/bin)
else()
  set_target_properties( ${target} PROPERTIES RUNTIME_OUTPUT_DIRECTORY                ${CMAKE_BINARY_DIR}/bin)
endif()


if(${ARGC} GREATER 1)
  legacy_name(${PROJECT_NAME} "${ARGV1}")
endif()


install(TARGETS ${target}
  RUNTIME DESTINATION bin
  BUNDLE  DESTINATION bin
  ARCHIVE DESTINATION bin
  )

endmacro()

