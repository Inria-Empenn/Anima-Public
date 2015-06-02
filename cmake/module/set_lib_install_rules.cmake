macro(set_lib_install_rules
  project_name
  )

################################################################################
#
# Usage: set_lib_install_rules(target)
# set rules for the library designed by the target
#
################################################################################

install(TARGETS ${project_name}
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  )

endmacro()

