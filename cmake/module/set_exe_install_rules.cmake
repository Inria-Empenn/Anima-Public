macro(set_exe_install_rules
  target
  )

################################################################################
#
# Usage: set_exe_install_rules(target)
# set rules for the executable designed by the target
#
################################################################################

install(TARGETS ${target}
  RUNTIME DESTINATION bin
  BUNDLE  DESTINATION bin
  ARCHIVE DESTINATION bin
  )

endmacro()

