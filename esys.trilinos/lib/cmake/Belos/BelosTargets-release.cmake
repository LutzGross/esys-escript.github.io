#----------------------------------------------------------------
# Generated CMake target import file for configuration "RELEASE".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "Belos::belos" for configuration "RELEASE"
set_property(TARGET Belos::belos APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Belos::belos PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libbelos.so.16.0.0"
  IMPORTED_SONAME_RELEASE "libbelos.so.16"
  )

list(APPEND _cmake_import_check_targets Belos::belos )
list(APPEND _cmake_import_check_files_for_Belos::belos "${_IMPORT_PREFIX}/lib/libbelos.so.16.0.0" )

# Import target "Belos::belosepetra" for configuration "RELEASE"
set_property(TARGET Belos::belosepetra APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Belos::belosepetra PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libbelosepetra.so.16.0.0"
  IMPORTED_SONAME_RELEASE "libbelosepetra.so.16"
  )

list(APPEND _cmake_import_check_targets Belos::belosepetra )
list(APPEND _cmake_import_check_files_for_Belos::belosepetra "${_IMPORT_PREFIX}/lib/libbelosepetra.so.16.0.0" )

# Import target "Belos::belostpetra" for configuration "RELEASE"
set_property(TARGET Belos::belostpetra APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Belos::belostpetra PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libbelostpetra.so.16.0.0"
  IMPORTED_SONAME_RELEASE "libbelostpetra.so.16"
  )

list(APPEND _cmake_import_check_targets Belos::belostpetra )
list(APPEND _cmake_import_check_files_for_Belos::belostpetra "${_IMPORT_PREFIX}/lib/libbelostpetra.so.16.0.0" )

# Import target "Belos::belosxpetra" for configuration "RELEASE"
set_property(TARGET Belos::belosxpetra APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Belos::belosxpetra PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libbelosxpetra.so.16.0.0"
  IMPORTED_SONAME_RELEASE "libbelosxpetra.so.16"
  )

list(APPEND _cmake_import_check_targets Belos::belosxpetra )
list(APPEND _cmake_import_check_files_for_Belos::belosxpetra "${_IMPORT_PREFIX}/lib/libbelosxpetra.so.16.0.0" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
