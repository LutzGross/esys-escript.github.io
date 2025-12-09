#----------------------------------------------------------------
# Generated CMake target import file for configuration "RELEASE".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "TpetraCore::tpetraclassic" for configuration "RELEASE"
set_property(TARGET TpetraCore::tpetraclassic APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(TpetraCore::tpetraclassic PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libtpetraclassic.so.16.0.0"
  IMPORTED_SONAME_RELEASE "libtpetraclassic.so.16"
  )

list(APPEND _cmake_import_check_targets TpetraCore::tpetraclassic )
list(APPEND _cmake_import_check_files_for_TpetraCore::tpetraclassic "${_IMPORT_PREFIX}/lib/libtpetraclassic.so.16.0.0" )

# Import target "TpetraCore::tpetra" for configuration "RELEASE"
set_property(TARGET TpetraCore::tpetra APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(TpetraCore::tpetra PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libtpetra.so.16.0.0"
  IMPORTED_SONAME_RELEASE "libtpetra.so.16"
  )

list(APPEND _cmake_import_check_targets TpetraCore::tpetra )
list(APPEND _cmake_import_check_files_for_TpetraCore::tpetra "${_IMPORT_PREFIX}/lib/libtpetra.so.16.0.0" )

# Import target "TpetraCore::tpetrainout" for configuration "RELEASE"
set_property(TARGET TpetraCore::tpetrainout APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(TpetraCore::tpetrainout PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libtpetrainout.so.16.0.0"
  IMPORTED_SONAME_RELEASE "libtpetrainout.so.16"
  )

list(APPEND _cmake_import_check_targets TpetraCore::tpetrainout )
list(APPEND _cmake_import_check_files_for_TpetraCore::tpetrainout "${_IMPORT_PREFIX}/lib/libtpetrainout.so.16.0.0" )

# Import target "TpetraCore::tpetraext" for configuration "RELEASE"
set_property(TARGET TpetraCore::tpetraext APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(TpetraCore::tpetraext PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libtpetraext.so.16.0.0"
  IMPORTED_SONAME_RELEASE "libtpetraext.so.16"
  )

list(APPEND _cmake_import_check_targets TpetraCore::tpetraext )
list(APPEND _cmake_import_check_files_for_TpetraCore::tpetraext "${_IMPORT_PREFIX}/lib/libtpetraext.so.16.0.0" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
