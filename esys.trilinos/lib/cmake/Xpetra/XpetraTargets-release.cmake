#----------------------------------------------------------------
# Generated CMake target import file for configuration "RELEASE".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "Xpetra::xpetra" for configuration "RELEASE"
set_property(TARGET Xpetra::xpetra APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Xpetra::xpetra PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libxpetra.so.16.0.0"
  IMPORTED_SONAME_RELEASE "libxpetra.so.16"
  )

list(APPEND _cmake_import_check_targets Xpetra::xpetra )
list(APPEND _cmake_import_check_files_for_Xpetra::xpetra "${_IMPORT_PREFIX}/lib/libxpetra.so.16.0.0" )

# Import target "Xpetra::xpetra-sup" for configuration "RELEASE"
set_property(TARGET Xpetra::xpetra-sup APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Xpetra::xpetra-sup PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libxpetra-sup.so.16.0.0"
  IMPORTED_SONAME_RELEASE "libxpetra-sup.so.16"
  )

list(APPEND _cmake_import_check_targets Xpetra::xpetra-sup )
list(APPEND _cmake_import_check_files_for_Xpetra::xpetra-sup "${_IMPORT_PREFIX}/lib/libxpetra-sup.so.16.0.0" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
