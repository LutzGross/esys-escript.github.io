#----------------------------------------------------------------
# Generated CMake target import file for configuration "RELEASE".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "Galeri::galeri-epetra" for configuration "RELEASE"
set_property(TARGET Galeri::galeri-epetra APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Galeri::galeri-epetra PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libgaleri-epetra.so.16.0.0"
  IMPORTED_SONAME_RELEASE "libgaleri-epetra.so.16"
  )

list(APPEND _cmake_import_check_targets Galeri::galeri-epetra )
list(APPEND _cmake_import_check_files_for_Galeri::galeri-epetra "${_IMPORT_PREFIX}/lib/libgaleri-epetra.so.16.0.0" )

# Import target "Galeri::galeri-xpetra" for configuration "RELEASE"
set_property(TARGET Galeri::galeri-xpetra APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Galeri::galeri-xpetra PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libgaleri-xpetra.so.16.0.0"
  IMPORTED_SONAME_RELEASE "libgaleri-xpetra.so.16"
  )

list(APPEND _cmake_import_check_targets Galeri::galeri-xpetra )
list(APPEND _cmake_import_check_files_for_Galeri::galeri-xpetra "${_IMPORT_PREFIX}/lib/libgaleri-xpetra.so.16.0.0" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
