#----------------------------------------------------------------
# Generated CMake target import file for configuration "RELEASE".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "EpetraExt::epetraext" for configuration "RELEASE"
set_property(TARGET EpetraExt::epetraext APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(EpetraExt::epetraext PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libepetraext.so.16.0.0"
  IMPORTED_SONAME_RELEASE "libepetraext.so.16"
  )

list(APPEND _cmake_import_check_targets EpetraExt::epetraext )
list(APPEND _cmake_import_check_files_for_EpetraExt::epetraext "${_IMPORT_PREFIX}/lib/libepetraext.so.16.0.0" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
