#----------------------------------------------------------------
# Generated CMake target import file for configuration "RELEASE".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "Ifpack2::ifpack2" for configuration "RELEASE"
set_property(TARGET Ifpack2::ifpack2 APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Ifpack2::ifpack2 PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libifpack2.so.16.0.0"
  IMPORTED_SONAME_RELEASE "libifpack2.so.16"
  )

list(APPEND _cmake_import_check_targets Ifpack2::ifpack2 )
list(APPEND _cmake_import_check_files_for_Ifpack2::ifpack2 "${_IMPORT_PREFIX}/lib/libifpack2.so.16.0.0" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
