#----------------------------------------------------------------
# Generated CMake target import file for configuration "RELEASE".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "Amesos2::amesos2" for configuration "RELEASE"
set_property(TARGET Amesos2::amesos2 APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Amesos2::amesos2 PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libamesos2.so.16.0.0"
  IMPORTED_SONAME_RELEASE "libamesos2.so.16"
  )

list(APPEND _cmake_import_check_targets Amesos2::amesos2 )
list(APPEND _cmake_import_check_files_for_Amesos2::amesos2 "${_IMPORT_PREFIX}/lib/libamesos2.so.16.0.0" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
