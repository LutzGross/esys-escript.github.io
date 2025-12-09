#----------------------------------------------------------------
# Generated CMake target import file for configuration "RELEASE".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "MueLu::muelu" for configuration "RELEASE"
set_property(TARGET MueLu::muelu APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(MueLu::muelu PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libmuelu.so.16.0.0"
  IMPORTED_SONAME_RELEASE "libmuelu.so.16"
  )

list(APPEND _cmake_import_check_targets MueLu::muelu )
list(APPEND _cmake_import_check_files_for_MueLu::muelu "${_IMPORT_PREFIX}/lib/libmuelu.so.16.0.0" )

# Import target "MueLu::muelu-adapters" for configuration "RELEASE"
set_property(TARGET MueLu::muelu-adapters APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(MueLu::muelu-adapters PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libmuelu-adapters.so.16.0.0"
  IMPORTED_SONAME_RELEASE "libmuelu-adapters.so.16"
  )

list(APPEND _cmake_import_check_targets MueLu::muelu-adapters )
list(APPEND _cmake_import_check_files_for_MueLu::muelu-adapters "${_IMPORT_PREFIX}/lib/libmuelu-adapters.so.16.0.0" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
