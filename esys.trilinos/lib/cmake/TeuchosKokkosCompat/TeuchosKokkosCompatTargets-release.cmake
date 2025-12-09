#----------------------------------------------------------------
# Generated CMake target import file for configuration "RELEASE".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "TeuchosKokkosCompat::teuchoskokkoscompat" for configuration "RELEASE"
set_property(TARGET TeuchosKokkosCompat::teuchoskokkoscompat APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(TeuchosKokkosCompat::teuchoskokkoscompat PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libteuchoskokkoscompat.so.16.0.0"
  IMPORTED_SONAME_RELEASE "libteuchoskokkoscompat.so.16"
  )

list(APPEND _cmake_import_check_targets TeuchosKokkosCompat::teuchoskokkoscompat )
list(APPEND _cmake_import_check_files_for_TeuchosKokkosCompat::teuchoskokkoscompat "${_IMPORT_PREFIX}/lib/libteuchoskokkoscompat.so.16.0.0" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
