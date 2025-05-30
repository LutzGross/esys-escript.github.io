#
# A) Package-specific configuration options
#

# Matain backward compatibility w.r.t. old Teuchos TPL dependencies.
# This list of TPLs is as of 11/01/2012.  Future TPL Dependencies
# should be listed and used as HAVE_TEUCHOS<SUBPACKAGE_NAME>_XXX which
# this is handled automatically by TriBITS once 
# #cmakedefine HAVE_TEUCHOS<SUBPACKAGE_NAME>_XXX is added to config.h file.
SET(OldTeuchos_SUBPACKAGE_NAMES BinUtils Boost MPI ARPREC QD QT)
FOREACH(OldTeuchos_SUBPACKAGE_NAME  ${OldTeuchos_SUBPACKAGE_NAMES})
  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_${OldTeuchos_SUBPACKAGE_NAME})
  IF (${PACKAGE_NAME}_ENABLE_${OldTeuchos_SUBPACKAGE_NAME})
    # Must set the preprocessor directives that C++ code will use
    STRING(TOUPPER "${OldTeuchos_SUBPACKAGE_NAME}" OldTeuchos_SUBPACKAGE_UNAME)
    SET(HAVE_TEUCHOS_${OldTeuchos_SUBPACKAGE_UNAME} TRUE)
    # Some client CMakeLists.txt code keys off of the old cmake variable names
    # so we have to set them as global variables.
    SET(Teuchos_ENABLE_${OldTeuchos_SUBPACKAGE_NAME} ON CACHE BOOL INTERNAL)
  ENDIF()
ENDFOREACH()

# Looking forward to TeuchosParameterList customers ...
ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_${PARENT_PACKAGE_NAME}ParameterList)
IF (${PROJECT_NAME}_ENABLE_${PARENT_PACKAGE_NAME}ParameterList)
  SET(HAVE_TEUCHOS_PARAMETERLIST TRUE)
ENDIF()

# FIXME: Remove this hard-coded setting.  We always have C++11 now.
SET(HAVE_TEUCHOSCORE_CXX11 ON)

SET(YAML_CPP_DEF 0)
IF(${PACKAGE_NAME}_ENABLE_yaml-cpp)
 SET(YAML_CPP_DEF 1)
ENDIF()

TRIBITS_CONFIGURE_FILE(${PARENT_PACKAGE_NAME}_config.h)
TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and include directories)
#

SET(HEADERS "")
SET(SOURCES "")

# Set the base Trilinos directory in order to find Trilinos_version.h
ASSERT_DEFINED(${${PACKAGE_NAME}_PARENT_REPOSITORY}_BINARY_DIR)
TRIBITS_INCLUDE_DIRECTORIES(${${${PACKAGE_NAME}_PARENT_REPOSITORY}_BINARY_DIR})

# Set headers (for IDEs and installing)

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

APPEND_SET(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/${PARENT_PACKAGE_NAME}_config.h
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

# NOTE: Above, we add the XXX_config.h file so tht MS Visual Studio
# will know to include this file in the dependency analysis

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.h)
APPEND_GLOB(HEADERS ${DIR}/*.c)
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)

LIST(REMOVE_ITEM  HEADERS
  # Not to be built yet
  ${CMAKE_CURRENT_SOURCE_DIR}/Teuchos_SharedAccessPolicy.hpp
  )

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Boost)
IF (NOT ${PACKAGE_NAME}_ENABLE_Boost)

  LIST(REMOVE_ITEM  HEADERS
      ${CMAKE_CURRENT_SOURCE_DIR}/Teuchos_RCPBoostSharedPtrConversions.hpp
      ${CMAKE_CURRENT_SOURCE_DIR}/Teuchos_RCPBoostSharedPtrConversionsDecl.hpp
      )

ENDIF()

ASSERT_DEFINED(TPL_ENABLE_Valgrind)
IF(TPL_ENABLE_Valgrind)
  TRIBITS_INCLUDE_DIRECTORIES(${TPL_Valgrind_INCLUDE_DIRS})
ENDIF()

#
# C) Define the targets for package's library(s)
#

IF (CMAKE_CXX_COMPILER_ID STREQUAL GNU)
  # Must suppress strict aliasing warnings for this file!
  SET_SOURCE_FILES_PROPERTIES(
    Teuchos_stacktrace.cpp
    PROPERTIES COMPILE_FLAGS -fno-strict-aliasing
    )
ENDIF()

TRIBITS_ADD_LIBRARY(
  teuchoscore
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEFINES -DTEUCHOSCORE_LIB_EXPORTS_MODE
  )

#
# Make a trivial change to this comment if you add a file either to
# this directory, or to the impl/ subdirectory.  Otherwise, users who 
# pull the latest version of Trilinos may get build errors, because
# their Makefile won't know to rerun CMake.
#
