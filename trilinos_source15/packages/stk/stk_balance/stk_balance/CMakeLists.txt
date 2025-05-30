# Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
# Solutions of Sandia, LLC (NTESS). Under the terms of Contract
# DE-NA0003525 with NTESS, the U.S. Government retains certain rights
# in this software.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
# 
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
# 
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
# 
#     * Neither the name of NTESS nor the names of its contributors
#       may be used to endorse or promote products derived from this
#       software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 


#
# A) Package-specific configuration options
#

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "")
SET(SOURCES "")

SET(HEADERS_IMPL "")
SET(SOURCES_IMPL "")

#
# Core headers
#

FILE(GLOB HEADERS *.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/mesh/*.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/setup/*.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/io/*.hpp)
FILE(GLOB SOURCES *.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/mesh/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/setup/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/io/*.cpp)

FILE(GLOB HEADERS_INTERNAL ${CMAKE_CURRENT_SOURCE_DIR}/internal/*.hpp)
FILE(GLOB SOURCES_INTERNAL ${CMAKE_CURRENT_SOURCE_DIR}/internal/*.cpp)

FILE(GLOB HEADERS_M2N ${CMAKE_CURRENT_SOURCE_DIR}/m2n/*.hpp)
FILE(GLOB SOURCES_M2N ${CMAKE_CURRENT_SOURCE_DIR}/m2n/*.cpp)

FILE(GLOB HEADERS_SETUP ${CMAKE_CURRENT_SOURCE_DIR}/setup/*.hpp)
FILE(GLOB SOURCES_SETUP ${CMAKE_CURRENT_SOURCE_DIR}/setup/*.cpp)

FILE(GLOB HEADERS_SEARCH_TOLERANCE ${CMAKE_CURRENT_SOURCE_DIR}/search_tolerance/*.hpp)

FILE(GLOB HEADERS_SEARCH_TOLERANCE_ALGS ${CMAKE_CURRENT_SOURCE_DIR}/search_tolerance_algs/*.hpp)
FILE(GLOB SOURCES_SEARCH_TOLERANCE_ALGS ${CMAKE_CURRENT_SOURCE_DIR}/search_tolerance_algs/*.cpp)

FILE(GLOB SOURCES_MAIN ${CMAKE_CURRENT_SOURCE_DIR}/balance_main/*.cpp)

FILE(GLOB SOURCES_M2N_MAIN ${CMAKE_CURRENT_SOURCE_DIR}/m2n_main/*.cpp)

IF (TPL_ENABLE_MPI)
ENDIF()

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  stk_balance_lib
  NOINSTALLHEADERS ${HEADERS} ${HEADERS_INTERNAL} ${HEADERS_M2N} ${HEADERS_SETUP} ${HEADERS_SEARCH_TOLERANCE} ${HEADERS_SEARCH_TOLERANCE_ALGS}
  SOURCES ${SOURCES} ${SOURCES_INTERNAL} ${SOURCES_M2N} ${SOURCES_SETUP} ${SOURCES_SEARCH_TOLERANCE_ALGS}
  )

target_include_directories(stk_balance_lib PUBLIC
  $<BUILD_INTERFACE:${STK_TOPLEVEL_SOURCE_DIR}/stk_balance>
  $<INSTALL_INTERFACE:include>
)

TRIBITS_ADD_EXECUTABLE(
    stk_balance
    SOURCES ${SOURCES_MAIN}
    NOEXEPREFIX INSTALLABLE
    )

TRIBITS_ADD_EXECUTABLE(
    stk_balance_m2n
    SOURCES ${SOURCES_M2N_MAIN}
    NOEXEPREFIX INSTALLABLE
    )

INSTALL(FILES ${HEADERS} DESTINATION
    ${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/stk_balance)
INSTALL(FILES ${HEADERS_INTERNAL} DESTINATION
    ${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/stk_balance/internal)
INSTALL(FILES ${HEADERS_M2N} DESTINATION
    ${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/stk_balance/m2n)
INSTALL(FILES ${HEADERS_SETUP} DESTINATION
    ${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/stk_balance/setup)
INSTALL(FILES ${HEADERS_SEARCH_TOLERANCE} DESTINATION
    ${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/stk_balance/search_tolerance)
INSTALL(FILES ${HEADERS_SEARCH_TOLERANCE_ALGS} DESTINATION
    ${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/stk_balance/search_tolerance_algs)

