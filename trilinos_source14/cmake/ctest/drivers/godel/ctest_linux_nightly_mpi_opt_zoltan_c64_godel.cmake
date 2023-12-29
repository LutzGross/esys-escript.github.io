# @HEADER
# ************************************************************************
#
#            Trilinos: An Object-Oriented Solver Framework
#                 Copyright (2001) Sandia Corporation
#
#
# Copyright (2001) Sandia Corporation. Under the terms of Contract
# DE-AC04-94AL85000, there is a non-exclusive license for use of this
# work by or on behalf of the U.S. Government.  Export of this program
# may require a license from the United States Government.
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# NOTICE:  The United States Government is granted for itself and others
# acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
# license in this data to reproduce, prepare derivative works, and
# perform publicly and display publicly.  Beginning five (5) years from
# July 25, 2001, the United States Government is granted for itself and
# others acting on its behalf a paid-up, nonexclusive, irrevocable
# worldwide license in this data to reproduce, prepare derivative works,
# distribute copies to the public, perform publicly and display
# publicly, and to permit others to do so.
#
# NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
# OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
# ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
# RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
# INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
# THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
#
# ************************************************************************
# @HEADER


INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.godel.gcc.cmake")

#
# Test only Zoltan with the C compiler with MPI.  Here the only reason
# that we are testing only Zoltan is that is the only package in
# Trilinos that can build with only C.
#   Use 64 bit global IDs in Zoltan.
#

SET(COMM_TYPE MPI)
SET(BUILD_TYPE RELEASE)
SET(BUILD_DIR_NAME ZOLTAN_C_64)
SET(CTEST_TEST_TYPE Nightly)
SET(CTEST_TEST_TIMEOUT 1200)

#SET(CTEST_DO_COVERAGE_TESTING TRUE)
#SET(CTEST_DO_MEMORY_TESTING TRUE)

SET(Trilinos_PACKAGES Zoltan)

SET( EXTRA_CONFIGURE_OPTIONS
  "-DMPI_EXEC_MAX_NUMPROCS:STRING=11"
  "-DTrilinos_ENABLE_CXX:BOOL=OFF"
  "-DTrilinos_ENABLE_ForTrilinos:BOOL=OFF"
  "-DTrilinos_ENABLE_Fortran:BOOL=OFF"
  "-DZoltan_ENABLE_ULONG_IDS:BOOL=ON"
  "-DZoltan_ENABLE_TESTS:BOOL=ON"
  "-DTrilinos_TEST_CATEGORIES:STRING=NIGHTLY"
  "-DTPL_ENABLE_ParMETIS:BOOL=ON"
  "-DParMETIS_LIBRARY_DIRS:PATH=/home/kddevin/code/ParMETIS3_1"
  "-DTPL_ENABLE_Scotch:BOOL=ON"
  "-DScotch_INCLUDE_DIRS:PATH=/home/kddevin/code/scotch_5.1/include"
  "-DScotch_LIBRARY_DIRS:PATH=/home/kddevin/code/scotch_5.1/lib"
  "-DTPL_ENABLE_BinUtils=OFF"
  )

#
# Set the rest of the system-specific options and run the dashboard build/test
#

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()
