
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../stk_util)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../stk_search_util)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../stk_transfer_util)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/)

FILE(GLOB SOURCES *.cpp)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
    transfer_util_unit_tests
    SOURCES ${SOURCES}
    TESTONLYLIBS stk_unit_main
    ARGS ""
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

