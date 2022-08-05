

#include <mpi.h>

#include <iostream>

namespace esys_pythonMPI
{

// forward declarations
void python_mpi_test_program_cxx(MPI_Comm comm);
void example_calculation_program(float x);

/////////////////////////////////////////////////////////////////
// wrapper that converts py_comm to an MPI_Comm 
static void pythonMPIWrapper(boost::python::object py_comm)
{
  PyObject* py_obj = py_comm.ptr();
  MPI_Comm *comm_p = PyMPIComm_Get(py_obj);
  if (comm_p == NULL) boost::python::throw_error_already_set();
  python_mpi_test_program_cxx(*comm_p);
}

/////////////////////////////////////////////////////////////////
// function that prints some MPI info to the console
void python_mpi_test_program_cxx(MPI_Comm comm)
{
	std::cout << "inside python_mpi_test_program_cxx" << std::endl;

	int size;
	int rank;

	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);

	std::cout << "Hello, World from process " << rank 
			  << " of " << size  << "." << std::endl;
}

} //end namespace esys_pythonMPI