
#include <boost/python.hpp>

#ifdef ESYS_HAVE_MPI4PY
#include <mpi4py/mpi4py.h>
#endif

#include <mpi.h>

#include <iostream>

namespace esys_pythonMPI
{

// forward declarations
void python_mpi_test_program_cxx(MPI_Comm comm);
float python_mpi_test_program_cxx_2(MPI_Comm comm, float x);
void example_calculation_program(float x);

/////////////////////////////////////////////////////////////////
// wrapper that converts py_comm to an MPI_Comm 
static void pythonMPIWrapper(boost::python::object py_comm)
{
	#ifdef ESYS_HAVE_MPI4PY
  PyObject* py_obj = py_comm.ptr();
  MPI_Comm *comm_p = PyMPIComm_Get(py_obj);
  if (comm_p == NULL) boost::python::throw_error_already_set();
  python_mpi_test_program_cxx(*comm_p);
  #endif
}

// wrapper for floating point number test function
static float pythonMPIWrapper2(boost::python::object py_comm, float x)
{
	#ifdef ESYS_HAVE_MPI4PY
  PyObject* py_obj = py_comm.ptr();
  MPI_Comm *comm_p = PyMPIComm_Get(py_obj);
  if (comm_p == NULL) boost::python::throw_error_already_set();
  return python_mpi_test_program_cxx_2(*comm_p, x);
  #else
  return 0.0;
  #endif
}

/////////////////////////////////////////////////////////////////
// function that prints some MPI info to the console
void python_mpi_test_program_cxx(MPI_Comm comm)
{
	int size;
	int rank;

	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);

	std::cout << "Hello, World from process " << rank 
			  << " of " << size  << "." << std::endl;
}

// function that takes in and modifies a float then returns it
float python_mpi_test_program_cxx_2(MPI_Comm comm, float x)
{
		int size;
		int rank;
		MPI_Comm_size(comm, &size);
		MPI_Comm_rank(comm, &rank);
		std::cout << "Process " << rank << " of " << size << " is returning x=" << x*x << std::endl;
		return x*x;
}

} //end namespace esys_pythonMPI