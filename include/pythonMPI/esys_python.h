
#include <boost/python.hpp>
#include <escript/EsysMPI.h>

#ifdef ESYS_HAVE_MPI4PY
#include <mpi4py/mpi4py.h>
#endif


#include <iostream>

namespace esys_pythonMPI
{

// forward declarations
void python_mpi_test_program_cxx(MPI_Comm comm);
float python_mpi_test_program_cxx_2(MPI_Comm comm, float x);
void example_calculation_program(float x);

/////////////////////////////////////////////////////////////////
// wrapper that converts py_comm to an MPI_Comm 
#ifdef ESYS_HAVE_MPI4PY


static void test_pythonMPIWrapper(boost::python::object py_comm)
{

  PyObject* py_obj = py_comm.ptr();
  MPI_Comm *comm_p = PyMPIComm_Get(py_obj);
  if (comm_p == NULL) boost::python::throw_error_already_set();
  python_mpi_test_program_cxx(*comm_p);

}

// wrapper for floating point number test function
static float test_pythonMPIWrapper2(boost::python::object py_comm, float x)
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
#endif
/////////////////////////////////////////////////////////////////
// function that prints some MPI info to the console
void python_mpi_test_program_cxx(MPI_Comm comm)
{

	int size=1;
	int rank=0;

    #ifdef ESYS_HAVE_MPI4PY
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);
    #endif

	std::cout << "Hello, World from process " << rank 
			  << " of " << size  << "." << std::endl;
}

// function that takes in and modifies a float then returns it
float python_mpi_test_program_cxx_2(MPI_Comm comm, float x)
{
		int size=1;
		int rank=0;
		#ifdef ESYS_HAVE_MPI4PY
		MPI_Comm_size(comm, &size);
		MPI_Comm_rank(comm, &rank);
		#endif

		std::cout << "Process " << rank << " of " << size << " is returning x=" << x*x << std::endl;
		return x*x;
}

} //end namespace esys_pythonMPI