/*****************************************************************************
*
* Copyright (c) 2003-2022 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#include <boost/python.hpp>

#include <pythonMPI/esys_python.h>


#ifdef ESYS_HAVE_MPI4PY
// Python interface
BOOST_PYTHON_MODULE(pythonMPIcpp)
{
	
	import_mpi4py();
	boost::python::def("test_mpi_program", esys_pythonMPI::test_pythonMPIWrapper);
	boost::python::def("test_mpi_program2", esys_pythonMPI::test_pythonMPIWrapper2, boost::python::args("comm", "x"));
	
}
#endif