/*****************************************************************************
*
* Copyright (c) 2003-2022 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014-2017 by Centre for Geoscience Computing (GeoComp)
* Development from 2019 by School of Earth and Environmental Sciences
**
*****************************************************************************/

#include <boost/python.hpp>
#include <mpi4py/mpi4py.h>

#include <pythonMPI/esys_python.h>


// Python interface
BOOST_PYTHON_MODULE(pythonMPIcpp)
{
	import_mpi4py();
	boost::python::def("test_mpi_program", esys_pythonMPI::pythonMPIWrapper);
}