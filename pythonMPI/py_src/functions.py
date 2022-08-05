from .pythonMPIcpp import *


def test_mpi(py_comm):
    test_mpi_program(py_comm)

def test_mpi_2(py_comm, x):
    return test_mpi_program2(py_comm, x)