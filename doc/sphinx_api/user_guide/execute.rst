.. _EXECUTION:

==============================
Execution of an escript Script
==============================

This chapter describes how to run **esys.escript** scripts using the ``run-escript`` launcher,
including options for parallel execution with OpenMP and MPI.

.. note::

   This chapter is being converted from LaTeX. For the complete content with all options
   and detailed examples, please refer to the `User Guide PDF <../../user/user.pdf>`_.

Overview
========

A typical way of starting your escript script ``myscript.py`` is with
the ``run-escript`` command. This command was renamed from ``escript``
(used in previous releases) to avoid clashing with an unrelated program
installed by default on some systems.

Basic Usage
-----------

To run your script:

.. code-block:: bash

   run-escript myscript.py

For interactive mode (useful for debugging):

.. code-block:: bash

   run-escript -i myscript.py

This will execute ``myscript.py`` and when it completes (or an error occurs),
a Python prompt will be provided. Press Ctrl-D to exit.

Parallel Execution
==================

Using OpenMP (Threading)
------------------------

To run the script using multiple threads (requires escript compiled with OpenMP support):

.. code-block:: bash

   run-escript -t 4 myscript.py

This sets the ``OMP_NUM_THREADS`` environment variable to control the number of threads.

Using MPI
---------

To run the script using MPI with multiple processes:

.. code-block:: bash

   run-escript -p 8 myscript.py

Hybrid OpenMP + MPI
-------------------

For multi-core processors or shared memory architectures, you can combine MPI and threading:

.. code-block:: bash

   run-escript -p 8 -t 4 myscript.py

This runs 8 MPI processes with 4 threads each.

Distributed Computing
---------------------

For supercomputers or clusters, distribute workload over multiple nodes:

.. code-block:: bash

   # 8 nodes with 4 MPI processes per node
   run-escript -n 8 -p 4 myscript.py

   # 8 nodes, 2 processes per node, 4 threads per process
   run-escript -n 8 -p 2 -t 4 myscript.py

Command Line Options
====================

The general form of the ``run-escript`` launcher is:

.. code-block:: text

   run-escript [options] [script.py] [ARGS]

**Common Options:**

``-n nn``
    Number of nodes to use

``-p np``
    Number of MPI processes (per node if -n is set)

``-t nt``
    Number of OpenMP threads

``-f hostfile``
    MPI host file

``-i``
    Interactive mode - drop to Python prompt after script

``-v``
    Print version information

``-h``
    Show help message

``-e``
    Show escript environment variables

``-b``
    Run in batch mode (no input)

``-m tool``
    Use specified memory debugging tool

For Complete Reference
======================

For the complete list of all command-line options and detailed usage examples,
please see the `User Guide PDF <../../user/user.pdf>`_.

See Also
========

* :doc:`mpi4py` - Using mpi4py with escript for sub-communicators
