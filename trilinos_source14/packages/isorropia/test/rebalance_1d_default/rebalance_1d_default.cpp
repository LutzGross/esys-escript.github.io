//@HEADER
//************************************************************************
//
//              Isorropia: Partitioning and Load Balancing Package
//                Copyright (2006) Sandia Corporation
//
//Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
//license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//************************************************************************
//@HEADER

#include <Isorropia_ConfigDefs.hpp>
#include <Isorropia_Exception.hpp>
#include <Isorropia_Epetra.hpp>
#include <Isorropia_EpetraPartitioner.hpp>
#include <Isorropia_EpetraRedistributor.hpp>

#include <Teuchos_CommandLineProcessor.hpp>

#include <ispatest_utils.hpp>
#include <ispatest_epetra_utils.hpp>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifdef HAVE_EPETRA
#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <Epetra_LinearProblem.h>
#endif

/* Isorropia now depends on Epetra, so I have replaced the
   below code with the lines below the commented block. JW
#if defined(HAVE_MPI) && defined(HAVE_EPETRA)
#include <Epetra_MpiComm.h>
#elif HAVE_EPETRA
#include <Epetra_SerialComm.h>
#endif*/

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#ifdef HAVE_EPETRA
Epetra_CrsGraph* create_epetra_test_graph_1(int numProcs,
                                            int localProc,
                                            bool verbose);

Epetra_CrsMatrix* create_epetra_test_matrix_1(int numProcs,
                                              int localProc,
                                              bool verbose);

bool test_rebalance_epetra_linproblem(int numProcs, int localProc, bool verbose);
bool test_rebalance_epetra_linproblem2(int numProcs, int localProc, bool verbose);
bool test_rebalance_epetra_crsmatrix(int numProcs, int localProc, bool verbose);
bool test_rebalance_epetra_rowmatrix(int numProcs, int localProc, bool verbose);
bool test_rebalance_epetra_graph(int numProcs, int localProc, bool verbose);
#endif

int main(int argc, char** argv) {
#ifdef HAVE_MPI

  bool verbose = false;
  int numProcs = 1;
  int localProc = 0;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &localProc);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  try {
    verbose = ispatest::set_verbose(localProc, argc, argv);
  }
  catch(std::exception& exc) {
    std::cout << "err, setting verbosity: " << exc.what() << std::endl;
    std::cout << "End Result: TEST FAILED" << std::endl;
    MPI_Finalize();
    return(-1);
  }

  // Only print on proc 0
  if (localProc != 0) {
    verbose = false;
  }

  bool test_passed = true;
#ifdef HAVE_EPETRA
  test_passed = test_passed &&
    test_rebalance_epetra_graph(numProcs, localProc, verbose);

  test_passed = test_passed &&
    test_rebalance_epetra_crsmatrix(numProcs, localProc, verbose);

  test_passed = test_passed &&
    test_rebalance_epetra_rowmatrix(numProcs, localProc, verbose);

  test_passed = test_passed &&
    test_rebalance_epetra_linproblem(numProcs, localProc, verbose);

  test_passed = test_passed &&
    test_rebalance_epetra_linproblem2(numProcs, localProc, verbose);

#else
  std::cout << "rebalance_1d_default main: currently can only test "
         << "rebalancing with Epetra enabled." << std::endl;
  test_passed = false;
#endif

  if (test_passed && verbose) {
    std::cout << "rebalance_1d_default main: tests passed."<<std::endl;
  }

  MPI_Finalize();

#else
  std::cout << "rebalance_1d_default: don't have MPI, can't run test."
            << std::endl;
#endif

  return(0);
}

#ifdef HAVE_EPETRA
//-------------------------------------------------------------------
bool test_rebalance_epetra_crsmatrix(int numProcs, int localProc, bool verbose)
{
  bool test_passed = false;
#ifndef HAVE_MPI
  return(test_passed);
#endif

  Epetra_CrsMatrix* input_matrix =
    create_epetra_test_matrix_1(numProcs, localProc, verbose);

  //Call the Isorropia rebalance method allowing default behavior. (I.e.,
  //we won't specify any parameters, weights, etc.) Default behavior should
  //be to balance the matrix so that the number of nonzeros on each processor
  //is roughly equal. i.e., by default, weights for each row are assumed to
  //be the number of nonzeros in that row.

  Teuchos::ParameterList paramlist;

  Epetra_CrsMatrix *balanced_matrix;

  try {
    if (verbose) {
      std::cout << " calling Isorropia::createBalancedCopy(Epetra_CrsMatrix)..."
                << std::endl;
    }

    balanced_matrix =
      Isorropia::Epetra::createBalancedCopy(*input_matrix,paramlist);
  }
  catch(std::exception& exc) {
    std::cout << "caught exception: " << exc.what() << std::endl;
    return(false);
  }

  //Now check the result matrix and make sure that the number of nonzeros
  //is indeed equal on each processor. (We constructed the input matrix
  //so that a correct rebalancing would result in the same number of
  //nonzeros being on each processor.)
  const Epetra_Map& bal_rowmap = balanced_matrix->RowMap();
  int bal_local_num_rows = bal_rowmap.NumMyElements();

  //count the local nonzeros.
  if (verbose) {
    std::cout << " counting local nnz for balanced matrix..." << std::endl;
  }

  int num_nonzeros = 0;
  for(int i=0; i<bal_local_num_rows; ++i) {
    num_nonzeros += balanced_matrix->NumMyEntries(i);
  }

  const Epetra_Comm& comm = input_matrix->Comm();

  int global_num_nonzeros;
  comm.SumAll(&num_nonzeros, &global_num_nonzeros, 1);

  int avg_nnz_per_proc = global_num_nonzeros/numProcs;

  if (verbose) {
    std::cout << " making sure local nnz ("<<num_nonzeros
             <<") is the same on every proc...\n" << std::endl;
  }

  if (num_nonzeros == avg_nnz_per_proc) test_passed = true;

  int local_int_result = test_passed ? 1 : 0;
  int global_int_result;
  comm.MinAll(&local_int_result, &global_int_result, 1);

  delete input_matrix;

  test_passed = global_int_result==1 ? true : false;

  if (!test_passed && verbose) {
    std::cout << "test FAILED!" << std::endl;
  }

  delete balanced_matrix;

  return(test_passed);
}

//-------------------------------------------------------------------
bool test_rebalance_epetra_rowmatrix(int numProcs, int localProc, bool verbose)
{
  bool test_passed = false;
#ifndef HAVE_MPI
  return(test_passed);
#endif

  Epetra_CrsMatrix* input_matrix =
    create_epetra_test_matrix_1(numProcs, localProc, verbose);

  Epetra_RowMatrix* input_rowmatrix = input_matrix;

  //Call the Isorropia rebalance method allowing default behavior. (I.e.,
  //we won't specify any parameters, weights, etc.) Default behavior should
  //be to balance the matrix so that the number of nonzeros on each processor
  //is roughly equal. i.e., by default, weights for each row are assumed to
  //be the number of nonzeros in that row.

  Epetra_RowMatrix *balanced_matrix;
  try {
    if (verbose) {
      std::cout << " calling Isorropia::createBalancedCopy(Epetra_RowMatrix)..."
                << std::endl;
    }

    balanced_matrix =
      Isorropia::Epetra::createBalancedCopy(*input_rowmatrix);
  }
  catch(std::exception& exc) {
    std::cout << "caught exception: " << exc.what() << std::endl;
    return(false);
  }

  //Now check the result matrix and make sure that the number of nonzeros
  //is indeed equal on each processor. (We constructed the input matrix
  //so that a correct rebalancing would result in the same number of
  //nonzeros being on each processor.)
  const Epetra_Map& bal_rowmap = balanced_matrix->RowMatrixRowMap();
  int bal_local_num_rows = bal_rowmap.NumMyElements();

  //count the local nonzeros.
  if (verbose) {
    std::cout << " counting local nnz for balanced matrix..." << std::endl;
  }

  int num_nonzeros = 0;
  int n;
  for(int i=0; i<bal_local_num_rows; ++i) {
    balanced_matrix->NumMyRowEntries(i, n);
    num_nonzeros += n;
  }

  const Epetra_Comm& comm = input_matrix->Comm();

  int global_num_nonzeros;
  comm.SumAll(&num_nonzeros, &global_num_nonzeros, 1);

  int avg_nnz_per_proc = global_num_nonzeros/numProcs;

  if (verbose) {
    std::cout << " making sure local nnz ("<<num_nonzeros
             <<") is the same on every proc...\n" << std::endl;
  }

  if (num_nonzeros == avg_nnz_per_proc) test_passed = true;

  int local_int_result = test_passed ? 1 : 0;
  int global_int_result;
  comm.MinAll(&local_int_result, &global_int_result, 1);

  delete input_matrix;

  test_passed = global_int_result==1 ? true : false;

  if (!test_passed && verbose) {
    std::cout << "test FAILED!" << std::endl;
  }

  delete balanced_matrix;

  return(test_passed);
}

//-------------------------------------------------------------------
bool test_rebalance_epetra_linproblem(int numProcs, int localProc, bool verbose)
{
  bool test_passed = false;
#ifndef HAVE_MPI
  return(test_passed);
#endif

  Epetra_CrsMatrix* input_matrix =
    create_epetra_test_matrix_1(numProcs, localProc, verbose);

  Epetra_Vector* x = new Epetra_Vector(input_matrix->RowMap());
  Epetra_Vector* b = new Epetra_Vector(input_matrix->RowMap());
  Epetra_LinearProblem problem(input_matrix, x, b);

  //Call the Isorropia rebalance method allowing default behavior. (I.e.,
  //we won't specify any parameters, weights, etc.) Default behavior should
  //be to balance the matrix so that the number of nonzeros on each processor
  //is roughly equal. i.e., by default, weights for each row are assumed to
  //be the number of nonzeros in that row.

  Epetra_LinearProblem *balanced_problem;

  try {
    if (verbose) {
      std::cout << " calling Isorropia::createBalancedCopy(Epetra_LinearProblem)..."
                << std::endl;
    }

    balanced_problem = Isorropia::Epetra::createBalancedCopy(problem);
  }
  catch(std::exception& exc) {
    std::cout << "caught exception: " << exc.what() << std::endl;
    return(false);
  }

  //Now check the result matrix and make sure that the number of nonzeros
  //is indeed equal on each processor. (We constructed the input matrix
  //so that a correct rebalancing would result in the same number of
  //nonzeros being on each processor.)
  const Epetra_Map& bal_rowmap = balanced_problem->GetMatrix()->RowMatrixRowMap();
  int bal_local_num_rows = bal_rowmap.NumMyElements();

  //count the local nonzeros.
  if (verbose) {
    std::cout << " counting local nnz for balanced matrix..." << std::endl;
  }

  int num_nonzeros = 0;
  for(int i=0; i<bal_local_num_rows; ++i) {
    int numrowentries = 0;
    balanced_problem->GetMatrix()->NumMyRowEntries(i,numrowentries);
    num_nonzeros += numrowentries;
  }

  delete balanced_problem;

  const Epetra_Comm& comm = input_matrix->Comm();

  int global_num_nonzeros;
  comm.SumAll(&num_nonzeros, &global_num_nonzeros, 1);

  int avg_nnz_per_proc = global_num_nonzeros/numProcs;

  if (verbose) {
    std::cout << " making sure local nnz ("<<num_nonzeros
             <<") is the same on every proc...\n" << std::endl;
  }

  if (num_nonzeros == avg_nnz_per_proc) test_passed = true;

  int local_int_result = test_passed ? 1 : 0;
  int global_int_result;
  comm.MinAll(&local_int_result, &global_int_result, 1);

  test_passed = global_int_result==1 ? true : false;

  if (!test_passed && verbose) {
    std::cout << "test FAILED!" << std::endl;
  }

  delete input_matrix;
  delete x;
  delete b;

  return(test_passed);
}

//-------------------------------------------------------------------
bool test_rebalance_epetra_linproblem2(int numProcs, int localProc, bool verbose)
{
  bool test_passed = false;
#ifndef HAVE_MPI
  return(test_passed);
#endif

  Epetra_CrsMatrix* input_matrix =
    create_epetra_test_matrix_1(numProcs, localProc, verbose);

  Epetra_Vector* x = new Epetra_Vector(input_matrix->RowMap());
  Epetra_Vector* b = new Epetra_Vector(input_matrix->RowMap());
  Epetra_LinearProblem problem(input_matrix, x, b);

  Teuchos::ParameterList paramlist;

  //Wrap a RefCountPtr around the matrix graph, and specify 'false', meaning
  //that the RefCountPtr will not take ownership of the graph (will not
  //delete it).
  Teuchos::RefCountPtr<const Epetra_CrsGraph> graph =
    Teuchos::rcp( &(input_matrix->Graph()), false);

  Teuchos::RefCountPtr<Isorropia::Epetra::Partitioner> partitioner =
    Teuchos::rcp(new Isorropia::Epetra::Partitioner(graph, paramlist));

  Isorropia::Epetra::Redistributor rd(partitioner);

  Teuchos::RefCountPtr<Epetra_CrsMatrix> bal_matrix =
    rd.redistribute(*(problem.GetMatrix()));

  Teuchos::RefCountPtr<Epetra_MultiVector> bal_x =
    rd.redistribute(*(problem.GetLHS()));

  Teuchos::RefCountPtr<Epetra_MultiVector> bal_b =
    rd.redistribute(*(problem.GetRHS()));

  Teuchos::RefCountPtr<Epetra_LinearProblem> balanced_problem =
    Teuchos::rcp(new Epetra_LinearProblem(bal_matrix.get(),
					  bal_x.get(), bal_b.get()));

  //Now check the result matrix and make sure that the number of nonzeros
  //is indeed equal on each processor. (We constructed the input matrix
  //so that a correct rebalancing would result in the same number of
  //nonzeros being on each processor.)
  const Epetra_Map& bal_rowmap = balanced_problem->GetMatrix()->RowMatrixRowMap();
  int bal_local_num_rows = bal_rowmap.NumMyElements();

  //count the local nonzeros.
  if (verbose) {
    std::cout << "test_rebalance_epetra_linproblem: " << std::endl;
    std::cout << " counting local nnz for balanced matrix..." << std::endl;
  }

  int num_nonzeros = 0;
  for(int i=0; i<bal_local_num_rows; ++i) {
    int numrowentries = 0;
    balanced_problem->GetMatrix()->NumMyRowEntries(i,numrowentries);
    num_nonzeros += numrowentries;
  }

  const Epetra_Comm& comm = input_matrix->Comm();

  int global_num_nonzeros;
  comm.SumAll(&num_nonzeros, &global_num_nonzeros, 1);

  int avg_nnz_per_proc = global_num_nonzeros/numProcs;

  if (verbose) {
    std::cout << " making sure local nnz ("<<num_nonzeros
             <<") is the same on every proc...\n" << std::endl;
  }

  if (num_nonzeros == avg_nnz_per_proc) test_passed = true;

  int local_int_result = test_passed ? 1 : 0;
  int global_int_result;
  comm.MinAll(&local_int_result, &global_int_result, 1);

  test_passed = global_int_result==1 ? true : false;

  if (!test_passed && verbose) {
    std::cout << "test FAILED!" << std::endl;
  }

  delete input_matrix;
  delete x;
  delete b;

  return(test_passed);
}

//-------------------------------------------------------------------
bool test_rebalance_epetra_graph(int numProcs, int localProc, bool verbose)
{
  bool test_passed = false;
#ifndef HAVE_MPI
  return(test_passed);
#endif

  Epetra_CrsGraph* input_graph =
    create_epetra_test_graph_1(numProcs, localProc, verbose);

  //Call the Isorropia rebalance method allowing default behavior. (I.e.,
  //we won't specify any parameters, weights, etc.) Default behavior should
  //be to balance the graph so that the number of nonzeros on each processor
  //is roughly equal. i.e., by default, weights for each row are assumed to
  //be the number of nonzeros in that row.

  Teuchos::ParameterList paramlist;
  Epetra_CrsGraph * balanced_graph;
  try {
    if (verbose) {
      std::cout << " calling Isorropia::createBalancedCopy(Epetra_CrsGraph)..."
                << std::endl;
    }

    balanced_graph =
      Isorropia::Epetra::createBalancedCopy(*input_graph, paramlist);
  }
  catch(std::exception& exc) {
    std::cout << "caught exception: " << exc.what() << std::endl;
    return(false);
  }

  //Now check the result graph and make sure that the number of nonzeros
  //is indeed equal on each processor. (We constructed the input graph
  //so that a correct rebalancing would result in the same number of
  //nonzeros being on each processor.)
  const Epetra_BlockMap& bal_rowmap = balanced_graph->RowMap();
  int bal_local_num_rows = bal_rowmap.NumMyElements();

  //count the local nonzeros.
  if (verbose) {
    std::cout << " counting local nnz for balanced graph..." << std::endl;
  }

  int num_nonzeros = 0;
  for(int i=0; i<bal_local_num_rows; ++i) {
    num_nonzeros += balanced_graph->NumMyIndices(i);
  }

  const Epetra_Comm& comm = input_graph->Comm();

  int global_num_nonzeros;
  comm.SumAll(&num_nonzeros, &global_num_nonzeros, 1);

  int avg_nnz_per_proc = global_num_nonzeros/numProcs;

  if (verbose) {
    std::cout << " making sure local nnz ("<<num_nonzeros
             <<") is the same on every proc...\n" << std::endl;
  }

  if (num_nonzeros == avg_nnz_per_proc) test_passed = true;

  int local_int_result = test_passed ? 1 : 0;
  int global_int_result;
  comm.MinAll(&local_int_result, &global_int_result, 1);

  delete input_graph;

  test_passed = global_int_result==1 ? true : false;

  if (!test_passed && verbose) {
    std::cout << "test FAILED!" << std::endl;
  }

  delete balanced_graph;

  return(test_passed);
}

Epetra_CrsMatrix* create_epetra_test_matrix_1(int numProcs,
                                              int localProc,
                                              bool verbose)
{
  if (verbose) {
    std::cout << " creating Epetra_CrsMatrix with un-even distribution..."
            << std::endl;
  }

  //create an Epetra_CrsMatrix with rows spread un-evenly over
  //processors.
#ifdef HAVE_MPI
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  int local_num_rows = 360;
  int nnz_per_row = local_num_rows/4+1;
  int global_num_rows = numProcs*local_num_rows;

  int mid_proc = numProcs/2;
  bool num_procs_even = numProcs%2==0 ? true : false;

  int adjustment = local_num_rows/2;

  //adjust local_num_rows so that it's not equal on all procs.
  if (localProc < mid_proc) {
    local_num_rows -= adjustment;
  }
  else {
    local_num_rows += adjustment;
  }

  //if numProcs is not an even number, undo the local_num_rows adjustment
  //on one proc so that the total will still be correct.
  if (localProc == numProcs-1) {
    if (num_procs_even == false) {
      local_num_rows -= adjustment;
    }
  }

  //now we're ready to create a row-map.
  Epetra_Map rowmap(global_num_rows, local_num_rows, 0, comm);

  //create a matrix
  Epetra_CrsMatrix* input_matrix =
    new Epetra_CrsMatrix(Copy, rowmap, nnz_per_row);

  int err = ispatest::fill_matrix(*input_matrix, nnz_per_row, verbose);
  if (err != 0) {
    delete input_matrix;
    return(0);
  }

  return(input_matrix);
}

Epetra_CrsGraph* create_epetra_test_graph_1(int numProcs,
                                              int localProc,
                                              bool verbose)
{
  if (verbose) {
    std::cout << " creating Epetra_CrsGraph with un-even distribution..."
            << std::endl;
  }

  //create an Epetra_CrsGraph with rows spread un-evenly over
  //processors.
#ifdef HAVE_MPI
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif
  int local_num_rows = 360;
  int nnz_per_row = local_num_rows/4+1;
  int global_num_rows = numProcs*local_num_rows;

  int mid_proc = numProcs/2;
  bool num_procs_even = numProcs%2==0 ? true : false;

  int adjustment = local_num_rows/2;

  //adjust local_num_rows so that it's not equal on all procs.
  if (localProc < mid_proc) {
    local_num_rows -= adjustment;
  }
  else {
    local_num_rows += adjustment;
  }

  //if numProcs is not an even number, undo the local_num_rows adjustment
  //on one proc so that the total will still be correct.
  if (localProc == numProcs-1) {
    if (num_procs_even == false) {
      local_num_rows -= adjustment;
    }
  }

  //now we're ready to create a row-map.
  Epetra_Map rowmap(global_num_rows, local_num_rows, 0, comm);

  //create a graph
  Epetra_CrsGraph* input_graph =
    new Epetra_CrsGraph(Copy, rowmap, nnz_per_row);

  int err = ispatest::fill_graph(*input_graph, nnz_per_row, verbose);
  if (err != 0) {
    delete input_graph;
    return(0);
  }

  return(input_graph);
}

#endif //HAVE_EPETRA

