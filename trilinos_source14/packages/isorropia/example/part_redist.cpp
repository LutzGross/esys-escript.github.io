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

//--------------------------------------------------------------------
//This file is a self-contained example of creating an Epetra_LinearProblem
//object, and using Isorropia to create a rebalanced copy of it.
//--------------------------------------------------------------------

//Include Isorropia_Exception.hpp only because the helper functions at
//the bottom of this file (which create the epetra objects) can
//potentially throw exceptions.
#include <Isorropia_Exception.hpp>

//The Isorropia symbols being demonstrated are declared
//in these headers:
#include <Isorropia_Epetra.hpp>
#include <Isorropia_EpetraRedistributor.hpp>
#include <Isorropia_EpetraPartitioner.hpp>
#include <Isorropia_EpetraCostDescriber.hpp>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifdef HAVE_EPETRA
#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>
#endif

#include "ispatest_lbeval_utils.hpp"

//Declaration for helper-function that creates epetra objects. This
//function is implemented at the bottom of this file.
#ifdef HAVE_EPETRA
Epetra_LinearProblem* create_epetra_problem(int numProcs,
                                            int localProc,
                                            int local_n);
#endif

int main(int argc, char** argv) {
#if defined(HAVE_MPI) && defined(HAVE_EPETRA)

  int numProcs = 1;
  int localProc = 0;

  //first, set up our MPI environment...
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &localProc);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  int local_n = 600;

  //Create a Epetra_LinearProblem object.

  Epetra_LinearProblem* linprob = 0;
  try {
    linprob = create_epetra_problem(numProcs, localProc, local_n);
  }
  catch(std::exception& exc) {
    std::cout << "linsys example: create_epetra_problem threw exception '"
          << exc.what() << "' on proc " << localProc << std::endl;
    MPI_Finalize();
    return(-1);
  }

  //We'll need a Teuchos::ParameterList object to pass to the
  //Isorropia::Epetra::Partitioner class.
  Teuchos::ParameterList paramlist;

  // If Zoltan is available, the Zoltan package will be used for
  // the partitioning operation. By default, Isorropia selects Zoltan's
  // Hypergraph partitioner. If a method other than Hypergraph is
  // desired, it can be specified by first creating a parameter sublist
  // named "Zoltan", and then setting appropriate Zoltan parameters in
  // that sublist. A sublist is created like this:
  //      Teuchos::ParameterList& sublist = paramlist.sublist("Zoltan");
  //

  // If Zoltan is not available, a simple linear partitioner will be
  // used to partition such that the number of nonzeros is equal (or
  // close to equal) on each processor.


  Epetra_RowMatrix* rowmatrix = linprob->GetMatrix();
  Teuchos::RCP<const Epetra_RowMatrix> rowmat =
    Teuchos::rcp(rowmatrix, false);


  //Now create the partitioner 

  Teuchos::RCP<Isorropia::Epetra::Partitioner> partitioner =
    Teuchos::rcp(new Isorropia::Epetra::Partitioner(rowmat, paramlist));

  //Next create a Redistributor object and use it to create balanced
  //copies of the objects in linprob.

  Isorropia::Epetra::Redistributor rd(partitioner);

  Teuchos::RCP<Epetra_CrsMatrix> bal_matrix;
  Teuchos::RCP<Epetra_MultiVector> bal_x;
  Teuchos::RCP<Epetra_MultiVector> bal_b;

  //Use a try-catch block because Isorropia will throw an exception
  //if it encounters an error.

  if (localProc == 0) {
    std::cout << " calling Isorropia::Epetra::Redistributor::redistribute..."
        << std::endl;
  }

  try {
    bal_matrix = rd.redistribute(*linprob->GetMatrix());
    bal_x = rd.redistribute(*linprob->GetLHS());
    bal_b = rd.redistribute(*linprob->GetRHS());
  }
  catch(std::exception& exc) {
    std::cout << "linsys example: Isorropia::Epetra::Redistributor threw "
         << "exception '" << exc.what() << "' on proc "
         << localProc << std::endl;
    MPI_Finalize();
    return(-1);
  }

  Epetra_LinearProblem balanced_problem(bal_matrix.get(),
                                        bal_x.get(), bal_b.get());


  // Results

  double bal0, bal1, cutn0, cutn1, cutl0, cutl1;
  Isorropia::Epetra::CostDescriber default_costs;

#if 1
  // Balance and cut quality before partitioning

  double goalWeight = 1.0 / (double)numProcs;

  ispatest::compute_hypergraph_metrics(*(linprob->GetMatrix()), default_costs, goalWeight,
                     bal0, cutn0, cutl0);

  // Balance and cut quality after partitioning

  ispatest::compute_hypergraph_metrics(*bal_matrix, default_costs, goalWeight,
                     bal1, cutn1, cutl1);
#else

  std::vector<double> bal(2), cutn(2), cutl(2);

  Epetra_Import &importer = rd.get_importer();

  default_costs.compareBeforeAndAfterHypergraph(*(linprob->GetMatrix()), *bal_matrix, importer,
             bal, cutn, cutl);

  bal0 = bal[0]; cutn0 = cutn[0]; cutl0 = cutl[0];
  bal1 = bal[1]; cutn1 = cutn[1]; cutl1 = cutl[1];
#endif


  if (localProc == 0){
    std::cout << "Before partitioning: ";
    std::cout << "Balance " << bal0 << " cutN " << cutn0 << " cutL " << cutl0;
    std::cout << std::endl;

    std::cout << "After partitioning:  ";
    std::cout << "Balance " << bal1 << " cutN " << cutn1 << " cutL " << cutl1;
    std::cout << std::endl;
  }

  //Finally, delete the pointer objects that we asked to be created.
  delete linprob->GetMatrix();
  delete linprob->GetLHS();
  delete linprob->GetRHS();
  delete linprob;

  if (localProc == 0) {
    std::cout << std::endl;
  }

  MPI_Finalize();

#else
  std::cout << "part_redist: must have both MPI and EPETRA. Make sure Trilinos "
    << "is configured with --enable-mpi and --enable-epetra." << std::endl;
#endif

  return(0);
}

//Below is the implementation of the helper-function that creates the
//poorly-balanced epetra objects for use in the above example program.

#if defined(HAVE_MPI) && defined(HAVE_EPETRA)

Epetra_LinearProblem* create_epetra_problem(int numProcs,
                                            int localProc,
                                            int local_n)
{
  if (localProc == 0) {
    std::cout << " creating Epetra_CrsMatrix with un-even distribution..."
            << std::endl;
  }

  //create an Epetra_CrsMatrix with rows spread un-evenly over
  //processors.
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int global_num_rows = numProcs*local_n;

  int mid_proc = numProcs/2;
  bool num_procs_even = numProcs%2==0 ? true : false;

  int adjustment = local_n/2;

  //adjust local_n so that it's not equal on all procs.
  if (localProc < mid_proc) {
    local_n -= adjustment;
  }
  else {
    local_n += adjustment;
  }

  //if numProcs is not an even number, undo the local_n adjustment
  //on one proc so that the total will still be correct.
  if (localProc == numProcs-1) {
    if (num_procs_even == false) {
      local_n -= adjustment;
    }
  }

  //now we're ready to create a row-map.
  Epetra_Map rowmap(global_num_rows, local_n, 0, comm);

  //create a matrix
  int nnz_per_row = 9;
  Epetra_CrsMatrix* matrix =
    new Epetra_CrsMatrix(Copy, rowmap, nnz_per_row);

  // Add  rows one-at-a-time
  double negOne = -1.0;
  double posTwo = 4.0;
  for (int i=0; i<local_n; i++) {
    int GlobalRow = matrix->GRID(i);
    int RowLess1 = GlobalRow - 1;
    int RowPlus1 = GlobalRow + 1;
    int RowLess2 = GlobalRow - 2;
    int RowPlus2 = GlobalRow + 2;
    int RowLess3 = GlobalRow - 3;
    int RowPlus3 = GlobalRow + 3;
    int RowLess4 = GlobalRow - 4;
    int RowPlus4 = GlobalRow + 4;

    if (RowLess4>=0) {
      matrix->InsertGlobalValues(GlobalRow, 1, &negOne, &RowLess4);
    }
    if (RowLess3>=0) {
      matrix->InsertGlobalValues(GlobalRow, 1, &negOne, &RowLess3);
    }
    if (RowLess2>=0) {
      matrix->InsertGlobalValues(GlobalRow, 1, &negOne, &RowLess2);
    }
    if (RowLess1>=0) {
      matrix->InsertGlobalValues(GlobalRow, 1, &negOne, &RowLess1);
    }
    if (RowPlus1<global_num_rows) {
      matrix->InsertGlobalValues(GlobalRow, 1, &negOne, &RowPlus1);
    }
    if (RowPlus2<global_num_rows) {
      matrix->InsertGlobalValues(GlobalRow, 1, &negOne, &RowPlus2);
    }
    if (RowPlus3<global_num_rows) {
      matrix->InsertGlobalValues(GlobalRow, 1, &negOne, &RowPlus3);
    }
    if (RowPlus4<global_num_rows) {
      matrix->InsertGlobalValues(GlobalRow, 1, &negOne, &RowPlus4);
    }

    matrix->InsertGlobalValues(GlobalRow, 1, &posTwo, &GlobalRow);
  }

  int err = matrix->FillComplete();
  if (err != 0) {
    throw Isorropia::Exception("create_epetra_matrix: error in matrix.FillComplete()");
  }

  Epetra_Vector* x = new Epetra_Vector(rowmap);
  Epetra_Vector* b = new Epetra_Vector(rowmap);
  return(new Epetra_LinearProblem(matrix, x, b));
}

#endif //HAVE_MPI && HAVE_EPETRA

