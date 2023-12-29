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

// We don't want warnings for using deprecated create_balanced_copy.
// TODO: remove this line when removing support of create_balanced_copy.
#ifdef __deprecated
#undef __deprecated
#endif
#define __deprecated

#include <Isorropia_Exception.hpp>
#include <Isorropia_Utils.hpp>
#include <Isorropia_Epetra.hpp>
#include <Isorropia_EpetraCostDescriber.hpp>
#include <Isorropia_EpetraPartitioner.hpp>
#include <Isorropia_EpetraRedistributor.hpp>

#ifdef HAVE_EPETRA
#include <Epetra_Map.h>
#include <Epetra_Import.h>
#include <Epetra_Vector.h>
#include <Epetra_Comm.h>
#include <Epetra_IntVector.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_RowMatrix.h>
#include <Epetra_LinearProblem.h>
#endif

#if defined(HAVE_EPETRA) && defined(HAVE_MPI)
#include <Epetra_MpiComm.h>
#endif

#ifdef HAVE_MPI
#include <mpi.h>
#endif

namespace Isorropia {

namespace Epetra {

#ifdef HAVE_EPETRA

Epetra_MultiVector* create_row_weights_nnz(const Epetra_RowMatrix& input_matrix)
{
  int stride;
  const Epetra_BlockMap& input_rowmap = input_matrix.RowMatrixRowMap();
  Epetra_MultiVector* weights = new Epetra_MultiVector(input_rowmap, 1);
  double* weights_ptr = 0;
  weights->ExtractView(&weights_ptr, &stride);
  int local_num_rows = input_rowmap.NumMyElements();

  for(int i=0; i<local_num_rows; ++i) {
    int nnz;
    int err = input_matrix.NumMyRowEntries(i, nnz);
    if (err != 0) {
      throw Isorropia::Exception("create_row_weights_nnz: error in input_matrix.NumMyRowEntries");
    }

    weights_ptr[i] = 1.0*nnz;
  }

  return( weights );
}

Epetra_MultiVector* create_row_weights_nnz(const Epetra_CrsGraph& input_graph)
{
  int stride;
  const Epetra_BlockMap& input_rowmap = input_graph.RowMap();
  Epetra_MultiVector* weights = new Epetra_MultiVector(input_rowmap, 1);
  double* weights_ptr = 0;
  weights->ExtractView(&weights_ptr, &stride);
  int local_num_rows = input_rowmap.NumMyElements();

  for(int i=0; i<local_num_rows; ++i) {
    int nnz = input_graph.NumMyIndices(i);

    weights_ptr[i] = 1.0*nnz;
  }

  return( weights );
}

Epetra_MultiVector* create_unit_weights(const Epetra_MultiVector& input_coords)
{
  int stride;
  double* weights_ptr = 0;

  Epetra_MultiVector* weights = new Epetra_MultiVector(input_coords.Map(), 1);
  weights->ExtractView(&weights_ptr, &stride);

  for(int i=0; i<input_coords.MyLength(); ++i) {
    weights_ptr[i] = 1.0;
  }

  return( weights );
}

double compute_imbalance(int nprocs, std::vector<int> &offsets, double *wgts, double target)
{
  double imbalance = 1.0;

  for (int p=0; p < nprocs; p++){

    double pweight = 0.0;

    for (int row=offsets[p]; row < offsets[p+1]; row++){
      pweight += wgts[row];
    }

    double ib = 1.0;

    if (pweight <= target)
      ib += ((target - pweight) / target);
    else
      ib += ((pweight - target) / target);

    if (ib > imbalance) imbalance = ib;
  }

  return imbalance;
}
int
repartition(const Epetra_BlockMap& input_map,
	    const Epetra_MultiVector& weights,
	    std::vector<int>& newPartition,    // new partition for each of my objects
	    int& exportsSize,                  // how many of my objects are exports
	    std::vector<int>& imports)         // list of gids I will import
{
  if (!input_map.PointSameAs(weights.Map())) {
    std::string str1("Epetra::repartition ERROR, input_map not ");
    std::string str2("equivalent size/layout to weights.Map()");
    throw Isorropia::Exception(str1+str2);
  }
  int stride=0;
  double *vals = NULL;

  const Epetra_Comm& comm = input_map.Comm();
  int base = input_map.IndexBase();
  int myPID = comm.MyPID();
  int numProcs = comm.NumProc();

  int globalSize = input_map.NumGlobalElements();
  int localSize = input_map.NumMyElements();
  int mySize = ((myPID==0) ? globalSize : 0);
  int myCountSize = ((myPID==0) ? numProcs: 0);

  // Create some maps.  The weight vector is not necessarily contiguous by
  // process rank, so we need to reorganize it.

  // Distributed map of objects that are contiguous by process rank, and those
  //  objects gathered onto process zero

  Epetra_BlockMap distMap(globalSize, localSize, 1, base, comm);
  Epetra_BlockMap procZeroMap(globalSize, mySize, 1, base, comm);

  // Map of counts for each process, and counts gathered on process 0

  Epetra_BlockMap distCountMap(numProcs, 1, 1, base, comm);
  Epetra_BlockMap procZeroCountMap(numProcs, myCountSize, 1, base, comm);

  // Some importers

  Epetra_Import DistToProcZero(procZeroMap, distMap);
  Epetra_Import ProcZeroToDist(distMap, procZeroMap);
  Epetra_Import CountsToProcZero(procZeroCountMap, distCountMap);
  Epetra_Import CountsToProcs(distCountMap, procZeroCountMap);

  // Move all weights to process 0

  weights.ExtractView(&vals, &stride);
  Epetra_Vector distWeights(Copy, distMap, vals);
  Epetra_Vector weightVal(procZeroMap);
  weightVal.Import(distWeights, DistToProcZero, Insert);

  // Move the count of weights for each process to process 0

  Epetra_Vector distCount(distCountMap);
  distCount[0] = localSize;
  Epetra_Vector weightCount(procZeroCountMap);
  weightCount.Import(distCount, CountsToProcZero, Insert);

  // Move all global IDs to process 0

  double *myGIDs = new double [localSize];
  for (int i=0; i<localSize; i++){
    myGIDs[i] = input_map.GID(i);
  }
  Epetra_Vector localGIDs(View, distMap, myGIDs);
  Epetra_Vector gids(procZeroMap);
  gids.Import(localGIDs, DistToProcZero, Insert);
  delete [] myGIDs;

  // Process zero computes a new balance

  double *numImports = new double [numProcs];
  int totalImports=0;
  double *importGIDs = NULL;
  double *ids = NULL;

  if (myPID == 0){

    double total_weight=0.0; 
    weightVal.ExtractView(&vals);
    for (int i=0; i<globalSize; i++){
      total_weight += vals[i];
    }

    double avg_weight = total_weight/numProcs;

    std::vector<int> all_proc_old_offsets(numProcs+1, globalSize);
    std::vector<int> all_proc_new_offsets(numProcs+1, globalSize);

    for (int i=numProcs-1; i >= 0 ; i--){
      all_proc_old_offsets[i] = all_proc_old_offsets[i+1] - static_cast<int>(weightCount[i]);
    }

    double old_imbalance =
      compute_imbalance(numProcs, all_proc_old_offsets, vals, avg_weight);

    int offset = 0;
    for(int p=0; p<numProcs; ++p) {
      all_proc_new_offsets[p] = offset;

      double tmp_weight = 0.0;

      while(offset < globalSize && tmp_weight < avg_weight) {
        tmp_weight += vals[offset++];
      }
    }

    double new_imbalance =
      compute_imbalance(numProcs, all_proc_new_offsets, vals, avg_weight);

    // Because this is a quick and dirty partitioning, it is possible that
    // if the balance was good to begin with, that we have just created a
    // slightly worse balance.  In that case, return to the old partitioning.

    if (new_imbalance >= old_imbalance){
      for (int proc=0; proc <= numProcs; proc++){
        all_proc_new_offsets[proc] = all_proc_old_offsets[proc];
      }
    }

    // Fill weight vector with the new process ID (partition ID) for each object

    for (int proc=0; proc < numProcs; proc++){
      int len = all_proc_new_offsets[proc+1] - all_proc_new_offsets[proc];
      for (int j=0; j<len; j++){
        *vals++ = static_cast<double>(proc);
      }
    }

    // Count the number of imports for each process, and create a list of
    // gids to be imported to each process

    totalImports = 0;
    importGIDs = new double [globalSize];
    double *curr = importGIDs;
    gids.ExtractView(&ids);

    for (int proc=0; proc < numProcs; proc++){
      numImports[proc] = 0;
      int from = all_proc_new_offsets[proc];
      int to =  all_proc_new_offsets[proc+1];
      int old_from = all_proc_old_offsets[proc];
      int old_to =  all_proc_old_offsets[proc+1];

      for (int i=from; i<to ; i++){
        if ((i< old_from) || (i >= old_to)){
          *curr++ = ids[i];
          numImports[proc]++;
        }
      }
      totalImports += static_cast<int>(numImports[proc]);
    }
  }

  // The new partition IDs are in the weightVal vector on process zero

  distWeights.Import(weightVal, ProcZeroToDist, Insert);

  newPartition.resize(localSize);
  for (int i=0; i<localSize; i++){
    newPartition[i] = static_cast<int>(distWeights[i]);
  }

  // Get the count of my imports - reuse the weightCount vector

  int *indices = new int [numProcs];
  for (int i=0; i<numProcs; i++){
    indices[i] = i + base;
  }

  weightCount.ReplaceGlobalValues(((myPID==0) ? numProcs : 0), numImports, indices);
  distCount.Import(weightCount, CountsToProcs, Insert);
  delete [] indices;
  if (numImports) delete [] numImports;

  // Get the list of my imports

  int numMyImports = static_cast<int>(distCount[0]);
  double n1;
  distCount.Norm1(&n1);
  totalImports = static_cast<int>(n1);
  int myImportSize = ((myPID==0) ? totalImports : 0);

  Epetra_BlockMap distImportMap(totalImports, numMyImports, 1, base, comm);
  Epetra_BlockMap proc0ImportMap(totalImports, myImportSize, 1, base, comm);
  Epetra_Import importer(distImportMap, proc0ImportMap);

  Epetra_Vector proc0imports(View, proc0ImportMap, importGIDs);
  Epetra_Vector distImports(distImportMap);

  distImports.Import(proc0imports, importer, Insert);

  if (importGIDs) delete [] importGIDs;

  imports.resize(numMyImports);
  for (int i=0; i<numMyImports; i++){
    imports[i] = static_cast<int>(distImports[i]);
  }

  // Finally, the number of exports

  exportsSize = 0;
  for (int i=0; i<localSize; i++){
    if (newPartition[i] != myPID){
      exportsSize++;
    }
  }

  return 0;
}

/* New createBalancedCopy functions */

Epetra_MultiVector * 
createBalancedCopy(const Epetra_MultiVector &coords)
{
  Teuchos::ParameterList paramlist;

  return createBalancedCopy(coords, paramlist);
}

Epetra_MultiVector * 
createBalancedCopy(const Epetra_MultiVector &coords,
		   const Teuchos::ParameterList& paramlist)
{
  Teuchos::RCP<const Epetra_MultiVector> coordRcp = Teuchos::rcp(&coords, false);

  Teuchos::RCP<Partitioner> partitioner =
    Teuchos::rcp(new Partitioner(coordRcp, paramlist));

  Redistributor rd(partitioner);

  Teuchos::RCP<Epetra_MultiVector> newVec = rd.redistribute(coords);

  newVec.release();

  return newVec.get();
}

Epetra_CrsGraph *
createBalancedCopy(const Epetra_CrsGraph& input_graph)
{
  Teuchos::ParameterList paramlist;

  return createBalancedCopy(input_graph, paramlist);
}

Epetra_CrsGraph *
createBalancedCopy(const Epetra_CrsGraph& input_graph,
		     const Teuchos::ParameterList& paramlist)
{
  Teuchos::RCP<const Epetra_CrsGraph> rcp_input_graph =
    Teuchos::rcp(&(input_graph), false);

  Teuchos::RCP<Partitioner> partitioner =
    Teuchos::rcp(new Partitioner(rcp_input_graph, paramlist));

  Redistributor rd(partitioner);

  Teuchos::RCP<Epetra_CrsGraph> balanced_graph =
    rd.redistribute(input_graph);

  balanced_graph.release();

  return balanced_graph.get();
}

Epetra_CrsMatrix *
createBalancedCopy(const Epetra_CrsMatrix& input_matrix)
{
  Teuchos::ParameterList paramlist;

  return createBalancedCopy(input_matrix, paramlist);
}

Epetra_CrsMatrix *
createBalancedCopy(const Epetra_CrsMatrix& input_matrix,
		     const Teuchos::ParameterList& paramlist)
{
  Teuchos::RCP<const Epetra_CrsGraph> input_graph =
    Teuchos::rcp(&(input_matrix.Graph()), false);

  Teuchos::RCP<Partitioner> partitioner =
    Teuchos::rcp(new Partitioner(input_graph, paramlist));

  Redistributor rd(partitioner);

  Teuchos::RCP<Epetra_CrsMatrix> balanced_matrix =
    rd.redistribute(input_matrix);

  balanced_matrix.release();
  return balanced_matrix.get();
}

Epetra_LinearProblem *
createBalancedCopy(const Epetra_LinearProblem& input_problem)
{
  Teuchos::ParameterList paramlist;

  return createBalancedCopy(input_problem, paramlist);
}

Epetra_LinearProblem *
createBalancedCopy(const Epetra_LinearProblem& input_problem,
		     const Teuchos::ParameterList& paramlist)
{
  Teuchos::RCP<const Epetra_RowMatrix> rowmat =
    Teuchos::rcp(input_problem.GetMatrix(), false);

  Teuchos::RCP<Partitioner> partitioner =
    Teuchos::rcp(new Partitioner(rowmat, paramlist));

  Redistributor rd(partitioner);

  Teuchos::RCP<Epetra_RowMatrix> balanced_matrix =
    rd.redistribute(*input_problem.GetMatrix());

  Teuchos::RCP<Epetra_MultiVector> balanced_rhs =
    rd.redistribute(*input_problem.GetRHS());

  Teuchos::RCP<Epetra_MultiVector> x=
    Teuchos::rcp(new Epetra_MultiVector(*input_problem.GetLHS()));

  // prevent these from being deallocated on return from createBalancedCopy
  balanced_matrix.release(); 
  balanced_rhs.release();
  x.release();

  Teuchos::RCP<Epetra_LinearProblem> linprob =
    Teuchos::rcp(new Epetra_LinearProblem(balanced_matrix.get(), x.get(), balanced_rhs.get()));

  linprob.release();

  return linprob.get();
}

/*
************* Deprecated methods below - remove later ********************
*/

Teuchos::RCP<Epetra_RowMatrix>
create_balanced_copy(const Epetra_RowMatrix& input_matrix)
{
  CostDescriber costs;
  Teuchos::ParameterList paramlist;

  Teuchos::RCP<Epetra_RowMatrix> balanced_matrix =
    create_balanced_copy(input_matrix, costs, paramlist);

  return balanced_matrix;
}

Teuchos::RCP<Epetra_RowMatrix>
create_balanced_copy(const Epetra_RowMatrix& input_matrix,
                     const Epetra_Vector &row_weights)
{
  CostDescriber costs; 
  Teuchos::ParameterList paramlist;

  Teuchos::RCP<const Epetra_Vector> vwgts = Teuchos::rcp(&row_weights);
  vwgts.release();
  costs.setVertexWeights(vwgts);

  Teuchos::RCP<Epetra_RowMatrix> balanced_matrix =
    create_balanced_copy(input_matrix, costs, paramlist);

  return balanced_matrix;
}

Teuchos::RCP<Epetra_RowMatrix>
create_balanced_copy(const Epetra_RowMatrix& input_matrix,
		     const Teuchos::ParameterList& paramlist)
{
  CostDescriber costs; 

  Teuchos::RCP<Epetra_RowMatrix> balanced_matrix =
    create_balanced_copy(input_matrix, costs, paramlist);

  return balanced_matrix;
}

Teuchos::RCP<Epetra_RowMatrix>
create_balanced_copy(const Epetra_RowMatrix& input_matrix,
                     CostDescriber &costs,
		     const Teuchos::ParameterList& paramlist)
{
  Teuchos::RCP<const Epetra_RowMatrix> matrixPtr=
    Teuchos::rcp(&(input_matrix), false);

  Teuchos::RCP<CostDescriber> costPtr =
    Teuchos::rcp(&(costs), false);

  Teuchos::RCP<Partitioner> partitioner =
    Teuchos::rcp(new Partitioner(matrixPtr, costPtr, paramlist));

  Redistributor rd(partitioner);

  Teuchos::RCP<Epetra_RowMatrix> balanced_matrix =
    rd.redistribute(input_matrix);

  return balanced_matrix;
}

Teuchos::RCP<Epetra_CrsMatrix>
create_balanced_copy(const Epetra_CrsMatrix& input_matrix)
{
  CostDescriber costs; 
  Teuchos::ParameterList paramlist;

  Teuchos::RCP<Epetra_CrsMatrix> balanced_matrix =
    create_balanced_copy(input_matrix, costs, paramlist);

  return balanced_matrix;
}

Teuchos::RCP<Epetra_CrsMatrix>
create_balanced_copy(const Epetra_CrsMatrix& input_matrix,
                     const Epetra_Vector &row_weights)
{
  CostDescriber costs; 
  Teuchos::ParameterList paramlist;

  Teuchos::RCP<const Epetra_Vector> vwgts = Teuchos::rcp(&row_weights);
  // release of ownership of vwgts, so destruction of stack objects
  // does not do an extra destruction of row_weights and cause
  // a fatal fault
  vwgts.release();  

  costs.setVertexWeights(vwgts);

  Teuchos::RCP<Epetra_CrsMatrix> balanced_matrix =
    create_balanced_copy(input_matrix, costs, paramlist);

  return balanced_matrix;
}

Teuchos::RCP<Epetra_CrsMatrix>
create_balanced_copy(const Epetra_CrsMatrix& input_matrix,
		     const Teuchos::ParameterList& paramlist)
{
  CostDescriber costs; 

  Teuchos::RCP<Epetra_CrsMatrix> balanced_matrix =
    create_balanced_copy(input_matrix, costs, paramlist);

  return balanced_matrix;
}

Teuchos::RCP<Epetra_CrsMatrix>
create_balanced_copy(const Epetra_CrsMatrix& input_matrix,
                     CostDescriber &costs,
		     const Teuchos::ParameterList& paramlist)
{
  Teuchos::RCP<const Epetra_CrsGraph> input_graph =
    Teuchos::rcp(&(input_matrix.Graph()), false);

  Teuchos::RCP<CostDescriber> costPtr =
    Teuchos::rcp(&(costs), false);

  Teuchos::RCP<Partitioner> partitioner =
    Teuchos::rcp(new Partitioner(input_graph, costPtr, paramlist));

  Redistributor rd(partitioner);

  Teuchos::RCP<Epetra_CrsMatrix> balanced_matrix =
    rd.redistribute(input_matrix);

  return balanced_matrix;
}

Teuchos::RCP<Epetra_CrsGraph>
create_balanced_copy(const Epetra_CrsGraph& input_graph)
{
  CostDescriber costs; 
  Teuchos::ParameterList paramlist;

  Teuchos::RCP<Epetra_CrsGraph> balanced_graph =
    create_balanced_copy(input_graph, costs, paramlist);

  return balanced_graph;
}

Teuchos::RCP<Epetra_CrsGraph>
create_balanced_copy(const Epetra_CrsGraph& input_graph,
                     const Epetra_Vector &row_weights)
{
  CostDescriber costs; 
  Teuchos::ParameterList paramlist;

  Teuchos::RCP<const Epetra_Vector> vwgts = Teuchos::rcp(&row_weights);
  // release of ownership of vwgts, so destruction of stack objects
  // does not do an extra destruction of row_weights and cause
  // a fatal fault
  vwgts.release();  

  costs.setVertexWeights(vwgts);

  Teuchos::RCP<Epetra_CrsGraph> balanced_graph =
    create_balanced_copy(input_graph, costs, paramlist);

  return balanced_graph;
}

Teuchos::RCP<Epetra_CrsGraph>
create_balanced_copy(const Epetra_CrsGraph& input_graph,
		     const Teuchos::ParameterList& paramlist)
{
  CostDescriber costs; 

  Teuchos::RCP<Epetra_CrsGraph> balanced_graph =
    create_balanced_copy(input_graph, costs, paramlist);

  return balanced_graph;
}

Teuchos::RCP<Epetra_CrsGraph>
create_balanced_copy(const Epetra_CrsGraph& input_graph,
                     CostDescriber &costs,
		     const Teuchos::ParameterList& paramlist)
{
  Teuchos::RCP<const Epetra_CrsGraph> graphPtr=
    Teuchos::rcp(&(input_graph), false);

  Teuchos::RCP<CostDescriber> costPtr =
    Teuchos::rcp(&(costs), false);

  Teuchos::RCP<Partitioner> partitioner =
    Teuchos::rcp(new Partitioner(graphPtr, costPtr, paramlist));

  Redistributor rd(partitioner);

  Teuchos::RCP<Epetra_CrsGraph> balanced_graph =
    rd.redistribute(input_graph);

  return balanced_graph;
}

Teuchos::RCP<Epetra_LinearProblem>
create_balanced_copy(const Epetra_LinearProblem& input_problem)
{
  CostDescriber costs; 
  Teuchos::ParameterList paramlist;

  Teuchos::RCP<Epetra_LinearProblem> linprob =
    create_balanced_copy(input_problem, costs, paramlist);

  return linprob;
}

Teuchos::RCP<Epetra_LinearProblem>
create_balanced_copy(const Epetra_LinearProblem& input_problem,
                     const Epetra_Vector &row_weights)
{
  CostDescriber costs; 
  Teuchos::ParameterList paramlist;

  Teuchos::RCP<const Epetra_Vector> vwgts = Teuchos::rcp(&row_weights);
  // release of ownership of vwgts, so destruction of stack objects
  // does not do an extra destruction of row_weights and cause
  // a fatal fault
  vwgts.release();  
  costs.setVertexWeights(vwgts);

  Teuchos::RCP<Epetra_LinearProblem> linprob =
    create_balanced_copy(input_problem, costs, paramlist);

  return linprob;
}

Teuchos::RCP<Epetra_LinearProblem>
create_balanced_copy(const Epetra_LinearProblem& input_problem,
		     const Teuchos::ParameterList& paramlist)
{
  CostDescriber costs; 

  Teuchos::RCP<Epetra_LinearProblem> linprob =
    create_balanced_copy(input_problem, costs, paramlist);

  return linprob;
}

Teuchos::RCP<Epetra_LinearProblem>
create_balanced_copy(const Epetra_LinearProblem& input_problem,
                     CostDescriber &costs,
		     const Teuchos::ParameterList& paramlist)
{
  Teuchos::RCP<const Epetra_RowMatrix> rowmat =
    Teuchos::rcp(input_problem.GetMatrix(), false);

  Teuchos::RCP<CostDescriber> costPtr =
    Teuchos::rcp(&(costs), false);

  Teuchos::RCP<Partitioner> partitioner =
    Teuchos::rcp(new Partitioner(rowmat, costPtr, paramlist));

  Redistributor rd(partitioner);

  Teuchos::RCP<Epetra_RowMatrix> balanced_matrix =
    rd.redistribute(*input_problem.GetMatrix());

  Teuchos::RCP<Epetra_MultiVector> balanced_rhs =
    rd.redistribute(*input_problem.GetRHS());

  Teuchos::RCP<Epetra_MultiVector> x=
    Teuchos::rcp(new Epetra_MultiVector(*input_problem.GetLHS()));

  // prevent these from being deallocated on return from create_balanced_copy
  balanced_matrix.release(); 
  balanced_rhs.release();
  x.release();

  Teuchos::RCP<Epetra_LinearProblem> linprob =
    Teuchos::rcp(new Epetra_LinearProblem(balanced_matrix.get(), x.get(), balanced_rhs.get()));

  return( linprob );
}

Teuchos::RCP<Epetra_MultiVector>
create_balanced_copy(const Epetra_MultiVector &coords,
		   const Teuchos::ParameterList& paramlist)
{
  // TODO make sure weights list with 0 vectors doesn't cause something
  // to crash
  Epetra_MultiVector noWeights(coords.Map(), 0);
  return create_balanced_copy(coords, noWeights, paramlist);
}

Teuchos::RCP<Epetra_MultiVector>
create_balanced_copy(const Epetra_MultiVector &coords,
                     const Epetra_MultiVector &weights,
		   const Teuchos::ParameterList& paramlist)
{
  Teuchos::RCP<const Epetra_MultiVector> coordRcp = Teuchos::rcp(&coords, false);
  Teuchos::RCP<const Epetra_MultiVector> weightRcp = Teuchos::rcp(&weights, false);

  Teuchos::RCP<Partitioner> partitioner =
    Teuchos::rcp(new Partitioner(coordRcp, weightRcp, paramlist));

  Redistributor rd(partitioner);

  Teuchos::RCP<Epetra_MultiVector> newVec = rd.redistribute(coords);

  return newVec;
}

#endif
}//namespace Epetra
}//namespace Isorropia
