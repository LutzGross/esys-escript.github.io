//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER

// Transform_Composite Test routine

#include <EpetraExt_ConfigDefs.h>
#include "EpetraExt_Version.h"

#include "Epetra_SerialComm.h"
#include "Epetra_Time.h"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#endif

#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"

#include "EpetraExt_Transform_Composite.h"
#include "EpetraExt_LPTrans_From_GraphTrans.h"
#include "EpetraExt_SymmRCM_CrsGraph.h"
#include "EpetraExt_Reindex_LinearProblem.h"
//#include "EpetraExt_CrsSingletonFilter_LinearProblem.h"

#ifdef EPETRA_MPI
#include "EpetraExt_Overlap_CrsGraph.h"
#endif

#include "../epetra_test_err.h"

int main(int argc, char *argv[]) {

  int i, ierr=0, returnierr=0;

#ifdef EPETRA_MPI

  // Initialize MPI

  MPI_Init(&argc,&argv);
  int size, rank; // Number of MPI processes, My process ID

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

#else

  int size = 1; // Serial case (not using MPI)
  int rank = 0;

#endif

  bool verbose = false;

  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;


#ifdef EPETRA_MPI
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif
  if (!verbose) Comm.SetTracebackMode(0); // This should shut down any error traceback reporting

  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  bool verbose1 = verbose;
  if( verbose ) verbose = (MyPID==0);

  if ( verbose )
    cout << EpetraExt::EpetraExt_Version() << endl << endl;

  Comm.Barrier();

  if( verbose1 ) cout << Comm << endl << flush;

  int NumMyElements = 3;
  int NumGlobalElements = NumProc*NumMyElements;
  int IndexBase = 0;
  
  Epetra_Map Map( NumGlobalElements, NumMyElements, 0, Comm );
  if( verbose1 ) cout << Map << endl << flush;

  Epetra_CrsGraph Graph( Copy, Map, 1 );

  int PIDFac = 10*MyPID;
  int index = PIDFac+2;
  Graph.InsertGlobalIndices( PIDFac+0, 1, &index );
  index = PIDFac+0;
  Graph.InsertGlobalIndices( PIDFac+1, 1, &index );
  index = PIDFac+1;
  Graph.InsertGlobalIndices( PIDFac+2, 1, &index );

  Graph.FillComplete();
  if( verbose1 ) cout << Graph << endl << flush;

  EpetraExt::Transform_Composite<Epetra_LinearProblem> CompTrans;

//  EpetraExt::LinearProblem_CrsSingletonFilter CSF_LPTrans;
//  EpetraExt::SameTypeTransform<Epetra_LinearProblem> * CSF_LPTransPtr = &CSF_LPTrans;
//  CompTrans.addTransform( CSF_LPTransPtr );

  EpetraExt::LinearProblem_Reindex * RI_Trans = new EpetraExt::LinearProblem_Reindex(0);
  EpetraExt::SameTypeTransform<Epetra_LinearProblem> * RI_LPTrans = RI_Trans;
  CompTrans.addTransform( RI_LPTrans );

  EpetraExt::CrsGraph_SymmRCM RCM_Trans;
  EpetraExt::SameTypeTransform<Epetra_LinearProblem> *
        RCM_LPTrans = new EpetraExt::LinearProblem_GraphTrans( RCM_Trans );
  CompTrans.addTransform( RCM_LPTrans );

#ifdef EPETRA_MPI
  EpetraExt::CrsGraph_Overlap Overlap_Trans(1);
  EpetraExt::SameTypeTransform<Epetra_LinearProblem> *
        Overlap_LPTrans = new EpetraExt::LinearProblem_GraphTrans( Overlap_Trans );
  CompTrans.addTransform( Overlap_LPTrans );
#endif

  Epetra_CrsMatrix Matrix( Copy, Graph );
  index = 2;
  double val = 2;
  Matrix.ReplaceMyValues( 0, 1, &val, &index );
  index = 0;
  val = 0;
  Matrix.ReplaceMyValues( 1, 1, &val, &index );
  index = 1;
  val = 1;
  Matrix.ReplaceMyValues( 2, 1, &val, &index);

  vector<double> valA(3);
  valA[0]=0; valA[1]=1; valA[2]=2;
  Epetra_BlockMap & MapRef = Map;
  Epetra_Vector LHS( Copy, MapRef, &valA[0] );
  Epetra_Vector RHS( Copy, MapRef, &valA[0] );

  Epetra_LinearProblem Prob( &Matrix, &LHS, &RHS );

  Epetra_LinearProblem & NewProb = CompTrans( Prob );

  CompTrans.fwd();
  CompTrans.rvs();

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return ierr;
}

