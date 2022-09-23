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

// Assumes we have Epetra and Zoltan support
#include <Isorropia_EpetraPartitioner2D.hpp>
#include <Isorropia_Zoltan_Repartition.hpp>
#include <Isorropia_EpetraZoltanLib.hpp>
#include <Isorropia_Exception.hpp>
#include <Isorropia_Epetra.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_Import.h>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>

#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
#include <ctype.h>

namespace Isorropia {

namespace Epetra {


  /** Constructor that accepts an Epetra_CrsGraph object, called by
        API function create_partitioner().

     \param input_graph Matrix-graph object for which a new partitioning
        is to be computed. A Teuchos::RCP is used here because a
        reference to the input object may be held by this object after
        this constructor completes and returns.

     \param paramlist Teuchos::ParameterList which will be copied to an
        internal ParameterList attribute. No reference to this input
        object is held after this constructor completes.<br>
  If the ParameterList object contains a sublist named "Zoltan", then
  the Zoltan library is used to perform the balancing. Also, any
  parameters in the "Zoltan" sublist will be relayed directly to Zoltan.
  Refer to the Zoltan users guide for specific parameters that Zoltan
  recognizes. A couple of important ones are "LB_METHOD" (valid values
  include "GRAPH", "HYPERGRAPH"), "DEBUG_LEVEL" (valid values are
  0 to 10, default is 1), etc.

     \param compute_partitioning_now Optional argument defaults to true.
        If true, the method compute_partitioning() will be called before
        this constructor returns.
  */
Partitioner2D::Partitioner2D(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (input_graph, paramlist, 0)
{
  if (compute_partitioning_now)
    partition(true);
}

  /**
     Constructor that accepts an Epetra_RowMatrix object, called by
       API function create_partitioner().

     \param input_matrix Matrix object for which a new partitioning is
        to be computed. A Teuchos::RCP is used here because a
        reference to the input object may be held by this object after
        this constructor completes and returns.

     \param paramlist Teuchos::ParameterList which will be copied to an
        internal ParameterList attribute. No reference to this input
        object is held after this constructor completes.<br>
  If the ParameterList object contains a sublist named "Zoltan", then
  the Zoltan library is used to perform the balancing. Also, any
  parameters in the "Zoltan" sublist will be relayed directly to Zoltan.
  Refer to the Zoltan users guide for specific parameters that Zoltan
  recognizes. A couple of important ones are "LB_METHOD" (valid values
  include "GRAPH", "HYPERGRAPH"), "DEBUG_LEVEL" (valid values are
  0 to 10, default is 1), etc.

     \param compute_partitioning_now Optional argument defaults to true.
        If true, the method compute_partitioning() will be called before
        this constructor returns.
  */
Partitioner2D::Partitioner2D(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
			 const Teuchos::ParameterList& paramlist,
			 bool compute_partitioning_now):
  Operator (input_matrix, paramlist, 0)
{
  if (compute_partitioning_now)
    partition(true);
}


  /** Destructor */
Partitioner2D::~Partitioner2D(){}

////////////////////////////////////////////////////////////////////////////////
  /** setParameters() is an internal Partitioner2D method which handles
      the parameters from a Teuchos::ParameterList object. 

      The input
      ParameterList object is copied into an internal ParameterList
      attribute, and no reference to the input object is held after
      this function returns. (Thus, the input paramlist object may be
      altered or destroyed as soon as this method returns.)<br>
  If the ParameterList object contains a sublist named "Zoltan", then
  the Zoltan library is used to perform the balancing. Also, any
  parameters in the "Zoltan" sublist will be relayed directly to Zoltan.
  Refer to the Zoltan users guide for specific parameters that Zoltan
  recognizes. A couple of important ones are "LB_METHOD" (valid values
  include "GRAPH", "HYPERGRAPH"), "DEBUG_LEVEL" (valid values are
  0 to 10, default is 1), etc.
   */

  /**  partition is a method that computes 
       a rebalanced partitioning for the data in the object
      that this class was constructed with.

      \param force_repartitioning Optional argument defaults to false. By
         default, compute_partitioning() only does anything the first time
         it is called, and subsequent repeated calls are no-ops. If the user's
         intent is to re-compute the partitioning (e.g., if parameters
         or other inputs have been changed), then setting this flag to
         true will force a new partitioning to be computed.
   */
////////////////////////////////////////////////////////////////////////////////
void Partitioner2D::
partition(bool force_repartitioning)
{
  std::string partitioning_method_str("PARTITIONING METHOD");
  std::string partitioning_method =
    paramlist_.get(partitioning_method_str, "UNSPECIFIED");

  if(partitioning_method == "UNSPECIFIED")
  {
    throw Isorropia::Exception("PARTITIONING_METHOD parameter must be specified.");
  }

  if(partitioning_method != "HYPERGRAPH2D")
  {
    throw Isorropia::Exception("PARTITIONING_METHOD parameter must be HYPERGRAPH2D.");
  }

  if (alreadyComputed() && !force_repartitioning)
    return;

  // Determine whether graph input or matrix input is used
  if (input_graph_.get() != 0)
    lib_ = Teuchos::rcp(new ZoltanLibClass(input_graph_, costs_, 
             Library::hgraph2d_finegrain_input_));
  else
    lib_ = Teuchos::rcp(new ZoltanLibClass(input_matrix_, costs_, 
             Library::hgraph2d_finegrain_input_));

  std::string zoltan("ZOLTAN");
  Teuchos::ParameterList &sublist = paramlist_.sublist(zoltan);

  if (paramlist_.isParameter("NUM PARTS")) 
  {
    sublist.set("NUM_GLOBAL_PARTS", paramlist_.get<std::string>("NUM PARTS"));
  }
  if (paramlist_.isParameter("IMBALANCE TOL")) 
  {
    sublist.set("IMBALANCE_TOL", paramlist_.get<std::string>("IMBALANCE TOL"));
  }

  lib_->repartition(sublist, properties_, exportsSize_, imports_);
  computeNumberOfProperties();
  operation_already_computed_ = true;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
void Partitioner2D::
compute(bool force_repartitioning)
{
  partition(force_repartitioning);
}
////////////////////////////////////////////////////////////////////////////////

//  /** An internal method which determines whether the 
//      method compute_partitioning() has already been
//      called on this class instance.
//  */
//bool Partitioner2D::partitioning_already_computed() const 
//{
//  return (alreadyComputed());
//}

//  /** An internal method which returns the new partition ID for a given element that
//     resided locally in the old partitioning.
//  */
//int Partitioner2D::newPartitionNumber(int myElem) const
//{
//  return ((*this)[myElem]);
//}

  /** An internal method which returns the number of elements in a given partition.

      (Currently only implemented for the case where 'partition' is local.)
  */
int Partitioner2D::numElemsInPart(int part) const
{
  return (numElemsWithProperty(part));
}

////////////////////////////////////////////////////////////////////////////////
  /** An internal method which fills caller-allocated list (of length len) with the
      global element ids to be located in the given partition.

      (Currently only implemented for the case where 'partition' is local.)
  */
////////////////////////////////////////////////////////////////////////////////
void Partitioner2D::elemsInPart(int partition, int* elementList, int len) const 
{
  //MMW
  std::cout << "MMW::NEED to reimplement" << std::endl;

  return (elemsWithProperty(partition, elementList, len));
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
  /** Returns indx for nonzero A(i,j) assuming this nonzero is located on this
      processor.  This is useful in getting a part number for this nonzero since
      the [] operator needs this index.  Returns -1 if nonzero not found locally.
  */
////////////////////////////////////////////////////////////////////////////////
int Partitioner2D::getNZIndx(int row, int column) const 
{


  if (input_graph_.get() != 0) // Graph
  {
    const Epetra_BlockMap *rowMap_ = &(input_graph_->RowMap());
    const Epetra_BlockMap *colMap_ = &(input_graph_->ColMap());
    int numRows = rowMap_->NumMyElements();

    int nzIndx = 0;

    for (int rowNum=0; rowNum<numRows; rowNum++)
    {
      int rowSize;
      rowSize = input_graph_->NumMyIndices(rowNum);

      if(row == rowNum)
      {
        int *tmprowCols = new int[rowSize];
        int numEntries;

        input_graph_->ExtractMyRowCopy (rowNum, rowSize, numEntries, tmprowCols);

        for(int colIndx=0;colIndx<rowSize;colIndx++)
	{

          if (column == colMap_->GID(tmprowCols[colIndx]))
	  {
            return nzIndx;
	  }
          else
	  {
	    nzIndx++;
          }
	}

        delete [] tmprowCols;
      }
      else
      {
        nzIndx += rowSize;
      }


    }



  }
  else  // matrix
  {

    const Epetra_BlockMap *rowMap_ = &(input_matrix_->RowMatrixRowMap());
    const Epetra_BlockMap *colMap_ = &(input_matrix_->RowMatrixColMap());
    int numRows = rowMap_->NumMyElements();

    int nzIndx = 0;
    for (int rowNum=0; rowNum<numRows; rowNum++)
    {
      int rowSize;
      input_matrix_->NumMyRowEntries(rowNum,rowSize);

      if(row == rowNum)
      {
        int *tmprowCols = new int[rowSize];
        double *tmprowVals = new double[rowSize];
        int numEntries;

        input_matrix_->ExtractMyRowCopy (rowNum, rowSize, numEntries, tmprowVals, tmprowCols);

        for(int colIndx=0; colIndx<rowSize; colIndx++)
	{
          if (column == colMap_->GID(tmprowCols[colIndx]))
	  {
            return nzIndx;
	  }
          else
	  {
	    nzIndx++;
          }
	}
	delete [] tmprowCols;
        delete [] tmprowVals;
      }
      else
      {
        nzIndx += rowSize;
      }
    }
  }

  return -1;
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
int Partitioner2D::createDomainAndRangeMaps(Epetra_Map *domainMap, 
			                    Epetra_Map *rangeMap) 
{
  std::cout << "MMW::NEED to implement" << std::endl;

  return 0;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
int Partitioner2D::createColumnMap(Epetra_Map *colMap) 
{
  std::cout << "MMW::NEED to implement" << std::endl;

  return 0;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
int Partitioner2D::createRowMap(Epetra_Map *rowMap) 
{
  std::cout << "MMW::NEED to implement" << std::endl;

  return 0;
}
////////////////////////////////////////////////////////////////////////////////

} // namespace EPETRA

}//namespace Isorropia

