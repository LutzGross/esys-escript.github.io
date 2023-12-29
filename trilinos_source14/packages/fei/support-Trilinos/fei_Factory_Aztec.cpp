/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/


#include "fei_trilinos_macros.hpp"

#include <fei_Factory_Aztec.hpp>

#include <fei_VectorReducer.hpp>
#include <fei_MatrixReducer.hpp>

Factory_Aztec::Factory_Aztec(MPI_Comm comm)
  : fei::Factory(comm),
    comm_(comm),
    reducer_(),
    blockEntryMatrix_(false),
    outputLevel_(0)
{
}

Factory_Aztec::~Factory_Aztec()
{
}

int Factory_Aztec::parameters(int numParams,
                                  const char* const* paramStrings)
{
  std::vector<std::string> stdstrings;
  fei::utils::char_ptrs_to_strings(numParams, paramStrings, stdstrings);

  fei::ParameterSet paramset;
  fei::utils::parse_strings(stdstrings, " ", paramset);

  parameters(paramset);
  return(0);
}

void Factory_Aztec::parameters(const fei::ParameterSet& parameterset)
{
  fei::Factory::parameters(parameterset);

  parameterset.getIntParamValue("outputLevel", outputLevel_);

  bool blkGraph = false;
  bool blkMatrix = false;

  parameterset.getBoolParamValue("BLOCK_GRAPH", blkGraph);
  parameterset.getBoolParamValue("BLOCK_MATRIX", blkMatrix);

  blockEntryMatrix_ = (blkGraph || blkMatrix);
}

fei::SharedPtr<fei::MatrixGraph>
Factory_Aztec::createMatrixGraph(fei::SharedPtr<fei::VectorSpace> rowSpace,
                      fei::SharedPtr<fei::VectorSpace> colSpace,
                      const char* name)
{
  static fei::MatrixGraph_Impl2::Factory factory2;
  return factory2.createMatrixGraph(rowSpace, colSpace, name);
}

fei::SharedPtr<fei::Vector>
Factory_Aztec::createVector(fei::SharedPtr<fei::VectorSpace> vecSpace,
                               bool isSolutionVector,
                               int numVectors)
{
  std::vector<int> indices;
  int err = 0, localSize = 0;
  if (reducer_.get() != NULL) {
    indices = reducer_->getLocalReducedEqns();
    localSize = indices.size();
  }
  else {
    if (blockEntryMatrix_) {
      localSize = vecSpace->getNumBlkIndices_Owned();
      indices.resize(localSize*2);
      err = vecSpace->getBlkIndices_Owned(localSize, &indices[0], &indices[localSize], localSize);
    }
    else {
      localSize = vecSpace->getNumIndices_Owned();
      err = vecSpace->getIndices_Owned(indices);
    }
  }
  if (err != 0) {
    throw std::runtime_error("fei::Factory_Aztec: error in vecSpace->getIndices_Owned");
  }

  fei::SharedPtr<fei::Vector> feivec, tmpvec;

  if (reducer_.get() != NULL) {
    feivec.reset(new fei::VectorReducer(reducer_,
                                        tmpvec, isSolutionVector));
  }
  else {
    feivec = tmpvec;
  }

  return(feivec);
}

fei::SharedPtr<fei::Vector>
Factory_Aztec::createVector(fei::SharedPtr<fei::VectorSpace> vecSpace,
                               int numVectors)
{
  bool isSolnVector = false;
  return(createVector(vecSpace, isSolnVector, numVectors));
}

fei::SharedPtr<fei::Vector>
Factory_Aztec::createVector(fei::SharedPtr<fei::MatrixGraph> matrixGraph,
                              int numVectors)
{
  bool isSolnVector = false;
  return(createVector(matrixGraph, isSolnVector, numVectors));
}

fei::SharedPtr<fei::Vector>
Factory_Aztec::createVector(fei::SharedPtr<fei::MatrixGraph> matrixGraph,
                               bool isSolutionVector,
                               int numVectors)
{
  int globalNumSlaves = matrixGraph->getGlobalNumSlaveConstraints();

  if (globalNumSlaves > 0 && reducer_.get()==NULL) {
    reducer_ = matrixGraph->getReducer();
  }

  fei::SharedPtr<fei::Vector> feivec, tmpvec;

  std::vector<int> indices;
  int err = 0, localSize;
  fei::SharedPtr<fei::VectorSpace> vecSpace = matrixGraph->getRowSpace();
  if (reducer_.get() != NULL) {
    indices = reducer_->getLocalReducedEqns();
    localSize = indices.size();
  }
  else {
    localSize = vecSpace->getNumIndices_Owned();
    indices.resize(localSize);
    err = vecSpace->getIndices_Owned(indices);
  }
  if (err != 0) {
    throw std::runtime_error("error in vecSpace->getIndices_Owned");
  }

  if (reducer_.get() != NULL) {
    feivec.reset(new fei::VectorReducer(reducer_, tmpvec, isSolutionVector));
  }
  else {
    feivec = tmpvec;
  }

  return(feivec);
}

fei::SharedPtr<fei::Matrix>
Factory_Aztec::createMatrix(fei::SharedPtr<fei::MatrixGraph> matrixGraph)
{
  fei::SharedPtr<fei::Matrix> feimat;
  int globalNumSlaves = matrixGraph->getGlobalNumSlaveConstraints();

  if (globalNumSlaves > 0 && reducer_.get()==NULL) {
    reducer_ = matrixGraph->getReducer();
  }

  return feimat;
}

fei::SharedPtr<fei::Solver>
Factory_Aztec::createSolver(const char* name)
{
  fei::SharedPtr<fei::Solver> solver;
  return(solver);
}

