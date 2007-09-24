
/* $Id$ */

/*******************************************************
 *
 *           Copyright 2003-2007 by ACceSS MNRF
 *       Copyright 2007 by University of Queensland
 *
 *                http://esscc.uq.edu.au
 *        Primary Business: Queensland, Australia
 *  Licensed under the Open Software License version 3.0
 *     http://www.opensource.org/licenses/osl-3.0.php
 *
 *******************************************************/

#include "MeshAdapter.h"
#include "escript/Data.h"
#include "escript/DataFactory.h"
extern "C" {
#include "escript/blocktimer.h"
}
#include <vector>

using namespace std;
using namespace escript;

namespace finley {

//
// define the static constants
MeshAdapter::FunctionSpaceNamesMapType MeshAdapter::m_functionSpaceTypeNames;
const int MeshAdapter::DegreesOfFreedom=FINLEY_DEGREES_OF_FREEDOM;
const int MeshAdapter::ReducedDegreesOfFreedom=FINLEY_REDUCED_DEGREES_OF_FREEDOM;
const int MeshAdapter::Nodes=FINLEY_NODES;
const int MeshAdapter::ReducedNodes=FINLEY_REDUCED_NODES;
const int MeshAdapter::Elements=FINLEY_ELEMENTS;
const int MeshAdapter::ReducedElements=FINLEY_REDUCED_ELEMENTS;
const int MeshAdapter::FaceElements=FINLEY_FACE_ELEMENTS;
const int MeshAdapter::ReducedFaceElements=FINLEY_REDUCED_FACE_ELEMENTS;
const int MeshAdapter::Points=FINLEY_POINTS;
const int MeshAdapter::ContactElementsZero=FINLEY_CONTACT_ELEMENTS_1;
const int MeshAdapter::ReducedContactElementsZero=FINLEY_REDUCED_CONTACT_ELEMENTS_1;
const int MeshAdapter::ContactElementsOne=FINLEY_CONTACT_ELEMENTS_2;
const int MeshAdapter::ReducedContactElementsOne=FINLEY_REDUCED_CONTACT_ELEMENTS_2;

MeshAdapter::MeshAdapter(Finley_Mesh* finleyMesh)
{
  setFunctionSpaceTypeNames();
  //
  // need to use a null_deleter as Finley_Mesh_free deletes the pointer
  // for us.
  m_finleyMesh.reset(finleyMesh,null_deleter());
}

//
// The copy constructor should just increment the use count
MeshAdapter::MeshAdapter(const MeshAdapter& in):
m_finleyMesh(in.m_finleyMesh)
{
  setFunctionSpaceTypeNames();
}

MeshAdapter::~MeshAdapter()
{
  //
  // I hope the case for the pointer being zero has been taken care of.
  //  cout << "In MeshAdapter destructor." << endl;
  if (m_finleyMesh.unique()) {
    Finley_Mesh_free(m_finleyMesh.get());
  }
}

int MeshAdapter::getMPISize() const
{
   return m_finleyMesh.get()->MPIInfo->size;
}
int MeshAdapter::getMPIRank() const
{
   return m_finleyMesh.get()->MPIInfo->rank;
}


Finley_Mesh* MeshAdapter::getFinley_Mesh() const {
   return m_finleyMesh.get();
}

void MeshAdapter::write(const std::string& fileName) const
{
  char *fName = (fileName.size()+1>0) ? TMPMEMALLOC(fileName.size()+1,char) : (char*)NULL;
  strcpy(fName,fileName.c_str());
  Finley_Mesh_write(m_finleyMesh.get(),fName);
  checkFinleyError();
  TMPMEMFREE(fName);
}

void MeshAdapter::dump(const std::string& fileName) const
{
  char *fName = (fileName.size()+1>0) ? TMPMEMALLOC(fileName.size()+1,char) : (char*)NULL;
  strcpy(fName,fileName.c_str());
  Finley_Mesh_dump(m_finleyMesh.get(),fName);
  checkFinleyError();
  TMPMEMFREE(fName);
}

string MeshAdapter::getDescription() const
{
  return "FinleyMesh";
}

string MeshAdapter::functionSpaceTypeAsString(int functionSpaceType) const
{
  FunctionSpaceNamesMapType::iterator loc;
  loc=m_functionSpaceTypeNames.find(functionSpaceType);
  if (loc==m_functionSpaceTypeNames.end()) {
    return "Invalid function space type code.";
  } else {
    return loc->second;
  }
}

bool MeshAdapter::isValidFunctionSpaceType(int functionSpaceType) const
{
  FunctionSpaceNamesMapType::iterator loc;
  loc=m_functionSpaceTypeNames.find(functionSpaceType);
  return (loc!=m_functionSpaceTypeNames.end());
}

void MeshAdapter::setFunctionSpaceTypeNames()
{
  m_functionSpaceTypeNames.insert
    (FunctionSpaceNamesMapType::value_type(DegreesOfFreedom,"Finley_DegreesOfFreedom"));
  m_functionSpaceTypeNames.insert
    (FunctionSpaceNamesMapType::value_type(ReducedDegreesOfFreedom,"Finley_ReducedDegreesOfFreedom"));
  m_functionSpaceTypeNames.insert
    (FunctionSpaceNamesMapType::value_type(Nodes,"Finley_Nodes"));
  m_functionSpaceTypeNames.insert
    (FunctionSpaceNamesMapType::value_type(ReducedNodes,"Finley_Reduced_Nodes"));
  m_functionSpaceTypeNames.insert
    (FunctionSpaceNamesMapType::value_type(Elements,"Finley_Elements"));
  m_functionSpaceTypeNames.insert
    (FunctionSpaceNamesMapType::value_type(ReducedElements,"Finley_Reduced_Elements"));
  m_functionSpaceTypeNames.insert
    (FunctionSpaceNamesMapType::value_type(FaceElements,"Finley_Face_Elements"));
  m_functionSpaceTypeNames.insert
    (FunctionSpaceNamesMapType::value_type(ReducedFaceElements,"Finley_Reduced_Face_Elements"));
  m_functionSpaceTypeNames.insert
    (FunctionSpaceNamesMapType::value_type(Points,"Finley_Points"));
  m_functionSpaceTypeNames.insert
    (FunctionSpaceNamesMapType::value_type(ContactElementsZero,"Finley_Contact_Elements_0"));
  m_functionSpaceTypeNames.insert
    (FunctionSpaceNamesMapType::value_type(ReducedContactElementsZero,"Finley_Reduced_Contact_Elements_0"));
  m_functionSpaceTypeNames.insert
    (FunctionSpaceNamesMapType::value_type(ContactElementsOne,"Finley_Contact_Elements_1"));
  m_functionSpaceTypeNames.insert
    (FunctionSpaceNamesMapType::value_type(ReducedContactElementsOne,"Finley_Reduced_Contact_Elements_1"));
}

int MeshAdapter::getContinuousFunctionCode() const
{
  return Nodes;
}
int MeshAdapter::getReducedContinuousFunctionCode() const
{
  return ReducedNodes;
}

int MeshAdapter::getFunctionCode() const
{
  return Elements;
}
int MeshAdapter::getReducedFunctionCode() const
{
  return ReducedElements;
}

int MeshAdapter::getFunctionOnBoundaryCode() const
{
  return FaceElements;
}
int MeshAdapter::getReducedFunctionOnBoundaryCode() const
{
  return ReducedFaceElements;
}

int MeshAdapter::getFunctionOnContactZeroCode() const
{
  return ContactElementsZero;
}
int MeshAdapter::getReducedFunctionOnContactZeroCode() const
{
  return ReducedContactElementsZero;
}

int MeshAdapter::getFunctionOnContactOneCode() const
{
  return ContactElementsOne;
}
int MeshAdapter::getReducedFunctionOnContactOneCode() const
{
  return ReducedContactElementsOne;
}

int MeshAdapter::getSolutionCode() const
{
  return DegreesOfFreedom;
}

int MeshAdapter::getReducedSolutionCode() const
{
  return ReducedDegreesOfFreedom;
}

int MeshAdapter::getDiracDeltaFunctionCode() const
{
  return Points;
}

//
// return the spatial dimension of the Mesh:
//
int MeshAdapter::getDim() const
{
  int numDim=Finley_Mesh_getDim(m_finleyMesh.get());
  checkFinleyError();
  return numDim;
}

//
// return the number of data points per sample and the number of samples
// needed to represent data on a parts of the mesh.
//
pair<int,int> MeshAdapter::getDataShape(int functionSpaceCode) const
{
   int numDataPointsPerSample=0;
   int numSamples=0;
   Finley_Mesh* mesh=m_finleyMesh.get();
   switch (functionSpaceCode) {
      case(Nodes):
           numDataPointsPerSample=1;
           numSamples=Finley_NodeFile_getNumNodes(mesh->Nodes);
           break;
      case(ReducedNodes):
           numDataPointsPerSample=1;
           numSamples=Finley_NodeFile_getNumReducedNodes(mesh->Nodes);
           break;
      case(Elements):
           if (mesh->Elements!=NULL) {
             numSamples=mesh->Elements->numElements;
             numDataPointsPerSample=mesh->Elements->ReferenceElement->numQuadNodes;
           }
           break;
      case(ReducedElements):
           if (mesh->Elements!=NULL) {
             numSamples=mesh->Elements->numElements;
             numDataPointsPerSample=mesh->Elements->ReferenceElementReducedOrder->numQuadNodes;
           }
           break;
      case(FaceElements):
           if (mesh->FaceElements!=NULL) {
                numDataPointsPerSample=mesh->FaceElements->ReferenceElement->numQuadNodes;
                numSamples=mesh->FaceElements->numElements;
           }
           break;
      case(ReducedFaceElements):
           if (mesh->FaceElements!=NULL) {
                numDataPointsPerSample=mesh->FaceElements->ReferenceElementReducedOrder->numQuadNodes;
                numSamples=mesh->FaceElements->numElements;
           }
           break;
      case(Points):
           if (mesh->Points!=NULL) {
             numDataPointsPerSample=1;
             numSamples=mesh->Points->numElements;
           }
           break;
      case(ContactElementsZero):
           if (mesh->ContactElements!=NULL) {
             numDataPointsPerSample=mesh->ContactElements->ReferenceElement->numQuadNodes;
             numSamples=mesh->ContactElements->numElements;
           }
           break;
      case(ReducedContactElementsZero):
           if (mesh->ContactElements!=NULL) {
             numDataPointsPerSample=mesh->ContactElements->ReferenceElementReducedOrder->numQuadNodes;
             numSamples=mesh->ContactElements->numElements;
           }
           break;
      case(ContactElementsOne):
           if (mesh->ContactElements!=NULL) {
             numDataPointsPerSample=mesh->ContactElements->ReferenceElement->numQuadNodes;
             numSamples=mesh->ContactElements->numElements;
           }
           break;
      case(ReducedContactElementsOne):
           if (mesh->ContactElements!=NULL) {
             numDataPointsPerSample=mesh->ContactElements->ReferenceElementReducedOrder->numQuadNodes;
             numSamples=mesh->ContactElements->numElements;
           }
           break;
      case(DegreesOfFreedom):
           if (mesh->Nodes!=NULL) {
             numDataPointsPerSample=1;
             numSamples=Finley_NodeFile_getNumDegreesOfFreedom(mesh->Nodes);
           }
           break;
      case(ReducedDegreesOfFreedom):
           if (mesh->Nodes!=NULL) {
             numDataPointsPerSample=1;
             numSamples=Finley_NodeFile_getNumReducedDegreesOfFreedom(mesh->Nodes);
           }
           break;
      default:
        stringstream temp;
        temp << "Error - Invalid function space type: " << functionSpaceCode << " for domain: " << getDescription();
        throw FinleyAdapterException(temp.str());
        break;
      }
      return pair<int,int>(numDataPointsPerSample,numSamples);
}

//
// adds linear PDE of second order into a given stiffness matrix and right hand side:
//
void MeshAdapter::addPDEToSystem(
                     SystemMatrixAdapter& mat, escript::Data& rhs,
                     const escript::Data& A, const escript::Data& B, const escript::Data& C,const  escript::Data& D,const  escript::Data& X,const  escript::Data& Y,
                     const escript::Data& d, const escript::Data& y, 
                     const escript::Data& d_contact,const escript::Data& y_contact) const
{
   escriptDataC _rhs=rhs.getDataC();
   escriptDataC _A  =A.getDataC();
   escriptDataC _B=B.getDataC();
   escriptDataC _C=C.getDataC();
   escriptDataC _D=D.getDataC();
   escriptDataC _X=X.getDataC();
   escriptDataC _Y=Y.getDataC();
   escriptDataC _d=d.getDataC();
   escriptDataC _y=y.getDataC();
   escriptDataC _d_contact=d_contact.getDataC();
   escriptDataC _y_contact=y_contact.getDataC();

   Finley_Mesh* mesh=m_finleyMesh.get();

   Finley_Assemble_PDE(mesh->Nodes,mesh->Elements,mat.getPaso_SystemMatrix(), &_rhs, &_A, &_B, &_C, &_D, &_X, &_Y );
   checkFinleyError();

   Finley_Assemble_PDE(mesh->Nodes,mesh->FaceElements, mat.getPaso_SystemMatrix(), &_rhs, 0, 0, 0, &_d, 0, &_y );
   checkFinleyError();

   Finley_Assemble_PDE(mesh->Nodes,mesh->ContactElements, mat.getPaso_SystemMatrix(), &_rhs , 0, 0, 0, &_d_contact, 0, &_y_contact );
   checkFinleyError();
}

void  MeshAdapter::addPDEToLumpedSystem(
                     escript::Data& mat,
                     const escript::Data& D,
                     const escript::Data& d) const
{
   escriptDataC _mat=mat.getDataC();
   escriptDataC _D=D.getDataC();
   escriptDataC _d=d.getDataC();

   Finley_Mesh* mesh=m_finleyMesh.get();

   Finley_Assemble_LumpedSystem(mesh->Nodes,mesh->Elements,&_mat, &_D);
   Finley_Assemble_LumpedSystem(mesh->Nodes,mesh->FaceElements,&_mat, &_d);

   checkFinleyError();
}


//
// adds linear PDE of second order into the right hand side only
//
void MeshAdapter::addPDEToRHS( escript::Data& rhs, const  escript::Data& X,const  escript::Data& Y, const escript::Data& y, const escript::Data& y_contact) const
{
   Finley_Mesh* mesh=m_finleyMesh.get();

   escriptDataC _rhs=rhs.getDataC();
   escriptDataC _X=X.getDataC();
   escriptDataC _Y=Y.getDataC();
   escriptDataC _y=y.getDataC();
   escriptDataC _y_contact=y_contact.getDataC();

   Finley_Assemble_PDE(mesh->Nodes,mesh->Elements, 0, &_rhs, 0, 0, 0, 0, &_X, &_Y );
   checkFinleyError();

   Finley_Assemble_PDE(mesh->Nodes,mesh->FaceElements, 0, &_rhs, 0, 0, 0, 0, 0, &_y );
   checkFinleyError();

   Finley_Assemble_PDE(mesh->Nodes,mesh->ContactElements, 0, &_rhs , 0, 0, 0, 0, 0, &_y_contact );
   checkFinleyError();
}

//
// interpolates data between different function spaces:
//
void MeshAdapter::interpolateOnDomain(escript::Data& target,const escript::Data& in) const
{
  const MeshAdapter& inDomain=dynamic_cast<const MeshAdapter&>(in.getFunctionSpace().getDomain());
  const MeshAdapter& targetDomain=dynamic_cast<const MeshAdapter&>(target.getFunctionSpace().getDomain());
  if (inDomain!=*this)  
    throw FinleyAdapterException("Error - Illegal domain of interpolant.");
  if (targetDomain!=*this) 
    throw FinleyAdapterException("Error - Illegal domain of interpolation target.");

  Finley_Mesh* mesh=m_finleyMesh.get();
  escriptDataC _target=target.getDataC();
  escriptDataC _in=in.getDataC();
  switch(in.getFunctionSpace().getTypeCode()) {
     case(Nodes):
        switch(target.getFunctionSpace().getTypeCode()) {
           case(Nodes):
           case(ReducedNodes):
           case(DegreesOfFreedom):
           case(ReducedDegreesOfFreedom):
               Finley_Assemble_CopyNodalData(mesh->Nodes,&_target,&_in);
               break;
           case(Elements):
           case(ReducedElements):
               Finley_Assemble_interpolate(mesh->Nodes,mesh->Elements,&_in,&_target);
               break;
           case(FaceElements):
           case(ReducedFaceElements):
               Finley_Assemble_interpolate(mesh->Nodes,mesh->FaceElements,&_in,&_target);
               break;
           case(Points):
               Finley_Assemble_interpolate(mesh->Nodes,mesh->Points,&_in,&_target);
               break;
           case(ContactElementsZero):
           case(ReducedContactElementsZero):
               Finley_Assemble_interpolate(mesh->Nodes,mesh->ContactElements,&_in,&_target);
               break;
           case(ContactElementsOne):
           case(ReducedContactElementsOne):
               Finley_Assemble_interpolate(mesh->Nodes,mesh->ContactElements,&_in,&_target);
               break;
           default:
               stringstream temp;
               temp << "Error - Interpolation on Domain: Finley does not know anything about function space type " << target.getFunctionSpace().getTypeCode();
               throw FinleyAdapterException(temp.str());
               break;
        }
        break;
     case(ReducedNodes):
        switch(target.getFunctionSpace().getTypeCode()) {
           case(Nodes):
           case(ReducedNodes):
           case(DegreesOfFreedom):
           case(ReducedDegreesOfFreedom):
               Finley_Assemble_CopyNodalData(mesh->Nodes,&_target,&_in);
               break;
           case(Elements):
           case(ReducedElements):
               Finley_Assemble_interpolate(mesh->Nodes,mesh->Elements,&_in,&_target);
               break;
           case(FaceElements):
           case(ReducedFaceElements):
               Finley_Assemble_interpolate(mesh->Nodes,mesh->FaceElements,&_in,&_target);
               break;
           case(Points):
               Finley_Assemble_interpolate(mesh->Nodes,mesh->Points,&_in,&_target);
               break;
           case(ContactElementsZero):
           case(ReducedContactElementsZero):
               Finley_Assemble_interpolate(mesh->Nodes,mesh->ContactElements,&_in,&_target);
               break;
           case(ContactElementsOne):
           case(ReducedContactElementsOne):
               Finley_Assemble_interpolate(mesh->Nodes,mesh->ContactElements,&_in,&_target);
               break;
           default:
               stringstream temp;
               temp << "Error - Interpolation on Domain: Finley does not know anything about function space type " << target.getFunctionSpace().getTypeCode();
               throw FinleyAdapterException(temp.str());
               break;
        }
        break;
     case(Elements):
        if (target.getFunctionSpace().getTypeCode()==Elements) {
           Finley_Assemble_CopyElementData(mesh->Elements,&_target,&_in);
        } else if (target.getFunctionSpace().getTypeCode()==ReducedElements) {
           Finley_Assemble_AverageElementData(mesh->Elements,&_target,&_in);
        } else {
           throw FinleyAdapterException("Error - No interpolation with data on elements possible.");
        }
        break;
     case(ReducedElements):
        if (target.getFunctionSpace().getTypeCode()==ReducedElements) {
           Finley_Assemble_CopyElementData(mesh->Elements,&_target,&_in);
        } else {
           throw FinleyAdapterException("Error - No interpolation with data on elements with reduced integration order possible.");
        }
        break;
     case(FaceElements):
        if (target.getFunctionSpace().getTypeCode()==FaceElements) {
           Finley_Assemble_CopyElementData(mesh->FaceElements,&_target,&_in);
        } else if (target.getFunctionSpace().getTypeCode()==ReducedFaceElements) {
           Finley_Assemble_AverageElementData(mesh->FaceElements,&_target,&_in);
        } else {
           throw FinleyAdapterException("Error - No interpolation with data on face elements possible.");
       }
       break;
     case(ReducedFaceElements):
        if (target.getFunctionSpace().getTypeCode()==FaceElements) {
           Finley_Assemble_CopyElementData(mesh->FaceElements,&_target,&_in);
        } else {
           throw FinleyAdapterException("Error - No interpolation with data on face elements with reduced integration order possible.");
       }
       break;
     case(Points):
        if (target.getFunctionSpace().getTypeCode()==Points) {
           Finley_Assemble_CopyElementData(mesh->Points,&_target,&_in);
        } else {
           throw FinleyAdapterException("Error - No interpolation with data on points possible.");
        }
        break;
     case(ContactElementsZero):
     case(ContactElementsOne):
        if (target.getFunctionSpace().getTypeCode()==ContactElementsZero || target.getFunctionSpace().getTypeCode()==ContactElementsOne) {
           Finley_Assemble_CopyElementData(mesh->ContactElements,&_target,&_in);
        } else if (target.getFunctionSpace().getTypeCode()==ReducedContactElementsZero || target.getFunctionSpace().getTypeCode()==ReducedContactElementsOne) {
           Finley_Assemble_AverageElementData(mesh->ContactElements,&_target,&_in);
        } else {
           throw FinleyAdapterException("Error - No interpolation with data on contact elements possible.");
           break;
        }
        break;
     case(ReducedContactElementsZero):
     case(ReducedContactElementsOne):
        if (target.getFunctionSpace().getTypeCode()==ReducedContactElementsZero || target.getFunctionSpace().getTypeCode()==ReducedContactElementsOne) {
           Finley_Assemble_CopyElementData(mesh->ContactElements,&_target,&_in);
        } else {
           throw FinleyAdapterException("Error - No interpolation with data on contact elements with reduced integration order possible.");
           break;
        }
        break;
     case(DegreesOfFreedom):      
        switch(target.getFunctionSpace().getTypeCode()) {
           case(ReducedDegreesOfFreedom):
           case(DegreesOfFreedom):
              Finley_Assemble_CopyNodalData(mesh->Nodes,&_target,&_in);
              break;

           case(Nodes):
           case(ReducedNodes):
              if (getMPISize()>1) {
                  escript::Data temp=escript::Data(in);
                  temp.expand();
                  escriptDataC _in2 = temp.getDataC();
                  Finley_Assemble_CopyNodalData(mesh->Nodes,&_target,&_in2);
              } else {
                  Finley_Assemble_CopyNodalData(mesh->Nodes,&_target,&_in);
              }
              break;
           case(Elements):
           case(ReducedElements):
              if (getMPISize()>1) {
                 escript::Data temp=escript::Data( in,  continuousFunction(asAbstractContinuousDomain()) );
                 escriptDataC _in2 = temp.getDataC();
                 Finley_Assemble_interpolate(mesh->Nodes,mesh->Elements,&_in2,&_target);
              } else {
                 Finley_Assemble_interpolate(mesh->Nodes,mesh->Elements,&_in,&_target);
              }
              break;
           case(FaceElements):
           case(ReducedFaceElements):
              if (getMPISize()>1) {
                 escript::Data temp=escript::Data( in,  continuousFunction(asAbstractContinuousDomain()) );
                 escriptDataC _in2 = temp.getDataC();
                 Finley_Assemble_interpolate(mesh->Nodes,mesh->FaceElements,&_in2,&_target);

              } else {
                 Finley_Assemble_interpolate(mesh->Nodes,mesh->FaceElements,&_in,&_target);
              }
              break;
           case(Points):
              if (getMPISize()>1) {
                 escript::Data temp=escript::Data( in,  continuousFunction(asAbstractContinuousDomain()) );
                 escriptDataC _in2 = temp.getDataC();
              } else {
                 Finley_Assemble_interpolate(mesh->Nodes,mesh->Points,&_in,&_target);
              }
              break;
           case(ContactElementsZero):
           case(ContactElementsOne):
           case(ReducedContactElementsZero):
           case(ReducedContactElementsOne):
              if (getMPISize()>1) {
                 escript::Data temp=escript::Data( in,  continuousFunction(asAbstractContinuousDomain()) );
                 escriptDataC _in2 = temp.getDataC();
                 Finley_Assemble_interpolate(mesh->Nodes,mesh->ContactElements,&_in2,&_target);
              } else {
                 Finley_Assemble_interpolate(mesh->Nodes,mesh->ContactElements,&_in,&_target);
              }
             break;
           default:
             stringstream temp;
             temp << "Error - Interpolation On Domain: Finley does not know anything about function space type " << target.getFunctionSpace().getTypeCode();
             throw FinleyAdapterException(temp.str());
             break;
        }
        break;
     case(ReducedDegreesOfFreedom):
       switch(target.getFunctionSpace().getTypeCode()) {
          case(Nodes):
             throw FinleyAdapterException("Error - Finley does not support interpolation from reduced degrees of freedom to mesh nodes.");
             break;
          case(ReducedNodes):
              if (getMPISize()>1) {
                  escript::Data temp=escript::Data(in);
                  temp.expand();
                  escriptDataC _in2 = temp.getDataC();
                  Finley_Assemble_CopyNodalData(mesh->Nodes,&_target,&_in2);
              } else {
                  Finley_Assemble_CopyNodalData(mesh->Nodes,&_target,&_in);
              }
              break;
          case(DegreesOfFreedom):
             throw FinleyAdapterException("Error - Finley does not support interpolation from reduced degrees of freedom to degrees of freedom");
             break;
          case(ReducedDegreesOfFreedom):
             Finley_Assemble_CopyNodalData(mesh->Nodes,&_target,&_in);
             break;
          case(Elements):
          case(ReducedElements):
              if (getMPISize()>1) {
                 escript::Data temp=escript::Data( in,  continuousFunction(asAbstractContinuousDomain()) );
                 escriptDataC _in2 = temp.getDataC();
                 Finley_Assemble_interpolate(mesh->Nodes,mesh->Elements,&_in2,&_target);
              } else {
                 Finley_Assemble_interpolate(mesh->Nodes,mesh->Elements,&_in,&_target);
             }
             break;
          case(FaceElements):
          case(ReducedFaceElements):
              if (getMPISize()>1) {
                 escript::Data temp=escript::Data( in,  continuousFunction(asAbstractContinuousDomain()) );
                 escriptDataC _in2 = temp.getDataC();
                 Finley_Assemble_interpolate(mesh->Nodes,mesh->FaceElements,&_in2,&_target);
              } else {
                 Finley_Assemble_interpolate(mesh->Nodes,mesh->FaceElements,&_in,&_target);
              }
             break;
          case(Points):
              if (getMPISize()>1) {
                 escript::Data temp=escript::Data( in,  continuousFunction(asAbstractContinuousDomain()) );
                 escriptDataC _in2 = temp.getDataC();
                 Finley_Assemble_interpolate(mesh->Nodes,mesh->Points,&_in2,&_target);
              } else {
                 Finley_Assemble_interpolate(mesh->Nodes,mesh->Points,&_in,&_target);
             }
             break;
          case(ContactElementsZero):
          case(ContactElementsOne):
          case(ReducedContactElementsZero):
          case(ReducedContactElementsOne):
              if (getMPISize()>1) {
                 escript::Data temp=escript::Data( in,  continuousFunction(asAbstractContinuousDomain()) );
                 escriptDataC _in2 = temp.getDataC();
                 Finley_Assemble_interpolate(mesh->Nodes,mesh->ContactElements,&_in2,&_target);
              } else {
                 Finley_Assemble_interpolate(mesh->Nodes,mesh->ContactElements,&_in,&_target);
             }
             break;
          default:
             stringstream temp;
             temp << "Error - Interpolation On Domain: Finley does not know anything about function space type " << target.getFunctionSpace().getTypeCode();
             throw FinleyAdapterException(temp.str());
             break;
       }
       break;
     default:
        stringstream temp;
        temp << "Error - Interpolation On Domain: Finley does not know anything about function space type %d" << in.getFunctionSpace().getTypeCode();
        throw FinleyAdapterException(temp.str());
        break;
  }
  checkFinleyError();
}

//
// copies the locations of sample points into x:
//
void MeshAdapter::setToX(escript::Data& arg) const
{
  const MeshAdapter& argDomain=dynamic_cast<const MeshAdapter&>(arg.getFunctionSpace().getDomain());
  if (argDomain!=*this) 
     throw FinleyAdapterException("Error - Illegal domain of data point locations");
  Finley_Mesh* mesh=m_finleyMesh.get();
  // in case of values node coordinates we can do the job directly:
  if (arg.getFunctionSpace().getTypeCode()==Nodes) {
     escriptDataC _arg=arg.getDataC();
     Finley_Assemble_NodeCoordinates(mesh->Nodes,&_arg);
  } else {
     escript::Data tmp_data=Vector(0.0,continuousFunction(asAbstractContinuousDomain()),true);
     escriptDataC _tmp_data=tmp_data.getDataC();
     Finley_Assemble_NodeCoordinates(mesh->Nodes,&_tmp_data);
     // this is then interpolated onto arg:
     interpolateOnDomain(arg,tmp_data);
  }
  checkFinleyError();
}

//
// return the normal vectors at the location of data points as a Data object:
//
void MeshAdapter::setToNormal(escript::Data& normal) const
{
  const MeshAdapter& normalDomain=dynamic_cast<const MeshAdapter&>(normal.getFunctionSpace().getDomain());
  if (normalDomain!=*this) 
     throw FinleyAdapterException("Error - Illegal domain of normal locations");
  Finley_Mesh* mesh=m_finleyMesh.get();
  escriptDataC _normal=normal.getDataC();
  switch(normal.getFunctionSpace().getTypeCode()) {
    case(Nodes):
      throw FinleyAdapterException("Error - Finley does not support surface normal vectors for nodes");
      break;
    case(ReducedNodes):
      throw FinleyAdapterException("Error - Finley does not support surface normal vectors for reduced nodes");
      break;
    case(Elements):
      throw FinleyAdapterException("Error - Finley does not support surface normal vectors for elements");
      break;
    case(ReducedElements):
      throw FinleyAdapterException("Error - Finley does not support surface normal vectors for elements with reduced integration order");
      break;
    case (FaceElements):
      Finley_Assemble_setNormal(mesh->Nodes,mesh->FaceElements,&_normal);
      break;
    case (ReducedFaceElements):
      Finley_Assemble_setNormal(mesh->Nodes,mesh->FaceElements,&_normal);
      break;
    case(Points):
      throw FinleyAdapterException("Error - Finley does not support surface normal vectors for point elements");
      break;
    case (ContactElementsOne):
    case (ContactElementsZero):
      Finley_Assemble_setNormal(mesh->Nodes,mesh->ContactElements,&_normal);
      break;
    case (ReducedContactElementsOne):
    case (ReducedContactElementsZero):
      Finley_Assemble_setNormal(mesh->Nodes,mesh->ContactElements,&_normal);
      break;
    case(DegreesOfFreedom):
      throw FinleyAdapterException("Error - Finley does not support surface normal vectors for degrees of freedom.");
      break;
    case(ReducedDegreesOfFreedom):
      throw FinleyAdapterException("Error - Finley does not support surface normal vectors for reduced degrees of freedom.");
      break;
    default:
      stringstream temp;
      temp << "Error - Normal Vectors: Finley does not know anything about function space type " << normal.getFunctionSpace().getTypeCode();
      throw FinleyAdapterException(temp.str());
      break;
  }
  checkFinleyError();
}

//
// interpolates data to other domain:
//
void MeshAdapter::interpolateACross(escript::Data& target,const escript::Data& source) const
{
  const MeshAdapter& targetDomain=dynamic_cast<const MeshAdapter&>(target.getFunctionSpace().getDomain());
  if (targetDomain!=*this) 
     throw FinleyAdapterException("Error - Illegal domain of interpolation target");

  throw FinleyAdapterException("Error - Finley does not allow interpolation across domains yet.");
}

//
// calculates the integral of a function defined of arg:
//
void MeshAdapter::setToIntegrals(std::vector<double>& integrals,const escript::Data& arg) const
{
  const MeshAdapter& argDomain=dynamic_cast<const MeshAdapter&>(arg.getFunctionSpace().getDomain());
  if (argDomain!=*this) 
     throw FinleyAdapterException("Error - Illegal domain of integration kernel");

  double blocktimer_start = blocktimer_time();
  Finley_Mesh* mesh=m_finleyMesh.get();
  escriptDataC _temp;
  escript::Data temp;
  escriptDataC _arg=arg.getDataC();
  switch(arg.getFunctionSpace().getTypeCode()) {
     case(Nodes):
        temp=escript::Data( arg, function(asAbstractContinuousDomain()) );
        _temp=temp.getDataC();
        Finley_Assemble_integrate(mesh->Nodes,mesh->Elements,&_temp,&integrals[0]);
        break;
     case(ReducedNodes):
        temp=escript::Data( arg, function(asAbstractContinuousDomain()) );
        _temp=temp.getDataC();
        Finley_Assemble_integrate(mesh->Nodes,mesh->Elements,&_temp,&integrals[0]);
        break;
     case(Elements):
        Finley_Assemble_integrate(mesh->Nodes,mesh->Elements,&_arg,&integrals[0]);
        break;
     case(ReducedElements):
        Finley_Assemble_integrate(mesh->Nodes,mesh->Elements,&_arg,&integrals[0]);
        break;
     case(FaceElements):
        Finley_Assemble_integrate(mesh->Nodes,mesh->FaceElements,&_arg,&integrals[0]);
        break;
     case(ReducedFaceElements):
        Finley_Assemble_integrate(mesh->Nodes,mesh->FaceElements,&_arg,&integrals[0]);
        break;
     case(Points):
        throw FinleyAdapterException("Error - Integral of data on points is not supported.");
        break;
     case(ContactElementsZero):
        Finley_Assemble_integrate(mesh->Nodes,mesh->ContactElements,&_arg,&integrals[0]);
        break;
     case(ReducedContactElementsZero):
        Finley_Assemble_integrate(mesh->Nodes,mesh->ContactElements,&_arg,&integrals[0]);
        break;
     case(ContactElementsOne):
        Finley_Assemble_integrate(mesh->Nodes,mesh->ContactElements,&_arg,&integrals[0]);
        break;
     case(ReducedContactElementsOne):
        Finley_Assemble_integrate(mesh->Nodes,mesh->ContactElements,&_arg,&integrals[0]);
        break;
     case(DegreesOfFreedom):
        temp=escript::Data( arg, function(asAbstractContinuousDomain()) );
        _temp=temp.getDataC();
        Finley_Assemble_integrate(mesh->Nodes,mesh->Elements,&_temp,&integrals[0]);
        break;
     case(ReducedDegreesOfFreedom):
        temp=escript::Data( arg, function(asAbstractContinuousDomain()) );
        _temp=temp.getDataC();
        Finley_Assemble_integrate(mesh->Nodes,mesh->Elements,&_temp,&integrals[0]);
        break;
     default:
        stringstream temp;
        temp << "Error - Integrals: Finley does not know anything about function space type " << arg.getFunctionSpace().getTypeCode();
        throw FinleyAdapterException(temp.str());
        break;
  }
  checkFinleyError();
  blocktimer_increment("integrate()", blocktimer_start);
}

//
// calculates the gradient of arg:
//
void MeshAdapter::setToGradient(escript::Data& grad,const escript::Data& arg) const
{
  const MeshAdapter& argDomain=dynamic_cast<const MeshAdapter&>(arg.getFunctionSpace().getDomain());
  if (argDomain!=*this)
    throw FinleyAdapterException("Error - Illegal domain of gradient argument");
  const MeshAdapter& gradDomain=dynamic_cast<const MeshAdapter&>(grad.getFunctionSpace().getDomain());
  if (gradDomain!=*this)
     throw FinleyAdapterException("Error - Illegal domain of gradient");

  Finley_Mesh* mesh=m_finleyMesh.get();
  escriptDataC _grad=grad.getDataC();
  escriptDataC nodeDataC;
  escript::Data temp;
  if (getMPISize()>1) {
      if( arg.getFunctionSpace().getTypeCode() == DegreesOfFreedom ) {
        temp=escript::Data( arg,  continuousFunction(asAbstractContinuousDomain()) );
        nodeDataC = temp.getDataC();
      } else if( arg.getFunctionSpace().getTypeCode() == ReducedDegreesOfFreedom ) {
        temp=escript::Data( arg,  reducedContinuousFunction(asAbstractContinuousDomain()) );
        nodeDataC = temp.getDataC();
      } else {
        nodeDataC = arg.getDataC();
      }
  } else {
     nodeDataC = arg.getDataC();
  }
  switch(grad.getFunctionSpace().getTypeCode()) {
       case(Nodes):
          throw FinleyAdapterException("Error - Gradient at nodes is not supported.");
          break;
       case(ReducedNodes):
          throw FinleyAdapterException("Error - Gradient at reduced nodes is not supported.");
          break;
       case(Elements):
          Finley_Assemble_gradient(mesh->Nodes,mesh->Elements,&_grad,&nodeDataC);
          break;
       case(ReducedElements):
          Finley_Assemble_gradient(mesh->Nodes,mesh->Elements,&_grad,&nodeDataC);
          break;
       case(FaceElements):
          Finley_Assemble_gradient(mesh->Nodes,mesh->FaceElements,&_grad,&nodeDataC);
          break;
       case(ReducedFaceElements):
          Finley_Assemble_gradient(mesh->Nodes,mesh->FaceElements,&_grad,&nodeDataC);
          break;
       case(Points):
          throw FinleyAdapterException("Error - Gradient at points is not supported.");
          break;
       case(ContactElementsZero):
          Finley_Assemble_gradient(mesh->Nodes,mesh->ContactElements,&_grad,&nodeDataC);
          break;
       case(ReducedContactElementsZero):
          Finley_Assemble_gradient(mesh->Nodes,mesh->ContactElements,&_grad,&nodeDataC);
          break;
       case(ContactElementsOne):
          Finley_Assemble_gradient(mesh->Nodes,mesh->ContactElements,&_grad,&nodeDataC);
          break;
       case(ReducedContactElementsOne):
          Finley_Assemble_gradient(mesh->Nodes,mesh->ContactElements,&_grad,&nodeDataC);
          break;
       case(DegreesOfFreedom):
          throw FinleyAdapterException("Error - Gradient at degrees of freedom is not supported.");
          break;
       case(ReducedDegreesOfFreedom):
          throw FinleyAdapterException("Error - Gradient at reduced degrees of freedom is not supported.");
          break;
       default:
          stringstream temp;
          temp << "Error - Gradient: Finley does not know anything about function space type " << arg.getFunctionSpace().getTypeCode();
          throw FinleyAdapterException(temp.str());
          break;
  }
  checkFinleyError();
}

//
// returns the size of elements:
//
void MeshAdapter::setToSize(escript::Data& size) const
{
  Finley_Mesh* mesh=m_finleyMesh.get();
  escriptDataC tmp=size.getDataC();
  switch(size.getFunctionSpace().getTypeCode()) {
       case(Nodes):
          throw FinleyAdapterException("Error - Size of nodes is not supported.");
          break;
       case(ReducedNodes):
          throw FinleyAdapterException("Error - Size of reduced nodes is not supported.");
          break;
       case(Elements):
          Finley_Assemble_getSize(mesh->Nodes,mesh->Elements,&tmp);
          break;
       case(ReducedElements):
          Finley_Assemble_getSize(mesh->Nodes,mesh->Elements,&tmp);
          break;
       case(FaceElements):
          Finley_Assemble_getSize(mesh->Nodes,mesh->FaceElements,&tmp);
          break;
       case(ReducedFaceElements):
          Finley_Assemble_getSize(mesh->Nodes,mesh->FaceElements,&tmp);
          break;
       case(Points):
          throw FinleyAdapterException("Error - Size of point elements is not supported.");
          break;
       case(ContactElementsZero):
       case(ContactElementsOne):
          Finley_Assemble_getSize(mesh->Nodes,mesh->ContactElements,&tmp);
          break;
       case(ReducedContactElementsZero):
       case(ReducedContactElementsOne):
          Finley_Assemble_getSize(mesh->Nodes,mesh->ContactElements,&tmp);
          break;
       case(DegreesOfFreedom):
          throw FinleyAdapterException("Error - Size of degrees of freedom is not supported.");
          break;
       case(ReducedDegreesOfFreedom):
          throw FinleyAdapterException("Error - Size of reduced degrees of freedom is not supported.");
          break;
       default:
          stringstream temp;
          temp << "Error - Element size: Finley does not know anything about function space type " << size.getFunctionSpace().getTypeCode();
          throw FinleyAdapterException(temp.str());
          break;
  }
  checkFinleyError();
}

// sets the location of nodes:
void MeshAdapter::setNewX(const escript::Data& new_x)
{
  Finley_Mesh* mesh=m_finleyMesh.get();
  escriptDataC tmp;
  const MeshAdapter& newDomain=dynamic_cast<const MeshAdapter&>(new_x.getFunctionSpace().getDomain());
  if (newDomain!=*this) 
     throw FinleyAdapterException("Error - Illegal domain of new point locations");
  tmp = new_x.getDataC();
  Finley_Mesh_setCoordinates(mesh,&tmp);
  checkFinleyError();
}

// saves a data array in openDX format:
void MeshAdapter::saveDX(const std::string& filename,const boost::python::dict& arg) const
{
    int MAX_namelength=256;
    const int num_data=boost::python::extract<int>(arg.attr("__len__")());
  /* win32 refactor */
  char* *names = (num_data>0) ? TMPMEMALLOC(num_data,char*) : (char**)NULL;
  for(int i=0;i<num_data;i++)
  {
  	names[i] = (MAX_namelength>0) ? TMPMEMALLOC(MAX_namelength,char) : (char*)NULL;
  }

  char* *c_names = (num_data>0) ? TMPMEMALLOC(num_data,char*) : (char**)NULL;
  escriptDataC *data = (num_data>0) ? TMPMEMALLOC(num_data,escriptDataC) : (escriptDataC*)NULL;
  escriptDataC* *ptr_data = (num_data>0) ? TMPMEMALLOC(num_data,escriptDataC*) : (escriptDataC**)NULL;

  boost::python::list keys=arg.keys();
  for (int i=0;i<num_data;++i) {
         std::string n=boost::python::extract<std::string>(keys[i]);
         escript::Data& d=boost::python::extract<escript::Data&>(arg[keys[i]]);
         if (dynamic_cast<const MeshAdapter&>(d.getFunctionSpace().getDomain()) !=*this) 
             throw FinleyAdapterException("Error  in saveVTK: Data must be defined on same Domain");
         data[i]=d.getDataC();
         ptr_data[i]=&(data[i]);
         c_names[i]=&(names[i][0]);
         if (n.length()>MAX_namelength-1) {
            strncpy(c_names[i],n.c_str(),MAX_namelength-1);
            c_names[i][MAX_namelength-1]='\0';
         } else {
            strcpy(c_names[i],n.c_str());
         }
    }
    Finley_Mesh_saveDX(filename.c_str(),m_finleyMesh.get(),num_data,c_names,ptr_data);
    checkFinleyError();
    
      /* win32 refactor */
  TMPMEMFREE(c_names);
  TMPMEMFREE(data);
  TMPMEMFREE(ptr_data);
  for(int i=0;i<num_data;i++)
  {
  	TMPMEMFREE(names[i]);
  }
  TMPMEMFREE(names);

    return;
}

// saves a data array in openVTK format:
void MeshAdapter::saveVTK(const std::string& filename,const boost::python::dict& arg) const
{
    int MAX_namelength=256;
    const int num_data=boost::python::extract<int>(arg.attr("__len__")());
  /* win32 refactor */
  char* *names = (num_data>0) ? TMPMEMALLOC(num_data,char*) : (char**)NULL;
  for(int i=0;i<num_data;i++)
  {
  	names[i] = (MAX_namelength>0) ? TMPMEMALLOC(MAX_namelength,char) : (char*)NULL;
  }

  char* *c_names = (num_data>0) ? TMPMEMALLOC(num_data,char*) : (char**)NULL;
  escriptDataC *data = (num_data>0) ? TMPMEMALLOC(num_data,escriptDataC) : (escriptDataC*)NULL;
  escriptDataC* *ptr_data = (num_data>0) ? TMPMEMALLOC(num_data,escriptDataC*) : (escriptDataC**)NULL;

    boost::python::list keys=arg.keys();
    for (int i=0;i<num_data;++i) {
         std::string n=boost::python::extract<std::string>(keys[i]);
         escript::Data& d=boost::python::extract<escript::Data&>(arg[keys[i]]);
         if (dynamic_cast<const MeshAdapter&>(d.getFunctionSpace().getDomain()) !=*this) 
             throw FinleyAdapterException("Error  in saveVTK: Data must be defined on same Domain");
         data[i]=d.getDataC();
         ptr_data[i]=&(data[i]);
         c_names[i]=&(names[i][0]);
         if (n.length()>MAX_namelength-1) {
            strncpy(c_names[i],n.c_str(),MAX_namelength-1);
            c_names[i][MAX_namelength-1]='\0';
         } else {
            strcpy(c_names[i],n.c_str());
         }
    }
    Finley_Mesh_saveVTK(filename.c_str(),m_finleyMesh.get(),num_data,c_names,ptr_data);

  checkFinleyError();
  /* win32 refactor */
  TMPMEMFREE(c_names);
  TMPMEMFREE(data);
  TMPMEMFREE(ptr_data);
  for(int i=0;i<num_data;i++)
  {
  	TMPMEMFREE(names[i]);
  }
  TMPMEMFREE(names);

    return;
}
                                                                                                                                                                        
                                                                                                                                                                        
// creates a SystemMatrixAdapter stiffness matrix an initializes it with zeros:
SystemMatrixAdapter MeshAdapter::newSystemMatrix(
                      const int row_blocksize,
                      const escript::FunctionSpace& row_functionspace,
                      const int column_blocksize,
                      const escript::FunctionSpace& column_functionspace,
                      const int type) const
{
    int reduceRowOrder=0;
    int reduceColOrder=0;
    // is the domain right?
    const MeshAdapter& row_domain=dynamic_cast<const MeshAdapter&>(row_functionspace.getDomain());
    if (row_domain!=*this) 
          throw FinleyAdapterException("Error - domain of row function space does not match the domain of matrix generator.");
    const MeshAdapter& col_domain=dynamic_cast<const MeshAdapter&>(column_functionspace.getDomain());
    if (col_domain!=*this) 
          throw FinleyAdapterException("Error - domain of columnn function space does not match the domain of matrix generator.");
    // is the function space type right 
    if (row_functionspace.getTypeCode()==DegreesOfFreedom) {
        reduceRowOrder=0;
    } else if (row_functionspace.getTypeCode()==ReducedDegreesOfFreedom) {
        reduceRowOrder=1;
    } else {
        throw FinleyAdapterException("Error - illegal function space type for system matrix rows.");
    }
    if (column_functionspace.getTypeCode()==DegreesOfFreedom) {
        reduceColOrder=0;
    } else if (column_functionspace.getTypeCode()==ReducedDegreesOfFreedom) {
        reduceColOrder=1;
    } else {
        throw FinleyAdapterException("Error - illegal function space type for system matrix columns.");
    }
    // generate matrix:
    
    Paso_SystemMatrixPattern* fsystemMatrixPattern=Finley_getPattern(getFinley_Mesh(),reduceRowOrder,reduceColOrder);
    checkFinleyError();
    Paso_SystemMatrix* fsystemMatrix;
    int trilinos = 0;
    if (trilinos) {
#ifdef TRILINOS
      /* Allocation Epetra_VrbMatrix here */
#endif
    }
    else {
      fsystemMatrix=Paso_SystemMatrix_alloc(type,fsystemMatrixPattern,row_blocksize,column_blocksize);
    }
    checkPasoError();
    Paso_SystemMatrixPattern_free(fsystemMatrixPattern);
    return SystemMatrixAdapter(fsystemMatrix,row_blocksize,row_functionspace,column_blocksize,column_functionspace);
}

//
// vtkObject MeshAdapter::createVtkObject() const
// TODO:
//

//
// returns true if data at the atom_type is considered as being cell centered:
bool MeshAdapter::isCellOriented(int functionSpaceCode) const
{
  switch(functionSpaceCode) {
       case(Nodes):
       case(DegreesOfFreedom):
       case(ReducedDegreesOfFreedom):
          return false;
          break;
       case(Elements):
       case(FaceElements):
       case(Points):
       case(ContactElementsZero):
       case(ContactElementsOne):
       case(ReducedElements):
       case(ReducedFaceElements):
       case(ReducedContactElementsZero):
       case(ReducedContactElementsOne):
          return true;
          break;
       default:
          stringstream temp;
          temp << "Error - Cell: Finley does not know anything about function space type " << functionSpaceCode;
          throw FinleyAdapterException(temp.str());
          break;
  }
  checkFinleyError();
  return false;
}

bool MeshAdapter::probeInterpolationOnDomain(int functionSpaceType_source,int functionSpaceType_target) const
{
  switch(functionSpaceType_source) {
     case(Nodes):
        switch(functionSpaceType_target) {
           case(Nodes):
           case(ReducedNodes):
           case(ReducedDegreesOfFreedom):
           case(DegreesOfFreedom):
           case(Elements):
           case(ReducedElements):
           case(FaceElements):
           case(ReducedFaceElements):
           case(Points):
           case(ContactElementsZero):
           case(ReducedContactElementsZero):
           case(ContactElementsOne):
           case(ReducedContactElementsOne):
               return true;
           default:
               stringstream temp;
               temp << "Error - Interpolation On Domain: Finley does not know anything about function space type " << functionSpaceType_target;
               throw FinleyAdapterException(temp.str());
        }
        break;
     case(ReducedNodes):
        switch(functionSpaceType_target) {
           case(ReducedNodes):
           case(ReducedDegreesOfFreedom):
           case(Elements):
           case(ReducedElements):
           case(FaceElements):
           case(ReducedFaceElements):
           case(Points):
           case(ContactElementsZero):
           case(ReducedContactElementsZero):
           case(ContactElementsOne):
           case(ReducedContactElementsOne):
               return true;
          case(Nodes):
          case(DegreesOfFreedom):
             return false;
           default:
               stringstream temp;
               temp << "Error - Interpolation On Domain: Finley does not know anything about function space type " << functionSpaceType_target;
               throw FinleyAdapterException(temp.str());
        }
        break;
     case(Elements):
        if (functionSpaceType_target==Elements) {
           return true;
        } else if (functionSpaceType_target==ReducedElements) {
           return true;
        } else {
           return false;
        }
     case(ReducedElements):
        if (functionSpaceType_target==ReducedElements) {
           return true;
        } else {
           return false;
        }
     case(FaceElements):
        if (functionSpaceType_target==FaceElements) {
           return true;
        } else if (functionSpaceType_target==ReducedFaceElements) {
           return true;
        } else {
           return false;
        }
     case(ReducedFaceElements):
        if (functionSpaceType_target==ReducedFaceElements) {
           return true;
        } else {
           return false;
        }
     case(Points):
        if (functionSpaceType_target==Points) {
           return true;
        } else {
           return false;
        }
     case(ContactElementsZero):
     case(ContactElementsOne):
        if (functionSpaceType_target==ContactElementsZero || functionSpaceType_target==ContactElementsOne) {
           return true;
        } else if (functionSpaceType_target==ReducedContactElementsZero || functionSpaceType_target==ReducedContactElementsOne) {
           return true;
        } else {
           return false;
        }
     case(ReducedContactElementsZero):
     case(ReducedContactElementsOne):
        if (functionSpaceType_target==ReducedContactElementsZero || functionSpaceType_target==ReducedContactElementsOne) {
           return true;
        } else {
           return false;
        }
     case(DegreesOfFreedom):
        switch(functionSpaceType_target) {
           case(ReducedDegreesOfFreedom):
           case(DegreesOfFreedom):
           case(Nodes):
           case(ReducedNodes):
           case(Elements):
           case(ReducedElements):
           case(Points):
           case(FaceElements):
           case(ReducedFaceElements):
           case(ContactElementsZero):
           case(ReducedContactElementsZero):
           case(ContactElementsOne):
           case(ReducedContactElementsOne):
              return true;
           default:
             stringstream temp;
             temp << "Error - Interpolation On Domain: Finley does not know anything about function space type " << functionSpaceType_target;
             throw FinleyAdapterException(temp.str());
        }
        break;
     case(ReducedDegreesOfFreedom):
       switch(functionSpaceType_target) {
          case(ReducedDegreesOfFreedom):
          case(ReducedNodes):
          case(Elements):
          case(ReducedElements):
          case(FaceElements):
          case(ReducedFaceElements):
          case(Points):
          case(ContactElementsZero):
          case(ReducedContactElementsZero):
          case(ContactElementsOne):
          case(ReducedContactElementsOne):
              return true;
          case(Nodes):
          case(DegreesOfFreedom):
             return false;
          default:
             stringstream temp;
             temp << "Error - Interpolation On Domain: Finley does not know anything about function space type " << functionSpaceType_target;
             throw FinleyAdapterException(temp.str());
       }
       break;
     default:
        stringstream temp;
        temp << "Error - Interpolation On Domain: Finley does not know anything about function space type " << functionSpaceType_source;
        throw FinleyAdapterException(temp.str());
        break;
  }
  checkFinleyError();
  return false;
}

bool MeshAdapter::probeInterpolationACross(int functionSpaceType_source,const AbstractDomain& targetDomain, int functionSpaceType_target) const
{
    return false;
}

bool MeshAdapter::operator==(const AbstractDomain& other) const
{
  const MeshAdapter* temp=dynamic_cast<const MeshAdapter*>(&other);
  if (temp!=0) {
    return (m_finleyMesh==temp->m_finleyMesh);
  } else {
    return false;
  }
}

bool MeshAdapter::operator!=(const AbstractDomain& other) const
{
  return !(operator==(other));
}

int MeshAdapter::getSystemMatrixTypeId(const int solver, const int package, const bool symmetry) const
{
   int out=Paso_SystemMatrix_getSystemMatrixTypeId(SystemMatrixAdapter::mapOptionToPaso(solver),SystemMatrixAdapter::mapOptionToPaso(package),symmetry?1:0);
   checkPasoError();
   return out;
}

escript::Data MeshAdapter::getX() const
{
  return continuousFunction(asAbstractContinuousDomain()).getX();
}

escript::Data MeshAdapter::getNormal() const
{
  return functionOnBoundary(asAbstractContinuousDomain()).getNormal();
}

escript::Data MeshAdapter::getSize() const
{
  return function(asAbstractContinuousDomain()).getSize();
}

int* MeshAdapter::borrowSampleReferenceIDs(int functionSpaceType) const
{
  int *out=0,i;
  Finley_Mesh* mesh=m_finleyMesh.get();
  switch (functionSpaceType) {
    case(Nodes):
      out=mesh->Nodes->Id;
      break;
    case(ReducedNodes):
      out=mesh->Nodes->reducedNodesId;
      break;
    case(Elements):
      out=mesh->Elements->Id;
      break;
    case(ReducedElements):
      out=mesh->Elements->Id;
      break;
    case(FaceElements):
      out=mesh->FaceElements->Id;
      break;
    case(ReducedFaceElements):
      out=mesh->FaceElements->Id;
      break;
    case(Points):
      out=mesh->Points->Id;
      break;
    case(ContactElementsZero):
      out=mesh->ContactElements->Id;
      break;
    case(ReducedContactElementsZero):
      out=mesh->ContactElements->Id;
      break;
    case(ContactElementsOne):
      out=mesh->ContactElements->Id;
      break;
    case(ReducedContactElementsOne):
      out=mesh->ContactElements->Id;
      break;
    case(DegreesOfFreedom):
      out=mesh->Nodes->degreesOfFreedomId;
      break;
    case(ReducedDegreesOfFreedom):
      out=mesh->Nodes->reducedDegreesOfFreedomId;
      break;
    default:
      stringstream temp;
      temp << "Error - Invalid function space type: " << functionSpaceType << " for domain: " << getDescription();
      throw FinleyAdapterException(temp.str());
      break;
  }
  return out;
}
int MeshAdapter::getTagFromSampleNo(int functionSpaceType, int sampleNo) const
{
  int out=0;
  Finley_Mesh* mesh=m_finleyMesh.get();
  switch (functionSpaceType) {
  case(Nodes):
    out=mesh->Nodes->Tag[sampleNo];
    break;
  case(ReducedNodes):
    throw FinleyAdapterException(" Error - ReducedNodes does not support tags.");
    break;
  case(Elements):
    out=mesh->Elements->Tag[sampleNo];
    break;
  case(ReducedElements):
    out=mesh->Elements->Tag[sampleNo];
    break;
  case(FaceElements):
    out=mesh->FaceElements->Tag[sampleNo];
    break;
  case(ReducedFaceElements):
    out=mesh->FaceElements->Tag[sampleNo];
    break;
  case(Points):
    out=mesh->Points->Tag[sampleNo];
    break;
  case(ContactElementsZero):
    out=mesh->ContactElements->Tag[sampleNo];
    break;
  case(ReducedContactElementsZero):
    out=mesh->ContactElements->Tag[sampleNo];
    break;
  case(ContactElementsOne):
    out=mesh->ContactElements->Tag[sampleNo];
    break;
  case(ReducedContactElementsOne):
    out=mesh->ContactElements->Tag[sampleNo];
    break;
  case(DegreesOfFreedom):
    throw FinleyAdapterException(" Error - DegreesOfFreedom does not support tags.");
    break;
  case(ReducedDegreesOfFreedom):
    throw FinleyAdapterException(" Error - ReducedDegreesOfFreedom does not support tags.");
    break;
  default:
    stringstream temp;
    temp << "Error - Invalid function space type: " << functionSpaceType << " for domain: " << getDescription();
    throw FinleyAdapterException(temp.str());
    break;
  }
  return out;
}


void MeshAdapter::setTags(const int functionSpaceType, const int newTag, const escript::Data& mask) const
{
  Finley_Mesh* mesh=m_finleyMesh.get();
  escriptDataC tmp=mask.getDataC();
  switch(functionSpaceType) {
       case(Nodes):
          Finley_NodeFile_setTags(mesh->Nodes,newTag,&tmp);
          break;
       case(ReducedNodes):
          throw FinleyAdapterException("Error - ReducedNodes does not support tags");
          break;
       case(DegreesOfFreedom):
          throw FinleyAdapterException("Error - DegreesOfFreedom does not support tags");
          break;
       case(ReducedDegreesOfFreedom):
          throw FinleyAdapterException("Error - ReducedDegreesOfFreedom does not support tags");
          break;
       case(Elements):
          Finley_ElementFile_setTags(mesh->Elements,newTag,&tmp);
          break;
       case(ReducedElements):
          Finley_ElementFile_setTags(mesh->Elements,newTag,&tmp);
          break;
       case(FaceElements):
          Finley_ElementFile_setTags(mesh->FaceElements,newTag,&tmp);
          break;
       case(ReducedFaceElements):
          Finley_ElementFile_setTags(mesh->FaceElements,newTag,&tmp);
          break;
       case(Points):
          Finley_ElementFile_setTags(mesh->Points,newTag,&tmp);
          break;
       case(ContactElementsZero):
          Finley_ElementFile_setTags(mesh->ContactElements,newTag,&tmp);
          break;
       case(ReducedContactElementsZero):
          Finley_ElementFile_setTags(mesh->ContactElements,newTag,&tmp);
          break;
       case(ContactElementsOne):
          Finley_ElementFile_setTags(mesh->ContactElements,newTag,&tmp);
          break;
       case(ReducedContactElementsOne):
          Finley_ElementFile_setTags(mesh->ContactElements,newTag,&tmp);
          break;
       default:
          stringstream temp;
          temp << "Error - Finley does not know anything about function space type " << functionSpaceType;
          throw FinleyAdapterException(temp.str());
  }
  checkFinleyError();
  return;
}

void MeshAdapter::setTagMap(const std::string& name,  int tag)
{
  Finley_Mesh* mesh=m_finleyMesh.get();
  Finley_Mesh_addTagMap(mesh, name.c_str(),tag);
  checkPasoError();
  // throwStandardException("MeshAdapter::set TagMap is not implemented.");
}

int MeshAdapter::getTag(const std::string& name) const
{
  Finley_Mesh* mesh=m_finleyMesh.get();
  int tag=0;
  tag=Finley_Mesh_getTag(mesh, name.c_str());
  checkPasoError();
  // throwStandardException("MeshAdapter::getTag is not implemented.");
  return tag;
}

bool MeshAdapter::isValidTagName(const std::string& name) const
{
  Finley_Mesh* mesh=m_finleyMesh.get();
  return Finley_Mesh_isValidTagName(mesh,name.c_str());
}

std::string MeshAdapter::showTagNames() const
{
  stringstream temp;
  Finley_Mesh* mesh=m_finleyMesh.get();
  Finley_TagMap* tag_map=mesh->TagMap;
  while (tag_map) {
     temp << tag_map->name;
     tag_map=tag_map->next;
     if (tag_map) temp << ", ";
  }
  return temp.str();
}

}  // end of namespace
