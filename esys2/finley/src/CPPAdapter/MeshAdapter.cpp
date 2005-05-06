// $Id$
/*
 ******************************************************************************
 *                                                                            *
 *       COPYRIGHT  ACcESS 2004 -  All Rights Reserved                        *
 *                                                                            *
 * This software is the property of ACcESS. No part of this code              *
 * may be copied in any form or by any means without the expressed written    *
 * consent of ACcESS.  Copying, use or modification of this software          *
 * by any unauthorised person is illegal unless that person has a software    *
 * license agreement with ACcESS.                                             *
 *                                                                            *
 ******************************************************************************
*/
extern "C" {
#include "finley/finleyC/Finley.h"
#include "finley/finleyC/Assemble.h"
#include "finley/finleyC/Mesh.h"
#include "finley/finleyC/Finley.h"
#include "finley/finleyC/System.h"
}
#include "finley/CPPAdapter/SystemMatrixAdapter.h"
#include "finley/CPPAdapter/MeshAdapter.h"
#include "finley/CPPAdapter/FinleyError.h"
#include "finley/CPPAdapter/FinleyAdapterException.h"
#include "escript/Data/FunctionSpaceFactory.h"
#include "escript/Data/Data.h"
#include "escript/Data/DataArrayView.h"
#include "escript/Data/FunctionSpace.h"
#include "escript/Data/DataFactory.h"
#include <iostream>
#include <vector>
#include <sstream>

using namespace std;
using namespace escript;

namespace finley {

struct null_deleter
{
  void operator()(void const *ptr) const
  {
  }
};

//
// define the statics
MeshAdapter::FunctionSpaceNamesMapType MeshAdapter::m_functionSpaceTypeNames;
const int MeshAdapter::DegreesOfFreedom=FINLEY_DEGREES_OF_FREEDOM;
const int MeshAdapter::ReducedDegreesOfFreedom=FINLEY_REDUCED_DEGREES_OF_FREEDOM;
const int MeshAdapter::Nodes=FINLEY_NODES;
const int MeshAdapter::Elements=FINLEY_ELEMENTS;
const int MeshAdapter::FaceElements=FINLEY_FACE_ELEMENTS;
const int MeshAdapter::Points=FINLEY_POINTS;
const int MeshAdapter::ContactElementsZero=FINLEY_CONTACT_ELEMENTS_1;
const int MeshAdapter::ContactElementsOne=FINLEY_CONTACT_ELEMENTS_2;

MeshAdapter::MeshAdapter(Finley_Mesh* finleyMesh)
{
  setFunctionSpaceTypeNames();
  //
  // need to use a null_deleter as Finley_Mesh_dealloc deletes the pointer
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
    //   cout << "Calling dealloc." << endl;
    Finley_Mesh_dealloc(m_finleyMesh.get());
    //   cout << "Finished dealloc." << endl;
  }
}

Finley_Mesh* MeshAdapter::getFinley_Mesh() const {
   return m_finleyMesh.get();
}

void MeshAdapter::write(const std::string& fileName) const
{
  char fName[fileName.size()+1];
  strcpy(fName,fileName.c_str());
  Finley_Mesh_write(m_finleyMesh.get(),fName);
  checkFinleyError();
}

// void MeshAdapter::getTagList(int functionSpaceType,
// 				    int* numTags) const
// {
//   Finley_Mesh_tagList(m_finleyMesh.get(),functionSpaceType,tagList,numTags);
//   return;
// }

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
    (FunctionSpaceNamesMapType::value_type(Elements,"Finley_Elements"));
  m_functionSpaceTypeNames.insert
    (FunctionSpaceNamesMapType::value_type(FaceElements,"Finley_Face_Elements"));
  m_functionSpaceTypeNames.insert
    (FunctionSpaceNamesMapType::value_type(Points,"Finley_Points"));
  m_functionSpaceTypeNames.insert
    (FunctionSpaceNamesMapType::value_type(ContactElementsZero,"Finley_Contact_Elements_0"));
  m_functionSpaceTypeNames.insert
    (FunctionSpaceNamesMapType::value_type(ContactElementsOne,"Finley_Contact_Elements_1"));
}

int MeshAdapter::getContinuousFunctionCode() const
{
  return Nodes;
}
int MeshAdapter::getFunctionCode() const
{
  return Elements;
}
int MeshAdapter::getFunctionOnBoundaryCode() const
{
  return FaceElements;
}
int MeshAdapter::getFunctionOnContactZeroCode() const
{
  return ContactElementsZero;
}
int MeshAdapter::getFunctionOnContactOneCode() const
{
  return ContactElementsOne;
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
// returns a pointer to the tag list of samples of functionSpaceType
//
void MeshAdapter::getTagList(int functionSpaceType, int** tagList,
			     int* numTags) const
{
  *tagList=NULL;
  *numTags=0;
  Finley_Mesh* mesh=m_finleyMesh.get();
  switch (functionSpaceType) {
  case(Nodes):
    if (mesh->Nodes!=NULL) {
      *tagList=mesh->Nodes->Tag;
      *numTags=mesh->Nodes->numNodes;
    }
    break;
  case(Elements):
    if (mesh->Elements!=NULL) {
      *tagList=mesh->Elements->Tag;
      *numTags=mesh->Elements->numElements;
    }
    break;
  case(FaceElements):
    if (mesh->FaceElements!=NULL) {
      *tagList=mesh->FaceElements->Tag;
      *numTags=mesh->FaceElements->numElements;
    }
    break;
  case(Points):
    if (mesh->Points!=NULL) {
      *tagList=mesh->Points->Tag;
      *numTags=mesh->Points->numElements;
    }
    break;
  case(ContactElementsZero):
    if (mesh->ContactElements!=NULL) {
      *tagList=mesh->ContactElements->Tag;
      *numTags=mesh->ContactElements->numElements;
    }
    break;
  case(ContactElementsOne):
    if (mesh->ContactElements!=NULL) {
      *tagList=mesh->ContactElements->Tag;
      *numTags=mesh->ContactElements->numElements;
    }
    break;
  case(DegreesOfFreedom):
    if (mesh->Nodes!=NULL) {
      *tagList=NULL;
      *numTags=0;
    }
    break;
  case(ReducedDegreesOfFreedom):
    if (mesh->Nodes!=NULL) {
      *tagList=NULL;
      *numTags=0;
    }
    break;
  default:
    stringstream temp;
    temp << "Error - Invalid function space type: "
	 << functionSpaceType << " for domain: " << getDescription();
    throw FinleyAdapterException(temp.str());
    break;
  }
  if (*tagList==NULL) {
    stringstream temp;
    temp << "Error - no tags available for " << functionSpaceType << " for domain: " << getDescription();
    throw FinleyAdapterException(temp.str());
  }
  return;
}
//
// returns a pointer to the reference no list of samples of functionSpaceType
//
void MeshAdapter::getReferenceNoList(int functionSpaceType, int** referenceNoList,
			     int* numReferenceNo) const
{
  *referenceNoList=NULL;
  *numReferenceNo=0;
  Finley_Mesh* mesh=m_finleyMesh.get();
  switch (functionSpaceType) {
  case(Nodes):
    if (mesh->Nodes!=NULL) {
      *referenceNoList=mesh->Nodes->Id;
      *numReferenceNo=mesh->Nodes->numNodes;
    }
    break;
  case(Elements):
    if (mesh->Elements!=NULL) {
      *referenceNoList=mesh->Elements->Id;
      *numReferenceNo=mesh->Elements->numElements;
    }
    break;
  case(FaceElements):
    if (mesh->FaceElements!=NULL) {
      *referenceNoList=mesh->FaceElements->Id;
      *numReferenceNo=mesh->FaceElements->numElements;
    }
    break;
  case(Points):
    if (mesh->Points!=NULL) {
      *referenceNoList=mesh->Points->Id;
      *numReferenceNo=mesh->Points->numElements;
    }
    break;
  case(ContactElementsZero):
    if (mesh->ContactElements!=NULL) {
      *referenceNoList=mesh->ContactElements->Id;
      *numReferenceNo=mesh->ContactElements->numElements;
    }
    break;
  case(ContactElementsOne):
    if (mesh->ContactElements!=NULL) {
      *referenceNoList=mesh->ContactElements->Id;
      *numReferenceNo=mesh->ContactElements->numElements;
    }
    break;
  case(DegreesOfFreedom):
    if (mesh->Nodes!=NULL) {
      *referenceNoList=NULL;
      *numReferenceNo=0;
    }
    break;
  case(ReducedDegreesOfFreedom):
    if (mesh->Nodes!=NULL) {
      *referenceNoList=NULL;
      *numReferenceNo=0;
    }
    break;
  default:
    stringstream temp;
    temp << "Error - Invalid function space type: "
	 << functionSpaceType << " for domain: " << getDescription();
    throw FinleyAdapterException(temp.str());
    break;
  }
  if (*referenceNoList==NULL) {
    stringstream temp;
    temp << "Error - reference number list available for " << functionSpaceType << " for domain: " << getDescription();
    throw FinleyAdapterException(temp.str());
  }
  return;
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
           if (mesh->Nodes!=NULL) numSamples=mesh->Nodes->numNodes;
           break;
      case(Elements):
           if (mesh->Elements!=NULL) {
             numSamples=mesh->Elements->numElements;
             numDataPointsPerSample=mesh->Elements->ReferenceElement->numQuadNodes;
           }
           break;
      case(FaceElements):
           if (mesh->FaceElements!=NULL) {
                numDataPointsPerSample=mesh->FaceElements->ReferenceElement->numQuadNodes;
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
      case(ContactElementsOne):
           if (mesh->ContactElements!=NULL) {
             numDataPointsPerSample=mesh->ContactElements->ReferenceElement->numQuadNodes;
             numSamples=mesh->ContactElements->numElements;
           }
           break;
      case(DegreesOfFreedom):
           if (mesh->Nodes!=NULL) {
             numDataPointsPerSample=1;
             numSamples=mesh->Nodes->numDegreesOfFreedom;
           }
           break;
      case(ReducedDegreesOfFreedom):
           if (mesh->Nodes!=NULL) {
             numDataPointsPerSample=1;
             numSamples=mesh->Nodes->reducedNumDegreesOfFreedom;
           }
           break;
      default:
        stringstream temp;
        temp << "Error - Invalid function space type: "
	     << functionSpaceCode << " for domain: " << getDescription();
        throw FinleyAdapterException(temp.str());
        break;
      }
      return pair<int,int>(numDataPointsPerSample,numSamples);
}
//
// adds linear PDE of second order into a given stiffness matrix and right hand side:
//
void MeshAdapter::addPDEToSystem(
                     SystemMatrixAdapter& mat, Data& rhs,
                     const Data& A, const Data& B, const Data& C,const  Data& D,const  Data& X,const  Data& Y,
                     const Data& d, const Data& y, 
                     const Data& d_contact,const Data& y_contact) const
{
   Finley_Mesh* mesh=m_finleyMesh.get();
   Finley_Assemble_PDE(mesh->Nodes,mesh->Elements,mat.getFinley_SystemMatrix(),&(rhs.getDataC()),
                       &(A.getDataC()),&(B.getDataC()),&(C.getDataC()),&(D.getDataC()),&(X.getDataC()),&(Y.getDataC()));
   checkFinleyError();
   Finley_Assemble_RobinCondition(mesh->Nodes,mesh->FaceElements,
				  mat.getFinley_SystemMatrix(),
				  &(rhs.getDataC()),
                                  &(d.getDataC()),&(y.getDataC()),
                                  Finley_Assemble_handelShapeMissMatch_Mean_out);
   checkFinleyError();
   Finley_Assemble_RobinCondition(mesh->Nodes,mesh->FaceElements,
				  mat.getFinley_SystemMatrix(),
				  &(rhs.getDataC()),
                                  &(d_contact.getDataC()),
				  &(y_contact.getDataC()),
                                  Finley_Assemble_handelShapeMissMatch_Step_out);
   checkFinleyError();
}
//
// adds linear PDE of second order into the right hand side only
//
void MeshAdapter::addPDEToRHS( Data& rhs,
                     const  Data& X,const  Data& Y, const Data& y, const Data& y_contact) const
{
   Finley_Mesh* mesh=m_finleyMesh.get();
   Finley_Assemble_PDE(mesh->Nodes,mesh->Elements,0,&(rhs.getDataC()),0,0,0,0,&(X.getDataC()),&(Y.getDataC()));
   checkFinleyError();
   Finley_Assemble_RobinCondition(mesh->Nodes,mesh->FaceElements,0,&(rhs.getDataC()),0,&(y.getDataC()),
                                  Finley_Assemble_handelShapeMissMatch_Mean_out);
// cout << "Calling :addPDEToRHS." << endl;
   checkFinleyError();
   Finley_Assemble_RobinCondition(mesh->Nodes,mesh->FaceElements,0,&(rhs.getDataC()),0,&(y_contact.getDataC()),
                                  Finley_Assemble_handelShapeMissMatch_Step_out);
// cout << "Calling :addPDEToRHS." << endl;
   checkFinleyError();
}
//
// interpolates data between different function spaces:
//
void MeshAdapter::interpolateOnDomain(Data& target,const Data& in) const
{
  const MeshAdapter& inDomain=dynamic_cast<const MeshAdapter&>(in.getFunctionSpace().getDomain());
  const MeshAdapter& targetDomain=dynamic_cast<const MeshAdapter&>(target.getFunctionSpace().getDomain());
  if (inDomain!=*this) 
     throw FinleyAdapterException("Error - Illegal domain of interpolant.");
  if (targetDomain!=*this) 
     throw FinleyAdapterException("Error - Illegal domain of interpolation target.");

  Finley_Mesh* mesh=m_finleyMesh.get();
  switch(in.getFunctionSpace().getTypeCode()) {
     case(Nodes):
        switch(target.getFunctionSpace().getTypeCode()) {
           case(Nodes):
           case(ReducedDegreesOfFreedom):
           case(DegreesOfFreedom):
               Finley_Assemble_CopyNodalData(mesh->Nodes,&(target.getDataC()),&(in.getDataC()));
               break;
           case(Elements):
               Finley_Assemble_interpolate(mesh->Nodes,mesh->Elements,&(in.getDataC()),&(target.getDataC()));
               break;
           case(FaceElements):
               Finley_Assemble_interpolate(mesh->Nodes,mesh->FaceElements,&(in.getDataC()),&(target.getDataC()));
               break;
           case(Points):
               Finley_Assemble_interpolate(mesh->Nodes,mesh->Points,&(in.getDataC()),&(target.getDataC()));
               break;
           case(ContactElementsZero):
               Finley_Assemble_interpolate(mesh->Nodes,mesh->ContactElements,&(in.getDataC()),&(target.getDataC()));
               break;
           case(ContactElementsOne):
               Finley_Assemble_interpolate(mesh->Nodes,mesh->ContactElements,&(in.getDataC()),&(target.getDataC()));
               break;
           default:
               Finley_ErrorCode=TYPE_ERROR;
               sprintf(Finley_ErrorMsg,"Interpolation on Domain: Finley does not know anything about function space type %d",target.getFunctionSpace().getTypeCode());
               break;
        }
        break;
     case(Elements):
        if (target.getFunctionSpace().getTypeCode()==Elements) {
           Finley_Assemble_CopyElementData(mesh->Elements,&(target.getDataC()),&(in.getDataC()));
        } else {
           Finley_ErrorCode=TYPE_ERROR;
           sprintf(Finley_ErrorMsg,"No interpolation with data on elements possible.");
        }
        break;
     case(FaceElements):
        if (target.getFunctionSpace().getTypeCode()==FaceElements) {
           Finley_Assemble_CopyElementData(mesh->FaceElements,&(target.getDataC()),&(in.getDataC()));
        } else {
           Finley_ErrorCode=TYPE_ERROR;
           sprintf(Finley_ErrorMsg,"No interpolation with data on face elements possible.");
           break;
       }
     case(Points):
        if (target.getFunctionSpace().getTypeCode()==Points) {
           Finley_Assemble_CopyElementData(mesh->Points,&(target.getDataC()),&(in.getDataC()));
        } else {
           Finley_ErrorCode=TYPE_ERROR;
           sprintf(Finley_ErrorMsg,"No interpolation with data on points possible.");
           break;
        }
        break;
     case(ContactElementsZero):
     case(ContactElementsOne):
        if (target.getFunctionSpace().getTypeCode()==ContactElementsZero || target.getFunctionSpace().getTypeCode()==ContactElementsOne) {
           Finley_Assemble_CopyElementData(mesh->ContactElements,&(target.getDataC()),&(in.getDataC()));
        } else {
           Finley_ErrorCode=TYPE_ERROR;
           sprintf(Finley_ErrorMsg,"No interpolation with data on contact elements possible.");
           break;
        }
        break;
     case(DegreesOfFreedom):
        switch(target.getFunctionSpace().getTypeCode()) {
           case(ReducedDegreesOfFreedom):
           case(DegreesOfFreedom):
           case(Nodes):
              Finley_Assemble_CopyNodalData(mesh->Nodes,&(target.getDataC()),&(in.getDataC()));
              break;
           case(Elements):
              Finley_Assemble_interpolate(mesh->Nodes,mesh->Elements,&(in.getDataC()),&(target.getDataC()));
              break;
           case(FaceElements):
              Finley_Assemble_interpolate(mesh->Nodes,mesh->FaceElements,&(in.getDataC()),&(target.getDataC()));
              break;
           case(Points):
              Finley_Assemble_interpolate(mesh->Nodes,mesh->Points,&(in.getDataC()),&(target.getDataC()));
              break;
           case(ContactElementsZero):
           case(ContactElementsOne):
              Finley_Assemble_interpolate(mesh->Nodes,mesh->ContactElements,&(in.getDataC()),&(target.getDataC()));
             break;
           default:
             Finley_ErrorCode=TYPE_ERROR;
             sprintf(Finley_ErrorMsg,"Interpolation On Domain: Finley does not know anything about function space type %d",target.getFunctionSpace().getTypeCode());
             break;
        }
        break;
     case(ReducedDegreesOfFreedom):
       switch(target.getFunctionSpace().getTypeCode()) {
          case(ReducedDegreesOfFreedom):
          case(Nodes):
             Finley_Assemble_CopyNodalData(mesh->Nodes,&(target.getDataC()),&(in.getDataC()));
             break;
          case(Elements):
             Finley_Assemble_interpolate(mesh->Nodes,mesh->Elements,&(in.getDataC()),&(target.getDataC()));
             break;
          case(FaceElements):
             Finley_Assemble_interpolate(mesh->Nodes,mesh->FaceElements,&(in.getDataC()),&(target.getDataC()));
             break;
          case(Points):
             Finley_Assemble_interpolate(mesh->Nodes,mesh->Points,&(in.getDataC()),&(target.getDataC()));
             break;
          case(ContactElementsZero):
          case(ContactElementsOne):
             Finley_Assemble_interpolate(mesh->Nodes,mesh->ContactElements,&(in.getDataC()),&(target.getDataC()));
             break;
          case(DegreesOfFreedom):
             Finley_ErrorCode=TYPE_ERROR;
             sprintf(Finley_ErrorMsg,"Finley does not support interpolation from reduced degrees of freedom to degrees of freedom");
             break;
          default:
             Finley_ErrorCode=TYPE_ERROR;
             sprintf(Finley_ErrorMsg,"Interpolation On Domain: Finley does not know anything about function space type %d",target.getFunctionSpace().getTypeCode());
             break;
       }
       break;
     default:
        Finley_ErrorCode=TYPE_ERROR;
        sprintf(Finley_ErrorMsg,"Interpolation On Domain: Finley does not know anything about function space type %d",in.getFunctionSpace().getTypeCode());
        break;
  }
  checkFinleyError();
}

//
// copies the locations of sample points into x:
//
void MeshAdapter::setToX(Data& arg) const
{
  const MeshAdapter& argDomain=dynamic_cast<const MeshAdapter&>(arg.getFunctionSpace().getDomain());
  if (argDomain!=*this) 
     throw FinleyAdapterException("Error - Illegal domain of data point locations");
  Finley_Mesh* mesh=m_finleyMesh.get();
  // in case of values node coordinates we can do the job directly:
  if (arg.getFunctionSpace().getTypeCode()==Nodes) {
     Finley_Assemble_NodeCoordinates(mesh->Nodes,&(arg.getDataC()));
  } else {
     Data tmp_data=Vector(0.0,continuousFunction(asAbstractContinuousDomain()),true);
     Finley_Assemble_NodeCoordinates(mesh->Nodes,&(tmp_data.getDataC()));
     // this is then interpolated onto arg:
     interpolateOnDomain(arg,tmp_data);
  }
  checkFinleyError();
}
//
// return the normal vectors at the location of data points as a Data object:
//
void MeshAdapter::setToNormal(Data& normal) const
{
  const MeshAdapter& normalDomain=dynamic_cast<const MeshAdapter&>(normal.getFunctionSpace().getDomain());
  if (normalDomain!=*this) 
     throw FinleyAdapterException("Error - Illegal domain of normal locations");
  Finley_Mesh* mesh=m_finleyMesh.get();
  switch(normal.getFunctionSpace().getTypeCode()) {
    case(Nodes):
      Finley_ErrorCode=TYPE_ERROR;
      sprintf(Finley_ErrorMsg,"Finley does not support surface normal vectors for nodes");
      break;
    case(Elements):
      Finley_ErrorCode=TYPE_ERROR;
      sprintf(Finley_ErrorMsg,"Finley does not support surface normal vectors for elements");
      break;
    case (FaceElements):
      Finley_Assemble_setNormal(mesh->Nodes,mesh->FaceElements,&(normal.getDataC()));
      break;
    case(Points):
      Finley_ErrorCode=TYPE_ERROR;
      sprintf(Finley_ErrorMsg,"Finley does not support surface normal vectors for point elements");
      break;
    case (ContactElementsOne):
    case (ContactElementsZero):
      Finley_Assemble_setNormal(mesh->Nodes,mesh->ContactElements,&(normal.getDataC()));
      break;
    case(DegreesOfFreedom):
      Finley_ErrorCode=TYPE_ERROR;
      sprintf(Finley_ErrorMsg,"Finley does not support surface normal vectors for degrees of freedom.");
      break;
    case(ReducedDegreesOfFreedom):
      Finley_ErrorCode=TYPE_ERROR;
      sprintf(Finley_ErrorMsg,"Finley does not support surface normal vectors for reduced degrees of freedom.");
      break;
    default:
      Finley_ErrorCode=TYPE_ERROR;
      sprintf(Finley_ErrorMsg,"Normal Vectors: Finley does not know anything about function space type %d",normal.getFunctionSpace().getTypeCode());
      break;
  }
  checkFinleyError();
}
//
// interpolates data to other domain:
//
void MeshAdapter::interpolateACross(Data& target,const Data& source) const
{
  const MeshAdapter& targetDomain=dynamic_cast<const MeshAdapter&>(target.getFunctionSpace().getDomain());
  if (targetDomain!=*this) 
     throw FinleyAdapterException("Error - Illegal domain of interpolation target");

  Finley_ErrorCode=SYSTEM_ERROR;
  sprintf(Finley_ErrorMsg,"Finley does not allow interpolation across domains yet.");
  checkFinleyError();
}
//
// calculates the integral of a function defined of arg:
//
void MeshAdapter::setToIntegrals(std::vector<double>& integrals,const Data& arg) const
{
  const MeshAdapter& argDomain=dynamic_cast<const MeshAdapter&>(arg.getFunctionSpace().getDomain());
  if (argDomain!=*this) 
     throw FinleyAdapterException("Error - Illegal domain of integration kernel");

  Finley_Mesh* mesh=m_finleyMesh.get();
  switch(arg.getFunctionSpace().getTypeCode()) {
     case(Nodes):
        Finley_ErrorCode=TYPE_ERROR;
        sprintf(Finley_ErrorMsg,"Integral of data on nodes is not supported.");
        break;
     case(Elements):
        Finley_Assemble_integrate(mesh->Nodes,mesh->Elements,&(arg.getDataC()),&integrals[0]);
        break;
     case(FaceElements):
        Finley_Assemble_integrate(mesh->Nodes,mesh->FaceElements,&(arg.getDataC()),&integrals[0]);
        break;
     case(Points):
        Finley_ErrorCode=TYPE_ERROR;
        sprintf(Finley_ErrorMsg,"Integral of data on points is not supported.");
        break;
     case(ContactElementsZero):
        Finley_Assemble_integrate(mesh->Nodes,mesh->ContactElements,&(arg.getDataC()),&integrals[0]);
        break;
     case(ContactElementsOne):
        Finley_Assemble_integrate(mesh->Nodes,mesh->ContactElements,&(arg.getDataC()),&integrals[0]);
        break;
     case(DegreesOfFreedom):
        Finley_ErrorCode=TYPE_ERROR;
        sprintf(Finley_ErrorMsg,"Integral of data on degrees of freedom is not supported.");
        break;
     case(ReducedDegreesOfFreedom):
        Finley_ErrorCode=TYPE_ERROR;
        sprintf(Finley_ErrorMsg,"Integral of data on reduced degrees of freedom is not supported.");
        break;
     default:
        Finley_ErrorCode=TYPE_ERROR;
        sprintf(Finley_ErrorMsg,"Integrals: Finley does not know anything about function space type %d",arg.getFunctionSpace().getTypeCode());
        break;
  }
  checkFinleyError();
}
//
// calculates the gradient of arg:
//
void MeshAdapter::setToGradient(Data& grad,const Data& arg) const
{
  const MeshAdapter& argDomain=dynamic_cast<const MeshAdapter&>(arg.getFunctionSpace().getDomain());
  if (argDomain!=*this)
    throw FinleyAdapterException("Error - Illegal domain of gradient argument");
  const MeshAdapter& gradDomain=dynamic_cast<const MeshAdapter&>(grad.getFunctionSpace().getDomain());
  if (gradDomain!=*this)
     throw FinleyAdapterException("Error - Illegal domain of gradient");

  Finley_Mesh* mesh=m_finleyMesh.get();
  switch(grad.getFunctionSpace().getTypeCode()) {
       case(Nodes):
          Finley_ErrorCode=TYPE_ERROR;
          sprintf(Finley_ErrorMsg,"Gradient at nodes is not supported.");
          break;
       case(Elements):
          Finley_Assemble_gradient(mesh->Nodes,mesh->Elements,&(grad.getDataC()),&(arg.getDataC()));
          break;
       case(FaceElements):
          Finley_Assemble_gradient(mesh->Nodes,mesh->FaceElements,&(grad.getDataC()),&(arg.getDataC()));
          break;
       case(Points):
          Finley_ErrorCode=TYPE_ERROR;
          sprintf(Finley_ErrorMsg,"Gradient at points is not supported.");
          break;
       case(ContactElementsZero):
          Finley_Assemble_gradient(mesh->Nodes,mesh->ContactElements,&(grad.getDataC()),&(arg.getDataC()));
          break;
       case(ContactElementsOne):
          Finley_Assemble_gradient(mesh->Nodes,mesh->ContactElements,&(grad.getDataC()),&(arg.getDataC()));
          break;
       case(DegreesOfFreedom):
          Finley_ErrorCode=TYPE_ERROR;
          sprintf(Finley_ErrorMsg,"Gradient at degrees of freedom is not supported.");
          break;
       case(ReducedDegreesOfFreedom):
          Finley_ErrorCode=TYPE_ERROR;
          sprintf(Finley_ErrorMsg,"Gradient at reduced degrees of freedom is not supported.");
          break;
       default:
          Finley_ErrorCode=TYPE_ERROR;
          sprintf(Finley_ErrorMsg,"Gradient: Finley does not know anything about function space type %d",arg.getFunctionSpace().getTypeCode());
          break;
  }
  checkFinleyError();
}
//
// returns the size of elements:
//
void MeshAdapter::setToSize(Data& size) const
{
  Finley_Mesh* mesh=m_finleyMesh.get();
  escriptDataC tmp=size.getDataC();
  switch(size.getFunctionSpace().getTypeCode()) {
       case(Nodes):
          Finley_ErrorCode=TYPE_ERROR;
          sprintf(Finley_ErrorMsg,"Size of nodes is not supported.");
          break;
       case(Elements):
          Finley_Assemble_getSize(mesh->Nodes,mesh->Elements,&tmp);
          break;
       case(FaceElements):
          Finley_Assemble_getSize(mesh->Nodes,mesh->FaceElements,&tmp);
          break;
       case(Points):
          Finley_ErrorCode=TYPE_ERROR;
          sprintf(Finley_ErrorMsg,"Size of point elements is not supported.");
          break;
       case(ContactElementsZero):
       case(ContactElementsOne):
          Finley_Assemble_getSize(mesh->Nodes,mesh->ContactElements,&tmp);
          break;
       case(DegreesOfFreedom):
          Finley_ErrorCode=TYPE_ERROR;
          sprintf(Finley_ErrorMsg,"Size of degrees of freedom is not supported.");
          break;
       case(ReducedDegreesOfFreedom):
          Finley_ErrorCode=TYPE_ERROR;
          sprintf(Finley_ErrorMsg,"Size of reduced degrees of freedom is not supported.");
          break;
       default:
          Finley_ErrorCode=TYPE_ERROR;
          sprintf(Finley_ErrorMsg,"Element size: Finley does not know anything about function space type %d",size.getFunctionSpace().getTypeCode());
          break;
  }
  checkFinleyError();
}
// sets the location of nodes:
void MeshAdapter::setNewX(const Data& new_x)
{
  Finley_Mesh* mesh=m_finleyMesh.get();
  const MeshAdapter& newDomain=dynamic_cast<const MeshAdapter&>(new_x.getFunctionSpace().getDomain());
  if (newDomain!=*this) 
     throw FinleyAdapterException("Error - Illegal domain of new point locations");
  Finley_NodeFile_setCoordinates(mesh->Nodes,&(new_x.getDataC()));
  checkFinleyError();
}
// saves a data array in openDX format:
void MeshAdapter::saveDX(const std::string& filename,const Data& arg) const
{ 
  Finley_Mesh_saveDX(filename.c_str(),m_finleyMesh.get(),&(arg.getDataC()));
  checkFinleyError();
}
// saves a data array in openVTK format:
void MeshAdapter::saveVTK(const std::string& filename,const Data& arg) const
{ 
  Finley_Mesh_saveVTK(filename.c_str(),m_finleyMesh.get(),&(arg.getDataC()));
  checkFinleyError();
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
    
    Finley_SystemMatrixPattern* fsystemMatrixPattern=Finley_getPattern(getFinley_Mesh(),reduceRowOrder,reduceColOrder);
    checkFinleyError();
    Finley_SystemMatrix* fsystemMatrix=Finley_SystemMatrix_alloc(type,fsystemMatrixPattern,row_blocksize,column_blocksize);
    checkFinleyError();
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
          return true;
          break;
       default:
          Finley_ErrorCode=TYPE_ERROR;
          sprintf(Finley_ErrorMsg,"Cell: Finley does not know anything about function space type %d",functionSpaceCode);
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
           case(ReducedDegreesOfFreedom):
           case(DegreesOfFreedom):
           case(Elements):
           case(FaceElements):
           case(Points):
           case(ContactElementsZero):
           case(ContactElementsOne):
               return true;
           default:
               Finley_ErrorCode=TYPE_ERROR;
               sprintf(Finley_ErrorMsg,"Interpolation On Domain: Finley does not know anything about function space type %d",functionSpaceType_target);
        }
        break;
     case(Elements):
        if (functionSpaceType_target==Elements) {
           return true;
        } else {
           return false;
        }
     case(FaceElements):
        if (functionSpaceType_target==FaceElements) {
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
        } else {
           return false;
        }
     case(DegreesOfFreedom):
        switch(functionSpaceType_target) {
           case(ReducedDegreesOfFreedom):
           case(DegreesOfFreedom):
           case(Nodes):
           case(Elements):
           case(FaceElements):
           case(Points):
           case(ContactElementsZero):
           case(ContactElementsOne):
              return true;
           default:
             Finley_ErrorCode=TYPE_ERROR;
             sprintf(Finley_ErrorMsg,"Interpolation On Domain: Finley does not know anything about function space type %d",functionSpaceType_target);
        }
        break;
     case(ReducedDegreesOfFreedom):
       switch(functionSpaceType_target) {
          case(ReducedDegreesOfFreedom):
          case(Nodes):
          case(Elements):
          case(FaceElements):
          case(Points):
          case(ContactElementsZero):
          case(ContactElementsOne):
              return true;
          case(DegreesOfFreedom):
             return false;
          default:
             Finley_ErrorCode=TYPE_ERROR;
             sprintf(Finley_ErrorMsg,"Interpolation On Domain: Finley does not know anything about function space type %d",functionSpaceType_target);
       }
       break;
     default:
        Finley_ErrorCode=TYPE_ERROR;
        sprintf(Finley_ErrorMsg,"Interpolation On Domain: Finley does not know anything about function space type %d",functionSpaceType_source);
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

int MeshAdapter::getSystemMatrixTypeId(const int solver, const bool symmetry) const
{
   int out=Finley_SystemMatrix_getSystemMatrixTypeId(solver,symmetry?1:0);
   checkFinleyError();
   return out;
}
Data MeshAdapter::getX() const
{
  return continuousFunction(asAbstractContinuousDomain()).getX();
}
Data MeshAdapter::getNormal() const
{
  return functionOnBoundary(asAbstractContinuousDomain()).getNormal();
}
Data MeshAdapter::getSize() const
{
  return function(asAbstractContinuousDomain()).getSize();
}

int MeshAdapter::getTagFromSampleNo(int functionSpaceType, int sampleNo) const
{
  int* tagList;
  int numTags;
  getTagList(functionSpaceType, &tagList, &numTags);
  return tagList[sampleNo];
}

int MeshAdapter::getReferenceNoFromSampleNo(int functionSpaceType, int sampleNo) const
{
  int* referenceNoList;
  int numReferenceNo;
  getReferenceNoList(functionSpaceType, &referenceNoList, &numReferenceNo);
  return referenceNoList[sampleNo];
}

}  // end of namespace
