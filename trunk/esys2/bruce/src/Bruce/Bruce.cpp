// $Id$
/*
 ******************************************************************************
 *                                                                            *
 *       COPYRIGHT  ACcESS 2005 -  All Rights Reserved                        *
 *                                                                            *
 * This software is the property of ACcESS. No part of this code              *
 * may be copied in any form or by any means without the expressed written    *
 * consent of ACcESS.  Copying, use or modification of this software          *
 * by any unauthorised person is illegal unless that person has a software    *
 * license agreement with ACcESS.                                             *
 *                                                                            *
 ******************************************************************************
*/

#include "bruce/Bruce/Bruce.h"
#include "bruce/Bruce/BruceException.h"

using namespace std;
using namespace escript;

namespace bruce {

const int Bruce::ContinuousFunction=0;
const int Bruce::Function=1;

Bruce::FunctionSpaceNamesMapType Bruce::m_functionSpaceTypeNames;

Bruce::Bruce()
{
  setFunctionSpaceTypeNames();
}

Bruce::Bruce(DimVec v0, DimVec v1, DimVec v2,
             int n0, int n1, int n2,
             DimVec origin):
  m_v0(v0), m_v1(v1), m_v2(v2),
  m_n0(n0), m_n1(n1), m_n2(n2),
  m_origin(origin)
{
  if (!checkParameters()) {
    stringstream temp;
    temp << "Error - Invalid parameters supplied to constructor.";
    throw BruceException(temp.str());
  }
  setFunctionSpaceTypeNames();
}

Bruce::Bruce(const Bruce& other):
  m_v0(other.m_v0), m_v1(other.m_v1), m_v2(other.m_v2),
  m_n0(other.m_n0), m_n1(other.m_n1), m_n2(other.m_n2),
  m_origin(other.m_origin)
{
  setFunctionSpaceTypeNames();
}

Bruce::~Bruce()
{
  m_n0=-1;
  m_n1=-1;
  m_n2=-1;
  m_origin.clear();
  m_v0.clear();
  m_v1.clear();
  m_v2.clear();
}

void
Bruce::setFunctionSpaceTypeNames()
{
  m_functionSpaceTypeNames.insert
    (FunctionSpaceNamesMapType::value_type(Function,"Bruce_Function"));
  m_functionSpaceTypeNames.insert
    (FunctionSpaceNamesMapType::value_type(ContinuousFunction,"Bruce_ContinuousFunction"));
}

bool
Bruce::checkParameters()
{
  if (m_origin.size()>3) {
    return false;
  }

  if (m_n0<1) {
    m_n0=1;
  }
  if (m_n1<1) {
    m_n1=1;
  }
  if (m_n2<1) {
    m_n2=1;
  }

  //
  // domains in 3d space
  if (m_origin.size()==3) {

    if (m_v0.size()==0) {

      // 0d domain in 3d space
      if ( (m_v1.size()!=0) || (m_v2.size()!=0) ) {
        return false;
      }

      m_v0.clear();
      m_v0.push_back(0);
      m_v0.push_back(0);
      m_v0.push_back(0);
      m_v1.clear();
      m_v1.push_back(0);
      m_v1.push_back(0);
      m_v1.push_back(0);
      m_v2.clear();
      m_v2.push_back(0);
      m_v2.push_back(0);
      m_v2.push_back(0);

      m_n0=1;
      m_n1=1;
      m_n2=1;

    } else if (m_v1.size()==0) {

      // 1d domain in 3d space
      if ( (m_v0.size()!=3) || (m_v2.size()!=0) ) {
        return false;
      }

      m_v1.clear();
      m_v1.push_back(0);
      m_v1.push_back(0);
      m_v1.push_back(0);
      m_v2.clear();
      m_v2.push_back(0);
      m_v2.push_back(0);
      m_v2.push_back(0);

      m_n1=1;
      m_n2=1;

    } else if (m_v2.size()==0) {

      // 2d domain in 3d space
      if ( (m_v0.size()!=3) || (m_v1.size()!=3) ) {
        return false;
      }

      m_v2.clear();
      m_v2.push_back(0);
      m_v2.push_back(0);
      m_v2.push_back(0);

      m_n2=1;

    } else {

      // 3d domain in 3d space
      if ( (m_v0.size()!=3) || (m_v1.size()!=3) || (m_v2.size()!=3) ) {
        return false;
      }

    }

  }

  //
  // domains in 2d space
  if (m_origin.size()==2) {

    if (m_v0.size()==0) {

      // 0d domain in 2d space
      if (m_v1.size()!=0) {
        return false;
      }

      m_v0.clear();
      m_v0.push_back(0);
      m_v0.push_back(0);
      m_v1.clear();
      m_v1.push_back(0);
      m_v1.push_back(0);
      m_v2.clear();

      m_n0=1;
      m_n1=1;
      m_n2=1;

    } else if (m_v1.size()==0) {

      // 1d domain in 2d space
      if (m_v0.size()!=2) {
        return false;
      }

      m_v1.clear();
      m_v1.push_back(0);
      m_v1.push_back(0);
      m_v2.clear();

      m_n1=1;
      m_n2=1;

    } else {

      // 2d domain in 2d space
      if ( (m_v0.size()!=2) || (m_v1.size()!=2) ) {
        return false;
      }

      m_v2.clear();

      m_n2=1;

    }

  }

  //
  // domains in 1d space
  if (m_origin.size()==1) {

    if (m_v0.size()==0) {

      // 0d domain in 1d space
      m_v0.clear();
      m_v0.push_back(0);
      m_v1.clear();
      m_v2.clear();

      m_n0=1;
      m_n1=1;
      m_n2=1;

    } else {

      // 1d domain in 1d space
      if (m_v0.size()!=1)  {
        return false;
      }

      m_v1.clear();
      m_v2.clear();

      m_n1=1;
      m_n2=1;

    }

  }

  //
  // domains in 0d space
  if (m_origin.size()==0) {

    // point (0d) domain in 0d space
    m_v0.clear();
    m_v1.clear();
    m_v2.clear();

    m_n0=1;
    m_n1=1;
    m_n2=1;

  }

  return true;
}

bool
Bruce::isValidFunctionSpaceType(int functionSpaceCode) const
{
  FunctionSpaceNamesMapType::iterator loc;
  loc=m_functionSpaceTypeNames.find(functionSpaceCode);
  return (loc!=m_functionSpaceTypeNames.end());
}

std::string
Bruce::functionSpaceTypeAsString(int functionSpaceCode) const
{
  FunctionSpaceNamesMapType::iterator loc;
  loc=m_functionSpaceTypeNames.find(functionSpaceCode);
  if (loc==m_functionSpaceTypeNames.end()) {
    stringstream temp;
    temp << "Error - Invalid function space type code.";
    throw BruceException(temp.str());
  } else {
    return loc->second;
  }
}

pair<int,int>
Bruce::getDataShape(int functionSpaceCode) const
{
  int numDataPointsPerSample=1;
  int numSamples;

  switch (functionSpaceCode) {

      //
      // Continuous functions
      case(ContinuousFunction):

          numSamples = m_n0 * m_n1 * m_n2;
          break;

      //
      // Functions
      case(Function):

          // 0d domains
          if (getDim()==0) {

            numSamples = 0;

          // 1d domains
          } else if (getDim()==1) {

            if (isZero(m_v0)) {
              numSamples = 0;
            } else {
              numSamples = m_n0-1;
            }

          // 2d domains
          } else if (getDim()==2) {

            if (isZero(m_v0)) {
              numSamples = 0;
            } else if (isZero(m_v1)) {
              numSamples = m_n0-1;
            } else {
              numSamples = (m_n0-1) * (m_n1-1);
            }

          // 3d domains
          } else {

            if (isZero(m_v0)) {
              numSamples = 0;
            } else if (isZero(m_v1)) {
              numSamples = m_n0-1;
            } else if (isZero(m_v2)) {
              numSamples = (m_n0-1) * (m_n1-1);
            } else {
              numSamples = (m_n0-1) * (m_n1-1) * (m_n2-1);
            }

          }

          break;

      default:
          stringstream temp;
          temp << "Error - Invalid function space type: "
               << functionSpaceCode << " for domain: " << getDescription();
          throw BruceException(temp.str());
          break;

  }

  return pair<int,int>(numDataPointsPerSample,numSamples);
}

Data
Bruce::getX() const
{
  FunctionSpace tempFunc = continuousFunction(asAbstractContinuousDomain());
  return tempFunc.getX();
}

void
Bruce::setToX(escript::Data& out) const
{

  //
  // determine functionSpace type expected by supplied Data object
  int functionSpaceCode = out.getFunctionSpace().getTypeCode();

  //
  // ensure supplied Data object has a matching number of data-points
  // for this Bruce domain object
  std::pair<int,int> domainShape = getDataShape(functionSpaceCode);
  if(domainShape.first!=out.getNumDataPointsPerSample() ||
     domainShape.second!=out.getNumSamples()) {
    stringstream temp;
    temp << "Error - Incompatible number of data-points Data object supplied to Bruce::setToX";
    throw BruceException(temp.str());
  }

  int dim = getDim();
  int numSamples = domainShape.second;

  //
  // ensure shape of data-points in supplied Data object matches the
  // shape needed to store the coordinates of each data-point in this
  // Bruce domain
  std::vector<int> dataShape = out.getDataPointShape();
  if (dim>0 && (dataShape.size()!=1 || dataShape[0]!=dim)) {
    stringstream temp;
    temp << "Error - Incompatible shape Data object supplied to Bruce::setToX";
    throw BruceException(temp.str());
  } else if (dim==0 && dataShape.size()!=0) {
    stringstream temp;
    temp << "Error - Incompatible shape Data object supplied to Bruce::setToX";
    throw BruceException(temp.str());
  }

  if (functionSpaceCode==ContinuousFunction) {

    //
    // calculate the spatial coordinates of data-points 
    // located on the nodes of this Bruce domain

    if (dim==0) {

      // Bruce domains in 0d space

      int sampleNo=0;
      DataAbstract::ValueType::value_type* sampleData = out.getSampleData(sampleNo);
      sampleNo++;
      if (sampleNo!=numSamples) {
        stringstream temp;
        temp << "Bruce::setToX: Didn't iterate across correct number of samples.";
        throw BruceException(temp.str());
      }

    } else if (dim==1) {

      // Bruce domains in 1d space

      int sampleNo=0;
      for (int i=0; i<m_n0; i++) {
        DataAbstract::ValueType::value_type* sampleData = out.getSampleData(sampleNo);
        sampleData[0] = m_origin[0] + m_v0[0]*i;
        sampleNo++;
      }
      if (sampleNo!=numSamples) {
        stringstream temp;
        temp << "Bruce::setToX: Didn't iterate across correct number of samples.";
        throw BruceException(temp.str());
      }

    } else if (dim==2) {

      // Bruce domains in 2d space

      int sampleNo=0;
      for (int i=0; i<m_n0; i++) {
        for (int j=0; j<m_n1; j++) {
          DataAbstract::ValueType::value_type* sampleData = out.getSampleData(sampleNo);
          for (int d=0; d<dim; d++) {
            sampleData[d] = m_origin[d] + m_v0[d]*i + m_v1[d]*j;
          }
          sampleNo++;
        }
      }
      if (sampleNo!=numSamples) {
        stringstream temp;
        temp << "Bruce::setToX: Didn't iterate across correct number of samples.";
        throw BruceException(temp.str());
      }

    } else if (dim==3) {

      // Bruce domains in 3d space

      int sampleNo=0;
      for (int i=0; i<m_n0; i++) {
        for (int j=0; j<m_n1; j++) {
          for (int k=0; k<m_n2; k++) {
            DataAbstract::ValueType::value_type* sampleData = out.getSampleData(sampleNo);
            for (int d=0; d<dim; d++) {
              sampleData[d] = m_origin[d] + m_v0[d]*i + m_v1[d]*j + m_v2[d]*k;
            }
            sampleNo++;
          }
        }
      }
      if (sampleNo!=numSamples) {
        stringstream temp;
        temp << "Bruce::setToX: Didn't iterate across correct number of samples.";
        throw BruceException(temp.str());
      }

    }

  } else if (functionSpaceCode==Function) {

    //
    // calculate the spatial coordinates of data-points 
    // located on the cell centres of this Bruce domain

    if (dim==0) {

      // Bruce domains in 0d space

      stringstream temp;
      temp << "Error - Invalid function space type: "
           << functionSpaceCode << " for Bruce::setToX";
      throw BruceException(temp.str());

    } else if (dim==1) {

      // Bruce domains in 1d space

      int n0max=m_n0-1;
      int sampleNo=0;
      if (isZero(m_v0)) {
        stringstream temp;
        temp << "Error - Invalid function space type: "
             << functionSpaceCode << " for Bruce::setToX";
        throw BruceException(temp.str());
      } else {
        for (int i=0; i<n0max; i++) {
          DataAbstract::ValueType::value_type* sampleData = out.getSampleData(sampleNo);
          sampleData[0] = m_origin[0] + m_v0[0]*(i + 0.5);
          sampleNo++;
        }
      }
      if (sampleNo!=numSamples) {
        stringstream temp;
        temp << "Bruce::setToX: Didn't iterate across correct number of samples.";
        throw BruceException(temp.str());
      }

    } else if (dim==2) {

      // Bruce domains in 2d space

      int n0max=m_n0-1;
      int n1max=m_n1-1;
      int sampleNo=0;
      if (isZero(m_v0)) {
        stringstream temp;
        temp << "Error - Invalid function space type: "
             << functionSpaceCode << " for Bruce::setToX";
        throw BruceException(temp.str());
      } else if (isZero(m_v1)) {
        for (int i=0; i<n0max; i++) {
          DataAbstract::ValueType::value_type* sampleData = out.getSampleData(sampleNo);
          for (int d=0; d<dim; d++) {
            sampleData[d] = m_origin[d] + m_v0[d]*(i + 0.5);
          }
          sampleNo++;
        }
      } else {
        for (int i=0; i<n0max; i++) {
          for (int j=0; j<n1max; j++) {
            DataAbstract::ValueType::value_type* sampleData = out.getSampleData(sampleNo);
            for (int d=0; d<dim; d++) {
              sampleData[d] = m_origin[d] + m_v0[d]*(i + 0.5) + m_v1[d]*(j + 0.5);
            }
            sampleNo++;
          }
        }
      }
      if (sampleNo!=numSamples) {
        stringstream temp;
        temp << "Bruce::setToX: Didn't iterate across correct number of samples.";
        throw BruceException(temp.str());
      }

    } else if (dim==3) {

      // Bruce domains in 3d space

      int n0max=m_n0-1;
      int n1max=m_n1-1;
      int n2max=m_n2-1;
      int sampleNo=0;
      if (isZero(m_v0)) {
        stringstream temp;
        temp << "Error - Invalid function space type: "
             << functionSpaceCode << " for Bruce::setToX";
        throw BruceException(temp.str());
      } else if (isZero(m_v1)) {
        for (int i=0; i<n0max; i++) {
          DataAbstract::ValueType::value_type* sampleData = out.getSampleData(sampleNo);
          for (int d=0; d<dim; d++) {
            sampleData[d] = m_origin[d] + m_v0[d]*(i + 0.5);
          }
          sampleNo++;
        }
      } else if (isZero(m_v2)) {
        for (int i=0; i<n0max; i++) {
          for (int j=0; j<n1max; j++) {
            DataAbstract::ValueType::value_type* sampleData = out.getSampleData(sampleNo);
            for (int d=0; d<dim; d++) {
              sampleData[d] = m_origin[d] + m_v0[d]*(i + 0.5) + m_v1[d]*(j + 0.5);
            }
            sampleNo++;
          }
        }
      } else {
        for (int i=0; i<n0max; i++) {
          for (int j=0; j<n1max; j++) {
            for (int k=0; k<n2max; k++) {
              DataAbstract::ValueType::value_type* sampleData = out.getSampleData(sampleNo);
              for (int d=0; d<dim; d++) {
                sampleData[d] = m_origin[d] + m_v0[d]*(i + 0.5) + m_v1[d]*(j + 0.5) + m_v2[d]*(k + 0.5);
              }
              sampleNo++;
            }
          }
        }
      }
      if (sampleNo!=numSamples) {
        stringstream temp;
        temp << "Bruce::setToX: Didn't iterate across correct number of samples.";
        throw BruceException(temp.str());
      }

    }

  } else {
    stringstream temp;
    temp << "Error - Invalid function space type: "
         << functionSpaceCode << " for domain: " << getDescription();
    throw BruceException(temp.str());
  }

}

Data
Bruce::getSize() const
{
  FunctionSpace tempFunc = function(asAbstractContinuousDomain());
  return tempFunc.getSize();
}

void
Bruce::setToSize(escript::Data& out) const
{

  //
  // determine functionSpace type expected by supplied Data object
  int functionSpaceCode = out.getFunctionSpace().getTypeCode();

  //
  // ensure supplied Data object has a matching number of data-points
  // for this Bruce domain object
  std::pair<int,int> domainShape = getDataShape(functionSpaceCode);
  if(domainShape.first!=out.getNumDataPointsPerSample() ||
     domainShape.second!=out.getNumSamples()) {
    stringstream temp;
    temp << "Error - Incompatible number of data-points Data object supplied to Bruce::setToSize";
    throw BruceException(temp.str());
  }

  int dim = getDim();
  int numSamples = domainShape.second;

  //
  // ensure shape of data-points in supplied Data object matches the
  // shape needed to store the size of each data-point in this Bruce domain
  std::vector<int> dataShape = out.getDataPointShape();
  if (dataShape.size()!=1 || dataShape[0]!=1) {
    stringstream temp;
    temp << "Error - Incompatible shape Data object supplied to Bruce::setToSize";
    throw BruceException(temp.str());
  }

  double dp_size;

  if (functionSpaceCode==ContinuousFunction) {

    //
    // calculate the size of data-points 
    // located on the nodes of this Bruce domain

    if (dim==0) {

      // Bruce domains in 0d space

      dp_size = 0.0;

      int sampleNo=0;
      DataAbstract::ValueType::value_type* sampleData = out.getSampleData(sampleNo);
      sampleData[0] = dp_size;
      sampleNo++;
      if (sampleNo!=numSamples) {
        stringstream temp;
        temp << "Bruce::setToSize: Didn't iterate across correct number of samples.";
        throw BruceException(temp.str());
      }

    } else if (dim==1) {

      // Bruce domains in 1d space

      dp_size = m_v0[0];

      int sampleNo=0;
      for (int i=0; i<m_n0; i++) {
        DataAbstract::ValueType::value_type* sampleData = out.getSampleData(sampleNo);
        sampleData[0] = dp_size;
        sampleNo++;
      }
      if (sampleNo!=numSamples) {
        stringstream temp;
        temp << "Bruce::setToSize: Didn't iterate across correct number of samples.";
        throw BruceException(temp.str());
      }

    } else if (dim==2) {

      // Bruce domains in 2d space

      double x = m_v0[0] + m_v1[0];
      double y = m_v0[1] + m_v1[1];
      dp_size = sqrt(pow(x,2)+pow(y,2));

      int sampleNo=0;
      for (int i=0; i<m_n0; i++) {
        for (int j=0; j<m_n1; j++) {
          DataAbstract::ValueType::value_type* sampleData = out.getSampleData(sampleNo);
          sampleData[0] = dp_size;
          sampleNo++;
        }
      }
      if (sampleNo!=numSamples) {
        stringstream temp;
        temp << "Bruce::setToSize: Didn't iterate across correct number of samples.";
        throw BruceException(temp.str());
      }

    } else if (dim==3) {

      // Bruce domains in 3d space

      double x = m_v0[0] + m_v1[0] + m_v2[0];
      double y = m_v0[1] + m_v1[1] + m_v2[1];
      double z = m_v0[2] + m_v1[2] + m_v2[2];
      dp_size = sqrt(pow(x,2)+pow(y,2)+pow(z,2));

      int sampleNo=0;
      for (int i=0; i<m_n0; i++) {
        for (int j=0; j<m_n1; j++) {
          for (int k=0; k<m_n2; k++) {
            DataAbstract::ValueType::value_type* sampleData = out.getSampleData(sampleNo);
            sampleData[0] = dp_size;
            sampleNo++;
          }
        }
      }
      if (sampleNo!=numSamples) {
        stringstream temp;
        temp << "Bruce::setToSize: Didn't iterate across correct number of samples.";
        throw BruceException(temp.str());
      }

    }

  } else if (functionSpaceCode==Function) {

    //
    // calculate the sizes of data-points 
    // located on the cell centres of this Bruce domain

    if (dim==0) {

      // Bruce domains in 0d space

      stringstream temp;
      temp << "Error - Invalid function space type: "
           << functionSpaceCode << " for Bruce::setToSize";
      throw BruceException(temp.str());

    } else if (dim==1) {

      // Bruce domains in 1d space

      dp_size = m_v0[0];

      int n0max=m_n0-1;
      int sampleNo=0;
      if (isZero(m_v0)) {
        stringstream temp;
        temp << "Error - Invalid function space type: "
             << functionSpaceCode << " for Bruce::setToSize";
        throw BruceException(temp.str());
      } else {
        for (int i=0; i<n0max; i++) {
          DataAbstract::ValueType::value_type* sampleData = out.getSampleData(sampleNo);
          sampleData[0] = dp_size;
          sampleNo++;
        }
      }
      if (sampleNo!=numSamples) {
        stringstream temp;
        temp << "Bruce::setToSize: Didn't iterate across correct number of samples.";
        throw BruceException(temp.str());
      }

    } else if (dim==2) {

      // Bruce domains in 2d space

      double x = m_v0[0] + m_v1[0];
      double y = m_v0[1] + m_v1[1];
      dp_size = sqrt(pow(x,2)+pow(y,2));

      int n0max=m_n0-1;
      int n1max=m_n1-1;
      int sampleNo=0;
      if (isZero(m_v0)) {
        stringstream temp;
        temp << "Error - Invalid function space type: "
             << functionSpaceCode << " for Bruce::setToSize";
        throw BruceException(temp.str());
      } else if (isZero(m_v1)) {
        for (int i=0; i<n0max; i++) {
          DataAbstract::ValueType::value_type* sampleData = out.getSampleData(sampleNo);
          sampleData[0] = dp_size;
          sampleNo++;
        }
      } else {
        for (int i=0; i<n0max; i++) {
          for (int j=0; j<n1max; j++) {
            DataAbstract::ValueType::value_type* sampleData = out.getSampleData(sampleNo);
            sampleData[0] = dp_size;
            sampleNo++;
          }
        }
      }
      if (sampleNo!=numSamples) {
        stringstream temp;
        temp << "Bruce::setToSize: Didn't iterate across correct number of samples.";
        throw BruceException(temp.str());
      }

    } else if (dim==3) {

      // Bruce domains in 3d space

      double x = m_v0[0] + m_v1[0] + m_v2[0];
      double y = m_v0[1] + m_v1[1] + m_v2[1];
      double z = m_v0[2] + m_v1[2] + m_v2[2];
      dp_size = sqrt(pow(x,2)+pow(y,2)+pow(z,2));

      int n0max=m_n0-1;
      int n1max=m_n1-1;
      int n2max=m_n2-1;
      int sampleNo=0;
      if (isZero(m_v0)) {
        stringstream temp;
        temp << "Error - Invalid function space type: "
             << functionSpaceCode << " for Bruce::setToSize";
        throw BruceException(temp.str());
      } else if (isZero(m_v1)) {
        for (int i=0; i<n0max; i++) {
          DataAbstract::ValueType::value_type* sampleData = out.getSampleData(sampleNo);
          sampleData[0] = dp_size;
          sampleNo++;
        }
      } else if (isZero(m_v2)) {
        for (int i=0; i<n0max; i++) {
          for (int j=0; j<n1max; j++) {
            DataAbstract::ValueType::value_type* sampleData = out.getSampleData(sampleNo);
            sampleData[0] = dp_size;
            sampleNo++;
          }
        }
      } else {
        for (int i=0; i<n0max; i++) {
          for (int j=0; j<n1max; j++) {
            for (int k=0; k<n2max; k++) {
              DataAbstract::ValueType::value_type* sampleData = out.getSampleData(sampleNo);
              sampleData[0] = dp_size;
              sampleNo++;
            }
          }
        }
      }
      if (sampleNo!=numSamples) {
        stringstream temp;
        temp << "Bruce::setToSize: Didn't iterate across correct number of samples.";
        throw BruceException(temp.str());
      }

    }

  } else {
    stringstream temp;
    temp << "Error - Invalid function space type: "
         << functionSpaceCode << " for domain: " << getDescription();
    throw BruceException(temp.str());
  }

}

bool
Bruce::operator==(const AbstractDomain& other) const
{
  const Bruce* temp=dynamic_cast<const Bruce*>(&other);
  if (temp!=0) {
    if ((m_v0 != temp->m_v0) || (m_v1 != temp->m_v1) || (m_v2 != temp->m_v2)) {
      return false;
    }
    if ((m_n0 != temp->m_n0) || (m_n1 != temp->m_n1) || (m_n2 != temp->m_n2)) {
      return false;
    }
    if (m_origin != temp->m_origin) {
      return false;
    }
    return true;
  } else {
    return false;
  }
}

bool
Bruce::operator!=(const AbstractDomain& other) const
{
  return !(operator==(other));
}

bool
Bruce::isZero(DimVec vec)
{
  for (int i=0; i<vec.size(); i++) {
    if (vec[i] != 0) {
      return false;
    }
  }
  return true;
}

}  // end of namespace
