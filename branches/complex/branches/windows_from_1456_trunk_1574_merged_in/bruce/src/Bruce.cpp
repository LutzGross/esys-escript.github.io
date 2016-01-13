
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

#include "Bruce.h"

#include "BruceException.h"

#include "escript/FunctionSpaceFactory.h"

#include <boost/python/extract.hpp>
#include <vector>
#include <string>
#include "vtkCellType.h"

using namespace std;
using namespace escript;

namespace bruce {

BRUCE_DLL_API
const int Bruce::ContinuousFunction=0;
BRUCE_DLL_API
const int Bruce::Function=1;

Bruce::FunctionSpaceNamesMapType Bruce::m_functionSpaceTypeNames;

Bruce::Bruce():
  m_n0(0), m_n1(0), m_n2(0)
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

  // reorder vectors and point counts according to point counts

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

int
Bruce::getNumSamples(int functionSpaceCode) const
{
  std::pair<int,int> domainShape = getDataShape(functionSpaceCode);
  return domainShape.second;
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

    } else if (dim==1) {

      // Bruce domains in 1d space

      for (int i=0; i<m_n0; i++) {
        int sampleNo=i;
        DataAbstract::ValueType::value_type* sampleData = out.getSampleData(sampleNo);
        sampleData[0] = m_origin[0] + m_v0[0]*i;
      }

    } else if (dim==2) {

      // Bruce domains in 2d space

      for (int i=0; i<m_n0; i++) {
        for (int j=0; j<m_n1; j++) {
          int sampleNo=(m_n1*i)+j;
          DataAbstract::ValueType::value_type* sampleData = out.getSampleData(sampleNo);
          for (int d=0; d<dim; d++) {
            sampleData[d] = m_origin[d] + m_v0[d]*i + m_v1[d]*j;
          }
        }
      }

    } else if (dim==3) {

      // Bruce domains in 3d space

      for (int i=0; i<m_n0; i++) {
        for (int j=0; j<m_n1; j++) {
          for (int k=0; k<m_n2; k++) {
            int sampleNo=(m_n1*m_n2*i)+(m_n2*j)+k;
            DataAbstract::ValueType::value_type* sampleData = out.getSampleData(sampleNo);
            for (int d=0; d<dim; d++) {
              sampleData[d] = m_origin[d] + m_v0[d]*i + m_v1[d]*j + m_v2[d]*k;
            }
          }
        }
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
      if (isZero(m_v0)) {
        stringstream temp;
        temp << "Error - Invalid function space type: "
             << functionSpaceCode << " for Bruce::setToX";
        throw BruceException(temp.str());
      } else {
        for (int i=0; i<n0max; i++) {
          int sampleNo=i;
          DataAbstract::ValueType::value_type* sampleData = out.getSampleData(sampleNo);
          sampleData[0] = m_origin[0] + m_v0[0]*(i + 0.5);
        }
      }

    } else if (dim==2) {

      // Bruce domains in 2d space

      int n0max=m_n0-1;
      int n1max=m_n1-1;
      if (isZero(m_v0)) {
        stringstream temp;
        temp << "Error - Invalid function space type: "
             << functionSpaceCode << " for Bruce::setToX";
        throw BruceException(temp.str());
      } else if (isZero(m_v1)) {
        for (int i=0; i<n0max; i++) {
          int sampleNo=i;
          DataAbstract::ValueType::value_type* sampleData = out.getSampleData(sampleNo);
          for (int d=0; d<dim; d++) {
            sampleData[d] = m_origin[d] + m_v0[d]*(i + 0.5);
          }
        }
      } else {
        for (int i=0; i<n0max; i++) {
          for (int j=0; j<n1max; j++) {
            int sampleNo=(n1max*i)+j;
            DataAbstract::ValueType::value_type* sampleData = out.getSampleData(sampleNo);
            for (int d=0; d<dim; d++) {
              sampleData[d] = m_origin[d] + m_v0[d]*(i + 0.5) + m_v1[d]*(j + 0.5);
            }
          }
        }
      }

    } else if (dim==3) {

      // Bruce domains in 3d space

      int n0max=m_n0-1;
      int n1max=m_n1-1;
      int n2max=m_n2-1;
      if (isZero(m_v0)) {
        stringstream temp;
        temp << "Error - Invalid function space type: "
             << functionSpaceCode << " for Bruce::setToX";
        throw BruceException(temp.str());
      } else if (isZero(m_v1)) {
        for (int i=0; i<n0max; i++) {
          int sampleNo=i;
          DataAbstract::ValueType::value_type* sampleData = out.getSampleData(sampleNo);
          for (int d=0; d<dim; d++) {
            sampleData[d] = m_origin[d] + m_v0[d]*(i + 0.5);
          }
        }
      } else if (isZero(m_v2)) {
        for (int i=0; i<n0max; i++) {
          for (int j=0; j<n1max; j++) {
            int sampleNo=(n1max*i)+j;
            DataAbstract::ValueType::value_type* sampleData = out.getSampleData(sampleNo);
            for (int d=0; d<dim; d++) {
              sampleData[d] = m_origin[d] + m_v0[d]*(i + 0.5) + m_v1[d]*(j + 0.5);
            }
          }
        }
      } else {
        for (int i=0; i<n0max; i++) {
          for (int j=0; j<n1max; j++) {
            for (int k=0; k<n2max; k++) {
              int sampleNo=(n1max*n2max*i)+(n2max*j)+k;
              DataAbstract::ValueType::value_type* sampleData = out.getSampleData(sampleNo);
              for (int d=0; d<dim; d++) {
                sampleData[d] = m_origin[d] + m_v0[d]*(i + 0.5) + m_v1[d]*(j + 0.5) + m_v2[d]*(k + 0.5);
              }
            }
          }
        }
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
  // this check should be satisfied by Data objects passed to setToSize, but
  // FunctionSpace::getSize() seems to create an object which is larger than
  // this... either way, this method can deal with this
/*
  if (dataShape.size()!=1 || dataShape[0]!=1) {
    stringstream temp;
    temp << "Error - Incompatible shape Data object supplied to Bruce::setToSize";
    throw BruceException(temp.str());
  }
*/

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

    } else if (dim==1) {

      // Bruce domains in 1d space

      dp_size = m_v0[0];

      for (int i=0; i<m_n0; i++) {
        int sampleNo=i;
        DataAbstract::ValueType::value_type* sampleData = out.getSampleData(sampleNo);
        sampleData[0] = dp_size;
      }

    } else if (dim==2) {

      // Bruce domains in 2d space

      double x = m_v0[0] + m_v1[0];
      double y = m_v0[1] + m_v1[1];
      dp_size = sqrt(pow(x,2)+pow(y,2));

      for (int i=0; i<m_n0; i++) {
        for (int j=0; j<m_n1; j++) {
          int sampleNo=(m_n1*i)+j;
          DataAbstract::ValueType::value_type* sampleData = out.getSampleData(sampleNo);
          sampleData[0] = dp_size;
        }
      }

    } else if (dim==3) {

      // Bruce domains in 3d space

      double x = m_v0[0] + m_v1[0] + m_v2[0];
      double y = m_v0[1] + m_v1[1] + m_v2[1];
      double z = m_v0[2] + m_v1[2] + m_v2[2];
      dp_size = sqrt(pow(x,2)+pow(y,2)+pow(z,2));

      for (int i=0; i<m_n0; i++) {
        for (int j=0; j<m_n1; j++) {
          for (int k=0; k<m_n2; k++) {
            int sampleNo=(m_n1*m_n2*i)+(m_n2*j)+k;
            DataAbstract::ValueType::value_type* sampleData = out.getSampleData(sampleNo);
            sampleData[0] = dp_size;
          }
        }
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
      if (isZero(m_v0)) {
        stringstream temp;
        temp << "Error - Invalid function space type: "
             << functionSpaceCode << " for Bruce::setToSize";
        throw BruceException(temp.str());
      } else {
        for (int i=0; i<n0max; i++) {
          int sampleNo=i;
          DataAbstract::ValueType::value_type* sampleData = out.getSampleData(sampleNo);
          sampleData[0] = dp_size;
        }
      }

    } else if (dim==2) {

      // Bruce domains in 2d space

      double x = m_v0[0] + m_v1[0];
      double y = m_v0[1] + m_v1[1];
      dp_size = sqrt(pow(x,2)+pow(y,2));

      int n0max=m_n0-1;
      int n1max=m_n1-1;
      if (isZero(m_v0)) {
        stringstream temp;
        temp << "Error - Invalid function space type: "
             << functionSpaceCode << " for Bruce::setToSize";
        throw BruceException(temp.str());
      } else if (isZero(m_v1)) {
        for (int i=0; i<n0max; i++) {
          int sampleNo=i;
          DataAbstract::ValueType::value_type* sampleData = out.getSampleData(sampleNo);
          sampleData[0] = dp_size;
        }
      } else {
        for (int i=0; i<n0max; i++) {
          for (int j=0; j<n1max; j++) {
            int sampleNo=(n1max*i)+j;
            DataAbstract::ValueType::value_type* sampleData = out.getSampleData(sampleNo);
            sampleData[0] = dp_size;
          }
        }
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
      if (isZero(m_v0)) {
        stringstream temp;
        temp << "Error - Invalid function space type: "
             << functionSpaceCode << " for Bruce::setToSize";
        throw BruceException(temp.str());
      } else if (isZero(m_v1)) {
        for (int i=0; i<n0max; i++) {
          int sampleNo=i;
          DataAbstract::ValueType::value_type* sampleData = out.getSampleData(sampleNo);
          sampleData[0] = dp_size;
        }
      } else if (isZero(m_v2)) {
        for (int i=0; i<n0max; i++) {
          for (int j=0; j<n1max; j++) {
            int sampleNo=(n1max*i)+j;
            DataAbstract::ValueType::value_type* sampleData = out.getSampleData(sampleNo);
            sampleData[0] = dp_size;
          }
        }
      } else {
        for (int i=0; i<n0max; i++) {
          for (int j=0; j<n1max; j++) {
            for (int k=0; k<n2max; k++) {
              int sampleNo=(n2max*n1max*i)+(n1max*j)+k;
              DataAbstract::ValueType::value_type* sampleData = out.getSampleData(sampleNo);
              sampleData[0] = dp_size;
            }
          }
        }
      }

    }

  } else {
    stringstream temp;
    temp << "Error - Invalid function space type: "
         << functionSpaceCode << " for domain: " << getDescription();
    throw BruceException(temp.str());
  }

}

void
Bruce::setToGradient(escript::Data& grad,
                     const escript::Data& arg) const
{
  // pre-conditions: arg must be rank 0->3 only, rank 4 cannot be accomodated
  //                 grad will be rank of arg, +1
  // grad is calculated by, for each data-point in this Bruce domain, retrieving
  // the associated value from arg, calculating the grad, and loading this value
  // into the corresponding data-point in grad
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

int
Bruce::getReferenceNoFromSampleNo(int functionSpaceCode,
                                  int sampleNo) const
{
  // ensure supplied sampleNo is valid
  std::pair<int,int> domainShape = getDataShape(functionSpaceCode);
  int numSamples = domainShape.second;
  if ( (sampleNo>=numSamples) || (sampleNo < 0)) {
    stringstream temp;
    temp << "Bruce::getReferenceNoFromSampleNo: Error - invalid sampleNo supplied";
    throw BruceException(temp.str());
  }
  // we just set the reference number to be the sample number for all samples
  return sampleNo;
}

void
Bruce::saveVTK(const std::string& filename,
               const boost::python::dict& dataDict) const
{

  const int num_data=boost::python::extract<int>(dataDict.attr("__len__")());
  boost::python::list keys=dataDict.keys();

  //
  // extract Data objects and associated names from dictionary object
  const std::string::size_type MAX_namelength=256;

  //
  // open archive file
  ofstream archiveFile;
  archiveFile.open(filename.data(), ios::out);
  if (!archiveFile.good()) {
      throw BruceException("Error in Bruce::saveVTK: File could not be opened for writing");
  }
  //
  // determine mesh type to be written
  bool write_celldata=false, write_pointdata=false;
  for (int i_data=0; i_data<num_data; i_data++) {
    Data& d=boost::python::extract<Data&>(dataDict[keys[i_data]]);
    if (!d.isEmpty()) {
      switch(d.getFunctionSpace().getTypeCode()) {
        case ContinuousFunction:
          write_pointdata=true;
          break;
        case Function:
          write_celldata=true;
          break;
      }
    }
  }

  //
  // determine number of points and cells
  int numPoints=getDataShape(ContinuousFunction).second;
  int numCells=getDataShape(Function).second;

  //
  // determine VTK element type
  int cellType;
  std::string elemTypeStr;
  int numVTKNodesPerElement;

  int nDim = getDim();
  switch (nDim) {

    case 0:
      cellType = VTK_VERTEX;
      numVTKNodesPerElement = 1;
      elemTypeStr = "VTK_VERTEX";
      break;

    case 1:
      cellType = VTK_LINE;
      numVTKNodesPerElement = 2;
      elemTypeStr = "VTK_LINE";
      break;

    case 2:
      cellType = VTK_QUAD;
      numVTKNodesPerElement = 4;
      elemTypeStr = "VTK_QUAD";
      break;

    case 3:
      cellType = VTK_HEXAHEDRON;
      numVTKNodesPerElement = 8;
      elemTypeStr = "VTK_HEXAHEDRON";
      break;

  }

  //
  // write xml header
  archiveFile << "<?xml version=\"1.0\"?>" << endl;

  //
  // determine grid extent

  // ??

  //
  // write grid type and extent
  archiveFile << "\t<VTKFile type=\"StructuredGrid\" version=\"0.1\">" << endl;
  archiveFile << "\t\t<StructuredGrid WholeExtent=\"x1 x2 y1 y2 z1 z2\">" << endl;
  archiveFile << "\t\t\t<Piece Extent=\"x1 x2 y1 y2 z1 z2\">" << endl;

  //
  // start to write out point definitions
  archiveFile << "\t\t\t\t<Points>" << endl;

  //
  // determine grid cooordinates

  // ??

  // vtk/mayavi doesn't like 2D data, it likes 3D data with
  // a degenerate third dimension to handle 2D data.
  // So if nDim is less than 3, must pad all empty dimensions,
  // so that the total number of dims is 3.

  archiveFile << "\t\t\t\t\t<DataArray NumberOfComponents=" << max(3,nDim) << " type=\"Float32\" format=\"ascii\">" << endl;
/*
  for (int i=0; i<numPoints; i++) {
    for (int j=0; j<nDim; j++)
      archiveFile << "%e "; //mesh_p->Nodes->Coordinates[INDEX2(j, i, nDim)])
    for (int k=0; !(k>=3-nDim); k++)
      archiveFile << "0. ";
    archiveFile << endl;
  }
*/
  archiveFile << "\t\t\t\t\t</DataArray>" << endl;
  archiveFile << "\t\t\t\t</Points>" << endl;

  //
  // cell data
  if (write_celldata) {
    //
    // mark the active cell-data data arrays
    bool set_scalar=false, set_vector=false, set_tensor=false;
    archiveFile << "\t\t\t\t<CellData";
    for (int i_data=0; i_data<num_data; i_data++) {
      Data& d = boost::python::extract<Data&>(dataDict[keys[i_data]]);
      const std::string &
      full_name = boost::python::extract<std::string>(keys[i_data]);
      std::string name(full_name, 0, MAX_namelength);
      
      if (!d.isEmpty() && d.getFunctionSpace().getTypeCode()==Function) {
        
        switch(d.getDataPointRank()) {

          case 0:
            if (!set_scalar) {
              archiveFile << " Scalars=" << name;
              set_scalar=true;
            }
            break;

          case 1:
            if (!set_vector) {
              archiveFile << " Vectors=" << name;
              set_vector=true;
            }
            break;

          case 2:
            if (!set_tensor) {
              archiveFile << " Tensors=" << name;
              set_tensor=true;
            }
            break;

          default:
            archiveFile.close();
            throw BruceException("saveVTK: VTK can't handle objects with rank greater than 2.");

        }
      }
    }
    archiveFile << ">" << endl;
    //
    // write the cell-data data arrays
    for (int i_data=0; i_data<num_data; i_data++) {
      Data& d = boost::python::extract<Data&>(dataDict[keys[i_data]]);
      const std::string &
      full_name = boost::python::extract<std::string>(keys[i_data]);
      std::string name(full_name, 0, MAX_namelength);
      
      if (!d.isEmpty() && d.getFunctionSpace().getTypeCode()==Function) {
        int numPointsPerSample=1;
        int rank = d.getDataPointRank();
        int nComp = d.getDataPointSize();
        int nCompReqd; // the number of components required by vtk
        int shape=0;

        switch (rank) {

          case 0:
            nCompReqd=1;
            break;

          case 1:
            shape = d.getDataPointShape()[0];
            if (shape>3) {
              archiveFile.close();
              throw BruceException("saveVTK: rank 1 object must have less than 4 components");
            }
            nCompReqd=3;
            break;

          case 2:
            shape=d.getDataPointShape()[0];
            if (shape>3 || shape!=d.getDataPointShape()[1]) {
              archiveFile.close();
              throw BruceException("saveVTK: rank 2 object must have less than 4x4 components and must have a square shape");
            }
            nCompReqd=9;
            break;

        }

        archiveFile << "\t\t\t\t\t<DataArray Name=" << name << " type=\"Float32\" NumberOfComponents=" << nCompReqd << " format=\"ascii\">" << endl;
/*
        //
        // write out the cell data

        // if the number of required components is more than the number
        // of actual components, pad with zeros
        // probably only need to get shape of first element
        // write the data different ways for scalar, vector and tensor

        for (int i=0; i<numCells; i++) {

          double *values=getSampleData(ptr_data[i_data], i);
          double sampleAvg[nComp];

          // average over the number of points in the sample (ie: 1)
          for (int k=0; k<nComp; k++) {
            sampleAvg[k] = values[INDEX2(k,0,nComp)];
          }

          if (nCompReqd == 1) {

            archiveFile << sampleAvg[0] << endl;

          } else if (nCompReqd == 3) {

            for (int m=0; m<shape; m++) {
              archiveFile << sampleAvg[m] << endl;
            }
            for (int m=0; m<nCompReqd-shape; m++) {
              archiveFile << "0." << endl;
            }

          } else if (nCompReqd == 9) {

            int count=0;
            for (int m=0; m<shape; m++) {
              for (int n=0; n<shape; n++) {
                archiveFile << sampleAvg[count] << endl;
                count++;
              }
              for (int n=0; n<3-shape; n++) {
                archiveFile << "0." << endl;
              }
            }
            for (int m=0; m<3-shape; m++) {
              for (int n=0; n<3; n++) {
                archiveFile << "0." << endl;
              }
            }

          }
        }
*/
        archiveFile << "\t\t\t\t\t</DataArray>" << endl;

      }
    }

    archiveFile << "\t\t\t\t</CellData>" << endl;
  }
  //
  // point data
  if (write_pointdata) {
    //
    // mark the active point-data data arrays 
    bool set_scalar=false, set_vector=false, set_tensor=false;
    archiveFile << "\t\t\t\t<PointData";
    for (int i_data=0; i_data<num_data; i_data++) {
      Data& d = boost::python::extract<Data&>(dataDict[keys[i_data]]);
      const std::string &
      full_name = boost::python::extract<std::string>(keys[i_data]);
      std::string name(full_name, 0, MAX_namelength);

      if (!d.isEmpty() && d.getFunctionSpace().getTypeCode()==ContinuousFunction) {

        switch(d.getDataPointRank()) {

          case 0:
            if (!set_scalar) {
              archiveFile << " Scalars=" << name;
              set_scalar=true;
            }
            break;

          case 1:
            if (!set_vector) {
              archiveFile << " Vectors=" << name;
              set_vector=true;
            }
            break;

          case 2:
            if (!set_tensor) {
              archiveFile << " Tensors=" << name;
              set_tensor=true;
            }
            break;

          default:
            archiveFile.close();
            throw BruceException("saveVTK: Vtk can't handle objects with rank greater than 2");
        }

      }
    }
    archiveFile << ">" << endl;
    //
    // write the point-data data arrays
    for (int i_data=0; i_data<num_data; i_data++) {
      Data& d = boost::python::extract<Data&>(dataDict[keys[i_data]]);
      const std::string &
      full_name = boost::python::extract<std::string>(keys[i_data]);
      std::string name(full_name, 0, MAX_namelength);

      if (!d.isEmpty() && 
          d.getFunctionSpace().getTypeCode()==ContinuousFunction) {

        int numPointsPerSample=1;
        int rank = d.getDataPointRank();
        int nComp = d.getDataPointSize();
        int nCompReqd; // the number of components required by vtk
        int shape=0;

        switch (rank) {

        case 0:
          nCompReqd=1;
          break;
        case 1:
          shape=d.getDataPointShape()[0];
          if (shape>3) {
            archiveFile.close();
            throw BruceException("saveVTK: rank 1 object must have less than 4 components");
          }
          nCompReqd=3;
          break;
        case 2:
          shape=d.getDataPointShape()[0];
          if (shape>3 || shape!=d.getDataPointShape()[1]) {
            archiveFile.close();
            throw BruceException("saveVTK: rank 2 object must have less than 4x4 components and must have a square shape");
          }
          nCompReqd=9;
          break;
        }

        archiveFile << "\t\t\t\t\t<DataArray Name=" << name << " type=\"Float32\" NumberOfComponents=" << nCompReqd << " format=\"ascii\">" << endl;
/*
        //
        // write out the point data

        // if the number of required components is more than
        // the number of actual components, pad with zeros
        // write the data different ways for scalar, vector and tensor

        for (int i=0; i<numNodes; i++) {

          values=getSampleData(data[i_data], i);

          if (nCompReqd==1) {

            archiveFile << values[0];

          } else if (nCompReqd==3) {

            for (int m=0; m<shape; m++) {
              archiveFile << values[m];
            }
            for (int m=0; m<nCompReqd-shape; m++) {
              archiveFile << "0.";
            }

          } else if (nCompReqd==9) {

             int count=0;
             for (int m=0; m<shape; m++) {
               for (int n=0; n<shape; n++) {
                 archiveFile << values[count];
                 count++;
               }
               for (int n=0; n<3-shape; n++) {
                 archiveFile << "0.";
               }
             }
             for (int m=0; m<3-shape; m++) {
               for (int n=0; n<3; n++) {
                 archiveFile << "0.";
               }
             }

          }

          archiveFile << endl;
        }
*/
        archiveFile << "\t\t\t\t\t</DataArray>" << endl;
      }
    }
    archiveFile << "\t\t\t\t</PointData>" << endl;
  }
  //
  // finish off the grid definition
  archiveFile << "\t\t\t</Piece>" << endl;
  archiveFile << "\t\t</StructuredGrid>" << endl;

  //
  // write the xml footer
  archiveFile << "\t</VTKFile>" << endl;

  //
  // Close archive file
  archiveFile.close();

  if (!archiveFile.good()) {
    throw BruceException("Error in Bruce::saveVTK: problem closing file");
  }

}

void
Bruce::interpolateOnDomain(escript::Data& target,
                           const escript::Data& source) const
{
}

bool
Bruce::probeInterpolationOnDomain(int functionSpaceType_source,
                                  int functionSpaceType_target) const
{
  return true;
}

void
Bruce::interpolateACross(escript::Data& target,
                         const escript::Data& source) const
{
}

bool
Bruce::probeInterpolationACross(int functionSpaceType_source,
                                const AbstractDomain& targetDomain,
                                int functionSpaceType_target) const
{
  return true;
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
