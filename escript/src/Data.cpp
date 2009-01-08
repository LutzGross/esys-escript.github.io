
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#include "Data.h"

#include "DataExpanded.h"
#include "DataConstant.h"
#include "DataTagged.h"
#include "DataEmpty.h"
#include "DataLazy.h"
#include "FunctionSpaceFactory.h"
#include "AbstractContinuousDomain.h"
#include "UnaryFuncs.h"
#include "FunctionSpaceException.h"
#include "EscriptParams.h"

extern "C" {
#include "esysUtils/blocktimer.h"
}

#include <fstream>
#include <algorithm>
#include <vector>
#include <functional>
#include <sstream>	// so we can throw messages about ranks

#include <boost/python/dict.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/long.hpp>

using namespace std;
using namespace boost::python;
using namespace boost;
using namespace escript;

// ensure the current object is not a DataLazy
// The idea was that we could add an optional warning whenever a resolve is forced
#define FORCERESOLVE if (isLazy()) {resolve();}
#define AUTOLAZYON escriptParams.getAUTOLAZY()
#define MAKELAZYOP(X)   if (isLazy() || (AUTOLAZYON && m_data->isExpanded())) \
  {\
	DataLazy* c=new DataLazy(borrowDataPtr(),X);\
	return Data(c);\
  }
#define MAKELAZYOPOFF(X,Y) if (isLazy() || (AUTOLAZYON && m_data->isExpanded())) \
  {\
	DataLazy* c=new DataLazy(borrowDataPtr(),X,Y);\
	return Data(c);\
  }

#define MAKELAZYBINSELF(R,X)   if (isLazy() || R.isLazy() || (AUTOLAZYON && (isExpanded() || R.isExpanded()))) \
  {\
	DataLazy* c=new DataLazy(m_data,R.borrowDataPtr(),X);\
        m_data=c->getPtr();\
	return (*this);\
  }

// like the above but returns a new data rather than *this
#define MAKELAZYBIN(R,X)   if (isLazy() || R.isLazy() || (AUTOLAZYON && (isExpanded() || R.isExpanded()))) \
  {\
	DataLazy* c=new DataLazy(m_data,R.borrowDataPtr(),X);\
	return Data(c);\
  }

#define MAKELAZYBIN2(L,R,X)   if (L.isLazy() || R.isLazy() || (AUTOLAZYON && (L.isExpanded() || R.isExpanded()))) \
  {\
	DataLazy* c=new DataLazy(L.borrowDataPtr(),R.borrowDataPtr(),X);\
	return Data(c);\
  }

Data::Data()
{
  //
  // Default data is type DataEmpty
  DataAbstract* temp=new DataEmpty();
  m_data=temp->getPtr();
  m_protected=false;
}

Data::Data(double value,
           const tuple& shape,
           const FunctionSpace& what,
           bool expanded)
{
  DataTypes::ShapeType dataPointShape;
  for (int i = 0; i < shape.attr("__len__")(); ++i) {
    dataPointShape.push_back(extract<const int>(shape[i]));
  }

  int len = DataTypes::noValues(dataPointShape);
  DataVector temp_data(len,value,len);
  initialise(temp_data, dataPointShape, what, expanded);
  m_protected=false;
}

Data::Data(double value,
	   const DataTypes::ShapeType& dataPointShape,
	   const FunctionSpace& what,
           bool expanded)
{
  int len = DataTypes::noValues(dataPointShape);

  DataVector temp_data(len,value,len);
//   DataArrayView temp_dataView(temp_data, dataPointShape);

//   initialise(temp_dataView, what, expanded);
  initialise(temp_data, dataPointShape, what, expanded);

  m_protected=false;
}

Data::Data(const Data& inData)
{
  m_data=inData.m_data;
  m_protected=inData.isProtected();
}


Data::Data(const Data& inData,
           const DataTypes::RegionType& region)
{
  DataAbstract_ptr dat=inData.m_data;
  if (inData.isLazy())
  {
	dat=inData.m_data->resolve();
  }
  else
  {
	dat=inData.m_data;
  }
  //
  // Create Data which is a slice of another Data
  DataAbstract* tmp = dat->getSlice(region);
  m_data=DataAbstract_ptr(tmp);
  m_protected=false;
}

Data::Data(const Data& inData,
           const FunctionSpace& functionspace)
{
  if (inData.isEmpty())
  {
    throw DataException("Error - will not interpolate for instances of DataEmpty.");
  }
  if (inData.getFunctionSpace()==functionspace) {
    m_data=inData.m_data;
  } 
  else 
  {

    if (inData.isConstant()) {	// for a constant function, we just need to use the new function space
      if (!inData.probeInterpolation(functionspace))
      {           // Even though this is constant, we still need to check whether interpolation is allowed
	throw FunctionSpaceException("Cannot interpolate across to the domain of the specified FunctionSpace. (DataConstant)");
      }
      // if the data is not lazy, this will just be a cast to DataReady
      DataReady_ptr dr=inData.m_data->resolve();
      DataConstant* dc=new DataConstant(functionspace,inData.m_data->getShape(),dr->getVector());	
      m_data=DataAbstract_ptr(dc); 
    } else {
      Data tmp(0,inData.getDataPointShape(),functionspace,true);
      // Note: Must use a reference or pointer to a derived object
      // in order to get polymorphic behaviour. Shouldn't really
      // be able to create an instance of AbstractDomain but that was done
      // as a boost:python work around which may no longer be required.
      /*const AbstractDomain& inDataDomain=inData.getDomain();*/
      const_Domain_ptr inDataDomain=inData.getDomain();
      if  (inDataDomain==functionspace.getDomain()) {
        inDataDomain->interpolateOnDomain(tmp,inData);
      } else {
        inDataDomain->interpolateACross(tmp,inData);
      }
      m_data=tmp.m_data;
    }
  }
  m_protected=false;
}

Data::Data(DataAbstract* underlyingdata)
{
// 	m_data=shared_ptr<DataAbstract>(underlyingdata);
	m_data=underlyingdata->getPtr();
	m_protected=false;
}

Data::Data(DataAbstract_ptr underlyingdata)
{
	m_data=underlyingdata;
	m_protected=false;
}


Data::Data(const numeric::array& value,
	   const FunctionSpace& what,
           bool expanded)
{
  initialise(value,what,expanded);
  m_protected=false;
}

Data::Data(const DataTypes::ValueType& value,
		 const DataTypes::ShapeType& shape,
                 const FunctionSpace& what,
                 bool expanded)
{
   initialise(value,shape,what,expanded);
   m_protected=false;
}


Data::Data(const object& value,
	   const FunctionSpace& what,
           bool expanded)
{
  numeric::array asNumArray(value);
  initialise(asNumArray,what,expanded);
  m_protected=false;
}


Data::Data(const object& value,
           const Data& other)
{
  numeric::array asNumArray(value);

  // extract the shape of the numarray
  DataTypes::ShapeType tempShape=DataTypes::shapeFromNumArray(asNumArray);
  if (DataTypes::getRank(tempShape)==0) {


    // get the space for the data vector
    int len1 = DataTypes::noValues(tempShape);
    DataVector temp_data(len1, 0.0, len1);
    temp_data.copyFromNumArray(asNumArray,1);

    int len = DataTypes::noValues(other.getDataPointShape());

    DataVector temp2_data(len, temp_data[0]/*temp_dataView()*/, len);
    DataConstant* t=new DataConstant(other.getFunctionSpace(),other.getDataPointShape(),temp2_data);
    m_data=DataAbstract_ptr(t);

  } else {
    //
    // Create a DataConstant with the same sample shape as other
    DataConstant* t=new DataConstant(asNumArray,other.getFunctionSpace());
//     boost::shared_ptr<DataAbstract> sp(t);
//     m_data=sp;
    m_data=DataAbstract_ptr(t);
  }
  m_protected=false;
}

Data::~Data()
{

}



void
Data::initialise(const boost::python::numeric::array& value,
                 const FunctionSpace& what,
                 bool expanded)
{
  //
  // Construct a Data object of the appropriate type.
  // Construct the object first as there seems to be a bug which causes
  // undefined behaviour if an exception is thrown during construction
  // within the shared_ptr constructor.
  if (expanded) {
    DataAbstract* temp=new DataExpanded(value, what);
//     boost::shared_ptr<DataAbstract> temp_data(temp);
//     m_data=temp_data;
    m_data=temp->getPtr();
  } else {
    DataAbstract* temp=new DataConstant(value, what);
//     boost::shared_ptr<DataAbstract> temp_data(temp);
//     m_data=temp_data;
    m_data=temp->getPtr();
  }
}


void
Data::initialise(const DataTypes::ValueType& value,
		 const DataTypes::ShapeType& shape,
                 const FunctionSpace& what,
                 bool expanded)
{
  //
  // Construct a Data object of the appropriate type.
  // Construct the object first as there seems to be a bug which causes
  // undefined behaviour if an exception is thrown during construction
  // within the shared_ptr constructor.
  if (expanded) {
    DataAbstract* temp=new DataExpanded(what, shape, value);
//     boost::shared_ptr<DataAbstract> temp_data(temp);
//     m_data=temp_data;
    m_data=temp->getPtr();
  } else {
    DataAbstract* temp=new DataConstant(what, shape, value);
//     boost::shared_ptr<DataAbstract> temp_data(temp);
//     m_data=temp_data;
    m_data=temp->getPtr();
  }
}


// void
// Data::CompareDebug(const Data& rd)
// {
// 	using namespace std;
// 	bool mismatch=false;
// 	std::cout << "Comparing left and right" << endl;
// 	const DataTagged* left=dynamic_cast<DataTagged*>(m_data.get());
// 	const DataTagged* right=dynamic_cast<DataTagged*>(rd.m_data.get());
// 	
// 	if (left==0)
// 	{
// 		cout << "left arg is not a DataTagged\n";
// 		return;
// 	}
// 	
// 	if (right==0)
// 	{
// 		cout << "right arg is not a DataTagged\n";
// 		return;
// 	}
// 	cout << "Num elements=" << left->getVector().size() << ":" << right->getVector().size() << std::endl;
// 	cout << "Shapes ";
// 	if (left->getShape()==right->getShape())
// 	{
// 		cout << "ok\n";
// 	}
// 	else
// 	{
// 		cout << "Problem: shapes do not match\n";
// 		mismatch=true;
// 	}
// 	int lim=left->getVector().size();
// 	if (right->getVector().size()) lim=right->getVector().size();
// 	for (int i=0;i<lim;++i)
// 	{
// 		if (left->getVector()[i]!=right->getVector()[i])
// 		{
// 			cout << "[" << i << "] value mismatch " << left->getVector()[i] << ":" << right->getVector()[i] << endl;
// 			mismatch=true;
// 		}
// 	}
// 
// 	// still need to check the tag map
// 	// also need to watch what is happening to function spaces, are they copied or what?
// 
// 	const DataTagged::DataMapType& mapleft=left->getTagLookup();
// 	const DataTagged::DataMapType& mapright=right->getTagLookup();
// 
// 	if (mapleft.size()!=mapright.size())
// 	{
// 		cout << "Maps are different sizes " << mapleft.size() << ":" << mapright.size() << endl;
// 		mismatch=true;
// 		cout << "Left map\n";
// 		DataTagged::DataMapType::const_iterator i,j;
// 		for (i=mapleft.begin();i!=mapleft.end();++i) {
// 			cout << "(" << i->first << "=>" << i->second << ")\n";
// 		}
// 		cout << "Right map\n";
// 		for (i=mapright.begin();i!=mapright.end();++i) {
// 			cout << "(" << i->first << "=>" << i->second << ")\n";
// 		}
// 		cout << "End map\n";
// 
// 	}
// 
// 	DataTagged::DataMapType::const_iterator i,j;
// 	for (i=mapleft.begin(),j=mapright.begin();i!=mapleft.end() && j!=mapright.end();++i,++j) {
// 	   if ((i->first!=j->first) || (i->second!=j->second))
// 	   {
// 		cout << "(" << i->first << "=>" << i->second << ")";
// 		cout << ":(" << j->first << "=>" << j->second << ") ";
// 		mismatch=true;
//            }
// 	}
// 	if (mismatch)
// 	{
// 		cout << "#Mismatch\n";
// 	}
// }

escriptDataC
Data::getDataC()
{
  escriptDataC temp;
  temp.m_dataPtr=(void*)this;
  return temp;
}

escriptDataC
Data::getDataC() const
{
  escriptDataC temp;
  temp.m_dataPtr=(void*)this;
  return temp;
}

const boost::python::tuple
Data::getShapeTuple() const
{
  const DataTypes::ShapeType& shape=getDataPointShape();
  switch(getDataPointRank()) {
     case 0:
        return make_tuple();
     case 1:
        return make_tuple(long_(shape[0]));
     case 2:
        return make_tuple(long_(shape[0]),long_(shape[1]));
     case 3:
        return make_tuple(long_(shape[0]),long_(shape[1]),long_(shape[2]));
     case 4:
        return make_tuple(long_(shape[0]),long_(shape[1]),long_(shape[2]),long_(shape[3]));
     default:
        throw DataException("Error - illegal Data rank.");
  }
}


// The different name is needed because boost has trouble with overloaded functions.
// It can't work out what type the function is based soley on its name.
// There are ways to fix this involving creating function pointer variables for each form
// but there doesn't seem to be a need given that the methods have the same name from the python point of view
Data
Data::copySelf()
{
   DataAbstract* temp=m_data->deepCopy();
   return Data(temp);
}

void
Data::copy(const Data& other)
{
  DataAbstract* temp=other.m_data->deepCopy();
  DataAbstract_ptr p=temp->getPtr();
  m_data=p;
}


Data
Data::delay()
{
  DataLazy* dl=new DataLazy(m_data);
  return Data(dl);
}

void
Data::delaySelf()
{
  if (!isLazy())
  {
	m_data=(new DataLazy(m_data))->getPtr();
  }
}

void
Data::setToZero()
{
  if (isEmpty())
  {
     throw DataException("Error - Operations not permitted on instances of DataEmpty.");
  }
  m_data->setToZero();
}

void
Data::copyWithMask(const Data& other,
                   const Data& mask)
{
  // 1. Interpolate if required so all Datas use the same FS as this
  // 2. Tag or Expand so that all Data's are the same type
  // 3. Iterate over the data vectors copying values where mask is >0
  if (other.isEmpty() || mask.isEmpty())
  {
	throw DataException("Error - copyWithMask not permitted using instances of DataEmpty.");
  }
  Data other2(other);
  Data mask2(mask);
  other2.resolve();
  mask2.resolve();
  this->resolve();
  FunctionSpace myFS=getFunctionSpace();
  FunctionSpace oFS=other2.getFunctionSpace();
  FunctionSpace mFS=mask2.getFunctionSpace();
  if (oFS!=myFS)
  {
     if (other2.probeInterpolation(myFS))
     {
	other2=other2.interpolate(myFS);
     }
     else
     {
	throw DataException("Error - copyWithMask: other FunctionSpace is not compatible with this one.");
     }
  }
  if (mFS!=myFS)
  {
     if (mask2.probeInterpolation(myFS))
     {
	mask2=mask2.interpolate(myFS);
     }
     else
     {
	throw DataException("Error - copyWithMask: mask FunctionSpace is not compatible with this one.");
     }
  }
  			// Ensure that all args have the same type
  if (this->isExpanded() || mask2.isExpanded() || other2.isExpanded())
  {
	this->expand();
	other2.expand();
	mask2.expand();
  }
  else if (this->isTagged() || mask2.isTagged() || other2.isTagged())
  {
	this->tag();
	other2.tag();
	mask2.tag();
  }
  else if (this->isConstant() && mask2.isConstant() && other2.isConstant())
  {
  }
  else
  {
	throw DataException("Error - Unknown DataAbstract passed to copyWithMask.");
  }
  // Now we iterate over the elements
  DataVector& self=getReadyPtr()->getVector();
  const DataVector& ovec=other2.getReadyPtr()->getVector();
  const DataVector& mvec=mask2.getReadyPtr()->getVector();
  if ((self.size()!=ovec.size()) || (self.size()!=mvec.size()))
  {
	throw DataException("Error - size mismatch in arguments to copyWithMask.");
  }
  size_t num_points=self.size();

  // OPENMP 3.0 allows unsigned loop vars.
#if defined(_OPENMP) && (_OPENMP < 200805)
  long i;
#else
  size_t i;
#endif
  #pragma omp parallel for private(i) schedule(static)
  for (i=0;i<num_points;++i)
  {
	if (mvec[i]>0)
	{
	   self[i]=ovec[i];
	}
  }
}



bool
Data::isExpanded() const
{
  DataExpanded* temp=dynamic_cast<DataExpanded*>(m_data.get());
  return (temp!=0);
}

bool
Data::isTagged() const
{
  DataTagged* temp=dynamic_cast<DataTagged*>(m_data.get());
  return (temp!=0);
}

bool
Data::isEmpty() const
{
  DataEmpty* temp=dynamic_cast<DataEmpty*>(m_data.get());
  return (temp!=0);
}

bool
Data::isConstant() const
{
  DataConstant* temp=dynamic_cast<DataConstant*>(m_data.get());
  return (temp!=0);
}

bool
Data::isLazy() const
{
  return m_data->isLazy();
}

// at the moment this is synonymous with !isLazy() but that could change
bool
Data::isReady() const
{
  return (dynamic_cast<DataReady*>(m_data.get())!=0);
}


void
Data::setProtection()
{
   m_protected=true;
}

bool
Data::isProtected() const
{
   return m_protected;
}



void
Data::expand()
{
  if (isConstant()) {
    DataConstant* tempDataConst=dynamic_cast<DataConstant*>(m_data.get());
    DataAbstract* temp=new DataExpanded(*tempDataConst);
//     shared_ptr<DataAbstract> temp_data(temp);
//     m_data=temp_data;
    m_data=temp->getPtr();
  } else if (isTagged()) {
    DataTagged* tempDataTag=dynamic_cast<DataTagged*>(m_data.get());
    DataAbstract* temp=new DataExpanded(*tempDataTag);
//     shared_ptr<DataAbstract> temp_data(temp);
//     m_data=temp_data;
    m_data=temp->getPtr();
  } else if (isExpanded()) {
    //
    // do nothing
  } else if (isEmpty()) {
    throw DataException("Error - Expansion of DataEmpty not possible.");
  } else if (isLazy()) {
    resolve();
    expand();		// resolve might not give us expanded data
  } else {
    throw DataException("Error - Expansion not implemented for this Data type.");
  }
}

void
Data::tag()
{
  if (isConstant()) {
    DataConstant* tempDataConst=dynamic_cast<DataConstant*>(m_data.get());
    DataAbstract* temp=new DataTagged(*tempDataConst);
//     shared_ptr<DataAbstract> temp_data(temp);
//     m_data=temp_data;
    m_data=temp->getPtr();
  } else if (isTagged()) {
    // do nothing
  } else if (isExpanded()) {
    throw DataException("Error - Creating tag data from DataExpanded not possible.");
  } else if (isEmpty()) {
    throw DataException("Error - Creating tag data from DataEmpty not possible.");
  } else if (isLazy()) {
     DataAbstract_ptr res=m_data->resolve();
     if (m_data->isExpanded())
     {
	throw DataException("Error - data would resolve to DataExpanded, tagging is not possible.");
     }
     m_data=res;	
     tag();
  } else {
    throw DataException("Error - Tagging not implemented for this Data type.");
  }
}

void
Data::resolve()
{
  if (isLazy())
  {
     m_data=m_data->resolve();
  }
}


Data
Data::oneOver() const
{
  MAKELAZYOP(RECIP)
  return C_TensorUnaryOperation(*this, bind1st(divides<double>(),1.));
}

Data
Data::wherePositive() const
{
  MAKELAZYOP(GZ)
  return C_TensorUnaryOperation(*this, bind2nd(greater<double>(),0.0));
}

Data
Data::whereNegative() const
{
  MAKELAZYOP(LZ)
  return C_TensorUnaryOperation(*this, bind2nd(less<double>(),0.0));
}

Data
Data::whereNonNegative() const
{
  MAKELAZYOP(GEZ)
  return C_TensorUnaryOperation(*this, bind2nd(greater_equal<double>(),0.0));
}

Data
Data::whereNonPositive() const
{
  MAKELAZYOP(LEZ)
  return C_TensorUnaryOperation(*this, bind2nd(less_equal<double>(),0.0));
}

Data
Data::whereZero(double tol) const
{
//   Data dataAbs=abs();
//   return C_TensorUnaryOperation(dataAbs, bind2nd(less_equal<double>(),tol));
   MAKELAZYOPOFF(EZ,tol)
   return C_TensorUnaryOperation(*this, bind2nd(AbsLTE(),tol));

}

Data
Data::whereNonZero(double tol) const
{
//   Data dataAbs=abs();
//   return C_TensorUnaryOperation(dataAbs, bind2nd(greater<double>(),tol));
  MAKELAZYOPOFF(NEZ,tol)
  return C_TensorUnaryOperation(*this, bind2nd(AbsGT(),tol));

}

Data
Data::interpolate(const FunctionSpace& functionspace) const
{
  return Data(*this,functionspace);
}

bool
Data::probeInterpolation(const FunctionSpace& functionspace) const
{
  return getFunctionSpace().probeInterpolation(functionspace);
}

Data
Data::gradOn(const FunctionSpace& functionspace) const
{
  if (isEmpty())
  {
	throw DataException("Error - operation not permitted on instances of DataEmpty.");
  }
  double blocktimer_start = blocktimer_time();
  if (functionspace.getDomain()!=getDomain())
    throw DataException("Error - gradient cannot be calculated on different domains.");
  DataTypes::ShapeType grad_shape=getDataPointShape();
  grad_shape.push_back(functionspace.getDim());
  Data out(0.0,grad_shape,functionspace,true);
  getDomain()->setToGradient(out,*this);
  blocktimer_increment("grad()", blocktimer_start);
  return out;
}

Data
Data::grad() const
{
  if (isEmpty())
  {
	throw DataException("Error - operation not permitted on instances of DataEmpty.");
  }
  return gradOn(escript::function(*getDomain()));
}

int
Data::getDataPointSize() const
{
  return m_data->getNoValues();
}

DataTypes::ValueType::size_type
Data::getLength() const
{
  return m_data->getLength();
}

const
boost::python::numeric::array
Data:: getValueOfDataPoint(int dataPointNo)
{
  int i, j, k, l;

  FORCERESOLVE;

  //
  // determine the rank and shape of each data point
  int dataPointRank = getDataPointRank();
  const DataTypes::ShapeType& dataPointShape = getDataPointShape();

  //
  // create the numeric array to be returned
  boost::python::numeric::array numArray(0.0);

  //
  // the shape of the returned numeric array will be the same
  // as that of the data point
  int arrayRank = dataPointRank;
  const DataTypes::ShapeType& arrayShape = dataPointShape;

  //
  // resize the numeric array to the shape just calculated
  if (arrayRank==0) {
    numArray.resize(1);
  }
  if (arrayRank==1) {
    numArray.resize(arrayShape[0]);
  }
  if (arrayRank==2) {
    numArray.resize(arrayShape[0],arrayShape[1]);
  }
  if (arrayRank==3) {
    numArray.resize(arrayShape[0],arrayShape[1],arrayShape[2]);
  }
  if (arrayRank==4) {
    numArray.resize(arrayShape[0],arrayShape[1],arrayShape[2],arrayShape[3]);
  }

  if (getNumDataPointsPerSample()>0) {
       int sampleNo = dataPointNo/getNumDataPointsPerSample();
       int dataPointNoInSample = dataPointNo - sampleNo * getNumDataPointsPerSample();
       //
       // Check a valid sample number has been supplied
       if ((sampleNo >= getNumSamples()) || (sampleNo < 0 )) {
           throw DataException("Error - Data::convertToNumArray: invalid sampleNo.");
       }

       //
       // Check a valid data point number has been supplied
       if ((dataPointNoInSample >= getNumDataPointsPerSample()) || (dataPointNoInSample < 0)) {
           throw DataException("Error - Data::convertToNumArray: invalid dataPointNoInSample.");
       }
       // TODO: global error handling
       // create a view of the data if it is stored locally
//       DataArrayView dataPointView = getDataPoint(sampleNo, dataPointNoInSample);
       DataTypes::ValueType::size_type offset=getDataOffset(sampleNo, dataPointNoInSample);


       switch( dataPointRank ){
			case 0 :
				numArray[0] = getDataAtOffset(offset);
				break;
			case 1 :
				for( i=0; i<dataPointShape[0]; i++ )
					numArray[i]=getDataAtOffset(offset+DataTypes::getRelIndex(dataPointShape, i));
				break;
			case 2 :
				for( i=0; i<dataPointShape[0]; i++ )
					for( j=0; j<dataPointShape[1]; j++)
						numArray[make_tuple(i,j)]=getDataAtOffset(offset+DataTypes::getRelIndex(dataPointShape, i,j));
				break;
			case 3 :
				for( i=0; i<dataPointShape[0]; i++ )
					for( j=0; j<dataPointShape[1]; j++ )
						for( k=0; k<dataPointShape[2]; k++)
							numArray[make_tuple(i,j,k)]=getDataAtOffset(offset+DataTypes::getRelIndex(dataPointShape, i,j,k));
				break;
			case 4 :
				for( i=0; i<dataPointShape[0]; i++ )
					for( j=0; j<dataPointShape[1]; j++ )
						for( k=0; k<dataPointShape[2]; k++ )
							for( l=0; l<dataPointShape[3]; l++)
								numArray[make_tuple(i,j,k,l)]=getDataAtOffset(offset+DataTypes::getRelIndex(dataPointShape, i,j,k,l));
				break;
	}
  }
  //
  // return the array
  return numArray;
}

void
Data::setValueOfDataPointToPyObject(int dataPointNo, const boost::python::object& py_object)
{
    // this will throw if the value cannot be represented
    boost::python::numeric::array num_array(py_object);
    setValueOfDataPointToArray(dataPointNo,num_array);
}

void
Data::setValueOfDataPointToArray(int dataPointNo, const boost::python::numeric::array& num_array)
{
  if (isProtected()) {
        throw DataException("Error - attempt to update protected Data object.");
  }
  FORCERESOLVE;
  //
  // check rank
  if (static_cast<unsigned int>(num_array.getrank())<getDataPointRank())
      throw DataException("Rank of numarray does not match Data object rank");

  //
  // check shape of num_array
  for (unsigned int i=0; i<getDataPointRank(); i++) {
    if (extract<int>(num_array.getshape()[i])!=getDataPointShape()[i])
       throw DataException("Shape of numarray does not match Data object rank");
  }
  //
  // make sure data is expanded:
  //
  if (!isExpanded()) {
    expand();
  }
  if (getNumDataPointsPerSample()>0) {
       int sampleNo = dataPointNo/getNumDataPointsPerSample();
       int dataPointNoInSample = dataPointNo - sampleNo * getNumDataPointsPerSample();
       m_data->copyToDataPoint(sampleNo, dataPointNoInSample,num_array);
  } else {
       m_data->copyToDataPoint(-1, 0,num_array);
  }
}

void
Data::setValueOfDataPoint(int dataPointNo, const double value)
{
  if (isProtected()) {
        throw DataException("Error - attempt to update protected Data object.");
  }
  //
  // make sure data is expanded:
  FORCERESOLVE;
  if (!isExpanded()) {
    expand();
  }
  if (getNumDataPointsPerSample()>0) {
       int sampleNo = dataPointNo/getNumDataPointsPerSample();
       int dataPointNoInSample = dataPointNo - sampleNo * getNumDataPointsPerSample();
       m_data->copyToDataPoint(sampleNo, dataPointNoInSample,value);
  } else {
       m_data->copyToDataPoint(-1, 0,value);
  }
}

const
boost::python::numeric::array
Data::getValueOfGlobalDataPoint(int procNo, int dataPointNo)
{
  size_t length=0;
  int i, j, k, l, pos;
  FORCERESOLVE;
  //
  // determine the rank and shape of each data point
  int dataPointRank = getDataPointRank();
  const DataTypes::ShapeType& dataPointShape = getDataPointShape();

  //
  // create the numeric array to be returned
  boost::python::numeric::array numArray(0.0);

  //
  // the shape of the returned numeric array will be the same
  // as that of the data point
  int arrayRank = dataPointRank;
  const DataTypes::ShapeType& arrayShape = dataPointShape;

  //
  // resize the numeric array to the shape just calculated
  if (arrayRank==0) {
    numArray.resize(1);
  }
  if (arrayRank==1) {
    numArray.resize(arrayShape[0]);
  }
  if (arrayRank==2) {
    numArray.resize(arrayShape[0],arrayShape[1]);
  }
  if (arrayRank==3) {
    numArray.resize(arrayShape[0],arrayShape[1],arrayShape[2]);
  }
  if (arrayRank==4) {
    numArray.resize(arrayShape[0],arrayShape[1],arrayShape[2],arrayShape[3]);
  }

  // added for the MPI communication
  length=1;
  for( i=0; i<arrayRank; i++ ) length *= arrayShape[i];
  double *tmpData = new double[length];

  //
  // load the values for the data point into the numeric array.

	// updated for the MPI case
	if( get_MPIRank()==procNo ){
             if (getNumDataPointsPerSample()>0) {
                int sampleNo = dataPointNo/getNumDataPointsPerSample();
                int dataPointNoInSample = dataPointNo - sampleNo * getNumDataPointsPerSample();
                //
                // Check a valid sample number has been supplied
                if ((sampleNo >= getNumSamples()) || (sampleNo < 0 )) {
                  throw DataException("Error - Data::convertToNumArray: invalid sampleNo.");
                }

                //
                // Check a valid data point number has been supplied
                if ((dataPointNoInSample >= getNumDataPointsPerSample()) || (dataPointNoInSample < 0)) {
                  throw DataException("Error - Data::convertToNumArray: invalid dataPointNoInSample.");
                }
                // TODO: global error handling
		// create a view of the data if it is stored locally
		//DataArrayView dataPointView = getDataPoint(sampleNo, dataPointNoInSample);
		DataTypes::ValueType::size_type offset=getDataOffset(sampleNo, dataPointNoInSample);

		// pack the data from the view into tmpData for MPI communication
		pos=0;
		switch( dataPointRank ){
			case 0 :
				tmpData[0] = getDataAtOffset(offset);
				break;
			case 1 :
				for( i=0; i<dataPointShape[0]; i++ )
					tmpData[i]=getDataAtOffset(offset+DataTypes::getRelIndex(dataPointShape, i));
				break;
			case 2 :
				for( i=0; i<dataPointShape[0]; i++ )
					for( j=0; j<dataPointShape[1]; j++, pos++ )
						tmpData[pos]=getDataAtOffset(offset+DataTypes::getRelIndex(dataPointShape, i,j));
				break;
			case 3 :
				for( i=0; i<dataPointShape[0]; i++ )
					for( j=0; j<dataPointShape[1]; j++ )
						for( k=0; k<dataPointShape[2]; k++, pos++ )
							tmpData[pos]=getDataAtOffset(offset+DataTypes::getRelIndex(dataPointShape, i,j,k));
				break;
			case 4 :
				for( i=0; i<dataPointShape[0]; i++ )
					for( j=0; j<dataPointShape[1]; j++ )
						for( k=0; k<dataPointShape[2]; k++ )
							for( l=0; l<dataPointShape[3]; l++, pos++ )
								tmpData[pos]=getDataAtOffset(offset+DataTypes::getRelIndex(dataPointShape, i,j,k,l));
				break;
		}
            }
	}
        #ifdef PASO_MPI
        // broadcast the data to all other processes
	MPI_Bcast( tmpData, length, MPI_DOUBLE, procNo, get_MPIComm() );
        #endif

	// unpack the data
	switch( dataPointRank ){
		case 0 :
			numArray[0]=tmpData[0];
			break;
		case 1 :
			for( i=0; i<dataPointShape[0]; i++ )
				numArray[i]=tmpData[i];
			break;
		case 2 :
			for( i=0; i<dataPointShape[0]; i++ )
				for( j=0; j<dataPointShape[1]; j++ )
				   numArray[make_tuple(i,j)]=tmpData[i+j*dataPointShape[0]];
			break;
		case 3 :
			for( i=0; i<dataPointShape[0]; i++ )
				for( j=0; j<dataPointShape[1]; j++ )
					for( k=0; k<dataPointShape[2]; k++ )
						numArray[make_tuple(i,j,k)]=tmpData[i+dataPointShape[0]*(j*+k*dataPointShape[1])];
			break;
		case 4 :
			for( i=0; i<dataPointShape[0]; i++ )
				for( j=0; j<dataPointShape[1]; j++ )
					for( k=0; k<dataPointShape[2]; k++ )
						for( l=0; l<dataPointShape[3]; l++ )
					        	numArray[make_tuple(i,j,k,l)]=tmpData[i+dataPointShape[0]*(j*+dataPointShape[1]*(k+l*dataPointShape[2]))];
			break;
	}

	delete [] tmpData;
  //
  // return the loaded array
  return numArray;
}


boost::python::numeric::array
Data::integrate_const() const
{
  if (isLazy())
  {
	throw DataException("Error - cannot integrate for constant lazy data.");
  }
  return integrateWorker();
}

boost::python::numeric::array
Data::integrate()
{
  if (isLazy())
  {
	expand(); 
  }
  return integrateWorker();
}



boost::python::numeric::array
Data::integrateWorker() const
{
  int index;
  int rank = getDataPointRank();
  DataTypes::ShapeType shape = getDataPointShape();
  int dataPointSize = getDataPointSize();

  //
  // calculate the integral values
  vector<double> integrals(dataPointSize);
  vector<double> integrals_local(dataPointSize);
  const AbstractContinuousDomain* dom=dynamic_cast<const AbstractContinuousDomain*>(getDomain().get());
  if (dom==0)
  {				
    throw DataException("Can not integrate over non-continuous domains.");
  }
#ifdef PASO_MPI
  dom->setToIntegrals(integrals_local,*this);
  // Global sum: use an array instead of a vector because elements of array are guaranteed to be contiguous in memory
  double *tmp = new double[dataPointSize];
  double *tmp_local = new double[dataPointSize];
  for (int i=0; i<dataPointSize; i++) { tmp_local[i] = integrals_local[i]; }
  MPI_Allreduce( &tmp_local[0], &tmp[0], dataPointSize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
  for (int i=0; i<dataPointSize; i++) { integrals[i] = tmp[i]; }
  delete[] tmp;
  delete[] tmp_local;
#else
  dom->setToIntegrals(integrals,*this);
#endif

  //
  // create the numeric array to be returned
  // and load the array with the integral values
  boost::python::numeric::array bp_array(1.0);
  if (rank==0) {
    bp_array.resize(1);
    index = 0;
    bp_array[0] = integrals[index];
  }
  if (rank==1) {
    bp_array.resize(shape[0]);
    for (int i=0; i<shape[0]; i++) {
      index = i;
      bp_array[i] = integrals[index];
    }
  }
  if (rank==2) {
       bp_array.resize(shape[0],shape[1]);
       for (int i=0; i<shape[0]; i++) {
         for (int j=0; j<shape[1]; j++) {
           index = i + shape[0] * j;
           bp_array[make_tuple(i,j)] = integrals[index];
         }
       }
  }
  if (rank==3) {
    bp_array.resize(shape[0],shape[1],shape[2]);
    for (int i=0; i<shape[0]; i++) {
      for (int j=0; j<shape[1]; j++) {
        for (int k=0; k<shape[2]; k++) {
          index = i + shape[0] * ( j + shape[1] * k );
          bp_array[make_tuple(i,j,k)] = integrals[index];
        }
      }
    }
  }
  if (rank==4) {
    bp_array.resize(shape[0],shape[1],shape[2],shape[3]);
    for (int i=0; i<shape[0]; i++) {
      for (int j=0; j<shape[1]; j++) {
        for (int k=0; k<shape[2]; k++) {
          for (int l=0; l<shape[3]; l++) {
            index = i + shape[0] * ( j + shape[1] * ( k + shape[2] * l ) );
            bp_array[make_tuple(i,j,k,l)] = integrals[index];
          }
        }
      }
    }
  }

  //
  // return the loaded array
  return bp_array;
}

Data
Data::sin() const
{
  MAKELAZYOP(SIN)
  return C_TensorUnaryOperation<double (*)(double)>(*this, ::sin);
}

Data
Data::cos() const
{
  MAKELAZYOP(COS)
  return C_TensorUnaryOperation<double (*)(double)>(*this, ::cos);
}

Data
Data::tan() const
{
  MAKELAZYOP(TAN)
  return C_TensorUnaryOperation<double (*)(double)>(*this, ::tan);
}

Data
Data::asin() const
{
  MAKELAZYOP(ASIN)
  return C_TensorUnaryOperation<double (*)(double)>(*this, ::asin);
}

Data
Data::acos() const
{
  MAKELAZYOP(ACOS)
  return C_TensorUnaryOperation<double (*)(double)>(*this, ::acos);
}


Data
Data::atan() const
{
  MAKELAZYOP(ATAN)
  return C_TensorUnaryOperation<double (*)(double)>(*this, ::atan);
}

Data
Data::sinh() const
{
  MAKELAZYOP(SINH)
  return C_TensorUnaryOperation<double (*)(double)>(*this, ::sinh);
}

Data
Data::cosh() const
{
  MAKELAZYOP(COSH)
  return C_TensorUnaryOperation<double (*)(double)>(*this, ::cosh);
}

Data
Data::tanh() const
{
  MAKELAZYOP(TANH)
  return C_TensorUnaryOperation<double (*)(double)>(*this, ::tanh);
}


Data
Data::erf() const
{
#if defined (_WIN32) && !defined(__INTEL_COMPILER)
  throw DataException("Error - Data:: erf function is not supported on _WIN32 platforms.");
#else
  MAKELAZYOP(ERF)
  return C_TensorUnaryOperation(*this, ::erf);
#endif
}

Data
Data::asinh() const
{
  MAKELAZYOP(ASINH)
#if defined (_WIN32) && !defined(__INTEL_COMPILER)
  return C_TensorUnaryOperation(*this, escript::asinh_substitute);
#else
  return C_TensorUnaryOperation(*this, ::asinh);
#endif
}

Data
Data::acosh() const
{
  MAKELAZYOP(ACOSH)
#if defined (_WIN32) && !defined(__INTEL_COMPILER)
  return C_TensorUnaryOperation(*this, escript::acosh_substitute);
#else
  return C_TensorUnaryOperation(*this, ::acosh);
#endif
}

Data
Data::atanh() const
{
  MAKELAZYOP(ATANH)
#if defined (_WIN32) && !defined(__INTEL_COMPILER)
  return C_TensorUnaryOperation(*this, escript::atanh_substitute);
#else
  return C_TensorUnaryOperation(*this, ::atanh);
#endif
}

Data
Data::log10() const
{
  MAKELAZYOP(LOG10)
  return C_TensorUnaryOperation<double (*)(double)>(*this, ::log10);
}

Data
Data::log() const
{
  MAKELAZYOP(LOG)
  return C_TensorUnaryOperation<double (*)(double)>(*this, ::log);
}

Data
Data::sign() const
{
  MAKELAZYOP(SIGN)
  return C_TensorUnaryOperation(*this, escript::fsign);
}

Data
Data::abs() const
{
  MAKELAZYOP(ABS)
  return C_TensorUnaryOperation<double (*)(double)>(*this, ::fabs);
}

Data
Data::neg() const
{
  MAKELAZYOP(NEG)
  return C_TensorUnaryOperation(*this, negate<double>());
}

Data
Data::pos() const
{
	// not doing lazy check here is deliberate.
	// since a deep copy of lazy data should be cheap, I'll just let it happen now
  Data result;
  // perform a deep copy
  result.copy(*this);
  return result;
}

Data
Data::exp() const
{
  MAKELAZYOP(EXP)
  return C_TensorUnaryOperation<double (*)(double)>(*this, ::exp);
}

Data
Data::sqrt() const
{
  MAKELAZYOP(SQRT)
  return C_TensorUnaryOperation<double (*)(double)>(*this, ::sqrt);
}

double
Data::Lsup_const() const
{
   if (isLazy())
   {
	throw DataException("Error - cannot compute Lsup for constant lazy data.");
   }
   return LsupWorker();
}

double
Data::Lsup() 
{
   if (isLazy())
   {
	resolve();
   }
   return LsupWorker();
}

double
Data::sup_const() const
{
   if (isLazy())
   {
	throw DataException("Error - cannot compute sup for constant lazy data.");
   }
   return supWorker();
}

double
Data::sup() 
{
   if (isLazy())
   {
	resolve();
   }
   return supWorker();
}

double
Data::inf_const() const
{
   if (isLazy())
   {
	throw DataException("Error - cannot compute inf for constant lazy data.");
   }
   return infWorker();
}

double
Data::inf() 
{
   if (isLazy())
   {
	resolve();
   }
   return infWorker();
}

double
Data::LsupWorker() const
{
  double localValue;
  //
  // set the initial absolute maximum value to zero

  AbsMax abs_max_func;
  localValue = algorithm(abs_max_func,0);
#ifdef PASO_MPI
  double globalValue;
  MPI_Allreduce( &localValue, &globalValue, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
  return globalValue;
#else
  return localValue;
#endif
}

double
Data::supWorker() const
{
  double localValue;
  //
  // set the initial maximum value to min possible double
  FMax fmax_func;
  localValue = algorithm(fmax_func,numeric_limits<double>::max()*-1);
#ifdef PASO_MPI
  double globalValue;
  MPI_Allreduce( &localValue, &globalValue, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
  return globalValue;
#else
  return localValue;
#endif
}

double
Data::infWorker() const
{
  double localValue;
  //
  // set the initial minimum value to max possible double
  FMin fmin_func;
  localValue = algorithm(fmin_func,numeric_limits<double>::max());
#ifdef PASO_MPI
  double globalValue;
  MPI_Allreduce( &localValue, &globalValue, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
  return globalValue;
#else
  return localValue;
#endif
}

/* TODO */
/* global reduction */
Data
Data::maxval() const
{
  if (isLazy())
  {
	Data temp(*this);	// to get around the fact that you can't resolve a const Data
	temp.resolve();
	return temp.maxval();
  }
  //
  // set the initial maximum value to min possible double
  FMax fmax_func;
  return dp_algorithm(fmax_func,numeric_limits<double>::max()*-1);
}

Data
Data::minval() const
{
  if (isLazy())
  {
	Data temp(*this);	// to get around the fact that you can't resolve a const Data
	temp.resolve();
	return temp.minval();
  }
  //
  // set the initial minimum value to max possible double
  FMin fmin_func;
  return dp_algorithm(fmin_func,numeric_limits<double>::max());
}

Data
Data::swapaxes(const int axis0, const int axis1) const
{
     if (isLazy())
     {
	Data temp(*this);
	temp.resolve();
	return temp.swapaxes(axis0,axis1);
     }
     int axis0_tmp,axis1_tmp;
     DataTypes::ShapeType s=getDataPointShape();
     DataTypes::ShapeType ev_shape;
     // Here's the equivalent of python s_out=s[axis_offset:]+s[:axis_offset]
     // which goes thru all shape vector elements starting with axis_offset (at index=rank wrap around to 0)
     int rank=getDataPointRank();
     if (rank<2) {
        throw DataException("Error - Data::swapaxes argument must have at least rank 2.");
     }
     if (axis0<0 || axis0>rank-1) {
        throw DataException("Error - Data::swapaxes: axis0 must be between 0 and rank-1=" + rank-1);
     }
     if (axis1<0 || axis1>rank-1) {
         throw DataException("Error - Data::swapaxes: axis1 must be between 0 and rank-1=" + rank-1);
     }
     if (axis0 == axis1) {
         throw DataException("Error - Data::swapaxes: axis indices must be different.");
     }
     if (axis0 > axis1) {
         axis0_tmp=axis1;
         axis1_tmp=axis0;
     } else {
         axis0_tmp=axis0;
         axis1_tmp=axis1;
     }
     for (int i=0; i<rank; i++) {
       if (i == axis0_tmp) {
          ev_shape.push_back(s[axis1_tmp]);
       } else if (i == axis1_tmp) {
          ev_shape.push_back(s[axis0_tmp]);
       } else {
          ev_shape.push_back(s[i]);
       }
     }
     Data ev(0.,ev_shape,getFunctionSpace());
     ev.typeMatchRight(*this);
     m_data->swapaxes(ev.m_data.get(), axis0_tmp, axis1_tmp);
     return ev;

}

Data
Data::symmetric() const
{
     // check input
     DataTypes::ShapeType s=getDataPointShape();
     if (getDataPointRank()==2) {
        if(s[0] != s[1])
           throw DataException("Error - Data::symmetric can only be calculated for rank 2 object with equal first and second dimension.");
     }
     else if (getDataPointRank()==4) {
        if(!(s[0] == s[2] && s[1] == s[3]))
           throw DataException("Error - Data::symmetric can only be calculated for rank 4 object with dim0==dim2 and dim1==dim3.");
     }
     else {
        throw DataException("Error - Data::symmetric can only be calculated for rank 2 or 4 object.");
     }
     MAKELAZYOP(SYM)
     Data ev(0.,getDataPointShape(),getFunctionSpace());
     ev.typeMatchRight(*this);
     m_data->symmetric(ev.m_data.get());
     return ev;
}

Data
Data::nonsymmetric() const
{
     MAKELAZYOP(NSYM)
     // check input
     DataTypes::ShapeType s=getDataPointShape();
     if (getDataPointRank()==2) {
        if(s[0] != s[1])
           throw DataException("Error - Data::nonsymmetric can only be calculated for rank 2 object with equal first and second dimension.");
        DataTypes::ShapeType ev_shape;
        ev_shape.push_back(s[0]);
        ev_shape.push_back(s[1]);
        Data ev(0.,ev_shape,getFunctionSpace());
        ev.typeMatchRight(*this);
        m_data->nonsymmetric(ev.m_data.get());
        return ev;
     }
     else if (getDataPointRank()==4) {
        if(!(s[0] == s[2] && s[1] == s[3]))
           throw DataException("Error - Data::nonsymmetric can only be calculated for rank 4 object with dim0==dim2 and dim1==dim3.");
        DataTypes::ShapeType ev_shape;
        ev_shape.push_back(s[0]);
        ev_shape.push_back(s[1]);
        ev_shape.push_back(s[2]);
        ev_shape.push_back(s[3]);
        Data ev(0.,ev_shape,getFunctionSpace());
        ev.typeMatchRight(*this);
        m_data->nonsymmetric(ev.m_data.get());
        return ev;
     }
     else {
        throw DataException("Error - Data::nonsymmetric can only be calculated for rank 2 or 4 object.");
     }
}

Data
Data::trace(int axis_offset) const
{     
     MAKELAZYOPOFF(TRACE,axis_offset)
     DataTypes::ShapeType s=getDataPointShape();
     if (getDataPointRank()==2) {
        DataTypes::ShapeType ev_shape;
        Data ev(0.,ev_shape,getFunctionSpace());
        ev.typeMatchRight(*this);
        m_data->trace(ev.m_data.get(), axis_offset);
        return ev;
     }
     if (getDataPointRank()==3) {
        DataTypes::ShapeType ev_shape;
        if (axis_offset==0) {
          int s2=s[2];
          ev_shape.push_back(s2);
        }
        else if (axis_offset==1) {
          int s0=s[0];
          ev_shape.push_back(s0);
        }
        Data ev(0.,ev_shape,getFunctionSpace());
        ev.typeMatchRight(*this);
        m_data->trace(ev.m_data.get(), axis_offset);
        return ev;
     }
     if (getDataPointRank()==4) {
        DataTypes::ShapeType ev_shape;
        if (axis_offset==0) {
          ev_shape.push_back(s[2]);
          ev_shape.push_back(s[3]);
        }
        else if (axis_offset==1) {
          ev_shape.push_back(s[0]);
          ev_shape.push_back(s[3]);
        }
	else if (axis_offset==2) {
	  ev_shape.push_back(s[0]);
	  ev_shape.push_back(s[1]);
	}
        Data ev(0.,ev_shape,getFunctionSpace());
        ev.typeMatchRight(*this);
	m_data->trace(ev.m_data.get(), axis_offset);
        return ev;
     }
     else {
        throw DataException("Error - Data::trace can only be calculated for rank 2, 3 or 4 object.");
     }
}

Data
Data::transpose(int axis_offset) const
{     
     MAKELAZYOPOFF(TRANS,axis_offset)
     DataTypes::ShapeType s=getDataPointShape();
     DataTypes::ShapeType ev_shape;
     // Here's the equivalent of python s_out=s[axis_offset:]+s[:axis_offset]
     // which goes thru all shape vector elements starting with axis_offset (at index=rank wrap around to 0)
     int rank=getDataPointRank();
     if (axis_offset<0 || axis_offset>rank) {
        throw DataException("Error - Data::transpose must have 0 <= axis_offset <= rank=" + rank);
     }
     for (int i=0; i<rank; i++) {
       int index = (axis_offset+i)%rank;
       ev_shape.push_back(s[index]); // Append to new shape
     }
     Data ev(0.,ev_shape,getFunctionSpace());
     ev.typeMatchRight(*this);
     m_data->transpose(ev.m_data.get(), axis_offset);
     return ev;
}

Data
Data::eigenvalues() const
{
     if (isLazy())
     {
	Data temp(*this);	// to get around the fact that you can't resolve a const Data
	temp.resolve();
	return temp.eigenvalues();
     }
     // check input
     DataTypes::ShapeType s=getDataPointShape();
     if (getDataPointRank()!=2)
        throw DataException("Error - Data::eigenvalues can only be calculated for rank 2 object.");
     if(s[0] != s[1])
        throw DataException("Error - Data::eigenvalues can only be calculated for object with equal first and second dimension.");
     // create return
     DataTypes::ShapeType ev_shape(1,s[0]);
     Data ev(0.,ev_shape,getFunctionSpace());
     ev.typeMatchRight(*this);
     m_data->eigenvalues(ev.m_data.get());
     return ev;
}

const boost::python::tuple
Data::eigenvalues_and_eigenvectors(const double tol) const
{
     if (isLazy())
     {
	Data temp(*this);	// to get around the fact that you can't resolve a const Data
	temp.resolve();
	return temp.eigenvalues_and_eigenvectors(tol);
     }
     DataTypes::ShapeType s=getDataPointShape();
     if (getDataPointRank()!=2)
        throw DataException("Error - Data::eigenvalues and eigenvectors can only be calculated for rank 2 object.");
     if(s[0] != s[1])
        throw DataException("Error - Data::eigenvalues and eigenvectors can only be calculated for object with equal first and second dimension.");
     // create return
     DataTypes::ShapeType ev_shape(1,s[0]);
     Data ev(0.,ev_shape,getFunctionSpace());
     ev.typeMatchRight(*this);
     DataTypes::ShapeType V_shape(2,s[0]);
     Data V(0.,V_shape,getFunctionSpace());
     V.typeMatchRight(*this);
     m_data->eigenvalues_and_eigenvectors(ev.m_data.get(),V.m_data.get(),tol);
     return make_tuple(boost::python::object(ev),boost::python::object(V));
}

const boost::python::tuple
Data::minGlobalDataPoint() const
{
  // NB: calc_minGlobalDataPoint( had to be split off from minGlobalDataPoint( as boost::make_tuple causes an
  // abort (for unknown reasons) if there are openmp directives with it in the
  // surrounding function

  int DataPointNo;
  int ProcNo;
  calc_minGlobalDataPoint(ProcNo,DataPointNo);
  return make_tuple(ProcNo,DataPointNo);
}

void
Data::calc_minGlobalDataPoint(int& ProcNo,
       	                int& DataPointNo) const
{
  if (isLazy())
  {
    Data temp(*this);	// to get around the fact that you can't resolve a const Data
    temp.resolve();
    return temp.calc_minGlobalDataPoint(ProcNo,DataPointNo);
  }
  int i,j;
  int lowi=0,lowj=0;
  double min=numeric_limits<double>::max();

  Data temp=minval();

  int numSamples=temp.getNumSamples();
  int numDPPSample=temp.getNumDataPointsPerSample();

  double next,local_min;
  int local_lowi=0,local_lowj=0;	

  #pragma omp parallel firstprivate(local_lowi,local_lowj) private(next,local_min)
  {
    local_min=min;
    #pragma omp for private(i,j) schedule(static)
    for (i=0; i<numSamples; i++) {
      for (j=0; j<numDPPSample; j++) {
        next=temp.getDataAtOffset(temp.getDataOffset(i,j));
        if (next<local_min) {
          local_min=next;
          local_lowi=i;
          local_lowj=j;
        }
      }
    }
    #pragma omp critical
    if (local_min<min) {
      min=local_min;
      lowi=local_lowi;
      lowj=local_lowj;
    }
  }

#ifdef PASO_MPI
	// determine the processor on which the minimum occurs
	next = temp.getDataPoint(lowi,lowj);
	int lowProc = 0;
	double *globalMins = new double[get_MPISize()+1];
	int error;
    error = MPI_Gather ( &next, 1, MPI_DOUBLE, globalMins, 1, MPI_DOUBLE, 0, get_MPIComm() );

	if( get_MPIRank()==0 ){
		next = globalMins[lowProc];
		for( i=1; i<get_MPISize(); i++ )
			if( next>globalMins[i] ){
				lowProc = i;
				next = globalMins[i];
			}
	}
	MPI_Bcast( &lowProc, 1, MPI_DOUBLE, 0, get_MPIComm() );

	delete [] globalMins;
	ProcNo = lowProc;
#else
	ProcNo = 0;
#endif
  DataPointNo = lowj + lowi * numDPPSample;
}

void
Data::saveDX(std::string fileName) const
{
  if (isEmpty())
  {
    throw DataException("Error - Operations not permitted on instances of DataEmpty.");
  }
  if (isLazy())
  {
     Data temp(*this);	// to get around the fact that you can't resolve a const Data
     temp.resolve();
     temp.saveDX(fileName);
     return;
  }
  boost::python::dict args;
  args["data"]=boost::python::object(this);
  getDomain()->saveDX(fileName,args);
  return;
}

void
Data::saveVTK(std::string fileName) const
{
  if (isEmpty())
  {
    throw DataException("Error - Operations not permitted on instances of DataEmpty.");
  }
  if (isLazy())
  {
     Data temp(*this);	// to get around the fact that you can't resolve a const Data
     temp.resolve();
     temp.saveVTK(fileName);
     return;
  }
  boost::python::dict args;
  args["data"]=boost::python::object(this);
  getDomain()->saveVTK(fileName,args);
  return;
}

Data&
Data::operator+=(const Data& right)
{
  if (isProtected()) {
        throw DataException("Error - attempt to update protected Data object.");
  }
  MAKELAZYBINSELF(right,ADD)	// for lazy + is equivalent to +=
  binaryOp(right,plus<double>());
  return (*this);
}

Data&
Data::operator+=(const boost::python::object& right)
{
  if (isProtected()) {
        throw DataException("Error - attempt to update protected Data object.");
  }
  Data tmp(right,getFunctionSpace(),false);
  MAKELAZYBINSELF(tmp,ADD)
  binaryOp(tmp,plus<double>());
  return (*this);
}

// Hmmm, operator= makes a deep copy but the copy constructor does not?
Data&
Data::operator=(const Data& other)
{
  copy(other);
  return (*this);
}

Data&
Data::operator-=(const Data& right)
{
  if (isProtected()) {
        throw DataException("Error - attempt to update protected Data object.");
  }
  MAKELAZYBINSELF(right,SUB)
  binaryOp(right,minus<double>());
  return (*this);
}

Data&
Data::operator-=(const boost::python::object& right)
{
  if (isProtected()) {
        throw DataException("Error - attempt to update protected Data object.");
  }
  Data tmp(right,getFunctionSpace(),false);
  MAKELAZYBINSELF(tmp,SUB)
  binaryOp(tmp,minus<double>());
  return (*this);
}

Data&
Data::operator*=(const Data& right)
{
  if (isProtected()) {
        throw DataException("Error - attempt to update protected Data object.");
  }
  MAKELAZYBINSELF(right,MUL)
  binaryOp(right,multiplies<double>());
  return (*this);
}

Data&
Data::operator*=(const boost::python::object& right)
{  
  if (isProtected()) {
        throw DataException("Error - attempt to update protected Data object.");
  }
  Data tmp(right,getFunctionSpace(),false);
  MAKELAZYBINSELF(tmp,MUL)
  binaryOp(tmp,multiplies<double>());
  return (*this);
}

Data&
Data::operator/=(const Data& right)
{
  if (isProtected()) {
        throw DataException("Error - attempt to update protected Data object.");
  }
  MAKELAZYBINSELF(right,DIV)
  binaryOp(right,divides<double>());
  return (*this);
}

Data&
Data::operator/=(const boost::python::object& right)
{
  if (isProtected()) {
        throw DataException("Error - attempt to update protected Data object.");
  }
  Data tmp(right,getFunctionSpace(),false);
  MAKELAZYBINSELF(tmp,DIV)
  binaryOp(tmp,divides<double>());
  return (*this);
}

Data
Data::rpowO(const boost::python::object& left) const
{
  Data left_d(left,*this);
  return left_d.powD(*this);
}

Data
Data::powO(const boost::python::object& right) const
{
  Data tmp(right,getFunctionSpace(),false);
  return powD(tmp);
}

Data
Data::powD(const Data& right) const
{
  MAKELAZYBIN(right,POW)
  return C_TensorBinaryOperation<double (*)(double, double)>(*this, right, ::pow);
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator+(const Data& left, const Data& right)
{
  MAKELAZYBIN2(left,right,ADD)
  return C_TensorBinaryOperation(left, right, plus<double>());
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator-(const Data& left, const Data& right)
{
  MAKELAZYBIN2(left,right,SUB)
  return C_TensorBinaryOperation(left, right, minus<double>());
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator*(const Data& left, const Data& right)
{
  MAKELAZYBIN2(left,right,MUL)
  return C_TensorBinaryOperation(left, right, multiplies<double>());
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator/(const Data& left, const Data& right)
{
  MAKELAZYBIN2(left,right,DIV)
  return C_TensorBinaryOperation(left, right, divides<double>());
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator+(const Data& left, const boost::python::object& right)
{
  Data tmp(right,left.getFunctionSpace(),false);
  MAKELAZYBIN2(left,tmp,ADD)
  return left+tmp;
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator-(const Data& left, const boost::python::object& right)
{
  Data tmp(right,left.getFunctionSpace(),false);
  MAKELAZYBIN2(left,tmp,SUB)
  return left-tmp;
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator*(const Data& left, const boost::python::object& right)
{
  Data tmp(right,left.getFunctionSpace(),false);
  MAKELAZYBIN2(left,tmp,MUL)
  return left*tmp;
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator/(const Data& left, const boost::python::object& right)
{
  Data tmp(right,left.getFunctionSpace(),false);
  MAKELAZYBIN2(left,tmp,DIV)
  return left/tmp;
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator+(const boost::python::object& left, const Data& right)
{
  Data tmp(left,right.getFunctionSpace(),false);
  MAKELAZYBIN2(tmp,right,ADD)
  return tmp+right;
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator-(const boost::python::object& left, const Data& right)
{
  Data tmp(left,right.getFunctionSpace(),false);
  MAKELAZYBIN2(tmp,right,SUB)
  return tmp-right;
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator*(const boost::python::object& left, const Data& right)
{
  Data tmp(left,right.getFunctionSpace(),false);
  MAKELAZYBIN2(tmp,right,MUL)
  return tmp*right;
}

//
// NOTE: It is essential to specify the namespace this operator belongs to
Data
escript::operator/(const boost::python::object& left, const Data& right)
{
  Data tmp(left,right.getFunctionSpace(),false);
  MAKELAZYBIN2(tmp,right,DIV)
  return tmp/right;
}


/* TODO */
/* global reduction */
Data
Data::getItem(const boost::python::object& key) const
{

  DataTypes::RegionType slice_region=DataTypes::getSliceRegion(getDataPointShape(),key);

  if (slice_region.size()!=getDataPointRank()) {
    throw DataException("Error - slice size does not match Data rank.");
  }

  return getSlice(slice_region);
}

/* TODO */
/* global reduction */
Data
Data::getSlice(const DataTypes::RegionType& region) const
{
  return Data(*this,region);
}

/* TODO */
/* global reduction */
void
Data::setItemO(const boost::python::object& key,
               const boost::python::object& value)
{
  Data tempData(value,getFunctionSpace());
  setItemD(key,tempData);
}

void
Data::setItemD(const boost::python::object& key,
               const Data& value)
{
//  const DataArrayView& view=getPointDataView();

  DataTypes::RegionType slice_region=DataTypes::getSliceRegion(getDataPointShape(),key);
  if (slice_region.size()!=getDataPointRank()) {
    throw DataException("Error - slice size does not match Data rank.");
  }
  if (getFunctionSpace()!=value.getFunctionSpace()) {
     setSlice(Data(value,getFunctionSpace()),slice_region);
  } else {
     setSlice(value,slice_region);
  }
}

void
Data::setSlice(const Data& value,
               const DataTypes::RegionType& region)
{
  if (isProtected()) {
        throw DataException("Error - attempt to update protected Data object.");
  }
  FORCERESOLVE;
/*  if (isLazy())
  {
	throw DataException("Error - setSlice not permitted on lazy data.");
  }*/
  Data tempValue(value);
  typeMatchLeft(tempValue);
  typeMatchRight(tempValue);
  getReady()->setSlice(tempValue.m_data.get(),region);
}

void
Data::typeMatchLeft(Data& right) const
{
  if (right.isLazy() && !isLazy())
  {
    right.resolve();
  }
  if (isExpanded()){
    right.expand();
  } else if (isTagged()) {
    if (right.isConstant()) {
      right.tag();
    }
  }
}

void
Data::typeMatchRight(const Data& right)
{
  if (isLazy() && !right.isLazy())
  {
    resolve();
  }
  if (isTagged()) {
    if (right.isExpanded()) {
      expand();
    }
  } else if (isConstant()) {
    if (right.isExpanded()) {
      expand();
    } else if (right.isTagged()) {
      tag();
    }
  }
}

void
Data::setTaggedValueByName(std::string name,
                           const boost::python::object& value)
{
     if (getFunctionSpace().getDomain()->isValidTagName(name)) {
	FORCERESOLVE;
        int tagKey=getFunctionSpace().getDomain()->getTag(name);
        setTaggedValue(tagKey,value);
     }
}
void
Data::setTaggedValue(int tagKey,
                     const boost::python::object& value)
{
  if (isProtected()) {
        throw DataException("Error - attempt to update protected Data object.");
  }
  //
  // Ensure underlying data object is of type DataTagged
  FORCERESOLVE;
  if (isConstant()) tag();
  numeric::array asNumArray(value);

  // extract the shape of the numarray
  DataTypes::ShapeType tempShape;
  for (int i=0; i < asNumArray.getrank(); i++) {
    tempShape.push_back(extract<int>(asNumArray.getshape()[i]));
  }

  DataVector temp_data2;
  temp_data2.copyFromNumArray(asNumArray,1);

  m_data->setTaggedValue(tagKey,tempShape, temp_data2);
}


void
Data::setTaggedValueFromCPP(int tagKey,
			    const DataTypes::ShapeType& pointshape,
                            const DataTypes::ValueType& value,
			    int dataOffset)
{
  if (isProtected()) {
        throw DataException("Error - attempt to update protected Data object.");
  }
  //
  // Ensure underlying data object is of type DataTagged
  FORCERESOLVE;
  if (isConstant()) tag();
  //
  // Call DataAbstract::setTaggedValue
  m_data->setTaggedValue(tagKey,pointshape, value, dataOffset);
}

int
Data::getTagNumber(int dpno)
{
  if (isEmpty())
  {
	throw DataException("Error - operation not permitted on instances of DataEmpty.");
  }
  return getFunctionSpace().getTagFromDataPointNo(dpno);
}


ostream& escript::operator<<(ostream& o, const Data& data)
{
  o << data.toString();
  return o;
}

Data
escript::C_GeneralTensorProduct(Data& arg_0,
                     Data& arg_1,
                     int axis_offset,
                     int transpose)
{
  // General tensor product: res(SL x SR) = arg_0(SL x SM) * arg_1(SM x SR)
  // SM is the product of the last axis_offset entries in arg_0.getShape().

  // deal with any lazy data
//   if (arg_0.isLazy()) {arg_0.resolve();}
//   if (arg_1.isLazy()) {arg_1.resolve();}
  if (arg_0.isLazy() || arg_1.isLazy() || (AUTOLAZYON && (arg_0.isExpanded() || arg_1.isExpanded())))
  {
	DataLazy* c=new DataLazy(arg_0.borrowDataPtr(), arg_1.borrowDataPtr(), PROD, axis_offset,transpose);
	return Data(c);
  }

  // Interpolate if necessary and find an appropriate function space
  Data arg_0_Z, arg_1_Z;
  if (arg_0.getFunctionSpace()!=arg_1.getFunctionSpace()) {
    if (arg_0.probeInterpolation(arg_1.getFunctionSpace())) {
      arg_0_Z = arg_0.interpolate(arg_1.getFunctionSpace());
      arg_1_Z = Data(arg_1);
    }
    else if (arg_1.probeInterpolation(arg_0.getFunctionSpace())) {
      arg_1_Z=arg_1.interpolate(arg_0.getFunctionSpace());
      arg_0_Z =Data(arg_0);
    }
    else {
      throw DataException("Error - C_GeneralTensorProduct: arguments have incompatible function spaces.");
    }
  } else {
      arg_0_Z = Data(arg_0);
      arg_1_Z = Data(arg_1);
  }
  // Get rank and shape of inputs
  int rank0 = arg_0_Z.getDataPointRank();
  int rank1 = arg_1_Z.getDataPointRank();
  const DataTypes::ShapeType& shape0 = arg_0_Z.getDataPointShape();
  const DataTypes::ShapeType& shape1 = arg_1_Z.getDataPointShape();

  // Prepare for the loops of the product and verify compatibility of shapes
  int start0=0, start1=0;
  if (transpose == 0)		{}
  else if (transpose == 1)	{ start0 = axis_offset; }
  else if (transpose == 2)	{ start1 = rank1-axis_offset; }
  else				{ throw DataException("C_GeneralTensorProduct: Error - transpose should be 0, 1 or 2"); }


  // Adjust the shapes for transpose
  DataTypes::ShapeType tmpShape0(rank0);	// pre-sizing the vectors rather
  DataTypes::ShapeType tmpShape1(rank1);	// than using push_back
  for (int i=0; i<rank0; i++)	{ tmpShape0[i]=shape0[(i+start0)%rank0]; }
  for (int i=0; i<rank1; i++)	{ tmpShape1[i]=shape1[(i+start1)%rank1]; }

#if 0
  // For debugging: show shape after transpose
  char tmp[100];
  std::string shapeStr;
  shapeStr = "(";
  for (int i=0; i<rank0; i++)	{ sprintf(tmp, "%d,", tmpShape0[i]); shapeStr += tmp; }
  shapeStr += ")";
  cout << "C_GeneralTensorProduct: Shape of arg0 is " << shapeStr << endl;
  shapeStr = "(";
  for (int i=0; i<rank1; i++)	{ sprintf(tmp, "%d,", tmpShape1[i]); shapeStr += tmp; }
  shapeStr += ")";
  cout << "C_GeneralTensorProduct: Shape of arg1 is " << shapeStr << endl;
#endif

  // Prepare for the loops of the product
  int SL=1, SM=1, SR=1;
  for (int i=0; i<rank0-axis_offset; i++)	{
    SL *= tmpShape0[i];
  }
  for (int i=rank0-axis_offset; i<rank0; i++)	{
    if (tmpShape0[i] != tmpShape1[i-(rank0-axis_offset)]) {
      throw DataException("C_GeneralTensorProduct: Error - incompatible shapes");
    }
    SM *= tmpShape0[i];
  }
  for (int i=axis_offset; i<rank1; i++)		{
    SR *= tmpShape1[i];
  }

  // Define the shape of the output (rank of shape is the sum of the loop ranges below)
  DataTypes::ShapeType shape2(rank0+rank1-2*axis_offset);	
  {			// block to limit the scope of out_index
     int out_index=0;
     for (int i=0; i<rank0-axis_offset; i++, ++out_index) { shape2[out_index]=tmpShape0[i]; } // First part of arg_0_Z
     for (int i=axis_offset; i<rank1; i++, ++out_index)   { shape2[out_index]=tmpShape1[i]; } // Last part of arg_1_Z
  }

  if (shape2.size()>ESCRIPT_MAX_DATA_RANK)
  {
     ostringstream os;
     os << "C_GeneralTensorProduct: Error - Attempt to create a rank " << shape2.size() << " object. The maximum rank is " << ESCRIPT_MAX_DATA_RANK << ".";
     throw DataException(os.str());
  }

  // Declare output Data object
  Data res;

  if      (arg_0_Z.isConstant()   && arg_1_Z.isConstant()) {
    res = Data(0.0, shape2, arg_1_Z.getFunctionSpace());	// DataConstant output
    double *ptr_0 = &(arg_0_Z.getDataAtOffset(0));
    double *ptr_1 = &(arg_1_Z.getDataAtOffset(0));
    double *ptr_2 = &(res.getDataAtOffset(0));
    matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
  }
  else if (arg_0_Z.isConstant()   && arg_1_Z.isTagged()) {

    // Prepare the DataConstant input
    DataConstant* tmp_0=dynamic_cast<DataConstant*>(arg_0_Z.borrowData());
    if (tmp_0==0) { throw DataException("GTP Programming error - casting to DataConstant."); }

    // Borrow DataTagged input from Data object
    DataTagged* tmp_1=dynamic_cast<DataTagged*>(arg_1_Z.borrowData());
    if (tmp_1==0) { throw DataException("GTP_1 Programming error - casting to DataTagged."); }

    // Prepare a DataTagged output 2
    res = Data(0.0, shape2, arg_1_Z.getFunctionSpace());	// DataTagged output
    res.tag();
    DataTagged* tmp_2=dynamic_cast<DataTagged*>(res.borrowData());
    if (tmp_2==0) { throw DataException("GTP Programming error - casting to DataTagged."); }

    // Prepare offset into DataConstant
    int offset_0 = tmp_0->getPointOffset(0,0);
    double *ptr_0 = &(arg_0_Z.getDataAtOffset(offset_0));
    // Get the views
//     DataArrayView view_1 = tmp_1->getDefaultValue();
//     DataArrayView view_2 = tmp_2->getDefaultValue();
//     // Get the pointers to the actual data
//     double *ptr_1 = &((view_1.getData())[0]);
//     double *ptr_2 = &((view_2.getData())[0]);

    double *ptr_1 = &(tmp_1->getDefaultValue(0));
    double *ptr_2 = &(tmp_2->getDefaultValue(0));


    // Compute an MVP for the default
    matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
    // Compute an MVP for each tag
    const DataTagged::DataMapType& lookup_1=tmp_1->getTagLookup();
    DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
    for (i=lookup_1.begin();i!=lookup_1.end();i++) {
      tmp_2->addTag(i->first);
//       DataArrayView view_1 = tmp_1->getDataPointByTag(i->first);
//       DataArrayView view_2 = tmp_2->getDataPointByTag(i->first);
//       double *ptr_1 = &view_1.getData(0);
//       double *ptr_2 = &view_2.getData(0);

      double *ptr_1 = &(tmp_1->getDataByTag(i->first,0));
      double *ptr_2 = &(tmp_2->getDataByTag(i->first,0));
	
      matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
    }

  }
  else if (arg_0_Z.isConstant()   && arg_1_Z.isExpanded()) {

    res = Data(0.0, shape2, arg_1_Z.getFunctionSpace(),true); // DataExpanded output
    DataConstant* tmp_0=dynamic_cast<DataConstant*>(arg_0_Z.borrowData());
    DataExpanded* tmp_1=dynamic_cast<DataExpanded*>(arg_1_Z.borrowData());
    DataExpanded* tmp_2=dynamic_cast<DataExpanded*>(res.borrowData());
    if (tmp_0==0) { throw DataException("GTP Programming error - casting to DataConstant."); }
    if (tmp_1==0) { throw DataException("GTP Programming error - casting to DataExpanded."); }
    if (tmp_2==0) { throw DataException("GTP Programming error - casting to DataExpanded."); }
    int sampleNo_1,dataPointNo_1;
    int numSamples_1 = arg_1_Z.getNumSamples();
    int numDataPointsPerSample_1 = arg_1_Z.getNumDataPointsPerSample();
    int offset_0 = tmp_0->getPointOffset(0,0);
    #pragma omp parallel for private(sampleNo_1,dataPointNo_1) schedule(static)
    for (sampleNo_1 = 0; sampleNo_1 < numSamples_1; sampleNo_1++) {
      for (dataPointNo_1 = 0; dataPointNo_1 < numDataPointsPerSample_1; dataPointNo_1++) {
        int offset_1 = tmp_1->getPointOffset(sampleNo_1,dataPointNo_1);
        int offset_2 = tmp_2->getPointOffset(sampleNo_1,dataPointNo_1);
        double *ptr_0 = &(arg_0_Z.getDataAtOffset(offset_0));
        double *ptr_1 = &(arg_1_Z.getDataAtOffset(offset_1));
        double *ptr_2 = &(res.getDataAtOffset(offset_2));
        matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
      }
    }

  }
  else if (arg_0_Z.isTagged()     && arg_1_Z.isConstant()) {

    // Borrow DataTagged input from Data object
    DataTagged* tmp_0=dynamic_cast<DataTagged*>(arg_0_Z.borrowData());
    if (tmp_0==0) { throw DataException("GTP_0 Programming error - casting to DataTagged."); }

    // Prepare the DataConstant input
    DataConstant* tmp_1=dynamic_cast<DataConstant*>(arg_1_Z.borrowData());
    if (tmp_1==0) { throw DataException("GTP Programming error - casting to DataConstant."); }

    // Prepare a DataTagged output 2
    res = Data(0.0, shape2, arg_0_Z.getFunctionSpace());	// DataTagged output
    res.tag();
    DataTagged* tmp_2=dynamic_cast<DataTagged*>(res.borrowData());
    if (tmp_2==0) { throw DataException("GTP Programming error - casting to DataTagged."); }

    // Prepare offset into DataConstant
    int offset_1 = tmp_1->getPointOffset(0,0);
    double *ptr_1 = &(arg_1_Z.getDataAtOffset(offset_1));
    // Get the views
//     DataArrayView view_0 = tmp_0->getDefaultValue();
//     DataArrayView view_2 = tmp_2->getDefaultValue();
//     // Get the pointers to the actual data
//     double *ptr_0 = &((view_0.getData())[0]);
//     double *ptr_2 = &((view_2.getData())[0]);

    double *ptr_0 = &(tmp_0->getDefaultValue(0));
    double *ptr_2 = &(tmp_2->getDefaultValue(0));

    // Compute an MVP for the default
    matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
    // Compute an MVP for each tag
    const DataTagged::DataMapType& lookup_0=tmp_0->getTagLookup();
    DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
    for (i=lookup_0.begin();i!=lookup_0.end();i++) {
//      tmp_2->addTaggedValue(i->first,tmp_2->getDefaultValue());
//       DataArrayView view_0 = tmp_0->getDataPointByTag(i->first);
//       DataArrayView view_2 = tmp_2->getDataPointByTag(i->first);
//       double *ptr_0 = &view_0.getData(0);
//       double *ptr_2 = &view_2.getData(0);

      tmp_2->addTag(i->first);
      double *ptr_0 = &(tmp_0->getDataByTag(i->first,0));
      double *ptr_2 = &(tmp_2->getDataByTag(i->first,0));
      matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
    }

  }
  else if (arg_0_Z.isTagged()     && arg_1_Z.isTagged()) {

    // Borrow DataTagged input from Data object
    DataTagged* tmp_0=dynamic_cast<DataTagged*>(arg_0_Z.borrowData());
    if (tmp_0==0) { throw DataException("GTP Programming error - casting to DataTagged."); }

    // Borrow DataTagged input from Data object
    DataTagged* tmp_1=dynamic_cast<DataTagged*>(arg_1_Z.borrowData());
    if (tmp_1==0) { throw DataException("GTP Programming error - casting to DataTagged."); }

    // Prepare a DataTagged output 2
    res = Data(0.0, shape2, arg_1_Z.getFunctionSpace());
    res.tag();	// DataTagged output
    DataTagged* tmp_2=dynamic_cast<DataTagged*>(res.borrowData());
    if (tmp_2==0) { throw DataException("GTP Programming error - casting to DataTagged."); }

//     // Get the views
//     DataArrayView view_0 = tmp_0->getDefaultValue();
//     DataArrayView view_1 = tmp_1->getDefaultValue();
//     DataArrayView view_2 = tmp_2->getDefaultValue();
//     // Get the pointers to the actual data
//     double *ptr_0 = &((view_0.getData())[0]);
//     double *ptr_1 = &((view_1.getData())[0]);
//     double *ptr_2 = &((view_2.getData())[0]);

    double *ptr_0 = &(tmp_0->getDefaultValue(0));
    double *ptr_1 = &(tmp_1->getDefaultValue(0));
    double *ptr_2 = &(tmp_2->getDefaultValue(0));


    // Compute an MVP for the default
    matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
    // Merge the tags
    DataTagged::DataMapType::const_iterator i; // i->first is a tag, i->second is an offset into memory
    const DataTagged::DataMapType& lookup_0=tmp_0->getTagLookup();
    const DataTagged::DataMapType& lookup_1=tmp_1->getTagLookup();
    for (i=lookup_0.begin();i!=lookup_0.end();i++) {
      tmp_2->addTag(i->first); // use tmp_2 to get correct shape
    }
    for (i=lookup_1.begin();i!=lookup_1.end();i++) {
      tmp_2->addTag(i->first);
    }
    // Compute an MVP for each tag
    const DataTagged::DataMapType& lookup_2=tmp_2->getTagLookup();
    for (i=lookup_2.begin();i!=lookup_2.end();i++) {
//       DataArrayView view_0 = tmp_0->getDataPointByTag(i->first);
//       DataArrayView view_1 = tmp_1->getDataPointByTag(i->first);
//       DataArrayView view_2 = tmp_2->getDataPointByTag(i->first);
//       double *ptr_0 = &view_0.getData(0);
//       double *ptr_1 = &view_1.getData(0);
//       double *ptr_2 = &view_2.getData(0);

      double *ptr_0 = &(tmp_0->getDataByTag(i->first,0));
      double *ptr_1 = &(tmp_1->getDataByTag(i->first,0));
      double *ptr_2 = &(tmp_2->getDataByTag(i->first,0));

      matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
    }

  }
  else if (arg_0_Z.isTagged()     && arg_1_Z.isExpanded()) {

    // After finding a common function space above the two inputs have the same numSamples and num DPPS
    res = Data(0.0, shape2, arg_1_Z.getFunctionSpace(),true); // DataExpanded output
    DataTagged*   tmp_0=dynamic_cast<DataTagged*>(arg_0_Z.borrowData());
    DataExpanded* tmp_1=dynamic_cast<DataExpanded*>(arg_1_Z.borrowData());
    DataExpanded* tmp_2=dynamic_cast<DataExpanded*>(res.borrowData());
    if (tmp_0==0) { throw DataException("GTP Programming error - casting to DataTagged."); }
    if (tmp_1==0) { throw DataException("GTP Programming error - casting to DataExpanded."); }
    if (tmp_2==0) { throw DataException("GTP Programming error - casting to DataExpanded."); }
    int sampleNo_0,dataPointNo_0;
    int numSamples_0 = arg_0_Z.getNumSamples();
    int numDataPointsPerSample_0 = arg_0_Z.getNumDataPointsPerSample();
    #pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
    for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
      int offset_0 = tmp_0->getPointOffset(sampleNo_0,0); // They're all the same, so just use #0
      double *ptr_0 = &(arg_0_Z.getDataAtOffset(offset_0));
      for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
        int offset_1 = tmp_1->getPointOffset(sampleNo_0,dataPointNo_0);
        int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
        double *ptr_1 = &(arg_1_Z.getDataAtOffset(offset_1));
        double *ptr_2 = &(res.getDataAtOffset(offset_2));
        matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
      }
    }

  }
  else if (arg_0_Z.isExpanded()   && arg_1_Z.isConstant()) {

    res = Data(0.0, shape2, arg_1_Z.getFunctionSpace(),true); // DataExpanded output
    DataExpanded* tmp_0=dynamic_cast<DataExpanded*>(arg_0_Z.borrowData());
    DataConstant* tmp_1=dynamic_cast<DataConstant*>(arg_1_Z.borrowData());
    DataExpanded* tmp_2=dynamic_cast<DataExpanded*>(res.borrowData());
    if (tmp_0==0) { throw DataException("GTP Programming error - casting to DataExpanded."); }
    if (tmp_1==0) { throw DataException("GTP Programming error - casting to DataConstant."); }
    if (tmp_2==0) { throw DataException("GTP Programming error - casting to DataExpanded."); }
    int sampleNo_0,dataPointNo_0;
    int numSamples_0 = arg_0_Z.getNumSamples();
    int numDataPointsPerSample_0 = arg_0_Z.getNumDataPointsPerSample();
    int offset_1 = tmp_1->getPointOffset(0,0);
    #pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
    for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
      for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
        int offset_0 = tmp_0->getPointOffset(sampleNo_0,dataPointNo_0);
        int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
        double *ptr_0 = &(arg_0_Z.getDataAtOffset(offset_0));
        double *ptr_1 = &(arg_1_Z.getDataAtOffset(offset_1));
        double *ptr_2 = &(res.getDataAtOffset(offset_2));
        matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
      }
    }


  }
  else if (arg_0_Z.isExpanded()   && arg_1_Z.isTagged()) {

    // After finding a common function space above the two inputs have the same numSamples and num DPPS
    res = Data(0.0, shape2, arg_1_Z.getFunctionSpace(),true); // DataExpanded output
    DataExpanded* tmp_0=dynamic_cast<DataExpanded*>(arg_0_Z.borrowData());
    DataTagged*   tmp_1=dynamic_cast<DataTagged*>(arg_1_Z.borrowData());
    DataExpanded* tmp_2=dynamic_cast<DataExpanded*>(res.borrowData());
    if (tmp_0==0) { throw DataException("GTP Programming error - casting to DataExpanded."); }
    if (tmp_1==0) { throw DataException("GTP Programming error - casting to DataTagged."); }
    if (tmp_2==0) { throw DataException("GTP Programming error - casting to DataExpanded."); }
    int sampleNo_0,dataPointNo_0;
    int numSamples_0 = arg_0_Z.getNumSamples();
    int numDataPointsPerSample_0 = arg_0_Z.getNumDataPointsPerSample();
    #pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
    for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
      int offset_1 = tmp_1->getPointOffset(sampleNo_0,0);
      double *ptr_1 = &(arg_1_Z.getDataAtOffset(offset_1));
      for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
        int offset_0 = tmp_0->getPointOffset(sampleNo_0,dataPointNo_0);
        int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
        double *ptr_0 = &(arg_0_Z.getDataAtOffset(offset_0));
        double *ptr_2 = &(res.getDataAtOffset(offset_2));
        matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
      }
    }

  }
  else if (arg_0_Z.isExpanded()   && arg_1_Z.isExpanded()) {

    // After finding a common function space above the two inputs have the same numSamples and num DPPS
    res = Data(0.0, shape2, arg_1_Z.getFunctionSpace(),true); // DataExpanded output
    DataExpanded* tmp_0=dynamic_cast<DataExpanded*>(arg_0_Z.borrowData());
    DataExpanded* tmp_1=dynamic_cast<DataExpanded*>(arg_1_Z.borrowData());
    DataExpanded* tmp_2=dynamic_cast<DataExpanded*>(res.borrowData());
    if (tmp_0==0) { throw DataException("GTP Programming error - casting to DataExpanded."); }
    if (tmp_1==0) { throw DataException("GTP Programming error - casting to DataExpanded."); }
    if (tmp_2==0) { throw DataException("GTP Programming error - casting to DataExpanded."); }
    int sampleNo_0,dataPointNo_0;
    int numSamples_0 = arg_0_Z.getNumSamples();
    int numDataPointsPerSample_0 = arg_0_Z.getNumDataPointsPerSample();
    #pragma omp parallel for private(sampleNo_0,dataPointNo_0) schedule(static)
    for (sampleNo_0 = 0; sampleNo_0 < numSamples_0; sampleNo_0++) {
      for (dataPointNo_0 = 0; dataPointNo_0 < numDataPointsPerSample_0; dataPointNo_0++) {
        int offset_0 = tmp_0->getPointOffset(sampleNo_0,dataPointNo_0);
        int offset_1 = tmp_1->getPointOffset(sampleNo_0,dataPointNo_0);
        int offset_2 = tmp_2->getPointOffset(sampleNo_0,dataPointNo_0);
        double *ptr_0 = &(arg_0_Z.getDataAtOffset(offset_0));
        double *ptr_1 = &(arg_1_Z.getDataAtOffset(offset_1));
        double *ptr_2 = &(res.getDataAtOffset(offset_2));
        matrix_matrix_product(SL, SM, SR, ptr_0, ptr_1, ptr_2, transpose);
      }
    }

  }
  else {
    throw DataException("Error - C_GeneralTensorProduct: unknown combination of inputs");
  }

  return res;
}

DataAbstract*
Data::borrowData() const
{
  return m_data.get();
}

// Not all that happy about returning a non-const from a const
DataAbstract_ptr
Data::borrowDataPtr() const
{
  return m_data;
}

// Not all that happy about returning a non-const from a const
DataReady_ptr
Data::borrowReadyPtr() const
{
   DataReady_ptr dr=dynamic_pointer_cast<DataReady>(m_data);
   EsysAssert((dr!=0), "Error - casting to DataReady.");
   return dr;
}

std::string
Data::toString() const
{
    if (!m_data->isEmpty() &&
	!m_data->isLazy() && 
	getLength()>escriptParams.getInt("TOO_MANY_LINES"))
    {
	stringstream temp;
	temp << "Summary: inf="<< inf_const() << " sup=" << sup_const() << " data points=" << getNumDataPoints();
	return  temp.str();
    }
    return m_data->toString();
}



DataTypes::ValueType::const_reference
Data::getDataAtOffset(DataTypes::ValueType::size_type i) const
{
    if (isLazy())
    {
	throw DataException("Programmer error - getDataAtOffset not permitted on lazy data (object is const which prevents resolving).");
    }
    return getReady()->getDataAtOffset(i);
}


DataTypes::ValueType::reference
Data::getDataAtOffset(DataTypes::ValueType::size_type i)
{
//     if (isLazy())
//     {
// 	throw DataException("getDataAtOffset not permitted on lazy data.");
//     }
    FORCERESOLVE;
    return getReady()->getDataAtOffset(i);
}

DataTypes::ValueType::const_reference
Data::getDataPoint(int sampleNo, int dataPointNo) const
{
  if (!isReady())
  {
	throw DataException("Programmer error - getDataPoint() not permitted on Lazy Data (object is const which prevents resolving).");
  }
  else
  {
	const DataReady* dr=getReady();
	return dr->getDataAtOffset(dr->getPointOffset(sampleNo, dataPointNo));
  }
}


DataTypes::ValueType::reference
Data::getDataPoint(int sampleNo, int dataPointNo)
{
  FORCERESOLVE;
  if (!isReady())
  {
	throw DataException("Programmer error - getDataPoint() not permitted on Lazy Data.");
  }
  else
  {
	DataReady* dr=getReady();
	return dr->getDataAtOffset(dr->getPointOffset(sampleNo, dataPointNo));
  }
}


/* Member functions specific to the MPI implementation */

void
Data::print()
{
  int i,j;

  printf( "Data is %dX%d\n", getNumSamples(), getNumDataPointsPerSample() );
  for( i=0; i<getNumSamples(); i++ )
  {
    printf( "[%6d]", i );
    for( j=0; j<getNumDataPointsPerSample(); j++ )
      printf( "\t%10.7g", (getSampleData(i))[j] );
    printf( "\n" );
  }
}
void
Data::dump(const std::string fileName) const
{
  try
  {
	if (isLazy())
	{
	  Data temp(*this);	// this is to get a non-const object which we can resolve
	  temp.resolve();
	  temp.dump(fileName);
	}
	else
	{
          return m_data->dump(fileName);
	}
  }
  catch (std::exception& e)
  {
        cout << e.what() << endl;
  }
}

int
Data::get_MPISize() const
{
	int size;
#ifdef PASO_MPI
	int error;
	error = MPI_Comm_size( get_MPIComm(), &size );
#else
	size = 1;
#endif
	return size;
}

int
Data::get_MPIRank() const
{
	int rank;
#ifdef PASO_MPI
	int error;
	error = MPI_Comm_rank( get_MPIComm(), &rank );
#else
	rank = 0;
#endif
	return rank;
}

MPI_Comm
Data::get_MPIComm() const
{
#ifdef PASO_MPI
	return MPI_COMM_WORLD;
#else
	return -1;
#endif
}


