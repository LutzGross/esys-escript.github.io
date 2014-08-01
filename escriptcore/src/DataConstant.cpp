
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


#include "Data.h"
#include "DataConstant.h"
#include "DataException.h"
#include "esysUtils/EsysAssert.h"
#include "esysUtils/Esys_MPI.h"

#include <iostream>
#include <boost/python/extract.hpp>
#include <boost/scoped_ptr.hpp>
#ifdef USE_NETCDF
#include <netcdfcpp.h>
#endif

#include "DataMaths.h"

// #define CHECK_FOR_EX_WRITE if (!checkNoSharing()) {throw DataException("Attempt to modify shared object");}

#define CHECK_FOR_EX_WRITE if (!checkNoSharing()) {std::ostringstream ss; ss << " Attempt to modify shared object. line " << __LINE__ << " of " << __FILE__; ss << m_owners.size(); cerr << ss << endl; /* *((int*)0)=17; */throw DataException(ss.str());}

using namespace std;
using namespace boost::python;

namespace escript {

DataConstant::DataConstant(const WrappedArray& value,
                           const FunctionSpace& what)
  : parent(what,value.getShape())
{
  m_data.copyFromArray(value,1);
}

DataConstant::DataConstant(const DataConstant& other)
  : parent(other.getFunctionSpace(),other.getShape())
{ 
  m_data=other.m_data;
}

DataConstant::DataConstant(const DataConstant& other,
                           const DataTypes::RegionType& region)
  : parent(other.getFunctionSpace(),DataTypes::getResultSliceShape(region))
{
  //
  // allocate space for this new DataConstant's data
  int len = getNoValues();
  m_data.resize(len,0.,len);
  //
  // create a view of the data with the correct shape
  DataTypes::RegionLoopRangeType region_loop_range=DataTypes::getSliceRegionLoopRange(region);
  //
  // load the view with the data from the slice
  DataTypes::copySlice(m_data,getShape(),0,other.getVectorRO(),other.getShape(),0,region_loop_range);
}

DataConstant::DataConstant(const FunctionSpace& what,
                           const DataTypes::ShapeType &shape,
                           const DataTypes::ValueType &data)
  : parent(what,shape)
{
  //
  // copy the data in the correct format
  m_data=data;
}

DataConstant::DataConstant(const FunctionSpace& what,
                           const DataTypes::ShapeType &shape,
                           const double v)
  : parent(what,shape), m_data(DataTypes::noValues(shape),v)
{
}


bool
DataConstant::hasNaN() const
{
  bool haveNaN=false;
  #pragma omp parallel for
	for (ValueType::size_type i=0;i<m_data.size();++i)
	{
		if (nancheck(m_data[i]))	
		{
			#pragma omp atomic write
            haveNaN=true;
		}
	}
	return haveNaN;
}

void
DataConstant::replaceNaN(double value)
{
  #pragma omp parallel for
  for (ValueType::size_type i=0;i<m_data.size();++i)
  {
    if (nancheck(m_data[i]))  
    {
      m_data[i] = value;
    } 
  }
}

string
DataConstant::toString() const
{
  return DataTypes::pointToString(m_data,getShape(),0,"");
}


DataAbstract*
DataConstant::deepCopy()
{
  return new DataConstant(*this);
}


DataTypes::ValueType::size_type
DataConstant::getPointOffset(int sampleNo,
                             int dataPointNo) const
{
// We avoid this check for constant data due to issues manipulating 
// data with no samples.

//  EsysAssert((validSamplePointNo(dataPointNo) && validSampleNo(sampleNo)),
//              "Invalid index, sampleNo: " << sampleNo << " dataPointNo: " << dataPointNo);
  //
  // Whatever the coord's always return the same value as this is constant data.
  return 0;
}

DataTypes::ValueType::size_type
DataConstant::getPointOffset(int sampleNo,
                             int dataPointNo)
{
// We avoid this check for constant data due to issues manipulating 
// data with no samples.
 
//  EsysAssert((validSamplePointNo(dataPointNo) && validSampleNo(sampleNo)),
//              "Invalid index, sampleNo: " << sampleNo << " dataPointNo: " << dataPointNo);
  //
  // Whatever the coord's always return the same value as this is constant data.
  //

  return 0;
}

DataTypes::ValueType::size_type
DataConstant::getLength() const
{
  return m_data.size();
}

DataAbstract*
DataConstant::getSlice(const DataTypes::RegionType& region) const
{
  return new DataConstant(*this,region);
}

void
DataConstant::setSlice(const DataAbstract* value,
                       const DataTypes::RegionType& region)
{
  const DataConstant* tempDataConst=dynamic_cast<const DataConstant*>(value);
  if (tempDataConst==0) {
    throw DataException("Programming error - casting to DataConstant.");
  }
  CHECK_FOR_EX_WRITE
  //
  DataTypes::ShapeType shape(DataTypes::getResultSliceShape(region));
  DataTypes::RegionLoopRangeType region_loop_range=DataTypes::getSliceRegionLoopRange(region);
  //
  // check shape:
  if (getRank()!=region.size()) {
    throw DataException("Error - Invalid slice region.");
  }
  if (getRank()>0 && !DataTypes::checkShape(value->getShape(),shape)) {
    throw DataException (DataTypes::createShapeErrorMessage(
                "Error - Couldn't copy slice due to shape mismatch.",shape,value->getShape()));
  }
  //   getPointDataView().copySliceFrom(tempDataConst->getPointDataView(),region_loop_range);
  DataTypes::copySliceFrom(m_data,getShape(),0,tempDataConst->getVectorRO(), tempDataConst->getShape(),0,region_loop_range);
}



void
DataConstant::symmetric(DataAbstract* ev)
{
  DataConstant* temp_ev=dynamic_cast<DataConstant*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataConstant::symmetric: casting to DataConstant failed (probably a programming error).");
  }
  DataMaths::symmetric(m_data,getShape(),0,temp_ev->getVectorRW(), temp_ev->getShape(),0);
}

void
DataConstant::nonsymmetric(DataAbstract* ev)
{
  DataConstant* temp_ev=dynamic_cast<DataConstant*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataConstant::nonsymmetric: casting to DataConstant failed (probably a programming error).");
  }
  DataMaths::nonsymmetric(m_data,getShape(),0,temp_ev->getVectorRW(), temp_ev->getShape(),0);
}

void
DataConstant::trace(DataAbstract* ev, int axis_offset)
{
  DataConstant* temp_ev=dynamic_cast<DataConstant*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataConstant::trace: casting to DataConstant failed (probably a programming error).");
  }
  ValueType& evVec=temp_ev->getVectorRW();
  const ShapeType& evShape=temp_ev->getShape();
  DataMaths::trace(m_data,getShape(),0,evVec,evShape,0,axis_offset);
}

void
DataConstant::swapaxes(DataAbstract* ev, int axis0, int axis1)
{
  DataConstant* temp_ev=dynamic_cast<DataConstant*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataConstant::swapaxes: casting to DataConstant failed (probably a programming error).");
  }
  DataMaths::swapaxes(m_data,getShape(),0,temp_ev->getVectorRW(), temp_ev->getShape(),0,axis0,axis1);
}

void
DataConstant::transpose(DataAbstract* ev, int axis_offset)
{
  DataConstant* temp_ev=dynamic_cast<DataConstant*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataConstant::transpose: casting to DataConstant failed (probably a programming error).");
  }
  DataMaths::transpose(m_data, getShape(),0, temp_ev->getVectorRW(),temp_ev->getShape(),0,axis_offset);
}

void
DataConstant::eigenvalues(DataAbstract* ev)
{
  DataConstant* temp_ev=dynamic_cast<DataConstant*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataConstant::eigenvalues: casting to DataConstant failed (probably a programming error).");
  }
  DataMaths::eigenvalues(m_data,getShape(),0,temp_ev->getVectorRW(), temp_ev->getShape(),0);
}
void
DataConstant::eigenvalues_and_eigenvectors(DataAbstract* ev,DataAbstract* V,const double tol)
{
  DataConstant* temp_ev=dynamic_cast<DataConstant*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataConstant::eigenvalues_and_eigenvectors: casting to DataConstant failed (probably a programming error).");
  }
  DataConstant* temp_V=dynamic_cast<DataConstant*>(V);
  if (temp_V==0) {
    throw DataException("Error - DataConstant::eigenvalues_and_eigenvectors: casting to DataConstant failed (probably a programming error).");
  }
  DataMaths::eigenvalues_and_eigenvectors(m_data, getShape(),0,temp_ev->getVectorRW(), temp_ev->getShape(),0,temp_V->getVectorRW(), temp_V->getShape(),0,tol);
}



int
DataConstant::matrixInverse(DataAbstract* out) const
{
  DataConstant* temp=dynamic_cast<DataConstant*>(out);
  if (temp==0)
  {
	throw DataException("Error - DataConstant::matrixInverse: casting to DataConstant failed (probably a programming error).");
  }
  if (getRank()!=2)
  {
	throw DataException("Error - DataExpanded::matrixInverse: input must be rank 2.");
  }
  LapackInverseHelper h(getShape()[0]);
  int res=DataMaths::matrix_inverse(m_data, getShape(), 0, temp->getVectorRW(), temp->getShape(), 0, 1, h);
  return res;
}

void
DataConstant::setToZero()
{
    CHECK_FOR_EX_WRITE
    DataTypes::ValueType::size_type n=m_data.size();
    for (int i=0; i<n ;++i) m_data[i]=0.;
}

void
DataConstant::dump(const std::string fileName) const
{
   #ifdef USE_NETCDF
   const NcDim* ncdims[DataTypes::maxRank];
   NcVar* var;
   int rank = getRank();
   int type=  getFunctionSpace().getTypeCode();
   int ndims =0;
   long dims[DataTypes::maxRank];
   const double* d_ptr=&(m_data[0]);
   DataTypes::ShapeType shape = getShape();
   int mpi_iam=getFunctionSpace().getDomain()->getMPIRank();
   int mpi_num=getFunctionSpace().getDomain()->getMPISize();
#ifdef ESYS_MPI
   MPI_Status status;
#endif

#ifdef ESYS_MPI
   /* Serialize NetCDF I/O */
   if (mpi_iam>0) MPI_Recv(&ndims, 0, MPI_INT, mpi_iam-1, 81802, MPI_COMM_WORLD, &status);
#endif

   // netCDF error handler
   NcError err(NcError::verbose_nonfatal);
   // Create the file.
   const std::string newFileName(esysUtils::appendRankToFileName(fileName,
                                                            mpi_num, mpi_iam));
   NcFile dataFile(newFileName.c_str(), NcFile::Replace);
   // check if writing was successful
   if (!dataFile.is_valid())
	throw DataException("Error - DataConstant:: opening of netCDF file for output failed.");
   if (!dataFile.add_att("type_id",0) )
	throw DataException("Error - DataConstant:: appending data type to netCDF file failed.");
   if (!dataFile.add_att("rank",rank) )
	throw DataException("Error - DataConstant:: appending rank attribute to netCDF file failed.");
   if (!dataFile.add_att("function_space_type",type))
	throw DataException("Error - DataConstant:: appending function space attribute to netCDF file failed.");

   if (rank == 0) {
      if( ! (ncdims[0] = dataFile.add_dim("l", 1)) )
		throw DataException("Error - DataConstant:: appending ncdimension 0 to netCDF file failed.");
      dims[0]=1,
      ndims=1;
   } else {
       ndims=rank;
       dims[0]=shape[0];
       if (! (ncdims[0] = dataFile.add_dim("d0",shape[0])) )
		throw DataException("Error - DataConstant:: appending ncdimension 0 to netCDF file failed.");
       if ( rank >1 ) {
           dims[1]=shape[1];
           if (! (ncdims[1] = dataFile.add_dim("d1",shape[1])) )
		throw DataException("Error - DataConstant:: appending ncdimension 1 to netCDF file failed.");
       }
       if ( rank >2 ) {
           dims[2]=shape[2];
           if (! (ncdims[2] = dataFile.add_dim("d2", shape[2])) )
		throw DataException("Error - DataConstant:: appending ncdimension 2 to netCDF file failed.");
       }
       if ( rank >3 ) {
           dims[3]=shape[3];
           if (! (ncdims[3] = dataFile.add_dim("d3", shape[3])) )
		throw DataException("Error - DataConstant:: appending ncdimension 3 to netCDF file failed.");
       }
   }

   if (! ( var = dataFile.add_var("data", ncDouble, ndims, ncdims)) )
	throw DataException("Error - DataConstant:: appending variable to netCDF file failed.");
   if (! (var->put(d_ptr,dims)) )
         throw DataException("Error - DataConstant:: copy data to netCDF buffer failed.");
#ifdef ESYS_MPI
   if (mpi_iam<mpi_num-1) MPI_Send(&ndims, 0, MPI_INT, mpi_iam+1, 81802, MPI_COMM_WORLD);
#endif
   #else
   throw DataException("Error - DataConstant:: dump is not configured with netCDF. Please contact your installation manager.");
   #endif
}

// These used to be marked as inline in DataConstant.
// But they are marked virtual in DataReady
DataTypes::ValueType&
DataConstant::getVectorRW()
{
  CHECK_FOR_EX_WRITE
  return m_data;
}

const DataTypes::ValueType&
DataConstant::getVectorRO() const
{
  return m_data;
}

}  // end of namespace
