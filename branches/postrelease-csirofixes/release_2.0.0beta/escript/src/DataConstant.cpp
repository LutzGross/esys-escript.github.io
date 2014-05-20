
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
#include "DataConstant.h"
#include "DataException.h"
#include "esysUtils/EsysAssert.h"

#include <iostream>
#include <boost/python/extract.hpp>
#ifdef USE_NETCDF
#include <netcdfcpp.h>
#endif
#ifdef PASO_MPI
#include <mpi.h>
#endif

#include <boost/python/extract.hpp>
#include "DataMaths.h"

using namespace std;
using namespace boost::python;

namespace escript {

DataConstant::DataConstant(const boost::python::numeric::array& value,
                           const FunctionSpace& what)
  : parent(what,DataTypes::shapeFromNumArray(value))
{
//   // extract the shape of the numarray
//   DataTypes::ShapeType tempShape;
//   for (int i=0; i < value.getrank(); i++) {
//     tempShape.push_back(extract<int>(value.getshape()[i]));
//   }

  // get the space for the data vector
//   int len = getNoValues();
//   DataVector temp_data(len, 0.0, len);
//   DataArrayView temp_dataView(temp_data, tempShape);
//   temp_dataView.copy(value);

  m_data.copyFromNumArray(value,1);
  //

  // copy the data in the correct format
//   m_data=temp_data;
  //
  // create the view of the data
//   DataArrayView tempView(m_data,temp_dataView.getShape());
//   setPointDataView(tempView);
}

// DataConstant::DataConstant(const DataArrayView& value,
//                            const FunctionSpace& what)
//   : DataAbstract(what)
// {
//   //
//   // copy the data in the correct format
//   m_data=value.getData();
//   //
//   // create the view of the data
//   DataArrayView tempView(m_data,value.getShape());
//   setPointDataView(tempView);
// }

DataConstant::DataConstant(const DataConstant& other)
  : parent(other.getFunctionSpace(),other.getShape())
{  //
  // copy the data in the correct format
  m_data=other.m_data;
  //
//   // create the view of the data
//   DataArrayView tempView(m_data,other.getPointDataView().getShape());
//   setPointDataView(tempView);
}

DataConstant::DataConstant(const DataConstant& other,
                           const DataTypes::RegionType& region)
  : parent(other.getFunctionSpace(),DataTypes::getResultSliceShape(region))
{
  //
  // get the shape of the slice to copy from
//   DataTypes::ShapeType shape(DataTypes::getResultSliceShape(region));
  //
  // allocate space for this new DataConstant's data
  int len = getNoValues();
  m_data.resize(len,0.,len);
  //
  // create a view of the data with the correct shape
//   DataArrayView tempView(m_data,shape);
  DataTypes::RegionLoopRangeType region_loop_range=DataTypes::getSliceRegionLoopRange(region);
  //
  // load the view with the data from the slice
//   tempView.copySlice(other.getPointDataView(),region_loop_range);
  DataTypes::copySlice(m_data,getShape(),0,other.getVector(),other.getShape(),0,region_loop_range);
//   setPointDataView(tempView);
}

DataConstant::DataConstant(const FunctionSpace& what,
                           const DataTypes::ShapeType &shape,
                           const DataTypes::ValueType &data)
  : parent(what,shape)
{
  //
  // copy the data in the correct format
  m_data=data;
  //
  // create the view of the data
//   DataArrayView tempView(m_data,shape);
//   setPointDataView(tempView);
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
  EsysAssert((validSamplePointNo(dataPointNo) && validSampleNo(sampleNo)),
              "Invalid index, sampleNo: " << sampleNo << " dataPointNo: " << dataPointNo);
  //
  // Whatever the coord's always return the same value as this is constant data.
  return 0;
}

DataTypes::ValueType::size_type
DataConstant::getPointOffset(int sampleNo,
                             int dataPointNo)
{
  EsysAssert((validSamplePointNo(dataPointNo) && validSampleNo(sampleNo)),
              "Invalid index, sampleNo: " << sampleNo << " dataPointNo: " << dataPointNo);
  //
  // Whatever the coord's always return the same value as this is constant data.
  return 0;
}

DataTypes::ValueType::size_type
DataConstant::getLength() const
{
  return m_data.size();
}

// DataArrayView
// DataConstant::getDataPoint(int sampleNo,
//                            int dataPointNo)
// {
//   EsysAssert((validSamplePointNo(dataPointNo) && validSampleNo(sampleNo)),
//              "Invalid index, sampleNo: " << sampleNo << " dataPointNo: " << dataPointNo);
//   //
//   // Whatever the coord's always return the same value as this is constant data.
//   return getPointDataView();
// }

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
  DataTypes::copySliceFrom(m_data,getShape(),0,tempDataConst->getVector(), tempDataConst->getShape(),0,region_loop_range);
}



void
DataConstant::symmetric(DataAbstract* ev)
{
  DataConstant* temp_ev=dynamic_cast<DataConstant*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataConstant::symmetric: casting to DataConstant failed (propably a programming error).");
  }
/*  DataArrayView& thisView=getPointDataView();
  DataArrayView& evView=ev->getPointDataView();*/
  DataMaths::symmetric(m_data,getShape(),0,temp_ev->getVector(), temp_ev->getShape(),0);
}

void
DataConstant::nonsymmetric(DataAbstract* ev)
{
  DataConstant* temp_ev=dynamic_cast<DataConstant*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataConstant::nonsymmetric: casting to DataConstant failed (propably a programming error).");
  }
/*  DataArrayView& thisView=getPointDataView();
  DataArrayView& evView=ev->getPointDataView();*/
  DataMaths::nonsymmetric(m_data,getShape(),0,temp_ev->getVector(), temp_ev->getShape(),0);
}

void
DataConstant::trace(DataAbstract* ev, int axis_offset)
{
  DataConstant* temp_ev=dynamic_cast<DataConstant*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataConstant::trace: casting to DataConstant failed (propably a programming error).");
  }
/*  DataArrayView& thisView=getPointDataView();
  DataArrayView& evView=ev->getPointDataView();*/
  ValueType& evVec=temp_ev->getVector();
  const ShapeType& evShape=temp_ev->getShape();
  DataMaths::trace(m_data,getShape(),0,evVec,evShape,0,axis_offset);
}

void
DataConstant::swapaxes(DataAbstract* ev, int axis0, int axis1)
{
  DataConstant* temp_ev=dynamic_cast<DataConstant*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataConstant::swapaxes: casting to DataConstant failed (propably a programming error).");
  }
//   DataArrayView& thisView=getPointDataView();
//   DataArrayView& evView=ev->getPointDataView();
  DataMaths::swapaxes(m_data,getShape(),0,temp_ev->getVector(), temp_ev->getShape(),0,axis0,axis1);
}

void
DataConstant::transpose(DataAbstract* ev, int axis_offset)
{
  DataConstant* temp_ev=dynamic_cast<DataConstant*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataConstant::transpose: casting to DataConstant failed (propably a programming error).");
  }
/*  DataArrayView& thisView=getPointDataView();
  DataArrayView& evView=ev->getPointDataView();*/
  DataMaths::transpose(m_data, getShape(),0, temp_ev->getVector(),temp_ev->getShape(),0,axis_offset);
}

void
DataConstant::eigenvalues(DataAbstract* ev)
{
  DataConstant* temp_ev=dynamic_cast<DataConstant*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataConstant::eigenvalues: casting to DataConstant failed (propably a programming error).");
  }
/*  DataArrayView& thisView=getPointDataView();
  DataArrayView& evView=ev->getPointDataView();*/
  DataMaths::eigenvalues(m_data,getShape(),0,temp_ev->getVector(), temp_ev->getShape(),0);
}
void
DataConstant::eigenvalues_and_eigenvectors(DataAbstract* ev,DataAbstract* V,const double tol)
{
  DataConstant* temp_ev=dynamic_cast<DataConstant*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataConstant::eigenvalues_and_eigenvectors: casting to DataConstant failed (propably a programming error).");
  }
  DataConstant* temp_V=dynamic_cast<DataConstant*>(V);
  if (temp_V==0) {
    throw DataException("Error - DataConstant::eigenvalues_and_eigenvectors: casting to DataConstant failed (propably a programming error).");
  }
//   DataArrayView thisView=getPointDataView();
//   DataArrayView evView=ev->getPointDataView();
//   DataArrayView VView=V->getPointDataView();

  DataMaths::eigenvalues_and_eigenvectors(m_data, getShape(),0,temp_ev->getVector(), temp_ev->getShape(),0,temp_V->getVector(), temp_V->getShape(),0,tol);
}

void
DataConstant::setToZero()
{
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
#ifdef PASO_MPI
   MPI_Status status;
#endif

#ifdef PASO_MPI
   /* Serialize NetCDF I/O */
   if (mpi_iam>0) MPI_Recv(&ndims, 0, MPI_INT, mpi_iam-1, 81802, MPI_COMM_WORLD, &status);
#endif

   // netCDF error handler
   NcError err(NcError::verbose_nonfatal);
   // Create the file.
   char *newFileName = Escript_MPI_appendRankToFileName(fileName.c_str(), mpi_num, mpi_iam);
   NcFile dataFile(newFileName, NcFile::Replace);
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
#ifdef PASO_MPI
   if (mpi_iam<mpi_num-1) MPI_Send(&ndims, 0, MPI_INT, mpi_iam+1, 81802, MPI_COMM_WORLD);
#endif
   #else
   throw DataException("Error - DataConstant:: dump is not configured with netCDF. Please contact your installation manager.");
   #endif
}

}  // end of namespace
