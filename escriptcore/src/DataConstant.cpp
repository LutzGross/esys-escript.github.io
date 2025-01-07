
/*****************************************************************************
*
* Copyright (c) 2003-2020 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014-2017 by Centre for Geoscience Computing (GeoComp)
* Development from 2019 by School of Earth and Environmental Sciences
**
*****************************************************************************/

#include "DataConstant.h"
#include "Data.h"
#include "DataException.h"
#include "DataVectorOps.h"

#include <iostream>

#ifdef ESYS_HAVE_NETCDF
 #ifdef NETCDF4
  #include <ncDim.h>
  #include <ncVar.h>
  #include <ncFile.h>
 #else
  #include <netcdfcpp.h>
 #endif
#endif

#ifdef ESYS_HAVE_HDF5
  #include <H5Cpp.h>
#endif

#ifdef SLOWSHARECHECK
  #define CHECK_FOR_EX_WRITE if (isShared()) {\
    std::ostringstream ss;\
    ss << "Attempt to modify shared object. line " << __LINE__ << " of " << __FILE__;\
    int nn=17;\
    try\
    {\
	nn=shared_from_this().use_count();\
	ss << "use count=" << nn << "\n";\
    } catch (...)\
    {\
	ss << "Failed to get a use count\n";\
    }\
    std::cerr << ss.str() << std::endl;\
    throw DataException(ss.str());}
#else
  #define CHECK_FOR_EX_WRITE 
#endif
    
    
using namespace std;
using namespace boost::python;
using namespace escript::DataTypes;

#ifdef NETCDF4
using namespace netCDF;
#endif

namespace escript {

DataConstant::DataConstant(const WrappedArray& value,
                           const FunctionSpace& what)
  : parent(what,value.getShape())
{
  if (value.isComplex())
  {
      m_data_c.copyFromArray(value,1);
      this->m_iscompl=true;
  }
  else
  {
      m_data_r.copyFromArray(value,1);
  }
}

DataConstant::DataConstant(const DataConstant& other)
  : parent(other.getFunctionSpace(),other.getShape())
{
  this->m_iscompl=other.m_iscompl;
  if (other.isComplex()) 
  {
      m_data_c=other.m_data_c;
  }
  else
  {
      m_data_r=other.m_data_r;
  }
}

DataConstant::DataConstant(const DataConstant& other,
                           const DataTypes::RegionType& region)
  : parent(other.getFunctionSpace(),DataTypes::getResultSliceShape(region))
{
    // create a view of the data with the correct shape
    DataTypes::RegionLoopRangeType region_loop_range=DataTypes::getSliceRegionLoopRange(region);
    int len = getNoValues();
    if (other.isComplex())
    {
        // allocate space for this new DataConstant's data
        m_data_c.resize(len,0.,len);
        // load the view with the data from the slice
        DataTypes::copySlice(m_data_c,getShape(),0,other.getVectorROC(),other.getShape(),0,region_loop_range);
	m_iscompl=true; 	
    }
    else
    {
        // allocate space for this new DataConstant's data
        m_data_r.resize(len,0.,len);
        // load the view with the data from the slice
        DataTypes::copySlice(m_data_r,getShape(),0,other.getVectorRO(),other.getShape(),0,region_loop_range);
	m_iscompl=false;
    }
}

DataConstant::DataConstant(const FunctionSpace& what,
                           const DataTypes::ShapeType &shape,
                           const DataTypes::RealVectorType &data)
  : parent(what,shape)
{
    // copy the data in the correct format
    m_data_r=data;
}

DataConstant::DataConstant(const FunctionSpace& what,
                           const DataTypes::ShapeType &shape,
                           const DataTypes::CplxVectorType &data)
  : parent(what,shape)
{
    // copy the data in the correct format
    m_data_c=data;
    m_iscompl=true;        
}

DataConstant::DataConstant(const FunctionSpace& what,
                           const DataTypes::ShapeType &shape,
                           const DataTypes::real_t v)
  : parent(what,shape), m_data_r(DataTypes::noValues(shape),v)
{
}

DataConstant::DataConstant(const FunctionSpace& what,
                           const DataTypes::ShapeType &shape,
                           const DataTypes::cplx_t v)
  : parent(what,shape), m_data_c(DataTypes::noValues(shape),v)
{
    m_iscompl=true;
}

bool DataConstant::hasNaN() const
{
    bool haveNaN=false;
    if (isComplex())
    {
        #pragma omp parallel for
        for (DataTypes::CplxVectorType::size_type i=0;i<m_data_c.size();++i)
        {
            if (std::isnan(m_data_c[i].real()) || std::isnan(m_data_c[i].imag()))
            {
                #pragma omp critical 
                {
                    haveNaN=true;
                }
            }
        }
    }
    else
    {
        #pragma omp parallel for
        for (DataTypes::RealVectorType::size_type i=0;i<m_data_r.size();++i)
        {
            if (std::isnan(m_data_r[i]))
            {
                #pragma omp critical 
                {
                    haveNaN=true;
                }
            }
        }
    }
    return haveNaN;
}

void
DataConstant::replaceNaN(double value)
{
    CHECK_FOR_EX_WRITE  
    if (isComplex())
    {
      #pragma omp parallel for
      for (DataTypes::CplxVectorType::size_type i=0;i<m_data_c.size();++i)
      {
        if (std::isnan(m_data_c[i].real()) || std::isnan(m_data_c[i].imag()))  
        {
          m_data_c[i] = value;
        }
      }
    }
    else
    {
      #pragma omp parallel for
      for (DataTypes::RealVectorType::size_type i=0;i<m_data_r.size();++i)
      {
        if (std::isnan(m_data_r[i]))  
        {
          m_data_r[i] = value;
        }
      }    
    }
}

void
DataConstant::replaceNaN(DataTypes::cplx_t value)
{
    CHECK_FOR_EX_WRITE  
    if (isComplex())
    {
      #pragma omp parallel for
      for (DataTypes::CplxVectorType::size_type i=0;i<m_data_c.size();++i)
      {
        if (std::isnan(m_data_c[i].real()) || std::isnan(m_data_c[i].imag())) 
        {
          m_data_c[i] = value;
        }
      }
    }
    else
    {
      complicate();
      replaceNaN(value);
    }
}

bool DataConstant::hasInf() const
{
    bool haveNaN=false;
    if (isComplex())
    {
        #pragma omp parallel for
        for (DataTypes::CplxVectorType::size_type i=0;i<m_data_c.size();++i)
        {
            if (std::isinf(m_data_c[i].real()) || std::isinf(m_data_c[i].imag()))
            {
                #pragma omp critical 
                {
                    haveNaN=true;
                }
            }
        }
    }
    else
    {
        #pragma omp parallel for
        for (DataTypes::RealVectorType::size_type i=0;i<m_data_r.size();++i)
        {
            if (std::isinf(m_data_r[i]))
            {
                #pragma omp critical 
                {
                    haveNaN=true;
                }
            }
        }
    }
    return haveNaN;
}

void
DataConstant::replaceInf(double value)
{
    CHECK_FOR_EX_WRITE  
    if (isComplex())
    {
      #pragma omp parallel for
      for (DataTypes::CplxVectorType::size_type i=0;i<m_data_c.size();++i)
      {
        if (std::isinf(m_data_c[i].real()) || std::isinf(m_data_c[i].imag()))  
        {
          m_data_c[i] = value;
        }
      }
    }
    else
    {
      #pragma omp parallel for
      for (DataTypes::RealVectorType::size_type i=0;i<m_data_r.size();++i)
      {
        if (std::isinf(m_data_r[i]))  
        {
          m_data_r[i] = value;
        }
      }    
    }
}

void
DataConstant::replaceInf(DataTypes::cplx_t value)
{
    CHECK_FOR_EX_WRITE  
    if (isComplex())
    {
      #pragma omp parallel for
      for (DataTypes::CplxVectorType::size_type i=0;i<m_data_c.size();++i)
      {
        if (std::isinf(m_data_c[i].real()) || std::isinf(m_data_c[i].imag())) 
        {
          m_data_c[i] = value;
        }
      }
    }
    else
    {
      complicate();
      replaceInf(value);
    }
}


string DataConstant::toString() const
{
    if (isComplex())
        return DataTypes::pointToString(m_data_c,getShape(),0,"");

    return DataTypes::pointToString(m_data_r,getShape(),0,"");
}


DataAbstract* DataConstant::deepCopy() const
{
    return new DataConstant(*this);
}


DataAbstract* DataConstant::zeroedCopy() const
{
    DataConstant* p=0;
    if (isComplex())
    {
        p=new DataConstant(this->getFunctionSpace(), this->getShape(), DataTypes::cplx_t(0));
    }
    else
    {
        p=new DataConstant(this->getFunctionSpace(), this->getShape(), DataTypes::real_t(0));
    }
    return p;
}



DataTypes::RealVectorType::size_type
DataConstant::getPointOffset(int sampleNo,
                             int dataPointNo) const
{
// We avoid this check for constant data due to issues manipulating 
// data with no samples.

//  ESYS_ASSERT((validSamplePointNo(dataPointNo) && validSampleNo(sampleNo)),
//              "Invalid index, sampleNo: " << sampleNo << " dataPointNo: " << dataPointNo);
  //
  // Whatever the coord's always return the same value as this is constant data.
  return 0;
}

DataTypes::RealVectorType::size_type
DataConstant::getPointOffset(int sampleNo,
                             int dataPointNo)
{
// We avoid this check for constant data due to issues manipulating 
// data with no samples.
 
//  ESYS_ASSERT((validSamplePointNo(dataPointNo) && validSampleNo(sampleNo)),
//              "Invalid index, sampleNo: " << sampleNo << " dataPointNo: " << dataPointNo);
  //
  // Whatever the coord's always return the same value as this is constant data.
  //

  return 0;
}

DataTypes::RealVectorType::size_type
DataConstant::getLength() const
{
    return std::max(m_data_c.size(), m_data_r.size());  
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
  
  if (isComplex()!=value->isComplex())
  {
    throw DataException("Error - cannot copy between slices of different complexity.");
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
  if (value->isComplex())
  {
      DataTypes::copySliceFrom(m_data_c,getShape(),0,tempDataConst->getVectorROC(), tempDataConst->getShape(),0,region_loop_range);
  }
  else
  {
      DataTypes::copySliceFrom(m_data_r,getShape(),0,tempDataConst->getVectorRO(), tempDataConst->getShape(),0,region_loop_range);
  }
}



void
DataConstant::symmetric(DataAbstract* ev)
{
  DataConstant* temp_ev=dynamic_cast<DataConstant*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataConstant::symmetric: casting to DataConstant failed (probably a programming error).");
  }
  if (isComplex())
  {
      escript::symmetric(m_data_c,getShape(),0,temp_ev->getVectorRWC(), temp_ev->getShape(),0);
  }
  else
  {
      escript::symmetric(m_data_r,getShape(),0,temp_ev->getVectorRW(), temp_ev->getShape(),0);
  }
}

void
DataConstant::antisymmetric(DataAbstract* ev)
{
  DataConstant* temp_ev=dynamic_cast<DataConstant*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataConstant::antisymmetric: casting to DataConstant failed (probably a programming error).");
  }
  if (isComplex())
  {
      escript::antisymmetric(m_data_c,getShape(),0,temp_ev->getVectorRWC(), temp_ev->getShape(),0);
  }
  else
  {
      escript::antisymmetric(m_data_r,getShape(),0,temp_ev->getVectorRW(), temp_ev->getShape(),0);
  }  
}

void
DataConstant::hermitian(DataAbstract* ev)
{
  DataConstant* temp_ev=dynamic_cast<DataConstant*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataConstant::hermitian: casting to DataConstant failed (probably a programming error).");
  }
  if (!isComplex() || !temp_ev->isComplex())
  {
      throw DataException("DataTagged::hermitian: do not call this method with real data");
  }  
  escript::hermitian(m_data_c,getShape(),0,temp_ev->getVectorRWC(), temp_ev->getShape(),0);
}

void
DataConstant::antihermitian(DataAbstract* ev)
{
  DataConstant* temp_ev=dynamic_cast<DataConstant*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataConstant::antihermitian: casting to DataConstant failed (probably a programming error).");
  }
  if (!isComplex() || !temp_ev->isComplex())
  {
      throw DataException("DataTagged::antihermitian: do not call this method with real data");
  }  
  escript::antihermitian(m_data_c,getShape(),0,temp_ev->getVectorRWC(), temp_ev->getShape(),0);
}


void
DataConstant::trace(DataAbstract* ev, int axis_offset)
{
  DataConstant* temp_ev=dynamic_cast<DataConstant*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataConstant::trace: casting to DataConstant failed (probably a programming error).");
  }
  const ShapeType& evShape=temp_ev->getShape();
  if (isComplex())
  {
      escript::trace(m_data_c,getShape(),0,temp_ev->getVectorRWC(),evShape,0,axis_offset);
  }
  else
  {
      escript::trace(m_data_r,getShape(),0,temp_ev->getVectorRW(),evShape,0,axis_offset);    
  }
}

void
DataConstant::swapaxes(DataAbstract* ev, int axis0, int axis1)
{
  DataConstant* temp_ev=dynamic_cast<DataConstant*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataConstant::swapaxes: casting to DataConstant failed (probably a programming error).");
  }
  if (isComplex())
  {
      escript::swapaxes(m_data_c,getShape(),0,temp_ev->getVectorRWC(), temp_ev->getShape(),0,axis0,axis1);
  }
  else
  {
      escript::swapaxes(m_data_r,getShape(),0,temp_ev->getVectorRW(), temp_ev->getShape(),0,axis0,axis1);
  }
}

void
DataConstant::transpose(DataAbstract* ev, int axis_offset)
{
  DataConstant* temp_ev=dynamic_cast<DataConstant*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataConstant::transpose: casting to DataConstant failed (probably a programming error).");
  }
  if (isComplex())
  {
      escript::transpose(m_data_c, getShape(),0, temp_ev->getVectorRWC(),temp_ev->getShape(),0,axis_offset);
  }
  else
  {
      escript::transpose(m_data_r, getShape(),0, temp_ev->getVectorRW(),temp_ev->getShape(),0,axis_offset);
  }
}

void
DataConstant::eigenvalues(DataAbstract* ev)
{
  DataConstant* temp_ev=dynamic_cast<DataConstant*>(ev);
  if (temp_ev==0) {
    throw DataException("Error - DataConstant::eigenvalues: casting to DataConstant failed (probably a programming error).");
  }
  if (isComplex())
  {
      escript::eigenvalues(m_data_c,getShape(),0,temp_ev->getVectorRWC(), temp_ev->getShape(),0);    
  }
  else
  {
      escript::eigenvalues(m_data_r,getShape(),0,temp_ev->getVectorRW(), temp_ev->getShape(),0);
  }
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
  escript::eigenvalues_and_eigenvectors(m_data_r, getShape(),0,temp_ev->getVectorRW(), temp_ev->getShape(),0,temp_V->getVectorRW(), temp_V->getShape(),0,tol);
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
  int res=escript::matrix_inverse(m_data_r, getShape(), 0, temp->getVectorRW(), temp->getShape(), 0, 1, h);
  return res;
}

void
DataConstant::setToZero()
{
    CHECK_FOR_EX_WRITE
    if (isComplex())
    {
        DataTypes::CplxVectorType::size_type n=m_data_c.size();
        for (int i=0; i<n ;++i) m_data_c[i]=0.;
    }
    else 
    {
        DataTypes::RealVectorType::size_type n=m_data_r.size();
        for (int i=0; i<n ;++i) m_data_r[i]=0.;
    }
}
void DataConstant::dump_hdf5(const std::string fileName) const
{
#ifdef ESYS_HAVE_HDF5
    int rank = getRank();
    int fs_type=  getFunctionSpace().getTypeCode();
    const DataTypes::ShapeType& shape = getShape();
    JMPI mpiInfo(getFunctionSpace().getDomain()->getMPI());
    if (isComplex())
    {
        throw DataException("Error - DataConstant::dump_hdf5: complex data are not supported. Split into real and imaginary part.");
    }

    const std::string newFileName(mpiInfo->appendRankToFileName(fileName));

#ifdef ESYS_MPI
   /* Serialize I/O */
   const int mpi_iam = mpiInfo->rank;
   const int mpi_num = mpiInfo->size;
   MPI_Status status;
   if (mpi_iam > 0)
       MPI_Recv(&ndims, 0, MPI_INT, mpi_iam-1, 81802, mpiInfo->comm, &status);
#endif
    try
    {
        H5::H5File h5_file(newFileName, H5F_ACC_TRUNC);
        // .... add meta data ............
        uint h5_shape[DataTypes::maxRank]; // dataset dimensions
        for (uint i = 0; i < rank; i++) {
            h5_shape[i]= shape[i];
        }
        hsize_t h5_shape_dims[1] = {rank};
        H5::DataSet h5_dmeta = h5_file.createDataSet("meta", H5::PredType::NATIVE_UINT, H5::DataSpace(1, h5_shape_dims ) );
        h5_dmeta.write( h5_shape, H5::PredType::NATIVE_UINT);
        // data type
        hsize_t h5_typeid_dims[1] = { 1 };
        uint h5_type_id[1] = { 0 };
        H5::Attribute h5_typeid_attr = h5_dmeta.createAttribute("type_id", H5::PredType::NATIVE_UINT, H5::DataSpace(1,h5_typeid_dims ) );
        h5_typeid_attr.write( H5::PredType::NATIVE_UINT , h5_type_id );

        uint h5_rank[1] = { rank };
        H5::Attribute h5_rank_attr = h5_dmeta.createAttribute("rank", H5::PredType::NATIVE_UINT, H5::DataSpace(1, h5_typeid_dims ) );
        h5_rank_attr.write( H5::PredType::NATIVE_UINT , h5_rank );

        uint dfs_type[1] = { fs_type };
        H5::Attribute h5_fs_type_attr = h5_dmeta.createAttribute("function_space_type", H5::PredType::NATIVE_UINT, H5::DataSpace(1,h5_typeid_dims) );
        h5_fs_type_attr.write( H5::PredType::NATIVE_UINT , dfs_type );
        // .... end meta data ............
        // ... add data ....
        hsize_t h5_data_length[1] = { m_data_r.size() };
        const double* d_ptr=&(m_data_r[0]);
        H5::DataSet h5_dataset_data = h5_file.createDataSet("data", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(1 , h5_data_length ) );
        h5_dataset_data.write(d_ptr, H5::PredType::NATIVE_DOUBLE);
        // .. we are done ....
        h5_file.close();
    }
    // catch failure caused by the H5File operations
    catch (H5::Exception& error)
    {
        error.printErrorStack();
        throw DataException("Error - DataConstant:: creating HDF5 file failed.");
    }
#ifdef ESYS_MPI
   if ( mpi_iam < mpi_num-1 ) MPI_Send(&ndims, 0, MPI_INT, mpi_iam+1, 81802, MPI_COMM_WORLD);
#endif
#else
    throw DataException("DataConstant::dump_hdf5: not configured with HDF5. Please contact your installation manager.");
#endif

}

#ifdef NETCDF4
void
DataConstant::dump(const std::string fileName) const
{
#ifdef ESYS_HAVE_NETCDF
    vector<NcDim> ncdims;
    int rank = getRank();
    int type=  getFunctionSpace().getTypeCode();

    //long dims[DataTypes::maxRank];
    const double* d_ptr=&(m_data_r[0]);
    DataTypes::ShapeType shape = getShape();
    JMPI mpiInfo(getFunctionSpace().getDomain()->getMPI());
#ifdef ESYS_MPI
    int ndims =0;    
    const int mpi_iam = mpiInfo->rank;
    const int mpi_num = mpiInfo->size;
    MPI_Status status;

    /* Serialize NetCDF I/O */
    if (mpi_iam > 0)
        MPI_Recv(&ndims, 0, MPI_INT, mpi_iam-1, 81802, mpiInfo->comm, &status);
#endif
    // Create the file.
    const std::string newFileName(mpiInfo->appendRankToFileName(fileName));
    NcFile dataFile;
    try
    {
        dataFile.open(newFileName.c_str(), NcFile::FileMode::replace,   NcFile::FileFormat::classic64);
    }
    catch (exceptions::NcException& e)
    {
        throw DataException("Error - DataConstant:: opening of netCDF file for output failed.");
    }
    int line=0;
    try
    {
        const NcInt ni;
        dataFile.putAtt("type_id", ni, 0);
        line++;
        dataFile.putAtt("rank", ni, rank);
        line++;
        dataFile.putAtt("function_space_type", ni, type);
    }
    catch (exceptions::NcException& e)
    {
        switch (line)
        {
        case 0: throw DataException("Error - DataConstant:: appending data type to netCDF file failed.");
        case 1: throw DataException("Error - DataConstant:: appending rank attribute to netCDF file failed.");
        case 2: throw DataException("Error - DataConstant:: appending function space attribute to netCDF file failed.");
        }
    }
    if (rank == 0) {
        try
        {
            ncdims.push_back(dataFile.addDim("l",  1));
        }
        catch (exceptions::NcException& e)
        {
            throw DataException("Error - DataConstant:: appending ncdimension 0 to netCDF file failed.");           
        }
//         dims[0]=1,
#ifdef ESYS_MPI
        ndims=1;
#endif        
    } else {
#ifdef ESYS_MPI        
        ndims=rank;
#endif        
//         dims[0]=shape[0];
        try
        {
            ncdims.push_back(dataFile.addDim("d0", shape[0]));
        }
        catch (exceptions::NcException& e)
        {
            throw DataException("Error - DataConstant:: appending ncdimension 0 to netCDF file failed.");
        }
        if ( rank >1 ) {
//             dims[1]=shape[1];
            try
            {
                ncdims.push_back(dataFile.addDim("d1", shape[1]));
            }
            catch (exceptions::NcException& e)
            {
                throw DataException("Error - DataConstant:: appending ncdimension 1 to netCDF file failed.");
            }
        }
        if ( rank >2 ) {
//             dims[2]=shape[2];
            try
            {
                ncdims.push_back(dataFile.addDim("d2",  shape[2]));
            } 
            catch (exceptions::NcException& e)
            {
                throw DataException("Error - DataConstant:: appending ncdimension 2 to netCDF file failed.");
            }
        }
        if ( rank >3 ) {
//             dims[3]=shape[3];
            try
            {
                ncdims.push_back(dataFile.addDim("d3",  shape[3]));
            } 
            catch (exceptions::NcException& e)
            {
                throw DataException("Error - DataConstant:: appending ncdimension 3 to netCDF file failed.");
            }
        }
    }

    try
    {
        line=0;     
        NcVar var = dataFile.addVar("data", ncDouble, ncdims);
        line++;
            // d_ptr is the data itself, dims is the shape
        // var.put(d_ptr,dims);
        var.putVar(d_ptr);
    }
    catch (exceptions::NcException& e)
    {
        if (line==0)
        {
            throw DataException("Error - DataConstant:: appending variable to netCDF file failed.");
        }
        throw DataException("Error - DataConstant:: copy data to netCDF buffer failed.");
    }
#ifdef ESYS_MPI
    if (mpi_iam<mpi_num-1) MPI_Send(&ndims, 0, MPI_INT, mpi_iam+1, 81802, MPI_COMM_WORLD);
#endif
    #else
    throw DataException("Error - DataConstant:: dump is not configured with netCDF. Please contact your installation manager.");
    #endif
}

#else

// old netcdf implementation
void
DataConstant::dump(const std::string fileName) const
{
#ifdef ESYS_HAVE_NETCDF
   const NcDim* ncdims[DataTypes::maxRank];
   NcVar* var;
   int rank = getRank();
   int type=  getFunctionSpace().getTypeCode();
   int ndims =0;
   long dims[DataTypes::maxRank];
   const double* d_ptr=&(m_data_r[0]);
   DataTypes::ShapeType shape = getShape();
   JMPI mpiInfo(getFunctionSpace().getDomain()->getMPI());
#ifdef ESYS_MPI
   const int mpi_iam = mpiInfo->rank;
   const int mpi_num = mpiInfo->size;
   MPI_Status status;

   /* Serialize NetCDF I/O */
   if (mpi_iam > 0)
       MPI_Recv(&ndims, 0, MPI_INT, mpi_iam-1, 81802, mpiInfo->comm, &status);
#endif

   // netCDF error handler
   NcError err(NcError::verbose_nonfatal);
   // Create the file.
   const std::string newFileName(mpiInfo->appendRankToFileName(fileName));
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

#endif

// These used to be marked as inline in DataConstant.
// But they are marked virtual in DataReady
DataTypes::RealVectorType&
DataConstant::getVectorRW()
{
    CHECK_FOR_EX_WRITE
    return m_data_r;
}

const DataTypes::RealVectorType&
DataConstant::getVectorRO() const
{
    return m_data_r;
}

DataTypes::CplxVectorType&
DataConstant::getVectorRWC()
{
    CHECK_FOR_EX_WRITE
    return m_data_c;
}

const DataTypes::CplxVectorType&
DataConstant::getVectorROC() const
{
    return m_data_c;
}

DataTypes::RealVectorType&
DataConstant::getTypedVectorRW(DataTypes::real_t dummy)
{
    CHECK_FOR_EX_WRITE
    return m_data_r;
}

const DataTypes::RealVectorType&
DataConstant::getTypedVectorRO(DataTypes::real_t dummy) const
{
    return m_data_r;
}

DataTypes::CplxVectorType&
DataConstant::getTypedVectorRW(DataTypes::cplx_t dummy)
{
    CHECK_FOR_EX_WRITE
    return m_data_c;
}

const DataTypes::CplxVectorType&
DataConstant::getTypedVectorRO(DataTypes::cplx_t dummy) const
{
    return m_data_c;
}

void DataConstant::complicate()
{
    if (!isComplex())
    {
        fillComplexFromReal(m_data_r, m_data_c);
        this->m_iscompl=true;
        m_data_r.resize(0,0,1);
    }
}

}  // end of namespace

