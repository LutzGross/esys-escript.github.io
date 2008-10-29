
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


#include "DataFactory.h"
#include "esysUtils/esys_malloc.h"

#include <boost/python/extract.hpp>
#include <iostream>
#include <exception>
#ifdef USE_NETCDF
#include <netcdfcpp.h>
#endif
#ifdef PASO_MPI
#include <mpi.h>
#endif

using namespace boost::python;

namespace escript {

Data
Scalar(double value,
       const FunctionSpace& what,
       bool expanded)
{
    //
    // an empty shape is a scalar
    DataTypes::ShapeType shape;
    return Data(value,shape,what,expanded);
}

Data
Vector(double value,
       const FunctionSpace& what,
       bool expanded)
{
    DataTypes::ShapeType shape(1,what.getDomain()->getDim());
    return Data(value,shape,what,expanded);
}

Data
Tensor(double value,
       const FunctionSpace& what,
       bool expanded)
{
    DataTypes::ShapeType shape(2,what.getDomain()->getDim());
    return Data(value,shape,what,expanded);
}

Data
Tensor3(double value,
        const FunctionSpace& what,
        bool expanded)
{
    DataTypes::ShapeType shape(3,what.getDomain()->getDim());
    return Data(value,shape,what,expanded);
}

Data
Tensor4(double value,
        const FunctionSpace& what,
        bool expanded)
{
    DataTypes::ShapeType shape(4,what.getDomain()->getDim());
    return Data(value,shape,what,expanded);
}

Data 
load(const std::string fileName,
     const AbstractDomain& domain)
{
   #ifdef USE_NETCDF
   NcAtt *type_att, *rank_att, *function_space_type_att;
   // netCDF error handler
   NcError err(NcError::silent_nonfatal);
   int mpi_iam=0, mpi_num=1;
   // Create the file.
#ifdef PASO_MPI
   MPI_Comm_rank(MPI_COMM_WORLD, &mpi_iam);
   MPI_Comm_size(MPI_COMM_WORLD, &mpi_num);
#endif
   char *newFileName = Escript_MPI_appendRankToFileName(fileName.c_str(), mpi_num, mpi_iam);
   NcFile dataFile(newFileName, NcFile::ReadOnly);
   if (!dataFile.is_valid())
        throw DataException("Error - load:: opening of netCDF file for input failed.");
   /* recover function space */
   if (! (function_space_type_att=dataFile.get_att("function_space_type")) )
        throw DataException("Error - load:: cannot recover function_space_type attribute from escript netCDF file.");
   int function_space_type = function_space_type_att->as_int(0);
   delete function_space_type_att;
   /* test if function space id is valid and create function space instance */
   if (! domain.isValidFunctionSpaceType(function_space_type) ) 
        throw DataException("Error - load:: function space type code in netCDF file is invalid for given domain.");
   FunctionSpace function_space=FunctionSpace(domain.getPtr(), function_space_type);
   /* recover rank */
   if (! (rank_att=dataFile.get_att("rank")) )
        throw DataException("Error - load:: cannot recover rank attribute from escript netCDF file.");
   int rank = rank_att->as_int(0);
   delete rank_att;
   if (rank<0 || rank>DataTypes::maxRank)
        throw DataException("Error - load:: rank in escript netCDF file is greater than maximum rank.");
   /* recover type attribute */
   int type=-1;
   if ((type_att=dataFile.get_att("type")) ) {
       char* type_str = type_att->as_string(0);
       if (strncmp(type_str, "constant", strlen("constant")) == 0 ) {
          type =0;
       } else if (strncmp(type_str, "tagged", strlen("tagged")) == 0 ) {
           type =1;
       } else if (strncmp(type_str, "expanded", strlen("expanded")) == 0 ) {
           type =2;
       }
       esysUtils::free(type_str);
   } else {
      if (! (type_att=dataFile.get_att("type_id")) )
  	throw DataException("Error - load:: cannot recover type attribute from escript netCDF file.");
      type=type_att->as_int(0);
   }
   delete type_att;

   /* recover dimension */
   int ndims=dataFile.num_dims();
   int ntags =0 , num_samples =0 , num_data_points_per_sample =0, d=0, len_data_point=1;
   NcDim *d_dim, *tags_dim, *num_samples_dim, *num_data_points_per_sample_dim;
   /* recover shape */
   DataTypes::ShapeType shape;
   long dims[DataTypes::maxRank+2];
   if (rank>0) {
     if (! (d_dim=dataFile.get_dim("d0")) )
          throw DataException("Error - load:: unable to recover d0 from netCDF file.");
      d=d_dim->size();
      shape.push_back(d);
      dims[0]=d;
      len_data_point*=d;
   }
   if (rank>1) {
     if (! (d_dim=dataFile.get_dim("d1")) )
          throw DataException("Error - load:: unable to recover d1 from netCDF file.");
      d=d_dim->size();
      shape.push_back(d);
      dims[1]=d;
      len_data_point*=d;
   }
   if (rank>2) {
     if (! (d_dim=dataFile.get_dim("d2")) )
          throw DataException("Error - load:: unable to recover d2 from netCDF file.");
      d=d_dim->size();
      shape.push_back(d);
      dims[2]=d;
      len_data_point*=d;
   }
   if (rank>3) {
     if (! (d_dim=dataFile.get_dim("d3")) )
          throw DataException("Error - load:: unable to recover d3 from netCDF file.");
      d=d_dim->size();
      shape.push_back(d);
      dims[3]=d;
      len_data_point*=d;
   }
   /* recover stuff */
   Data out;
   NcVar *var, *ids_var, *tags_var;
   if (type == 0) {
      /* constant data */
      if ( ! ( (ndims == rank && rank >0) || ( ndims ==1 && rank == 0 ) ) )
          throw DataException("Error - load:: illegal number of dimensions for constant data in netCDF file.");
      if (rank == 0) {
          if (! (d_dim=dataFile.get_dim("l")) )
              throw DataException("Error - load:: unable to recover d0 for scalar constant data in netCDF file.");
          int d0 = d_dim->size();
          if (! d0 == 1) 
              throw DataException("Error - load:: d0 is expected to be one for scalar constant data in netCDF file.");
          dims[0]=1;
      }
      out=Data(0,shape,function_space);
      if (!(var = dataFile.get_var("data")))
              throw DataException("Error - load:: unable to find data in netCDF file.");
      if (! var->get(&(out.getDataAtOffset(out.getDataOffset(0,0))), dims) ) 
              throw DataException("Error - load:: unable to recover data from netCDF file.");
   } else if (type == 1) { 
      /* tagged data */
      if ( ! (ndims == rank + 1) )
         throw DataException("Error - load:: illegal number of dimensions for tagged data in netCDF file.");
      if (! (tags_dim=dataFile.get_dim("num_tags")) )
         throw DataException("Error - load:: unable to recover number of tags from netCDF file.");
      ntags=tags_dim->size();
      dims[rank]=ntags;
      int *tags = (int *) esysUtils::malloc(ntags*sizeof(int));
      if (! ( tags_var = dataFile.get_var("tags")) )
      {
         esysUtils::free(tags);
         throw DataException("Error - load:: unable to find tags in netCDF file.");
      }
      if (! tags_var->get(tags, ntags) ) 
      {
         esysUtils::free(tags);
         throw DataException("Error - load:: unable to recover tags from netCDF file.");
      }

// Current Version
/*      DataVector data(len_data_point * ntags, 0., len_data_point * ntags);
      if (!(var = dataFile.get_var("data")))
      {
         esysUtils::free(tags);
         throw DataException("Error - load:: unable to find data in netCDF file.");
      }
      if (! var->get(&(data[0]), dims) ) 
      {
         esysUtils::free(tags);
         throw DataException("Error - load:: unable to recover data from netCDF file.");
      }
      out=Data(DataArrayView(data,shape,0),function_space);
      for (int t=1; t<ntags; ++t) {
	 out.setTaggedValueFromCPP(tags[t],shape, data, t*len_data_point);
//         out.setTaggedValueFromCPP(tags[t],DataArrayView(data,shape,t*len_data_point));
      }*/
// End current version
	
// New version

	// A) create a DataTagged dt
	// B) Read data from file
	// C) copy default value into dt
	// D) copy tagged values into dt
	// E) create a new Data based on dt

      NcVar* var1;
      DataVector data1(len_data_point * ntags, 0., len_data_point * ntags);
      if (!(var1 = dataFile.get_var("data")))
      {
         esysUtils::free(tags);
         throw DataException("Error - load:: unable to find data in netCDF file.");
      }
      if (! var1->get(&(data1[0]), dims) ) 
      {
         esysUtils::free(tags);
         throw DataException("Error - load:: unable to recover data from netCDF file.");
      }
      DataTagged* dt=new DataTagged(function_space, shape, tags,data1);
      out=Data(dt);
      esysUtils::free(tags);
   } else if (type == 2) {
      /* expanded data */
      if ( ! (ndims == rank + 2) )
          throw DataException("Error - load:: illegal number of dimensions for exanded data in netCDF file.");
      if ( ! (num_samples_dim = dataFile.get_dim("num_samples") ) )
          throw DataException("Error - load:: unable to recover number of samples from netCDF file.");
      num_samples = num_samples_dim->size();
      if ( ! (num_data_points_per_sample_dim = dataFile.get_dim("num_data_points_per_sample") ) )
          throw DataException("Error - load:: unable to recover number of data points per sample from netCDF file.");
      num_data_points_per_sample=num_data_points_per_sample_dim->size();
      // check shape:
      if ( ! (num_samples == function_space.getNumSamples() && num_data_points_per_sample == function_space.getNumDataPointsPerSample()) )
          throw DataException("Error - load:: data sample layout of file does not match data layout of function space.");
      if (num_samples==0) {
	out = Data(0,shape,function_space,true);
      }
      else {
	// get ids
	if (! ( ids_var = dataFile.get_var("id")) )
		throw DataException("Error - load:: unable to find reference ids in netCDF file.");
	const int* ids_p=function_space.borrowSampleReferenceIDs();
	int *ids_of_nc = (int *)esysUtils::malloc(num_samples*sizeof(int));
	if (! ids_var->get(ids_of_nc, (long) num_samples) ) 
	{
		esysUtils::free(ids_of_nc);
		throw DataException("Error - load:: unable to recover ids from netCDF file.");
	}
	// check order:
	int failed=-1, local_failed=-1, i;
	#pragma omp parallel private(local_failed)
	{
		local_failed=-1;
		#pragma omp for private(i) schedule(static)
		for (i=0;i < num_samples; ++i) {
		if (ids_of_nc[i]!=ids_p[i]) local_failed=i;
		}
		#pragma omp critical
		if (local_failed>=0) failed = local_failed;
	}
	/* if (failed>=0) 
	{
		esysUtils::free(ids_of_nc);
		throw DataException("Error - load:: data ordering in netCDF file does not match ordering of FunctionSpace.");
	} */
	// get the data:
	dims[rank]=num_data_points_per_sample;
	dims[rank+1]=num_samples;
	out=Data(0,shape,function_space,true);
	if (!(var = dataFile.get_var("data")))
	{
		esysUtils::free(ids_of_nc);
		throw DataException("Error - load:: unable to find data in netCDF file.");
	}
	if (! var->get(&(out.getDataAtOffset(out.getDataOffset(0,0))), dims) ) 
	{
		esysUtils::free(ids_of_nc);
		throw DataException("Error - load:: unable to recover data from netCDF file.");
	}
	if (failed>=0) {
		try {
		std::cout << "Information - load: start reordering data from netCDF file " << fileName << std::endl;
		out.borrowData()->reorderByReferenceIDs(ids_of_nc);
		} 
		catch (std::exception&) {
		esysUtils::free(ids_of_nc);
		throw DataException("Error - load:: unable to reorder data in netCDF file.");
		}
	}
      }
   } else {
       throw DataException("Error - load:: unknown escript data type in netCDF file.");
   }
   return out;
   #else
   throw DataException("Error - load:: is not compiled with netCFD. Please contact your insatllation manager.");
   #endif
}

bool 
loadConfigured()
{
   #ifdef USE_NETCDF
   return true;
   #else
   return false;
   #endif
}

Data
convertToData(const boost::python::object& value,
              const FunctionSpace& what) 
{
     // first we try to extract a Data object from value 
     extract<Data> value_data(value);
     if (value_data.check()) {
         Data extracted_data=value_data();
         if (extracted_data.isEmpty()) {
            return extracted_data;
         } else {
            return Data(extracted_data,what);
         }
     } else {
        return Data(value,what);
     }
}

}  // end of namespace
