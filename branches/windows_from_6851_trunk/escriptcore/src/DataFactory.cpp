
/*****************************************************************************
*
* Copyright (c) 2003-2018 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

#include "DataFactory.h"

#include <boost/python/extract.hpp>
#include <boost/scoped_array.hpp>

#include <exception>
#include <iostream>
#include <fstream>


#ifdef ESYS_HAVE_NETCDF
 #ifdef NETCDF4
  #include <ncDim.h>
  #include <ncVar.h>
  #include <ncFile.h>
  
 #include "NCHelper.h"  
  
 #else
  #include <netcdfcpp.h>
 #endif
#endif

namespace bp = boost::python;
#ifdef NETCDF4
using namespace netCDF;
#endif


namespace escript {

Data Scalar(double value, const FunctionSpace& what, bool expanded)
{
    // an empty shape is a scalar
    DataTypes::ShapeType shape;
    return Data(value, shape, what, expanded);
}

Data Scalar(DataTypes::cplx_t value, const FunctionSpace& what, bool expanded)
{
    // an empty shape is a scalar
    DataTypes::ShapeType shape;
    return Data(value, shape, what, expanded);
}

Data
ScalarFromObj(boost::python::object o,
	const FunctionSpace& what,
	bool expanded)
{
    // check for real first
    try {
        double v = bp::extract<double>(o);
        return Scalar(v, what, expanded);
    } catch(...) {
        PyErr_Clear();
    }    
    // check for real first
    try {
        DataTypes::cplx_t v = bp::extract<DataTypes::cplx_t>(o);
        return Scalar(v, what, expanded);
    } catch(...) {
        PyErr_Clear();
    }    
    throw DataException("Can not make a Scalar from a non-scalar value.");
}

Data Vector(double value, const FunctionSpace& what, bool expanded)
{
    DataTypes::ShapeType shape(1, what.getDomain()->getDim());
    return Data(value, shape, what, expanded);
}

Data VectorFromObj(bp::object o, const FunctionSpace& what, bool expanded)
{    
    // first try to get a double and route it to the other method
    try {
        double v = bp::extract<double>(o);
        return Vector(v, what, expanded);
    } catch(...) {
        PyErr_Clear();
    }
    DataTypes::ShapeType shape(1, what.getDomain()->getDim());
    Data d(o, what, expanded);
    if (d.getDataPointShape() != shape) {
        throw DataException("VectorFromObj: Shape of vector passed to function"
               " does not match the dimension of the domain. ");
    }
    return d;
}

Data Tensor(double value, const FunctionSpace& what, bool expanded)
{
    DataTypes::ShapeType shape(2, what.getDomain()->getDim());
    return Data(value, shape, what, expanded);
}

Data TensorC(DataTypes::cplx_t value, const FunctionSpace& what, bool expanded)
{
    DataTypes::ShapeType shape(2, what.getDomain()->getDim());
    return Data(value, shape, what, expanded);
}

// We need to take some care here because this signature trumps the other one from boost's point of view
Data TensorFromObj(bp::object o, const FunctionSpace& what, bool expanded)
{
    // first try to get a double and route it to the other method
    try {
        double v = bp::extract<double>(o);
        return Tensor(v, what, expanded);
    } catch(...) {
        PyErr_Clear();
    }
    // now try to get a complex and route to scalar factory
    try {
        DataTypes::cplx_t v = bp::extract<DataTypes::cplx_t>(o);
        return TensorC(v, what, expanded);
    } catch(...) {
        PyErr_Clear();
    }    
    DataTypes::ShapeType shape(2, what.getDomain()->getDim());
    Data d(o, what, expanded);
    if (d.getDataPointShape() != shape) {
        throw DataException("TensorFromObj: Shape of tensor passed to function"
               " does not match the dimension of the domain.");
    }
    return d;
}

Data Tensor3(double value, const FunctionSpace& what, bool expanded)
{
    DataTypes::ShapeType shape(3, what.getDomain()->getDim());
    return Data(value, shape, what, expanded);
}

Data Tensor3C(DataTypes::cplx_t  value, const FunctionSpace& what, bool expanded)
{
    DataTypes::ShapeType shape(3, what.getDomain()->getDim());
    return Data(value, shape, what, expanded);
}


Data Tensor3FromObj(bp::object o, const FunctionSpace& what, bool expanded)
{
    // first try to get a double and route it to the other method
    try {
        double v = bp::extract<double>(o);
        return Tensor3(v, what, expanded);
    } catch(...) {
        PyErr_Clear();
    }
    // first try to get a complex and route it to the other method
    try {
        DataTypes::cplx_t v = bp::extract<DataTypes::cplx_t>(o);
        return Tensor3C(v, what, expanded);
    } catch(...) {
        PyErr_Clear();
    }    
    DataTypes::ShapeType shape(3, what.getDomain()->getDim());
    Data d(o, what, expanded);
    if (d.getDataPointShape() != shape) {
        throw DataException("Tensor3FromObj: Shape of tensor passed to "
                "function does not match the dimension of the domain.");
    }
    return d;
}

Data Tensor4(double value, const FunctionSpace& what, bool expanded)
{
    DataTypes::ShapeType shape(4, what.getDomain()->getDim());
    return Data(value, shape, what, expanded);
}

Data Tensor4C(DataTypes::cplx_t value, const FunctionSpace& what, bool expanded)
{
    DataTypes::ShapeType shape(4, what.getDomain()->getDim());
    return Data(value, shape, what, expanded);
}


Data Tensor4FromObj(bp::object o, const FunctionSpace& what, bool expanded)
{
    // first try to get a double and route it to the other method
    try {
        double v = bp::extract<double>(o);
        return Tensor4(v, what, expanded);
    } catch(...) {
        PyErr_Clear();
    }
    // first try to get a double and route it to the other method
    try {
        DataTypes::cplx_t v = bp::extract<DataTypes::cplx_t>(o);
        return Tensor4C(v, what, expanded);
    } catch(...) {
        PyErr_Clear();
    }    
    DataTypes::ShapeType shape(4, what.getDomain()->getDim());
    Data d(o, what, expanded);
    if (d.getDataPointShape() != shape) {
        throw DataException("VectorFromObj: Shape of tensor passed to function"
               " does not match the dimension of the domain.");
    }
    return d;
}

#ifdef NETCDF4

// Warning: at time of writing, calls to ncVar::getVar are not range checked as was done under the previous API
// we want to see if it works first, but at some stage that should be tightened up

Data load(const std::string fileName, const AbstractDomain& domain)
{
#ifdef ESYS_HAVE_NETCDF
    JMPI mpiInfo(domain.getMPI());
    const std::string newFileName(mpiInfo->appendRankToFileName(fileName));
    NcFile dataFile;
    if (!openNcFile(dataFile, newFileName))
    {
        throw DataException("load: opening of netCDF file for input failed.");
    }
    Data out;
    int error = 0;
    std::string msg;
    try {
        int line=0;
        int rank=-1;
        int type=-1;
        int function_space_type=0;  // =0 should not actually be used but we keep compiler happy
        try
        {
            // recover function space            
            NcGroupAtt fst=dataFile.getAtt("function_space_type");
            if (fst.getAttLength()!=1)
            {
                throw DataException("load: oversize attribute function_space_type");
            }
            fst.getValues(&function_space_type);
        }
        catch (exceptions::NcException* e)
        {
                throw DataException("load: cannot recover function_space_type attribute from escript netCDF file.");    
        }
        line=0;
        if (!domain.isValidFunctionSpaceType(function_space_type))
        {       // jump to outer catch
            throw DataException("load: function space type code in netCDF file is invalid for given domain.");
        }
        FunctionSpace function_space=FunctionSpace(domain.getPtr(), function_space_type);
        try
        {
            line++;
            NcGroupAtt rt=dataFile.getAtt("rank");
            if (rt.getAttLength()!=1)
            {
                throw DataException("load: oversize attribute rank");
            }
            rt.getValues(&rank);
            line++;
            if (rank<0 || rank>DataTypes::maxRank)
                throw DataException("load: rank in escript netCDF file is greater than maximum rank.");            
            // recover type attribute
            //   looks like we can have either "type" or "typeid"
            NcGroupAtt tatt=dataFile.getAtt("type");
            if (!tatt.isNull())
            {
                std::string type_str;
                tatt.getValues(&type_str);
                if (type_str=="constant")
                {
                    type = 0;
                }
                else if (type_str=="tagged")
                {
                    type = 1;
                }
                else // if (type_str=="expanded")
                {
                    type = 2;
                }
            }
            else 
            {
                tatt=dataFile.getAtt("type_id");
                if (!tatt.isNull())
                {
                    if (tatt.getAttLength()>1)
                    {
                        throw DataException("load: oversize attribute type_id");
                    }
                    //type=dataFile.getAtt("type_id").as_int(0);
                    tatt.getValues(&type);
                }
                else
                {
                    throw DataException("load: cannot recover type attribute from escript netCDF file.");
                }
            }
        }
        catch (exceptions::NcException* e)
        {
            switch (line)
            {
            case 1: throw DataException("load: cannot recover rank attribute from escript netCDF file.");
            default:
                throw DataException("load: unspecified error.");
            }
        }
        // recover dimension
        int ndims=dataFile.getDimCount();
        int ntags =0 , num_samples =0 , num_data_points_per_sample =0, d=0, len_data_point=1;
        NcDim d_dim, tags_dim, num_samples_dim, num_data_points_per_sample_dim;
        /* recover shape */
        DataTypes::ShapeType shape;
//        long dims[DataTypes::maxRank+2];
        if (rank>0) {
            if ((d_dim=dataFile.getDim("d0")).isNull() )
                throw DataException("load: unable to recover d0 from netCDF file.");
            d=d_dim.getSize();
            shape.push_back(d);
//             dims[0]=d;
            len_data_point*=d;
        }
        if (rank>1) {
            if ((d_dim=dataFile.getDim("d1")).isNull() )
                throw DataException("load: unable to recover d1 from netCDF file.");
            d=d_dim.getSize();
            shape.push_back(d);
//             dims[1]=d;
            len_data_point*=d;
        }
        if (rank>2) {
            if ((d_dim=dataFile.getDim("d2")).isNull() )
                throw DataException("load: unable to recover d2 from netCDF file.");
            d=d_dim.getSize();
            shape.push_back(d);
//             dims[2]=d;
            len_data_point*=d;
        }
        if (rank>3) {
            if ((d_dim=dataFile.getDim("d3")).isNull() )
                throw DataException("load: unable to recover d3 from netCDF file.");
            d=d_dim.getSize();
            shape.push_back(d);
//             dims[3]=d;
            len_data_point*=d;
        }

        NcVar var, ids_var, tags_var;
        if (type == 0) {
            // constant data
            if ( ! ( (ndims == rank && rank >0) || ( ndims ==1 && rank == 0 ) ) )
                throw DataException("load: illegal number of dimensions for constant data in netCDF file.");
            if (rank == 0) {
                if ((d_dim=dataFile.getDim("l")).isNull() )
                    throw DataException("load: unable to recover d0 for scalar constant data in netCDF file.");
                int d0 = d_dim.getSize();
                if (d0 != 1)
                    throw DataException("load: d0 is expected to be one for scalar constant data in netCDF file.");
//                 dims[0]=1;
            }
            out=Data(0,shape,function_space,false);
            if ((var = dataFile.getVar("data")).isNull())
                throw DataException("load: unable to find data in netCDF file.");
            var.getVar(&(out.getDataAtOffsetRW(out.getDataOffset(0,0), static_cast<DataTypes::real_t>(0))));
        } else if (type == 1) { 
            // tagged data
            if ( ! (ndims == rank + 1) )
                throw DataException("load: illegal number of dimensions for tagged data in netCDF file.");
            if ((tags_dim=dataFile.getDim("num_tags")).isNull() )
                throw DataException("load: unable to recover number of tags from netCDF file.");
            ntags=tags_dim.getSize();
//             dims[rank]=ntags;
            std::vector<int> tags(ntags);
            if (( tags_var = dataFile.getVar("tags")).isNull() )
                throw DataException("load: unable to find tags in netCDF file.");
            // oversize could be a problem here?
            tags_var.getVar(&tags[0]);

            // A) create a DataTagged dt
            // B) Read data from file
            // C) copy default value into dt
            // D) copy tagged values into dt
            // E) create a new Data based on dt

            NcVar var1;
            DataTypes::RealVectorType data1(len_data_point * ntags, 0., len_data_point * ntags);
            if ((var1 = dataFile.getVar("data")).isNull())
                throw DataException("load: unable to find data in netCDF file.");
            var1.getVar(&(data1[0])); 
            DataTagged* dt=new DataTagged(function_space, shape, &tags[0], data1);
            out=Data(dt);
        } else if (type == 2) {
            // expanded data
            if ( ! (ndims == rank + 2) )
                throw DataException("load: illegal number of dimensions for expanded data in netCDF file.");
            if ((num_samples_dim = dataFile.getDim("num_samples")).isNull() )
                throw DataException("load: unable to recover number of samples from netCDF file.");
            num_samples = num_samples_dim.getSize();
            if ((num_data_points_per_sample_dim = dataFile.getDim("num_data_points_per_sample")).isNull() )
                throw DataException("load: unable to recover number of data points per sample from netCDF file.");
            num_data_points_per_sample=num_data_points_per_sample_dim.getSize();
            // check shape:
            if ( ! (num_samples == function_space.getNumSamples() && num_data_points_per_sample == function_space.getNumDataPointsPerSample()) )
                throw DataException("load: data sample layout of file does not match data layout of function space.");
            if (num_samples==0) {
                out = Data(0,shape,function_space,true);
            } else {
                // get ids
                if (( ids_var = dataFile.getVar("id")).isNull() )
                    throw DataException("load: unable to find reference ids in netCDF file.");
                const DataTypes::dim_t* ids_p=function_space.borrowSampleReferenceIDs();
                std::vector<DataTypes::dim_t> ids_of_nc(num_samples);
                // oversize could be a problem here
                ids_var.getVar(&ids_of_nc[0]);
                // check order:
                int failed=-1, local_failed=-1, i;
#pragma omp parallel private(local_failed)
                {
                    local_failed=-1;
#pragma omp for private(i) schedule(static)
                    for (i=0; i < num_samples; ++i) {
                        if (ids_of_nc[i]!=ids_p[i]) local_failed=i;
                    }
#pragma omp critical
                    if (local_failed>=0) failed = local_failed;
                }
                // get the data:
//                 dims[rank]=num_data_points_per_sample;
//                 dims[rank+1]=num_samples;
                out=Data(0,shape,function_space,true);
                if ((var = dataFile.getVar("data")).isNull())
                    throw DataException("load: unable to find data in netCDF file.");
                var.getVar(&(out.getDataAtOffsetRW(out.getDataOffset(0,0), static_cast<DataTypes::real_t>(0))));
                if (failed >= 0) {
                    try {
                        std::cout << "Information - load: start reordering data from netCDF file " << fileName << std::endl;
                        out.borrowData()->reorderByReferenceIDs(&ids_of_nc[0]);
                    } catch (std::exception&) {
                        throw DataException("load: unable to reorder data in netCDF file.");
                    }
                }
            }
        } else {
            throw DataException("load: unknown escript data type in netCDF file.");
        }
    } catch (DataException& e) {
        error=1;
        msg=e.what();
    }
    int gerror = error;
    checkResult(error, gerror, mpiInfo);
    if (gerror > 0) {
        char* gmsg;
        shipString(msg.c_str(), &gmsg, mpiInfo->comm);
        throw DataException(gmsg);
    }
    return out;
#else
    throw DataException("load: not compiled with netCDF. Please contact your"
                        " installation manager.");
#endif // ESYS_HAVE_NETCDF
}

#else

Data load(const std::string fileName, const AbstractDomain& domain)
{
#ifdef ESYS_HAVE_NETCDF
    NcAtt *type_att, *rank_att, *function_space_type_att;
    // netCDF error handler
    NcError err(NcError::silent_nonfatal);
    JMPI mpiInfo(domain.getMPI());
    const std::string newFileName(mpiInfo->appendRankToFileName(fileName));
    NcFile dataFile(newFileName.c_str(), NcFile::ReadOnly);
    Data out;
    int error = 0;
    std::string msg;
    try {
        if (!dataFile.is_valid())
            throw DataException("load: opening of netCDF file for input failed.");
       // recover function space
        if (! (function_space_type_att=dataFile.get_att("function_space_type")) )
            throw DataException("load: cannot recover function_space_type attribute from escript netCDF file.");
        int function_space_type = function_space_type_att->as_int(0);
        delete function_space_type_att;
        // test if function space id is valid and create function space instance
        if (!domain.isValidFunctionSpaceType(function_space_type)) 
            throw DataException("load: function space type code in netCDF file is invalid for given domain.");
        FunctionSpace function_space=FunctionSpace(domain.getPtr(), function_space_type);
        // recover rank
        if (! (rank_att=dataFile.get_att("rank")) )
            throw DataException("load: cannot recover rank attribute from escript netCDF file.");
        int rank = rank_att->as_int(0);
        delete rank_att;
        if (rank<0 || rank>DataTypes::maxRank)
            throw DataException("load: rank in escript netCDF file is greater than maximum rank.");
        // recover type attribute
        int type=-1;
        if ((type_att=dataFile.get_att("type")) ) {
            boost::scoped_array<char> type_str(type_att->as_string(0));
            if (strncmp(type_str.get(), "constant", strlen("constant")) == 0 ) {
                type = 0;
            } else if (strncmp(type_str.get(), "tagged", strlen("tagged")) == 0 ) {
                type = 1;
            } else if (strncmp(type_str.get(), "expanded", strlen("expanded")) == 0 ) {
                type = 2;
            }
        } else {
            if (! (type_att=dataFile.get_att("type_id")) )
                throw DataException("load: cannot recover type attribute from escript netCDF file.");
            type=type_att->as_int(0);
        }
        delete type_att;

        // recover dimension
        int ndims=dataFile.num_dims();
        int ntags =0 , num_samples =0 , num_data_points_per_sample =0, d=0, len_data_point=1;
        NcDim *d_dim, *tags_dim, *num_samples_dim, *num_data_points_per_sample_dim;
        /* recover shape */
        DataTypes::ShapeType shape;
        long dims[DataTypes::maxRank+2];
        if (rank>0) {
            if (! (d_dim=dataFile.get_dim("d0")) )
                throw DataException("load: unable to recover d0 from netCDF file.");
            d=d_dim->size();
            shape.push_back(d);
            dims[0]=d;
            len_data_point*=d;
        }
        if (rank>1) {
            if (! (d_dim=dataFile.get_dim("d1")) )
                throw DataException("load: unable to recover d1 from netCDF file.");
            d=d_dim->size();
            shape.push_back(d);
            dims[1]=d;
            len_data_point*=d;
        }
        if (rank>2) {
            if (! (d_dim=dataFile.get_dim("d2")) )
                throw DataException("load: unable to recover d2 from netCDF file.");
            d=d_dim->size();
            shape.push_back(d);
            dims[2]=d;
            len_data_point*=d;
        }
        if (rank>3) {
            if (! (d_dim=dataFile.get_dim("d3")) )
                throw DataException("load: unable to recover d3 from netCDF file.");
            d=d_dim->size();
            shape.push_back(d);
            dims[3]=d;
            len_data_point*=d;
        }

        NcVar *var, *ids_var, *tags_var;
        if (type == 0) {
            // constant data
            if ( ! ( (ndims == rank && rank >0) || ( ndims ==1 && rank == 0 ) ) )
                throw DataException("load: illegal number of dimensions for constant data in netCDF file.");
            if (rank == 0) {
                if (! (d_dim=dataFile.get_dim("l")) )
                    throw DataException("load: unable to recover d0 for scalar constant data in netCDF file.");
                int d0 = d_dim->size();
                if (d0 != 1)
                    throw DataException("load: d0 is expected to be one for scalar constant data in netCDF file.");
                dims[0]=1;
            }
            out=Data(0,shape,function_space,false);
            if (!(var = dataFile.get_var("data")))
                throw DataException("load: unable to find data in netCDF file.");
            if (! var->get(&(out.getDataAtOffsetRW(out.getDataOffset(0,0), static_cast<DataTypes::real_t>(0))), dims) ) 
                throw DataException("load: unable to recover data from netCDF file.");
        } else if (type == 1) { 
            // tagged data
            if ( ! (ndims == rank + 1) )
                throw DataException("load: illegal number of dimensions for tagged data in netCDF file.");
            if (! (tags_dim=dataFile.get_dim("num_tags")) )
                throw DataException("load: unable to recover number of tags from netCDF file.");
            ntags=tags_dim->size();
            dims[rank]=ntags;
            std::vector<int> tags(ntags);
            if (! ( tags_var = dataFile.get_var("tags")) )
                throw DataException("load: unable to find tags in netCDF file.");
            if (! tags_var->get(&tags[0], ntags) ) 
                throw DataException("load: unable to recover tags from netCDF file.");

            // A) create a DataTagged dt
            // B) Read data from file
            // C) copy default value into dt
            // D) copy tagged values into dt
            // E) create a new Data based on dt

            NcVar* var1;
            DataTypes::RealVectorType data1(len_data_point * ntags, 0., len_data_point * ntags);
            if (!(var1 = dataFile.get_var("data")))
                throw DataException("load: unable to find data in netCDF file.");
            if (! var1->get(&(data1[0]), dims) ) 
                throw DataException("load: unable to recover data from netCDF file.");
            DataTagged* dt=new DataTagged(function_space, shape, &tags[0], data1);
            out=Data(dt);
        } else if (type == 2) {
            // expanded data
            if ( ! (ndims == rank + 2) )
                throw DataException("load: illegal number of dimensions for expanded data in netCDF file.");
            if ( ! (num_samples_dim = dataFile.get_dim("num_samples") ) )
                throw DataException("load: unable to recover number of samples from netCDF file.");
            num_samples = num_samples_dim->size();
            if ( ! (num_data_points_per_sample_dim = dataFile.get_dim("num_data_points_per_sample") ) )
                throw DataException("load: unable to recover number of data points per sample from netCDF file.");
            num_data_points_per_sample=num_data_points_per_sample_dim->size();
            // check shape:
            if ( ! (num_samples == function_space.getNumSamples() && num_data_points_per_sample == function_space.getNumDataPointsPerSample()) )
                throw DataException("load: data sample layout of file does not match data layout of function space.");
            if (num_samples==0) {
                out = Data(0,shape,function_space,true);
            } else {
                // get ids
                if (! ( ids_var = dataFile.get_var("id")) )
                    throw DataException("load: unable to find reference ids in netCDF file.");
                const DataTypes::dim_t* ids_p=function_space.borrowSampleReferenceIDs();
                std::vector<DataTypes::dim_t> ids_of_nc(num_samples);
                if (! ids_var->get(&ids_of_nc[0], (long) num_samples) ) 
                    throw DataException("load: unable to recover ids from netCDF file.");
                // check order:
                int failed=-1, local_failed=-1, i;
#pragma omp parallel private(local_failed)
                {
                    local_failed=-1;
#pragma omp for private(i) schedule(static)
                    for (i=0; i < num_samples; ++i) {
                        if (ids_of_nc[i]!=ids_p[i]) local_failed=i;
                    }
#pragma omp critical
                    if (local_failed>=0) failed = local_failed;
                }
                // get the data:
                dims[rank]=num_data_points_per_sample;
                dims[rank+1]=num_samples;
                out=Data(0,shape,function_space,true);
                if (!(var = dataFile.get_var("data")))
                    throw DataException("load: unable to find data in netCDF file.");
                if (! var->get(&(out.getDataAtOffsetRW(out.getDataOffset(0,0), static_cast<DataTypes::real_t>(0))), dims)) 
                    throw DataException("load: unable to recover data from netCDF file.");
                if (failed >= 0) {
                    try {
                        std::cout << "Information - load: start reordering data from netCDF file " << fileName << std::endl;
                        out.borrowData()->reorderByReferenceIDs(&ids_of_nc[0]);
                    } catch (std::exception&) {
                        throw DataException("load: unable to reorder data in netCDF file.");
                    }
                }
            }
        } else {
            throw DataException("load: unknown escript data type in netCDF file.");
        }
    } catch (DataException& e) {
        error=1;
        msg=e.what();
    }
    int gerror = error;
    checkResult(error, gerror, mpiInfo);
    if (gerror > 0) {
        char* gmsg;
        shipString(msg.c_str(), &gmsg, mpiInfo->comm);
        throw DataException(gmsg);
    }
    return out;
#else
    throw DataException("load: not compiled with netCDF. Please contact your"
                        " installation manager.");
#endif // ESYS_HAVE_NETCDF
}

#endif

bool loadConfigured()
{
#ifdef ESYS_HAVE_NETCDF
    return true;
#else
    return false;
#endif
}

Data convertToData(const bp::object& value, const FunctionSpace& what) 
{
    // first we try to extract a Data object from value 
    bp::extract<Data> value_data(value);
    if (value_data.check()) {
        Data extracted_data=value_data();
        if (extracted_data.isEmpty()) {
            return extracted_data;
        } else {
            return Data(extracted_data,what);
        }
    } else {
        return Data(value,what,false);
    }
}

}  // end of namespace

