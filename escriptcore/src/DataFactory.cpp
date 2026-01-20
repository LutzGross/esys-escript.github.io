
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#include "DataFactory.h"

#include <boost/python/extract.hpp>
#include <boost/scoped_array.hpp>

#include <exception>
#include <iostream>
#include <fstream>

#ifdef ESYS_HAVE_HDF5
  #include <H5Cpp.h>
#endif

namespace bp = boost::python;

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

Data ComplexScalar(double value, const FunctionSpace& what, bool expanded)
{
    // an empty shape is a scalar
    DataTypes::ShapeType shape;
    escript::Data newdata = Data(value, shape, what, expanded);
    newdata.complicate();
    return newdata;
}

Data ComplexScalar(DataTypes::cplx_t value, const FunctionSpace& what, bool expanded)
{
    // an empty shape is a scalar
    DataTypes::ShapeType shape;
    escript::Data newdata = Data(value, shape, what, expanded);
    newdata.complicate();
    return newdata;
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

Data
ComplexScalarFromObj(boost::python::object o,
	const FunctionSpace& what,
	bool expanded)
{
    // check for real first
    try {
        double v = bp::extract<double>(o);
        return ComplexScalar(v, what, expanded);
    } catch(...) {
        PyErr_Clear();
    }
    // check for real first
    try {
        DataTypes::cplx_t v = bp::extract<DataTypes::cplx_t>(o);
        return ComplexScalar(v, what, expanded);
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

Data ComplexVector(double value, const FunctionSpace& what, bool expanded)
{
    DataTypes::ShapeType shape(1, what.getDomain()->getDim());
    escript::Data newdata = Data(value, shape, what, expanded);
    newdata.complicate();
    return newdata;
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

Data ComplexVectorFromObj(bp::object o, const FunctionSpace& what, bool expanded)
{
    // first try to get a double and route it to the other method
    try {
        double v = bp::extract<double>(o);
        return ComplexVector(v, what, expanded);
    } catch(...) {
        PyErr_Clear();
    }
    DataTypes::ShapeType shape(1, what.getDomain()->getDim());
    Data d(o, what, expanded);
    d.complicate();
    if (d.getDataPointShape() != shape) {
        throw DataException("ComplexVectorFromObj: Shape of vector passed to function"
               " does not match the dimension of the domain. ");
    }
    return d;
}

Data Tensor(double value, const FunctionSpace& what, bool expanded)
{
    DataTypes::ShapeType shape(2, what.getDomain()->getDim());
    return Data(value, shape, what, expanded);
}

Data ComplexTensor(double value, const FunctionSpace& what, bool expanded)
{
    DataTypes::ShapeType shape(2, what.getDomain()->getDim());
    escript::Data newdata = Data(value, shape, what, expanded);
    newdata.complicate();
    return newdata;
}

Data TensorC(DataTypes::cplx_t value, const FunctionSpace& what, bool expanded)
{
    DataTypes::ShapeType shape(2, what.getDomain()->getDim());
    return Data(value, shape, what, expanded);
}

Data ComplexTensorC(DataTypes::cplx_t value, const FunctionSpace& what, bool expanded)
{
    DataTypes::ShapeType shape(2, what.getDomain()->getDim());
    escript::Data newdata = Data(value, shape, what, expanded);
    newdata.complicate();
    return newdata;
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

Data ComplexTensorFromObj(bp::object o, const FunctionSpace& what, bool expanded)
{
    // first try to get a double and route it to the other method
    try {
        double v = bp::extract<double>(o);
        return ComplexTensor(v, what, expanded);
    } catch(...) {
        PyErr_Clear();
    }
    // now try to get a complex and route to scalar factory
    try {
        DataTypes::cplx_t v = bp::extract<DataTypes::cplx_t>(o);
        return ComplexTensorC(v, what, expanded);
    } catch(...) {
        PyErr_Clear();
    }
    DataTypes::ShapeType shape(2, what.getDomain()->getDim());
    Data d(o, what, expanded);
    d.complicate();
    if (d.getDataPointShape() != shape) {
        throw DataException("ComplexTensorFromObj: Shape of tensor passed to function"
               " does not match the dimension of the domain.");
    }
    return d;
}

Data Tensor3(double value, const FunctionSpace& what, bool expanded)
{
    DataTypes::ShapeType shape(3, what.getDomain()->getDim());
    return Data(value, shape, what, expanded);
}

Data ComplexTensor3(double value, const FunctionSpace& what, bool expanded)
{
    DataTypes::ShapeType shape(3, what.getDomain()->getDim());
    escript::Data newdata = Data(value, shape, what, expanded);
    newdata.complicate();
    return newdata;
}

Data Tensor3C(DataTypes::cplx_t  value, const FunctionSpace& what, bool expanded)
{
    DataTypes::ShapeType shape(3, what.getDomain()->getDim());
    return Data(value, shape, what, expanded);
}

Data ComplexTensor3C(DataTypes::cplx_t  value, const FunctionSpace& what, bool expanded)
{
    DataTypes::ShapeType shape(3, what.getDomain()->getDim());
    escript::Data newdata = Data(value, shape, what, expanded);
    newdata.complicate();
    return newdata;
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

Data ComplexTensor3FromObj(bp::object o, const FunctionSpace& what, bool expanded)
{
    // first try to get a double and route it to the other method
    try {
        double v = bp::extract<double>(o);
        return ComplexTensor3(v, what, expanded);
    } catch(...) {
        PyErr_Clear();
    }
    // first try to get a complex and route it to the other method
    try {
        DataTypes::cplx_t v = bp::extract<DataTypes::cplx_t>(o);
        return ComplexTensor3C(v, what, expanded);
    } catch(...) {
        PyErr_Clear();
    }
    DataTypes::ShapeType shape(3, what.getDomain()->getDim());
    Data d(o, what, expanded);
    d.complicate();
    if (d.getDataPointShape() != shape) {
        throw DataException("ComplexTensor3FromObj: Shape of tensor passed to "
                "function does not match the dimension of the domain.");
    }
    return d;
}

Data Tensor4(double value, const FunctionSpace& what, bool expanded)
{
    DataTypes::ShapeType shape(4, what.getDomain()->getDim());
    return Data(value, shape, what, expanded);
}

Data ComplexTensor4(double value, const FunctionSpace& what, bool expanded)
{
    DataTypes::ShapeType shape(4, what.getDomain()->getDim());
    escript::Data newdata = Data(value, shape, what, expanded);
    newdata.complicate();
    return newdata;
}

Data Tensor4C(DataTypes::cplx_t value, const FunctionSpace& what, bool expanded)
{
    DataTypes::ShapeType shape(4, what.getDomain()->getDim());
    return Data(value, shape, what, expanded);
}

Data ComplexTensor4C(DataTypes::cplx_t value, const FunctionSpace& what, bool expanded)
{
    DataTypes::ShapeType shape(4, what.getDomain()->getDim());
    escript::Data newdata = Data(value, shape, what, expanded);
    newdata.complicate();
    return newdata;
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

Data ComplexTensor4FromObj(bp::object o, const FunctionSpace& what, bool expanded)
{
    // first try to get a double and route it to the other method
    try {
        double v = bp::extract<double>(o);
        return ComplexTensor4(v, what, expanded);
    } catch(...) {
        PyErr_Clear();
    }
    // first try to get a double and route it to the other method
    try {
        DataTypes::cplx_t v = bp::extract<DataTypes::cplx_t>(o);
        return ComplexTensor4C(v, what, expanded);
    } catch(...) {
        PyErr_Clear();
    }
    DataTypes::ShapeType shape(4, what.getDomain()->getDim());
    Data d(o, what, expanded);
    d.complicate();
    if (d.getDataPointShape() != shape) {
        throw DataException("ComplexTensor4FromObj: Shape of tensor passed to function"
               " does not match the dimension of the domain.");
    }
    return d;
}

Data ComplexData(boost::python::object o, const FunctionSpace& what, bool expanded)
{
    // first try to get a double and route it to the other method
    try {
        double v = bp::extract<double>(o);
        std::complex<double> c(v,0);
        DataTypes::ShapeType shape;
        escript::Data newdata = Data(c, shape, what, expanded);
        newdata.complicate();
        return newdata;
    } catch(...) {
        PyErr_Clear();
    }
    // first try to get a double and route it to the other method
    try {
        DataTypes::cplx_t v = bp::extract<DataTypes::cplx_t>(o);
        DataTypes::ShapeType shape;
        escript::Data newdata = Data(v, shape, what, expanded);
        newdata.complicate();
        return newdata;
    } catch(...) {
        PyErr_Clear();
    }
    DataTypes::ShapeType shape(1, what.getDomain()->getDim());
    Data d(o, what, expanded);
    d.complicate();
    if (d.getDataPointShape() != shape) {
        throw DataException("ComplexData: Shape of tensor passed to function"
               " does not match the dimension of the domain.");
    }
    return d;
}

#ifdef ESYS_HAVE_HDF5
Data load_hdf5grp(const H5::Group h5_grp, const AbstractDomain& domain)
{
    Data out;
    int error = 0;
    std::string msg;
    JMPI mpiInfo(domain.getMPI());
    /* .. read meta data ... */
    try {
            uint rank = 0;
            int type=-1;
            int function_space_type=-1;
            uint h5_shape[DataTypes::maxRank];

            // .... read meta data ...
            H5::DataSet h5_meta_data=h5_grp.openDataSet("Meta");
            // .. rank:
            H5::Attribute h5_attr_rank(h5_meta_data.openAttribute("rank"));
            H5::DataType h5_type_rank(h5_attr_rank.getDataType());
            if ( h5_type_rank != H5::PredType::NATIVE_UINT ) {
                 throw DataException("Error - load_hdf5: illegal rank data type in HDF5 file.");
            }
            if ( h5_attr_rank.getStorageSize() != 1 * h5_type_rank.getSize() ) {
                 throw DataException("Error - load_hdf5: rank in HDF5 file needs to be a single value.");
            }
            h5_attr_rank.read(h5_type_rank, &rank);
            // .. type id:
            H5::Attribute h5_attr_type(h5_meta_data.openAttribute("type_id"));
            H5::DataType h5_type_type(h5_attr_type.getDataType());
            if ( h5_type_type != H5::PredType::NATIVE_INT ) {
                 throw DataException("Error - load_hdf5: illegal type_id data type in HDF5 file.");
            }
            if ( h5_attr_type.getStorageSize() !=  1 * h5_type_type.getSize() ) {
                 throw DataException("Error - load_hdf5: type_id  in HDF5 file  needs to be a single value.");
            }
            h5_attr_type.read(h5_type_type, &type);
            // .. functionspace type:
            H5::Attribute h5_attr_fstype(h5_meta_data.openAttribute("function_space_type"));
            H5::DataType h5_type_fstype(h5_attr_fstype.getDataType());
            if ( h5_type_fstype != H5::PredType::NATIVE_INT ) {
                 throw DataException("Error - load_hdf5: illegal function_space_type data type in HDF5 file.");
            }
            if ( h5_attr_fstype.getStorageSize() != 1 * h5_type_fstype.getSize() ) {
                 throw DataException("Error - load_hdf5: function_space_type in HDF5 file  needs to be a single value.");
            }
            h5_attr_fstype.read(h5_type_fstype, &function_space_type);
            // ... shape
            H5::DataType h5_type_shape(h5_meta_data.getDataType());
            if ( h5_type_shape != H5::PredType::NATIVE_UINT ) {
                 throw DataException("Error - load_hdf5: illegal shape data type in HDF5 file.");
            }
            if ( h5_meta_data.getStorageSize() != h5_type_shape.getSize()  * rank ) {
                 throw DataException("Error - load_hdf5: shape length  in HDF5 file needs to be equal to rank.");
            }
            DataTypes::ShapeType shape;
            shape.resize(rank);
            h5_meta_data.read(&h5_shape, h5_type_shape);
            int num_values_per_data_point = 1;
            for (uint i=0; i < rank; ++i) {
                shape[i] = h5_shape[i];
                num_values_per_data_point*=shape[i];
            }
            if (! ( num_values_per_data_point > 0) )
            {       // jump to outer catch
                throw DataException("Error - load_hdf5: illegal number of values per data point.");
            }
            // ... recover function space:
            if (!domain.isValidFunctionSpaceType(function_space_type))
            {       // jump to outer catch
                throw DataException("Error - load_hdf5: function space type code in HDF5 file is invalid for given domain.");
            }
            FunctionSpace function_space=FunctionSpace(domain.getPtr(), function_space_type);

            // ... recover data
            if (type == 0) {
                // constant data
               H5::DataSet h5_data_data =h5_grp.openDataSet("data");
               H5::DataType h5_type_data(h5_data_data.getDataType());
               out=Data(0, shape, function_space, false);
               if ( h5_type_data != H5::PredType::NATIVE_DOUBLE ) {
                 throw DataException("Error - load_hdf5: illegal value data type for constant data in HDF5 file.");
               }
               if ( out.getDataPointSize() * sizeof(DataTypes::real_t) != h5_data_data.getStorageSize() )
               {
                     throw DataException("Error - load_hdf5: insufficient data values for constant data  in HDF5 file.");
               }
               h5_data_data.read(&(out.getDataAtOffsetRW(out.getDataOffset(0,0), static_cast<DataTypes::real_t>(0))), h5_type_data);
            } else if (type == 1) {
               // tagged data

               // .... check the sample id order ...
               H5::DataSet h5_data_tags =h5_grp.openDataSet("tags");
               H5::DataType h5_type_tags(h5_data_tags.getDataType());

               const int ntags = h5_data_tags.getStorageSize() / h5_type_tags.getSize();
               if ( h5_type_tags != H5::PredType::NATIVE_INT ) {
                 throw DataException("Error - load_hdf5: illegal value data type in HDF5 file.");
               }

               if (  ! ( ntags > 0 ) )
               {
                     throw DataException("Error - load_hdf5: no tags found for tagged data in HDF5 file.");
               }
               std::vector<int> h5_tags(ntags);
               h5_data_tags.read(&h5_tags[0], h5_type_tags);

               H5::DataSet h5_data_data =h5_grp.openDataSet("data");
               H5::DataType h5_type_data(h5_data_data.getDataType());
               if ( h5_type_data != H5::PredType::NATIVE_DOUBLE ) {
                 throw DataException("Error - load_hdf5: illegal value data type for tagged data in HDF5 file.");
               }
               if ( ntags * num_values_per_data_point * sizeof(DataTypes::real_t) != h5_data_data.getStorageSize() )
               {
                     throw DataException("Error - load_hdf5: insufficient data values for tagged data in HDF5 file.");
               }
                DataTypes::RealVectorType values(num_values_per_data_point * ntags, 0., num_values_per_data_point * ntags);
                h5_data_data.read(&(values[0]), h5_type_data);
                DataTagged* taggedData=new DataTagged(function_space, shape, &h5_tags[0], values);
                out=Data(taggedData);
            } else if (type == 2) {
                // ... expanded data
                const int num_samples = function_space.getNumSamples();
                const int num_data_points_per_sample = function_space.getNumDataPointsPerSample();
                const DataTypes::dim_t* ids_p=function_space.borrowSampleReferenceIDs();
                // .... check the sample id order ...
                H5::DataSet h5_data_sample_id =h5_grp.openDataSet("sample_id");
                H5::DataType h5_type_sample_id(h5_data_sample_id.getDataType());
                #ifdef ESYS_INDEXTYPE_LONG
                    if ( h5_type_sample_id != H5::PredType::NATIVE_LONG ) {
                        throw DataException("Error - load_hdf5: illegal value data type for sample id for expanded data in HDF5 file.");
                    }
                #else
                    if ( h5_type_sample_id != H5::PredType::NATIVE_INT ) {
                        throw DataException("Error - load_hdf5: illegal value data type for sample id for expanded data in HDF5 file.");
                    }
                #endif
               if ( num_samples * sizeof(DataTypes::dim_t) != h5_data_sample_id.getStorageSize() )
               {
                     throw DataException("Error - load_hdf5: insufficient sample id values for expanded data in HDF5 file.");
               }
               std::vector<DataTypes::dim_t> H5_ids(num_samples);
               h5_data_sample_id.read(&H5_ids[0], h5_type_sample_id);
               // ..... check order .....
               int wrong_position=-1, local_wrong_position=-1, i;
               #pragma omp parallel private(local_wrong_position)
               {
                    local_wrong_position=-1;
                    #pragma omp for private(i) schedule(static)
                    for (i=0; i < num_samples; ++i) {
                           if (H5_ids[i]!=ids_p[i]) local_wrong_position=i;
                    }
                    #pragma omp critical
                    if (local_wrong_position>=0) wrong_position = local_wrong_position;
               }
               //... load data ....
               H5::DataSet h5_data_data =h5_grp.openDataSet("data");
               H5::DataType h5_type_data(h5_data_data.getDataType());
               out=Data(0, shape, function_space, true);
               if ( h5_type_data != H5::PredType::NATIVE_DOUBLE ) {
                 throw DataException("Error - load_hdf5: illegal value data type for expanded data in HDF5 file.");
               }
               if ( out.getDataPointSize() * num_samples * num_data_points_per_sample * sizeof(DataTypes::real_t) != h5_data_data.getStorageSize() )
               {
                     throw DataException("Error - load_hdf5: insufficient data values for expanded data in HDF5 file.");
               }
               h5_data_data.read(&(out.getDataAtOffsetRW(out.getDataOffset(0,0), static_cast<DataTypes::real_t>(0))), h5_type_data);
               // if order is not the same we try to reorder the data:
               if (wrong_position >= 0)
               {
                    try {
                        std::cout << "Information - load: start reordering data from HDF5 file" << std::endl;
                        out.borrowData()->reorderByReferenceIDs(&H5_ids[0]);
                    }
                    catch (std::exception&) {
                        throw DataException("load: unable to reorder data in HDF5 file.");
                    }
               }
            } else {
                throw DataException("Error - load_hdf5:: unknown data type in HDF5 file detected.");
            }
    }
    // catch failure caused by the H5File operations
    catch (H5::Exception& e)
    {
        error=1;
        e.printErrorStack();
        msg=e.getCDetailMsg();
    }
    catch (DataException& e) {
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
}
#endif  // ESYS_HAVE_HDF5


Data load_hdf5(const std::string fileName, const AbstractDomain& domain)
{
    Data out;
#ifdef ESYS_HAVE_HDF5
    JMPI mpiInfo(domain.getMPI());
    const std::string newFileName(mpiInfo->appendRankToFileName(fileName));
    H5::H5File h5_file(newFileName, H5F_ACC_RDONLY);

    if (h5_file.nameExists("Data"))
    {
        return load_hdf5grp(h5_file.openGroup("Data"), domain);

    }
    else if ( h5_file.nameExists("Data_Re") and h5_file.nameExists("Data_Im") )
    {
            Data data_re = load_hdf5grp(h5_file.openGroup("Data_Re"), domain);
            Data data_im = load_hdf5grp(h5_file.openGroup("Data_Im"), domain);
            return data_re +  data_im * Scalar(std::complex(0.0, 1.0), data_im.getFunctionSpace(), false );
    }
    else
    {
        throw DataException("load_hdf5: unable to load HDF neither as real nor complex data as expected groups are missing.");
    }
#else
    throw DataException("load_hdf5: not compiled with HDF5 (serial). Please contact your installation manager.");
#endif // ESYS_HAVE_HDF5
    return out;
}

bool loadConfigured()
{
#ifdef ESYS_HAVE_HDF5
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
