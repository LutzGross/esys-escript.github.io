//$Id$
/* 
 ************************************************************
 *          Copyright 2006 by ACcESS MNRF                   *
 *                                                          *
 *              http://www.access.edu.au                    *
 *       Primary Business: Queensland, Australia            *
 *  Licensed under the Open Software License version 3.0    *
 *     http://www.opensource.org/licenses/osl-3.0.php       *
 *                                                          *
 ************************************************************
*/

#if !defined escript_DataConstant_20040323_H
#define escript_DataConstant_20040323_H
#include "system_dep.h"

#include "DataAbstract.h"
#include "DataArrayView.h"

#include <boost/python/numeric.hpp>

namespace escript {

/**
   \brief
   DataConstant stores a single data point which represents the entire
   function space.

   Description:
   DataConstant stores a single data point which represents the entire
   function space.
*/
class DataConstant : public DataAbstract  {

 public:

  /**
     \brief
     Constructor for DataConstant objects.

     Description:
     Constructor for DataConstant objects.

     \param value - Input - Data value for a single point.
     \param what - Input - A description of what this data object represents.
  */
  ESCRIPT_DLL_API
  DataConstant(const boost::python::numeric::array& value,
               const FunctionSpace& what);

  /**
     \brief
     Copy constructor. Performs a deep copy.
  */
  ESCRIPT_DLL_API
  DataConstant(const DataConstant& other);

  /**
     \brief
     Alternative constructor for DataConstant objects.

     Description:
     Alternative Constructor for DataConstant objects.
     \param value - Input - Data value for a single point.
     \param what - Input - A description of what this data object represents.
  */
  ESCRIPT_DLL_API
  DataConstant(const DataArrayView& value,
               const FunctionSpace& what);

  /**
     \brief
     Alternative constructor for DataConstant objects.

     Description:
     Alternative Constructor for DataConstant objects.
     \param other - Input - Data object to copy from.
     \param region - Input - region to copy.
  */
  ESCRIPT_DLL_API
  DataConstant(const DataConstant& other,
               const DataArrayView::RegionType& region);

  /**
     \brief
     Alternative constructor for DataConstant objects.

     Description:
     Alternative Constructor for DataConstant objects.
     \param what - Input - A description of what this data object represents.
     \param shape - Input - the shape of each data-point.
     \param data - the data values for each data-point.
  */
  ESCRIPT_DLL_API
  DataConstant(const FunctionSpace& what,
               const DataArrayView::ShapeType &shape,
               const DataArrayView::ValueType &data);

  /**
     \brief
     Write the data as a string.
  */
  ESCRIPT_DLL_API
  std::string
  toString() const;

  /**
     \brief
     Return the offset for the given sample. This is a somewhat artificial notion
     but returns the offset in bytes for the given point into the container
     holding the point data. Only really necessary to avoid many DataArrayView
     objects.
     \param sampleNo - Input - sample number.
     \param dataPointNo - Input - data point number for the sample.
   */
  ESCRIPT_DLL_API
  virtual
  DataArrayView::ValueType::size_type
  getPointOffset(int sampleNo,
                 int dataPointNo) const;

  /**
     \brief
     Return a view into the data for the data point specified.
     \param sampleNo - Input - sample number.
     \param dataPointNo - Input - data point number for the sample.
  */
  ESCRIPT_DLL_API
  virtual
  DataArrayView
  getDataPoint(int sampleNo,
               int dataPointNo);

  /**
     \brief
     Return the number of doubles stored for the Data object.
  */
  ESCRIPT_DLL_API
  virtual
  DataArrayView::ValueType::size_type
  getLength() const;

  /**
     \brief
     Factory method that returns a newly created DataConstant object
     sliced from the specified region of this object.
     The caller is reponsible for managing the object created.
     \param region - Input - region to slice from this object.
  */
  ESCRIPT_DLL_API
  virtual
  DataAbstract*
  getSlice(const DataArrayView::RegionType& region) const;

  /**
     \brief
     Copy the specified region from the given value.
     \param value - Input - Data object to copy from.
     \param region - Input - Region to copy.
  */
  ESCRIPT_DLL_API
  virtual
  void
  setSlice(const DataAbstract* value,
           const DataArrayView::RegionType& region);

  /**
     \brief
     Reshape the data point if the data point is currently rank 0.
     The original data point value is used for all values of the new
     data point.
  */
  ESCRIPT_DLL_API
  void
  reshapeDataPoint(const DataArrayView::ShapeType& shape);

  /**
    \brief
    Archive the underlying data values to the file referenced
    by ofstream. A count of the number of values expected to be written
    is provided as a cross-check.

    The return value indicates success (0) or otherwise (1).
  */
  ESCRIPT_DLL_API
  int
  archiveData(std::ofstream& archiveFile,
              const DataArrayView::ValueType::size_type noValues) const;

  /**
    \brief
    Extract the number of values specified by noValues from the file
    referenced by ifstream to the underlying data structure.

    The return value indicates success (0) or otherwise (1).
  */
  ESCRIPT_DLL_API
  int
  extractData(std::ifstream& archiveFile,
              const DataArrayView::ValueType::size_type noValues);

  /**
     \brief
     Computes a symmetric matrix (A + AT) / 2

     \param ev - Output - symmetric matrix

  */
  ESCRIPT_DLL_API
  virtual void
  symmetric(DataAbstract* ev);

  /**
     \brief
     Computes a nonsymmetric matrix (A - AT) / 2

     \param ev - Output - nonsymmetric matrix

  */
  ESCRIPT_DLL_API
  virtual void
  nonsymmetric(DataAbstract* ev);

  /**
     \brief
     Computes the trace of a matrix

     \param ev - Output - trace of matrix

  */
  ESCRIPT_DLL_API
  virtual void
  trace(DataAbstract* ev, int axis_offset);

  /**
     \brief
     Transpose each data point of this Data object around the given axis.

     \param ev - Output - transpose of matrix

  */
  ESCRIPT_DLL_API
  virtual void
  transpose(DataAbstract* ev, int axis_offset);

  /**
     \brief
     swaps components axis_offset and axis_offset+1

     \param ev - Output - swapped components

  */
  ESCRIPT_DLL_API
  virtual void
  swap(DataAbstract* ev, int axis_offset);


  /**
     \brief
     solves the eigenvalue problem this*V=ev*V for the eigenvalues ev

     \param ev - Output - eigenvalues in increasing order at each data point

  */
  ESCRIPT_DLL_API
  virtual void
  eigenvalues(DataAbstract* ev);

  /**
     \brief
     solves the eigenvalue problem this*V=ev*V for the eigenvalues ev and eigenvectors V

     \param ev - Output - eigenvalues in increasing order at each data point
     \param V - Output - corresponding eigenvectors. They are normalized such that their length is one
                         and the first nonzero component is positive.
     \param tol - Input - eigenvalue with relative distance tol are treated as equal.

  */

  ESCRIPT_DLL_API
  virtual void
  eigenvalues_and_eigenvectors(DataAbstract* ev,DataAbstract* V,const double tol=1.e-13);


 protected:

 private:
  //
  // the actual data
  DataArrayView::ValueType m_data;

};

} // end of namespace
#endif
