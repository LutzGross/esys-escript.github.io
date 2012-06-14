
/*******************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


#if !defined escript_UnaryOp_20040315_H
#define escript_UnaryOp_20040315_H
#include "system_dep.h"

#include "DataConstant.h"
#include "DataTagged.h"
#include "DataExpanded.h"
#include "DataTypes.h"

namespace escript {

/**
   \brief
   Perform the given unary operation on each data point of the given Data object.
   Called by Data::unaryOp.
   Calls DataArrayView::unaryOp.
   For DataExpanded objects, operation is done in parallel.
   \param data Input/Output - The data.
   \param operation Input - The operation to perform.
*/

template <class UnaryFunction>
inline
void
unaryOp(DataExpanded& data,
        UnaryFunction operation)
{
  int i,j;
  DataTypes::ValueType::size_type numDPPSample=data.getNumDPPSample();
  DataTypes::ValueType::size_type numSamples=data.getNumSamples();
  DataTypes::ValueType& left=data.getVectorRW();
  const DataTypes::ShapeType& shape=data.getShape();
  #pragma omp parallel for private(i,j) schedule(static)
  for (i=0;i<numSamples;i++) {
    for (j=0;j<numDPPSample;j++) {
      DataMaths::unaryOp(left,shape,data.getPointOffset(i,j),operation);
    }
  }
}

template <class UnaryFunction>
inline
void
unaryOp(DataTagged& data,
        UnaryFunction operation)
{
  // perform the operation on each tagged value
  const DataTagged::DataMapType& lookup=data.getTagLookup();
  DataTagged::DataMapType::const_iterator i;
  DataTagged::DataMapType::const_iterator lookupEnd=lookup.end();
  DataTypes::ValueType& left=data.getVectorRW();
  const DataTypes::ShapeType& shape=data.getShape();
  for (i=lookup.begin();i!=lookupEnd;i++) {
    DataMaths::unaryOp(left,shape,i->second,operation);
  }
  // perform the operation on the default value
  DataMaths::unaryOp(left,shape,data.getDefaultOffset(),operation);
}

template <class UnaryFunction>
inline
void
unaryOp(DataConstant& data,
        UnaryFunction operation)
{
  DataMaths::unaryOp(data.getVectorRW(),data.getShape(),0,operation);
}

} // end of namespace
#endif
