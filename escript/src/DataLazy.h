
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


#if !defined escript_DataLazy_20081008_H
#define escript_DataLazy_20081008_H
#include "system_dep.h"

#include "DataAbstract.h"
//#include "DataTypes.h"
//#include "FunctionSpace.h"

#include <string>

namespace escript {

enum ES_optype
{
	UNKNOWNOP=0,
	IDENTITY=1
};

const std::string&
opToString(ES_optype op);

/**
\class escript::DataLazy
\brief Wraps an expression tree of other DataObjects.
The values of DataPoints are computed when requested rather than all at once.
*/

class DataLazy : public DataAbstract
{

typedef DataAbstract parent;
typedef DataTypes::ValueType ValueType;
typedef DataTypes::ShapeType ShapeType;

public:
  ESCRIPT_DLL_API
  DataLazy(DataAbstract_ptr p);

  ESCRIPT_DLL_API
  DataLazy(DataAbstract_ptr left, DataAbstract_ptr right, ES_optype op);

  ESCRIPT_DLL_API
  ~DataLazy();



  /**
  \brief Compute all data points in the expression tree
  */
  ESCRIPT_DLL_API
  DataReady_ptr resolve();

  ESCRIPT_DLL_API
  std::string
  toString() const;

  ESCRIPT_DLL_API
  DataAbstract* 
  deepCopy();


  /**
     \brief
     Return the number of doubles that would be stored for this Data object if it were resolved.
  */
  ESCRIPT_DLL_API
  ValueType::size_type
  getLength() const;


  ESCRIPT_DLL_API
  DataAbstract*
  getSlice(const DataTypes::RegionType& region) const;


  DataTypes::ValueType::size_type 
  getPointOffset(int sampleNo,
                 int dataPointNo) const;

private:
  DataAbstract_ptr m_left, m_right;
  ES_optype m_op;
  size_t length;	// number of values represented by the operation
};

}
#endif
