
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

#ifndef __ESCRIPT_BINARYOP_H__
#define __ESCRIPT_BINARYOP_H__

#include "system_dep.h"
#include "DataTypes.h"
#include "DataConstant.h"
#include "DataExpanded.h"
#include "DataVectorOps.h"
#include "DataTagged.h"

/**
\file BinaryDataReadyOps.h 
\brief Describes binary operations performed on instances of DataAbstract.

For operations on DataVector see DataMaths.h.
For operations on double* see LocalOps.h.
*/

namespace escript {

void binaryOpDataCCC(DataConstant& result, const DataConstant& left, const DataConstant& right, 
		     escript::ES_optype operation);
void binaryOpDataTCT(DataTagged& result, const DataConstant& left, const DataTagged& right, 
		     escript::ES_optype operation);
void binaryOpDataTTC(DataTagged& result, const DataTagged& left, const DataConstant& right, 
		     escript::ES_optype operation);
void binaryOpDataTTT(DataTagged& result, const DataTagged& left, const DataTagged& right, 
		     escript::ES_optype operation);
void binaryOpDataEEC(DataExpanded& result, const DataExpanded& left, const DataConstant& right, 
		     escript::ES_optype operation);
void binaryOpDataECE(DataExpanded& result, const DataConstant& left, const DataExpanded& right, 
		     escript::ES_optype operation);
void binaryOpDataEEE(DataExpanded& result, const DataExpanded& left, const DataExpanded& right, 
		     escript::ES_optype operation);
void binaryOpDataETE(DataExpanded& result, const DataTagged& left, const DataExpanded& right, 
		     escript::ES_optype operation);
void binaryOpDataEET(DataExpanded& result, const DataExpanded& left, const DataTagged& right, 
 		     escript::ES_optype operation);

} // end of namespace

#endif // __ESCRIPT_BINARYOP_H__

