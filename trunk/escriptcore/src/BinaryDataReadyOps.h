
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
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
		     escript::ESFunction operation);
void binaryOpDataTCT(DataTagged& result, const DataConstant& left, const DataTagged& right, 
		     escript::ESFunction operation);
void binaryOpDataTTC(DataTagged& result, const DataTagged& left, const DataConstant& right, 
		     escript::ESFunction operation);
void binaryOpDataTTT(DataTagged& result, const DataTagged& left, const DataTagged& right, 
		     escript::ESFunction operation);
void binaryOpDataEEC(DataExpanded& result, const DataExpanded& left, const DataConstant& right, 
		     escript::ESFunction operation);
void binaryOpDataECE(DataExpanded& result, const DataConstant& left, const DataExpanded& right, 
		     escript::ESFunction operation);
void binaryOpDataEEE(DataExpanded& result, const DataExpanded& left, const DataExpanded& right, 
		     escript::ESFunction operation);
void binaryOpDataETE(DataExpanded& result, const DataTagged& left, const DataExpanded& right, 
		     escript::ESFunction operation);
void binaryOpDataEET(DataExpanded& result, const DataExpanded& left, const DataTagged& right, 
 		     escript::ESFunction operation);

} // end of namespace

#endif // __ESCRIPT_BINARYOP_H__

