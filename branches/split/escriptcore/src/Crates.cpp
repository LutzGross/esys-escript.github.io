/*****************************************************************************
*
* Copyright (c) 2014 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

#include "Crates.h"
#include "SplitWorldException.h"

using namespace escript;

AbstractCrate::AbstractCrate(const std::string& name, const std::string& opname)
  : label(name), op(opname)
{
}

AbstractCrate::~AbstractCrate()
{}

const std::string& AbstractCrate::getLabel()
{
    return label;
}

DataCrate::DataCrate(const std::string& name, const std::string& opname)
  : AbstractCrate(name, opname)
{
    if (!OpLegal(opname))
    {
	std::string res=std::string("<")+opname+"> is not a legal merge operation for DataCrate."; 
	throw SplitWorldException(res.c_str());
    }
  
}


DataCrate::~DataCrate()
{
}

bool DataCrate::OpLegal(const std::string& opname)
{
    if (opname=="readonly")
    {
	return true;
    }
    return false;
}
