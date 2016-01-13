
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/

#include "DataReady.h"

namespace escript
{

DataReady::DataReady(const FunctionSpace& what, const ShapeType& shape, bool isDataEmpty)
	:parent(what,shape,isDataEmpty)
{
}


DataReady_ptr 
DataReady::resolve()
{
	return boost::dynamic_pointer_cast<DataReady>(this->getPtr());
}



}
