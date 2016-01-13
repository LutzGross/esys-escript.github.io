
/*******************************************************
*
* Copyright (c) 2008-2009 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/** \file BufferGroup.cpp */

#include "BufferGroup.h"

using namespace escript;

BufferGroup::BufferGroup(size_t buffersize,size_t numbuffs)
	: m_vec(numbuffs*buffersize,0,1),m_numbuffs(numbuffs),m_step(buffersize)
{
}

BufferGroup::~BufferGroup()
{	// DataVector can deallocate itself
}

