
/*******************************************************
*
* Copyright (c) 2008 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/** \file BufferGroup.h */

#ifndef BUFFERGROUP_H
#define BUFFERGROUP_H

#include "DataVector.h"
#ifdef _OPENMP
#include <omp.h>
#endif

namespace escript
{

class BufferGroup
{
public:
	BufferGroup(size_t buffersize,size_t numbuffs);
	~BufferGroup();
	DataVector& getBuffer(size_t buffnum);
	size_t getOffset(size_t buffnum);

private:
	DataVector m_vec;
	size_t m_numbuffs;
	size_t m_step;
};

inline DataVector& BufferGroup::getBuffer(size_t buffnum)
{
	return m_vec;
}

inline size_t BufferGroup::getOffset(size_t buffnum)
{
	return m_step*buffnum;
}

} // end namespace

#endif	// BUFFERGROUP_H