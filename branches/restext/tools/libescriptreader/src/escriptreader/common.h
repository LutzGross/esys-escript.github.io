
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

//
// common.h
//
#ifndef __ESCRIPTREADER_COMMON_H__
#define __ESCRIPTREADER_COMMON_H__

#include <string>
#include <vector>
#include <map>

namespace EscriptReader {

typedef std::vector<std::string> StringVec;
typedef std::vector<float> FloatVec;
typedef std::vector<int> IntVec;
typedef std::vector<float*> CoordArray;
typedef std::map<int, size_t> IndexMap;

} // namespace EscriptReader

#endif // __ESCRIPTREADER_COMMON_H__

