
/*****************************************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/

#ifndef __WEIPA_H__
#define __WEIPA_H__

#include <string>
#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>

#ifdef ESYS_MPI
#define HAVE_MPI 1
#endif

#if HAVE_MPI
#include <mpi.h>
#endif

#define WEIPA_DLL_API

#ifdef _WIN32
#   ifndef WEIPA_STATIC_LIB
#       undef WEIPA_DLL_API
#       ifdef WEIPA_EXPORTS
#           define WEIPA_DLL_API __declspec(dllexport)
#       else
#           define WEIPA_DLL_API __declspec(dllimport)
#       endif
#   endif
#endif

namespace weipa {

class DataVar;
class DomainChunk;
class ElementData;
class EscriptDataset;
class NodeData;

typedef std::vector<float>       FloatVec;
typedef std::vector<int>         IntVec;
typedef std::vector<std::string> StringVec;
typedef std::vector<float*>      CoordArray;
typedef std::map<int, size_t>    IndexMap;

typedef boost::shared_ptr<DataVar>        DataVar_ptr;
typedef boost::shared_ptr<DomainChunk>    DomainChunk_ptr;
typedef boost::shared_ptr<ElementData>    ElementData_ptr;
typedef boost::shared_ptr<EscriptDataset> EscriptDataset_ptr;
typedef boost::shared_ptr<NodeData>       NodeData_ptr;

typedef boost::shared_ptr<const DomainChunk>    const_DomainChunk_ptr;
typedef boost::shared_ptr<const EscriptDataset> const_EscriptDataset_ptr;

} // namespace weipa

#endif // __WEIPA_H__

