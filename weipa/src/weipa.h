
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#ifndef __WEIPA_H__
#define __WEIPA_H__

#ifndef VISIT_PLUGIN
#include <escript/DataTypes.h>
#endif

#include <string>
#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>

#ifdef ESYS_MPI
#define WEIPA_HAVE_MPI 1
#endif

#if WEIPA_HAVE_MPI
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

