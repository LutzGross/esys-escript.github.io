
/*******************************************************
*
* Copyright (c) 2003-2010 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/

#ifndef __WEIPA_H__
#define __WEIPA_H__

#include <string>
#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>

#ifdef PASO_MPI
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
class ElementData;
class EscriptDataset;
class FinleyMesh;
class NodeData;

typedef std::vector<std::string> StringVec;
typedef std::vector<float> FloatVec;
typedef std::vector<int> IntVec;
typedef std::vector<float*> CoordArray;
typedef std::map<int, size_t> IndexMap;

typedef boost::shared_ptr<DataVar> DataVar_ptr;
typedef boost::shared_ptr<ElementData> ElementData_ptr;
typedef boost::shared_ptr<FinleyMesh> FinleyMesh_ptr;
typedef boost::shared_ptr<NodeData> NodeData_ptr;
typedef boost::shared_ptr<EscriptDataset> EscriptDataset_ptr;

} // namespace weipa

#endif // __WEIPA_H__

