
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

#include <weipa/FinleyElements.h>
#include <weipa/NodeData.h>

#ifndef VISIT_PLUGIN
#include <dudley/CppAdapter/MeshAdapter.h>
#include <finley/CppAdapter/MeshAdapter.h>
#else
#define ABS(X) ((X)>0?(X):-(X))
#endif

#include <iostream>

#if USE_NETCDF
#include <netcdfcpp.h>
#endif

#if USE_SILO
#include <silo.h>
#endif


using namespace std;

// The following arrays contain indices to convert unsupported element
// types into supported ones (lines, triangles, quads, hexahedrons)

static const size_t line3indices[2*2] = {
    0, 2,
    2, 1
};
static const size_t tri6indices[4*3] = {
    0, 3, 5,
    5, 4, 2,
    3, 1, 4,
    4, 5, 3
};
static const size_t rec8indices[6*3] = {
    0, 4, 7,
    4, 1, 5,
    5, 2, 6,
    6, 3, 7,
    7, 5, 6,
    7, 4, 5
};
static const size_t rec9indices[4*4] = {
    0, 4, 8, 7,
    4, 1, 5, 8,
    7, 8, 6, 3,
    8, 5, 2, 6
};
static const size_t tet10indices[8*4] = {
    6, 4, 0, 7,
    6, 5, 4, 8,
    5, 1, 4, 8,
    9, 8, 7, 3,
    2, 5, 6, 9,
    8, 9, 5, 6,
    6, 7, 9, 8,
    6, 4, 7, 8
};
static const size_t hex20indices[36*3] = {
     0,  8, 12,   8,  1, 13,  13,  5, 16,
    16,  4, 12,   8, 13, 16,   8, 16, 12,
     1,  9, 13,   9,  2, 14,  14,  6, 17,
    17,  5, 13,   9, 14, 17,   9, 17, 13,
     2, 10, 14,  10,  3, 15,  15,  7, 18,
    18, 14,  6,  10, 15, 18,  10, 18, 14,
     3, 11, 15,  11,  0, 12,  12,  4, 19,
    19,  7, 15,  11, 12, 19,  11, 19, 15,
     4, 16, 19,  16,  5, 17,  17,  6, 18,
    18,  7, 19,  16, 17, 18,  16, 18, 19,
     3, 10, 11,  10,  2,  9,   9,  1,  8,
     8,  0, 11,  10,  9,  8,  10,  8, 11
};
static const size_t hex27indices[8*8] = {
     0,  8, 20, 11, 12, 21, 26, 24,
     8,  1,  9, 20, 21, 13, 22, 26,
    11, 20, 10,  3, 24, 26, 23, 15,
    20,  9,  2, 10, 26, 22, 14, 23,
    12, 21, 26, 24,  4, 16, 25, 19,
    21, 13, 22, 26, 16,  5, 17, 25,
    24, 26, 23, 15, 19, 25, 18,  7,
    26, 22, 14, 23, 25, 17,  6, 18
};

namespace weipa {
    
//
// Constructor
//
FinleyElements::FinleyElements(const string& elementName, FinleyNodes_ptr nodeData)
    : originalMesh(nodeData), name(elementName), numElements(0),
      numGhostElements(0), nodesPerElement(0),
      type(ZONETYPE_UNKNOWN), finleyTypeId(Finley_NoRef), elementFactor(1)
{
    nodeMesh.reset(new FinleyNodes(name));
}

//
// Copy constructor
//
FinleyElements::FinleyElements(const FinleyElements& e)
{
    name = e.name;
    numElements = e.numElements;
    numGhostElements = e.numGhostElements;
    type = e.type;
    finleyTypeId = e.finleyTypeId;
    nodesPerElement = e.nodesPerElement;
    elementFactor = e.elementFactor;
    originalMesh = e.originalMesh;
    if (e.nodeMesh)
        nodeMesh.reset(new FinleyNodes(*e.nodeMesh));
    else
        nodeMesh.reset(new FinleyNodes(name));

    nodes = e.nodes;
    ID = e.ID;
    color = e.color;
    tag = e.tag;
    owner = e.owner;

    if (e.reducedElements)
        reducedElements = FinleyElements_ptr(
                new FinleyElements(*e.reducedElements));
}

//
//
//
bool FinleyElements::initFromDudley(const Dudley_ElementFile* dudleyFile)
{
#ifndef VISIT_PLUGIN
    numElements = dudleyFile->numElements;

    if (numElements > 0) {
        nodesPerElement = dudleyFile->numNodes;

        int* iPtr;
   
        iPtr = dudleyFile->Nodes;
        nodes.clear();
        nodes.insert(nodes.end(), numElements*nodesPerElement, 0);
        copy(iPtr, iPtr+numElements*nodesPerElement, nodes.begin());

        iPtr = dudleyFile->Color;
        color.clear();
        color.insert(color.end(), numElements, 0);
        copy(iPtr, iPtr+numElements, color.begin());

        iPtr = dudleyFile->Id;
        ID.clear();
        ID.insert(ID.end(), numElements, 0);
        copy(iPtr, iPtr+numElements, ID.begin());

        iPtr = dudleyFile->Owner;
        owner.clear();
        owner.insert(owner.end(), numElements, 0);
        copy(iPtr, iPtr+numElements, owner.begin());

        iPtr = dudleyFile->Tag;
        tag.clear();
        tag.insert(tag.end(), numElements, 0);
        copy(iPtr, iPtr+numElements, tag.begin());

        FinleyElementInfo f = getDudleyTypeInfo(dudleyFile->etype);
        type = f.elementType;
        elementFactor = f.elementFactor;
        if (elementFactor > 1 || f.reducedElementSize != nodesPerElement)
            buildReducedElements(f);

        buildMeshes();
    }
    return true;

#else // VISIT_PLUGIN
    return false;
#endif
}

//
//
//
bool FinleyElements::initFromFinley(const Finley_ElementFile* finleyFile)
{
#ifndef VISIT_PLUGIN
    numElements = finleyFile->numElements;

    if (numElements > 0) {
        nodesPerElement = finleyFile->numNodes;

        int* iPtr;
   
        iPtr = finleyFile->Nodes;
        nodes.clear();
        nodes.insert(nodes.end(), numElements*nodesPerElement, 0);
        copy(iPtr, iPtr+numElements*nodesPerElement, nodes.begin());

        iPtr = finleyFile->Color;
        color.clear();
        color.insert(color.end(), numElements, 0);
        copy(iPtr, iPtr+numElements, color.begin());

        iPtr = finleyFile->Id;
        ID.clear();
        ID.insert(ID.end(), numElements, 0);
        copy(iPtr, iPtr+numElements, ID.begin());

        iPtr = finleyFile->Owner;
        owner.clear();
        owner.insert(owner.end(), numElements, 0);
        copy(iPtr, iPtr+numElements, owner.begin());

        iPtr = finleyFile->Tag;
        tag.clear();
        tag.insert(tag.end(), numElements, 0);
        copy(iPtr, iPtr+numElements, tag.begin());

        finleyTypeId = finleyFile->referenceElementSet->referenceElement
            ->Type->TypeId;
        FinleyElementInfo f = getFinleyTypeInfo(finleyTypeId);
        type = f.elementType;
        elementFactor = f.elementFactor;
        if (elementFactor > 1 || f.reducedElementSize != nodesPerElement)
            buildReducedElements(f);

        if (f.useQuadNodes) {
            CoordArray quadNodes;
            int numQuadNodes;
            Finley_ShapeFunction* sf = finleyFile->referenceElementSet
                ->referenceElement->Parametrization;
            numQuadNodes = sf->numQuadNodes;
            for (int i=0; i<f.quadDim; i++) {
                double* srcPtr = sf->QuadNodes+i;
                float* c = new float[numQuadNodes];
                quadNodes.push_back(c);
                for (int j=0; j<numQuadNodes; j++, srcPtr+=f.quadDim) {
                    *c++ = (float) *srcPtr;
                }
            }
            quadMask = buildQuadMask(quadNodes, numQuadNodes);
            for (int i=0; i<f.quadDim; i++)
                delete[] quadNodes[i];
            quadNodes.clear();


            // now the reduced quadrature
            sf = finleyFile->referenceElementSet
                ->referenceElementReducedQuadrature->Parametrization;
            numQuadNodes = sf->numQuadNodes;
            for (int i=0; i<f.quadDim; i++) {
                double* srcPtr = sf->QuadNodes+i;
                float* c = new float[numQuadNodes];
                quadNodes.push_back(c);
                for (int j=0; j<numQuadNodes; j++, srcPtr+=f.quadDim) {
                    *c++ = (float) *srcPtr;
                }
            }
            reducedQuadMask = buildQuadMask(quadNodes, numQuadNodes);
            for (int i=0; i<f.quadDim; i++)
                delete[] quadNodes[i];
            quadNodes.clear();
        }

        buildMeshes();
    }
    return true;

#else // VISIT_PLUGIN
    return false;
#endif
}

//
// Reads element data from given NetCDF file
//
bool FinleyElements::readFromNc(NcFile* ncfile)
{
#if USE_NETCDF
    string num_str("num_");
    num_str += name;

    NcAtt* att = ncfile->get_att(num_str.c_str());
    numElements = att->as_int(0);

    // Only attempt to read further if there are any elements.
    // Having no elements is not an error.
    if (numElements > 0) {
        att = ncfile->get_att((num_str + string("_numNodes")).c_str());
        nodesPerElement = att->as_int(0);

        nodes.insert(nodes.end(), numElements*nodesPerElement, 0);
        NcVar* var = ncfile->get_var((name + string("_Nodes")).c_str());
        var->get(&nodes[0], numElements, nodesPerElement);

        color.insert(color.end(), numElements, 0);
        var = ncfile->get_var((name + string("_Color")).c_str());
        var->get(&color[0], numElements);

        ID.insert(ID.end(), numElements, 0);
        var = ncfile->get_var((name + string("_Id")).c_str());
        var->get(&ID[0], numElements);

        owner.insert(owner.end(), numElements, 0);
        var = ncfile->get_var((name + string("_Owner")).c_str());
        var->get(&owner[0], numElements);

        tag.insert(tag.end(), numElements, 0);
        var = ncfile->get_var((name + string("_Tag")).c_str());
        var->get(&tag[0], numElements);

        att = ncfile->get_att((name + string("_TypeId")).c_str());
        finleyTypeId = (Finley_ElementTypeId)att->as_int(0);
        FinleyElementInfo f = getFinleyTypeInfo(finleyTypeId);
        type = f.elementType;
        elementFactor = f.elementFactor;
        if (f.elementFactor > 1 || f.reducedElementSize != nodesPerElement)
            buildReducedElements(f);

        // if we don't link with finley we can't get the quadrature nodes
        // and hence cannot interpolate data properly
#ifndef VISIT_PLUGIN
        if (f.useQuadNodes) {
            att = ncfile->get_att("order");
            int order = att->as_int(0);
            att = ncfile->get_att("reduced_order");
            int reduced_order = att->as_int(0);
            finley::ReferenceElementSetWrapper wrapper(finleyTypeId, order,
                    reduced_order);
            Finley_ReferenceElementSet* refElements = wrapper.getElementSet();

            CoordArray quadNodes;
            int numQuadNodes;
            Finley_ShapeFunction* sf = refElements->referenceElement
                ->Parametrization;
            numQuadNodes = sf->numQuadNodes;
            for (int i=0; i<f.quadDim; i++) {
                double* srcPtr = sf->QuadNodes+i;
                float* c = new float[numQuadNodes];
                quadNodes.push_back(c);
                for (int j=0; j<numQuadNodes; j++, srcPtr+=f.quadDim) {
                    *c++ = (float) *srcPtr;
                }
            }
            quadMask = buildQuadMask(quadNodes, numQuadNodes);
            for (int i=0; i<f.quadDim; i++)
                delete[] quadNodes[i];
            quadNodes.clear();

            // now the reduced quadrature
            sf = refElements->referenceElementReducedQuadrature->Parametrization;
            numQuadNodes = sf->numQuadNodes;
            for (int i=0; i<f.quadDim; i++) {
                double* srcPtr = sf->QuadNodes+i;
                float* c = new float[numQuadNodes];
                quadNodes.push_back(c);
                for (int j=0; j<numQuadNodes; j++, srcPtr+=f.quadDim) {
                    *c++ = (float) *srcPtr;
                }
            }
            reducedQuadMask = buildQuadMask(quadNodes, numQuadNodes);
            for (int i=0; i<f.quadDim; i++)
                delete[] quadNodes[i];
            quadNodes.clear();
        }
#endif // VISIT_PLUGIN

        buildMeshes();
    }

    return true; 
#else // !USE_NETCDF
    return false;
#endif
}

//
//
//
StringVec FinleyElements::getMeshNames() const
{
    StringVec res;
    if (nodeMesh)
        res.push_back(nodeMesh->getName());
    if (reducedElements) {
        StringVec rNames = reducedElements->getMeshNames();
        if (!rNames.empty())
            res.insert(res.end(), rNames.begin(), rNames.end());
    }
    return res;
}

//
//
//
StringVec FinleyElements::getVarNames() const
{
    StringVec res;
    res.push_back(name + string("_Color"));
    res.push_back(name + string("_Id"));
    res.push_back(name + string("_Owner"));
    res.push_back(name + string("_Tag"));
    return res;
}

//
//
//
const IntVec& FinleyElements::getVarDataByName(const string varName) const
{
    if (varName == name+string("_Color"))
        return color;
    else if (varName == name+string("_Id"))
        return ID;
    else if (varName == name+string("_Owner")) {
        return owner;
    } else if (varName == name+string("_Tag"))
        return tag;
    else if (reducedElements)
        return reducedElements->getVarDataByName(varName);
    else
        throw "Invalid variable name";
}

//
//
//
const QuadMaskInfo& FinleyElements::getQuadMask(int fsCode) const
{
    if (fsCode == FINLEY_REDUCED_ELEMENTS ||
            fsCode == FINLEY_REDUCED_FACE_ELEMENTS ||
            fsCode == FINLEY_REDUCED_CONTACT_ELEMENTS_1) {
        return reducedQuadMask;
    } else {
        return quadMask;
    }
}

//
//
//
void FinleyElements::reorderArray(IntVec& v, const IntVec& idx,
                               int elementsPerIndex)
{
    IntVec newArray(v.size());
    IntVec::iterator arrIt = newArray.begin();
    IntVec::const_iterator idxIt;
    if (elementsPerIndex == 1) {
        for (idxIt=idx.begin(); idxIt!=idx.end(); idxIt++) {
            *arrIt++ = v[*idxIt];
        }
    } else {
        for (idxIt=idx.begin(); idxIt!=idx.end(); idxIt++) {
            int i = *idxIt;
            copy(&v[i*elementsPerIndex], &v[(i+1)*elementsPerIndex], arrIt);
            arrIt += elementsPerIndex;
        }
    }
    v.swap(newArray);
}

//
//
//
void FinleyElements::buildReducedElements(const FinleyElementInfo& f)
{
    // create the node list for the new element type
    IntVec reducedNodes(f.reducedElementSize*numElements, 0);

    IntVec::iterator reducedIt = reducedNodes.begin();
    IntVec::const_iterator origIt;
    for (origIt=nodes.begin(); origIt!=nodes.end();
         origIt+=nodesPerElement)
    {
        copy(origIt, origIt+f.reducedElementSize, reducedIt);
        reducedIt += f.reducedElementSize;
    }

    if (f.elementFactor > 1) {
        // replace each element by multiple smaller ones which will be the
        // new 'full' elements, whereas the original ones are the reduced
        // elements, e.g.:
        //
        // new reduced:         new full:
        //   _________           _________
        //  |         |         |    |    |
        //  |         |    >    |____|____|
        //  |         |    >    |    |    |
        //  |_________|         |____|____|
        //

        // create the reduced elements which are basically a copy of the
        // current elements
        reducedElements = FinleyElements_ptr(new FinleyElements(
                    "Reduced"+name, originalMesh));
        reducedElements->nodes = reducedNodes;
        reducedElements->numElements = numElements;
        reducedElements->type = f.reducedElementType;
        reducedElements->nodesPerElement = f.reducedElementSize;
        reducedElements->owner = owner;
        reducedElements->color = color;
        reducedElements->ID = ID;
        reducedElements->tag = tag;

        // now update full element data
        IntVec fullNodes(f.elementSize*f.elementFactor*numElements);
        IntVec::iterator cellIt = fullNodes.begin();

        owner.clear();
        color.clear();
        ID.clear();
        tag.clear();
        for (int i=0; i < numElements; i++) {
            owner.insert(owner.end(), f.elementFactor, reducedElements->owner[i]);
            color.insert(color.end(), f.elementFactor, reducedElements->color[i]);
            ID.insert(ID.end(), f.elementFactor, reducedElements->ID[i]);
            tag.insert(tag.end(), f.elementFactor, reducedElements->tag[i]);
            for (int j=0; j < f.elementFactor*f.elementSize; j++)
                *cellIt++ = nodes[nodesPerElement*i+f.multiCellIndices[j]];
        }

        nodes.swap(fullNodes);
        nodesPerElement = f.elementSize;
        numElements *= f.elementFactor;

    } else {
        // we merely converted element types and don't need reduced elements
        // so just replace node list and type
        nodes.swap(reducedNodes);
        nodesPerElement = f.reducedElementSize;
        type = f.reducedElementType;
    }
}

//
//
//
IntVec FinleyElements::prepareGhostIndices(int ownIndex)
{
    IntVec indexArray;
    numGhostElements = 0;
    
    // move indices of "ghost zones" to the end to be able to reorder
    // data accordingly
    for (int i=0; i<numElements; i++) {
        if (owner[i] == ownIndex)
            indexArray.push_back(i);
    }

    for (int i=0; i<numElements; i++) {
        if (owner[i] != ownIndex) {
            numGhostElements++;
            indexArray.push_back(i);
        }
    }
    return indexArray;
}

//
//
//
void FinleyElements::reorderGhostZones(int ownIndex)
{
    IntVec indexArray = prepareGhostIndices(ownIndex);

    // move "ghost data" to the end of the arrays
    if (numGhostElements > 0) {
        reorderArray(nodes, indexArray, nodesPerElement);
        reorderArray(owner, indexArray, 1);
        reorderArray(color, indexArray, 1);
        reorderArray(ID, indexArray, 1);
        reorderArray(tag, indexArray, 1);
    }

    if (reducedElements)
        reducedElements->reorderGhostZones(ownIndex);
}

//
//
//
void FinleyElements::removeGhostZones(int ownIndex)
{
    reorderGhostZones(ownIndex);

    if (numGhostElements > 0) {
        numElements -= numGhostElements;
        nodes.resize(numElements*nodesPerElement);
        owner.resize(numElements);
        color.resize(numElements);
        ID.resize(numElements);
        tag.resize(numElements);
        numGhostElements = 0;
    }

    if (reducedElements)
        reducedElements->removeGhostZones(ownIndex);
}

//
//
//
void FinleyElements::buildMeshes()
{
    // build a new mesh containing only the required nodes
    if (numElements > 0) {
        if (nodeMesh && nodeMesh->getNumNodes() > 0) {
            FinleyNodes_ptr newMesh(new FinleyNodes(nodeMesh, nodes, name));
            nodeMesh.swap(newMesh);
        } else {
            nodeMesh.reset(new FinleyNodes(originalMesh, nodes, name));
        }
#ifdef _DEBUG
        cout << nodeMesh->getName() << " has " << nodeMesh->getNumNodes()
            << " nodes and " << numElements << " elements" << endl;
#endif
    }

    if (reducedElements)
        reducedElements->buildMeshes();
}

//
//
//
void FinleyElements::writeConnectivityVTK(ostream& os)
{
    if (numElements > 0) {
        const IntVec& gNI = nodeMesh->getGlobalNodeIndices();
        IntVec::const_iterator it;
        int count = 1;
        for (it=nodes.begin(); it!=nodes.end(); it++, count++) {
            os << gNI[*it];
            if (count % nodesPerElement == 0)
                os << endl;
            else
                os << " ";
        }
    }
}

#if USE_SILO
//
//
//
inline int toSiloElementType(int type)
{
    switch (type) {
        case ZONETYPE_BEAM: return DB_ZONETYPE_BEAM;
        case ZONETYPE_HEX: return DB_ZONETYPE_HEX;
        case ZONETYPE_POLYGON: return DB_ZONETYPE_POLYGON;
        case ZONETYPE_QUAD: return DB_ZONETYPE_QUAD;
        case ZONETYPE_TET: return DB_ZONETYPE_TET;
        case ZONETYPE_TRIANGLE: return DB_ZONETYPE_TRIANGLE;
    }
    return 0;
}
#endif

//
//
//
bool FinleyElements::writeToSilo(DBfile* dbfile, const string& siloPath,
                                 const StringVec& labels,
                                 const StringVec& units)
{
#if USE_SILO
    if (numElements == 0)
        return true;

    int ret;

    if (siloPath != "") {
        ret = DBSetDir(dbfile, siloPath.c_str());
        if (ret != 0)
            return false;
    }

    // write out the full mesh in any case
    nodeMesh->setSiloPath(siloPath);
    string siloMeshNameStr = nodeMesh->getFullSiloName();
    const char* siloMeshName = siloMeshNameStr.c_str();
    int arraylen = numElements * nodesPerElement;
    int eltype = toSiloElementType(type);

    string varName = name + string("_zones");
    ret = DBPutZonelist2(dbfile, varName.c_str(), numElements,
            nodeMesh->getNumDims(), &nodes[0], arraylen, 0, 0,
            numGhostElements, &eltype, &nodesPerElement, &numElements, 1, NULL);
    if (ret == 0) {
        CoordArray& coordbase = const_cast<CoordArray&>(nodeMesh->getCoords());
        DBoptlist* optList = NULL;
        int nOpts = labels.size()+units.size();
        if (nOpts>0) {
            optList = DBMakeOptlist(nOpts);
            if (labels.size()>0)
                DBAddOption(optList, DBOPT_XLABEL, (void*)labels[0].c_str());
            if (labels.size()>1)
                DBAddOption(optList, DBOPT_YLABEL, (void*)labels[1].c_str());
            if (labels.size()>2)
                DBAddOption(optList, DBOPT_ZLABEL, (void*)labels[2].c_str());
            if (units.size()>0)
                DBAddOption(optList, DBOPT_XUNITS, (void*)units[0].c_str());
            if (units.size()>1)
                DBAddOption(optList, DBOPT_YUNITS, (void*)units[1].c_str());
            if (units.size()>2)
                DBAddOption(optList, DBOPT_ZUNITS, (void*)units[2].c_str());
        }
        ret = DBPutUcdmesh(dbfile, siloMeshName,
                nodeMesh->getNumDims(), NULL, &coordbase[0],
                nodeMesh->getNumNodes(), numElements, varName.c_str(),
                /*"facelist"*/NULL, DB_FLOAT, optList);
        if (optList) {
            DBFreeOptlist(optList);
        }
    }
    
    // Point mesh is useful for debugging
    if (0) {
        CoordArray& coordbase = const_cast<CoordArray&>(nodeMesh->getCoords());
        DBPutPointmesh(dbfile, "/pointmesh", nodeMesh->getNumDims(),
                &coordbase[0], nodeMesh->getNumNodes(), DB_FLOAT, NULL);
    }

    if (ret != 0)
        return false;

    // write out the element-centered variables
    varName = name + string("_Color");
    ret = DBPutUcdvar1(dbfile, varName.c_str(), siloMeshName,
            (float*)&color[0], numElements, NULL, 0, DB_INT, DB_ZONECENT,
            NULL);
    if (ret == 0) {
        varName = name + string("_Id");
        ret = DBPutUcdvar1(dbfile, varName.c_str(), siloMeshName,
            (float*)&ID[0], numElements, NULL, 0, DB_INT, DB_ZONECENT,
            NULL);
    }
    if (ret == 0) {
        varName = name + string("_Owner");
        ret = DBPutUcdvar1(dbfile, varName.c_str(), siloMeshName,
            (float*)&owner[0], numElements, NULL, 0, DB_INT, DB_ZONECENT,
            NULL);
    }
    if (ret == 0) {
        varName = name + string("_Tag");
        ret = DBPutUcdvar1(dbfile, varName.c_str(), siloMeshName,
            (float*)&tag[0], numElements, NULL, 0, DB_INT, DB_ZONECENT,
            NULL);
    }

    if (reducedElements) {
        reducedElements->writeToSilo(dbfile, siloPath, labels, units);
    }

    // "Elements" is a special case
    if (name == "Elements") {
        nodeMesh->writeToSilo(dbfile);
    }

    return (ret == 0);

#else // !USE_SILO
    return false;
#endif
}

//
//
//
FinleyElementInfo FinleyElements::getDudleyTypeInfo(Dudley_ElementTypeId typeId)
{
    FinleyElementInfo ret;
    ret.multiCellIndices = NULL;
    ret.elementFactor = 1;
    ret.useQuadNodes = false;
    ret.quadDim = 0;

    switch (typeId) {
        case Dudley_Line2Face://untested
        case Dudley_Point1://untested
            cerr << "WARNING: Dudley type " <<typeId<< " is untested!" << endl;
            ret.elementSize = 1;
            ret.elementType = ZONETYPE_POLYGON;
            break;

        case Dudley_Tri3Face://untested
            cerr << "WARNING: Dudley type " <<typeId<< " is untested!" << endl;
        case Dudley_Line2:
            ret.elementSize = ret.reducedElementSize = 2;
            ret.elementType = ret.reducedElementType = ZONETYPE_BEAM;
            break;

        case Dudley_Tet4Face://untested
            cerr << "WARNING: Dudley type " <<typeId<< " is untested!" << endl;
        case Dudley_Tri3:
            ret.elementSize = ret.reducedElementSize = 3;
            ret.elementType = ret.reducedElementType = ZONETYPE_TRIANGLE;
            break;

        case Dudley_Tet4:
            ret.elementSize = ret.reducedElementSize = 4;
            ret.elementType = ret.reducedElementType = ZONETYPE_TET;
            break;

        default:
            cerr << "WARNING: Unknown Dudley Type " << typeId << endl;
            break;
    }
    return ret;
}

//
//
//
FinleyElementInfo FinleyElements::getFinleyTypeInfo(Finley_ElementTypeId typeId)
{
    FinleyElementInfo ret;
    ret.multiCellIndices = NULL;
    ret.elementFactor = 1;
    ret.useQuadNodes = false;
    ret.quadDim = 0;

    switch (typeId) {
        case Finley_Point1_Contact://untested
        case Finley_Line2Face_Contact://untested
        case Finley_Line3Face_Contact://untested
        case Finley_Line2Face://untested
        case Finley_Line3Face://untested
        case Finley_Point1://untested
            cerr << "WARNING: Finley type " <<typeId<< " is untested!" << endl;
            ret.elementSize = 1;
            ret.elementType = ZONETYPE_POLYGON;
            break;

        case Finley_Tri3Face://untested
            cerr << "WARNING: Finley type " <<typeId<< " is untested!" << endl;
        case Finley_Tri3Face_Contact:
        case Finley_Line2_Contact:
        case Finley_Rec4Face_Contact:
        case Finley_Rec4Face:
        case Finley_Line2:
            ret.elementSize = ret.reducedElementSize = 2;
            ret.elementType = ret.reducedElementType = ZONETYPE_BEAM;
            break;

        case Finley_Line3Macro:
            ret.useQuadNodes = true;
            ret.quadDim = 1;
            // fall through
        case Finley_Line3:
            ret.multiCellIndices = line3indices;
            ret.elementFactor = 2;
            // fall through
        case Finley_Line3_Contact:
        case Finley_Tri6Face_Contact://untested
        case Finley_Rec8Face_Contact:
        case Finley_Tri6Face://untested
        case Finley_Rec8Face:
            //VTK_QUADRATIC_EDGE
            ret.elementSize = ret.reducedElementSize = 2;
            ret.elementType = ret.reducedElementType = ZONETYPE_BEAM;
            break;

        case Finley_Tet4Face_Contact://untested
        case Finley_Tet4Face://untested
            cerr << "WARNING: Finley type " <<typeId<< " is untested!" << endl;
        case Finley_Tri3_Contact:
        case Finley_Tri3:
            ret.elementSize = ret.reducedElementSize = 3;
            ret.elementType = ret.reducedElementType = ZONETYPE_TRIANGLE;
            break;

        case Finley_Rec4_Contact:
        case Finley_Hex8Face_Contact:
        case Finley_Hex8Face:
        case Finley_Rec4:
            ret.elementSize = ret.reducedElementSize = 4;
            ret.elementType = ret.reducedElementType = ZONETYPE_QUAD;
            break;

        case Finley_Rec9:
        case Finley_Rec9Macro:
            ret.useQuadNodes = true;
            ret.quadDim = 2;
            ret.multiCellIndices = rec9indices;
            ret.elementFactor = 4;
            // fall through
        case Finley_Rec9_Contact:
            ret.elementSize = ret.reducedElementSize = 4;
            ret.elementType = ret.reducedElementType = ZONETYPE_QUAD;
            break;

        case Finley_Tet4:
            ret.elementSize = ret.reducedElementSize = 4;
            ret.elementType = ret.reducedElementType = ZONETYPE_TET;
            break;

        case Finley_Tri6:
        case Finley_Tri6Macro:
            ret.useQuadNodes = true;
            ret.quadDim = 2;
            ret.multiCellIndices = tri6indices;
            ret.elementFactor = 4;
            // fall through
        case Finley_Tri6_Contact:
        case Finley_Tet10Face_Contact://untested
        case Finley_Tet10Face://untested
            //VTK_QUADRATIC_TRIANGLE
            ret.elementSize = ret.reducedElementSize = 3;
            ret.elementType = ret.reducedElementType = ZONETYPE_TRIANGLE;
            break;

        case Finley_Rec8:
            ret.multiCellIndices = rec8indices;
            ret.elementFactor = 6;
            // fall through
        case Finley_Hex20Face:
        case Finley_Rec8_Contact:
        case Finley_Hex20Face_Contact:
            //VTK_QUADRATIC_QUAD
            ret.elementSize = 3;
            ret.elementType = ZONETYPE_TRIANGLE; 
            ret.reducedElementSize = 4;
            ret.reducedElementType = ZONETYPE_QUAD;
            break;

        case Finley_Tet10:
        case Finley_Tet10Macro:
            //VTK_QUADRATIC_TETRA
            ret.useQuadNodes = true;
            ret.quadDim = 3;
            ret.multiCellIndices = tet10indices;
            ret.elementFactor = 8;
            ret.elementSize = ret.reducedElementSize = 4;
            ret.elementType = ret.reducedElementType = ZONETYPE_TET;
            break;

        case Finley_Hex20:
            //VTK_QUADRATIC_HEXAHEDRON
            ret.multiCellIndices = hex20indices;
            ret.elementFactor = 36;
            ret.elementSize = 3;
            ret.elementType = ZONETYPE_TRIANGLE;
            ret.reducedElementSize = 8;
            ret.reducedElementType = ZONETYPE_HEX;
            break;

        case Finley_Hex27:
        case Finley_Hex27Macro:
            ret.useQuadNodes = true;
            ret.quadDim = 3;
            ret.multiCellIndices = hex27indices;
            ret.elementFactor = 8;
            // fall through
        case Finley_Hex8:
            ret.elementSize = ret.reducedElementSize = 8;
            ret.elementType = ret.reducedElementType = ZONETYPE_HEX;
            break;

        default:
            cerr << "WARNING: Unknown Finley Type " << typeId << endl;
            break;
    }
    return ret;
}

/////////////////////////////////
// Helpers for buildQuadMask() //
/////////////////////////////////

// returns true if |x-c| <= r, false otherwise
inline bool inside1D(float x, float c, float r)
{
    return (ABS(x-c) <= r);
}

// returns true if |x-cx| <= r and |y-cy| <= r, false otherwise
inline bool inside2D(float x, float y, float cx, float cy, float r)
{
    return (inside1D(x, cx, r) && inside1D(y, cy, r));
}

// returns true if |x-cx| <= r and |y-cy| <= r and |z-cz| <= r, false otherwise
inline bool inside3D(float x, float y, float z,
                     float cx, float cy, float cz, float r)
{
    return (inside2D(x, y, cx, cy, r) && inside1D(z, cz, r));
}

// returns true if d1 and d2 have the same sign or at least one of them is
// close to 0, false otherwise
inline bool sameSide(float d1, float d2)
{
    const float TOL = 1.e-8f;
    return (ABS(d1) < TOL || ABS(d2) < TOL || d1*d2>=0.);
}

// computes the determinant of the 4x4 matrix given by its elements m_ij
static float det4x4(float m_00, float m_01, float m_02, float m_03,
                    float m_10, float m_11, float m_12, float m_13,
                    float m_20, float m_21, float m_22, float m_23,
                    float m_30, float m_31, float m_32, float m_33)
{
    float det1 = m_12 * m_23 - m_22 * m_13;
    float det2 = m_11 * m_23 - m_21 * m_13;
    float det3 = m_11 * m_22 - m_21 * m_12;
    float det4 = m_10 * m_23 - m_20 * m_13;
    float det5 = m_10 * m_22 - m_20 * m_12;
    float det6 = m_10 * m_21 - m_20 * m_11;
    return -m_30 * (m_01 * det1 - m_02 * det2 + m_03 * det3) +
            m_31 * (m_00 * det1 - m_02 * det4 + m_03 * det5) -
            m_32 * (m_00 * det2 - m_01 * det4 + m_03 * det6) +
            m_33 * (m_00 * det3 - m_01 * det5 + m_02 * det6);
}

// returns true if point (x,y,z) is in or on the tetrahedron given by its
// corner points p0, p1, p2 and p3, false otherwise.
static bool pointInTet(float x, float y, float z,
                       const float* p0, const float* p1,
                       const float* p2, const float* p3)
{
    float d0 = det4x4(
            p0[0], p0[1], p0[2], 1.f,
            p1[0], p1[1], p1[2], 1.f,
            p2[0], p2[1], p2[2], 1.f,
            p3[0], p3[1], p3[2], 1.f);
    float d1 = det4x4(
                x,     y,     z, 1.f,
            p1[0], p1[1], p1[2], 1.f,
            p2[0], p2[1], p2[2], 1.f,
            p3[0], p3[1], p3[2], 1.f);
    float d2 = det4x4(
            p0[0], p0[1], p0[2], 1.f,
                x,     y,     z, 1.f,
            p2[0], p2[1], p2[2], 1.f,
            p3[0], p3[1], p3[2], 1.f);
    float d3 = det4x4(
            p0[0], p0[1], p0[2], 1.f,
            p1[0], p1[1], p1[2], 1.f,
                x,     y,     z, 1.f,
            p3[0], p3[1], p3[2], 1.f);
    float d4 = det4x4(
            p0[0], p0[1], p0[2], 1.f,
            p1[0], p1[1], p1[2], 1.f,
            p2[0], p2[1], p2[2], 1.f,
                x,     y,     z, 1.f);

    return (sameSide(d0,d1) && sameSide(d1,d2) &&
            sameSide(d2,d3) && sameSide(d3,d4));
}

// returns true if point (x,y) is in or on the triangle given by its corner
// points p0, p1 and p2, false otherwise
static bool pointInTri(float x, float y,
                       const float* p0, const float* p1, const float* p2)
{
    const float TOL = 1.e-8f;
    float v0[2] = { p2[0]-p0[0], p2[1]-p0[1] };
    float v1[2] = { p1[0]-p0[0], p1[1]-p0[1] };
    float v2[2] = { x - p0[0],   y - p0[1] };

    float dot00 = v0[0]*v0[0]+v0[1]*v0[1];
    float dot01 = v0[0]*v1[0]+v0[1]*v1[1];
    float dot02 = v0[0]*v2[0]+v0[1]*v2[1];
    float dot11 = v1[0]*v1[0]+v1[1]*v1[1];
    float dot12 = v1[0]*v2[0]+v1[1]*v2[1];
    float invDenom = dot00*dot11 - dot01*dot01;
    if (ABS(invDenom) < TOL) invDenom = TOL;
    invDenom = 1.f/invDenom;
    float u = (dot11*dot02 - dot01*dot12) * invDenom;
    float v = (dot00*dot12 - dot01*dot02) * invDenom;
    return (u>=0.f) && (v>=0.f) && (u+v<=1.f);
}

//
//
//
QuadMaskInfo FinleyElements::buildQuadMask(const CoordArray& qnodes, int numQNodes)
{
    QuadMaskInfo qmi;
    if (numQNodes == 0)
        return qmi;

    if (finleyTypeId == Finley_Line3Macro) {
        for (int i=0; i<elementFactor; i++) {
            const float bounds[] = { .25, .75 };
            IntVec m(numQNodes, 0);
            int hits = 0;
            for (size_t j=0; j<numQNodes; j++) {
                if (inside1D(qnodes[0][j], bounds[i], .25)) {
                    m[j] = 1;
                    hits++;
                }
            }
            qmi.mask.push_back(m);
            if (hits == 0)
                qmi.factor.push_back(1);
            else
                qmi.factor.push_back(hits);
        }
    } else if ((finleyTypeId == Finley_Tri6) || (finleyTypeId == Finley_Tri6Macro)) {
        for (int i=0; i<elementFactor; i++) {
            const float bounds[][2] = { { 0., 0. }, { 1., 0. },
                                        { 0., 1. }, { .5, 0. },
                                        { .5, .5 }, { 0., .5 } };
            const size_t* nodeIdx = &tri6indices[i*nodesPerElement];
            IntVec m(numQNodes, 0);
            int hits = 0;
            for (size_t j=0; j<numQNodes; j++) {
                // check if point j is in triangle i
                if (pointInTri(qnodes[0][j], qnodes[1][j],
                        bounds[nodeIdx[0]], bounds[nodeIdx[1]],
                        bounds[nodeIdx[2]])) {
                    m[j] = 1;
                    hits++;
                }
            }
            if (hits == 0) {
                // if an element does not contain any quadrature points we
                // simply average over all data points within that element
                m = IntVec(numQNodes, 1);
                qmi.factor.push_back(numQNodes);
            } else {
                qmi.factor.push_back(hits);
            }
            qmi.mask.push_back(m);
        }
    } else if ((finleyTypeId == Finley_Tet10) || (finleyTypeId == Finley_Tet10Macro)) {
        for (int i=0; i<elementFactor; i++) {
            const float bounds[][3] = {
                { 0., 0., 0. },
                { 1., 0., 0. },
                { 0., 0., 1. },
                { 0., 1., 0. },
                { .5, 0., 0. },
                { .5, 0., .5 },
                { 0., 0., .5 },
                { 0., .5, 0. },
                { .5, .5, 0. },
                { 0., .5, .5 }
            };
            // need to reorder the elements
            const size_t elNumIdx[] = { 0,1,2,4,3,5,6,7 };
            const size_t* nodeIdx = &tet10indices[elNumIdx[i]*nodesPerElement];
            IntVec m(numQNodes, 0);
            int hits = 0;
            for (size_t j=0; j<numQNodes; j++) {
                // check if point j is in tetrahedron i
                if (pointInTet(qnodes[0][j], qnodes[1][j], qnodes[2][j],
                               bounds[nodeIdx[0]], bounds[nodeIdx[1]],
                               bounds[nodeIdx[2]], bounds[nodeIdx[3]])) {
                    m[j] = 1;
                    hits++;
                }
            }
            if (hits == 0) {
                // if an element does not contain any quadrature points we
                // simply average over all data points within that element
                m = IntVec(numQNodes, 1);
                qmi.factor.push_back(numQNodes);
            } else {
                qmi.factor.push_back(hits);
            }
            qmi.mask.push_back(m);
        }
    } else if ((finleyTypeId == Finley_Rec9) || (finleyTypeId == Finley_Rec9Macro)) {
        for (int i=0; i<elementFactor; i++) {
            const float bounds[][2] = { { 0.25, 0.25 }, { 0.75, 0.25 },
                                        { 0.25, 0.75 }, { 0.75, 0.75 } };
            IntVec m(numQNodes, 0);
            int hits = 0;
            for (size_t j=0; j<numQNodes; j++) {
                if (inside2D(qnodes[0][j], qnodes[1][j],
                        bounds[i][0], bounds[i][1], .25)) {
                    m[j] = 1;
                    hits++;
                }
            }
            qmi.mask.push_back(m);
            if (hits == 0)
                qmi.factor.push_back(1);
            else
                qmi.factor.push_back(hits);
        }
    } else if ((finleyTypeId == Finley_Hex27) || (finleyTypeId == Finley_Hex27Macro) ){
        for (int i=0; i<elementFactor; i++) {
            const float bounds[][3] = {
                { 0.25, 0.25, 0.25 }, { 0.75, 0.25, 0.25 },
                { 0.25, 0.75, 0.25 }, { 0.75, 0.75, 0.25 },
                { 0.25, 0.25, 0.75 }, { 0.75, 0.25, 0.75 },
                { 0.25, 0.75, 0.75 }, { 0.75, 0.75, 0.75 } };
            IntVec m(numQNodes, 0);
            int hits = 0;
            for (size_t j=0; j<numQNodes; j++) {
                if (inside3D(qnodes[0][j], qnodes[1][j], qnodes[2][j],
                        bounds[i][0], bounds[i][1], bounds[i][2], .25)) {
                    m[j] = 1;
                    hits++;
                }
            }
            qmi.mask.push_back(m);
            if (hits == 0)
                qmi.factor.push_back(1);
            else
                qmi.factor.push_back(hits);
        }
    }
    return qmi;
}

} // namespace weipa

