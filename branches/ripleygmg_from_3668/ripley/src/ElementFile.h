
/*******************************************************
*
* Copyright (c) 2003-2011 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/

#ifndef __RIPLEY_ELEMENTFILE_H__
#define __RIPLEY_ELEMENTFILE_H__

#include <ripley/Ripley.h>
#include <ripley/IndexList.h>

#ifdef USE_NETCDF
class NcFile;
#endif

namespace escript {
class Data;
}

namespace ripley {

class RipleyDomain;
class ElementFile;
typedef boost::shared_ptr<ElementFile> ElementFile_ptr;

typedef enum {
    InvalidElementType = -1,
    Point1 = 0,
    Line2 = 1,
    Rec4 = 2,
    Hex8 = 3,
    Line2Face = 4,
    Rec4Face = 5,
    Hex8Face = 6
} ElementTypeId;

/**
   \brief
   An ElementFile stores information about one of the element types of a
   ripley domain.
*/
class ElementFile
{
public:
    /**
       \brief
       Returns the element type ID corresponding to the type name s.
    */
    static ElementTypeId ripleyTypeFromString(const std::string &s);

    /**
       \brief
       Constructor with dimensionality and MPI info structure.
    */
    RIPLEY_DLL_API
    ElementFile(ElementTypeId type, Esys_MPIInfo *mpiInfo);

    /**
       \brief
       Destructor.
    */
    RIPLEY_DLL_API
    ~ElementFile();

#ifdef USE_NETCDF
    /**
       \brief
       Reads ID, tags, color etc. from a netCDF file.
    */
    void readFromNetCDF(NcFile &dataFile, dim_t numElements, const std::string &name);

    /**
       \brief
       Stores ID, tags, color etc. into a netCDF file.
       'name' should be unique among the element files in order to be able
       to recover the data using readFromNetCDF.
    */
    void dumpToNetCDF(NcFile &dataFile, const std::string &name);
#endif

    /**
       \brief
       Returns a textual representation of the element type.
    */
    std::string getName() const;

    /**
       \brief
       Returns the type identifier for this element file.
    */
    ElementTypeId getTypeId() const { return m_type; }

    /**
       \brief
       Returns the number of nodes per element.
    */
    dim_t getNumNodes() const { return 1 << m_localDim; } /*Line=2,Rec=4,Hex=8*/

    /**
       \brief
       Returns the element dimensionality, e.g. Line2 has dim 2 in 2D and 3D.
    */
    dim_t getNumLocalDim() const { return m_localDim; }

    /**
       \brief
       Returns the number of spatial dimensions.
    */
    dim_t getNumDim() const { return m_numDim; }

    /**
       \brief
       Returns the active MPI info structure.
    */
    Esys_MPIInfo* getMPIInfo() const { return m_mpiInfo; }

    /**
       \brief
       Returns a reference to the vector of element IDs.
    */
    const IndexVector &getIdVector() const { return m_id; }

    /**
       \brief
       Returns a reference to the vector of element tags.
    */
    const IndexVector &getTagVector() const { return m_tag; }

    /**
       \brief
       Returns a reference to the vector of element colors.
    */
    const IndexVector &getColorVector() const { return m_color; }

    /**
       \brief
       Returns a reference to the vector of element owners.
    */
    const RankVector &getOwnerVector() const { return m_owner; }

    /**
       \brief
       Returns a reference to the vector of element nodes.
    */
    const IndexVector &getNodes() const { return m_nodes; }

    /**
       \brief
       Returns a reference to the vector of tags that are in use.
    */
    const IndexVector &getTagsInUse() const { return m_tagsInUse; }

    /**
       \brief
       Returns the number of elements in this element file.
    */
    dim_t getNumElements() const { return m_id.size(); }

    /**
       \brief
       Sets tags to newTag where mask>0.
    */
    void setTags(int newTag, const escript::Data &mask);

    /**
       \brief
       Refreshes the list of tags that are in use.
    */
    void updateTagsInUse();

    /**
       \brief
       Sets tags to newTag where mask>0.
    */
    dim_t getNumberOfTagsInUse() const { return m_tagsInUse.size(); }

    /**
       \brief
       Swaps the contents of the vectors with the node data vectors.
    */
    void swapEntries(IndexVector &idIn, IndexVector &tagIn,
                     IndexVector &colorIn, RankVector &ownerIn,
                     IndexVector &nodesIn);

    /**
       \brief
       Marks the used nodes with offset, i.e. sets mask[n-offset]=1 for all
       nodes n used by this element file.
    */
    void markNodes(IndexVector &mask, index_t offset);

    /**
       \brief
       Assigns new node reference numbers to elements.
       If k is the old node, the new node is newNode[k-offset].
    */
    void relabelNodes(const IndexVector &newNode, index_t offset);

    /**
       \brief
       Returns the minimum and maximum node IDs used in this element file.
    */
    IndexPair getNodeRange() const;

    /**
       \brief
       Reorders the elements so that the elements are stored close to the
       nodes.
    */
    void optimizeOrdering();

    /**
       \brief
       Tries to reduce the number of colors used to color the elements.
    */
    void createColoring(const IndexVector &degreesOfFreedom);

    /**
       \brief
       Redistributes the elements including overlap.
    */
    void distributeByRankOfDOF(const RankVector &mpiRankOfDOF, const IndexVector &id);

    /**
       \brief
       Inserts elements into an index matrix.
    */
    void insertIntoIndexMatrix(IndexMatrix &matrix, const IndexVector &rowMap,
                               const IndexVector &colMap, index_t firstRow=-1,
                               index_t lastRow=-1);

    /**
       \brief
       Inserts elements into an index matrix omitting the main diagonal.
    */
    void insertIntoIndexMatrixNoMainDiagonal(IndexMatrix &matrix,
            const IndexVector &rowMap, const IndexVector &colMap,
            index_t firstRow, index_t lastRow);

private:
    // MPI information
    Esys_MPIInfo *m_mpiInfo;

    // element type identifier
    ElementTypeId m_type;

    // number of spatial dimensions in the domain
    dim_t m_numDim;

    // number of dimensions of the element type
    dim_t m_localDim;

    // m_nodes[INDEX2(k,i,numNodes)] is the k-th node in element i
    IndexVector m_nodes;

    // m_owner[i] is the process ID that owns element i
    RankVector m_owner;

    // m_id[i] is the unique ID of element i
    IndexVector m_id;

    // m_tag[i] is the tag of element i
    IndexVector m_tag;

    // m_color[i] is the color of element i. Elements with the same color don't
    // share a node so they can be processed simultaneously.
    IndexVector m_color;

    // minimum color number
    index_t m_minColor;

    // maximum color number
    index_t m_maxColor;

    // array of tags which are actually used
    IndexVector m_tagsInUse;
};

} // end of namespace ripley

#endif // __RIPLEY_ELEMENTFILE_H__

