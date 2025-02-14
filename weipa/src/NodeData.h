
/*****************************************************************************
*
* Copyright (c) 2003-2020 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014-2017 by Centre for Geoscience Computing (GeoComp)
* Development from 2019 by School of Earth and Environmental Sciences
**
*****************************************************************************/

#ifndef __WEIPA_NODEDATA_H__
#define __WEIPA_NODEDATA_H__

#include <weipa/weipa.h>
#include <ostream>

namespace weipa {

/// \brief
class NodeData
{
public:
    /// \brief Writes coordinates to a stream in VTK text format.
    virtual void writeCoordinatesVTK(std::ostream& os, int ownIndex) = 0;

    /// \brief Returns a vector with the mesh variable names.
    virtual StringVec getVarNames() const = 0;

    /// \brief Returns the name of this node mesh.
    virtual std::string getName() const = 0;

    /// \brief Returns full Silo mesh name, e.g. "/block0000/Nodes".
    virtual std::string getFullSiloName() const = 0;

    /// \brief Returns the node ID array.
    virtual const IntVec& getNodeIDs() const = 0;

    /// \brief Returns the node distribution array
    virtual const IntVec& getNodeDistribution() const = 0;

    /// \brief Returns the global node index array.
    virtual const IntVec& getGlobalNodeIndices() const = 0;

    /// \brief Returns the coordinates of the mesh nodes.
    virtual const CoordArray& getCoords() const = 0;

    /// \brief Returns the dimensionality of this mesh (2 or 3).
    virtual int getNumDims() const = 0;

    /// \brief Returns the number of mesh nodes.
    virtual int getNumNodes() const = 0;

    /// \brief Returns the total number of mesh nodes for a distributed mesh.
    virtual int getGlobalNumNodes() const = 0;

protected:
    /// \brief Virtual destructor
    virtual ~NodeData() {}
};

} // namespace weipa

#endif // __WEIPA_NODEDATA_H__

