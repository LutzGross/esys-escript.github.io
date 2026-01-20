
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

#ifndef __WEIPA_DOMAINCHUNK_H__
#define __WEIPA_DOMAINCHUNK_H__

#include <weipa/weipa.h>

class DBfile;

namespace escript {
    class AbstractDomain;
}

namespace weipa {

typedef enum {
    NODE_CENTERED=0,
    ZONE_CENTERED
} Centering;


/// \brief Abstract base class for weipa's interface to an Escript domain or
///        one chunk thereof if domain decomposition was used.
///
/// Implementations of this class load or convert data from an Escript domain
/// with all meshes for weipa.
///
/// Note that this class is not MPI aware, that is if domain decomposition
/// was used then one instance of this class will hold one 'chunk' of the
/// full domain. See the EscriptDataset class for how to process full domains.
class DomainChunk
{
public:
    /// \brief Initialises the domain using an escript domain instance.
    virtual bool initFromEscript(const escript::AbstractDomain* domain) = 0;

    /// \brief Reads the domain from a dump file.
    virtual bool initFromFile(const std::string& filename) = 0;

    /// \brief Writes the domain to a Silo file.
    virtual bool writeToSilo(DBfile* dbfile, const std::string& pathInSilo,
                             const StringVec& labels,
                             const StringVec& units,
                             bool writeMeshData) = 0;

    /// \brief Reorders elements so that 'ghost' elements (i.e. those that
    ///        do not belong to ownIndex) appear last.
    virtual void reorderGhostZones(int ownIndex) = 0;

    /// \brief Removes 'ghost' elements and nodes.
    virtual void removeGhostZones(int ownIndex) = 0;

    /// \brief Returns the names of all meshes within this domain.
    virtual StringVec getMeshNames() const = 0;

    /// \brief Returns the names of all 'special' domain variables.
    virtual StringVec getVarNames() const = 0;

    /// \brief Returns element data with given name.
    virtual ElementData_ptr getElementsByName(const std::string& name) const=0;

    /// \brief Returns the node mesh with given name.
    virtual NodeData_ptr getMeshByName(const std::string& name) const = 0;

    /// \brief Creates and returns a variable with domain data.
    virtual DataVar_ptr getDataVarByName(const std::string& name) const = 0;

    /// \brief Returns whether data on given function space is node or cell
    ///        centered
    virtual Centering getCenteringForFunctionSpace(int fsCode) const = 0;

    /// \brief Returns the node mesh for given function space code.
    virtual NodeData_ptr getMeshForFunctionSpace(int fsCode) const = 0;

    /// \brief Returns the element data for given function space code.
    virtual ElementData_ptr getElementsForFunctionSpace(int fsCode) const = 0;

    /// \brief Returns a pointer to the full nodes.
    virtual NodeData_ptr getNodes() const = 0;

    /// \brief Returns the absolute path within Silo file if writeToSilo()
    ///        or setSiloPath() was called before, the empty string otherwise.
    virtual std::string getSiloPath() const = 0;

    /// \brief Sets the silo path to be used when saving to a Silo file.
    virtual void setSiloPath(const std::string& path) = 0;

protected:
    /// \brief Destructor.
    virtual ~DomainChunk() {}
};

} // namespace weipa

#endif // __WEIPA_DOMAINCHUNK_H__

