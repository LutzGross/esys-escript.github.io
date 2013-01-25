
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
// MeshWithElements.h
//
#ifndef __MESHWITHELEMENTS_H__
#define __MESHWITHELEMENTS_H__

#include <escriptreader/Mesh.h>

class DBfile;

namespace EscriptReader {

class ElementData;

/// \brief A full escript domain including elements.
///
/// This class represents a mesh including cells, faces, contact elements and
/// points if applicable.
class MeshWithElements : public Mesh
{
public:
    /// \brief Default constructor.
    MeshWithElements();

    /// \brief Copy constructor.
    MeshWithElements(const MeshWithElements& m);
    
    /// \brief Virtual destructor.
    virtual ~MeshWithElements();

    virtual bool readFromNc(const std::string& filename);
    virtual void handleGhostZones(int ownIndex);
    virtual void removeGhostZones();
    virtual bool writeToSilo(DBfile* dbfile, const std::string& pathInSilo);

    /// \brief
    StringVec getMeshNames() const;

    /// \brief
    StringVec getVarNames() const;

    /// \brief
    ElementData* getElementsByName(const std::string name) const;

    /// \brief
    Mesh* getMeshByName(const std::string name) const;

    /// \brief
    const IntVec& getVarDataByName(const std::string name) const;

    /// \brief Returns a pointer to the elements.
    ElementData* getElements() { return cells; }

    /// \brief Returns a pointer to the face elements.
    ElementData* getFaceElements() { return faces; }

    /// \brief Returns a pointer to the contact elements.
    ElementData* getContactElements() { return contacts; }

    /// \brief Returns a pointer to the point "elements".
    ElementData* getPoints() { return points; }

private:
    IntVec nodeTag, nodeGDOF, nodeGNI, nodeGRDFI, nodeGRNI;
    ElementData* cells;
    ElementData* faces;
    ElementData* contacts;
    ElementData* points;
};

} // namespace EscriptReader

#endif // __MESHWITHELEMENTS_H__

