
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
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

//
//
//
class MeshWithElements : public Mesh
{
public:
    /// Default constructor
    MeshWithElements();

    /// Copy constructor
    MeshWithElements(const MeshWithElements& m);
    
    /// Virtual destructor
    virtual ~MeshWithElements();

    virtual bool readFromNc(const std::string& filename);
    virtual void handleGhostZones(int ownIndex);
    virtual void removeGhostZones();
    virtual bool writeToSilo(DBfile* dbfile, const std::string& pathInSilo);
    StringVec getMeshNames() const;
    StringVec getVarNames() const;

    ElementData* getElementsByName(const std::string name) const;
    Mesh* getMeshByName(const std::string name) const;
    const IntVec& getVarDataByName(const std::string name) const;

    ElementData* getElements() { return cells; }
    ElementData* getFaceElements() { return faces; }
    ElementData* getContactElements() { return contacts; }
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

