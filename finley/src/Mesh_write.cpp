
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


/*****************************************************************************/

/*   Finley: write Mesh in finley file format */

/*****************************************************************************/

#define ESNEEDPYTHON
#include <esysUtils/first.h>

#include "Mesh.h"

#include <iomanip>

using std::cout;
using std::endl;
using std::ios;
using std::setw;
using std::string;

namespace finley {

// private
void Mesh::writeElementInfo(std::ostream& stream, const ElementFile* e,
                            const string defaultType) const
{
    if (e != NULL) {
        stream << e->referenceElementSet->referenceElement->Type->Name
          << " " << e->numElements << endl;
        const int NN = e->numNodes;
        for (index_t i=0; i < e->numElements; i++) {
            stream << e->Id[i] << " " << e->Tag[i];
            for (int j=0; j<NN; j++)
                stream << " " << Nodes->Id[e->Nodes[INDEX2(j,i,NN)]];
            stream << endl;
        }
    } else {
        stream << defaultType << " 0" << endl;
    }
}

// private
void Mesh::printElementInfo(const ElementFile* e, const string title,
                            const string defaultType, bool full) const
{
    if (e != NULL) {
        dim_t mine=0, overlap=0;
        for (index_t i=0; i < e->numElements; i++) {
            if (e->Owner[i] == MPIInfo->rank)
                mine++;
            else
                overlap++;
        }
        cout << "\t" << title << ": "
            << e->referenceElementSet->referenceElement->Type->Name
            << " " << e->numElements << " (TypeId="
            << e->referenceElementSet->referenceElement->Type->TypeId
            << ") owner=" << mine << " overlap=" << overlap << endl;
        if (full) {
            const int NN = e->numNodes;
            cout << "\t     Id   Tag Owner Color:  Nodes" << endl;
            for (index_t i=0; i < e->numElements; i++) {
                cout << "\t" << setw(7) << e->Id[i]
                     << setw(6) << e->Tag[i]
                     << setw(6) << e->Owner[i]
                     << setw(6) << e->Color[i] << ": ";
                for (int j=0; j<NN; j++)
                    cout << setw(6) << Nodes->Id[e->Nodes[INDEX2(j,i,NN)]];
                cout << endl;
            }
        }
    } else {
        cout << "\t" << title << ": " << defaultType << " 0" << endl;
    }
}

                             
/// writes the mesh to an external file using the 'fly' file format
void Mesh::write(const string filename) const
{
    if (MPIInfo->size >1) {
        throw escript::NotImplementedError("Mesh::write: only single rank runs are supported.");
    }

    std::ofstream f(filename.c_str());
    if (!f.is_open()) {
        std::stringstream ss;
        ss << "Mesh::write: Opening file " << filename << " for writing failred.";
        const string err(ss.str());
        throw escript::IOError(err);
    }

    // write header
    f << m_name << endl;

    // write nodes
    if (Nodes != NULL) {
        const int numDim = getDim();
        f << numDim << "D-Nodes " << Nodes->numNodes << endl;
        for (index_t i=0; i<Nodes->numNodes; i++) {
            f << Nodes->Id[i] << " " << Nodes->globalDegreesOfFreedom[i]
              << " " << Nodes->Tag[i];
            f.setf(ios::scientific, ios::floatfield);
            f.precision(15);
            for (int j=0; j<numDim; j++)
                f << " " << Nodes->Coordinates[INDEX2(j,i,numDim)];
            f << endl;
        }
    } else {
        f << "0D-Nodes 0" << endl;
    }

    // write elements
    writeElementInfo(f, Elements, "Tet4");

    // write face elements
    writeElementInfo(f, FaceElements, "Tri3");

    // write contact elements
    writeElementInfo(f, ContactElements, "Tri3_Contact");

    // write points
    writeElementInfo(f, Points, "Point1");

    // write tags
    if (tagMap.size() > 0) {
        f <<  "Tags" << endl;
        TagMap::const_iterator it;
        for (it=tagMap.begin(); it!=tagMap.end(); it++) {
            f << it->first << " " << it->second << endl;
        }
    }
    f.close();
#ifdef Finley_TRACE
    cout << "mesh " << m_name << " has been written to file " << filename << endl;
#endif
}

void Mesh::printInfo(bool full)
{
    cout << "PrintMesh_Info running on CPU " << MPIInfo->rank << " of "
              << MPIInfo->size << endl;
    cout << "\tMesh name '" << m_name << "'\n";
    cout << "\tApproximation order " << approximationOrder << endl;
    cout << "\tReduced Approximation order " <<reducedApproximationOrder << endl;
    cout << "\tIntegration order " << integrationOrder << endl;
    cout << "\tReduced Integration order " << reducedIntegrationOrder << endl;

    // write nodes
    if (Nodes != NULL) {
        const int numDim = getDim();
        cout << "\tNodes: " << numDim << "D-Nodes " << Nodes->numNodes << endl;
        if (full) {
            cout << "\t     Id   Tag  gDOF   gNI grDfI  grNI:  Coordinates\n";
            for (index_t i=0; i < Nodes->numNodes; i++) {
                cout << "\t" << setw(7) << Nodes->Id[i]
                     << setw(6) << Nodes->Tag[i]
                     << setw(6) << Nodes->globalDegreesOfFreedom[i]
                     << setw(6) << Nodes->globalNodesIndex[i]
                     << setw(6) << Nodes->globalReducedDOFIndex[i]
                     << setw(6) << Nodes->globalReducedNodesIndex[i] << ": ";
                cout.setf(ios::scientific, ios::floatfield);
                cout.precision(15);
                for (int j=0; j<numDim; j++)
                    cout << " " << Nodes->Coordinates[INDEX2(j,i,numDim)];
                cout << endl;
            }
        }
    } else {
        cout << "\tNodes: 0D-Nodes 0\n";
    }

    // write elements
    printElementInfo(Elements, "Elements", "Tet4", full);

    // write face elements
    printElementInfo(FaceElements, "Face elements", "Tri3", full);

    // write contact elements
    printElementInfo(ContactElements, "Contact elements", "Tri3_Contact", full);

    // write points
    printElementInfo(Points, "Points", "Point1", full);

    // write tags
    if (tagMap.size() > 0) {
        cout << "\tTags:\n";
        TagMap::const_iterator it;
        for (it=tagMap.begin(); it!=tagMap.end(); it++) {
            cout << "\t" << setw(7) << it->second << " " << it->first << endl;
        }
    }
}

} // namespace finley

