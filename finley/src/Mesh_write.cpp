
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

#include "FinleyDomain.h"

#include <escript/index.h>

#include <iomanip>
#include <iostream>

using std::cout;
using std::endl;
using std::ios;
using std::setw;
using std::string;

namespace finley {

// private
void FinleyDomain::writeElementInfo(std::ostream& stream, const ElementFile* e,
                                    const string& defaultType) const
{
    if (e != NULL) {
        stream << e->referenceElementSet->referenceElement->Type->Name
          << " " << e->numElements << endl;
        const int NN = e->numNodes;
        for (index_t i = 0; i < e->numElements; i++) {
            stream << e->Id[i] << " " << e->Tag[i];
            for (int j = 0; j < NN; j++)
                stream << " " << m_nodes->Id[e->Nodes[INDEX2(j,i,NN)]];
            stream << endl;
        }
    } else {
        stream << defaultType << " 0" << endl;
    }
}

// private
void FinleyDomain::printElementInfo(const ElementFile* e, const string& title,
                                    const string& defaultType, bool full) const
{
    if (e != NULL) {
        dim_t mine = 0, overlap = 0;
        for (index_t i = 0; i < e->numElements; i++) {
            if (e->Owner[i] == m_mpiInfo->rank)
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
            for (index_t i = 0; i < e->numElements; i++) {
                cout << "\t" << setw(7) << e->Id[i]
                     << setw(6) << e->Tag[i]
                     << setw(6) << e->Owner[i]
                     << setw(6) << e->Color[i] << ": ";
                for (int j = 0; j < NN; j++)
                    cout << setw(6) << m_nodes->Id[e->Nodes[INDEX2(j,i,NN)]];
                cout << endl;
            }
        }
    } else {
        cout << "\t" << title << ": " << defaultType << " 0" << endl;
    }
}


void FinleyDomain::write(const string& filename) const
{
    if (m_mpiInfo->size > 1)
        throw escript::NotImplementedError("FinleyDomain::write: only single rank "
                                           "runs are supported.");

    std::ofstream f(filename.c_str());
    if (!f.is_open()) {
        std::stringstream ss;
        ss << "FinleyDomain::write: Opening file " << filename << " for writing failed";
        throw escript::IOError(ss.str());
    }

    // write header
    f << m_name << endl;

    // write nodes
    if (m_nodes != NULL) {
        const int numDim = getDim();
        f << numDim << "D-Nodes " << m_nodes->getNumNodes() << endl;
        for (index_t i = 0; i < m_nodes->getNumNodes(); i++) {
            f << m_nodes->Id[i] << " " << m_nodes->globalDegreesOfFreedom[i]
              << " " << m_nodes->Tag[i];
            f.setf(ios::scientific, ios::floatfield);
            f.precision(15);
            for (int j = 0; j < numDim; j++)
                f << " " << m_nodes->Coordinates[INDEX2(j,i,numDim)];
            f << endl;
        }
    } else {
        f << "0D-Nodes 0" << endl;
    }

    // write elements
    writeElementInfo(f, m_elements, "Tet4");

    // write face elements
    writeElementInfo(f, m_faceElements, "Tri3");

    // write contact elements
    writeElementInfo(f, m_contactElements, "Tri3_Contact");

    // write points
    writeElementInfo(f, m_points, "Point1");

    // write tags
    if (m_tagMap.size() > 0) {
        f <<  "Tags" << endl;
        TagMap::const_iterator it;
        for (it = m_tagMap.begin(); it != m_tagMap.end(); it++) {
            f << it->first << " " << it->second << endl;
        }
    }
    f << endl;
    f.close();
#ifdef Finley_TRACE
    cout << "mesh " << m_name << " has been written to file " << filename << endl;
#endif
}

void FinleyDomain::Print_Mesh_Info(bool full) const
{
    cout << "PrintMeshInfo running on CPU " << m_mpiInfo->rank << " of "
              << m_mpiInfo->size << endl;
    cout << "\tMesh name '" << m_name << "'\n";
    cout << "\tApproximation order " << approximationOrder << endl;
    cout << "\tReduced Approximation order " <<reducedApproximationOrder << endl;
    cout << "\tIntegration order " << integrationOrder << endl;
    cout << "\tReduced Integration order " << reducedIntegrationOrder << endl;

    // write nodes
    if (m_nodes != NULL) {
        const int numDim = getDim();
        cout << "\tNodes: " << numDim << "D-Nodes " << m_nodes->getNumNodes() << endl;
        if (full) {
            cout << "\t     Id   Tag  gDOF   gNI grDfI  grNI:  Coordinates\n";
            for (index_t i = 0; i < m_nodes->getNumNodes(); i++) {
                cout << "\t" << setw(7) << m_nodes->Id[i]
                     << setw(6) << m_nodes->Tag[i]
                     << setw(6) << m_nodes->globalDegreesOfFreedom[i]
                     << setw(6) << m_nodes->globalNodesIndex[i]
                     << setw(6) << m_nodes->globalReducedDOFIndex[i]
                     << setw(6) << m_nodes->globalReducedNodesIndex[i] << ": ";
                cout.setf(ios::scientific, ios::floatfield);
                cout.precision(15);
                for (int j = 0; j < numDim; j++)
                    cout << " " << m_nodes->Coordinates[INDEX2(j,i,numDim)];
                cout << endl;
            }
        }
    } else {
        cout << "\tNodes: 0D-Nodes 0\n";
    }

    // write elements
    printElementInfo(m_elements, "Elements", "Tet4", full);

    // write face elements
    printElementInfo(m_faceElements, "Face elements", "Tri3", full);

    // write contact elements
    printElementInfo(m_contactElements, "Contact elements", "Tri3_Contact", full);

    // write points
    printElementInfo(m_points, "Points", "Point1", full);

    // write tags
    if (m_tagMap.size() > 0) {
        cout << "\tTags:\n";
        TagMap::const_iterator it;
        for (it = m_tagMap.begin(); it != m_tagMap.end(); it++) {
            cout << "\t" << setw(7) << it->second << " " << it->first << endl;
        }
    }
}

} // namespace finley

