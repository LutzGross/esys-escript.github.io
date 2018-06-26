
/*****************************************************************************
*
* Copyright (c) 2003-2018 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

#include "DudleyDomain.h"

#include <escript/index.h>

#include <iomanip>

using std::cout;
using std::endl;
using std::ios;
using std::setw;
using std::string;

namespace dudley {

// private
void DudleyDomain::writeElementInfo(std::ostream& stream, const ElementFile* e,
                                    const string& defaultType) const
{
    if (e != NULL) {
        stream << e->ename << " " << e->numElements << endl;
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
void DudleyDomain::printElementInfo(const ElementFile* e, const string& title,
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
            << e->ename << " " << e->numElements << " (TypeId=" << e->etype
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


void DudleyDomain::write(const std::string& filename) const
{
    if (m_mpiInfo->size > 1)
        throw escript::NotImplementedError("DudleyDomain::write: only single rank "
                                           "runs are supported.");

    std::ofstream f(filename.c_str());
    if (!f.is_open()) {
        std::stringstream ss;
        ss << "DudleyDomain::write: Opening file " << filename << " for writing failed";
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
    f.close();
#ifdef Dudley_TRACE
    cout << "mesh " << m_name << " has been written to file " << filename << endl;
#endif
}

void DudleyDomain::Print_Mesh_Info(bool full) const
{
    cout << "PrintMeshInfo running on CPU " << m_mpiInfo->rank << " of "
              << m_mpiInfo->size << endl;
    cout << "\tMesh name '" << m_name << "'\n";
    cout << "\tApproximation order " << 1 << endl;
    cout << "\tIntegration order " << 2 << endl;
    cout << "\tReduced Integration order " << 0 << endl;

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
                     << setw(6) << m_nodes->globalDegreesOfFreedom[i]
                     << setw(6) << m_nodes->globalNodesIndex[i] << ": ";
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

} // namespace dudley

