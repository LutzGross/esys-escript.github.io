
/*****************************************************************************
*
* Copyright (c) 2003-2015 by The University of Queensland
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

namespace finley {

/// writes the mesh to the external file fname using the Finley file format
void Mesh::write(const std::string fname) const
{
    if (MPIInfo->size >1) {
        setError(IO_ERROR, "Mesh_write: only single processor runs are supported.");
        return;
    }
    // open file
    FILE* f=fopen(fname.c_str(), "w");
    if (f==NULL) {
        char error_msg[LenErrorMsg_MAX];
        sprintf(error_msg, "Mesh_write: Opening file %s for writing failed.",fname.c_str());
        setError(IO_ERROR,error_msg);
        return;
    }

    // write header
    fprintf(f, "%s\n", m_name.c_str());

    // write nodes:
    if (Nodes != NULL) {
        const int numDim = getDim();
        fprintf(f,"%1dD-Nodes %d\n", numDim, Nodes->numNodes);
        for (index_t i=0; i<Nodes->numNodes; i++) {
            fprintf(f,"%d %d %d",Nodes->Id[i],Nodes->globalDegreesOfFreedom[i],Nodes->Tag[i]);
            for (int j=0; j<numDim; j++)
                fprintf(f," %20.15e",Nodes->Coordinates[INDEX2(j,i,numDim)]);
            fprintf(f,"\n");
        }
    } else {
        fprintf(f,"0D-Nodes 0\n");
    }

    // write elements:
    if (Elements != NULL) {
        fprintf(f, "%s %d\n",Elements->referenceElementSet->referenceElement->Type->Name,Elements->numElements);
        const int NN=Elements->numNodes;
        for (index_t i=0; i<Elements->numElements; i++) {
            fprintf(f,"%d %d",Elements->Id[i],Elements->Tag[i]);
            for (int j=0; j<NN; j++)
                fprintf(f," %d",Nodes->Id[Elements->Nodes[INDEX2(j,i,NN)]]);
            fprintf(f,"\n");
        }
    } else {
        fprintf(f,"Tet4 0\n");
    }

    // write face elements:
    if (FaceElements != NULL) {
        fprintf(f, "%s %d\n", FaceElements->referenceElementSet->referenceElement->Type->Name,FaceElements->numElements);
        const int NN=FaceElements->numNodes;
        for (index_t i=0; i<FaceElements->numElements; i++) {
            fprintf(f,"%d %d",FaceElements->Id[i],FaceElements->Tag[i]);
            for (int j=0; j<NN; j++)
                fprintf(f," %d",Nodes->Id[FaceElements->Nodes[INDEX2(j,i,NN)]]);
            fprintf(f,"\n");
        }
    } else {
        fprintf(f,"Tri3 0\n");
    }

    // write contact elements:
    if (ContactElements != NULL) {
        fprintf(f, "%s %d\n",ContactElements->referenceElementSet->referenceElement->Type->Name,ContactElements->numElements);
        const int NN=ContactElements->numNodes;
        for (index_t i=0; i<ContactElements->numElements; i++) {
            fprintf(f,"%d %d",ContactElements->Id[i],ContactElements->Tag[i]);
            for (int j=0; j<NN; j++)
                fprintf(f," %d",Nodes->Id[ContactElements->Nodes[INDEX2(j,i,NN)]]);
            fprintf(f,"\n");
        }
    } else {
        fprintf(f,"Tri3_Contact 0\n");
    }

    // write points:
    if (Points != NULL) {
        fprintf(f, "%s %d\n",Points->referenceElementSet->referenceElement->Type->Name,Points->numElements);
        for (index_t i=0; i<Points->numElements; i++) {
            fprintf(f,"%d %d %d\n",Points->Id[i],Points->Tag[i],Nodes->Id[Points->Nodes[INDEX2(0,i,1)]]);
        }
    } else {
        fprintf(f,"Point1 0\n");
    }

    // write tags:
    if (tagMap.size() > 0) {
        fprintf(f, "Tags\n");
        TagMap::const_iterator it;
        for (it=tagMap.begin(); it!=tagMap.end(); it++) {
            fprintf(f, "%s %d\n", it->first.c_str(), it->second);
        }
    }
    fclose(f);
#ifdef Finley_TRACE
    cout << "mesh " << m_name << " has been written to file " << fname << endl;
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

    // write nodes:
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

    // write elements:
    if (Elements != NULL) {
        dim_t mine=0, overlap=0;
        for (index_t i=0; i<Elements->numElements; i++) {
            if (Elements->Owner[i] == MPIInfo->rank)
                mine++;
            else
                overlap++;
        }
        cout << "\tElements: "
            << Elements->referenceElementSet->referenceElement->Type->Name
            << " " << Elements->numElements << " (TypeId="
            << Elements->referenceElementSet->referenceElement->Type->TypeId
            << ") owner=" << mine << " overlap=" << overlap << endl;
        if (full) {
            const int NN=Elements->numNodes;
            cout << "\t     Id   Tag Owner Color:  Nodes\n";
            for (index_t i=0; i<Elements->numElements; i++) {
                cout << "\t" << setw(7) << Elements->Id[i]
                     << setw(6) << Elements->Tag[i]
                     << setw(6) << Elements->Owner[i]
                     << setw(6) << Elements->Color[i] << ": ";
                for (int j=0; j<NN; j++)
                    cout << setw(6) << Nodes->Id[Elements->Nodes[INDEX2(j,i,NN)]];
                cout << endl;
            }
        }
    } else {
        cout << "\tElements: Tet4 0\n";
    }

    // write face elements:
    if (FaceElements != NULL) {
        dim_t mine=0, overlap=0;
        for (index_t i=0; i < FaceElements->numElements; i++) {
            if (FaceElements->Owner[i] == MPIInfo->rank)
                mine++;
            else
                overlap++;
        }
        cout << "\tFace elements: "
            << FaceElements->referenceElementSet->referenceElement->Type->Name
            << " " << FaceElements->numElements << " (TypeId="
            << FaceElements->referenceElementSet->referenceElement->Type->TypeId
            << ") owner=" << mine << " overlap=" << overlap << endl;
        if (full) {
            const int NN=FaceElements->numNodes;
            cout << "\t     Id   Tag Owner Color:  Nodes\n";
            for (index_t i=0; i<FaceElements->numElements; i++) {
                cout << "\t" << setw(7) << FaceElements->Id[i]
                     << setw(6) << FaceElements->Tag[i]
                     << setw(6) << FaceElements->Owner[i]
                     << setw(6) << FaceElements->Color[i] << ": ";
                for (int j=0; j<NN; j++)
                    cout << setw(6) << Nodes->Id[FaceElements->Nodes[INDEX2(j,i,NN)]];
                cout << endl;
            }
        }
    } else {
        cout << "\tFace elements: Tri3 0\n";
    }

    // write Contact elements:
    if (ContactElements != NULL) {
        dim_t mine=0, overlap=0;
        for (index_t i=0; i<ContactElements->numElements; i++) {
            if (ContactElements->Owner[i] == MPIInfo->rank)
                mine++;
            else
                overlap++;
        }
        cout << "\tContact elements: "
            << ContactElements->referenceElementSet->referenceElement->Type->Name
            << " " << ContactElements->numElements << " (TypeId="
            << ContactElements->referenceElementSet->referenceElement->Type->TypeId
            << ") owner=" << mine << " overlap=" << overlap << endl;
        if (full) {
            const int NN=ContactElements->numNodes;
            cout << "\t     Id   Tag Owner Color:  Nodes\n";
            for (index_t i=0; i<ContactElements->numElements; i++) {
                cout << "\t" << setw(7) << ContactElements->Id[i]
                     << setw(6) << ContactElements->Tag[i]
                     << setw(6) << ContactElements->Owner[i]
                     << setw(6) << ContactElements->Color[i] << ": ";
                for (int j=0; j<NN; j++)
                    cout << setw(6) << Nodes->Id[ContactElements->Nodes[INDEX2(j,i,NN)]];
                cout << endl;
            }
        }
    } else {
        cout << "\tContact elements: Tri3_Contact 0\n";
    }

    // write points:
    if (Points != NULL) {
        dim_t mine=0, overlap=0;
        for (index_t i=0; i<Points->numElements; i++) {
            if (Points->Owner[i] == MPIInfo->rank)
                mine++;
            else
                overlap++;
        }
        cout << "\tPoints: "
            << Points->referenceElementSet->referenceElement->Type->Name
            << " " << Points->numElements << " (TypeId="
            << Points->referenceElementSet->referenceElement->Type->TypeId
            << ") owner=" << mine << " overlap=" << overlap << endl;
        if (full) {
            cout << "\t     Id   Tag Owner Color:  Nodes\n";
            for (index_t i=0; i<Points->numElements; i++) {
                cout << "\t" << setw(7) << Points->Id[i]
                     << setw(6) << Points->Tag[i]
                     << setw(6) << Points->Owner[i]
                     << setw(6) << Points->Color[i]
                     << setw(8) << Nodes->Id[Points->Nodes[INDEX2(0,i,1)]]
                     << endl;
            }
        }
    } else {
        cout << "\tPoints: Point1 0\n";
    }

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

