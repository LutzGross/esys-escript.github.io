
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


/****************************************************************************/

/*   Paso: pattern: a simple algorithm to create a new labeling */
/*   of input/output with a minimum bandwidth. The algorithm    */
/*   starts from a vertex with a minimum number of neighbours   */
/*   (connections in the pattern). From this root it searches   */
/*   recursively the neighbours. In the last level a new root   */
/*   is picked (again the vertex with the minimum number of     */
/*   neighbours) and the process is repeated. The algorithm is  */
/*   repeated until the maximum number of vertices in a level   */
/*   cannot be reduced.                                         */

/****************************************************************************/

/*   Author: Lutz Gross, l.gross@uq.edu.au */

/****************************************************************************/

#include "Pattern.h"

namespace paso {

// calculate initial bandwidth for a given labeling
dim_t Pattern::getBandwidth(index_t* label) const
{
    dim_t bandwidth = 0;

#pragma omp parallel
    {
        dim_t local_bandwidth=0;
#pragma omp for
        for (dim_t i=0; i < numOutput; ++i) {
            const index_t k = label[i];
            for (index_t iptr=ptr[i]; iptr < ptr[i+1]; ++iptr) {
                local_bandwidth = std::max(local_bandwidth,std::abs(k-label[index[iptr]]));
            }
        }
#pragma omp critical
        {
            bandwidth=std::max(local_bandwidth, bandwidth);
        }
    }
    return bandwidth;
}

// structure used to create an ordering by increasing degree
struct DegreeAndIdx
{
    dim_t deg;
    index_t idx;
};

int comparDegreeAndIdx(const void *arg1, const void *arg2)
{
    DegreeAndIdx a1=*(DegreeAndIdx*)arg1;
    DegreeAndIdx a2=*(DegreeAndIdx*)arg2;
    if (a1.deg < a2.deg) {
        return -1;
    } else if (a1.deg > a2.deg) {
        return 1;
    } else if (a1.idx < a2.idx) {
        return -1;
    } else if (a1.idx > a2.idx) {
        return 1;
    }
    return 0;
}

/*  dropTree() drops a tree in pattern from root
 *  root                 - on input the starting point of the tree.
 *  AssignedLevel        - array of length numOutput indicating the level
 *                         assigned to vertex
 *  VerticesInTree       - on output contains the vertices used in tree
 *                         (array of length numOutput)
 *  numLevels            - on output the number of levels used
 *  firstVertexInLevel   - on output firstVertexInLevel[i] points to the first
 *                         vertex of level i in VerticesInTree
 *                         firstVertexInLevel[i+1]-firstVertexInLevel[i] is the
 *                         number of vertices in level i (array of length
 *                         numOutput+1)
 *  max_LevelWidth_abort - input param which triggers early return if
 *                         LevelWidth becomes >= max_LevelWidth_abort
 */
bool dropTree(index_t root, const Pattern* pattern, index_t* AssignedLevel,
              index_t* VerticesInTree, dim_t* numLevels,
              index_t* firstVertexInLevel, dim_t max_LevelWidth_abort,
              dim_t N)
{
#pragma omp parallel for
    for (dim_t i=0; i < pattern->numInput; ++i)
        AssignedLevel[i]=-1;

    dim_t nlvls = 0;
    AssignedLevel[root] = 0;
    VerticesInTree[0] = root;
    firstVertexInLevel[0]=0;
    dim_t level_top = firstVertexInLevel[0]+1;

    while (firstVertexInLevel[nlvls] < level_top) {
        nlvls++;
        firstVertexInLevel[nlvls] = level_top;
        if (firstVertexInLevel[nlvls]-firstVertexInLevel[nlvls-1] >= max_LevelWidth_abort)
            return false;

        for (dim_t i=firstVertexInLevel[nlvls-1]; i < firstVertexInLevel[nlvls]; ++i) {
            const index_t k = VerticesInTree[i];
            for (index_t j=pattern->ptr[k]; j<pattern->ptr[k+1]; ++j) {
                const index_t itest = pattern->index[j];
                if (AssignedLevel[itest] < 0) {
#ifdef BOUNDS_CHECK
                    if (itest < 0 || itest >= N) {
                        printf("BOUNDS_CHECK %s %d itest=%d\n",
                                __FILE__, __LINE__, itest);
                        exit(1);
                    }
                    if (level_top < 0 || level_top >= N) {
                        printf("BOUNDS_CHECK %s %d level_top=%d\n",
                                __FILE__, __LINE__, level_top);
                        exit(1);
                    }
#endif
                    AssignedLevel[itest] = nlvls;
                    VerticesInTree[level_top] = itest;
                    level_top++;
                }
            }
        }
    }
    *numLevels=nlvls;
    return true;
}

/* the driver */
void Pattern::reduceBandwidth(index_t* oldToNew)
{
    if (numOutput != numInput) {
        Esys_setError(VALUE_ERROR, "Pattern::reduceBandwidth: pattern needs to be for a square matrix.");
        return;
    } else if (numOutput == 0) {
        return;
    }

    dim_t i, N=numOutput, numLevels, deg, bandwidth;
    index_t k;

    /* printf("relabeling of %d DOFs started.\n",N); */
    DegreeAndIdx* degAndIdx = new DegreeAndIdx[N];
    index_t* oldLabel = new index_t[N];
    index_t* AssignedLevel = new index_t[N];
    index_t* VerticesInTree = new index_t[N];
    index_t* firstVertexInLevel = new index_t[N+1];

    // get the initial bandwidth
#pragma omp parallel for private(i)
    for (i=0; i<N; ++i) oldToNew[i]=i;

    dim_t initial_bandwidth = getBandwidth(oldToNew);
    // printf("initial bandwidth = %d\n",initial_bandwidth);

    #pragma omp parallel for private(i)
    for (i=0; i<N; ++i) {
        oldToNew[i] = -1;
        degAndIdx[i].idx = i;
        degAndIdx[i].deg = ptr[i+1] - ptr[i];
    }

    // create an ordering with increasing degree
    qsort(degAndIdx, (size_t)N, sizeof(DegreeAndIdx), comparDegreeAndIdx);
    index_t root = degAndIdx[0].idx;
    dim_t numLabledVertices = 0;

    while (root >= 0) {
        dim_t max_LevelWidth = N+1;
        dim_t numVerticesInTree = 0;
        while (dropTree(root, this, AssignedLevel, VerticesInTree,
                        &numLevels, firstVertexInLevel, max_LevelWidth, N)) {
            // find new maximum level width
            max_LevelWidth=0;
            for (i=0; i<numLevels; ++i) {
#ifdef BOUNDS_CHECK
                if (i >= N+1) {
                    printf("BOUNDS_CHECK %s %d i=%d N=%d\n", __FILE__,
                            __LINE__, i, N); exit(1);
                }
#endif
                max_LevelWidth=std::max(max_LevelWidth, firstVertexInLevel[i+1]-firstVertexInLevel[i]);
            }
            // find a vertex in the last level which has minimum degree
            dim_t min_deg = N+1;
            root=-1;
            for (i = firstVertexInLevel[numLevels-1]; i<firstVertexInLevel[numLevels];++i) {
                k=VerticesInTree[i];
                deg=ptr[k+1]-ptr[k];
                if (deg<min_deg) {
                    min_deg=deg;
                    root=k;
                }
            }
            // save the vertices in the current tree
            numVerticesInTree=firstVertexInLevel[numLevels];
            for (i=0;i<firstVertexInLevel[numLevels];++i) {
#ifdef BOUNDS_CHECK
                if (numLabledVertices+i < 0 || numLabledVertices+i >= N) { printf("BOUNDS_CHECK %s %d i=%d numLabeledVertices=%d root=%d N=%d firstVertexInLevel[numLevels]=%d\n", __FILE__, __LINE__, i, numLabledVertices, root, N, firstVertexInLevel[numLevels]); exit(1); }
#endif
                oldLabel[numLabledVertices+i]=VerticesInTree[i];
            }
        }
        // now the vertices in the current tree
        for (i=0; i<numVerticesInTree; ++i) {
#ifdef BOUNDS_CHECK
            if (numLabledVertices+i < 0 || numLabledVertices+i >= N) { printf("BOUNDS_CHECK %s %d i=%d numLabeledVertices=%d root=%d N=%d\n", __FILE__, __LINE__, i, numLabledVertices, root, N); exit(1); }
#endif
            oldToNew[oldLabel[numLabledVertices+i]]=numLabledVertices+i;
        }
        numLabledVertices+=numVerticesInTree;
        // new search for a vertex which is not labeled yet
        root=-1;
        for (i=0; i<N; ++i) {
            if (oldToNew[degAndIdx[i].idx] < 0) {
                root=degAndIdx[i].idx;
                break;
            }
        }
    } // end of while root loop
    bandwidth=getBandwidth(oldToNew);
    // printf("bandwidth after DOF relabeling= %d\n",bandwidth);
    if (bandwidth >= initial_bandwidth) {
        //printf("initial labeling used.\n");
#pragma omp parallel for private(i)
        for (i=0; i<N; ++i) oldToNew[i]=i;
    }
    delete[] degAndIdx;
    delete[] oldLabel;
    delete[] AssignedLevel;
    delete[] VerticesInTree;
    delete[] firstVertexInLevel;
}

} // namespace paso

