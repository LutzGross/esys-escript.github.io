
/*****************************************************************************
*
* Copyright (c) 2003-2013 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/


/************************************************************************************/

/*   Paso: pattern: a simple algorithm to create a new labeling */
/*   of input/output with a minimum bandwidth. The algorithm    */
/*   starts from a vertex with a minimum number of neighbours   */
/*   (connections in the pattern). From this root it searches   */
/*   recursively the neighbours. In the last level a new root   */
/*   is picked (again the vertex with the minimum number of     */
/*   neighbours) and the process is repeated. The algorithm is  */
/*   repeated until the maximum number of vertices in a level   */
/*   cannot be reduced.                                         */

/************************************************************************************/

/*   Author: Lutz Gross, l.gross@uq.edu.au */

/************************************************************************************/


#include "Paso.h"
#include "Pattern.h"
#include "Common.h"

/*   calculate initial bandwidth for a given labeling */
dim_t Paso_Pattern_getBandwidth(Paso_Pattern* pattern, index_t* label) {
      register index_t k;
      index_t iptr;
      dim_t bandwidth = 0, local_bandwidth=0,i ;

      #pragma omp parallel private(local_bandwidth)
      {
         local_bandwidth=0;
         #pragma omp for private(i,iptr,k)
         for (i=0;i<pattern->numOutput;++i) {
              k=label[i];
              for (iptr=pattern->ptr[i];iptr<pattern->ptr[i+1];++iptr) {
                   local_bandwidth = MAX(local_bandwidth, ABS(k - label[pattern->index[iptr]]));
              }
          }
          #pragma omp critical
          {
              bandwidth=MAX(local_bandwidth,bandwidth);
          }
     }
     return bandwidth;
}
/*    structures used to create an ordering by increasing degree */
typedef struct Paso_DegreeAndIdx {
   dim_t deg;
   index_t idx;
} Paso_DegreeAndIdx;

int Paso_comparDegreeAndIdx(const void *arg1,const void *arg2){
   Paso_DegreeAndIdx a1=*(Paso_DegreeAndIdx*)arg1;
   Paso_DegreeAndIdx a2=*(Paso_DegreeAndIdx*)arg2;
   if (a1.deg<a2.deg) {
      return -1;
   } else if (a1.deg>a2.deg) {
      return 1;
   } else if (a1.idx<a2.idx) {
       return -1;
   } else if (a1.idx>a2.idx) {
      return 1;
   } else {
      return 0;
   }
}

/*  Paso_Pattern_dropTree drops a tree in pattern from root */

/*  root - on input the starting point of the tree. */
/*  AssignedLevel- array of length pattern->numOutput indicating the level assigned to vertex */
/*  VerticesInTree - on output contains the vertices used in tree  (array of length pattern->numOutput) */
/*  numLevels-  on output the number of levels used. */
/*  firstVertexInLevel - on output firstVertexInLevel[i] points to the first vertex of level i in VerticesInTree. */
/*                       firstVertexInLevel[i+1]-firstVertexInLevel[i] is the number of vertices in level i. */
/*                         (array of length pattern->numOutput+1) */
/*  max_LevelWidth_abort-  input param which triggers early return if LevelWidth becomes >= max_LevelWidth_abort */
bool Paso_Pattern_dropTree(index_t root, 
                             Paso_Pattern *pattern, 
                             index_t *AssignedLevel, 
                             index_t *VerticesInTree,
                             dim_t* numLevels, 
                             index_t* firstVertexInLevel, 
                             dim_t max_LevelWidth_abort,
			     dim_t N) 
{
  dim_t nlvls,i;
  index_t level_top,k ,itest,j;

  #pragma omp parallel for private(i)
  for (i=0;i<pattern->numInput;++i) AssignedLevel[i]=-1;

  nlvls=0;
  AssignedLevel[root] = 0;
  VerticesInTree[0] = root;
  firstVertexInLevel[0]=0;
  level_top=firstVertexInLevel[0]+1;
  while (firstVertexInLevel[nlvls]<level_top) {
      nlvls++;
      firstVertexInLevel[nlvls]=level_top;
      if (firstVertexInLevel[nlvls]-firstVertexInLevel[nlvls-1]>=max_LevelWidth_abort) return FALSE;
      for (i=firstVertexInLevel[nlvls-1];i<firstVertexInLevel[nlvls];++i) {
         k = VerticesInTree[i];
         for (j=pattern->ptr[k];j<pattern->ptr[k+1];++j) {
           itest = pattern->index[j];
           if (AssignedLevel[itest]<0) {
#ifdef BOUNDS_CHECK
		       if (itest < 0 || itest >= N) { printf("BOUNDS_CHECK %s %d itest=%d\n", __FILE__, __LINE__, itest); exit(1); }
		       if (level_top < 0 || level_top >= N) { printf("BOUNDS_CHECK %s %d level_top=%d\n", __FILE__, __LINE__, level_top); exit(1); }
#endif
              AssignedLevel[itest] = nlvls;
              VerticesInTree[level_top] = itest;
              level_top++;
           }
         }
      }
   }
   *numLevels=nlvls;
   return TRUE;
}

/* the driver */
void Paso_Pattern_reduceBandwidth(Paso_Pattern* pattern,index_t* oldToNew) {
   dim_t i, initial_bandwidth, N=pattern->numOutput, numLabledVertices, numLevels, max_LevelWidth, min_deg,deg, numVerticesInTree=0, bandwidth;
   Paso_DegreeAndIdx* degAndIdx=NULL;
   index_t root, *AssignedLevel=NULL, *VerticesInTree=NULL, *firstVertexInLevel=NULL,k, *oldLabel=NULL;
   /* check input */
   if (N != pattern->numInput) {
      Esys_setError(VALUE_ERROR,"Paso_Pattern_reduceBandwidth: pattern needs to be for a square matrix.");
   } else if (N > 0) {
/* printf("relabeling of %d DOFs started.\n",N); */
      degAndIdx=new Paso_DegreeAndIdx[N];
      oldLabel=new index_t[N];
      AssignedLevel=new index_t[N];
      VerticesInTree=new index_t[N];
      firstVertexInLevel=new index_t[N+1];
      if (! ( Esys_checkPtr(degAndIdx) || Esys_checkPtr(oldLabel) || Esys_checkPtr(AssignedLevel) || Esys_checkPtr(VerticesInTree) || Esys_checkPtr(firstVertexInLevel) ) ) {
         /* get the initial bandwidth */
         #pragma omp parallel for private(i)
         for (i=0;i<N;++i) oldToNew[i]=i; 
         initial_bandwidth=Paso_Pattern_getBandwidth(pattern,oldToNew);
/* printf("initial bandwidth = %d\n",initial_bandwidth); */
         /* get the initial bandwidth */
         #pragma omp parallel for private(i)
         for (i=0;i<N;++i) {
            oldToNew[i]=-1; 
            degAndIdx[i].idx=i; 
            degAndIdx[i].deg=pattern->ptr[i+1]-pattern->ptr[i];
         }
         /* create an ordering with increasing degree */
         #ifdef USE_QSORTG
            qsortG(degAndIdx,(size_t)N,sizeof(Paso_DegreeAndIdx),Paso_comparDegreeAndIdx);
         #else
             qsort(degAndIdx,(size_t)N,sizeof(Paso_DegreeAndIdx),Paso_comparDegreeAndIdx);
         #endif

         root=degAndIdx[0].idx;
         numLabledVertices=0;

         while (root>=0) {
             max_LevelWidth=N+1;
                 
             while (Paso_Pattern_dropTree(root,pattern,AssignedLevel,VerticesInTree,
                                          &numLevels,firstVertexInLevel,max_LevelWidth,N) ) {

                 /* find new maximum level width */
                 max_LevelWidth=0;
                 for (i=0;i<numLevels;++i) {
#ifdef BOUNDS_CHECK
		      if (i < 0 || i >= N+1) { printf("BOUNDS_CHECK %s %d i=%d N=%d\n", __FILE__, __LINE__, i, N); exit(1); }
#endif
                      max_LevelWidth=MAX(max_LevelWidth,firstVertexInLevel[i+1]-firstVertexInLevel[i]);
                 }
                 /* find a vertex in the last level which has minimum degree */
                 min_deg=N+1;
                 root=-1;
                 for (i = firstVertexInLevel[numLevels-1]; i<firstVertexInLevel[numLevels];++i) {
                      k=VerticesInTree[i];
                      deg=pattern->ptr[k+1]-pattern->ptr[k];
                      if (deg<min_deg) {
                         min_deg=deg;
                         root=k;
                      } 
                 }
                 /* save the vertices in the current tree */
                 numVerticesInTree=firstVertexInLevel[numLevels];
                 for (i=0;i<firstVertexInLevel[numLevels];++i) {
#ifdef BOUNDS_CHECK
		       if (numLabledVertices+i < 0 || numLabledVertices+i >= N) { printf("BOUNDS_CHECK %s %d i=%d numLabeledVertices=%d root=%d N=%d firstVertexInLevel[numLevels]=%d\n", __FILE__, __LINE__, i, numLabledVertices, root, N, firstVertexInLevel[numLevels]); exit(1); }
#endif
                       oldLabel[numLabledVertices+i]=VerticesInTree[i];
                 }
             }
             /* now the vertices in the current tree */
             for (i=0;i<numVerticesInTree;++i) {
#ifdef BOUNDS_CHECK
		 if (numLabledVertices+i < 0 || numLabledVertices+i >= N) { printf("BOUNDS_CHECK %s %d i=%d numLabeledVertices=%d root=%d N=%d\n", __FILE__, __LINE__, i, numLabledVertices, root, N); exit(1); }
#endif
                 oldToNew[oldLabel[numLabledVertices+i]]=numLabledVertices+i;
             }
             numLabledVertices+=numVerticesInTree;
             /* new search for a vertex which is not labeled yet */
             root=-1;
             for (i=0;i<N;++i) {
                 if (oldToNew[degAndIdx[i].idx] < 0) {
                     root=degAndIdx[i].idx;
                     break;
                 }
             }
        } /* end of while root loop */
        bandwidth=Paso_Pattern_getBandwidth(pattern,oldToNew);
/* printf("bandwidth after DOF relabeling= %d\n",bandwidth); */
        if (bandwidth>=initial_bandwidth) {
/* printf("initial labeling used.\n"); */
           #pragma omp parallel for private(i)
           for (i=0;i<N;++i) oldToNew[i]=i; 
        }
      }      
      delete[] degAndIdx;
      delete[] oldLabel;
      delete[] AssignedLevel;
      delete[] VerticesInTree;
      delete[] firstVertexInLevel;
   }
}

