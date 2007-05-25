/*
 ************************************************************
 *          Copyright 2006, 2007 by ACcESS MNRF             *
 *                                                          *
 *              http://www.access.edu.au                    *
 *       Primary Business: Queensland, Australia            *
 *  Licensed under the Open Software License version 3.0    *
 *     http://www.opensource.org/licenses/osl-3.0.php       *
 *                                                          *
 ************************************************************
*/

/**************************************************************/

/*   Finley: Mesh : optimizes the labeling of DOFs */

/**************************************************************/

/*   Author: gross@access.edu.au */
/*   Version: $Id:$ */

/**************************************************************/

#include "Mesh.h"

/**************************************************************/

void Finley_Mesh_optimizeNodeLabeling(Finley_Mesh* mesh_p) {

#if 0


   index_t i, n, initial_bandwidth, first_available_DOF, first_available_level, max_level_size, num_levels;
   bool_t next_iteration,
   Paso_SystemMatrixPattern* pattern_p=NULL;

   /* get access to the matrix pattern */
   pattern_p=Finley_getPattern(mesh_p,FALSE,FALSE);
   if (Finley_noError()) {
   if no error {
       n=pattern_p->n_ptr;
old_to_new=MEMALLOC(n, index_t);
availbale_DOF=MEMALLOC(n, bool_t);
       #pragma omp for private(i)
       for (i=0;i<n;++i) {
           old_to_new[i]=i;
           availbale_DOF=TRUE;
       }
       /* get initial bandwidth */
       initial_bandwidth=Finley_Mesh_getDegree(pattern_p,old_to_new);
printf("initial_bandwidth = %d\n",initial_bandwidth);
       /* make sure that all connection components are processed */
       first_available_DOF=0;
       first_available_level=0;
       while (root=Finley_Mesh_FindMinDegreeNode(pattern_p,availbale_DOF, TRUE) >-1) {
            // get an available DOF with minimum number of naighbours 
            max_level_size=n;
            next_iteration=TRUE;
            num_levels=n;
            #pragma omp parallel for private(i)
            for (i=0;i<n;++i) {
                if (availbale_DOF[i]) 
                    level_mask[i]=num_levels-1;
                } else {
                    level_mask[i]=-1;
                }
            }
            /* run through graph diameter and creates levels until the maximum level size cannot be improved anymore */
            while (root>-1) {
               /*   num_levels_new is the number of levels spanned from root.
                    The DOFs ptr_DOFs_in_level_new[i] to ptr_DOFs_in_level_new[i+1]-1 in DOFs_in_level_new are 
                    maked for level i = 0,...,num_levels_new-1 */
               Finley_Mesh_MarkLevels(root,pattern_p,level_mask,num_levels-1,
                                      &num_levels_new, ptr_DOFs_in_level_new, DOFs_in_level_new);
               /* now we calculate the maximum level size */
               max_level_size_new=0;
               for (i=0;i<num_levels_new;++i) {
                   max_level_size_new=MAX(ptr_DOFs_in_level_new[i+1]-ptr_DOFs_in_level_new[i], max_level_size_new);
               }
               /* if there is a reduction in the maximum level size the we accept this leveling  */
               /* otherwise we give up                                                           */
               if (max_level_size_new < max_level_size) {
                  max_level_size=max_level_size_new;
                  num_levels=num_levels_new;
                  for (i=0;i<num_levels_new+1;++i) {
                      ptr_DOFs_in_level[i+first_available_level]=ptr_DOFs_in_level_new[i];
                  }
                  for (i=0;i<ptr_DOFs_in_level_new[num_levels_new];++i) {
                      DOFs_in_level[i+first_available_DOF]=DOFs_in_level_new[i]+first_available_level;
                  }
                  root=Finley_Mesh_FindMinDegreeNodeFromList(pattern_p,
                                                             ptr_DOFs_in_level_new[num_levels_new]-ptr_DOFs_in_level_new[num_levels_new-1],
                                                             DOFs_in_level_new[ptr_DOFs_in_level_new[num_levels_new-1]]);
                                                             
               } else { 
                   root=-1;
               }
            }
            #pragma omp parallel for private(i) 
            for (i=first_available_DOF;i<ptr_DOFs_in_level_new[num_levels];++i) {
                 availbale_DOF[DOFs_in_level[i]]=FALSE;
            }
            first_available_DOF+=ptr_DOFs_in_level_new[num_levels];
            first_available_level+=num_levels;
       }
       /* create new_to_old labeling (we don't do anything here) */
       #pragma omp parallel for private(i)
       for (i=0;i<n;++i) new_to_old[i]=DOFs_in_level[i];
       /* invert the new_to_old labeling */
       #pragma omp parallel for private(i)
       for (i=0;i<n;++i) old_to_new[new_to_old[i]]=i;
       /* now we can start to assign new labels to DOFs */
       new_bandwidth=Finley_Mesh_getDegree(pattern_p,old_to_new);
       if (new_bandwidth < initial_bandwidth) {


       }
#endif
}
#if 0
dim_t Finley_Mesh_FindMinDegreeNode(Paso_SystemMatrixPattern* pattern_p,index_t* available,index_t indicator) {
      index_t min_deg_local, min_deg_node_local, min_deg, min_deg_node;
      register index_t deg;
      min_deg=pattern_p->n_ptr;
      min_deg_node=-1;
      #pragma omp parallel private(min_deg_local, min_deg_node_local)
      {
        min_deg_local=pattern_p->n_ptr;
        min_deg_node=-1;
        #pragma omp for private(i,iptr,deg)
        for (i=0,i<pattern->n_ptr,++i) {
           if ( available[i] == indicator) {
              deg=0;
              for (iptr = pattern->ptr[i+1]; iptr<pattern->ptr[i+1]; ++iptr) {
                   if (available[pattern->index[iptr]]==indicator) ++deg;
              }
              if (deg < min_deg_local) {
                  min_deg_local=deg;
                  min_deg_node_local=i;
              }
           }
         }
         #pragma omp critical
         {
            if ((min_deg_local<min_deg) && (min_deg_local<min_deg)) {
                min_deg=min_deg_local;
                min_deg_node=min_deg_node_local;
            }
         }
      }
      return min_deg_node;
}

index_t Finley_Mesh_getDegree(Paso_SystemMatrixPattern* pattern_p, index_t *label) {
      index_t bandwidth, i, iptr, local_bandwidth;
      bandwidth = 0;
      #pragma omp parallel private(local_bandwidth)
      {
        local_bandwidth=0;
        #pragma omp for private(i,i_ptr)
        for (i=0,i<pattern->n_ptr,++i) {
           for (iptr=pattern->ptr[i],iptr<pattern->ptr[i+1],++iptr) {
               local_bandwidth = MAX(local_bandwidth,ABS(label[i] - label(pattern->index[iptr])));
           }
         }
         #pragma omp critical
         {
            bandwidth=MAX(bandwidth,local_bandwidth);
         }
      }
      return bandwidth;
}

void Finley_Mesh_MarkLevels(index_t root,
                            Paso_SystemMatrixPattern* pattern_p,
                            dim_t *num_levels, dim_t *ptr_DOFs_in_level, index_t *DOFs_in_level)

{
  DOFs_in_level[0]=root;
  available[root]=FALSE;
  ptr_DOFs_in_level[0]=0;
  level_count=1;
  nn=1;
  while (ptr_DOFs_in_level[level_count-1] < nn ) {
      ptr_DOFs_in_level[level_count] = nn;
      ++level_count;
      for (i = ptr_DOFs_in_level[level_count-1]; i < ptr_DOFs_in_level[level_count]; ++i) {
           dof=DOFs_in_level[i]
           for (iptr=pattern->ptr[dof],iptr<pattern->ptr[dof+1],++iptr) {
              itmp=pattern->index[iptr];
              if (! available[itmp]) {
                   available[itmp]=FALSE;
                   DOFs_in_level[nn]=itmp;
                   ++nn;
              }



      }
  }
  *num_levels=level_count;
  ptr_DOFs_in_level[level_count]=  ;
}
#endif
