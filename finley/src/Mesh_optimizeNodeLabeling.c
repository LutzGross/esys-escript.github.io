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

/*   Finley: Mesh : optimizes the labeling of nodes */

/**************************************************************/

/*   Author: gross@access.edu.au */
/*   Version: $Id:$ */

/**************************************************************/

#include "Mesh.h"

/**************************************************************/

void Finley_Mesh_optimizeNodeLabeling(Finley_Mesh* mesh_p) {

   index_t *DOF_permutation=NULL, i;
/*
   Paso_SystemMatrixPattern* pattern_p=NULL;
   pattern_p=Finley_getPattern(mesh_p,FALSE,FALSE);
   if no error {
       XXX=pattern_p. ;
       DOF_permutation=MEMALLOC(XXX, index_t);
       availbale=MEMALLOC(XXX, index_t);
       #pragma omp for private(i)
       for (i=0;i< ;++i) {
           DOF_permutation[i]=i;
       }
       first_available_node=0;
       while (first_available_node <= XXX) {
            // get an available node with minimum number of naighbours 
            max_level_size=XXX;
            root=...
            // get the leveling string from root 
            num_levels_tmp=1;
            num_nodes_in_level_tmp[0]=0;
            levels_tmp[];
            // get maximum level size 
            max_level_size_tmp=MAX(num_nodes_in_level_tmp[i+1]-num_nodes_in_level_tmp[i], max_level_size_tmp);
            // use new leveling if 
            if (max_level_size_tmp<max_level_size) {
                max_level_size=max_level_size_tmp;
                num_levels=num_levels_tmp;
                num_nodes_in_level=num_nodes_in_level_tmp;
            }
            
            
       
       }
   }
*/
}
