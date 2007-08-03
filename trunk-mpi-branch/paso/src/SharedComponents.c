/* $Id:$ */

/*
********************************************************************************
*               Copyright   2006, 2007 by ACcESS MNRF                          *
*                                                                              * 
*                 http://www.access.edu.au                                     *
*           Primary Business: Queensland, Australia                            *
*     Licensed under the Open Software License version 3.0 		       *
*        http://www.opensource.org/licenses/osl-3.0.php                        *
********************************************************************************
*/

/**************************************************************/

/* Paso: SharedComponents organizes the coupling with in a pattern/matrix   */
/*       across processors                                         */

/**************************************************************/
 
/* Author: gross@access.edu.au */

/**************************************************************/

#include "SharedComponents.h"

/**************************************************************/

/* allocates a SharedComponents  */


/**************************************************************/

typedef struct Paso_SharedComponents_Container {
  index_t i;
  index_t n;
  Paso_MPI_rank p;
  struct Finley_IndexList *extension;
} Paso_SharedComponents_Container;

int Paso_SharedComponents_compar(const void *arg1 , const void *arg2 ) {
   Paso_SharedComponents_Container* e1=(Paso_SharedComponents_Container*) arg1;
   Paso_SharedComponents_Container* e2=(Paso_SharedComponents_Container*) arg2;
   if ( e1->p < e2->p ) {
       return -1;
   } else if ( e1->p > e2->p ) {
       return 1;
   } else {
       if ( e1->i < e2->i ) {
           return -1;
       } else if ( e1->i > e2->i ) {
           return 1;
       } else {
           return 0;
       }
   }
}

Paso_SharedComponents* Paso_SharedComponents_alloc(
                                                   dim_t numComponents,
                                                   index_t* globalComponent, dim_t m, index_t b,
                                                   Paso_Distribution *distribution)
{
  Paso_SharedComponents_Container* comps=NULL;
  Paso_SharedComponents*out=NULL;
  dim_t numComponents2=numComponents*m;
  Paso_MPI_rank p_min, p_max;
  dim_t i,j;
  register index_t itmp, offset;
  Paso_resetError();
  #ifdef XXXXX
  comps=MEMALLOC(MAX(numComponents2,Paso_SharedComponents_Container);
  if (!Paso_checkPtr(comps)) {
         #pragma omp parallel private(offset,i,j,itmp, local_minglobalComponent, local_maxglobalComponent, p)
         {
            if (globalComponent == NULL) {
               offset=Paso_Distribution_getFirstComponent(distribution)+b;
               #pragma omp for 
               for (i=0;i<numComponents;++i){
                   itmp=m*i+offset;
                   for (j=0;j<m;++j) comps[m*i+j].n=itmp+j;
               }
            } else {
               #pragma omp for
               for (i=0;i<numComponents;++i){
                   itmp=m*globalComponent[i]+b;
                   for (j=0;j<m;++j) comps[m*i+j].n=itmp+j;
               }
            }
            /* find p_min and p_max to minize the serach for the processor of a component */
            #pragma omp master
            {
               minglobalComponent=Paso_Distribution_getMaxGlobalComponents(distribution);
               maxglobalComponent=Paso_Distribution_getMinGlobalComponents(distribution);
            }
            local_minglobalComponent=minglobalComponent;
            local_maxglobalComponent=maxglobalComponent;
            #pragma omp for
            for (i=0;i<numComponents2;++i) {
                 local_minglobalComponent=MIN(local_minglobalComponent,comps[i].n);
                 local_maxglobalComponent=MAX(local_maxglobalComponent,comps[i].n);
                 comps[i].i=i;
                 comps[i].p=-1;
            }
            #pragma omp critical
            {
               minglobalComponent=MIN(minglobalComponent,local_minglobalComponent);
               maxglobalComponent=MAX(maxglobalComponent,local_maxglobalComponent);
            }
            #pragma omp master
            {
                p_min=-1;
                p_max=-1;
                for (p=0; p<distribution->mpi_info->size, ++p) {
                   if (distribution->firstComponent[p]<=minglobalComponent) p_min=p;
                   if (maxglobalComponent<distribution->firstComponent[p+1]) p_max=p;
                }
                fprintf("Lutz: detected processor range %d %d\n",p_min,p_max);
             }
             #pragma omp for
             for (i=0;i<numComponents2;++i) {
                 itmp=comps[i].n;
                 for (p=p_min; p<=p_max;++p) {
                    if ( p != distribution->mpi_info->rank) {
                       if ( (distribution->firstComponent[p]<=itmp) &&  (itmp<distribution->firstComponent[p+1])) {
                          comps[i].p=p;
                       }
                    }
                 }
             }
         } /* end parallel region */
         /* sort the comps py processor: */
         qsort(comps,comps,sizeof(Paso_SharedComponents_Container),Paso_SharedComponents_compar);

         /* find the first remote processor */
         nnnn=0;
         while (nnn<numComponents2) {
           if (comps[nnn].p != -1) break;
           nnn++;
         }
         /* the send case looks for the components needed elsewhere 
            these are the entries 0,...,nnn */

         /* the receive case looks for the component stored elsewhere 
            these are the entries nnn,...,numComponents2              */
     Paso_SharedComponents_createAndSortContainer(r_numComponents,comps,m,p,distribution);
     Paso_SharedComponents_createAndSortContainer(s_numComponents,comps,m,p,distribution);

           
    neighbor=NULL;
    offsetInShared
         numNeighbors=0;
         p=-1;
         for (i=nnn;i<numComponents2;++i) {
            if (comps[i].p != p ) {
              offsetInShared[numNeighbors]=i-nnn;
              neighbor[numNeighbors]=comps[i].p;
              fprintf("Lutz: detected neighbor %d %d\n",neighbor[numNeighbors], offsetInShared[numNeighbors]);
              numNeighbors++;
            }
         }
         offsetInShared[numNeighbors]=numComponents2-nnn+1;
   }
         
  MEMFREE(comps);
  
for (i=0; i<numComponents;++i) {
n=globalCompontents[i];
for (p=p_min; p<=p_max;++p) {
   if (distribution->firstComponent[p]<=n N<distribution->firstComponent[p+1]) {
      recv_processor[i].i=i;
      recv_processor[i].p=p;
   }
}
}

sort(recv_processor) by processor 
set recv_component_ordering[i]=recv_processor[i].i
find recv_numNeighbours
set recv_neighbours, offsetInRecvBuffer
#endif

  out=MEMALLOC(1,Paso_SharedComponents);
  if (!Paso_checkPtr(out)) {
      out->numComponents=numComponents*m;
      out->mpi_info = Paso_MPIInfo_getReference(distribution->mpi_info);
      out->globalComponent=NULL;
      out->numNeighbors=0;
      out->neighbor=NULL;
      out->shared=NULL;         
      out->offsetInShared=NULL; 
      out->distribution=Paso_Distribution_getReference(distribution);
      out->myNumComponents=Paso_Distribution_getMyNumComponents(out->distribution);
      out->numSharedComponents=out->numComponents-out->myNumComponents;
      out->reference_counter=1;

      out->globalComponent=MEMALLOC(out->numComponents,index_t);
      out->shared=MEMALLOC(out->numSharedComponents,index_t);
      if (! (Paso_checkPtr(out->globalComponent) || Paso_checkPtr(out->shared) ) ) {
         if (globalComponent == NULL) {
            offset=Paso_Distribution_getFirstComponent(out->distribution)+b;
            #pragma omp parallel for private(i,j,itmp)
            for (i=0;i<numComponents;++i){
                itmp=m*i+offset;
                for (j=0;j<m;++j)out->globalComponent[m*i+j]=itmp+j;
            }
         } else {
            #pragma omp parallel for private(i,j,itmp)
            for (i=0;i<numComponents;++i){
                itmp=m*globalComponent[i]+b;
                for (j=0;j<m;++j)out->globalComponent[m*i+j]=itmp+j;
            }
         }
      }
/*
  dim_t ordering;
  dim_t numNeighbors;      
  Paso_MPI_rank* neighbor; 
  index_t* offsetInShared; 
*/
     Paso_setError(SYSTEM_ERROR,"Paso_SharedComponents_alloc: not implemented yet.");

  }
  if (Paso_noError()) {
     return out;
  } else {
     Paso_SharedComponents_free(out);
     return NULL;
  }
}

/* returns a reference to in */

Paso_SharedComponents* Paso_SharedComponents_getReference(Paso_SharedComponents* in) {
     if (in!=NULL) {
        ++(in->reference_counter);
     }
     return in;
}
  
/* deallocates a SharedComponents: */

void Paso_SharedComponents_free(Paso_SharedComponents* in) {
  if (in!=NULL) {
     in->reference_counter--;
     if (in->reference_counter<=0) {
        MEMFREE(in->globalComponent);
        MEMFREE(in->ordering);
        Paso_Distribution_free(in->distribution);
        MEMFREE(in->neighbor);
        MEMFREE(in->shared);
        MEMFREE(in->offsetInShared);
        Paso_MPIInfo_free(in->mpi_info);
        MEMFREE(in);
        #ifdef Paso_TRACE
        printf("Paso_SharedComponents_dealloc: system matrix pattern as been deallocated.\n");
        #endif
     }
   }
}
