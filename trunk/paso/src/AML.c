
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

/* Paso: AML preconditioner:

This is just a collection of older code. This does not compile.

*/

/************************************************************************************/

/* Author: artak@uq.edu.au                                */

/************************************************************************************/
   
#include "Paso.h"
#include "Preconditioner.h"
#include "Options.h"
#include "PasoUtil.h"
#include "UMFPACK.h"
#include "MKL.h"
#include "SystemMatrix.h"
#include "Coarsening.h"
#include "BlockOps.h"

#ifndef INC_COARSE
#define INC_COARSE

#include "SystemMatrix.h"


/* Remove:
#define PASO_COARSENING_IN_F TRUE
#define PASO_COARSENING_IN_C FALSE
*/

void Paso_Coarsening_Local(index_t* marker_F, Paso_SparseMatrix* A,  Paso_Options* options);
void Paso_Coarsening_Local_Aggregation(Paso_SparseMatrix* A, index_t* marker_F, const double theta);
void Paso_Coarsening_Local_Aggregation_blk(Paso_SparseMatrix* A, index_t* marker_F, const double theta);
void Paso_Coarsening_Local_YS(Paso_SparseMatrix* A, index_t* marker_F, const double theta);
void Paso_Coarsening_Local_YS_blk(Paso_SparseMatrix* A, index_t* marker_F, const double theta);

void Paso_Coarsening_Local_RS(Paso_SparseMatrix* A, index_t* marker_F, double theta);

void Paso_Coarsening_Local_Partition(Paso_Pattern* pattern,index_t* marker_F);
void Paso_Coarsening_Local_greedy(Paso_Pattern* pattern, index_t* marker_F);



/*=== REVISE ============*/
void Paso_Coarsening_Local_greedy_color(Paso_Pattern* pattern, index_t* marker_F);
void Paso_Coarsening_Local_greedy_diag(Paso_SparseMatrix* A, index_t* marker_F, double thershold);

void Paso_Coarsening_Local_YS_plus(Paso_SparseMatrix* A, index_t* marker_F, double alpha, double taw, double delta);
void Paso_Coarsening_Local_Standard(Paso_SparseMatrix* A, index_t* marker_F, double theta);
void Paso_Coarsening_Local_greedy_RS(Paso_SparseMatrix* A, index_t* marker_F, double theta);
void Paso_Coarsening_Local_greedy_Agg(Paso_SparseMatrix* A, index_t* marker_F, double theta);
/*dim_t how_many(dim_t n,dim_t* S_i, int value1, dim_t* addedSet, int value2);*/
void Paso_Coarsening_Local_Standard_Block(Paso_SparseMatrix* A, index_t* marker_F, double theta);

dim_t how_many(dim_t i,Paso_Pattern * S, bool_t transpose);
dim_t arg_max(dim_t n, dim_t* lambda, dim_t mask);
Paso_Pattern* Paso_Coarsening_Local_getTranspose(Paso_Pattern* P);

void Paso_Coarsening_Local_getReport(dim_t n,index_t* marker_F);
void Paso_Coarsening_Local_Read(char *fileName,dim_t n,index_t* marker_F);
void Paso_Coarsening_Local_Write(char *fileName,dim_t n,index_t* marker_F);


#endif


/***********************************************************************************,amli->b_F*/

/* free all memory used by AMLI                                */

void Paso_Solver_AMLI_System_free(Paso_Solver_AMLI_System * in) {
     dim_t i;
     if (in!=NULL) {
        for (i=0;i<in->block_size;++i) {
          Paso_Solver_AMLI_free(in->amliblock[i]);
          Paso_SparseMatrix_free(in->block[i]);
        }
        MEMFREE(in);
     }
}


/* free all memory used by AMLI                                */

void Paso_Solver_AMLI_free(Paso_Solver_AMLI * in) {
     if (in!=NULL) {
	Paso_Preconditioner_LocalSmoother_free(in->Smoother);
	
        MEMFREE(in->inv_A_FF);
        MEMFREE(in->A_FF_pivot);
        Paso_SparseMatrix_free(in->A_FC);
        Paso_SparseMatrix_free(in->A_CF);
        Paso_SparseMatrix_free(in->A);
        if(in->coarsest_level==TRUE) {
        #ifdef MKL
          Paso_MKL_free1(in->AOffset1);
          Paso_SparseMatrix_free(in->AOffset1);
        #else
          #ifdef UMFPACK
          Paso_UMFPACK1_free((Paso_UMFPACK_Handler*)(in->solver));
          #endif
        #endif
        }
        MEMFREE(in->rows_in_F);
        MEMFREE(in->rows_in_C);
        MEMFREE(in->mask_F);
        MEMFREE(in->mask_C);
        MEMFREE(in->x_F);
        MEMFREE(in->b_F);
        MEMFREE(in->x_C);
        MEMFREE(in->b_C);
        in->solver=NULL;
        Paso_Solver_AMLI_free(in->AMLI_of_Schur);
        MEMFREE(in->b_C);
        MEMFREE(in);
     }
}

/************************************************************************************/

/*   constructs the block-block factorization of 

        [ A_FF A_FC ]
   A_p= 
        [ A_CF A_FF ]

to 

  [      I         0  ]  [ A_FF 0 ] [ I    invA_FF*A_FF ]  
  [ A_CF*invA_FF   I  ]  [   0  S ] [ 0          I      ] 


   where S=A_FF-ACF*invA_FF*A_FC within the shape of S

   then AMLI is applied to S again until S becomes empty 

*/
Paso_Solver_AMLI* Paso_Solver_getAMLI(Paso_SparseMatrix *A_p,dim_t level,Paso_Options* options) {
  Paso_Solver_AMLI* out=NULL;
  Paso_Pattern* temp1=NULL;
  Paso_Pattern* temp2=NULL;
  bool_t verbose=options->verbose;
  dim_t n=A_p->numRows;
  dim_t n_block=A_p->row_block_size;
  index_t* mis_marker=NULL;  
  index_t* counter=NULL;
  index_t iPtr,*index, *where_p, iPtr_s;
  dim_t i,j;
  Paso_SparseMatrix * schur=NULL;
  Paso_SparseMatrix * schur_withFillIn=NULL;
  double S=0;
  
  
 /* char filename[8];
  sprintf(filename,"AMLILevel%d",level);
  
 Paso_SparseMatrix_saveMM(A_p,filename);
  */
  
  /*Make sure we have block sizes 1*/
  if (A_p->col_block_size>1) {
     Esys_setError(TYPE_ERROR,"Paso_Solver_getAMLI: AMLI requires column block size 1.");
     return NULL;
  }
  if (A_p->row_block_size>1) {
     Esys_setError(TYPE_ERROR,"Paso_Solver_getAMLI: AMLI requires row block size 1.");
     return NULL;
  }
  out=MEMALLOC(1,Paso_Solver_AMLI);
  /* identify independend set of rows/columns */
  mis_marker=TMPMEMALLOC(n,index_t);
  counter=TMPMEMALLOC(n,index_t);
  if ( !( Esys_checkPtr(mis_marker) || Esys_checkPtr(counter) || Esys_checkPtr(out)) ) {
     out->AMLI_of_Schur=NULL;
     out->inv_A_FF=NULL;
     out->A_FF_pivot=NULL;
     out->A_FC=NULL;
     out->A_CF=NULL;
     out->rows_in_F=NULL;
     out->rows_in_C=NULL;
     out->mask_F=NULL;
     out->mask_C=NULL;
     out->x_F=NULL;
     out->b_F=NULL;
     out->x_C=NULL;
     out->b_C=NULL;
     out->A=Paso_SparseMatrix_getReference(A_p);
     out->Smoother=NULL;
     out->solver=NULL;
     /*out->GS=Paso_Solver_getGS(A_p,verbose);*/
     out->level=level;
     out->n=n;
     out->n_F=n+1;
     out->n_block=n_block;
     out->post_sweeps=options->post_sweeps;
     out->pre_sweeps=options->pre_sweeps;
     
     if (level==0 || n<=options->min_coarse_matrix_size) {
         out->coarsest_level=TRUE;
         #ifdef MKL
                  out->AOffset1=Paso_SparseMatrix_alloc(MATRIX_FORMAT_BLK1 + MATRIX_FORMAT_OFFSET1, out->A->pattern,1,1, FALSE);
                  #pragma omp parallel for private(i) schedule(static)
                  for (i=0;i<out->A->len;++i) {
                       out->AOffset1->val[i]=out->A->val[i];
                  }
         #else
            #ifdef UMFPACK 
            #else
            out->Smoother=Paso_Preconditioner_LocalSmoother_alloc(A_p,TRUE,verbose);
            #endif
         #endif
     } else {
         out->coarsest_level=FALSE;
	 out->Smoother=Paso_Preconditioner_LocalSmoother_alloc(A_p,TRUE,verbose);
 
	 Paso_Coarsening_Local(mis_marker,A_p, options->coarsening_threshold, options->coarsening_method);
	         
        #pragma omp parallel for private(i) schedule(static)
        for (i = 0; i < n; ++i) {
	   mis_marker[i]=(mis_marker[i]== PASO_COARSENING_IN_F);
	   counter[i]=mis_marker[i];
	    
        }
        
        out->n_F=Paso_Util_cumsum(n,counter);
        
        if (out->n_F==0) {
           out->coarsest_level=TRUE;
           level=0;
           if (verbose) { 
               /*printf("AMLI coarsening eliminates all unknowns, switching to Jacobi preconditioner.\n");*/
               printf("AMLI coarsening does not eliminate any of the unknowns, switching to Jacobi preconditioner.\n");
           }
        }
        else if (out->n_F==n) {
          out->coarsest_level=TRUE;
           level=0;
           if (verbose) { 
               /*printf("AMLI coarsening eliminates all unknowns, switching to Jacobi preconditioner.\n");*/
               printf("AMLI coarsening eliminates all of the unknowns, switching to Jacobi preconditioner.\n");
          
            }
        } else {
     
              if (Esys_noError()) {
                 
                 /*#pragma omp parallel for private(i) schedule(static)
                 for (i = 0; i < n; ++i) counter[i]=mis_marker[i];
                 out->n_F=Paso_Util_cumsum(n,counter);
                 */
                 
                 /*if(level==3) {
                   printf("##TOTAL: %d, ELIMINATED: %d\n",n,out->n_F);
                   for (i = 0; i < n; ++i) {
                    printf("##%d %d\n",i,mis_marker[i]);
                   }
                 }*/
                 
                 out->mask_F=MEMALLOC(n,index_t);
                 out->rows_in_F=MEMALLOC(out->n_F,index_t);
                 out->inv_A_FF=MEMALLOC(n_block*n_block*out->n_F,double);
                 out->A_FF_pivot=NULL; /* later use for block size>3 */
                 if (! (Esys_checkPtr(out->mask_F) || Esys_checkPtr(out->inv_A_FF) || Esys_checkPtr(out->rows_in_F) ) ) {
                    /* creates an index for F from mask */
                    #pragma omp parallel for private(i) schedule(static)
                    for (i = 0; i < out->n_F; ++i) out->rows_in_F[i]=-1;
                    #pragma omp parallel for private(i) schedule(static)
                    for (i = 0; i < n; ++i) {
                       if  (mis_marker[i]) {
                              out->rows_in_F[counter[i]]=i;
                              out->mask_F[i]=counter[i];
                       } else {
                              out->mask_F[i]=-1;
                       }
                    }
                    
                    /* Compute row-sum for getting rs(A_FF)^-1*/
                    #pragma omp parallel for private(i,iPtr,j,S) schedule(static)
                    for (i = 0; i < out->n_F; ++i) {
                      S=0;
      /*printf("[%d ]: [%d] -> ",i, out->rows_in_F[i]);*/
                      for (iPtr=A_p->pattern->ptr[out->rows_in_F[i]];iPtr<A_p->pattern->ptr[out->rows_in_F[i] + 1]; ++iPtr) {
                       j=A_p->pattern->index[iPtr];
      /*if (j==out->rows_in_F[i]) printf("diagonal %e",A_p->val[iPtr]);*/
                       if (mis_marker[j])
                           S+=A_p->val[iPtr];
                      }
      /*printf("-> %e \n",S);*/
                      out->inv_A_FF[i]=1./S;
                    }
                 }
              }
      
              /*check whether coarsening process actually makes sense to continue.
              if coarse matrix at least smaller by 30% then continue, otherwise we stop.*/
              if ((out->n_F*100/n)<30) {
                    level=1;
               }
             
              if ( Esys_noError()) {
                    /* if there are no nodes in the coarse level there is no more work to do */
                    out->n_C=n-out->n_F;
                   
                   /*if (out->n_F>500) */
                    out->rows_in_C=MEMALLOC(out->n_C,index_t);
                    out->mask_C=MEMALLOC(n,index_t);
                    if (! (Esys_checkPtr(out->mask_C) || Esys_checkPtr(out->rows_in_C) ) ) {
                         /* creates an index for C from mask */
                         #pragma omp parallel for private(i) schedule(static)
                         for (i = 0; i < n; ++i) counter[i]=! mis_marker[i];
                         Paso_Util_cumsum(n,counter);
                         #pragma omp parallel for private(i) schedule(static)
                         for (i = 0; i < out->n_C; ++i) out->rows_in_C[i]=-1;
                         #pragma omp parallel for private(i) schedule(static)
                         for (i = 0; i < n; ++i) {
                                  if  (! mis_marker[i]) {
                                      out->rows_in_C[counter[i]]=i;
                                      out->mask_C[i]=counter[i];
                                   } else {
                                      out->mask_C[i]=-1;
                                   }
                         }
                    }
              } 
              if ( Esys_noError()) {
                      /* get A_CF block: */
                      out->A_CF=Paso_SparseMatrix_getSubmatrix(A_p,out->n_C,out->n_F,out->rows_in_C,out->mask_F);
                      /* get A_FC block: */
                      out->A_FC=Paso_SparseMatrix_getSubmatrix(A_p,out->n_F,out->n_C,out->rows_in_F,out->mask_C);
                      /* get A_CC block: */
                      schur=Paso_SparseMatrix_getSubmatrix(A_p,out->n_C,out->n_C,out->rows_in_C,out->mask_C);
      
              }
              if ( Esys_noError()) {
                     /*find the pattern of the schur complement with fill in*/
                    temp1=Paso_Pattern_multiply(PATTERN_FORMAT_DEFAULT,out->A_CF->pattern,out->A_FC->pattern);
                    temp2=Paso_Pattern_binop(PATTERN_FORMAT_DEFAULT, schur->pattern, temp1);
                    schur_withFillIn=Paso_SparseMatrix_alloc(A_p->type,temp2,1,1, TRUE);
                    Paso_Pattern_free(temp1);
                    Paso_Pattern_free(temp2);
              }
              if ( Esys_noError()) {
                    /* copy values over*/ 
                    #pragma omp parallel for private(i,iPtr,j,iPtr_s,index,where_p) schedule(static)
                    for (i = 0; i < schur_withFillIn->numRows; ++i) {
                      for (iPtr=schur_withFillIn->pattern->ptr[i];iPtr<schur_withFillIn->pattern->ptr[i + 1]; ++iPtr) {
                         j=schur_withFillIn->pattern->index[iPtr];
                         iPtr_s=schur->pattern->ptr[i];
                         index=&(schur->pattern->index[iPtr_s]);
                         where_p=(index_t*)bsearch(&j,
                                              index,
                                              schur->pattern->ptr[i + 1]-schur->pattern->ptr[i],
                                              sizeof(index_t),
                                              Paso_comparIndex);
                         if (where_p!=NULL) {
                                schur_withFillIn->val[iPtr]=schur->val[iPtr_s+(index_t)(where_p-index)];
                         }
                       }
                    }
                    Paso_Solver_updateIncompleteSchurComplement(schur_withFillIn,out->A_CF,out->inv_A_FF,out->A_FF_pivot,out->A_FC);
                    out->AMLI_of_Schur=Paso_Solver_getAMLI(schur_withFillIn,level-1,options);
              }
              /* allocate work arrays for AMLI application */
              if (Esys_noError()) {
                         out->x_F=MEMALLOC(n_block*out->n_F,double);
                         out->b_F=MEMALLOC(n_block*out->n_F,double);
                         out->x_C=MEMALLOC(n_block*out->n_C,double);
                         out->b_C=MEMALLOC(n_block*out->n_C,double);
      
                         if (! (Esys_checkPtr(out->x_F) || Esys_checkPtr(out->b_F) || Esys_checkPtr(out->x_C) || Esys_checkPtr(out->b_C) ) ) {
                             #pragma omp parallel for private(i) schedule(static)
                             for (i = 0; i < out->n_F; ++i) {
                                         out->x_F[i]=0.;
                                         out->b_F[i]=0.;
                              }
                              #pragma omp parallel for private(i) schedule(static)
                              for (i = 0; i < out->n_C; ++i) {
                                     out->x_C[i]=0.;
                                     out->b_C[i]=0.;
                              }
                         }
              }
            Paso_SparseMatrix_free(schur);
            Paso_SparseMatrix_free(schur_withFillIn);
         }
     }
  }
  TMPMEMFREE(mis_marker);
  TMPMEMFREE(counter);

  if (Esys_noError()) {
      if (verbose && level>0 && !out->coarsest_level) {
         printf("AMLI: level: %d: %d unknowns eliminated. %d left.\n",level, out->n_F,out->n_C);
     }
     return out;
  } else  {
     Paso_Solver_AMLI_free(out);
     return NULL;
  }
}

/************************************************************************************/

/* apply AMLI precondition b-> x                               

     in fact it solves 

  [      I         0  ]  [ A_FF 0 ] [ I    invA_FF*A_FC ]  [ x_F ]  = [b_F]
  [ A_CF*invA_FF   I  ]  [   0  S ] [ 0          I      ]  [ x_C ]  = [b_C]

 in the form 

   b->[b_F,b_C] 
   x_F=invA_FF*b_F
   b_C=b_C-A_CF*x_F
   x_C=AMLI(b_C)
   b_F=b_F-A_FC*x_C
   x_F=invA_FF*b_F
   x<-[x_F,x_C]

 should be called within a parallel region                                              
 barrier synconization should be performed to make sure that the input vector available 

*/

void Paso_Solver_solveAMLI(Paso_Solver_AMLI * amli, double * x, double * b) {
     dim_t i;
     double time0=0;
     double *r=NULL, *x0=NULL,*x_F_temp=NULL;
     bool_t verbose=0;
     
     dim_t post_sweeps=amli->post_sweeps;
     dim_t pre_sweeps=amli->pre_sweeps;
     
     #ifdef UMFPACK 
          Paso_UMFPACK_Handler * ptr=NULL;
     #endif
     
     r=MEMALLOC(amli->n,double);
     x0=MEMALLOC(amli->n,double);
     x_F_temp=MEMALLOC(amli->n_F,double);
     
     if (amli->coarsest_level) {
      
      time0=Esys_timer();
      /*If all unknown are eliminated then Jacobi is the best preconditioner*/
      if (amli->n_F==0 || amli->n_F==amli->n) {
	 Paso_Preconditioner_LocalSmoother_solve(amli->A, amli->Smoother,x,b,1,FALSE);
      }
       else {
       #ifdef MKL
          Paso_MKL1(amli->AOffset1,x,b,verbose);
       #else
          #ifdef UMFPACK
             ptr=(Paso_UMFPACK_Handler *)(amli->solver);
             Paso_UMFPACK1(&ptr,amli->A,x,b,verbose);
             amli->solver=(void*) ptr;
          #else      
          Paso_Preconditioner_LocalSmoother_solve(amli->A,amli->Smoother,x,b,1,FALSE);
         #endif
       #endif
       }
       time0=Esys_timer()-time0;
       if (verbose) fprintf(stderr,"timing: DIRECT SOLVER: %e\n",time0);
       
     } else {
        /* presmoothing */
         time0=Esys_timer();
	 Paso_Preconditioner_LocalSmoother_solve(amli->A,amli->Smoother,x,b,pre_sweeps,FALSE);
         time0=Esys_timer()-time0;
         if (verbose) fprintf(stderr,"timing: Presmooting: %e\n",time0);
        /* end of presmoothing */
        
        
         time0=Esys_timer();
         #pragma omp parallel for private(i) schedule(static)
         for (i=0;i<amli->n;++i) r[i]=b[i];
         
         /*r=b-Ax*/ 
         Paso_SparseMatrix_MatrixVector_CSR_OFFSET0(-1.,amli->A,x,1.,r);
       
        /* r->[b_F,b_C]     */
        #pragma omp parallel for private(i) schedule(static)
        for (i=0;i<amli->n_F;++i) amli->b_F[i]=r[amli->rows_in_F[i]];
        
        #pragma omp parallel for private(i) schedule(static)
        for (i=0;i<amli->n_C;++i) amli->b_C[i]=r[amli->rows_in_C[i]];

        /* x_F=invA_FF*b_F  */
	Paso_Copy(amli->n_F, amli->x_F,amli->b_F);
        Paso_BlockOps_solveAll(1,amli->n_F,amli->inv_A_FF,amli->A_FF_pivot,amli->x_F);
        
        /* b_C=b_C-A_CF*x_F */
        Paso_SparseMatrix_MatrixVector_CSR_OFFSET0(-1.,amli->A_CF,amli->x_F,1.,amli->b_C);
        
        time0=Esys_timer()-time0;
        if (verbose) fprintf(stderr,"timing: Before next level: %e\n",time0);
        
        /* x_C=AMLI(b_C)     */
        Paso_Solver_solveAMLI(amli->AMLI_of_Schur,amli->x_C,amli->b_C);
        
        time0=Esys_timer();
        
        /* b_F=-A_FC*x_C */ 
        Paso_SparseMatrix_MatrixVector_CSR_OFFSET0(-1.,amli->A_FC,amli->x_C,0.,amli->b_F);
        /* x_F_temp=invA_FF*b_F  */
	Paso_Copy(amli->n_F, x_F_temp,amli->b_F);
	Paso_BlockOps_solveAll(1,amli->n_F,amli->inv_A_FF,amli->A_FF_pivot,x_F_temp);
        
        #pragma omp parallel for private(i) schedule(static)
        for (i=0;i<amli->n_F;++i) {
                 amli->x_F[i]+=x_F_temp[i];
        }
        
        /* x<-[x_F,x_C]     */
        #pragma omp parallel for private(i) schedule(static)
        for (i=0;i<amli->n;++i) {
            if (amli->mask_C[i]>-1) {
                 x[i]+=amli->x_C[amli->mask_C[i]];
            } else {
                 x[i]+=amli->x_F[amli->mask_F[i]];
            }
        }
        
        time0=Esys_timer()-time0;
        if (verbose) fprintf(stderr,"timing: After next level: %e\n",time0);

     /*postsmoothing*/
     time0=Esys_timer();
     Paso_Preconditioner_LocalSmoother_solve(amli->A,amli->Smoother,x,b,post_sweeps,TRUE);     
     time0=Esys_timer()-time0;
     if (verbose) fprintf(stderr,"timing: Postsmoothing: %e\n",time0);

     /*end of postsmoothing*/
     
     }
     MEMFREE(r);
     MEMFREE(x0);
     MEMFREE(x_F_temp);
     return;
}

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


/***********************************************************************************

Paso: Coarsening strategies (no MPI)   

the given matrix A is splitted in to the form 


[ A_FF A_FC ]
A  =  [           ]
[ A_CF A_CC ]

such that unknowns/equations in set C are weakly connected via A_CC and strongly connected 
to at least one unknowns/equation in F via A_CF and A_FC. The unknowns/equations in C and F
are marked in the marker_F by FALSE and TRUE respectively.

The weak/strong connection is controlled by coarsening_threshold.

three strategies are implemented:

a) YAIR_SHAPIRA (YS): |a_ij| >= theta * |a_ii|
b) Ruge-Stueben (RS): |a_ij| >= theta * max_(k<>i) |a_ik|
c) Aggregation :     |a_ij|^2 >= theta**2 * |a_ii||a_jj|

where theta = coarsening_threshold/maximum_pattern_degree 

Remark:

- a strong connection in YAIR_SHAPIRA is a strong connection in Aggregation

*/
/************************************************************************************

Author: artak@uq.edu.au, l.gross@uq.edu.au                   

************************************************************************************/

#include "Coarsening.h"
#include "PasoUtil.h"
#include <limits.h>



/*Used in Paso_Coarsening_Local_RS_MI*/

/*Computes how many connections unknown i have in S.
bool_t transpose - TRUE if we want to compute how many strong connection of i in S^T, FALSE otherwise.
Note that if we first transpose S and then call method on S^T then "transpose" should be set to FALSE.
*/
dim_t how_many(dim_t i,Paso_Pattern * S, bool_t transpose) {
   dim_t j,n;
   dim_t total,ltotal;
   index_t *index,*where_p;
   
   /*index_t iptr;*/
   total=0;
   ltotal=0;
   
   n=S->numOutput;
   
   if(transpose) {
      #pragma omp parallel 
      {
	 ltotal=0;
	 #pragma omp for private(j,index,where_p,ltotal) schedule(static)
	 for (j=0;j<n;++j) {
	    index=&(S->index[S->ptr[j]]);
	    where_p=(index_t*)bsearch(&i,
				      index,
				      S->ptr[j + 1]-S->ptr[j],
				      sizeof(index_t),
				      Paso_comparIndex);
				      if (where_p!=NULL) {
					 ltotal++;
				      }
	 }
      }
      #pragma omp critical
      {
	 total+=ltotal;
      }
      
   }
   else {
      total=S->ptr[i+1]-S->ptr[i];
      /*#pragma omp parallel for private(iptr) schedule(static)*/
      /*for (iptr=S->ptr[i];iptr<S->ptr[i+1]; ++iptr) {
      if(S->index[iptr]!=i && marker_F[S->index[iptr]]==IS_AVAILABLE)
	 total++;
   }*/
      
   }
   
   if (total==0) total=IS_NOT_AVAILABLE;
   
   return total; 
}





Paso_Pattern* Paso_Coarsening_Local_getTranspose(Paso_Pattern* P){
   
   Paso_Pattern *outpattern=NULL;
   
   dim_t C=P->numInput;
   dim_t F=P->numOutput-C;
   dim_t n=C+F;
   dim_t i,j;
   index_t iptr;
   Paso_IndexListArray* index_list = Paso_IndexListArray_alloc(C);
   
   
   /*#pragma omp parallel for private(i,iptr,j) schedule(static)*/
   for (i=0;i<n;++i) {
      for (iptr=P->ptr[i];iptr<P->ptr[i+1]; ++iptr) {
	 j=P->index[iptr];
	 Paso_IndexListArray_insertIndex(index_list,i,j);
      }
   }
   
   outpattern=Paso_Pattern_fromIndexListArray(0, index_list,0,n,0);
   
   /* clean up */
   Paso_IndexListArray_free(index_list);
   
   return outpattern;
}


/************** BLOCK COARSENENING *********************/

void Paso_Coarsening_Local_Standard_Block(Paso_SparseMatrix* A, index_t* marker_F, double theta)
{
   const dim_t n=A->numRows;
   
   dim_t i,j,k;
   index_t iptr,jptr;
   /*index_t *index,*where_p;*/
   double threshold,max_offdiagonal;
   dim_t *lambda;   /*mesure of importance */
   dim_t maxlambda=0;
   index_t index_maxlambda=0;
   double time0=0;
   bool_t verbose=0;
   dim_t n_block=A->row_block_size;
   
   double fnorm=0;
   dim_t bi;
   
   Paso_Pattern *S=NULL;
   Paso_Pattern *ST=NULL;
   Paso_IndexListArray* index_list = Paso_IndexListArray_alloc(n);
   
   time0=Esys_timer();
   k=0;
   /*Paso_Coarsening_Local_getReport(n,marker_F);*/
   /*printf("Blocks %d %d\n",n_block,A->len);*/
   
   /*S_i={j \in N_i; i strongly coupled to j}*/
   #pragma omp parallel for private(i,bi,fnorm,iptr,max_offdiagonal,threshold,j) schedule(static)
   for (i=0;i<n;++i) {
      if(marker_F[i]==IS_AVAILABLE) {
	 max_offdiagonal = DBL_MIN;
	 for (iptr=A->pattern->ptr[i];iptr<A->pattern->ptr[i+1]; ++iptr) {
	    j=A->pattern->index[iptr];
	    if( j != i){
	       fnorm=0;
	       for(bi=0;bi<n_block*n_block;++bi)
	       {
		  fnorm+=A->val[iptr*n_block*n_block+bi]*A->val[iptr*n_block*n_block+bi];
	       }
	       fnorm=sqrt(fnorm);
	       max_offdiagonal = MAX(max_offdiagonal,fnorm);
	    }
	 }
	 threshold = theta*max_offdiagonal;
	 for (iptr=A->pattern->ptr[i];iptr<A->pattern->ptr[i+1]; ++iptr) {
	    j=A->pattern->index[iptr];
	    if( j != i){
	       fnorm=0;
	       for(bi=0;bi<n_block*n_block;++bi)
	       {
		  fnorm+=A->val[iptr*n_block*n_block+bi]*A->val[iptr*n_block*n_block+bi];
	       }
	       fnorm=sqrt(fnorm);
	       if(fnorm>=threshold) {
		  Paso_IndexListArray_insertIndex(index_list,i,j);
	       }
	    }
	 }
      }
   }
   
   S=Paso_Pattern_fromIndexListArray(0,index_list,0,A->pattern->numInput,0);
   ST=Paso_Coarsening_Local_getTranspose(S);
   
   /*printf("Patterns len %d %d\n",S->len,ST->len);*/
   
   time0=Esys_timer()-time0;
   if (verbose) fprintf(stdout,"timing: RS filtering and pattern creation: %e\n",time0);
   
   lambda=TMPMEMALLOC(n,dim_t);
   
   #pragma omp parallel for private(i) schedule(static)
   for (i=0;i<n;++i) { lambda[i]=IS_NOT_AVAILABLE; }
   
   /*S_i={j \in N_i; i strongly coupled to j}*/
   
   /*
   #pragma omp parallel for private(i,iptr,lk) schedule(static)
   for (i=0;i<n;++i) {
      lk=0;
      for (iptr=A->pattern->ptr[i];iptr<A->pattern->ptr[i+1]; ++iptr) {
	 if(ABS(A->val[iptr])>1.e-15 && A->pattern->index[iptr]!=i )
	    lk++;
}
#pragma omp critical
k+=lk;
if(k==0) {
   marker_F[i]=TRUE;
}
}
*/
   
   
   k=0;
   maxlambda=0;
   
   time0=Esys_timer();
   
   for (i=0;i<n;++i) {
      if(marker_F[i]==IS_AVAILABLE) {
	 lambda[i]=how_many(k,ST,FALSE);   /* if every row is available then k and i are the same.*/
	 /*lambda[i]=how_many(i,S,TRUE);*/
	 /*printf("lambda[%d]=%d, k=%d \n",i,lambda[i],k);*/
	 k++;
	 if(maxlambda<lambda[i]) {
	    maxlambda=lambda[i];
	    index_maxlambda=i;
	 }
      }
   }
   
   k=0;
   time0=Esys_timer()-time0;
   if (verbose) fprintf(stdout,"timing: Lambdas computations at the begining: %e\n",time0);
   
   time0=Esys_timer();
   
   /*Paso_Coarsening_Local_getReport(n,marker_F);*/
   
   while (Paso_Util_isAny(n,marker_F,IS_AVAILABLE)) {
      
      if(index_maxlambda<0) {
	 break;
      }
      
      i=index_maxlambda;
      if(marker_F[i]==IS_AVAILABLE) {
	 marker_F[index_maxlambda]=FALSE;
	 lambda[index_maxlambda]=IS_NOT_AVAILABLE;
	 for (iptr=ST->ptr[i];iptr<ST->ptr[i+1]; ++iptr) {
	    j=ST->index[iptr];
	    if(marker_F[j]==IS_AVAILABLE) {
	       marker_F[j]=TRUE;
	       lambda[j]=IS_NOT_AVAILABLE;
	       for (jptr=S->ptr[j];jptr<S->ptr[j+1]; ++jptr) {
		  k=S->index[jptr];
		  if(marker_F[k]==IS_AVAILABLE) {
		     lambda[k]++; 
		  }
	       }
	    }
	 }
	 for (iptr=S->ptr[i];iptr<S->ptr[i+1]; ++iptr) {
	    j=S->index[iptr];
	    if(marker_F[j]==IS_AVAILABLE) {
	       lambda[j]--; 
	    }
	 }
      }
      
      /* Used when transpose of S is not available */
      /*
      for (i=0;i<n;++i) {
	 if(marker_F[i]==IS_AVAILABLE) {
	    if (i==index_maxlambda) {
	       marker_F[index_maxlambda]=FALSE;
	       lambda[index_maxlambda]=-1;
	       for (j=0;j<n;++j) {
		  if(marker_F[j]==IS_AVAILABLE) {
		     index=&(S->index[S->ptr[j]]);
		     where_p=(index_t*)bsearch(&i,
		     index,
		     S->ptr[j + 1]-S->ptr[j],
		     sizeof(index_t),
		     Paso_comparIndex);
		     if (where_p!=NULL) {
			marker_F[j]=TRUE;
			lambda[j]=-1;
			for (iptr=S->ptr[j];iptr<S->ptr[j+1]; ++iptr) {
			   k=S->index[iptr];
			   if(marker_F[k]==IS_AVAILABLE) {
			      lambda[k]++; 
   }
   }
   }
   
   }
   }
   
   }
   }
   }
   */
      index_maxlambda=arg_max(n,lambda, IS_NOT_AVAILABLE);
   }
   
   time0=Esys_timer()-time0;
   if (verbose) fprintf(stdout,"timing: Loop : %e\n",time0);
   
   /*Paso_Coarsening_Local_getReport(n,marker_F);*/
   
   #pragma omp parallel for private(i) schedule(static)
   for (i=0;i<n;++i) 
      if(marker_F[i]==IS_AVAILABLE) {
	 marker_F[i]=TRUE;
      }
      
      /*Paso_Coarsening_Local_getReport(n,marker_F);*/
      
      TMPMEMFREE(lambda);
   
   /* clean up */
   Paso_IndexListArray_free(index_list);
   Paso_Pattern_free(S);
   
   /* swap to TRUE/FALSE in marker_F */
   #pragma omp parallel for private(i) schedule(static)
   for (i=0;i<n;i++) marker_F[i]=(marker_F[i]==TRUE)? TRUE : FALSE;
   
}



#undef IS_AVAILABLE
#undef IS_NOT_AVAILABLE
#undef TRUE 
#undef FALSE

#undef IS_UNDECIDED
#undef IS_STRONG 
#undef IS_WEAK 

#undef TRUEB 
#undef TRUEB 

