
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/


/**************************************************************/

/* Paso: SparseMatrix_AMGcomponents */

/**************************************************************/

/* Author: Artak Amirbekyan, artak@uq.edu.au */

/**************************************************************/

#include "Paso.h"
#include "SparseMatrix.h"
#include "PasoUtil.h"

/**************************************************************

    Methods nessecary for AMG preconditioner

*/

/* Prolongation matrix is constracted by W_FC as P<- [W_FC]
                                                     [I_CC]
*/

Paso_SparseMatrix* Paso_SparseMatrix_getProlongation(Paso_SparseMatrix* W, index_t* mis_marker){
 
  Paso_Pattern *outpattern=NULL;
  Paso_SparseMatrix *out=NULL;
  index_t iptr,wptr;
  
  Paso_IndexList* index_list=NULL;
  dim_t n=W->numRows+W->numCols;
  dim_t i,j,k=0;

  index_list=TMPMEMALLOC(n,Paso_IndexList);
   if (! Paso_checkPtr(index_list)) {
        #pragma omp parallel for private(i) schedule(static)
        for(i=0;i<n;++i) {
             index_list[i].extension=NULL;
             index_list[i].n=0;
        }
    }


    for (i=0;i<n;++i) {
      if (mis_marker[i]) {
          for (iptr=W->pattern->ptr[k];iptr<W->pattern->ptr[k+1]; ++iptr) {
            j=W->pattern->index[iptr];
            Paso_IndexList_insertIndex(&(index_list[i]),j);
          }
          k++;
      }
      else {
          Paso_IndexList_insertIndex(&(index_list[i]),i-k);
      }
    }
    
    outpattern=Paso_IndexList_createPattern(0, n,index_list,0,W->numCols,0);
    out=Paso_SparseMatrix_alloc(W->type,outpattern,1,1,TRUE);
    
    k=0;

    for (i=0;i<n;++i) {
      if (mis_marker[i]) {
        wptr=W->pattern->ptr[k];
        for (iptr=out->pattern->ptr[i];iptr<out->pattern->ptr[i+1]; ++iptr) {
            out->val[iptr]=W->val[wptr];
            wptr++;
           }
        k++;   
      }
      else {
           for (iptr=out->pattern->ptr[i];iptr<out->pattern->ptr[i+1]; ++iptr){
             out->val[iptr]=1;
           }
      }
    }
    
     /* clean up */
   if (index_list!=NULL) {
        #pragma omp parallel for private(i) schedule(static)
        for(i=0;i<n;++i) Paso_IndexList_free(index_list[i].extension);
     }
    TMPMEMFREE(index_list);
    Paso_Pattern_free(outpattern);
    return out;
}


/* Restriction matrix R=P^T */

Paso_SparseMatrix* Paso_SparseMatrix_getRestriction(Paso_SparseMatrix* P){
 
  Paso_Pattern *outpattern=NULL;
  Paso_SparseMatrix *out=NULL;
  
  Paso_IndexList* index_list=NULL;
  dim_t C=P->numCols;
  dim_t F=P->numRows-C;
  dim_t n=C+F;
  dim_t i,j,k=0;
  index_t iptr,jptr;

  index_list=TMPMEMALLOC(C,Paso_IndexList);
   if (! Paso_checkPtr(index_list)) {
        #pragma omp parallel for private(i) schedule(static)
        for(i=0;i<C;++i) {
             index_list[i].extension=NULL;
             index_list[i].n=0;
        }
    }
  

    for (i=0;i<n;++i) {
          for (iptr=P->pattern->ptr[i];iptr<P->pattern->ptr[i+1]; ++iptr) {
             j=P->pattern->index[iptr];
             Paso_IndexList_insertIndex(&(index_list[j]),i);
        }
    }
   
    outpattern=Paso_IndexList_createPattern(0, C,index_list,0,C+F,0);
    out=Paso_SparseMatrix_alloc(P->type,outpattern,1,1,TRUE);
    
    for (i=0;i<out->numRows;++i) {
           for (iptr=out->pattern->ptr[i];iptr<out->pattern->ptr[i+1]; ++iptr) {
                 j=out->pattern->index[iptr];
                  /*This for later will be replaced by bsearch!!*/
                  for (jptr=P->pattern->ptr[j];jptr<P->pattern->ptr[j+1]; ++jptr) {
                        k=P->pattern->index[jptr];
                        if(k==i) {
                             out->val[iptr]=P->val[jptr];
                        }
                  }
           }
      }
    
     /* clean up */
   if (index_list!=NULL) {
        #pragma omp parallel for private(i) schedule(static)
        for(i=0;i<C;++i) Paso_IndexList_free(index_list[i].extension);
     }
    TMPMEMFREE(index_list);
    Paso_Pattern_free(outpattern);
    return out;
}


void Paso_SparseMatrix_updateWeights(Paso_SparseMatrix* A,Paso_SparseMatrix* W_FC, index_t* mis_marker){
 
  double *alpha;
  double *beta;
  double a_ii=0;
  double sum_all_neg,sum_all_pos,sum_strong_neg,sum_strong_pos;
  
  dim_t n=A->numRows;
  dim_t F=W_FC->numRows;
  dim_t i,j,k;
  
  index_t iPtr,dptr=0;
  /*index_t *index, *where_p;*/
 
  alpha=TMPMEMALLOC(F,double);
  beta=TMPMEMALLOC(F,double);


  k=0;
  for (i = 0; i < n; ++i) {
      if(mis_marker[i]) {
            alpha[k]=0;
            beta[k]=0;
            sum_all_neg=0;
            sum_all_pos=0;
            sum_strong_neg=0;
            sum_strong_pos=0;
            /*compute beta_i=sum_{j\inN_i}a^{+}_ij and alpha_i=sum_{j\inN_i}a^{-}_ij and denominators for both alpha_i and beta_i*/
            for (iPtr=A->pattern->ptr[i];iPtr<A->pattern->ptr[i + 1]; ++iPtr) {
                 j=A->pattern->index[iPtr];
                 if(j!=i) {
                        if(A->val[iPtr]<0) {
                          sum_all_neg+=A->val[iPtr];
                          if(!mis_marker[j]) {
                            sum_strong_neg+=A->val[iPtr];
                          }
                        }
                        else if(A->val[iPtr]>0) {
                          sum_all_pos+=A->val[iPtr];
                          if(!mis_marker[j]) {
                            sum_strong_pos+=A->val[iPtr];
                          }
                        }
                 }
                 else {
                      a_ii=A->val[iPtr];
                      dptr=iPtr;
                 }
             }

            if(sum_strong_neg!=0) {
                 alpha[k]=sum_all_neg/(sum_strong_neg);
            }
            if(sum_strong_pos!=0) {
                 beta[k]=sum_all_pos/(sum_strong_pos);
            }
/*            else {
              a_ii+=beta[k];
              A->val[dptr]=a_ii;
              beta[k]=0;
            }
*/           alpha[k]=alpha[k]/(a_ii);
            beta[k]=beta[k]/(a_ii);
          k++;
       /*printf("Got in row=%d, alpha[%d]=%e, beta[%d]=%e, a_den=%e, b_den=%e \n",i,k-1,alpha[k-1],k-1,beta[k-1],alpha_den[k-1],beta_den[k-1]);*/
      }
   }
      #pragma omp parallel for private(i,iPtr,j) schedule(static)
      for (i = 0; i < W_FC->numRows; ++i) {
            for (iPtr=W_FC->pattern->ptr[i];iPtr<W_FC->pattern->ptr[i + 1]; ++iPtr) {
                 j=W_FC->pattern->index[iPtr];
                   if(W_FC->val[iPtr]<0) {
                      W_FC->val[iPtr]=-alpha[i]*W_FC->val[iPtr];
                    }
                   else if(W_FC->val[iPtr]>0) {
                      W_FC->val[iPtr]=-beta[i]*W_FC->val[iPtr];
                   }
            }
         }
      
      TMPMEMFREE(alpha);
      TMPMEMFREE(beta);
}

  
/*A_c=R*A*P taking into account coarseneing */
/*Paso_SparseMatrix* Paso_Solver_getCoarseMatrix(Paso_SparseMatrix* A, Paso_SparseMatrix *R, Paso_SparseMatrix *P) {

Paso_SparseMatrix *temp=NULL;
Paso_SparseMatrix *out=NULL;

    temp=Paso_SparseMatrix_MatrixMatrix(A,P);
    out=Paso_SparseMatrix_MatrixMatrix(R,temp);
    
    Paso_SparseMatrix_free(temp);
      
return out;      
}
 */

Paso_SparseMatrix* Paso_SparseMatrix_MatrixMatrix(Paso_SparseMatrix* A, Paso_SparseMatrix* B) {
  index_t iptrA,iptrB,iptrC;
  dim_t i,j=0,k,l;
  Paso_Pattern* outpattern=NULL;
  Paso_SparseMatrix *out=NULL;
  double sum,b_lj=0;
  
  double time0=0;
  bool_t verbose=0;
  
  time0=Paso_timer();
  
  outpattern=Paso_Pattern_multiply(PATTERN_FORMAT_DEFAULT,A->pattern,B->pattern);
  out=Paso_SparseMatrix_alloc(A->type,outpattern,1,1, TRUE);
  
  time0=Paso_timer()-time0;
  if (verbose) fprintf(stdout,"timing: Paso_SparseMatrix_MatrixMatrix: Pattern creation: %e\n",time0);
  
  time0=Paso_timer();
  
  #pragma omp parallel for private(i,iptrC,j,sum,iptrA,k,b_lj,iptrB,l) schedule(static)
  for(i = 0; i < out->numRows; i++) {
    for(iptrC = out->pattern->ptr[i]; iptrC < out->pattern->ptr[i+1]; ++iptrC) {
      j = out->pattern->index[iptrC];
      sum=0;
      for(iptrA = A->pattern->ptr[i]; iptrA < A->pattern->ptr[i+1]; ++iptrA) {
            k=A->pattern->index[iptrA];
            /*can be replaced by bsearch */
            b_lj=0;
            for(iptrB = B->pattern->ptr[k]; iptrB < B->pattern->ptr[k+1]; ++iptrB) {
               l=B->pattern->index[iptrB];
               if(l==j) {
                  b_lj=B->val[iptrB];
                  break;
               }
            }
            sum+=A->val[iptrA]*b_lj;
      }
      out->val[iptrC]=sum;
   }
  }
  
  time0=Paso_timer()-time0;
  if (verbose) fprintf(stdout,"timing: Paso_SparseMatrix_MatrixMatrix: Matrix multiplication: %e\n",time0);
      
  Paso_Pattern_free(outpattern);

return out;      
}


Paso_SparseMatrix* Paso_Solver_getCoarseMatrix(Paso_SparseMatrix* A,Paso_SparseMatrix *R,Paso_SparseMatrix *P) {

  index_t iptrA_c,iptrR,iptrP,iptrA;
  dim_t i,j,k,l,m;
  double first_sum=0,second_sum=0,p_lj=0,a_kl=0,r_ik=0;
  Paso_Pattern* temp=NULL;
  Paso_Pattern* outpattern=NULL;
  Paso_SparseMatrix *A_c=NULL;
  
  double time0=0;
  bool_t verbose=0;
  
  time0=Paso_timer();
  
  temp=Paso_Pattern_multiply(PATTERN_FORMAT_DEFAULT,A->pattern,P->pattern);
  outpattern=Paso_Pattern_multiply(PATTERN_FORMAT_DEFAULT,R->pattern,temp);
  A_c=Paso_SparseMatrix_alloc(A->type,outpattern,1,1, TRUE);
  
  time0=Paso_timer()-time0;
  if (verbose) fprintf(stdout,"timing: Paso_Solver_getCoarseMatrix: Pattern creation: %e\n",time0);
  
  /*a^c_ij=sum_k^n(r_ik)sum_l^n(a_kl*P_lj)*/

  time0=Paso_timer();
  
  #pragma omp parallel for private(i,iptrA_c,j,second_sum,iptrR,k,first_sum,p_lj,iptrP,m,a_kl,r_ik) schedule(static)
  for(i = 0; i < A_c->numRows; i++) {
     for(iptrA_c = A_c->pattern->ptr[i]; iptrA_c < A_c->pattern->ptr[i+1]; ++iptrA_c) {
      j = A_c->pattern->index[iptrA_c];
      second_sum=0;
      for(iptrR = R->pattern->ptr[i]; iptrR < R->pattern->ptr[i+1]; ++iptrR) {
    	k = R->pattern->index[iptrR];
        first_sum=0;
        for(iptrA = A->pattern->ptr[k]; iptrA < A->pattern->ptr[k+1]; ++iptrA) {
         l= A->pattern->index[iptrA];
         p_lj=0;
        /*This loop can be replaced by bsearch */
         for(iptrP = P->pattern->ptr[l]; iptrP < P->pattern->ptr[l+1]; ++iptrP) {
           m=P->pattern->index[iptrP];
           if(m==j) {
            p_lj=P->val[iptrP];
            break;
           }
         }
         a_kl=A->val[iptrA];
         first_sum+=a_kl*p_lj;
        }
        r_ik=R->val[iptrR];
        second_sum+=r_ik*first_sum;
      }
      A_c->val[iptrA_c]=second_sum;
     }
  }
  time0=Paso_timer()-time0;
  if (verbose) fprintf(stdout,"timing: Paso_Solver_getCoarseMatrix: Matrix multiplication: %e\n",time0);
    
  Paso_Pattern_free(outpattern);
  Paso_Pattern_free(temp);
    
return A_c;
}

