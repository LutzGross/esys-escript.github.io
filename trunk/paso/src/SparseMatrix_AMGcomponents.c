
/*******************************************************
*
* Copyright (c) 2003-2010 by University of Queensland
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
  dim_t block_size=W->row_block_size;

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
    out=Paso_SparseMatrix_alloc(W->type,outpattern,block_size,block_size,FALSE);
    
    k=0;
    
    if (block_size==1) {
            for (i=0;i<n;++i) {
              if (mis_marker[i]) {
                wptr=W->pattern->ptr[k];
                for (iptr=out->pattern->ptr[i];iptr<out->pattern->ptr[i+1]; ++iptr) {
                    out->val[iptr*block_size*block_size]=W->val[wptr*block_size*block_size];
                    wptr++;
                   }
                k++;   
              }
              else {
                   iptr=out->pattern->ptr[i];
                   out->val[iptr*block_size*block_size]=1;
              }
            }
    } else if (block_size==2) {
            for (i=0;i<n;++i) {
              if (mis_marker[i]) {
                wptr=W->pattern->ptr[k];
                for (iptr=out->pattern->ptr[i];iptr<out->pattern->ptr[i+1]; ++iptr) {
                    out->val[iptr*block_size*block_size]=W->val[wptr*block_size*block_size];
                    out->val[iptr*block_size*block_size+1]=0;
                    out->val[iptr*block_size*block_size+2]=0;
                    out->val[iptr*block_size*block_size+3]=W->val[wptr*block_size*block_size+3];
                    wptr++;
                   }
                k++;   
              }
              else {
                   iptr=out->pattern->ptr[i];
                   out->val[iptr*block_size*block_size]=1;
                   out->val[iptr*block_size*block_size+1]=0;
                   out->val[iptr*block_size*block_size+2]=0;
                   out->val[iptr*block_size*block_size+3]=1;
              }
            }
    } else if (block_size==3) {
            for (i=0;i<n;++i) {
              if (mis_marker[i]) {
                wptr=W->pattern->ptr[k];
                for (iptr=out->pattern->ptr[i];iptr<out->pattern->ptr[i+1]; ++iptr) {
                    out->val[iptr*block_size*block_size]=W->val[wptr*block_size*block_size];
                    out->val[iptr*block_size*block_size+1]=0;
                    out->val[iptr*block_size*block_size+2]=0;
                    out->val[iptr*block_size*block_size+3]=0;
                    out->val[iptr*block_size*block_size+4]=W->val[wptr*block_size*block_size+4];
                    out->val[iptr*block_size*block_size+5]=0;
                    out->val[iptr*block_size*block_size+6]=0;
                    out->val[iptr*block_size*block_size+7]=0;
                    out->val[iptr*block_size*block_size+8]=W->val[wptr*block_size*block_size+8];
                    wptr++;
                   }
                k++;   
              }
              else {
                   iptr=out->pattern->ptr[i];
                   out->val[iptr*block_size*block_size]=1;
                   out->val[iptr*block_size*block_size+1]=0;
                   out->val[iptr*block_size*block_size+2]=0;
                   out->val[iptr*block_size*block_size+3]=0;
                   out->val[iptr*block_size*block_size+4]=1;
                   out->val[iptr*block_size*block_size+5]=0;
                   out->val[iptr*block_size*block_size+6]=0;
                   out->val[iptr*block_size*block_size+7]=0;
                   out->val[iptr*block_size*block_size+8]=1;
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
  dim_t block_size=P->row_block_size;
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
    out=Paso_SparseMatrix_alloc(P->type,outpattern,block_size,block_size,FALSE);
    
    
    if (block_size==1) {
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
    } else if (block_size==2) {
           for (i=0;i<out->numRows;++i) {
                 for (iptr=out->pattern->ptr[i];iptr<out->pattern->ptr[i+1]; ++iptr) {
                       j=out->pattern->index[iptr];
                        /*This for later will be replaced by bsearch!!*/
                        for (jptr=P->pattern->ptr[j];jptr<P->pattern->ptr[j+1]; ++jptr) {
                              k=P->pattern->index[jptr];
                              if(k==i) {
                                   out->val[iptr*block_size*block_size]=P->val[jptr*block_size*block_size];
                                   out->val[iptr*block_size*block_size+1]=P->val[jptr*block_size*block_size+2];
                                   out->val[iptr*block_size*block_size+2]=P->val[jptr*block_size*block_size+1];
                                   out->val[iptr*block_size*block_size+3]=P->val[jptr*block_size*block_size+3];
                              }
                        }
                 }
            }
    } else if (block_size==3) {
           for (i=0;i<out->numRows;++i) {
                 for (iptr=out->pattern->ptr[i];iptr<out->pattern->ptr[i+1]; ++iptr) {
                       j=out->pattern->index[iptr];
                        /*This for later will be replaced by bsearch!!*/
                        for (jptr=P->pattern->ptr[j];jptr<P->pattern->ptr[j+1]; ++jptr) {
                              k=P->pattern->index[jptr];
                              if(k==i) {
                                   out->val[iptr*block_size*block_size]=P->val[jptr*block_size*block_size];
                                   out->val[iptr*block_size*block_size+1]=P->val[jptr*block_size*block_size+3];
                                   out->val[iptr*block_size*block_size+2]=P->val[jptr*block_size*block_size+6];
                                   out->val[iptr*block_size*block_size+3]=P->val[jptr*block_size*block_size+1];
                                   out->val[iptr*block_size*block_size+4]=P->val[jptr*block_size*block_size+4];
                                   out->val[iptr*block_size*block_size+5]=P->val[jptr*block_size*block_size+7];
                                   out->val[iptr*block_size*block_size+6]=P->val[jptr*block_size*block_size+2];
                                   out->val[iptr*block_size*block_size+7]=P->val[jptr*block_size*block_size+5];
                                   out->val[iptr*block_size*block_size+8]=P->val[jptr*block_size*block_size+8];
                              }
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
  double sum_all_neg1,sum_all_pos1,sum_strong_neg1,sum_strong_pos1,sum_all_neg2,sum_all_pos2,sum_strong_neg2,sum_strong_pos2,sum_all_neg3,sum_all_pos3,sum_strong_neg3,sum_strong_pos3;
  double a_ii1=0;
  double a_ii2=0;
  double a_ii3=0;
  
  dim_t n=A->numRows;
  dim_t F=W_FC->numRows;
  dim_t i,j,k;
  
  dim_t block_size=A->row_block_size;
  
  index_t iPtr,dptr=0;
  /*index_t *index, *where_p;*/
 
  alpha=TMPMEMALLOC(F*block_size,double);
  beta=TMPMEMALLOC(F*block_size,double);


  if (block_size==1) {
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
                  else {
                    a_ii+=sum_all_pos;
                    beta[k]=0;
                  }
                  alpha[k]=alpha[k]/(a_ii);
                  beta[k]=beta[k]/(a_ii);
                k++;
             /*printf("Got in row=%d, alpha[%d]=%e, beta[%d]=%e, a_den=%e, b_den=%e \n",i,k-1,alpha[k-1],k-1,beta[k-1],alpha_den[k-1],beta_den[k-1]);*/
            }
         }
  } else if (block_size==2) {
            k=0;
        for (i = 0; i < n; ++i) {
            if(mis_marker[i]) {
                  alpha[k*block_size]=0;
                  alpha[k*block_size+1]=0;
                  beta[k*block_size]=0;
                  beta[k*block_size+1]=0;
                  sum_all_neg1=0;
                  sum_all_pos1=0;
                  sum_strong_neg1=0;
                  sum_strong_pos1=0;
                  sum_all_neg2=0;
                  sum_all_pos2=0;
                  sum_strong_neg2=0;
                  sum_strong_pos2=0;
                  /*compute beta_i=sum_{j\inN_i}a^{+}_ij and alpha_i=sum_{j\inN_i}a^{-}_ij and denominators for both alpha_i and beta_i*/
                  for (iPtr=A->pattern->ptr[i];iPtr<A->pattern->ptr[i + 1]; ++iPtr) {
                       j=A->pattern->index[iPtr];
                       if(j!=i) {
                              if(A->val[iPtr*block_size*block_size]<0) {
                                sum_all_neg1+=A->val[iPtr*block_size*block_size];
                                if(!mis_marker[j]) {
                                  sum_strong_neg1+=A->val[iPtr*block_size*block_size];
                                }
                              }
                              else if(A->val[iPtr*block_size*block_size]>0) {
                                sum_all_pos1+=A->val[iPtr*block_size*block_size];
                                if(!mis_marker[j]) {
                                  sum_strong_pos1+=A->val[iPtr*block_size*block_size];
                                }
                              }
                              if(A->val[iPtr*block_size*block_size+3]<0) {
                                sum_all_neg2+=A->val[iPtr*block_size*block_size+3];
                                if(!mis_marker[j]) {
                                  sum_strong_neg2+=A->val[iPtr*block_size*block_size+3];
                                }
                              } else if(A->val[iPtr*block_size*block_size+3]>0) {
                                sum_all_pos2+=A->val[iPtr*block_size*block_size+3];
                                if(!mis_marker[j]) {
                                  sum_strong_pos2+=A->val[iPtr*block_size*block_size+3];
                                }
                              }
                              
                       }
                       else {
                            a_ii1=A->val[iPtr*block_size*block_size];
                            a_ii2=A->val[iPtr*block_size*block_size+3];
                            dptr=iPtr;
                       }
                   }
      
                  if(sum_strong_neg1!=0) {
                       alpha[k*block_size]=sum_all_neg1/(sum_strong_neg1);
                  }
                  if(sum_strong_neg2!=0) {
                       alpha[k*block_size+1]=sum_all_neg2/(sum_strong_neg2);
                  }
                  
                  if(sum_strong_pos1!=0) {
                       beta[k*block_size]=sum_all_pos1/(sum_strong_pos1);
                  }
                  else {
                    a_ii1+=sum_all_pos1;
                    beta[k*block_size]=0;
                  }
                  if(sum_strong_pos2!=0) {
                       beta[k*block_size+1]=sum_all_pos2/(sum_strong_pos2);
                  }
                  else {
                    a_ii2+=sum_all_pos2;
                    beta[k*block_size+1]=0;
                  }
                  
                  alpha[k*block_size]=alpha[k*block_size]/(a_ii1);
                  alpha[k*block_size+1]=alpha[k*block_size+1]/(a_ii2);
                  beta[k*block_size]=beta[k*block_size]/(a_ii1);
                  beta[k*block_size+1]=beta[k*block_size+1]/(a_ii2);
                k++;
             /*printf("Got in row=%d, alpha[%d]=%e, beta[%d]=%e, a_den=%e, b_den=%e \n",i,k-1,alpha[k-1],k-1,beta[k-1],alpha_den[k-1],beta_den[k-1]);*/
            }
         }
  } else if (block_size==3) {
            k=0;
        for (i = 0; i < n; ++i) {
            if(mis_marker[i]) {
                  alpha[k*block_size]=0;
                  alpha[k*block_size+1]=0;
                  alpha[k*block_size+2]=0;
                  beta[k*block_size]=0;
                  beta[k*block_size+1]=0;
                  beta[k*block_size+2]=0;
                  sum_all_neg1=0;
                  sum_all_pos1=0;
                  sum_strong_neg1=0;
                  sum_strong_pos1=0;
                  sum_all_neg2=0;
                  sum_all_pos2=0;
                  sum_strong_neg2=0;
                  sum_strong_pos2=0;
                  sum_all_neg3=0;
                  sum_all_pos3=0;
                  sum_strong_neg3=0;
                  sum_strong_pos3=0;
                  /*compute beta_i=sum_{j\inN_i}a^{+}_ij and alpha_i=sum_{j\inN_i}a^{-}_ij and denominators for both alpha_i and beta_i*/
                  for (iPtr=A->pattern->ptr[i];iPtr<A->pattern->ptr[i + 1]; ++iPtr) {
                       j=A->pattern->index[iPtr];
                       if(j!=i) {
                              if(A->val[iPtr*block_size*block_size]<0) {
                                sum_all_neg1+=A->val[iPtr*block_size*block_size];
                                if(!mis_marker[j]) {
                                  sum_strong_neg1+=A->val[iPtr*block_size*block_size];
                                }
                              }
                              else if(A->val[iPtr*block_size*block_size]>0) {
                                sum_all_pos1+=A->val[iPtr*block_size*block_size];
                                if(!mis_marker[j]) {
                                  sum_strong_pos1+=A->val[iPtr*block_size*block_size];
                                }
                              }
                              if(A->val[iPtr*block_size*block_size+4]<0) {
                                sum_all_neg2+=A->val[iPtr*block_size*block_size+4];
                                if(!mis_marker[j]) {
                                  sum_strong_neg2+=A->val[iPtr*block_size*block_size+4];
                                }
                              } else if(A->val[iPtr*block_size*block_size+4]>0) {
                                sum_all_pos2+=A->val[iPtr*block_size*block_size+4];
                                if(!mis_marker[j]) {
                                  sum_strong_pos2+=A->val[iPtr*block_size*block_size+4];
                                }
                              }
                              if(A->val[iPtr*block_size*block_size+8]<0) {
                                sum_all_neg3+=A->val[iPtr*block_size*block_size+8];
                                if(!mis_marker[j]) {
                                  sum_strong_neg3+=A->val[iPtr*block_size*block_size+8];
                                }
                              } else if(A->val[iPtr*block_size*block_size+8]>0) {
                                sum_all_pos3+=A->val[iPtr*block_size*block_size+8];
                                if(!mis_marker[j]) {
                                  sum_strong_pos3+=A->val[iPtr*block_size*block_size+8];
                                }
                              }
 
                              
                       }
                       else {
                            a_ii1=A->val[iPtr*block_size*block_size];
                            a_ii2=A->val[iPtr*block_size*block_size+4];
                            a_ii3=A->val[iPtr*block_size*block_size+8];
                            dptr=iPtr;
                       }
                   }
      
                  if(sum_strong_neg1!=0) {
                       alpha[k*block_size]=sum_all_neg1/(sum_strong_neg1);
                  }
                  if(sum_strong_neg2!=0) {
                       alpha[k*block_size+1]=sum_all_neg2/(sum_strong_neg2);
                  }
                  if(sum_strong_neg3!=0) {
                       alpha[k*block_size+2]=sum_all_neg3/(sum_strong_neg3);
                  }
                  
                  if(sum_strong_pos1!=0) {
                       beta[k*block_size]=sum_all_pos1/(sum_strong_pos1);
                  }
                  else {
                    a_ii1+=sum_all_pos1;
                    beta[k*block_size]=0;
                  }
                  if(sum_strong_pos2!=0) {
                       beta[k*block_size+1]=sum_all_pos2/(sum_strong_pos2);
                  }
                  else {
                    a_ii2+=sum_all_pos2;
                    beta[k*block_size+1]=0;
                  }
                  if(sum_strong_pos3!=0) {
                       beta[k*block_size+2]=sum_all_pos3/(sum_strong_pos3);
                  }
                  else {
                    a_ii3+=sum_all_pos3;
                    beta[k*block_size+2]=0;
                  }
                  
                  alpha[k*block_size]=alpha[k*block_size]/(a_ii1);
                  alpha[k*block_size+1]=alpha[k*block_size+1]/(a_ii2);
                  alpha[k*block_size+2]=alpha[k*block_size+2]/(a_ii3);
                  beta[k*block_size]=beta[k*block_size]/(a_ii1);
                  beta[k*block_size+1]=beta[k*block_size+1]/(a_ii2);
                  beta[k*block_size+2]=beta[k*block_size+2]/(a_ii3);
                k++;
             /*printf("Got in row=%d, alpha[%d]=%e, beta[%d]=%e, a_den=%e, b_den=%e \n",i,k-1,alpha[k-1],k-1,beta[k-1],alpha_den[k-1],beta_den[k-1]);*/
            }
         }
  }

  if (block_size==1) {
        #pragma omp parallel for private(i,iPtr) schedule(static)
        for (i = 0; i < W_FC->numRows; ++i) {
              for (iPtr=W_FC->pattern->ptr[i];iPtr<W_FC->pattern->ptr[i + 1]; ++iPtr) {
                     if(W_FC->val[iPtr]<0) {
                        W_FC->val[iPtr]=-alpha[i]*W_FC->val[iPtr];
                      }
                     else if(W_FC->val[iPtr]>0) {
                        W_FC->val[iPtr]=-beta[i]*W_FC->val[iPtr];
                     }
              }
        }
  } else if (block_size==2) {
        #pragma omp parallel for private(i,iPtr) schedule(static)
        for (i = 0; i < W_FC->numRows; ++i) {
              for (iPtr=W_FC->pattern->ptr[i];iPtr<W_FC->pattern->ptr[i + 1]; ++iPtr) {
                     if(W_FC->val[iPtr*block_size*block_size]<0) {
                        W_FC->val[iPtr*block_size*block_size]=-alpha[i*block_size]*W_FC->val[iPtr*block_size*block_size];
                      }
                     else if(W_FC->val[iPtr*block_size*block_size]>0) {
                        W_FC->val[iPtr*block_size*block_size]=-beta[i*block_size]*W_FC->val[iPtr*block_size*block_size];
                     }
                     if(W_FC->val[iPtr*block_size*block_size+3]<0) {
                        W_FC->val[iPtr*block_size*block_size+3]=-alpha[i*block_size+1]*W_FC->val[iPtr*block_size*block_size+3];
                      }
                     else if(W_FC->val[iPtr*block_size*block_size+3]>0) {
                        W_FC->val[iPtr*block_size*block_size+3]=-beta[i*block_size+1]*W_FC->val[iPtr*block_size*block_size+3];
                     }
                     W_FC->val[iPtr*block_size*block_size+1]=0;
                     W_FC->val[iPtr*block_size*block_size+2]=0;
              }
        }
  } else if (block_size==3) {
       #pragma omp parallel for private(i,iPtr) schedule(static)
       for (i = 0; i < W_FC->numRows; ++i) {
             for (iPtr=W_FC->pattern->ptr[i];iPtr<W_FC->pattern->ptr[i + 1]; ++iPtr) {
                    if(W_FC->val[iPtr*block_size*block_size]<0) {
                       W_FC->val[iPtr*block_size*block_size]=-alpha[i*block_size]*W_FC->val[iPtr*block_size*block_size];
                     }
                    else if(W_FC->val[iPtr*block_size*block_size]>0) {
                       W_FC->val[iPtr*block_size*block_size]=-beta[i*block_size]*W_FC->val[iPtr*block_size*block_size];
                    }
                    if(W_FC->val[iPtr*block_size*block_size+4]<0) {
                       W_FC->val[iPtr*block_size*block_size+4]=-alpha[i*block_size+1]*W_FC->val[iPtr*block_size*block_size+4];
                     }
                    else if(W_FC->val[iPtr*block_size*block_size+4]>0) {
                       W_FC->val[iPtr*block_size*block_size+4]=-beta[i*block_size+1]*W_FC->val[iPtr*block_size*block_size+4];
                    }
                    if(W_FC->val[iPtr*block_size*block_size+8]<0) {
                       W_FC->val[iPtr*block_size*block_size+8]=-alpha[i*block_size+2]*W_FC->val[iPtr*block_size*block_size+8];
                     }
                    else if(W_FC->val[iPtr*block_size*block_size+8]>0) {
                       W_FC->val[iPtr*block_size*block_size+8]=-beta[i*block_size+2]*W_FC->val[iPtr*block_size*block_size+8];
                    }
                    W_FC->val[iPtr*block_size*block_size+1]=0;
                    W_FC->val[iPtr*block_size*block_size+2]=0;
                    W_FC->val[iPtr*block_size*block_size+3]=0;
                    W_FC->val[iPtr*block_size*block_size+5]=0;
                    W_FC->val[iPtr*block_size*block_size+6]=0;
                    W_FC->val[iPtr*block_size*block_size+7]=0;
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
  double sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9,b_lj1=0,b_lj2=0,b_lj3=0,b_lj4=0,b_lj5=0,b_lj6=0,b_lj7=0,b_lj8=0,b_lj9=0;
  
  dim_t block_size=A->row_block_size;
  
  double time0=0;
  bool_t verbose=0;
  
  time0=Paso_timer();
  
  outpattern=Paso_Pattern_multiply(PATTERN_FORMAT_DEFAULT,A->pattern,B->pattern);
  out=Paso_SparseMatrix_alloc(A->type, outpattern, block_size, block_size, FALSE);
  
  time0=Paso_timer()-time0;
  if (verbose) fprintf(stdout,"timing: Paso_SparseMatrix_MatrixMatrix: Pattern creation: %e\n",time0);
  
  time0=Paso_timer();
  
  if(block_size==1) {
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
  } else if (block_size==2) {
          #pragma omp parallel for private(i,iptrC,j,sum1,sum2,sum3,sum4,iptrA,k,b_lj1,b_lj2,b_lj3,b_lj4,iptrB,l) schedule(static)
        for(i = 0; i < out->numRows; i++) {
          for(iptrC = out->pattern->ptr[i]; iptrC < out->pattern->ptr[i+1]; ++iptrC) {
            j = out->pattern->index[iptrC];
            sum1=0;
            sum2=0;
            sum3=0;
            sum4=0;
            for(iptrA = A->pattern->ptr[i]; iptrA < A->pattern->ptr[i+1]; ++iptrA) {
                  k=A->pattern->index[iptrA];
                  /*can be replaced by bsearch */
                  b_lj1=0;
                  b_lj2=0;
                  b_lj3=0;
                  b_lj4=0;
                  for(iptrB = B->pattern->ptr[k]; iptrB < B->pattern->ptr[k+1]; ++iptrB) {
                     l=B->pattern->index[iptrB];
                     if(l==j) {
                        b_lj1=B->val[iptrB*block_size*block_size];
                        b_lj2=B->val[iptrB*block_size*block_size+1];
                        b_lj3=B->val[iptrB*block_size*block_size+2];
                        b_lj4=B->val[iptrB*block_size*block_size+3];
                        break;
                     }
                  }
                  sum1+=A->val[iptrA*block_size*block_size]*b_lj1+A->val[iptrA*block_size*block_size+1]*b_lj3;
                  sum2+=A->val[iptrA*block_size*block_size]*b_lj2+A->val[iptrA*block_size*block_size+1]*b_lj4;
                  sum3+=A->val[iptrA*block_size*block_size+2]*b_lj1+A->val[iptrA*block_size*block_size+3]*b_lj3;
                  sum4+=A->val[iptrA*block_size*block_size+2]*b_lj2+A->val[iptrA*block_size*block_size+3]*b_lj4;
            }
            out->val[iptrC*block_size*block_size]=sum1;
            out->val[iptrC*block_size*block_size+1]=sum2;
            out->val[iptrC*block_size*block_size+2]=sum3;
            out->val[iptrC*block_size*block_size+3]=sum4;
         }
        }
  } else if (block_size==3) {
          #pragma omp parallel for private(i,iptrC,j,sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9,iptrA,k,b_lj1,b_lj2,b_lj3,b_lj4,b_lj5,b_lj6,b_lj7,b_lj8,b_lj9,iptrB,l) schedule(static)
        for(i = 0; i < out->numRows; i++) {
          for(iptrC = out->pattern->ptr[i]; iptrC < out->pattern->ptr[i+1]; ++iptrC) {
            j = out->pattern->index[iptrC];
            sum1=0;
            sum2=0;
            sum3=0;
            sum4=0;
            sum5=0;
            sum6=0;
            sum7=0;
            sum8=0;
            sum9=0;
            for(iptrA = A->pattern->ptr[i]; iptrA < A->pattern->ptr[i+1]; ++iptrA) {
                  k=A->pattern->index[iptrA];
                  /*can be replaced by bsearch */
                  b_lj1=0;
                  b_lj2=0;
                  b_lj3=0;
                  b_lj4=0;
                  b_lj5=0;
                  b_lj6=0;
                  b_lj7=0;
                  b_lj8=0;
                  b_lj9=0;
                  for(iptrB = B->pattern->ptr[k]; iptrB < B->pattern->ptr[k+1]; ++iptrB) {
                     l=B->pattern->index[iptrB];
                     if(l==j) {
                        b_lj1=B->val[iptrB*block_size*block_size];
                        b_lj2=B->val[iptrB*block_size*block_size+1];
                        b_lj3=B->val[iptrB*block_size*block_size+2];
                        b_lj4=B->val[iptrB*block_size*block_size+3];
                        b_lj5=B->val[iptrB*block_size*block_size+4];
                        b_lj6=B->val[iptrB*block_size*block_size+5];
                        b_lj7=B->val[iptrB*block_size*block_size+6];
                        b_lj8=B->val[iptrB*block_size*block_size+7];
                        b_lj9=B->val[iptrB*block_size*block_size+8];
                        break;
                     }
                  }
                  sum1+=A->val[iptrA*block_size*block_size]*b_lj1+A->val[iptrA*block_size*block_size+1]*b_lj4+A->val[iptrA*block_size*block_size+2]*b_lj7;
                  sum2+=A->val[iptrA*block_size*block_size]*b_lj2+A->val[iptrA*block_size*block_size+1]*b_lj5+A->val[iptrA*block_size*block_size+2]*b_lj8;
                  sum3+=A->val[iptrA*block_size*block_size]*b_lj3+A->val[iptrA*block_size*block_size+1]*b_lj6+A->val[iptrA*block_size*block_size+2]*b_lj9;
                  sum4+=A->val[iptrA*block_size*block_size+3]*b_lj1+A->val[iptrA*block_size*block_size+4]*b_lj4+A->val[iptrA*block_size*block_size+5]*b_lj7;
                  sum5+=A->val[iptrA*block_size*block_size+3]*b_lj2+A->val[iptrA*block_size*block_size+4]*b_lj5+A->val[iptrA*block_size*block_size+5]*b_lj8;
                  sum6+=A->val[iptrA*block_size*block_size+3]*b_lj3+A->val[iptrA*block_size*block_size+4]*b_lj6+A->val[iptrA*block_size*block_size+5]*b_lj9;
                  sum7+=A->val[iptrA*block_size*block_size+6]*b_lj1+A->val[iptrA*block_size*block_size+7]*b_lj4+A->val[iptrA*block_size*block_size+8]*b_lj7;
                  sum8+=A->val[iptrA*block_size*block_size+6]*b_lj2+A->val[iptrA*block_size*block_size+7]*b_lj5+A->val[iptrA*block_size*block_size+8]*b_lj8;
                  sum9+=A->val[iptrA*block_size*block_size+6]*b_lj3+A->val[iptrA*block_size*block_size+7]*b_lj6+A->val[iptrA*block_size*block_size+8]*b_lj9;
            }
            out->val[iptrC*block_size*block_size]=sum1;
            out->val[iptrC*block_size*block_size+1]=sum2;
            out->val[iptrC*block_size*block_size+2]=sum3;
            out->val[iptrC*block_size*block_size+3]=sum4;
            out->val[iptrC*block_size*block_size+4]=sum5;
            out->val[iptrC*block_size*block_size+5]=sum6;
            out->val[iptrC*block_size*block_size+6]=sum7;
            out->val[iptrC*block_size*block_size+7]=sum8;
            out->val[iptrC*block_size*block_size+8]=sum9;
         }
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


Paso_SparseMatrix* Paso_SparseMatrix_RemovePositiveOffdiagonals(Paso_SparseMatrix* P){
 
  Paso_Pattern *outpattern=NULL;
  Paso_SparseMatrix *out=NULL;
  
  Paso_IndexList* index_list=NULL;
  dim_t n=P->numRows;
  dim_t i,j,k=0;
  index_t iptr,jptr;

  index_list=TMPMEMALLOC(n,Paso_IndexList);
   if (! Paso_checkPtr(index_list)) {
        #pragma omp parallel for private(i) schedule(static)
        for(i=0;i<n;++i) {
             index_list[i].extension=NULL;
             index_list[i].n=0;
        }
    }
  

    for (i=0;i<n;++i) {
          for (iptr=P->pattern->ptr[i];iptr<P->pattern->ptr[i+1]; ++iptr) {
             j=P->pattern->index[iptr];
             if(P->val[iptr]<0) {
             Paso_IndexList_insertIndex(&(index_list[j]),i);
             }
             if(i==j) {
              Paso_IndexList_insertIndex(&(index_list[j]),i);
             }
        }
    }
   
    outpattern=Paso_IndexList_createPattern(0, n,index_list,0,n,0);
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
        for(i=0;i<n;++i) Paso_IndexList_free(index_list[i].extension);
     }
    TMPMEMFREE(index_list);
    Paso_Pattern_free(outpattern);
    return out;
}


