
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


/**********************************************************************/

/* Paso: Pattern: Paso_Pattern_coupling 

   searches for a maximal independent set MIS in the matrix pattern 
   vertices in the maximal independent set are marked in mis_marker
   nodes to be considered are marked by -1 on the input in mis_marker

*/
/**********************************************************************/

/* Copyrights by ACcESS Australia 2003,2004,2005              */
/* Author: artak@uq.edu.au                                */

/**************************************************************/

#include "PasoUtil.h"
#include "Pattern_coupling.h"
#include <limits.h>


/***************************************************************/
 
#define IS_AVAILABLE -1
#define IS_IN_F -3   /* in F (strong) */
#define IS_IN_C -4  /* in C (weak) */ 

void Paso_Pattern_YS(Paso_SparseMatrix* A, index_t* mis_marker, double threshold) {

  dim_t i,j;
  /*double sum;*/
  index_t iptr,*index,*where_p,*diagptr;
  bool_t passed=FALSE;
  dim_t n=A->numRows;
  diagptr=MEMALLOC(n,index_t);

  if (A->pattern->type & PATTERN_FORMAT_SYM) {
    Paso_setError(TYPE_ERROR,"Paso_Pattern_coup: symmetric matrix pattern is not supported yet");
    return;
  }
   
   #pragma omp parallel for private(i) schedule(static)
   for (i=0;i<n;++i)
        if(mis_marker[i]==IS_AVAILABLE)
                    mis_marker[i]=IS_IN_C;

    /*#pragma omp parallel for private(i,index,where_p) schedule(static)*/
    for (i=0;i<n;++i) {
         diagptr[i]=A->pattern->ptr[i];
         index=&(A->pattern->index[A->pattern->ptr[i]]);
         where_p=(index_t*)bsearch(&i,
                                index,
                                A->pattern->ptr[i + 1]-A->pattern->ptr[i],
                                sizeof(index_t),
                                Paso_comparIndex);
        if (where_p==NULL) {
            Paso_setError(VALUE_ERROR, "Paso_Pattern_coup: main diagonal element missing.");
        } else {
                diagptr[i]+=(index_t)(where_p-index);
        }
    }
    


    /*This loop cannot be parallelized, as order matters here.*/ 
    for (i=0;i<n;++i) {
      if (mis_marker[i]==IS_IN_C) {
        for (iptr=A->pattern->ptr[i];iptr<A->pattern->ptr[i+1]; ++iptr) {
             j=A->pattern->index[iptr];
             if (j!=i && ABS(A->val[iptr])>=threshold*ABS(A->val[diagptr[i]])) {
                mis_marker[j]=IS_IN_F;
             }
        }
      }
    }
    
    
     
      /*This loop cannot be parallelized, as order matters here.*/ 
    for (i=0;i<n;i++) {
        if (mis_marker[i]==IS_IN_F) {
           passed=TRUE;
           for (iptr=A->pattern->ptr[i];iptr<A->pattern->ptr[i+1]; ++iptr) {
              j=A->pattern->index[iptr];
              if (mis_marker[j]==IS_IN_C) {
                if ((A->val[iptr]/A->val[diagptr[i]])>=-threshold) {
                    passed=TRUE;
                }
                else {
                    passed=FALSE;
                    break;
                }
              } 
           }
           if (passed) mis_marker[i]=IS_IN_C;
        }
    }

     /* swap to TRUE/FALSE in mis_marker */
     #pragma omp parallel for private(i) schedule(static)
     for (i=0;i<n;i++) mis_marker[i]=(mis_marker[i]==IS_IN_F);
     
     MEMFREE(diagptr);
}


/*
 * Ruge-Stueben strength of connection mask.
 *
 */
void Paso_Pattern_RS(Paso_SparseMatrix* A, index_t* mis_marker, double theta)
{
  dim_t i,n,j;
  index_t iptr;
  double threshold,max_offdiagonal;
  
  Paso_Pattern *out=NULL;
  
  Paso_IndexList* index_list=NULL;

 index_list=TMPMEMALLOC(A->pattern->numOutput,Paso_IndexList);
   if (! Paso_checkPtr(index_list)) {
        #pragma omp parallel for private(i) schedule(static)
        for(i=0;i<A->pattern->numOutput;++i) {
             index_list[i].extension=NULL;
             index_list[i].n=0;
        }
    }
  
  
  n=A->numRows;
  if (A->pattern->type & PATTERN_FORMAT_SYM) {
    Paso_setError(TYPE_ERROR,"Paso_Pattern_RS: symmetric matrix pattern is not supported yet");
    return;
  }
    /*#pragma omp parallel for private(i,iptr,max_offdiagonal,threshold,j) schedule(static)*/
    for (i=0;i<n;++i) {
      if(mis_marker[i]==IS_AVAILABLE) {
        max_offdiagonal = DBL_MIN;
        for (iptr=A->pattern->ptr[i];iptr<A->pattern->ptr[i+1]; ++iptr) {
            if(A->pattern->index[iptr] != i){
                max_offdiagonal = MAX(max_offdiagonal,-A->val[iptr]);
            }
        }
        
        threshold = theta*max_offdiagonal;
        for (iptr=A->pattern->ptr[i];iptr<A->pattern->ptr[i+1]; ++iptr) {
            j=A->pattern->index[iptr];
            if((-A->val[iptr])>=threshold) {
                Paso_IndexList_insertIndex(&(index_list[i]),j);
                Paso_IndexList_insertIndex(&(index_list[j]),i);
            }
        }
       }
      }
    
   
    out=Paso_IndexList_createPattern(0, A->pattern->numOutput,index_list,0,A->pattern->numInput,0);
    
     /* clean up */
   if (index_list!=NULL) {
        #pragma omp parallel for private(i) schedule(static)
        for(i=0;i<A->pattern->numOutput;++i) Paso_IndexList_free(index_list[i].extension);
     }
    TMPMEMFREE(index_list);

    /*Paso_Pattern_mis(out,mis_marker);*/
    Paso_Pattern_greedy(out,mis_marker);
    Paso_Pattern_free(out);
}

void Paso_Pattern_Aggregiation(Paso_SparseMatrix* A, index_t* mis_marker, double theta)
{
  dim_t i,j,n;
  index_t iptr;
  double diag,eps_Aii,val;
  double* diags;


  Paso_Pattern *out=NULL;
  Paso_IndexList* index_list=NULL;

  n=A->numRows;  
  diags=MEMALLOC(n,double);

  index_list=TMPMEMALLOC(A->pattern->numOutput,Paso_IndexList);
   if (! Paso_checkPtr(index_list)) {
        #pragma omp parallel for private(i) schedule(static)
        for(i=0;i<A->pattern->numOutput;++i) {
             index_list[i].extension=NULL;
             index_list[i].n=0;
        }
    }
    
  if (A->pattern->type & PATTERN_FORMAT_SYM) {
    Paso_setError(TYPE_ERROR,"Paso_Pattern_Aggregiation: symmetric matrix pattern is not supported yet");
    return;
  }


    #pragma omp parallel for private(i,iptr,diag) schedule(static)
      for (i=0;i<n;++i) {
        diag = 0;
        for (iptr=A->pattern->ptr[i];iptr<A->pattern->ptr[i+1]; ++iptr) {
            if(A->pattern->index[iptr] == i){
                diag+=A->val[iptr];
            }
        }
        diags[i]=ABS(diag);
      }


    #pragma omp parallel for private(i,iptr,j,val,eps_Aii) schedule(static)
     for (i=0;i<n;++i) {
       if (mis_marker[i]==IS_AVAILABLE) {
        eps_Aii = theta*theta*diags[i];
        for (iptr=A->pattern->ptr[i];iptr<A->pattern->ptr[i+1]; ++iptr) {
            j=A->pattern->index[iptr];
            val=A->val[iptr];
              if((val*val)>=(eps_Aii*diags[j])) {
               Paso_IndexList_insertIndex(&(index_list[i]),j);
              }
        }
       }
     }

    out=Paso_IndexList_createPattern(0, A->pattern->numOutput,index_list,0,A->pattern->numInput,0);
    
     /* clean up */
    if (index_list!=NULL) {
        #pragma omp parallel for private(i) schedule(static)
        for(i=0;i<A->pattern->numOutput;++i) Paso_IndexList_free(index_list[i].extension);
     }

    TMPMEMFREE(index_list);
    MEMFREE(diags);
    
    
    /*Paso_Pattern_mis(out,mis_marker);*/
    Paso_Pattern_greedy(out,mis_marker);
    Paso_Pattern_free(out);

}

/* Greedy algorithm */
void Paso_Pattern_greedy(Paso_Pattern* pattern, index_t* mis_marker) {

  dim_t i,j;
  /*double sum;*/
  index_t iptr;
  bool_t passed=FALSE;
  dim_t n=pattern->numOutput;

  if (pattern->type & PATTERN_FORMAT_SYM) {
    Paso_setError(TYPE_ERROR,"Paso_Pattern_greedy: symmetric matrix pattern is not supported yet");
    return;
  }
   
   #pragma omp parallel for private(i) schedule(static)
   for (i=0;i<n;++i)
        if(mis_marker[i]==IS_AVAILABLE)
                    mis_marker[i]=IS_IN_C;


    for (i=0;i<n;++i) {
      if (mis_marker[i]==IS_IN_C) {
        for (iptr=pattern->ptr[i];iptr<pattern->ptr[i+1]; ++iptr) {
             j=pattern->index[iptr];
             mis_marker[j]=IS_IN_F;
        }
      }
    }
    
    
     
    for (i=0;i<n;i++) {
        if (mis_marker[i]==IS_IN_F) {
           passed=TRUE;
           for (iptr=pattern->ptr[i];iptr<pattern->ptr[i+1]; ++iptr) {
              j=pattern->index[iptr];
                if (mis_marker[j]==IS_IN_F) {
                    passed=TRUE;
                }
                else {
                    passed=FALSE;
                    break;
                }
              } 
           if (passed) mis_marker[i]=IS_IN_C;
           }
        }

     /* swap to TRUE/FALSE in mis_marker */
     #pragma omp parallel for private(i) schedule(static)
     for (i=0;i<n;i++) mis_marker[i]=(mis_marker[i]==IS_IN_F);
     
}


void Paso_Pattern_greedy_color(Paso_Pattern* pattern, index_t* mis_marker) {

  dim_t i,j;
  /*double sum;*/
  index_t iptr;
  index_t num_colors;
  index_t* colorOf;
  register index_t color;
  bool_t passed=FALSE;
  dim_t n=pattern->numOutput;

  
  colorOf=MEMALLOC(n,index_t);

  if (pattern->type & PATTERN_FORMAT_SYM) {
    Paso_setError(TYPE_ERROR,"Paso_Pattern_greedy: symmetric matrix pattern is not supported yet");
    return;
  }
   
   Paso_Pattern_color(pattern,&num_colors,colorOf);
   
   /* We do not need this loop if we set IS_IN_MIS=IS_AVAILABLE. */
   #pragma omp parallel for private(i) schedule(static)
   for (i=0;i<n;++i)
        if(mis_marker[i]==IS_AVAILABLE)
                    mis_marker[i]=IS_IN_F;

   #pragma omp barrier
   for (color=0;color<num_colors;++color) {
    #pragma omp parallel for schedule(static) private(i,iptr,j)
    for (i=0;i<n;++i) {
     if (colorOf[i]==color) {  
      if (mis_marker[i]==IS_IN_F) {
        for (iptr=pattern->ptr[i];iptr<pattern->ptr[i+1]; ++iptr) {
             j=pattern->index[iptr];
             if (colorOf[j]<color)
              mis_marker[j]=IS_IN_C;
        }
      }
     }
    }
   }
    
    
   #pragma omp barrier
   for (color=0;color<num_colors;++color) {
   #pragma omp parallel for schedule(static) private(i,iptr,j) 
    for (i=0;i<n;i++) {
      if (colorOf[i]==color) {  
        if (mis_marker[i]==IS_IN_C) {
           passed=TRUE;
           for (iptr=pattern->ptr[i];iptr<pattern->ptr[i+1]; ++iptr) {
              j=pattern->index[iptr];
               if (colorOf[j]<color && passed) {
                if (mis_marker[j]==IS_IN_C) {
                    passed=TRUE;
                }
                else {
                    passed=FALSE;
                    /*break;*/
                }
              }
           }
           if (passed) mis_marker[i]=IS_IN_F;
           }
           
        }
    }
   }

     /* swap to TRUE/FALSE in mis_marker */
     #pragma omp parallel for private(i) schedule(static)
     for (i=0;i<n;i++) mis_marker[i]=(mis_marker[i]==IS_IN_F);
    
    MEMFREE(colorOf); 
}

/*For testing */
void Paso_Pattern_greedy_diag(Paso_SparseMatrix* A, index_t* mis_marker, double threshold) {

  dim_t i,j=0,k;
  double *theta;
  index_t iptr;
  dim_t n=A->numRows;
  double rsum,diag=0;
  index_t *AvADJ;
  theta=MEMALLOC(n,double);
  AvADJ=MEMALLOC(n,index_t);


  

  if (A->pattern->type & PATTERN_FORMAT_SYM) {
    Paso_setError(TYPE_ERROR,"Paso_Pattern_coup: symmetric matrix pattern is not supported yet");
    return;
  }
   

    #pragma omp parallel for private(i,iptr,j,rsum) schedule(static) 
    for (i=0;i<n;++i) {
        rsum=0;
        for (iptr=A->pattern->ptr[i];iptr<A->pattern->ptr[i+1]; ++iptr) {
            j=A->pattern->index[iptr];
            if(j!=i) {
              rsum+=ABS(A->val[iptr]);    
            }
            else {
                diag=ABS(A->val[iptr]);
            }
        }
        theta[i]=diag/rsum;
        if(theta[i]>threshold) {
            mis_marker[i]=IS_IN_F;
        }
    }
    
    while (Paso_Util_isAny(n,mis_marker,IS_AVAILABLE)) {
         k=0;
         #pragma omp parallel for private(i,j,k) schedule(static) 
         for (i=0;i<n;++i) {
           if(mis_marker[i]==IS_AVAILABLE) {
                if(k==0) {
                    j=i;
                    k++;
                }
                if(theta[j]>theta[i]) {
                    j=i;
                }
            }
         }
         mis_marker[j]=IS_IN_C;
         
         for (iptr=A->pattern->ptr[j];iptr<A->pattern->ptr[j+1]; ++iptr) {
            k=A->pattern->index[iptr];
            if(mis_marker[k]==IS_AVAILABLE) {
               AvADJ[k]=1; 
            }
            else {
                AvADJ[k]=-1;
            }
            
         }
            
         #pragma omp parallel for private(i,iptr,j,rsum) schedule(static) 
        for (i=0;i<n;++i) {
            if(AvADJ[i]) {
                rsum=0;
                for (iptr=A->pattern->ptr[i];iptr<A->pattern->ptr[i+1]; ++iptr) {
                    k=A->pattern->index[iptr];
                    if(k!=i && mis_marker[k]!=IS_IN_C ) {
                      rsum+=ABS(A->val[iptr]);    
                    }
                    if(j==i) {
                        diag=ABS(A->val[iptr]);
                    }
                }
                theta[i]=diag/rsum;
                if(theta[i]>threshold) {
                   mis_marker[i]=IS_IN_F;
                }
            }
        }
         
        
    }

     /* swap to TRUE/FALSE in mis_marker */
     #pragma omp parallel for private(i) schedule(static)
     for (i=0;i<n;i++) mis_marker[i]=(mis_marker[i]==IS_IN_F);
     
     MEMFREE(AvADJ);
     MEMFREE(theta);
}

#undef IS_AVAILABLE 
#undef IS_IN_F 
#undef IS_IN_C










