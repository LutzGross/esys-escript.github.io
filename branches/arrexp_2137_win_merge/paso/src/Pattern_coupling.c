
/*******************************************************
*
* Copyright (c) 2003-2008 by University of Queensland
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


/***************************************************************/
 
#define IS_AVAILABLE -1
#define IS_IN_SET -3
#define IS_REMOVED -4

void Paso_Pattern_coup(Paso_SparseMatrix* A, index_t* mis_marker, double threshold) {

  dim_t i,j;
  /*double threshold=0.05;*/
  double sum;
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
                    mis_marker[i]=IS_IN_SET;

    #pragma omp parallel for private(i,index,where_p) schedule(static) 
    for (i=0;i<n;++i) {
         diagptr[i]=A->pattern->ptr[i];
         index=&(A->pattern->index[diagptr[i]]);
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
      if (mis_marker[i]==IS_IN_SET) {
        for (iptr=A->pattern->ptr[i];iptr<A->pattern->ptr[i+1]; ++iptr) {
             j=A->pattern->index[iptr];
             if (j!=i && ABS(A->val[iptr])>=threshold*ABS(A->val[diagptr[i]])) {
                mis_marker[j]=IS_REMOVED;
             }
        }
      }
    }
    
    
     
      /*This loop cannot be parallelized, as order matters here.*/ 
    for (i=0;i<n;i++) {
        if (mis_marker[i]==IS_REMOVED) {
           for (iptr=A->pattern->ptr[i];iptr<A->pattern->ptr[i+1]; ++iptr) {
              j=A->pattern->index[iptr];
              if (mis_marker[j]==IS_IN_SET) {
                if ((A->val[iptr]/A->val[diagptr[i]])>=-threshold) {
                    passed=TRUE;
                }
                else {
                    passed=FALSE;
                    break;
                }
              }
            }
            if (passed) {
               mis_marker[i]=IS_IN_SET;
               passed=FALSE;
            }
        }
    }

    /* This check is to make sure we dont get some nusty rows which were not removed durring coarsing process.*/
    /* TODO: we have to mechanism that this does not happend at all, and get rid of this 'If'. */
     /*This loop cannot be parallelized, as order matters here.*/ 
    for (i=0;i<n;i++) {
        if (mis_marker[i]==IS_REMOVED) {
           sum=0;
           for (iptr=A->pattern->ptr[i];iptr<A->pattern->ptr[i+1]; ++iptr) {
             j=A->pattern->index[iptr];
             if (mis_marker[j]==IS_REMOVED)
                sum+=A->val[iptr];
           }
           if(ABS(sum)<1.e-12)
             mis_marker[i]=IS_IN_SET;
        }
    }


     /* swap to TRUE/FALSE in mis_marker */
     #pragma omp parallel for private(i) schedule(static)
     for (i=0;i<n;i++) mis_marker[i]=(mis_marker[i]!=IS_IN_SET);
     
     MEMFREE(diagptr);
}

/*
 *
 * Return a strength of connection mask using the classical 
 * strength of connection measure by Ruge and Stuben.
 *
 * Specifically, an off-diagonal entry A[i.j] is a strong 
 * connection if:
 *  
 *      -A[i,j] >= theta * max( -A[i,k] )   where k != i
 * 
 * Otherwise, the connection is weak.
 *  
 */  
void Paso_Pattern_RS(Paso_SparseMatrix* A, index_t* mis_marker, double theta) 
{
  dim_t i;
  index_t iptr;
  double threshold,max_offdiagonal;
  dim_t n=A->numRows;
  if (A->pattern->type & PATTERN_FORMAT_SYM) {
    Paso_setError(TYPE_ERROR,"Paso_Pattern_RS: symmetric matrix pattern is not supported yet");
    return;
  }

    #pragma omp parallel for private(i,iptr,max_offdiagonal,threshold) schedule(static) 
    for (i=0;i<n;++i) {
      if(mis_marker[i]==IS_AVAILABLE) {
        /*min_offdiagonal = A->val[A->pattern->ptr[i]-index_offset];*/
        max_offdiagonal = -1.e+15;
        for (iptr=A->pattern->ptr[i];iptr<A->pattern->ptr[i+1]; ++iptr) {
            if(A->pattern->index[iptr] != i){
                max_offdiagonal = MAX(max_offdiagonal,-(A->val[iptr]));
            }
        }

        threshold = theta*max_offdiagonal;
        for (iptr=A->pattern->ptr[i];iptr<A->pattern->ptr[i+1]; ++iptr) {
            if(-(A->val[iptr])<threshold && A->pattern->index[iptr]!=i){
                    A->val[iptr]=0.;
            }
        }
       }
      }

    Paso_Pattern_coup(A,mis_marker,0.05);
}




void Paso_Pattern_Aggregiation(Paso_SparseMatrix* A, index_t* mis_marker, double theta) 
{
  dim_t i,j;
  index_t iptr;
  double diag,eps_Aii,val;
  dim_t n=A->numRows;
  double* diags=MEMALLOC(n,double);
  
  if (A->pattern->type & PATTERN_FORMAT_SYM) {
    Paso_setError(TYPE_ERROR,"Paso_Pattern_RS: symmetric matrix pattern is not supported yet");
    return;
  }
    

    #pragma omp parallel for private(i,iptr,diag) schedule(static) 
      for (i=0;i<n;++i) {
        diag = 0;
        for (iptr=A->pattern->ptr[i];iptr<A->pattern->ptr[i+1]; ++iptr) {
            if(A->pattern->index[iptr] != i){
                diag+=A->val[iptr];
            }
        }
        diags[i]=ABS(diag);
      }
    
   
    #pragma omp parallel for private(i,iptr,j,val,eps_Aii) schedule(static) 
     for (i=0;i<n;++i) {
       if (mis_marker[i]==IS_AVAILABLE) {
        eps_Aii = theta*theta*diags[i];
        val=0.; 
        for (iptr=A->pattern->ptr[i];iptr<A->pattern->ptr[i+1]; ++iptr) {
            j=A->pattern->index[iptr];
            val=A->val[iptr];
            if(j!= i) {
              if(val*val>=eps_Aii * diags[j]) {
               A->val[iptr]=val;
              }
              else {
               A->val[iptr]=0.;
              }
            }
        }
       }
     }
    
    Paso_Pattern_coup(A,mis_marker,0.05);
     MEMFREE(diags);
}


#undef IS_AVAILABLE 
#undef IS_IN_SET 
#undef IS_REMOVED

void Paso_Pattern_color1(Paso_SparseMatrix* A, index_t* num_colors, index_t* colorOf) {
  index_t out=0, *mis_marker=NULL,i;
  dim_t n=A->numRows;
  mis_marker=TMPMEMALLOC(n,index_t);
  if ( !Paso_checkPtr(mis_marker) ) {
    /* get coloring */
    #pragma omp parallel for private(i) schedule(static)
    for (i = 0; i < n; ++i) {
        colorOf[i]=-1;
        mis_marker[i]=-1;
    }

    while (Paso_Util_isAny(n,colorOf,-1) && Paso_noError()) {
       Paso_Pattern_coup(A,mis_marker,0.05);

       #pragma omp parallel for private(i) schedule(static)
       for (i = 0; i < n; ++i) {
        if (mis_marker[i]) colorOf[i]=out;
        mis_marker[i]=colorOf[i];
       }
       ++out;
    }
    TMPMEMFREE(mis_marker);
    *num_colors=out;
  }
}
