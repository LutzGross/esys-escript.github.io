
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


/*************************************************************

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
/**************************************************************

  Author: artak@uq.edu.au, l.gross@uq.edu.au                   

**************************************************************/

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

