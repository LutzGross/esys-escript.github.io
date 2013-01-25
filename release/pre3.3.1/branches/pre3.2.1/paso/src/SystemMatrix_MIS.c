
#include "SystemMatrix_MIS.h"


#define MISIN 0
#define MISOUT 100
/* This macro is less complicated than it could have been because it is only used on available nodes */
#define ISLESS(x,y) (x<y)
#define ISAVAILABLE(x) ((x!=MISIN) && (x!=MISOUT))
#define IMAX(x,y) (x>y?x:y)

#define MISSTRING(x) ((x==MISIN)?"IN":((x==MISOUT)?"OUT":"UNKNOWN"))


/* returns the nodes in this system matrix which connect to things outside the main block
   the reference parameter count will be set to the number of nodes in the return value.
   Cleanup of return value is the callers responsibility.
   
   This routine assumes that the System matrix is in default format with a CSR col_coupleBlock 
*/
index_t* Paso_SparseMatrix_getBorderNodes(Paso_SystemMatrix* A, index_t* count) {
   const index_t MAXNEIGHBOURS=A->col_coupleBlock->len;
   index_t* border=MEMALLOC(MAXNEIGHBOURS, index_t);	
   index_t len=0;
   int i;
   /* A node is in the border if it has a non-empty row in the coupling block */
   #pragma omp parallel for schedule(static) private(i) 
   for (i=0;i<A->col_coupleBlock->pattern->numOutput;++i) {
	if (A->col_coupleBlock->pattern->ptr[i]!=A->col_coupleBlock->pattern->ptr[i+1]) {	/* non-empty row */
	   #pragma omp critical   
	   {
		border[len++]=i; 
	   }  
	}  
   } 
   *count=len;
   return border;
}


/* takes in a list of border nodes and weights and computes the MIS for all border ndoes
   on all ranks. It does this one node at a time. Later this should use rank colouring
   to do a number of ranks at once.
*/
void Paso_SystemMatrix_CalcBorderMIS(Paso_SystemMatrix* A, index_t* border, index_t bordercount, double* weights, index_t n) {
    index_t i=0;
    index_t j=0;
    index_t k=0;
    int mpi_iam = 0;	/* rank within world */
    int mpi_num = 1;	/* size of the world */
    #ifdef ESYS_MPI
    double *remote_values=NULL;
    #endif
    
    if (A->type!=MATRIX_FORMAT_DEFAULT) {		/* We only support CSR matricies here */
        Esys_setError(TYPE_ERROR,"Paso_SystemMatrix_CalcBorderMIS: Symmetric matrix patterns are not supported.");      
    }
    
    #ifdef ESYS_MPI
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_iam);
    #endif

    for (i=0;i<mpi_num;++i) {
	if (i==mpi_iam) {	/* mark my border in preparation for sending to other ranks */
	    for (j=0;j<bordercount;++j) {
		index_t bnode=border[j];
		if (ISAVAILABLE(weights[bnode])) {
		    Paso_Pattern* p=A->mainBlock->pattern;
		    weights[bnode]=MISIN;
		    /* Now walk the neighbours and mark them unavailable */
		    for (k=p->ptr[bnode];k<p->ptr[bnode+1];++k) {	/* Walk along the row */
			if (p->index[k]!=bnode) {	/* ignore diagonal link to self */
			    weights[p->index[k]]=MISOUT;	/* node can't be in the set */
		        }     
		    } /* walk row */
	        }     /* if avail */
	     }   /* walk border */
	} /* if not me */
	#ifdef ESYS_MPI
	/* Now we set the weights using the coupler */
  	remote_values=NULL;
	  /* start exchange */
	Paso_SystemMatrix_startCollect(A,weights);
	  /* finish exchange */
	remote_values=Paso_SystemMatrix_finishCollect(A);
	if (i!=mpi_iam) {
	    /* walk over nodes and see if they have a marked neighbour */
	    for (j=0;j<A->col_coupleBlock->pattern->numOutput;++j) {
	        for (k=A->col_coupleBlock->pattern->ptr[j];k<A->col_coupleBlock->pattern->ptr[j+1];++k) {
		    index_t con=A->col_coupleBlock->pattern->index[k];
		    if (remote_values[con]==MISIN) {
		        weights[j]=MISOUT;
			break;		/* done with this row */
		    }
		}
	    }
	}
	#endif
    }   /* End loop over ranks */
}


/* used to generate pseudo random numbers: */

static double Paso_Pattern_mis_seed=.4142135623730951;

/* Return a list of nodes which belong the a maximal independent set.
   Note: Only nodes local to this rank will be returned.

   Caller is responsible for cleaning up return value.
*/
index_t Paso_SystemMatrix_getMIS(Paso_SystemMatrix* A, index_t** set) {
/* identify the nodes on my border */
    index_t count=0;
    index_t missize=0;
    index_t* border=NULL;
    index_t n=A->mainBlock->numRows;
    char* inborder=NULL;
    double* weights=NULL;
    index_t* mis=NULL;
    Paso_Pattern* pat=A->mainBlock->pattern;
    int i,j,k, retry;
    char done=1;
    double seed=Paso_Pattern_mis_seed;
    Esys_resetError();
    border=Paso_SparseMatrix_getBorderNodes(A, &count);
    if (!Esys_noError()) {
        *set=NULL;
	return 0;
    }

    inborder=MEMALLOC(n,char);
    weights=MEMALLOC(n,double);
    mis=MEMALLOC(n,index_t);

    #pragma omp parallel for schedule(static) private(i)
    for (i=0;i<n;++i) {
	inborder[i]=0;
	weights[i]=0.5;
    }
   /* This loop will use a different memory access pattern to the setup
      ie a NUMA problem. However I'm gambling that the border is small and
      that the memory won't be relocated.
   */
    #pragma omp parallel for schedule(static) private(i)
    for (i=0;i<count;++i) {
	inborder[border[i]]=1;
    }

    /* compute set membership for border nodes */
    Paso_SystemMatrix_CalcBorderMIS(A, border, count, weights, n);
    
    do {
	done=1;
	retry=0;
	j=0;
	/* Randomise weights of unmarked nodes */
	#pragma omp parallel for schedule(static) private(i)
	for (i=0;i<n;++i) {
		if (ISAVAILABLE(weights[i])) {
		    weights[i]=fmod(seed*(i+1),1.)+0.1;
		    ++j;
		}
	}
	k=0;
	/* Now we evaluate nodes relative to their neighbours */
	#pragma omp parallel for schedule(static) private(i)
	for (i=0;i<n;++i) {
		char ok=1;
		for (j=pat->ptr[i];j<pat->ptr[i+1];++j) {
		    if (i!=pat->index[j] && !ISLESS(weights[i],weights[pat->index[j]])) {
			ok=0;
			break;
		    }
		}
		if (ok) {
		    weights[i]=MISIN;
		    k++;
		    retry=1;	/* multiple threads writing same value should be fine */
		}
	}

        /* Go through and mark all the neighbours of nodes definitly in */
	#pragma omp parallel for schedule(static) private(i)
	for (i=0;i<n;++i) {
		if (weights[i]==MISIN) {
		    for (j=pat->ptr[i];j<pat->ptr[i+1];++j) {
			if (pat->index[j]!=i) {
			    weights[pat->index[j]]=MISOUT;
			}
		    }
		}
	}
	
	/* Are we finished yet? */
	#pragma omp parallel for schedule(static) private(i)
	for (i=0;i<n;++i) {
		if ((weights[i]!=MISIN) && (weights[i]!=MISOUT)) {	/* neither in or neighbour of in */
			done=0;		/* multiple threads writing same value should be fine */
		}
	}
	seed+=0.014159; /* Arbitrary change to seed */
    } while (!done && retry);

    for (i=0;i<n;++i) {
        if (weights[i]==MISIN) {
	    mis[missize++]=i;
	}
    }
    MEMFREE(border);
    MEMFREE(weights);
    MEMFREE(inborder);
    if (done==0) {
        Esys_setError(NO_PROGRESS_ERROR,"Error in MIS - no progress.");
	MEMFREE(mis);
	*set=NULL;
	return 0;
    }
    *set=mis;
    return missize;
}