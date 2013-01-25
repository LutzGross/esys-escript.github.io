
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

/* Paso: SparseMatrix */

/**************************************************************/

/* Author: Lutz Gross, l.gross@uq.edu.au */

/**************************************************************/

#include "Paso.h"
#include "SparseMatrix.h"
#include "MKL.h"
#include "Preconditioner.h"
#include "UMFPACK.h"
#include "TRILINOS.h"
#include "mmio.h"

/**************************************************************/
static void swap( index_t*, index_t*, double*, int, int );
static void q_sort( index_t*, index_t*, double*, int, int );
/*static void print_entries( index_t*, index_t*, double* );*/

static int M, N, nz;


/* debug: print the entries */
/*
void print_entries( index_t *r, index_t *c, double *v )
{
	int i;

	for( i=0; i<nz; i++ )
	{
		printf( "(%ld, %ld) == %e\n", (long)r[i], (long)c[i], v[i] );
	}
}
*/

/* swap function */
void swap( index_t *r, index_t *c, double *v, int left, int right )
{
	double v_temp;
	index_t temp;

	temp = r[left];
	r[left] = r[right];
	r[right] = temp;

	temp = c[left];
	c[left] = c[right];
	c[right] = temp;

	v_temp = v[left];
	v[left] = v[right];
	v[right] = v_temp;
}

void q_sort( index_t *row, index_t *col, double *val, int begin, int end )
{
	int l, r;
	index_t pivot, lval;

	if( end > begin )
	{
		pivot = N * row[begin] + col[begin];
		l = begin + 1;
		r = end;

		while( l < r )
		{
			lval = N * row[l] + col[l];
			if( lval < pivot )
				l++;
			else
			{
				r--;
				swap( row, col, val, l, r );
			}
		}
		l--;
		swap( row, col, val, begin, l );
		q_sort( row, col, val, begin, l );
		q_sort( row, col, val, r, end );
	}
}


/* allocates a SparseMatrix of type type using the given matrix pattern 
   Values are initialized by zero. 
   if patternIsUnrolled and type & MATRIX_FORMAT_BLK1, it is assumed that the pattern is allready unrolled to match the requested block size
   and offsets otherwise unrolling and offset adjustment will be performed. 
*/
Paso_SparseMatrix* Paso_SparseMatrix_alloc(Paso_SparseMatrixType type,Paso_Pattern *pattern, int row_block_size, int col_block_size, const bool_t patternIsUnrolled) {

  Paso_SparseMatrix*out=NULL;
  Paso_SparseMatrixType pattern_format_out;
  bool_t unroll=FALSE;

  if (patternIsUnrolled) {
     if (! XNOR(type & MATRIX_FORMAT_OFFSET1, pattern->type & MATRIX_FORMAT_OFFSET1) ) {
         Esys_setError(TYPE_ERROR,"Paso_SparseMatrix_alloc: requested offset and pattern offset does not match.");
         return NULL;
     }
  }
  /* do we need to apply unrolling ? */
  unroll
        /* we don't like non-square blocks */
    =   (row_block_size!=col_block_size)
        /* or any block size bigger than 3 */
    ||  (col_block_size>3)
        /* or if lock size one requested and the block size is not 1 */
    ||  ((type & MATRIX_FORMAT_BLK1) &&  (col_block_size>1) ) 
        /* offsets don't match */
    || ( (type & MATRIX_FORMAT_OFFSET1) != ( pattern->type & MATRIX_FORMAT_OFFSET1) ) ;

  pattern_format_out= (type & MATRIX_FORMAT_OFFSET1)? MATRIX_FORMAT_OFFSET1:  MATRIX_FORMAT_DEFAULT;

  Esys_resetError();
  out=MEMALLOC(1,Paso_SparseMatrix);
  if (! Esys_checkPtr(out)) {  
     out->pattern=NULL;  
     out->val=NULL;  
     out->reference_counter=1;
     out->type=type;
     out->solver_package=PASO_PASO;  
     out->solver_p=NULL;  
     
     /* ====== compressed sparse columns === */
     if (type & MATRIX_FORMAT_CSC) {
           if (unroll) {
              if (patternIsUnrolled) {
                    out->pattern=Paso_Pattern_getReference(pattern); 
              } else {
                    out->pattern=Paso_Pattern_unrollBlocks(pattern,pattern_format_out,col_block_size,row_block_size); 
              }
              out->row_block_size=1;
              out->col_block_size=1;
           } else {
             out->pattern=Paso_Pattern_unrollBlocks(pattern,pattern_format_out,1,1);
              out->row_block_size=row_block_size;
              out->col_block_size=col_block_size;
           }
           if (Esys_noError()) {
              out->numRows = out->pattern->numInput;
              out->numCols = out->pattern->numOutput;
           }
 
     } else {
     /* ====== compressed sparse row === */
           if (unroll) {
              if (patternIsUnrolled) {
                   out->pattern=Paso_Pattern_getReference(pattern); 
              } else {
                   out->pattern=Paso_Pattern_unrollBlocks(pattern,pattern_format_out,row_block_size,col_block_size);
              }
              out->row_block_size=1;
              out->col_block_size=1;
           } else {
              out->pattern=Paso_Pattern_unrollBlocks(pattern,pattern_format_out,1,1);
              out->row_block_size=row_block_size;
              out->col_block_size=col_block_size;
           }
           if (Esys_noError()) {
               out->numRows = out->pattern->numOutput;
               out->numCols = out->pattern->numInput;
           }
     }
     if (Esys_noError()) {
	 if (type & MATRIX_FORMAT_DIAGONAL_BLOCK) {
	    out->block_size=MIN(out->row_block_size,out->col_block_size);
	 } else {
	    out->block_size=out->row_block_size*out->col_block_size;
	 }
         out->len=(size_t)(out->pattern->len)*(size_t)(out->block_size);
    
         out->val=MEMALLOC(out->len,double);
         if (! Esys_checkPtr(out->val)) Paso_SparseMatrix_setValues(out,DBLE(0));
     }
  }
  /* all done: */
  if (Esys_noError()) {
    return out;
  } else {
    Paso_SparseMatrix_free(out);
    return NULL;
  }
}

/* returns a reference to Paso_SparseMatrix in */

Paso_SparseMatrix* Paso_SparseMatrix_getReference(Paso_SparseMatrix* in) {
   if (in!=NULL) ++(in->reference_counter);
   return in;
}

/* deallocates a SparseMatrix: */

void Paso_SparseMatrix_free(Paso_SparseMatrix* in) {
  if (in!=NULL) {
     in->reference_counter--;
     if (in->reference_counter<=0) {
	switch(in->solver_package) {
	   
	    case PASO_SMOOTHER:
	       Paso_Preconditioner_LocalSmoother_free((Paso_Preconditioner_LocalSmoother*) in->solver_p);
	       break;
	   
	    case PASO_MKL:
	       Paso_MKL_free(in); /* releases solver_p */
	       break;
	   
	    case PASO_UMFPACK:
	       Paso_UMFPACK_free(in); /* releases solver_p */
	       break;
	}
	MEMFREE(in->val);
	Paso_Pattern_free(in->pattern);
        MEMFREE(in);
     }
   }
}

Paso_SparseMatrix* Paso_SparseMatrix_loadMM_toCSR( char *fileName_p )
{
	index_t *col_ind = NULL;
	index_t *row_ind = NULL;
	index_t *row_ptr = NULL;
	double *val = NULL;
	FILE *fileHandle_p = NULL;
        Paso_Pattern* mainPattern=NULL;
	Paso_SparseMatrix *out = NULL;
	int i, curr_row, scan_ret;
	MM_typecode matrixCode;
        Esys_resetError();

	/* open the file */
	fileHandle_p = fopen( fileName_p, "r" );
	if( fileHandle_p == NULL )
	{
		Esys_setError(IO_ERROR, "Paso_SparseMatrix_loadMM_toCSR: Cannot read file for reading.");
		return NULL;
	}

	/* process banner */
	if( mm_read_banner(fileHandle_p, &matrixCode) != 0 )
	{
		Esys_setError(IO_ERROR, "Paso_SparseMatrix_loadMM_toCSR: Error processing MM banner.");
		fclose( fileHandle_p );
		return NULL;
	}
	if( !(mm_is_real(matrixCode) && mm_is_sparse(matrixCode) && mm_is_general(matrixCode)) )
	{

		Esys_setError(TYPE_ERROR,"Paso_SparseMatrix_loadMM_toCSR: found Matrix Market type is not supported.");
		fclose( fileHandle_p );
		return NULL;
	}

	/* get matrix size */
	if( mm_read_mtx_crd_size(fileHandle_p, &M, &N, &nz) != 0 )
	{
		Esys_setError(IO_ERROR, "Paso_SparseMatrix_loadMM_toCSR: Could not parse matrix size");
		fclose( fileHandle_p );
		return NULL;
	}

	/* prepare storage */
	col_ind = MEMALLOC( nz, index_t );
	row_ind = MEMALLOC( nz, index_t );
	val = MEMALLOC( nz, double );

	row_ptr = MEMALLOC( (M+1), index_t );

	if( col_ind == NULL || row_ind == NULL || val == NULL || row_ptr == NULL )
	{
		Esys_setError(MEMORY_ERROR, "Paso_SparseMatrix_loadMM_toCSR: Could not allocate memory" );
		fclose( fileHandle_p );
		return NULL;
	}

	/* perform actual read of elements */
	for( i=0; i<nz; i++ )
	{
		scan_ret = fscanf( fileHandle_p, "%d %d %le\n", &row_ind[i], &col_ind[i], &val[i] );
		if (scan_ret!=3)
		{
			MEMFREE( val );
			MEMFREE( row_ind );
			MEMFREE( col_ind );
			MEMFREE( row_ptr );
			fclose(fileHandle_p);
			return NULL;
		}
		row_ind[i]--;
		col_ind[i]--;
	}
	fclose( fileHandle_p );

	/* sort the entries */
	q_sort( row_ind, col_ind, val, 0, nz );

	/* setup row_ptr */
	curr_row = 0;
	for( i=0; (i<nz && curr_row<M); curr_row++ )
	{
		while( row_ind[i] != curr_row )
			i++;
		row_ptr[curr_row] = i;
	}
	row_ptr[M] = nz;

        mainPattern=Paso_Pattern_alloc(MATRIX_FORMAT_DEFAULT,M,N,row_ptr,col_ind);
	out  = Paso_SparseMatrix_alloc(MATRIX_FORMAT_DEFAULT, mainPattern, 1, 1, TRUE);
	/* copy values and cleanup temps */
	for( i=0; i<nz; i++ ) out->val[i] = val[i];
	
	Paso_Pattern_free(mainPattern);
	MEMFREE( val );
	MEMFREE( row_ind );
	return out;
}


void Paso_SparseMatrix_saveMM(Paso_SparseMatrix * A_p, char * fileName_p) {
   FILE * fileHandle_p = NULL;
   dim_t N,M,i,j, irow, icol, ib, irb, icb, iptr;
   MM_typecode matcode;  
   const dim_t col_block_size = A_p->col_block_size;
   const dim_t row_block_size = A_p->row_block_size;
   const dim_t block_size = A_p->block_size;
   
   if (col_block_size !=row_block_size) {
      Esys_setError(TYPE_ERROR, "Paso_SparseMatrix_saveMM: currently only square blocks are supported.");
      return;
   }
   
   /* open the file */
   fileHandle_p = fopen(fileName_p, "w");
   if (fileHandle_p==NULL) {
      Esys_setError(IO_ERROR,"file could not be opened for writing");
      return;
   }
   if (A_p->type & MATRIX_FORMAT_CSC) {
      Esys_setError(TYPE_ERROR,"Paso_SparseMatrix_saveMM does not support CSC yet.");
   } else {
      mm_initialize_typecode(&matcode);
      mm_set_matrix(&matcode);
      mm_set_coordinate(&matcode);
      mm_set_real(&matcode);
      
      N=Paso_SparseMatrix_getNumRows(A_p);
      M=Paso_SparseMatrix_getNumCols(A_p);
      mm_write_banner(fileHandle_p, matcode); 
      mm_write_mtx_crd_size(fileHandle_p, N*row_block_size, M*col_block_size, A_p->pattern->ptr[N]*block_size);
      
      if (A_p->type & MATRIX_FORMAT_DIAGONAL_BLOCK) {
	 for (i=0; i<N; i++) {
	    for (iptr = A_p->pattern->ptr[i];iptr<A_p->pattern->ptr[i+1]; ++iptr) {
	       j=A_p->pattern->index[iptr];
	       for (ib=0;ib<block_size;ib++) {
		  irow=ib+row_block_size*i;
		  icol=ib+col_block_size*j;
		  fprintf(fileHandle_p, "%d %d %25.15e\n", irow+1, icol+1, A_p->val[iptr*block_size+ib]);
	       }
	    }
	 }
      } else {
	 for (i=0; i<N; i++) {
	    for (iptr = A_p->pattern->ptr[i];iptr<A_p->pattern->ptr[i+1]; ++iptr) {
	       j=A_p->pattern->index[iptr]; 
	       for (irb=0;irb<row_block_size;irb++) {
		  irow=irb+row_block_size*i;
		  for (icb=0;icb<col_block_size;icb++) {
		     icol=icb+col_block_size*j;
		     fprintf(fileHandle_p, "%d %d %25.15e\n", irow+1, icol+1, A_p->val[iptr*block_size+irb+row_block_size*icb]);
		  }
	       }
	    }
	 }
      }
   }
   /* close the file */
   fclose(fileHandle_p);
   return;
}

index_t* Paso_SparseMatrix_borrowMainDiagonalPointer(Paso_SparseMatrix * A_p) 
{
    return Paso_Pattern_borrowMainDiagonalPointer(A_p->pattern);
}

dim_t Paso_SparseMatrix_getNumColors(Paso_SparseMatrix* A_p) 
{
   return Paso_Pattern_getNumColors(A_p->pattern);
}
index_t* Paso_SparseMatrix_borrowColoringPointer(Paso_SparseMatrix* A_p)
{
   return Paso_Pattern_borrowColoringPointer(A_p->pattern);
}
dim_t Paso_SparseMatrix_maxDeg(Paso_SparseMatrix * A_p)
{
   return Paso_Pattern_maxDeg(A_p->pattern);
}
dim_t Paso_SparseMatrix_getTotalNumRows(const Paso_SparseMatrix* A){
   return A->numRows * A->row_block_size;
}

dim_t Paso_SparseMatrix_getTotalNumCols(const Paso_SparseMatrix* A){
   return A->numCols * A->col_block_size;
}
dim_t Paso_SparseMatrix_getNumRows(Paso_SparseMatrix* A) {
   return A->numRows;
}
dim_t Paso_SparseMatrix_getNumCols(Paso_SparseMatrix* A) {
   return A->numCols;  
}
