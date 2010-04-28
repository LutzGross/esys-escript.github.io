
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

  if (type & MATRIX_FORMAT_SYM) {
     Paso_setError(TYPE_ERROR,"Paso_SparseMatrix_alloc: symmetric matrix pattern are not supported.");
     return NULL;
  }

  if (patternIsUnrolled) {
     if (! XNOR(type & MATRIX_FORMAT_OFFSET1, pattern->type & PATTERN_FORMAT_OFFSET1) ) {
         Paso_setError(TYPE_ERROR,"Paso_SparseMatrix_alloc: requested offset and pattern offset does not match.");
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
    || ( (type & MATRIX_FORMAT_OFFSET1) != ( pattern->type & PATTERN_FORMAT_OFFSET1) ) ;

  pattern_format_out= (type & MATRIX_FORMAT_OFFSET1)? PATTERN_FORMAT_OFFSET1:  PATTERN_FORMAT_DEFAULT;

  Paso_resetError();
  out=MEMALLOC(1,Paso_SparseMatrix);
  if (! Paso_checkPtr(out)) {  
     out->pattern=NULL;  
     out->val=NULL;  
     out->reference_counter=1;
     out->type=type;
     out->solver=NULL; 

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
           if (Paso_noError()) {
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
           if (Paso_noError()) {
               out->numRows = out->pattern->numOutput;
               out->numCols = out->pattern->numInput;
           }
     }
     if (Paso_noError()) {
         out->block_size=out->row_block_size*out->col_block_size;
         out->len=(size_t)(out->pattern->len)*(size_t)(out->block_size);
    
         out->val=MEMALLOC(out->len,double);
         if (! Paso_checkPtr(out->val)) Paso_SparseMatrix_setValues(out,DBLE(0));
     }
  }
  /* all done: */
  if (Paso_noError()) {
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
        Paso_Pattern_free(in->pattern);
        MEMFREE(in->val);
	if(in->solver!=NULL) {
		MEMFREE(in->solver);
	}
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
        Paso_resetError();

	/* open the file */
	fileHandle_p = fopen( fileName_p, "r" );
	if( fileHandle_p == NULL )
	{
		Paso_setError(IO_ERROR, "Paso_SparseMatrix_loadMM_toCSR: Cannot read file for reading.");
		return NULL;
	}

	/* process banner */
	if( mm_read_banner(fileHandle_p, &matrixCode) != 0 )
	{
		Paso_setError(IO_ERROR, "Paso_SparseMatrix_loadMM_toCSR: Error processing MM banner.");
		fclose( fileHandle_p );
		return NULL;
	}
	if( !(mm_is_real(matrixCode) && mm_is_sparse(matrixCode) && mm_is_general(matrixCode)) )
	{

		Paso_setError(TYPE_ERROR,"Paso_SparseMatrix_loadMM_toCSR: found Matrix Market type is not supported.");
		fclose( fileHandle_p );
		return NULL;
	}

	/* get matrix size */
	if( mm_read_mtx_crd_size(fileHandle_p, &M, &N, &nz) != 0 )
	{
		Paso_setError(IO_ERROR, "Paso_SparseMatrix_loadMM_toCSR: Could not parse matrix size");
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
		Paso_setError(MEMORY_ERROR, "Paso_SparseMatrix_loadMM_toCSR: Could not allocate memory" );
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

        mainPattern=Paso_Pattern_alloc(PATTERN_FORMAT_DEFAULT,M,N,row_ptr,col_ind);
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
  dim_t N,M,i, iptr_ij;
  MM_typecode matcode;                        

  if (A_p->col_block_size !=A_p->row_block_size) {
    Paso_setError(TYPE_ERROR, "Paso_SparseMatrix_saveMM: currently only square blocks are supported.");
    return;
  }
  if (A_p->row_block_size>3) {
       Paso_setError(TYPE_ERROR,"Paso_SparseMatrix_saveMM: currently only block size 3 is supported.\n");
       return;
  }

  if (A_p->type & MATRIX_FORMAT_SYM) {
    Paso_setError(TYPE_ERROR,"Paso_SparseMatrix_saveMM does not support symmetric storage scheme");
    return;
  }
  
  /* open the file */
  fileHandle_p = fopen(fileName_p, "w");
  if (fileHandle_p==NULL) {
    Paso_setError(IO_ERROR,"file could not be opened for writing");
    return;
  }

  if (A_p->type & MATRIX_FORMAT_CSC) {
    Paso_setError(TYPE_ERROR,"Paso_SparseMatrix_saveMM does not support CSC yet.");
  } else {
    mm_initialize_typecode(&matcode);
    mm_set_matrix(&matcode);
    mm_set_coordinate(&matcode);
    mm_set_real(&matcode);

    N= A_p->numRows; 
    M= A_p->numCols; 
    mm_write_banner(fileHandle_p, matcode); 
    mm_write_mtx_crd_size(fileHandle_p, N, M, A_p->pattern->ptr[N]);
    if(A_p->row_block_size==1)
      for (i=0; i<N; i++) 
        for (iptr_ij=A_p->pattern->ptr[i];iptr_ij<A_p->pattern->ptr[i+1]; ++iptr_ij) {
          fprintf(fileHandle_p, "%d %d %25.15e\n", i+1, A_p->pattern->index[iptr_ij]+1, A_p->val[iptr_ij]);
        }

    if(A_p->row_block_size==2)
      for (i=0; i<N; i++) {
        for (iptr_ij=A_p->pattern->ptr[i];iptr_ij<A_p->pattern->ptr[i+1]; ++iptr_ij) {
          fprintf(fileHandle_p, "%d %d %25.15e\n", i*2+1, 2*(A_p->pattern->index[iptr_ij])+1, A_p->val[iptr_ij*4]);
          fprintf(fileHandle_p, "%d %d %25.15e\n", i*2+1, 2*(A_p->pattern->index[iptr_ij])+1+1, A_p->val[iptr_ij*4+2]);
          fprintf(fileHandle_p, "%d %d %25.15e\n", i*2+1+1, 2*(A_p->pattern->index[iptr_ij])+1, A_p->val[iptr_ij*4+1]);
          fprintf(fileHandle_p, "%d %d %25.15e\n", i*2+1+1, 2*(A_p->pattern->index[iptr_ij])+1+1, A_p->val[iptr_ij*4+3]);
        }
      }

    if(A_p->row_block_size==3)
      for (i=0; i<N; i++) {
       for (iptr_ij=A_p->pattern->ptr[i];iptr_ij<A_p->pattern->ptr[i+1]; ++iptr_ij) {
          fprintf(fileHandle_p, "%d %d %25.15e\n", i*3+1, 3*(A_p->pattern->index[iptr_ij])+1, A_p->val[iptr_ij*9]);
          fprintf(fileHandle_p, "%d %d %25.15e\n", i*3+1, 3*(A_p->pattern->index[iptr_ij])+1+1, A_p->val[iptr_ij*9+3]);
          fprintf(fileHandle_p, "%d %d %25.15e\n", i*3+1, 3*(A_p->pattern->index[iptr_ij])+2+1, A_p->val[iptr_ij*9+6]);
          fprintf(fileHandle_p, "%d %d %25.15e\n", i*3+1+1, 3*(A_p->pattern->index[iptr_ij])+1, A_p->val[iptr_ij*9+1]);
          fprintf(fileHandle_p, "%d %d %25.15e\n", i*3+1+1, 3*(A_p->pattern->index[iptr_ij])+1+1, A_p->val[iptr_ij*9+4]);
          fprintf(fileHandle_p, "%d %d %25.15e\n", i*3+1+1, 3*(A_p->pattern->index[iptr_ij])+2+1, A_p->val[iptr_ij*9+7]);
          fprintf(fileHandle_p, "%d %d %25.15e\n", i*3+2+1, 3*(A_p->pattern->index[iptr_ij])+1, A_p->val[iptr_ij*9+2]);
          fprintf(fileHandle_p, "%d %d %25.15e\n", i*3+2+1, 3*(A_p->pattern->index[iptr_ij])+1+1, A_p->val[iptr_ij*9+5]);
          fprintf(fileHandle_p, "%d %d %25.15e\n", i*3+2+1, 3*(A_p->pattern->index[iptr_ij])+2+1, A_p->val[iptr_ij*9+8]);
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

Paso_SparseMatrix* Paso_SparseMatrix_unroll(Paso_SparseMatrix* A) {
	Paso_SparseMatrix *out = NULL;
	/*Paso_Pattern* mainPattern=NULL;*/
	int i;
	dim_t output_block_size=A->row_block_size;
	dim_t input_block_size=A->col_block_size;
	dim_t block_size=A->row_block_size;
	index_t iptr, optr1=0, optr2=0, optr3=0;
	
	/*mainPattern= Paso_Pattern_unrollBlocks(A->pattern, PATTERN_FORMAT_DEFAULT, output_block_size,input_block_size);*/
	out  = Paso_SparseMatrix_alloc(MATRIX_FORMAT_BLK1, A->pattern, output_block_size,input_block_size, FALSE);
	
	#pragma omp parallel for private(i,optr1,optr2,optr3,iptr) schedule(static)
	for (i=0;i<A->numRows;++i) {
		if (block_size==1) {
		    optr1=out->pattern->ptr[i];
		} else if (block_size==2) {
		    optr1=out->pattern->ptr[i*block_size];
		    optr2=out->pattern->ptr[i*block_size+1];
		} else if (block_size==3) {
		    optr1=out->pattern->ptr[i*block_size];
		    optr2=out->pattern->ptr[i*block_size+1];
		    optr3=out->pattern->ptr[i*block_size+2];
		}
            for (iptr=A->pattern->ptr[i];iptr<A->pattern->ptr[i+1]; ++iptr) {
                if (block_size==1) {
			out->val[optr1]=A->val[iptr];
			optr1++;
		} else if (block_size==2) {
			out->val[optr1]=A->val[iptr*block_size*block_size];
			optr1++;
			out->val[optr1]=A->val[iptr*block_size*block_size+1];
			optr1++;
			out->val[optr2]=A->val[iptr*block_size*block_size+2];
			optr2++;
			out->val[optr2]=A->val[iptr*block_size*block_size+3];
			optr2++;
		} else if (block_size==3) {
			out->val[optr1]=A->val[iptr*block_size*block_size];
			optr1++;
			out->val[optr1]=A->val[iptr*block_size*block_size+1];
			optr1++;
			out->val[optr1]=A->val[iptr*block_size*block_size+2];
			optr1++;
			out->val[optr2]=A->val[iptr*block_size*block_size+3];
			optr2++;
			out->val[optr2]=A->val[iptr*block_size*block_size+4];
			optr2++;
			out->val[optr2]=A->val[iptr*block_size*block_size+5];
			optr2++;
			out->val[optr3]=A->val[iptr*block_size*block_size+6];
			optr3++;
			out->val[optr3]=A->val[iptr*block_size*block_size+7];
			optr3++;
			out->val[optr3]=A->val[iptr*block_size*block_size+8];
			optr3++;
		}
            }
	}

 return out;
}