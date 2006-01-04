/* $Id$ */

/**************************************************************/

/* Paso: Matrix Market format is loaded to a SystemMatrix   */

/**************************************************************/

/* Copyrights by ACcESS Australia 2003,2004,2005 */
/* Author: imran@access.edu.au */

/**************************************************************/

#include "Paso.h"
#include "mmio.h"
#include "SystemMatrix.h"

static void swap( index_t*, index_t*, double*, int, int );
static void q_sort( index_t*, index_t*, double*, int, int );
static void print_entries( index_t*, index_t*, double* );

static int M, N, nz;


/* debug: print the entries */
void print_entries( index_t *r, index_t *c, double *v )
{
	int i;

	for( i=0; i<nz; i++ )
	{
		printf( "(%ld, %ld) == %e\n", (long)r[i], (long)c[i], v[i] );
	}
}

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

Paso_SystemMatrix* Paso_SystemMatrix_loadMM_toCSR( char *fileName_p )
{
	int i, curr_row;
	MM_typecode matrixCode;
        Paso_resetError();

	index_t *col_ind = NULL;
	index_t *row_ind = NULL;
	index_t *row_ptr = NULL;
	double *val = NULL;

	Paso_SystemMatrixPattern *loc_pattern = NULL;
	Paso_SystemMatrix *out = NULL;

	/* open the file */
	FILE *fileHandle_p = fopen( fileName_p, "r" );
	if( fileHandle_p == NULL )
	{
		Paso_setError(IO_ERROR, "Paso_SystemMatrix_loadMM_toCSR: Cannot read file for reading.");
		return NULL;
	}

	/* process banner */
	if( mm_read_banner(fileHandle_p, &matrixCode) != 0 )
	{
		Paso_setError(IO_ERROR, "Paso_SystemMatrix_loadMM_toCSR: Error processing MM banner.");
		fclose( fileHandle_p );
		return NULL;
	}
	if( !(mm_is_real(matrixCode) && mm_is_sparse(matrixCode) && mm_is_general(matrixCode)) )
	{

		Paso_setError(TYPE_ERROR,"Paso_SystemMatrix_loadMM_toCSR: found Matrix Market type is not supported.");
		fclose( fileHandle_p );
		return NULL;
	}

	/* get matrix size */
	if( mm_read_mtx_crd_size(fileHandle_p, &M, &N, &nz) != 0 )
	{
		Paso_setError(IO_ERROR, "Paso_SystemMatrix_loadMM_toCSR: Could not parse matrix size");
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
		Paso_setError(MEMORY_ERROR, "Paso_SystemMatrix_loadMM_toCSR: Could not allocate memory" );

		fclose( fileHandle_p );
		return NULL;
	}

	/* perform actual read of elements */
	for( i=0; i<nz; i++ )
	{
		fscanf( fileHandle_p, "%d %d %le\n", &row_ind[i], &col_ind[i], &val[i] );
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

	/* create F_SMP and F_SM */
	loc_pattern = Paso_SystemMatrixPattern_alloc(PATTERN_FORMAT_DEFAULT, M, row_ptr, col_ind );
	if(! Paso_noError() )
		return NULL;

 	out = Paso_SystemMatrix_alloc(MATRIX_FORMAT_DEFAULT, loc_pattern, 1, 1 );
 	if(! Paso_noError() )
 		return NULL;

	/* copy values and cleanup temps */
	for( i=0; i<nz; i++ )
		out->val[i] = val[i];

	MEMFREE( val );
	MEMFREE( row_ind );

	return out;
}

Paso_SystemMatrix* Paso_SystemMatrix_loadMM_toCSC( char *fileName_p )
{
	int i;
//	int curr_row;
	int curr_col;
	MM_typecode matrixCode;

	index_t *col_ind = NULL;
	index_t *row_ind = NULL;
//	index_t *row_ptr = NULL;
	index_t *col_ptr = NULL;
	double *val = NULL;

	Paso_SystemMatrixPattern *loc_pattern = NULL;
	Paso_SystemMatrix *out = NULL;

	Paso_resetError();

	/* open the file */
	FILE *fileHandle_p = fopen( fileName_p, "r" );
	if( fileHandle_p == NULL )
	{
		Paso_setError(IO_ERROR,"Paso_SystemMatrix_loadMM_toCSC: File could not be opened for reading");
		return NULL;
	}

	/* process banner */
	if( mm_read_banner(fileHandle_p, &matrixCode) != 0 )
	{
		Paso_setError(IO_ERROR,"Paso_SystemMatrix_loadMM_toCSC: Error processing MM banner");
		fclose( fileHandle_p );
		return NULL;
	}
	if( !(mm_is_real(matrixCode) && mm_is_sparse(matrixCode) && mm_is_general(matrixCode)) )
	{
		Paso_setError(TYPE_ERROR,"Paso_SystemMatrix_loadMM_toCSC: found Matrix Market type is not supported.");
		fclose( fileHandle_p );
		return NULL;
	}

	/* get matrix size */
	if( mm_read_mtx_crd_size(fileHandle_p, &M, &N, &nz) != 0 )
	{
		Paso_setError(TYPE_ERROR,"Paso_SystemMatrix_loadMM_toCSC: found Matrix Market type is not supported.");
		fclose( fileHandle_p );
		return NULL;
	}

	/* prepare storage */
	col_ind = MEMALLOC( nz, index_t );
	row_ind = MEMALLOC( nz, index_t );
	val = MEMALLOC( nz, double );

//	row_ptr = MEMALLOC( (M+1), index_t );
	col_ptr = MEMALLOC( (N+1), index_t );


//	if( col_ind == NULL || row_ind == NULL || val == NULL || row_ptr == NULL )
	if( col_ind == NULL || row_ind == NULL || val == NULL || col_ptr == NULL )
	{
		Paso_setError(MEMORY_ERROR,"Could not allocate memory" );
		fclose( fileHandle_p );
		return NULL;
	}

	/* perform actual read of elements */
	for( i=0; i<nz; i++ )
	{
		fscanf( fileHandle_p, "%d %d %le\n", &row_ind[i], &col_ind[i], &val[i] );
		row_ind[i]--;
		col_ind[i]--;
	}
	fclose( fileHandle_p );

	/* sort the entries */
//	q_sort( row_ind, col_ind, val, 0, nz );
	q_sort( col_ind, row_ind, val, 0, nz );

	/* setup row_ptr */
//	curr_row = 0;
	curr_col = 0;
//	for( i=0; (i<nz && curr_row<M); curr_row++ )
	for( i=0; (i<nz && curr_col<N); curr_col++ )
	{
//		while( row_ind[i] != curr_row )
		while( col_ind[i] != curr_col )
			i++;
//		row_ptr[curr_row] = i;
		col_ptr[curr_col] = i;
	}
//	row_ptr[M] = nz;
	col_ptr[N] = nz;

	/* create F_SMP and F_SM */
//	loc_pattern = Paso_SystemMatrixPattern_alloc(PATTERN_FORMAT_DEFAULT, M, row_ptr, col_ind );
	loc_pattern = Paso_SystemMatrixPattern_alloc(PATTERN_FORMAT_DEFAULT, N, col_ptr, row_ind );
	if(! Paso_noError() )
		return NULL;

//      out = Paso_SystemMatrix_alloc( MATRIX_FORMAT_DEFAULT, loc_pattern, 1, 1 );
 	out = Paso_SystemMatrix_alloc(  MATRIX_FORMAT_CSC, loc_pattern, 1, 1 );
 	if(! Paso_noError() )
 		return NULL;

	/* copy values and cleanup temps */
	for( i=0; i<nz; i++ )
		out->val[i] = val[i];

	MEMFREE( val );
//	MEMFREE( row_ind );
	MEMFREE( col_ind );

	return out;
}


/*
 * $Log$
 * Revision 1.2  2005/09/15 03:44:39  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-15
 *
 * Revision 1.1.2.1  2005/09/05 06:29:48  gross
 * These files have been extracted from finley to define a stand alone libray for iterative
 * linear solvers on the ALTIX. main entry through Paso_solve. this version compiles but
 * has not been tested yet.
 *
 *
 */
