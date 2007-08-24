/* $Id$ */


/*
********************************************************************************
*               Copyright   2006 by ACcESS MNRF                                *
*                                                                              * 
*                 http://www.access.edu.au                                     *
*           Primary Business: Queensland, Australia                            *
*     Licensed under the Open Software License version 3.0 		       *
*        http://www.opensource.org/licenses/osl-3.0.php                        *
********************************************************************************
*/

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
        index_t dist[2];
        Paso_Distribution* input_dist=NULL, *output_dist=NULL;
	index_t *col_ind = NULL;
	index_t *row_ind = NULL;
	index_t *row_ptr = NULL;
	double *val = NULL;
	FILE *fileHandle_p = NULL;
        Paso_Pattern* mainPattern=NULL, *couplePattern=NULL;
	Paso_SystemMatrixPattern *pattern = NULL;
	Paso_SystemMatrix *out = NULL;
        Paso_SharedComponents *send =NULL;
        Paso_Coupler *coupler=NULL;
	int i, curr_row;
	MM_typecode matrixCode;
        Paso_MPIInfo* mpi_info=Paso_MPIInfo_alloc( MPI_COMM_WORLD);
        Paso_resetError();
        if (mpi_info->size >1) {
		Paso_setError(IO_ERROR, "Paso_SystemMatrix_loadMM_toCSR: support single processor only");
		return NULL;
        }
	/* open the file */
	fileHandle_p = fopen( fileName_p, "r" );
	if( fileHandle_p == NULL )
	{
		Paso_setError(IO_ERROR, "Paso_SystemMatrix_loadMM_toCSR: Cannot read file for reading.");
                Paso_MPIInfo_free(mpi_info);
		return NULL;
	}

	/* process banner */
	if( mm_read_banner(fileHandle_p, &matrixCode) != 0 )
	{
		Paso_setError(IO_ERROR, "Paso_SystemMatrix_loadMM_toCSR: Error processing MM banner.");
                Paso_MPIInfo_free(mpi_info);
		fclose( fileHandle_p );
		return NULL;
	}
	if( !(mm_is_real(matrixCode) && mm_is_sparse(matrixCode) && mm_is_general(matrixCode)) )
	{

		Paso_setError(TYPE_ERROR,"Paso_SystemMatrix_loadMM_toCSR: found Matrix Market type is not supported.");
                Paso_MPIInfo_free(mpi_info);
		fclose( fileHandle_p );
		return NULL;
	}

	/* get matrix size */
	if( mm_read_mtx_crd_size(fileHandle_p, &M, &N, &nz) != 0 )
	{
		Paso_setError(IO_ERROR, "Paso_SystemMatrix_loadMM_toCSR: Could not parse matrix size");
                Paso_MPIInfo_free(mpi_info);
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

                Paso_MPIInfo_free(mpi_info);
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

        /* create return value */
	/* create F_SMP and F_SM */
        dist[0]=0;
        dist[1]=M;
        output_dist=Paso_Distribution_alloc(mpi_info, dist,1,0);
        dist[1]=N;
        input_dist=Paso_Distribution_alloc(mpi_info, dist,1,0);
        mainPattern=Paso_Pattern_alloc(PATTERN_FORMAT_DEFAULT,M,row_ptr,row_ind);
        couplePattern=Paso_Pattern_alloc(PATTERN_FORMAT_DEFAULT,M,NULL,NULL);
        send=Paso_SharedComponents_alloc(0,NULL,NULL,NULL,1,0,mpi_info);
        coupler=Paso_Coupler_alloc(send,send);
        pattern=Paso_SystemMatrixPattern_alloc(PATTERN_FORMAT_DEFAULT,output_dist,input_dist,
                                               mainPattern,couplePattern,coupler);

 	out = Paso_SystemMatrix_alloc(MATRIX_FORMAT_DEFAULT, pattern, 1, 1);
	/* copy values and cleanup temps */
	for( i=0; i<nz; i++ ) out->mainBlock->val[i] = val[i];

        Paso_SystemMatrixPattern_free(pattern);
        Paso_Pattern_free(mainPattern);
        Paso_Pattern_free(couplePattern);
        Paso_Coupler_free(coupler);
        Paso_Distribution_free(output_dist);
        Paso_Distribution_free(input_dist);
	Paso_SharedComponents_free(send);
        Paso_MPIInfo_free(mpi_info);
	MEMFREE( val );
	MEMFREE( row_ind );

	return out;
}

Paso_SystemMatrix* Paso_SystemMatrix_loadMM_toCSC( char *fileName_p )
{
        index_t dist[2];
        Paso_Distribution* input_dist=NULL, *output_dist=NULL;
	FILE *fileHandle_p = NULL;
        Paso_Pattern* mainPattern=NULL, *couplePattern=NULL;
	Paso_SystemMatrixPattern *pattern = NULL;
	Paso_SystemMatrix *out = NULL;
        Paso_SharedComponents *send =NULL;
        Paso_Coupler *coupler=NULL;
	index_t *col_ind = NULL;
	index_t *row_ind = NULL;
	index_t *col_ptr = NULL;
	double *val = NULL;
	int i, curr_col=0;
	MM_typecode matrixCode;
        Paso_MPIInfo* mpi_info=Paso_MPIInfo_alloc( MPI_COMM_WORLD);
        if (mpi_info->size >1) {
		Paso_setError(IO_ERROR, "Paso_SystemMatrix_loadMM_toCSC: support single processor only");
		return NULL;
        }

	Paso_resetError();

	/* open the file */
	fileHandle_p = fopen( fileName_p, "r" );
	if( fileHandle_p == NULL )
	{
		Paso_setError(IO_ERROR,"Paso_SystemMatrix_loadMM_toCSC: File could not be opened for reading");
                Paso_MPIInfo_free(mpi_info);
		return NULL;
	}

	/* process banner */
	if( mm_read_banner(fileHandle_p, &matrixCode) != 0 )
	{
		Paso_setError(IO_ERROR,"Paso_SystemMatrix_loadMM_toCSC: Error processing MM banner");
		fclose( fileHandle_p );
                Paso_MPIInfo_free(mpi_info);
		return NULL;
	}
	if( !(mm_is_real(matrixCode) && mm_is_sparse(matrixCode) && mm_is_general(matrixCode)) )
	{
		Paso_setError(TYPE_ERROR,"Paso_SystemMatrix_loadMM_toCSC: found Matrix Market type is not supported.");
		fclose( fileHandle_p );
                Paso_MPIInfo_free(mpi_info);
		return NULL;
	}

	/* get matrix size */
	if( mm_read_mtx_crd_size(fileHandle_p, &M, &N, &nz) != 0 )
	{
		Paso_setError(TYPE_ERROR,"Paso_SystemMatrix_loadMM_toCSC: found Matrix Market type is not supported.");
		fclose( fileHandle_p );
                Paso_MPIInfo_free(mpi_info);
		return NULL;
	}

	/* prepare storage */
	col_ind = MEMALLOC( nz, index_t );
	row_ind = MEMALLOC( nz, index_t );
	val = MEMALLOC( nz, double );

	col_ptr = MEMALLOC( (N+1), index_t );


	/* perform actual read of elements */
	for( i=0; i<nz; i++ )
	{
		fscanf( fileHandle_p, "%d %d %le\n", &row_ind[i], &col_ind[i], &val[i] );
		row_ind[i]--;
		col_ind[i]--;
	}
	fclose( fileHandle_p );

	/* sort the entries */
	q_sort( col_ind, row_ind, val, 0, nz );

	/* setup row_ptr */
	for( i=0; (i<nz && curr_col<N); curr_col++ )
	{
		while( col_ind[i] != curr_col )
			i++;
		col_ptr[curr_col] = i;
	}
	col_ptr[N] = nz;

	/* create F_SMP and F_SM */
        dist[0]=0;
        dist[1]=N;
        output_dist=Paso_Distribution_alloc(mpi_info, dist,1,0);
        dist[1]=M;
        input_dist=Paso_Distribution_alloc(mpi_info, dist,1,0);
        mainPattern=Paso_Pattern_alloc(PATTERN_FORMAT_DEFAULT,N,col_ptr,col_ind);
        couplePattern=Paso_Pattern_alloc(PATTERN_FORMAT_DEFAULT,N,NULL,NULL);
        send=Paso_SharedComponents_alloc(0,NULL,NULL,NULL,1,0,mpi_info);
        coupler=Paso_Coupler_alloc(send,send);
        pattern=Paso_SystemMatrixPattern_alloc(PATTERN_FORMAT_DEFAULT,output_dist,input_dist,
                                               mainPattern,couplePattern,coupler);
 	out = Paso_SystemMatrix_alloc(MATRIX_FORMAT_CSC, pattern, 1, 1);
	/* copy values and cleanup temps */
	for( i=0; i<nz; i++ )
		out->mainBlock->val[i] = val[i];

        Paso_SystemMatrixPattern_free(pattern);
        Paso_Pattern_free(mainPattern);
        Paso_Pattern_free(couplePattern);
        Paso_Coupler_free(coupler);
        Paso_Distribution_free(output_dist);
        Paso_Distribution_free(input_dist);
	Paso_SharedComponents_free(send);
        Paso_MPIInfo_free(mpi_info);
	MEMFREE( val );
	MEMFREE( col_ind );
	return out;
}
