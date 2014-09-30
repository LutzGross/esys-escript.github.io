
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/


/****************************************************************************/

/* Paso: Matrix Market format is loaded to a SystemMatrix     */

/****************************************************************************/

/* Copyrights by ACcESS Australia 2003,2004,2005 */
/* Author: imran@access.edu.au */

/****************************************************************************/

#include "Paso.h"
#include "mmio.h"
#include "SystemMatrix.h"

#include "limits.h"

namespace paso {

static void swap(index_t*, index_t*, double*, int, int);
static void q_sort(index_t*, index_t*, double*, int, int);

static int M, N, nz;

void swap(index_t *r, index_t *c, double *v, int left, int right)
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

void q_sort(index_t *row, index_t *col, double *val, int begin, int end)
{
    int l, r;
    int flag;

    if (end > begin) {
        l = begin + 1;
        r = end;

        while(l < r) {
            /* This whole section is for checking lval<pivot, where
            pivot=N*row[begin]+col[begin] and lval=N*row[l]+col[l]. */
            if (row[l]<row[begin]) {
                if (ABS(row[l]-row[begin])==1 && ABS(col[l]-col[begin])==N)
                    flag=0;
                else
                    flag=1;
            } else if (row[l]==row[begin]) {
                if (col[l]<col[begin])
                    flag=1;
                else
                    flag=0;
            } else {
                if (ABS(row[l]-row[begin])==1 && ABS(col[l]-col[begin])==N)
                    flag=1;
                else
                    flag=0;
            }

            if (flag==1) {
                l++;
            } else {
                r--;
                swap(row, col, val, l, r);
            }
        }
        l--;
        swap(row, col, val, begin, l);
        q_sort(row, col, val, begin, l);
        q_sort(row, col, val, r, end);
    }
}

SystemMatrix_ptr SystemMatrix::loadMM_toCSR(const char *filename)
{
    index_t *col_ind = NULL;
    index_t *row_ind = NULL;
    index_t *row_ptr = NULL;
    double *val = NULL;
    SystemMatrix_ptr out;
    int curr_row;
    MM_typecode matrixCode;
    esysUtils::JMPI mpi_info=esysUtils::makeInfo(MPI_COMM_WORLD);
    Esys_resetError();
    if (mpi_info->size > 1) {
        Esys_setError(IO_ERROR, "SystemMatrix::loadMM_toCSR: supports single processor only");
        return out;
    }

    // open the file
    std::ifstream f(filename);
    if (!f.good()) {
        Esys_setError(IO_ERROR, "SystemMatrix::loadMM_toCSR: Cannot open file for reading.");
        return out;
    }

    // process banner
    if (mm_read_banner(f, &matrixCode) != 0) {
        Esys_setError(IO_ERROR, "SystemMatrix::loadMM_toCSR: Error processing MM banner.");
        f.close();
        return out;
    }
    if ( !(mm_is_real(matrixCode) && mm_is_sparse(matrixCode) && mm_is_general(matrixCode)) ) {
        Esys_setError(TYPE_ERROR, "SystemMatrix::loadMM_toCSR: found Matrix Market type is not supported.");
        f.close();
        return out;
    }

    // get matrix size
    if (mm_read_mtx_crd_size(f, &M, &N, &nz) != 0) {
        Esys_setError(IO_ERROR, "SystemMatrix::loadMM_toCSR: Could not read sparse matrix size.");
        f.close();
        return out;
    }

    // prepare storage
    col_ind = new index_t[nz];
    row_ind = new index_t[nz];
    val = new double[nz];

    row_ptr = new index_t[M+1];

    // perform actual read of elements
    for (int i=0; i<nz; i++) {
        f >> row_ind[i] >> col_ind[i] >> val[i];
        if (!f.good()) {
            delete[] val;
            delete[] row_ind;
            delete[] col_ind;
            delete[] row_ptr;
            f.close();
            return out;
        }
        row_ind[i]--;
        col_ind[i]--;
    }
    f.close();

    // sort the entries
    q_sort(row_ind, col_ind, val, 0, nz);

    // setup row_ptr
    curr_row = 0;
    for (int i=0; (i<nz && curr_row<M); curr_row++) {
        while (row_ind[i] != curr_row) {
            i++;
        }
        row_ptr[curr_row] = i;
    }
    row_ptr[M] = nz;

    // create return value
    index_t dist[2];
    dist[0]=0;
    dist[1]=M;
    Distribution_ptr output_dist(new Distribution(mpi_info, dist, 1, 0));
    dist[1]=N;
    Distribution_ptr input_dist(new Distribution(mpi_info, dist, 1, 0));
    Pattern_ptr mainPattern(new Pattern(MATRIX_FORMAT_DEFAULT, M, N, row_ptr, col_ind));
    Pattern_ptr couplePattern(new Pattern(MATRIX_FORMAT_DEFAULT, M, N, NULL, NULL));
    dist[0]=M;
    SharedComponents_ptr send(new SharedComponents(
                                    M, 0, NULL, NULL, dist, 1, 0, mpi_info));
    Connector_ptr connector(new Connector(send, send));
    SystemMatrixPattern_ptr pattern(new SystemMatrixPattern(
                MATRIX_FORMAT_DEFAULT, output_dist, input_dist, mainPattern,
                couplePattern, couplePattern, connector, connector));
    out.reset(new SystemMatrix(MATRIX_FORMAT_DEFAULT, pattern, 1, 1, true));

    // copy values
#pragma omp parallel for
    for (int i=0; i<nz; i++)
        out->mainBlock->val[i] = val[i];

    delete[] val;
    delete[] row_ind;
    return out;
}

SystemMatrix_ptr SystemMatrix::loadMM_toCSC(const char* filename)
{
    Pattern_ptr mainPattern, couplePattern;
    SystemMatrixPattern_ptr pattern;
    SystemMatrix_ptr out;
    Connector_ptr connector;
    index_t *col_ind = NULL;
    index_t *row_ind = NULL;
    index_t *col_ptr = NULL;
    double *val = NULL;
    int curr_col=0;
    MM_typecode matrixCode;
    esysUtils::JMPI mpi_info=esysUtils::makeInfo( MPI_COMM_WORLD);
    if (mpi_info->size > 1) {
        Esys_setError(IO_ERROR, "SystemMatrix::loadMM_toCSC: supports single processor only");
        return out;
    }

    Esys_resetError();

    // open the file
    std::ifstream f(filename);
    if (!f.good()) {
        Esys_setError(IO_ERROR, "SystemMatrix::loadMM_toCSC: File could not be opened for reading.");
        return out;
    }

    // process banner
    if (mm_read_banner(f, &matrixCode) != 0) {
        Esys_setError(IO_ERROR,"SystemMatrix::loadMM_toCSC: Error processing MM banner.");
        f.close();
        return out;
    }
    if( !(mm_is_real(matrixCode) && mm_is_sparse(matrixCode) && mm_is_general(matrixCode)) ) {
        Esys_setError(TYPE_ERROR, "SystemMatrix::loadMM_toCSC: found Matrix Market type is not supported.");
        f.close();
        return out;
    }

    // get matrix size
    if (mm_read_mtx_crd_size(f, &M, &N, &nz) != 0) {
        Esys_setError(TYPE_ERROR, "SystemMatrix::loadMM_toCSC: found Matrix Market type is not supported.");
        f.close();
        return out;
    }

    // prepare storage
    col_ind = new index_t[nz];
    row_ind = new index_t[nz];
    val = new double [nz];
    col_ptr = new index_t[N+1];

    // perform actual read of elements
    for (int i=0; i<nz; i++) {
        f >> row_ind[i] >> col_ind[i] >> val[i];
        if (!f.good()) {
            delete[] val;
            delete[] row_ind;
            delete[] col_ind;
            delete[] col_ptr;
            f.close();
            return out;
        }
        row_ind[i]--;
        col_ind[i]--;
    }
    f.close();

    // sort the entries
    q_sort( col_ind, row_ind, val, 0, nz );

    // setup row_ptr
    for(int i=0; (i<nz && curr_col<N); curr_col++) {
        while( col_ind[i] != curr_col )
            i++;
        col_ptr[curr_col] = i;
    }
    col_ptr[N] = nz;

    index_t dist[2];
    dist[0]=0;
    dist[1]=N;
    Distribution_ptr output_dist(new Distribution(mpi_info, dist,1,0));
    dist[1]=M;
    Distribution_ptr input_dist(new Distribution(mpi_info, dist,1,0));
    mainPattern.reset(new Pattern(MATRIX_FORMAT_DEFAULT,N,M,col_ptr,col_ind));
    couplePattern.reset(new Pattern(MATRIX_FORMAT_DEFAULT,N,M,NULL,NULL));
    SharedComponents_ptr send(new SharedComponents(
                    N, 0, NULL, NULL, NULL, 1, 0, mpi_info));
    connector.reset(new Connector(send,send));
    pattern.reset(new SystemMatrixPattern(MATRIX_FORMAT_DEFAULT,
                output_dist, input_dist, mainPattern, couplePattern,
                couplePattern, connector, connector));
    out.reset(new SystemMatrix(MATRIX_FORMAT_CSC, pattern, 1, 1, true));

    // copy values
#pragma omp parallel for
    for (int i=0; i<nz; i++)
        out->mainBlock->val[i] = val[i];

    delete[] val;
    delete[] row_ind;
    return out;
}

void RHS_loadMM_toCSR(const char *filename, double *b, dim_t size)
{
    MM_typecode matrixCode;
    Esys_resetError();
    // open the file
    std::ifstream f(filename);
    if (!f.good()) {
        Esys_setError(IO_ERROR, "RHS_loadMM_toCSR: Cannot open file for reading.");
    }

    // process banner
    if (mm_read_banner(f, &matrixCode) != 0) {
        Esys_setError(IO_ERROR, "RHS_loadMM_toCSR: Error processing MM banner.");
    }
    if( !(mm_is_real(matrixCode) && mm_is_general(matrixCode) && mm_is_array(matrixCode)) ) {
        Esys_setError(TYPE_ERROR,"RHS_loadMM_toCSR: found Matrix Market type is not supported.");
    }

    // get matrix size
    if (mm_read_mtx_array_size(f, &M, &N) != 0) {
        Esys_setError(IO_ERROR, "RHS_loadMM_toCSR: Could not read sparse matrix size.");
    }

    if (M != size) {
        Esys_setError(IO_ERROR, "RHS_loadMM_toCSR: Actual and provided sizes do not match.");
    }

    if (Esys_noError()) {
        nz=M;
        // perform actual read of elements
        for (int i=0; i<nz; i++) {
            f >> b[i];
            if (!f.good()) {
                f.close();
                Esys_setError(IO_ERROR, "RHS_loadMM_toCSR: Could not read some of the values.");
            }
        }
    }
    f.close();
}

} // namespace paso

