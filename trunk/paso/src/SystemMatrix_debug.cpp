
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


/****************************************************************************

 Paso: SystemMatrix: debugging tools

*****************************************************************************

 Author: Lutz Gross, l.gross@uq.edu.au 

*****************************************************************************/

#include "SystemMatrix.h"

#include <cstring> // strcat

namespace paso {

// fills the matrix with values i+f1*j where i and j are the global row
// and column indices of the matrix entry
void SystemMatrix::fillWithGlobalCoordinates(double f1)
{
    const dim_t n = getNumRows();
    const dim_t m = getNumCols();
    const index_t me = mpi_info->rank;
    const index_t row_offset = row_distribution->first_component[me];
    const index_t col_offset = col_distribution->first_component[me];
    double* cols = new double[m];
    double* rows = new double[n];
    Coupler_ptr col_couple(new Coupler(col_coupler->connector, 1));
    Coupler_ptr row_couple(new Coupler(col_coupler->connector, 1));

#pragma omp parallel for
    for (dim_t i=0; i<n; ++i)
        rows[i]=row_offset+i;

    col_couple->startCollect(rows);
   
#pragma omp parallel for
    for (dim_t i=0; i<m; ++i)
        cols[i]=col_offset+i;

    row_couple->startCollect(cols);

    // main block
    for (dim_t q=0; q<n; ++q) {
        for (dim_t iPtr = mainBlock->pattern->ptr[q];
                iPtr < mainBlock->pattern->ptr[q+1]; ++iPtr) {
            const dim_t p = mainBlock->pattern->index[iPtr];
            for (dim_t ib=0; ib<block_size; ib++)
                mainBlock->val[iPtr*block_size+ib]=f1*rows[q]+cols[p];
        }
    }

    col_couple->finishCollect();
    if (col_coupleBlock != NULL) {
        for (dim_t q=0; q<col_coupleBlock->pattern->numOutput; ++q) {
            for (dim_t iPtr = col_coupleBlock->pattern->ptr[q];
                    iPtr < col_coupleBlock->pattern->ptr[q+1]; ++iPtr) {
                const dim_t p = col_coupleBlock->pattern->index[iPtr];
                for (dim_t ib=0; ib<block_size; ib++)
                    col_coupleBlock->val[iPtr*block_size+ib] =
                            f1*rows[q]+col_couple->recv_buffer[p];
            }
        }
    }

    row_couple->finishCollect();
    if (row_coupleBlock != NULL) {
        for (dim_t p=0; p<row_coupleBlock->pattern->numOutput; ++p) {
            for (dim_t iPtr = row_coupleBlock->pattern->ptr[p];
                    iPtr < row_coupleBlock->pattern->ptr[p+1]; ++iPtr) {
                const dim_t q = row_coupleBlock->pattern->index[iPtr];
                for (dim_t ib=0; ib<block_size; ib++)
                    row_coupleBlock->val[iPtr*block_size+ib] =
                        row_couple->recv_buffer[p]*f1+cols[q];
            }
        }      
    }
   
    delete[] cols;
    delete[] rows;
}

void SystemMatrix::print() const
{
    dim_t iPtr, q, p, ib;
    const dim_t n = getNumRows();
    const index_t rank = mpi_info->rank;
    char *str1, *str2;
    str1 = new char[n*n*block_size*30+100];
    str2 = new char[30];
   
    sprintf(str1, "rank %d Main Block:\n-----------\n", rank);
    for (q=0; q<n; ++q){
        sprintf(str2, "Row %d: ",q);
        strcat(str1, str2);
        for (iPtr=mainBlock->pattern->ptr[q]; iPtr<mainBlock->pattern->ptr[q+1]; ++iPtr) {
            sprintf(str2, "(%d ",mainBlock->pattern->index[iPtr]);
            strcat(str1, str2);
            for (ib=0; ib<block_size; ib++){
                sprintf(str2, "%f ", mainBlock->val[iPtr*block_size+ib]);
                strcat(str1, str2);
            }
            sprintf(str2, "),");
            strcat(str1, str2);
        }
        sprintf(str1, "%s\n", str1);
    }
    fprintf(stderr, "%s", str1);
    if (col_coupleBlock != NULL && mpi_info->size>1) {
        sprintf(str1, "rank %d Column Couple Block:\n------------------\n", rank);
        for (q=0; q<col_coupleBlock->pattern->numOutput; ++q) {
            sprintf(str2, "Row %d: ", q);
            strcat(str1, str2);
            for (iPtr=col_coupleBlock->pattern->ptr[q]; iPtr<col_coupleBlock->pattern->ptr[q+1]; ++iPtr) {
                if (global_id) 
                    sprintf(str2, "(%d %f),", global_id[col_coupleBlock->pattern->index[iPtr]], col_coupleBlock->val[iPtr*block_size]);
                else 
                    sprintf(str2, "(%d %f),", col_coupleBlock->pattern->index[iPtr], col_coupleBlock->val[iPtr*block_size]);
                strcat(str1, str2);
            }
            sprintf(str1, "%s\n", str1);
        }
        fprintf(stderr, "%s", str1);
    }
    if (row_coupleBlock != NULL && mpi_info->size>1) {
        sprintf(str1, "rank %d Row Couple Block:\n--------------------\n", rank);
        for (p=0; p<row_coupleBlock->pattern->numOutput; ++p) {
            sprintf(str2, "Row %d:", p);
            strcat(str1, str2);
            for (iPtr=row_coupleBlock->pattern->ptr[p]; iPtr<row_coupleBlock->pattern->ptr[p+1]; ++iPtr) {
                sprintf(str2, "(%d %f),", row_coupleBlock->pattern->index[iPtr], row_coupleBlock->val[iPtr*block_size]);
                strcat(str1, str2);
            }
            sprintf(str1, "%s\n", str1);
        }
        fprintf(stderr, "%s", str1);
    }
    if (remote_coupleBlock != NULL && mpi_info->size>1) {
        sprintf(str1, "rank %d Remote Couple Block:\n--------------------\n", rank);
        for (p=0; p<remote_coupleBlock->pattern->numOutput; ++p) {
            sprintf(str2, "Row %d:", p);
            strcat(str1, str2);
            for (iPtr=remote_coupleBlock->pattern->ptr[p]; iPtr<remote_coupleBlock->pattern->ptr[p+1]; ++iPtr) {
                sprintf(str2, "(%d %f),", remote_coupleBlock->pattern->index[iPtr], remote_coupleBlock->val[iPtr*block_size]);
                strcat(str1, str2);
            }
            sprintf(str1, "%s\n", str1);
        }
        fprintf(stderr, "%s", str1);
    }
    delete[] str1;
    delete[] str2;
}

} // namespace paso

