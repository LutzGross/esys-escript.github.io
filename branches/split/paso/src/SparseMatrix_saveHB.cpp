
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

/* Paso: SparseMatrix saving to Harwell-Boeing format         */

/****************************************************************************/

/* Copyright: ACcESS Australia 2005                           */
/* Author: imran@esscc.uq.edu.au                              */

/****************************************************************************/

#include "Paso.h"
#include "SparseMatrix.h"

namespace paso {

/* TODO: Refactor the stuff in here into hbio, just like mmio! */

static dim_t M, N, nz;

static int calc_digits(int);
static void fmt_str(int, int, int*, int*, int*, char*, char*);
static void print_data(FILE*, int, int, int, char*, const void*, int, int);
static void generate_HB(FILE*, dim_t*, dim_t*, const double*);

/* function to get number of digits in an integer */
int calc_digits(int var)
{
    int digits = 1;
    while( (var/=10) )
        digits++;

    return digits;
}

/* function to generate the format string.
 *
 * use maxlen to determine no. of entries per line
 * use nvalues to determine no. of lines
 */
void fmt_str(int nvalues, int integer, int *width, int *nlines, int *nperline, char *pfmt, char *fmt)
{
    int per_line;
    int maxlen = *width;

    if( integer && maxlen < 10 )
        maxlen = 10;
    else
        maxlen = 13;

    per_line = 80 / maxlen;
    *nlines = nvalues / per_line;
    if( nvalues % per_line )
        (*nlines)++;

    *nperline = per_line;
    if (integer) {
        sprintf( pfmt, "(%dI%d)", per_line, maxlen );
        sprintf( fmt, "%%%dd", maxlen );
    } else {
        sprintf( pfmt, "(1P%dE%d.6)", per_line, maxlen );
        sprintf( fmt, "%%%d.6E", maxlen );
    }
    *width = maxlen;
}

/* function to print the actual data in the right format */
void print_data(FILE *fp, int n_perline, int width, int nval, char *fmt,
                const void *ptr, int integer, int adjust)
{
    int entries_done = 0;
    int padding, i;
    char pad_fmt[10];

    padding = 80 - n_perline*width;
    sprintf(pad_fmt, "%%%dc", padding);

    if (adjust != 1)
        adjust = 0;

    if (integer) {
        const dim_t *data = reinterpret_cast<const dim_t*>(ptr);
        for(i=0; i<nval; i++) {
            fprintf(fp, fmt, data[i]+adjust);
            entries_done++;
            if (entries_done == n_perline) {
                if (padding)
                    fprintf(fp, pad_fmt, ' ');
                    fprintf(fp, "\n");
                    entries_done = 0;
            }
        }
    } else {
        const double *data = reinterpret_cast<const double*>(ptr);
        for (i=0; i<nval; i++) {
            fprintf(fp, fmt, data[i]);
            entries_done++;
            if (entries_done == n_perline) {
                if (padding)
                    fprintf(fp, pad_fmt, ' ');
                fprintf(fp, "\n");
                entries_done = 0;
            }
        }
    }
    if ( entries_done ) {
        sprintf(pad_fmt, "%%%dc\n", (80 - entries_done*width));
        fprintf(fp, pad_fmt, ' ');
    }
}

void generate_HB(FILE *fp, dim_t *col_ptr, dim_t *row_ind, const double *val)
{
    char buffer[81];

    int val_lines, ind_lines, ptr_lines;
    int val_perline, ind_perline, ptr_perline;
    int val_width, ind_width, ptr_width;
    char ptr_pfmt[7], ind_pfmt[7], val_pfmt[11];
    char ptr_fmt[10], ind_fmt[10], val_fmt[10];

    /* line 1 */
    sprintf( buffer, "%-72s%-8s", "Matrix Title", "Key" );
    buffer[80] = '\0';
    fprintf( fp, "%s\n", buffer );

    /* line 2 */
    ptr_width = calc_digits( nz+1 );
    fmt_str( N+1, 1, &ptr_width, &ptr_lines, &ptr_perline, ptr_pfmt, ptr_fmt );
    ind_width = calc_digits( N );
    fmt_str( nz, 1, &ind_width, &ind_lines, &ind_perline, ind_pfmt, ind_fmt );
    val_width = 13;
    fmt_str( nz, 0, &val_width, &val_lines, &val_perline, val_pfmt, val_fmt );
    sprintf( buffer, "%14d%14d%14d%14d%14d%10c", (ptr_lines+ind_lines+val_lines), ptr_lines, ind_lines, val_lines, 0, ' ' );
    buffer[80] = '\0';
    fprintf( fp, "%s\n", buffer );

    /* line 3 */
    sprintf( buffer, "%c%c%c%11c%14d%14d%14d%14d%10c", 'R', 'U', 'A', ' ', M, N, nz, 0, ' ' );
    buffer[80] = '\0';
    fprintf( fp, "%s\n", buffer );

    /* line 4 */
    sprintf( buffer, "%16s%16s%20s%28c", ptr_pfmt, ind_pfmt, val_pfmt, ' ');
    buffer[80]='\0';
    fprintf( fp, "%s\n", buffer );

    /* line 5 */
    /* NOT PRESENT */

    /* write the actual data */
    print_data(fp, ptr_perline, ptr_width, (N+1), ptr_fmt, col_ptr, 1, 1);
    print_data(fp, ind_perline, ind_width, nz, ind_fmt, row_ind, 1, 1);
    print_data(fp, val_perline, val_width, nz, val_fmt, val, 0, 0);
}

void SparseMatrix::saveHB_CSC(const char* filename) const
{
    FILE* fileHandle = fopen(filename, "w");
    if (!fileHandle) {
        Esys_setError(IO_ERROR, "SparseMatrix::saveHB_CSC: File could not be opened for writing.");
        return;
    }

    int i, curr_col,j ;
    int iPtr, iCol, ir, ic;
    const index_t index_offset=(type & MATRIX_FORMAT_OFFSET1 ? 1:0);
    dim_t *row_ind=NULL, *col_ind = NULL, *col_ptr=NULL;

    if (row_block_size == 1 && col_block_size == 1) {
        M = numRows;
        N = numCols;
        generate_HB(fileHandle, pattern->ptr, pattern->index, val);
    } else {
        M = numRows*row_block_size;
        N = numCols*col_block_size;

        row_ind = new dim_t[len];
        col_ind = new dim_t[len];

        i = 0;
        for (iCol=0; iCol<pattern->numOutput; iCol++)
            for (ic=0; ic<col_block_size; ic++)
                for (iPtr=pattern->ptr[iCol]-index_offset; iPtr<pattern->ptr[iCol+1]-index_offset; iPtr++)
                    for (ir=0; ir<row_block_size; ir++) {
                        row_ind[i] = (pattern->index[iPtr]-index_offset)*row_block_size+ir+1;
                        col_ind[i] = iCol*col_block_size+ic+1;
                        i++;
                    }
        /* get the col_ptr */
        col_ptr = new dim_t[(N+1)];

        curr_col = 0;
        for (j=0; (j<len && curr_col<N); curr_col++) {
            while( col_ind[j] != curr_col) j++;
            col_ptr[curr_col] = j;
        }
        col_ptr[N] = len;

        /* generate the HB file */
        generate_HB(fileHandle, col_ptr, row_ind, val);

        /* free the allocated memory */
        delete[] col_ptr;
        delete[] col_ind;
        delete[] row_ind;
    }
    fclose(fileHandle);
}

} // namespace paso

