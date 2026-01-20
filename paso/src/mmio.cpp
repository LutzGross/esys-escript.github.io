
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
*
*****************************************************************************/


/*
*   Matrix Market I/O library for ANSI C
*
*   See http://math.nist.gov/MatrixMarket for details.
*
*   (Version 1.01, 5/2003)
*/

#include "mmio.h"

#include <cstring>
#include <iostream>


int mm_read_unsymmetric_sparse(const char* fname, int* M_, int* N_, int* nz_,
                double** val_, int** I_, int** J_)
{
    MM_typecode matcode;
    int M, N, nz;
    double *val;
    int *Ip, *Jp;

    std::ifstream f(fname);
    if (!f.good())
        return -1;

    if (mm_read_banner(f, &matcode) != 0) {
        std::cerr << "mm_read_unsymmetric_sparse: Could not process "
            "Matrix Market banner in file " << fname << std::endl;
        return -1;
    }

    if ( !(mm_is_real(matcode) && mm_is_matrix(matcode) &&
            mm_is_sparse(matcode))) {
        std::cerr << "Sorry, this application does not support "
                     "Matrix Market type: " << mm_typecode_to_str(matcode)
                     << std::endl;
        return -1;
    }

    // find out size of sparse matrix: M, N, nz ....
    if (mm_read_mtx_crd_size(f, &M, &N, &nz) != 0) {
        std::cerr << "mm_read_unsymmetric_sparse: Could not parse matrix size."
            << std::endl;
        return -1;
    }

    // reserve memory for matrices
    Ip = new int[nz];
    Jp = new int[nz];
    val = new double[nz];

    for (int i=0; i<nz; i++) {
        f >> Ip[i] >> Jp[i] >> val[i];
        if (!f.good()) {
            delete[] Ip;
            delete[] Jp;
            delete[] val;
            f.close();
            return -1;
        }
        Ip[i]--;  /* adjust from 1-based to 0-based */
        Jp[i]--;
    }

    f.close();
    *M_ = M;
    *N_ = N;
    *nz_ = nz;
    *val_ = val;
    *I_ = Ip;
    *J_ = Jp;
    return 0;
}

int mm_is_valid(MM_typecode matcode)
{
    if (!mm_is_matrix(matcode)) return 0;
    if (mm_is_dense(matcode) && mm_is_pattern(matcode)) return 0;
    if (mm_is_real(matcode) && mm_is_hermitian(matcode)) return 0;
    if (mm_is_pattern(matcode) && (mm_is_hermitian(matcode) ||
                mm_is_skew(matcode))) return 0;
    return 1;
}

int mm_read_banner(std::istream& f, MM_typecode* matcode)
{
    char line[MM_MAX_LINE_LENGTH];
    char banner[MM_MAX_TOKEN_LENGTH];
    char mtx[MM_MAX_TOKEN_LENGTH];
    char crd[MM_MAX_TOKEN_LENGTH];
    char data_type[MM_MAX_TOKEN_LENGTH];
    char storage_scheme[MM_MAX_TOKEN_LENGTH];
    char *p;

    mm_clear_typecode(matcode);

    f.get(line, MM_MAX_LINE_LENGTH);
    if (!f.good())
        return MM_PREMATURE_EOF;

    if (sscanf(line, "%s %s %s %s %s", banner, mtx, crd, data_type, storage_scheme) != 5)
        return MM_PREMATURE_EOF;

    // convert to lower case
    for (p=mtx; *p!='\0'; *p=(char) tolower(*p), p++);
    for (p=crd; *p!='\0'; *p=(char) tolower(*p), p++);
    for (p=data_type; *p!='\0'; *p=(char) tolower(*p), p++);
    for (p=storage_scheme; *p!='\0'; *p=(char) tolower(*p), p++);

    // check for banner
    if (strncmp(banner, MatrixMarketBanner, strlen(MatrixMarketBanner)) != 0)
        return MM_NO_HEADER;

    // first field should be "mtx"
    if (strcmp(mtx, MM_MTX_STR) != 0)
        return MM_UNSUPPORTED_TYPE;
    mm_set_matrix(matcode);

    // second field describes whether this is a sparse matrix (in coordinate
    // storage) or a dense array
    if (strcmp(crd, MM_SPARSE_STR) == 0)
        mm_set_sparse(matcode);
    else if (strcmp(crd, MM_DENSE_STR) == 0)
        mm_set_dense(matcode);
    else
        return MM_UNSUPPORTED_TYPE;

    // third field
    if (strcmp(data_type, MM_REAL_STR) == 0)
        mm_set_real(matcode);
    else if (strcmp(data_type, MM_COMPLEX_STR) == 0)
        mm_set_complex(matcode);
    else if (strcmp(data_type, MM_PATTERN_STR) == 0)
        mm_set_pattern(matcode);
    else if (strcmp(data_type, MM_INT_STR) == 0)
        mm_set_integer(matcode);
    else
        return MM_UNSUPPORTED_TYPE;

    // fourth field
    if (strcmp(storage_scheme, MM_GENERAL_STR) == 0)
        mm_set_general(matcode);
    else if (strcmp(storage_scheme, MM_SYMM_STR) == 0)
        mm_set_symmetric(matcode);
    else if (strcmp(storage_scheme, MM_HERM_STR) == 0)
        mm_set_hermitian(matcode);
    else if (strcmp(storage_scheme, MM_SKEW_STR) == 0)
        mm_set_skew(matcode);
    else
        return MM_UNSUPPORTED_TYPE;

    return 0;
}

int mm_write_mtx_crd_size(std::ostream& f, int M, int N, int nz)
{
    f << M << " " << N << " " << nz << std::endl;
    if (!f.good())
        return MM_COULD_NOT_WRITE_FILE;
    return 0;
}

int mm_read_mtx_crd_size(std::istream& f, int* M, int* N, int* nz)
{
    char line[MM_MAX_LINE_LENGTH];

    // initialize return values, in case we exit with errors
    *M = *N = *nz = 0;

    // now continue scanning until end of comments
    do {
        f.get(line, MM_MAX_LINE_LENGTH);
        if (!f.good())
            return MM_PREMATURE_EOF;
    } while (line[0] == '%');

    // line[] is either blank or has M, N, nz
    if (sscanf(line, "%d %d %d", M, N, nz) == 3)
        return 0;
    else {
        int num_items_read;
        do {
            f.get(line, MM_MAX_LINE_LENGTH);
            if (!f.good()) return MM_PREMATURE_EOF;
            num_items_read = sscanf(line, "%d %d %d", M, N, nz);
        } while (num_items_read != 3);
    }

    return 0;
}

int mm_read_mtx_array_size(std::istream& f, int* M, int* N)
{
    char line[MM_MAX_LINE_LENGTH];

    // initialize return values, in case we exit with errors
    *M = *N = 0;

    // now continue scanning until end of comments
    do {
        f.get(line, MM_MAX_LINE_LENGTH);
        if (!f.good())
            return MM_PREMATURE_EOF;
    } while (line[0] == '%');

    // line[] is either blank or has M, N
    if (sscanf(line, "%d %d", M, N) == 2)
        return 0;
    else { // we have a blank line
        int num_items_read;
        do {
            f.get(line, MM_MAX_LINE_LENGTH);
            if (!f.good()) return MM_PREMATURE_EOF;
            num_items_read = sscanf(line, "%d %d", M, N);
        } while (num_items_read != 2);
    }

    return 0;
}

int mm_write_mtx_array_size(std::ostream& f, int M, int N)
{
    f << M << " " << N << std::endl;
    if (!f.good())
        return MM_COULD_NOT_WRITE_FILE;
    return 0;
}


/*-------------------------------------------------------------------------*/

/*****************************************************************************/
/* use when Ip[], Jp[], and val[] are already allocated                      */
/*****************************************************************************/

int mm_read_mtx_crd_data(std::istream& f, int M, int N, int nz, int* Ip,
                         int* Jp, double* val, MM_typecode matcode)
{
    int i;
    if (mm_is_complex(matcode)) {
        for (i=0; i<nz; i++) {
            f >> Ip[i] >> Jp[i] >> val[2*i] >> val[2*i+1];
            if (!f.good())
                return MM_PREMATURE_EOF;
        }
    } else if (mm_is_real(matcode)) {
        for (i=0; i<nz; i++) {
            f >> Ip[i] >> Jp[i] >> val[i];
            if (!f.good())
                return MM_PREMATURE_EOF;
        }
    } else if (mm_is_pattern(matcode)) {
        for (i=0; i<nz; i++) {
            f >> Ip[i] >> Jp[i];
            if (!f.good())
                return MM_PREMATURE_EOF;
        }
    } else
        return MM_UNSUPPORTED_TYPE;

    return 0;
}

int mm_read_mtx_crd_entry(std::istream& f, int* Ip, int* Jp, double* real,
                          double* imag, MM_typecode matcode)
{
    if (mm_is_complex(matcode)) {
        f >> *Ip >> *Jp >> *real >> *imag;
        if (!f.good())
            return MM_PREMATURE_EOF;
    } else if (mm_is_real(matcode)) {
        f >> *Ip >> *Jp >> *real;
        if (!f.good())
            return MM_PREMATURE_EOF;
    } else if (mm_is_pattern(matcode)) {
        f >> *Ip >> *Jp;
        if (!f.good())
            return MM_PREMATURE_EOF;
    } else
        return MM_UNSUPPORTED_TYPE;

    return 0;
}


/****************************************************************************
    mm_read_mtx_crd()  fills M, N, nz, array of values, and return
                       type code, e.g. 'MCRS'

                       if matrix is complex, values[] is of size 2*nz,
                           (nz pairs of real/imaginary values)
*****************************************************************************/

int mm_read_mtx_crd(char *fname, int *M, int *N, int *nz, int **Ip, int **Jp,
        double **val, MM_typecode *matcode)
{
    int ret_code;

    std::ifstream f(fname);

    if (!f.good())
        return MM_COULD_NOT_READ_FILE;

    if ((ret_code = mm_read_banner(f, matcode)) != 0)
        return ret_code;

    if (!(mm_is_valid(*matcode) && mm_is_sparse(*matcode) &&
            mm_is_matrix(*matcode)))
        return MM_UNSUPPORTED_TYPE;

    if ((ret_code = mm_read_mtx_crd_size(f, M, N, nz)) != 0)
        return ret_code;

    *Ip = new int[*nz];
    *Jp = new int[*nz];
    *val = NULL;

    if (mm_is_complex(*matcode)) {
        *val = new double[*nz * 2];
        ret_code = mm_read_mtx_crd_data(f, *M, *N, *nz, *Ip, *Jp, *val,
                *matcode);
        if (ret_code != 0) return ret_code;
    } else if (mm_is_real(*matcode)) {
        *val = new  double[*nz];
        ret_code = mm_read_mtx_crd_data(f, *M, *N, *nz, *Ip, *Jp, *val,
                *matcode);
        if (ret_code != 0) return ret_code;
    } else if (mm_is_pattern(*matcode)) {
        ret_code = mm_read_mtx_crd_data(f, *M, *N, *nz, *Ip, *Jp, *val,
                *matcode);
        if (ret_code != 0) return ret_code;
    }

    f.close();
    return 0;
}

int mm_write_banner(std::ostream& f, MM_typecode matcode)
{
    f << MatrixMarketBanner << " " << mm_typecode_to_str(matcode) << std::endl;
    if (!f.good())
        return MM_COULD_NOT_WRITE_FILE;
    return 0;
}

int mm_write_mtx_crd(char* fname, int M, int N, int nz, int* Ip, int* Jp,
                     double* val, MM_typecode matcode)
{
    FILE *f;
    int i;

    if (strcmp(fname, "stdout") == 0)
        f = stdout;
    else
    if ((f = fopen(fname, "w")) == NULL)
        return MM_COULD_NOT_WRITE_FILE;

    /* print banner followed by typecode */
    fprintf(f, "%s ", MatrixMarketBanner);
    fprintf(f, "%s\n", mm_typecode_to_str(matcode));

    /* print matrix sizes and nonzeros */
    fprintf(f, "%d %d %d\n", M, N, nz);

    /* print values */
    if (mm_is_pattern(matcode))
        for (i=0; i<nz; i++)
            fprintf(f, "%d %d\n", Ip[i], Jp[i]);
    else if (mm_is_real(matcode))
        for (i=0; i<nz; i++)
            fprintf(f, "%d %d %20.16g\n", Ip[i], Jp[i], val[i]);
    else if (mm_is_complex(matcode))
        for (i=0; i<nz; i++)
            fprintf(f, "%d %d %20.16g %20.16g\n", Ip[i], Jp[i], val[2*i],
                        val[2*i+1]);
    else {
        if (f != stdout) fclose(f);
        return MM_UNSUPPORTED_TYPE;
    }

    if (f !=stdout) fclose(f);

    return 0;
}


char *mm_typecode_to_str(MM_typecode matcode)
{
    static char buffer[MM_MAX_LINE_LENGTH];
    const char *types[4];

    /* check for MTX type */
    if (mm_is_matrix(matcode))
        types[0] = MM_MTX_STR;
    else
        return NULL;

    /* check for CRD or ARR matrix */
    if (mm_is_sparse(matcode))
        types[1] = MM_SPARSE_STR;
    else
    if (mm_is_dense(matcode))
        types[1] = MM_DENSE_STR;
    else
        return NULL;

    /* check for element data type */
    if (mm_is_real(matcode))
        types[2] = MM_REAL_STR;
    else
    if (mm_is_complex(matcode))
        types[2] = MM_COMPLEX_STR;
    else
    if (mm_is_pattern(matcode))
        types[2] = MM_PATTERN_STR;
    else
    if (mm_is_integer(matcode))
        types[2] = MM_INT_STR;
    else
        return NULL;

    /* check for symmetry type */
    if (mm_is_general(matcode))
        types[3] = MM_GENERAL_STR;
    else
    if (mm_is_symmetric(matcode))
        types[3] = MM_SYMM_STR;
    else
    if (mm_is_hermitian(matcode))
        types[3] = MM_HERM_STR;
    else
    if (mm_is_skew(matcode))
        types[3] = MM_SKEW_STR;
    else
        return NULL;

    snprintf(buffer, MM_MAX_LINE_LENGTH, "%s %s %s %s", types[0], types[1], types[2], types[3]);
    return &buffer[0];
}

/*
 * $Log$
 * Revision 1.1  2004/10/26 06:53:59  jgs
 * Initial revision
 *
 * Revision 1.1  2004/07/02 00:48:35  gross
 * matrix market io function added
 *
 *
 */
