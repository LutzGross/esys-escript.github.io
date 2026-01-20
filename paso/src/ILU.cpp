
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


/****************************************************************************/

/* Paso: ILU preconditioner with reordering                 */

/****************************************************************************/

/* Copyrights by ACcESS Australia 2003,2004,2005              */
/* Author: Lutz Gross, l.gross@uq.edu.au                      */

/****************************************************************************/

#include "Paso.h"
#include "PasoUtil.h"
#include "Preconditioner.h"

namespace paso {

void Solver_ILU_free(Solver_ILU * in)
{
    if (in!=NULL) {
        delete[] in->factors;
        delete in;
    }
}

/// constructs the incomplete block factorization
Solver_ILU* Solver_getILU(SparseMatrix_ptr<double> A, bool verbose)
{
    const dim_t n=A->numRows;
    const dim_t n_block=A->row_block_size;
    const index_t* colorOf = A->pattern->borrowColoringPointer();
    const dim_t num_colors = A->pattern->getNumColors();
    const index_t *ptr_main = A->borrowMainDiagonalPointer();
    double A11,A12,A13,A21,A22,A23,A31,A32,A33,D;
    double S11,S12,S13,S21,S22,S23,S31,S32,S33;
    index_t i,iptr_main,iptr_ik,k,iptr_kj,j,iptr_ij,color,color2, iptr;
    Solver_ILU* out=new Solver_ILU;
    out->factors=new double[A->len];

    double time0 = escript::gettime();

#pragma omp parallel for schedule(static) private(i,iptr,k)
    for (i = 0; i < n; ++i) {
        for (iptr=A->pattern->ptr[i];iptr<A->pattern->ptr[i+1]; iptr++) {
            for (k=0;k<n_block*n_block;++k)
                out->factors[n_block*n_block*iptr+k]=A->val[n_block*n_block*iptr+k];
        }
    }

    // start factorization
    for (color=0; color<num_colors; ++color) {
        if (n_block==1) {
#pragma omp parallel for schedule(static) private(i,color2,iptr_ik,k,iptr_kj,S11,j,iptr_ij,A11,iptr_main,D)
            for (i = 0; i < n; ++i) {
                if (colorOf[i]==color) {
                    for (color2=0;color2<color;++color2) {
                        for (iptr_ik=A->pattern->ptr[i];iptr_ik<A->pattern->ptr[i+1]; ++iptr_ik) {
                            k=A->pattern->index[iptr_ik];
                            if (colorOf[k]==color2) {
                                A11=out->factors[iptr_ik];
                                /* a_ij=a_ij-a_ik*a_kj */
                                for (iptr_kj=A->pattern->ptr[k];iptr_kj<A->pattern->ptr[k+1]; iptr_kj++) {
                                    j=A->pattern->index[iptr_kj];
                                    if (colorOf[j]>color2) {
                                        S11=out->factors[iptr_kj];
                                        for (iptr_ij=A->pattern->ptr[i];iptr_ij<A->pattern->ptr[i+1]; iptr_ij++) {
                                            if (j==A->pattern->index[iptr_ij]) {
                                                out->factors[iptr_ij]-=A11*S11;
                                                break;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    iptr_main=ptr_main[i];
                    D=out->factors[iptr_main];
                    if (std::abs(D)>0.) {
                        D=1./D;
                        out->factors[iptr_main]=D;
                        /* a_ik=a_ii^{-1}*a_ik */
                        for (iptr_ik=A->pattern->ptr[i];iptr_ik<A->pattern->ptr[i+1]; ++iptr_ik) {
                            k=A->pattern->index[iptr_ik];
                            if (colorOf[k]>color) {
                                A11=out->factors[iptr_ik];
                                out->factors[iptr_ik]=A11*D;
                            }
                        }
                    } else {
                        throw PasoException("Solver_getILU: non-regular main diagonal block.");
                    }
                }
            }
        } else if (n_block==2) {
#pragma omp parallel for schedule(static) private(i,color2,iptr_ik,k,iptr_kj,S11,S21,S12,S22,j,iptr_ij,A11,A21,A12,A22,iptr_main,D)
            for (i = 0; i < n; ++i) {
                if (colorOf[i]==color) {
                    for (color2=0;color2<color;++color2) {
                        for (iptr_ik=A->pattern->ptr[i];iptr_ik<A->pattern->ptr[i+1]; ++iptr_ik) {
                            k=A->pattern->index[iptr_ik];
                            if (colorOf[k]==color2) {
                                A11=out->factors[iptr_ik*4  ];
                                A21=out->factors[iptr_ik*4+1];
                                A12=out->factors[iptr_ik*4+2];
                                A22=out->factors[iptr_ik*4+3];
                                /* a_ij=a_ij-a_ik*a_kj */
                                for (iptr_kj=A->pattern->ptr[k];iptr_kj<A->pattern->ptr[k+1]; iptr_kj++) {
                                    j=A->pattern->index[iptr_kj];
                                    if (colorOf[j]>color2) {
                                        S11=out->factors[iptr_kj*4];
                                        S21=out->factors[iptr_kj*4+1];
                                        S12=out->factors[iptr_kj*4+2];
                                        S22=out->factors[iptr_kj*4+3];
                                        for (iptr_ij=A->pattern->ptr[i];iptr_ij<A->pattern->ptr[i+1]; iptr_ij++) {
                                            if (j==A->pattern->index[iptr_ij]) {
                                                out->factors[4*iptr_ij  ]-=A11*S11+A12*S21;
                                                out->factors[4*iptr_ij+1]-=A21*S11+A22*S21;
                                                out->factors[4*iptr_ij+2]-=A11*S12+A12*S22;
                                                out->factors[4*iptr_ij+3]-=A21*S12+A22*S22;
                                                break;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    iptr_main=ptr_main[i];
                    A11=out->factors[iptr_main*4];
                    A21=out->factors[iptr_main*4+1];
                    A12=out->factors[iptr_main*4+2];
                    A22=out->factors[iptr_main*4+3];
                    D = A11*A22-A12*A21;
                    if (std::abs(D)>0.) {
                        D=1./D;
                        S11= A22*D;
                        S21=-A21*D;
                        S12=-A12*D;
                        S22= A11*D;
                        out->factors[iptr_main*4]  = S11;
                        out->factors[iptr_main*4+1]= S21;
                        out->factors[iptr_main*4+2]= S12;
                        out->factors[iptr_main*4+3]= S22;
                        /* a_ik=a_ii^{-1}*a_ik */
                        for (iptr_ik=A->pattern->ptr[i];iptr_ik<A->pattern->ptr[i+1]; ++iptr_ik) {
                            k=A->pattern->index[iptr_ik];
                            if (colorOf[k]>color) {
                                A11=out->factors[iptr_ik*4  ];
                                A21=out->factors[iptr_ik*4+1];
                                A12=out->factors[iptr_ik*4+2];
                                A22=out->factors[iptr_ik*4+3];
                                out->factors[4*iptr_ik  ]=S11*A11+S12*A21;
                                out->factors[4*iptr_ik+1]=S21*A11+S22*A21;
                                out->factors[4*iptr_ik+2]=S11*A12+S12*A22;
                                out->factors[4*iptr_ik+3]=S21*A12+S22*A22;
                            }
                        }
                    } else {
                        throw PasoException("Solver_getILU: non-regular main diagonal block.");
                    }
                }
            }
        } else if (n_block==3) {
#pragma omp parallel for schedule(static) private(i,color2,iptr_ik,k,iptr_kj,S11,S21,S31,S12,S22,S32,S13,S23,S33,j,iptr_ij,A11,A21,A31,A12,A22,A32,A13,A23,A33,iptr_main,D)
            for (i = 0; i < n; ++i) {
                if (colorOf[i]==color) {
                    for (color2=0;color2<color;++color2) {
                        for (iptr_ik=A->pattern->ptr[i];iptr_ik<A->pattern->ptr[i+1]; ++iptr_ik) {
                            k=A->pattern->index[iptr_ik];
                            if (colorOf[k]==color2) {
                                A11=out->factors[iptr_ik*9  ];
                                A21=out->factors[iptr_ik*9+1];
                                A31=out->factors[iptr_ik*9+2];
                                A12=out->factors[iptr_ik*9+3];
                                A22=out->factors[iptr_ik*9+4];
                                A32=out->factors[iptr_ik*9+5];
                                A13=out->factors[iptr_ik*9+6];
                                A23=out->factors[iptr_ik*9+7];
                                A33=out->factors[iptr_ik*9+8];
                                /* a_ij=a_ij-a_ik*a_kj */
                                for (iptr_kj=A->pattern->ptr[k];iptr_kj<A->pattern->ptr[k+1]; iptr_kj++) {
                                    j=A->pattern->index[iptr_kj];
                                    if (colorOf[j]>color2) {
                                        S11=out->factors[iptr_kj*9  ];
                                        S21=out->factors[iptr_kj*9+1];
                                        S31=out->factors[iptr_kj*9+2];
                                        S12=out->factors[iptr_kj*9+3];
                                        S22=out->factors[iptr_kj*9+4];
                                        S32=out->factors[iptr_kj*9+5];
                                        S13=out->factors[iptr_kj*9+6];
                                        S23=out->factors[iptr_kj*9+7];
                                        S33=out->factors[iptr_kj*9+8];
                                        for (iptr_ij=A->pattern->ptr[i];iptr_ij<A->pattern->ptr[i+1]; iptr_ij++) {
                                            if (j==A->pattern->index[iptr_ij]) {
                                                out->factors[iptr_ij*9  ]-=A11*S11+A12*S21+A13*S31;
                                                out->factors[iptr_ij*9+1]-=A21*S11+A22*S21+A23*S31;
                                                out->factors[iptr_ij*9+2]-=A31*S11+A32*S21+A33*S31;
                                                out->factors[iptr_ij*9+3]-=A11*S12+A12*S22+A13*S32;
                                                out->factors[iptr_ij*9+4]-=A21*S12+A22*S22+A23*S32;
                                                out->factors[iptr_ij*9+5]-=A31*S12+A32*S22+A33*S32;
                                                out->factors[iptr_ij*9+6]-=A11*S13+A12*S23+A13*S33;
                                                out->factors[iptr_ij*9+7]-=A21*S13+A22*S23+A23*S33;
                                                out->factors[iptr_ij*9+8]-=A31*S13+A32*S23+A33*S33;
                                                break;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    iptr_main=ptr_main[i];
                    A11=out->factors[iptr_main*9  ];
                    A21=out->factors[iptr_main*9+1];
                    A31=out->factors[iptr_main*9+2];
                    A12=out->factors[iptr_main*9+3];
                    A22=out->factors[iptr_main*9+4];
                    A32=out->factors[iptr_main*9+5];
                    A13=out->factors[iptr_main*9+6];
                    A23=out->factors[iptr_main*9+7];
                    A33=out->factors[iptr_main*9+8];
                    D = A11*(A22*A33-A23*A32)+ A12*(A31*A23-A21*A33)+A13*(A21*A32-A31*A22);
                    if (std::abs(D)>0.) {
                        D=1./D;
                        S11=(A22*A33-A23*A32)*D;
                        S21=(A31*A23-A21*A33)*D;
                        S31=(A21*A32-A31*A22)*D;
                        S12=(A13*A32-A12*A33)*D;
                        S22=(A11*A33-A31*A13)*D;
                        S32=(A12*A31-A11*A32)*D;
                        S13=(A12*A23-A13*A22)*D;
                        S23=(A13*A21-A11*A23)*D;
                        S33=(A11*A22-A12*A21)*D;

                        out->factors[iptr_main*9  ]=S11;
                        out->factors[iptr_main*9+1]=S21;
                        out->factors[iptr_main*9+2]=S31;
                        out->factors[iptr_main*9+3]=S12;
                        out->factors[iptr_main*9+4]=S22;
                        out->factors[iptr_main*9+5]=S32;
                        out->factors[iptr_main*9+6]=S13;
                        out->factors[iptr_main*9+7]=S23;
                        out->factors[iptr_main*9+8]=S33;

                        /* a_ik=a_ii^{-1}*a_ik */
                        for (iptr_ik=A->pattern->ptr[i];iptr_ik<A->pattern->ptr[i+1]; ++iptr_ik) {
                            k=A->pattern->index[iptr_ik];
                            if (colorOf[k]>color) {
                                A11=out->factors[iptr_ik*9  ];
                                A21=out->factors[iptr_ik*9+1];
                                A31=out->factors[iptr_ik*9+2];
                                A12=out->factors[iptr_ik*9+3];
                                A22=out->factors[iptr_ik*9+4];
                                A32=out->factors[iptr_ik*9+5];
                                A13=out->factors[iptr_ik*9+6];
                                A23=out->factors[iptr_ik*9+7];
                                A33=out->factors[iptr_ik*9+8];
                                out->factors[iptr_ik*9  ]=S11*A11+S12*A21+S13*A31;
                                out->factors[iptr_ik*9+1]=S21*A11+S22*A21+S23*A31;
                                out->factors[iptr_ik*9+2]=S31*A11+S32*A21+S33*A31;
                                out->factors[iptr_ik*9+3]=S11*A12+S12*A22+S13*A32;
                                out->factors[iptr_ik*9+4]=S21*A12+S22*A22+S23*A32;
                                out->factors[iptr_ik*9+5]=S31*A12+S32*A22+S33*A32;
                                out->factors[iptr_ik*9+6]=S11*A13+S12*A23+S13*A33;
                                out->factors[iptr_ik*9+7]=S21*A13+S22*A23+S23*A33;
                                out->factors[iptr_ik*9+8]=S31*A13+S32*A23+S33*A33;
                            }
                        }
                    } else {
                        throw PasoException("Solver_getILU: non-regular main diagonal block.");
                    }
                }
            }
        } else {
            throw PasoException("Solver_getILU: block size greater than 3 is not supported.");
        }
#pragma omp barrier
    }

    if (verbose) {
        const double time_fac=escript::gettime()-time0;
        printf("timing: ILU: coloring/elimination: %e sec\n",time_fac);
    }
    return out;
}

/****************************************************************************/

/* Applies ILU precondition b-> x

   In fact it solves LUx=b in the form x= U^{-1} L^{-1}b

   Should be called within a parallel region.
   Barrier synchronization should be performed to make sure that the input
   vector is available.
*/

void Solver_solveILU(SparseMatrix_ptr<double> A, Solver_ILU* ilu, double* x,
                     const double* b)
{
    dim_t i,k;
    index_t color,iptr_ik,iptr_main;
    double S1,S2,S3,R1,R2,R3;
    const dim_t n=A->numRows;
    const dim_t n_block=A->row_block_size;
    const index_t* colorOf = A->pattern->borrowColoringPointer();
    const dim_t num_colors = A->pattern->getNumColors();
    const index_t *ptr_main = A->borrowMainDiagonalPointer();

    /* copy x into b */
#pragma omp parallel for private(i) schedule(static)
    for (i=0;i<n*n_block;++i)
        x[i]=b[i];

    /* forward substitution */
    for (color=0;color<num_colors;++color) {
        if (n_block==1) {
#pragma omp parallel for schedule(static) private(i,iptr_ik,k,S1,R1,iptr_main)
            for (i = 0; i < n; ++i) {
                if (colorOf[i]==color) {
                    /* x_i=x_i-a_ik*x_k */
                    S1=x[i];
                    for (iptr_ik=A->pattern->ptr[i];iptr_ik<A->pattern->ptr[i+1]; ++iptr_ik) {
                        k=A->pattern->index[iptr_ik];
                        if (colorOf[k]<color) {
                            R1=x[k];
                            S1-=ilu->factors[iptr_ik]*R1;
                        }
                    }
                    iptr_main=ptr_main[i];
                    x[i]=ilu->factors[iptr_main]*S1;
                }
            }
        } else if (n_block==2) {
#pragma omp parallel for schedule(static) private(i,iptr_ik,k,iptr_main,S1,S2,R1,R2)
            for (i = 0; i < n; ++i) {
                if (colorOf[i]==color) {
                    /* x_i=x_i-a_ik*x_k */
                    S1=x[2*i];
                    S2=x[2*i+1];
                    for (iptr_ik=A->pattern->ptr[i];iptr_ik<A->pattern->ptr[i+1]; ++iptr_ik) {
                        k=A->pattern->index[iptr_ik];
                        if (colorOf[k]<color) {
                            R1=x[2*k];
                            R2=x[2*k+1];
                            S1-=ilu->factors[4*iptr_ik  ]*R1+ilu->factors[4*iptr_ik+2]*R2;
                            S2-=ilu->factors[4*iptr_ik+1]*R1+ilu->factors[4*iptr_ik+3]*R2;
                        }
                    }
                    iptr_main=ptr_main[i];
                    x[2*i  ]=ilu->factors[4*iptr_main  ]*S1+ilu->factors[4*iptr_main+2]*S2;
                    x[2*i+1]=ilu->factors[4*iptr_main+1]*S1+ilu->factors[4*iptr_main+3]*S2;
                }
            }
        } else if (n_block==3) {
#pragma omp parallel for schedule(static) private(i,iptr_ik,iptr_main,k,S1,S2,S3,R1,R2,R3)
            for (i = 0; i < n; ++i) {
                if (colorOf[i]==color) {
                    /* x_i=x_i-a_ik*x_k */
                    S1=x[3*i];
                    S2=x[3*i+1];
                    S3=x[3*i+2];
                    for (iptr_ik=A->pattern->ptr[i];iptr_ik<A->pattern->ptr[i+1]; ++iptr_ik) {
                        k=A->pattern->index[iptr_ik];
                        if (colorOf[k]<color) {
                            R1=x[3*k];
                            R2=x[3*k+1];
                            R3=x[3*k+2];
                            S1-=ilu->factors[9*iptr_ik  ]*R1+ilu->factors[9*iptr_ik+3]*R2+ilu->factors[9*iptr_ik+6]*R3;
                            S2-=ilu->factors[9*iptr_ik+1]*R1+ilu->factors[9*iptr_ik+4]*R2+ilu->factors[9*iptr_ik+7]*R3;
                            S3-=ilu->factors[9*iptr_ik+2]*R1+ilu->factors[9*iptr_ik+5]*R2+ilu->factors[9*iptr_ik+8]*R3;
                        }
                    }
                    iptr_main=ptr_main[i];
                    x[3*i  ]=ilu->factors[9*iptr_main  ]*S1+ilu->factors[9*iptr_main+3]*S2+ilu->factors[9*iptr_main+6]*S3;
                    x[3*i+1]=ilu->factors[9*iptr_main+1]*S1+ilu->factors[9*iptr_main+4]*S2+ilu->factors[9*iptr_main+7]*S3;
                    x[3*i+2]=ilu->factors[9*iptr_main+2]*S1+ilu->factors[9*iptr_main+5]*S2+ilu->factors[9*iptr_main+8]*S3;
                }
            }
        }
    }
    /* backward substitution */
    for (color=num_colors-1; color>-1; --color) {
        if (n_block==1) {
#pragma omp parallel for schedule(static) private(i,iptr_ik,k,S1,R1)
            for (i = 0; i < n; ++i) {
                if (colorOf[i]==color) {
                    /* x_i=x_i-a_ik*x_k */
                    S1=x[i];
                    for (iptr_ik=A->pattern->ptr[i];iptr_ik<A->pattern->ptr[i+1]; ++iptr_ik) {
                        k=A->pattern->index[iptr_ik];
                        if (colorOf[k]>color) {
                            R1=x[k];
                            S1-=ilu->factors[iptr_ik]*R1;
                        }
                    }
                    x[i]=S1;
                }
            }
        } else if (n_block==2) {
#pragma omp parallel for schedule(static) private(i,iptr_ik,k,S1,S2,R1,R2)
            for (i = 0; i < n; ++i) {
                if (colorOf[i]==color) {
                    /* x_i=x_i-a_ik*x_k */
                    S1=x[2*i];
                    S2=x[2*i+1];
                    for (iptr_ik=A->pattern->ptr[i];iptr_ik<A->pattern->ptr[i+1]; ++iptr_ik) {
                        k=A->pattern->index[iptr_ik];
                        if (colorOf[k]>color) {
                            R1=x[2*k];
                            R2=x[2*k+1];
                            S1-=ilu->factors[4*iptr_ik  ]*R1+ilu->factors[4*iptr_ik+2]*R2;
                            S2-=ilu->factors[4*iptr_ik+1]*R1+ilu->factors[4*iptr_ik+3]*R2;
                        }
                    }
                    x[2*i]=S1;
                    x[2*i+1]=S2;
                }
            }
        } else if (n_block==3) {
#pragma omp parallel for schedule(static) private(i,iptr_ik,k,S1,S2,S3,R1,R2,R3)
            for (i = 0; i < n; ++i) {
                if (colorOf[i]==color) {
                    /* x_i=x_i-a_ik*x_k */
                    S1=x[3*i  ];
                    S2=x[3*i+1];
                    S3=x[3*i+2];
                    for (iptr_ik=A->pattern->ptr[i];iptr_ik<A->pattern->ptr[i+1]; ++iptr_ik) {
                        k=A->pattern->index[iptr_ik];
                        if (colorOf[k]>color) {
                            R1=x[3*k];
                            R2=x[3*k+1];
                            R3=x[3*k+2];
                            S1-=ilu->factors[9*iptr_ik  ]*R1+ilu->factors[9*iptr_ik+3]*R2+ilu->factors[9*iptr_ik+6]*R3;
                            S2-=ilu->factors[9*iptr_ik+1]*R1+ilu->factors[9*iptr_ik+4]*R2+ilu->factors[9*iptr_ik+7]*R3;
                            S3-=ilu->factors[9*iptr_ik+2]*R1+ilu->factors[9*iptr_ik+5]*R2+ilu->factors[9*iptr_ik+8]*R3;
                        }
                    }
                    x[3*i]=S1;
                    x[3*i+1]=S2;
                    x[3*i+2]=S3;
                }
            }
        }
#pragma omp barrier
    }
}

} // namespace paso

