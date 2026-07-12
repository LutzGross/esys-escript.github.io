
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#include <oxley/DefaultAssembler3D.h>
#include <oxley/domainhelpers.h>

#include <escript/DataTypes.h>
#include <escript/index.h>

using namespace std;
using escript::AbstractSystemMatrix;
using escript::Data;

namespace oxley {

template<class Scalar>
void DefaultAssembler3D<Scalar>::collateFunctionSpaceTypes(
                             vector<int>& fsTypes, const DataMap& coefs) const
{
    if (isNotEmpty("A", coefs))
        fsTypes.push_back(coefs.find("A")->second.getFunctionSpace().getTypeCode());
    if (isNotEmpty("B", coefs))
        fsTypes.push_back(coefs.find("B")->second.getFunctionSpace().getTypeCode());
    if (isNotEmpty("C", coefs))
        fsTypes.push_back(coefs.find("C")->second.getFunctionSpace().getTypeCode());
    if (isNotEmpty("D", coefs))
        fsTypes.push_back(coefs.find("D")->second.getFunctionSpace().getTypeCode());
    if (isNotEmpty("X", coefs))
        fsTypes.push_back(coefs.find("X")->second.getFunctionSpace().getTypeCode());
    if (isNotEmpty("Y", coefs))
        fsTypes.push_back(coefs.find("Y")->second.getFunctionSpace().getTypeCode());
}

/****************************************************************************/
// wrappers
/****************************************************************************/

template<class Scalar>
void DefaultAssembler3D<Scalar>::assemblePDESingle(AbstractSystemMatrix* mat,
                                        Data& rhs, const DataMap& coefs) const
{
    const Data& A = unpackData("A", coefs);
    const Data& B = unpackData("B", coefs);
    const Data& C = unpackData("C", coefs);
    const Data& D = unpackData("D", coefs);
    const Data& X = unpackData("X", coefs);
    const Data& Y = unpackData("Y", coefs);
    assemblePDESingle(mat, rhs, A, B, C, D, X, Y);

}

template<class Scalar>
void DefaultAssembler3D<Scalar>::assemblePDEBoundarySingle(
                                        AbstractSystemMatrix* mat, Data& rhs,
                                        const DataMap& coefs) const 
{
    const Data& d = unpackData("d", coefs);
    const Data& y = unpackData("y", coefs);
    assemblePDEBoundarySingle(mat, rhs, d, y);
}

template<class Scalar>
void DefaultAssembler3D<Scalar>::assemblePDESingleReduced(
                                        AbstractSystemMatrix* mat, Data& rhs,
                                        const DataMap& coefs) const
{
    const Data& A = unpackData("A", coefs);
    const Data& B = unpackData("B", coefs);
    const Data& C = unpackData("C", coefs);
    const Data& D = unpackData("D", coefs);
    const Data& X = unpackData("X", coefs);
    const Data& Y = unpackData("Y", coefs);
    assemblePDESingleReduced(mat, rhs, A, B, C, D, X, Y);
}

template<class Scalar>
void DefaultAssembler3D<Scalar>::assemblePDEBoundarySingleReduced(
                                        AbstractSystemMatrix* mat, Data& rhs,
                                        const DataMap& coefs) const
{
    const Data& d = unpackData("d", coefs);
    const Data& y = unpackData("y", coefs);
    assemblePDEBoundarySingleReduced(mat, rhs, d, y);
}

template<class Scalar>
void DefaultAssembler3D<Scalar>::assemblePDESystem(AbstractSystemMatrix* mat,
                                        Data& rhs, const DataMap& coefs) const
{
    const Data& A = unpackData("A", coefs);
    const Data& B = unpackData("B", coefs);
    const Data& C = unpackData("C", coefs);
    const Data& D = unpackData("D", coefs);
    const Data& X = unpackData("X", coefs);
    const Data& Y = unpackData("Y", coefs);
    assemblePDESystem(mat, rhs, A, B, C, D, X, Y);
}

template<class Scalar>
void DefaultAssembler3D<Scalar>::assemblePDEBoundarySystem(
                                        AbstractSystemMatrix* mat, Data& rhs,
                                        const DataMap& coefs) const
{
    const Data& d = unpackData("d", coefs);
    const Data& y = unpackData("y", coefs);
    assemblePDEBoundarySystem(mat, rhs, d, y);
}

template<class Scalar>
void DefaultAssembler3D<Scalar>::assemblePDESystemReduced(
                                        AbstractSystemMatrix* mat, Data& rhs,
                                        const DataMap& coefs) const
{
    const Data& A = unpackData("A", coefs);
    const Data& B = unpackData("B", coefs);
    const Data& C = unpackData("C", coefs);
    const Data& D = unpackData("D", coefs);
    const Data& X = unpackData("X", coefs);
    const Data& Y = unpackData("Y", coefs);
    assemblePDESystemReduced(mat, rhs, A, B, C, D, X, Y);
}

template<class Scalar>
void DefaultAssembler3D<Scalar>::assemblePDEBoundarySystemReduced(
                                        AbstractSystemMatrix* mat,
                                        Data& rhs, const DataMap& coefs) const
{
    const Data& d = unpackData("d", coefs);
    const Data& y = unpackData("y", coefs);
    assemblePDEBoundarySystemReduced(mat, rhs, d, y);
}

#ifdef ESYS_HAVE_TRILINOS
template<class Scalar>
void DefaultAssembler3D<Scalar>::assemblePDEHanging(Teuchos::RCP<Tpetra::CrsMatrix<double,esys_trilinos::LO,esys_trilinos::GO,esys_trilinos::NT>>* mat) const
{

}
#endif

/****************************************************************************/
// PDE SINGLE
/****************************************************************************/

template<class Scalar>
void DefaultAssembler3D<Scalar>::assemblePDESingle(AbstractSystemMatrix* mat,
                                       Data& rhs, const Data& A, const Data& B,
                                       const Data& C, const Data& D,
                                       const Data& X, const Data& Y) const
{
    // A single PDE is numEq=numComp=1; for that block size the escript coefficient
    // layouts collapse exactly to the single-PDE layouts, so reuse the generic
    // (Gauss-order-correct) system assembler rather than the old codegen kernel.
    assemblePDESystem(mat, rhs, A, B, C, D, X, Y);
}

/****************************************************************************/
// PDE SINGLE BOUNDARY
/****************************************************************************/

template<class Scalar>
void DefaultAssembler3D<Scalar>::assemblePDEBoundarySingle(
                                        AbstractSystemMatrix* mat, Data& rhs,
                                        const Data& d, const Data& y) const
{
    // List-based boundary (face) assembly. Each boundary face element is a quad
    // with 4 nodes (borderNodeInfo.neighbours[0..3], ordered (a0b0,a1b0,a0b1,a1b1)
    // in the face's two in-plane axes). d and y live on FunctionOnBoundary and
    // are indexed by the face sample id m_faceOffset[face]+k -- NOT the volume
    // element index. Uses 2x2 Gauss with the standard bilinear shape functions.
    const bool addEM_S = !d.isEmpty();
    const bool addEM_F = !y.isEmpty();
    const Scalar zero = static_cast<Scalar>(0);
    rhs.requireWrite();

    // 1D shape value: node p at Gauss point g (both in {0,1}) is c_hi if p==g
    // else c_lo, with c_hi,c_lo = (1 +/- 1/sqrt(3))/2. Gauss point g=(gi,gj),
    // node a=(pa,qa) matching the interpolateNodesOnFaces ordering.
    const double SQ = 1.0/std::sqrt(3.0);
    const double chi = (1.0 + SQ)/2.0, clo = (1.0 - SQ)/2.0;
    double N[4][4];
    for(int a=0;a<4;++a){ int pa=a&1, qa=(a>>1)&1;
        for(int g=0;g<4;++g){ int gi=g&1, gj=(g>>1)&1;
            N[a][g] = ((pa==gi)?chi:clo) * ((qa==gj)?chi:clo); } }

    const std::vector<borderNodeInfo>* lists[6] = {
        &domain->NodeIDsLeft, &domain->NodeIDsRight, &domain->NodeIDsBottom,
        &domain->NodeIDsTop, &domain->NodeIDsAbove, &domain->NodeIDsBelow };
    static const int planeAxes[6][2] = {{1,2},{1,2},{0,2},{0,2},{0,1},{0,1}};

    std::vector<Scalar> EM_S(4*4, zero), EM_F(4, zero);
    for(int fc=0; fc<6; ++fc) {
        if(domain->m_faceOffset[fc] < 0) continue;
        const std::vector<borderNodeInfo>& L = *lists[fc];
        for(size_t k=0; k<L.size(); ++k) {
            if(addEM_S) std::fill(EM_S.begin(), EM_S.end(), zero);
            if(addEM_F) std::fill(EM_F.begin(), EM_F.end(), zero);
            const double h = (double)(1 << L[k].level);
            const double A = (domain->m_NX[planeAxes[fc][0]]/h)
                           * (domain->m_NX[planeAxes[fc][1]]/h);
            const index_t sample = domain->m_faceOffset[fc] + (index_t)k;

            if(addEM_S) {
                const Scalar* d_p = d.getSampleDataRO(sample, zero);
                if(d.actsExpanded()) {
                    for(int a=0;a<4;++a) for(int b=0;b<4;++b) {
                        Scalar s = zero;
                        for(int g=0;g<4;++g) s += d_p[g]*N[a][g]*N[b][g];
                        EM_S[INDEX2(a,b,4)] = s*(A/4.);
                    }
                } else {
                    for(int a=0;a<4;++a) for(int b=0;b<4;++b) {
                        double s = 0.;
                        for(int g=0;g<4;++g) s += N[a][g]*N[b][g];
                        EM_S[INDEX2(a,b,4)] = d_p[0]*(A/4.)*s;
                    }
                }
            }
            if(addEM_F) {
                const Scalar* y_p = y.getSampleDataRO(sample, zero);
                if(y.actsExpanded()) {
                    for(int a=0;a<4;++a) {
                        Scalar s = zero;
                        for(int g=0;g<4;++g) s += N[a][g]*y_p[g];
                        EM_F[a] = s*(A/4.);
                    }
                } else {
                    for(int a=0;a<4;++a) EM_F[a] = y_p[0]*(A/4.);  // sum_g N[a][g]=1
                }
            }
            domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, addEM_S, addEM_F, L[k]);
        }
    }
}

/****************************************************************************/
// PDE SINGLE REDUCED
/****************************************************************************/

template<class Scalar>
void DefaultAssembler3D<Scalar>::assemblePDESingleReduced(
                                       AbstractSystemMatrix* mat, Data& rhs,
                                       const Data& A, const Data& B,
                                       const Data& C, const Data& D,
                                       const Data& X, const Data& Y) const
{
    assemblePDESystemReduced(mat, rhs, A, B, C, D, X, Y);
}

/****************************************************************************/
// PDE SINGLE REDUCED BOUNDARY
/****************************************************************************/

template<class Scalar>
void DefaultAssembler3D<Scalar>::assemblePDEBoundarySingleReduced(
                                        AbstractSystemMatrix* mat, Data& rhs,
                                        const Data& d, const Data& y) const
{
    // Reduced (1-point) boundary assembly: the single centre point has weight A
    // and every shape function equals 1/4 there.
    const bool addEM_S = !d.isEmpty();
    const bool addEM_F = !y.isEmpty();
    const Scalar zero = static_cast<Scalar>(0);
    rhs.requireWrite();

    const std::vector<borderNodeInfo>* lists[6] = {
        &domain->NodeIDsLeft, &domain->NodeIDsRight, &domain->NodeIDsBottom,
        &domain->NodeIDsTop, &domain->NodeIDsAbove, &domain->NodeIDsBelow };
    static const int planeAxes[6][2] = {{1,2},{1,2},{0,2},{0,2},{0,1},{0,1}};

    std::vector<Scalar> EM_S(4*4, zero), EM_F(4, zero);
    for(int fc=0; fc<6; ++fc) {
        if(domain->m_faceOffset[fc] < 0) continue;
        const std::vector<borderNodeInfo>& L = *lists[fc];
        for(size_t k=0; k<L.size(); ++k) {
            if(addEM_S) std::fill(EM_S.begin(), EM_S.end(), zero);
            if(addEM_F) std::fill(EM_F.begin(), EM_F.end(), zero);
            const double h = (double)(1 << L[k].level);
            const double A = (domain->m_NX[planeAxes[fc][0]]/h)
                           * (domain->m_NX[planeAxes[fc][1]]/h);
            const index_t sample = domain->m_faceOffset[fc] + (index_t)k;
            if(addEM_S) {
                const Scalar* d_p = d.getSampleDataRO(sample, zero);
                for(int a=0;a<4;++a) for(int b=0;b<4;++b)
                    EM_S[INDEX2(a,b,4)] = d_p[0]*(A/16.);
            }
            if(addEM_F) {
                const Scalar* y_p = y.getSampleDataRO(sample, zero);
                for(int a=0;a<4;++a) EM_F[a] = y_p[0]*(A/4.);
            }
            domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, addEM_S, addEM_F, L[k]);
        }
    }
}

/****************************************************************************/
// PDE SYSTEM
/****************************************************************************/

template<class Scalar>
void DefaultAssembler3D<Scalar>::assemblePDESystem(AbstractSystemMatrix* mat,
                                        Data& rhs, const Data& A, const Data& B,
                                        const Data& C, const Data& D,
                                        const Data& X, const Data& Y) const
{
    // Generic vector/system interior assembly by 2x2x2 Gauss quadrature on the
    // trilinear hex (axis-aligned, so the Jacobian is diag(hx,hy,hz)).
    const int DIM = 3;
    const bool addEM_S = (!A.isEmpty() || !B.isEmpty() || !C.isEmpty() || !D.isEmpty());
    const bool addEM_F = (!X.isEmpty() || !Y.isEmpty());
    const Scalar zero = static_cast<Scalar>(0);
    rhs.requireWrite();

    dim_t numEq, numComp;
    if (!mat) {
        numEq = numComp = (rhs.isEmpty() ? 1 : rhs.getDataPointSize());
    } else {
        numEq = mat->getRowBlockSize();
        numComp = mat->getColumnBlockSize();
    }

    // reference shape values / gradients at the 8 Gauss points.
    const double SQ = 1.0/std::sqrt(3.0);
    const double chi = (1.0 + SQ)/2.0, clo = (1.0 - SQ)/2.0;
    // N1d(p,gi) = chi if p==gi else clo ; sign(p) = p? +1 : -1
    double N[8][8];            // N[a][g]
    double gradRef[8][8][3];   // dN_a/dxi at g, reference (no 1/h)
    for(int a=0;a<8;++a){
        const int pa=a&1, qa=(a>>1)&1, ra=(a>>2)&1;
        for(int g=0;g<8;++g){
            const int gi=(g>>2)&1, gj=(g>>1)&1, gk=g&1;  // interp gauss order (x slowest)
            const double Xa=(pa==gi)?chi:clo, Ya=(qa==gj)?chi:clo, Za=(ra==gk)?chi:clo;
            N[a][g] = Xa*Ya*Za;
            const double sx=pa?1.0:-1.0, sy=qa?1.0:-1.0, sz=ra?1.0:-1.0;
            gradRef[a][g][0] = sx*Ya*Za;
            gradRef[a][g][1] = Xa*sy*Za;
            gradRef[a][g][2] = Xa*Ya*sz;
        }
    }

    const bool Aex=A.actsExpanded(), Bex=B.actsExpanded(), Cex=C.actsExpanded();
    const bool Dex=D.actsExpanded(), Xex=X.actsExpanded(), Yex=Y.actsExpanded();

    std::vector<Scalar> EM_S(8*8*numEq*numComp, zero), EM_F(8*numEq, zero);

    long idbase = 0; // running local leaf index
    for (p8est_topidx_t t = domain->p8est->first_local_tree; t <= domain->p8est->last_local_tree; t++)
    {
        p8est_tree_t * currenttree = p8est_tree_array_index(domain->p8est->trees, t);
        sc_array_t * tquadrants = &currenttree->quadrants;
        p8est_locidx_t Q = (p8est_locidx_t) tquadrants->elem_count;
        for (int q = 0; q < Q; ++q, ++idbase)
        {
            p8est_quadrant_t * quad = p8est_quadrant_array_index(tquadrants, q);
            const int l = quad->level;
            const double hh = (double)(1 << l);
            const double h[3] = { domain->m_NX[0]/hh, domain->m_NX[1]/hh, domain->m_NX[2]/hh };
            const double V = h[0]*h[1]*h[2];
            const double w = V/8.0;                 // Gauss weight
            const long id = idbase;                 // element coefficient sample

            if(addEM_S) std::fill(EM_S.begin(), EM_S.end(), zero);
            if(addEM_F) std::fill(EM_F.begin(), EM_F.end(), zero);

            const Scalar* A_p = A.isEmpty()? nullptr : A.getSampleDataRO(id, zero);
            const Scalar* B_p = B.isEmpty()? nullptr : B.getSampleDataRO(id, zero);
            const Scalar* C_p = C.isEmpty()? nullptr : C.getSampleDataRO(id, zero);
            const Scalar* D_p = D.isEmpty()? nullptr : D.getSampleDataRO(id, zero);
            const Scalar* X_p = X.isEmpty()? nullptr : X.getSampleDataRO(id, zero);
            const Scalar* Y_p = Y.isEmpty()? nullptr : Y.getSampleDataRO(id, zero);

            // physical gradients per (a,g)
            for(int g=0; g<8; ++g) {
                double gr[8][3];
                for(int a=0;a<8;++a) for(int dd=0;dd<3;++dd) gr[a][dd]=gradRef[a][g][dd]/h[dd];

                if(addEM_S) {
                    for(index_t k=0;k<numEq;++k) for(index_t m=0;m<numComp;++m) {
                        for(int a=0;a<8;++a) for(int b=0;b<8;++b) {
                            Scalar v = zero;
                            if(A_p) {
                                for(int i=0;i<DIM;++i) for(int j=0;j<DIM;++j) {
                                    const Scalar Aij = Aex ? A_p[INDEX5(k,i,m,j,g,numEq,DIM,numComp,DIM)]
                                                           : A_p[INDEX4(k,i,m,j,numEq,DIM,numComp)];
                                    v += Aij * gr[a][i]*gr[b][j];
                                }
                            }
                            if(B_p) {
                                for(int i=0;i<DIM;++i) {
                                    const Scalar Bim = Bex ? B_p[INDEX4(k,i,m,g,numEq,DIM,numComp)]
                                                           : B_p[INDEX3(k,i,m,numEq,DIM)];
                                    v += Bim * gr[a][i]*N[b][g];
                                }
                            }
                            if(C_p) {
                                for(int j=0;j<DIM;++j) {
                                    const Scalar Cmj = Cex ? C_p[INDEX4(k,m,j,g,numEq,numComp,DIM)]
                                                           : C_p[INDEX3(k,m,j,numEq,numComp)];
                                    v += Cmj * N[a][g]*gr[b][j];
                                }
                            }
                            if(D_p) {
                                const Scalar Dkm = Dex ? D_p[INDEX3(k,m,g,numEq,numComp)]
                                                       : D_p[INDEX2(k,m,numEq)];
                                v += Dkm * N[a][g]*N[b][g];
                            }
                            EM_S[INDEX4(k,m,a,b,numEq,numComp,8)] += w*v;
                        }
                    }
                }
                if(addEM_F) {
                    for(index_t k=0;k<numEq;++k) for(int a=0;a<8;++a) {
                        Scalar v = zero;
                        if(X_p) {
                            for(int i=0;i<DIM;++i) {
                                const Scalar Xki = Xex ? X_p[INDEX3(k,i,g,numEq,DIM)]
                                                       : X_p[INDEX2(k,i,numEq)];
                                v += Xki * gr[a][i];
                            }
                        }
                        if(Y_p) {
                            const Scalar Yk = Yex ? Y_p[INDEX2(k,g,numEq)] : Y_p[k];
                            v += Yk * N[a][g];
                        }
                        EM_F[INDEX2(k,a,numEq)] += w*v;
                    }
                }
            }
            domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, addEM_S, addEM_F, q, t, numEq, numComp);
        }
    }
}

/****************************************************************************/
// PDE SYSTEM BOUNDARY
/****************************************************************************/

template<class Scalar>
void DefaultAssembler3D<Scalar>::assemblePDEBoundarySystem(
                                        AbstractSystemMatrix* mat, Data& rhs,
                                        const Data& d, const Data& y) const
{
    // Vector/system boundary assembly, list-based (see assemblePDEBoundarySingle).
    // Block layout: EM_S[INDEX4(k,m,a,b,numEq,numComp,4)], EM_F[INDEX2(k,a,numEq)].
    const bool addEM_S = !d.isEmpty();
    const bool addEM_F = !y.isEmpty();
    const Scalar zero = static_cast<Scalar>(0);
    rhs.requireWrite();

    dim_t numEq, numComp;
    if (!mat) {
        numEq = numComp = (rhs.isEmpty() ? 1 : rhs.getDataPointSize());
    } else {
        numEq = mat->getRowBlockSize();
        numComp = mat->getColumnBlockSize();
    }

    const double SQ = 1.0/std::sqrt(3.0);
    const double chi = (1.0 + SQ)/2.0, clo = (1.0 - SQ)/2.0;
    double N[4][4];
    for(int a=0;a<4;++a){ int pa=a&1, qa=(a>>1)&1;
        for(int g=0;g<4;++g){ int gi=g&1, gj=(g>>1)&1;
            N[a][g] = ((pa==gi)?chi:clo) * ((qa==gj)?chi:clo); } }

    const std::vector<borderNodeInfo>* lists[6] = {
        &domain->NodeIDsLeft, &domain->NodeIDsRight, &domain->NodeIDsBottom,
        &domain->NodeIDsTop, &domain->NodeIDsAbove, &domain->NodeIDsBelow };
    static const int planeAxes[6][2] = {{1,2},{1,2},{0,2},{0,2},{0,1},{0,1}};

    std::vector<Scalar> EM_S(4*4*numEq*numComp, zero), EM_F(4*numEq, zero);
    for(int fc=0; fc<6; ++fc) {
        if(domain->m_faceOffset[fc] < 0) continue;
        const std::vector<borderNodeInfo>& L = *lists[fc];
        for(size_t kk=0; kk<L.size(); ++kk) {
            if(addEM_S) std::fill(EM_S.begin(), EM_S.end(), zero);
            if(addEM_F) std::fill(EM_F.begin(), EM_F.end(), zero);
            const double h = (double)(1 << L[kk].level);
            const double A = (domain->m_NX[planeAxes[fc][0]]/h)
                           * (domain->m_NX[planeAxes[fc][1]]/h);
            const index_t sample = domain->m_faceOffset[fc] + (index_t)kk;

            if(addEM_S) {
                const Scalar* d_p = d.getSampleDataRO(sample, zero);
                const bool ex = d.actsExpanded();
                for(index_t k=0;k<numEq;++k) for(index_t m=0;m<numComp;++m)
                    for(int a=0;a<4;++a) for(int b=0;b<4;++b) {
                        Scalar s = zero;
                        if(ex) {
                            for(int g=0;g<4;++g)
                                s += d_p[INDEX3(k,m,g,numEq,numComp)]*N[a][g]*N[b][g];
                            s *= (A/4.);
                        } else {
                            double gg=0.; for(int g=0;g<4;++g) gg += N[a][g]*N[b][g];
                            s = d_p[INDEX2(k,m,numEq)]*(A/4.)*gg;
                        }
                        EM_S[INDEX4(k,m,a,b,numEq,numComp,4)] = s;
                    }
            }
            if(addEM_F) {
                const Scalar* y_p = y.getSampleDataRO(sample, zero);
                const bool ex = y.actsExpanded();
                for(index_t k=0;k<numEq;++k) for(int a=0;a<4;++a) {
                    Scalar s = zero;
                    if(ex) { for(int g=0;g<4;++g) s += N[a][g]*y_p[INDEX2(k,g,numEq)]; s *= (A/4.); }
                    else s = y_p[k]*(A/4.);
                    EM_F[INDEX2(k,a,numEq)] = s;
                }
            }
            domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, addEM_S, addEM_F, L[kk], numEq, numComp);
        }
    }
}

/****************************************************************************/
// PDE SYSTEM REDUCED
/****************************************************************************/

template<class Scalar>
void DefaultAssembler3D<Scalar>::assemblePDESystemReduced(
                                        AbstractSystemMatrix* mat, Data& rhs,
                                        const Data& A, const Data& B,
                                        const Data& C, const Data& D,
                                        const Data& X, const Data& Y) const
{
    // Reduced (1-point centre) vector/system interior assembly on the trilinear
    // hex. At the centre N_a = 1/8 and dN_a/dx_d = sign(pa_d)/(4 h_d); weight V.
    const int DIM = 3;
    const bool addEM_S = (!A.isEmpty() || !B.isEmpty() || !C.isEmpty() || !D.isEmpty());
    const bool addEM_F = (!X.isEmpty() || !Y.isEmpty());
    const Scalar zero = static_cast<Scalar>(0);
    rhs.requireWrite();

    dim_t numEq, numComp;
    if (!mat) {
        numEq = numComp = (rhs.isEmpty() ? 1 : rhs.getDataPointSize());
    } else {
        numEq = mat->getRowBlockSize();
        numComp = mat->getColumnBlockSize();
    }

    double N[8];
    double gsign[8][3];
    for(int a=0;a<8;++a){
        const int pa=a&1, qa=(a>>1)&1, ra=(a>>2)&1;
        N[a] = 1.0/8.0;
        gsign[a][0] = (pa?1.0:-1.0)*0.25;
        gsign[a][1] = (qa?1.0:-1.0)*0.25;
        gsign[a][2] = (ra?1.0:-1.0)*0.25;
    }

    std::vector<Scalar> EM_S(8*8*numEq*numComp, zero), EM_F(8*numEq, zero);

    long idbase = 0;
    for (p8est_topidx_t t = domain->p8est->first_local_tree; t <= domain->p8est->last_local_tree; t++)
    {
        p8est_tree_t * currenttree = p8est_tree_array_index(domain->p8est->trees, t);
        sc_array_t * tquadrants = &currenttree->quadrants;
        p8est_locidx_t Q = (p8est_locidx_t) tquadrants->elem_count;
        for (int q = 0; q < Q; ++q, ++idbase)
        {
            p8est_quadrant_t * quad = p8est_quadrant_array_index(tquadrants, q);
            const double hh = (double)(1 << quad->level);
            const double h[3] = { domain->m_NX[0]/hh, domain->m_NX[1]/hh, domain->m_NX[2]/hh };
            const double V = h[0]*h[1]*h[2];
            const long id = idbase;

            if(addEM_S) std::fill(EM_S.begin(), EM_S.end(), zero);
            if(addEM_F) std::fill(EM_F.begin(), EM_F.end(), zero);

            const Scalar* A_p = A.isEmpty()? nullptr : A.getSampleDataRO(id, zero);
            const Scalar* B_p = B.isEmpty()? nullptr : B.getSampleDataRO(id, zero);
            const Scalar* C_p = C.isEmpty()? nullptr : C.getSampleDataRO(id, zero);
            const Scalar* D_p = D.isEmpty()? nullptr : D.getSampleDataRO(id, zero);
            const Scalar* X_p = X.isEmpty()? nullptr : X.getSampleDataRO(id, zero);
            const Scalar* Y_p = Y.isEmpty()? nullptr : Y.getSampleDataRO(id, zero);

            double gr[8][3];
            for(int a=0;a<8;++a) for(int dd=0;dd<3;++dd) gr[a][dd]=gsign[a][dd]/h[dd];

            if(addEM_S) {
                for(index_t k=0;k<numEq;++k) for(index_t m=0;m<numComp;++m)
                    for(int a=0;a<8;++a) for(int b=0;b<8;++b) {
                        Scalar v = zero;
                        if(A_p) for(int i=0;i<DIM;++i) for(int j=0;j<DIM;++j)
                            v += A_p[INDEX4(k,i,m,j,numEq,DIM,numComp)]*gr[a][i]*gr[b][j];
                        if(B_p) for(int i=0;i<DIM;++i)
                            v += B_p[INDEX3(k,i,m,numEq,DIM)]*gr[a][i]*N[b];
                        if(C_p) for(int j=0;j<DIM;++j)
                            v += C_p[INDEX3(k,m,j,numEq,numComp)]*N[a]*gr[b][j];
                        if(D_p) v += D_p[INDEX2(k,m,numEq)]*N[a]*N[b];
                        EM_S[INDEX4(k,m,a,b,numEq,numComp,8)] = V*v;
                    }
            }
            if(addEM_F) {
                for(index_t k=0;k<numEq;++k) for(int a=0;a<8;++a) {
                    Scalar v = zero;
                    if(X_p) for(int i=0;i<DIM;++i) v += X_p[INDEX2(k,i,numEq)]*gr[a][i];
                    if(Y_p) v += Y_p[k]*N[a];
                    EM_F[INDEX2(k,a,numEq)] = V*v;
                }
            }
            domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, addEM_S, addEM_F, q, t, numEq, numComp);
        }
    }
}

/****************************************************************************/
// PDE SYSTEM REDUCED BOUNDARY
/****************************************************************************/

template<class Scalar>
void DefaultAssembler3D<Scalar>::assemblePDEBoundarySystemReduced(
                                        AbstractSystemMatrix* mat, Data& rhs,
                                        const Data& d, const Data& y) const
{
    // Reduced vector/system boundary assembly (1-point rule, N=1/4, weight A).
    const bool addEM_S = !d.isEmpty();
    const bool addEM_F = !y.isEmpty();
    const Scalar zero = static_cast<Scalar>(0);
    rhs.requireWrite();

    dim_t numEq, numComp;
    if (!mat) {
        numEq = numComp = (rhs.isEmpty() ? 1 : rhs.getDataPointSize());
    } else {
        numEq = mat->getRowBlockSize();
        numComp = mat->getColumnBlockSize();
    }

    const std::vector<borderNodeInfo>* lists[6] = {
        &domain->NodeIDsLeft, &domain->NodeIDsRight, &domain->NodeIDsBottom,
        &domain->NodeIDsTop, &domain->NodeIDsAbove, &domain->NodeIDsBelow };
    static const int planeAxes[6][2] = {{1,2},{1,2},{0,2},{0,2},{0,1},{0,1}};

    std::vector<Scalar> EM_S(4*4*numEq*numComp, zero), EM_F(4*numEq, zero);
    for(int fc=0; fc<6; ++fc) {
        if(domain->m_faceOffset[fc] < 0) continue;
        const std::vector<borderNodeInfo>& L = *lists[fc];
        for(size_t kk=0; kk<L.size(); ++kk) {
            if(addEM_S) std::fill(EM_S.begin(), EM_S.end(), zero);
            if(addEM_F) std::fill(EM_F.begin(), EM_F.end(), zero);
            const double h = (double)(1 << L[kk].level);
            const double A = (domain->m_NX[planeAxes[fc][0]]/h)
                           * (domain->m_NX[planeAxes[fc][1]]/h);
            const index_t sample = domain->m_faceOffset[fc] + (index_t)kk;
            if(addEM_S) {
                const Scalar* d_p = d.getSampleDataRO(sample, zero);
                for(index_t k=0;k<numEq;++k) for(index_t m=0;m<numComp;++m)
                    for(int a=0;a<4;++a) for(int b=0;b<4;++b)
                        EM_S[INDEX4(k,m,a,b,numEq,numComp,4)] = d_p[INDEX2(k,m,numEq)]*(A/16.);
            }
            if(addEM_F) {
                const Scalar* y_p = y.getSampleDataRO(sample, zero);
                for(index_t k=0;k<numEq;++k) for(int a=0;a<4;++a)
                    EM_F[INDEX2(k,a,numEq)] = y_p[k]*(A/4.);
            }
            domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, addEM_S, addEM_F, L[kk], numEq, numComp);
        }
    }
}

// instantiate our two supported versions
template class DefaultAssembler3D<escript::DataTypes::real_t>;
template class DefaultAssembler3D<escript::DataTypes::cplx_t>;

} // namespace oxley

