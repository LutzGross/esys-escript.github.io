
/*****************************************************************************
*
* Copyright (c) 2003-2020 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#include <ripley/LameAssembler3D.h>
#include <ripley/domainhelpers.h>

#include <escript/index.h>

using namespace std;

using escript::AbstractSystemMatrix;
using escript::Data;

namespace ripley {

void LameAssembler3D::collateFunctionSpaceTypes(vector<int>& fsTypes,
                                                const DataMap& coefs) const
{
    if (isNotEmpty("lame_mu", coefs))
        fsTypes.push_back(coefs.find("lame_mu")->second.getFunctionSpace().getTypeCode());
    if (isNotEmpty("lame_lambda", coefs))
        fsTypes.push_back(coefs.find("lame_lambda")->second.getFunctionSpace().getTypeCode());
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


void LameAssembler3D::assemblePDESingle(AbstractSystemMatrix* mat, Data& rhs,
                                        const DataMap& coefs) const
{
    throw RipleyException("assemblePDESingle not implemented in LameAssembler3D");
}

void LameAssembler3D::assemblePDEBoundarySingle(AbstractSystemMatrix* mat,
                                        Data& rhs, const DataMap& coefs) const 
{
    throw RipleyException("assemblePDESingle not implemented in LameAssembler3D");
}

void LameAssembler3D::assemblePDESingleReduced(AbstractSystemMatrix* mat,
                                        Data& rhs, const DataMap& coefs) const
{
    throw RipleyException("assemblePDESingle not implemented in LameAssembler3D");
}

void LameAssembler3D::assemblePDEBoundarySingleReduced(
                                        AbstractSystemMatrix* mat,
                                        Data& rhs, const DataMap& coefs) const
{
    throw RipleyException("assemblePDESingle not implemented in LameAssembler3D");
}

void LameAssembler3D::assemblePDESystemReduced(AbstractSystemMatrix* mat,
                                        Data& rhs, const DataMap& coefs) const
{
    throw RipleyException("assemblePDESystemReduced not implemented in LameAssembler3D");
}

void LameAssembler3D::assemblePDEBoundarySystemReduced(
                                        AbstractSystemMatrix* mat,
                                        Data& rhs, const DataMap& coefs) const
{
    throw RipleyException("assemblePDEBoundarySystemReduced not implemented in LameAssembler3D");
}

void LameAssembler3D::assemblePDEBoundarySystem(AbstractSystemMatrix* mat,
                                        Data& rhs, const DataMap& coefs) const
{
    const Data& d = unpackData("d", coefs);
    const Data& y = unpackData("y", coefs);
    dim_t numEq, numComp;
    if (!mat)
        numEq=numComp=(rhs.isEmpty() ? 1 : rhs.getDataPointSize());
    else {
        numEq=mat->getRowBlockSize();
        numComp=mat->getColumnBlockSize();
    }
    const double SQRT3 = 1.73205080756887719318;
    const double w12 = m_dx[0]*m_dx[1]/144;
    const double w10 = w12*(-SQRT3 + 2);
    const double w11 = w12*(SQRT3 + 2);
    const double w13 = w12*(-4*SQRT3 + 7);
    const double w14 = w12*(4*SQRT3 + 7);
    const double w7 = m_dx[0]*m_dx[2]/144;
    const double w5 = w7*(-SQRT3 + 2);
    const double w6 = w7*(SQRT3 + 2);
    const double w8 = w7*(-4*SQRT3 + 7);
    const double w9 = w7*(4*SQRT3 + 7);
    const double w2 = m_dx[1]*m_dx[2]/144;
    const double w0 = w2*(-SQRT3 + 2);
    const double w1 = w2*(SQRT3 + 2);
    const double w3 = w2*(-4*SQRT3 + 7);
    const double w4 = w2*(4*SQRT3 + 7);
    const bool add_EM_S = !d.isEmpty();
    const bool add_EM_F = !y.isEmpty();
    rhs.requireWrite();

#pragma omp parallel
    {
        vector<double> EM_S(8*8*numEq*numComp, 0);
        vector<double> EM_F(8*numEq, 0);

        if (domain->m_faceOffset[0] > -1) {
            if (add_EM_S)
                fill(EM_S.begin(), EM_S.end(), 0);
            if (add_EM_F)
                fill(EM_F.begin(), EM_F.end(), 0);

            for (index_t k2_0=0; k2_0<2; k2_0++) { // colouring
#pragma omp for
                for (index_t k2=k2_0; k2<m_NE[2]; k2+=2) {
                    for (index_t k1=0; k1<m_NE[1]; ++k1) {
                        const index_t e = INDEX2(k1,k2,m_NE[1]);
                        ///////////////
                        // process d //
                        ///////////////
                        if (add_EM_S) {
                            const double* d_p=d.getSampleDataRO(e);
                            if (d.actsExpanded()) {
                                for (index_t k=0; k<numEq; k++) {
                                    for (index_t m=0; m<numComp; m++) {
                                        const double d_0 = d_p[INDEX3(k,m,0,numEq,numComp)];
                                        const double d_1 = d_p[INDEX3(k,m,1,numEq,numComp)];
                                        const double d_2 = d_p[INDEX3(k,m,2,numEq,numComp)];
                                        const double d_3 = d_p[INDEX3(k,m,3,numEq,numComp)];
                                        const double tmp0 = w0*(d_0 + d_1);
                                        const double tmp1 = w1*(d_2 + d_3);
                                        const double tmp2 = w0*(d_0 + d_2);
                                        const double tmp3 = w1*(d_1 + d_3);
                                        const double tmp4 = w0*(d_1 + d_3);
                                        const double tmp5 = w1*(d_0 + d_2);
                                        const double tmp6 = w0*(d_2 + d_3);
                                        const double tmp7 = w1*(d_0 + d_1);
                                        const double tmp8 = w2*(d_0 + d_3);
                                        const double tmp9 = w2*(d_1 + d_2);
                                        const double tmp10 = w2*(d_0 + d_1 + d_2 + d_3);
                                        EM_S[INDEX4(k,m,0,0,numEq,numComp,8)] = d_0*w4 + d_3*w3 + tmp9;
                                        EM_S[INDEX4(k,m,0,2,numEq,numComp,8)] = tmp6 + tmp7;
                                        EM_S[INDEX4(k,m,0,4,numEq,numComp,8)] = tmp4 + tmp5;
                                        EM_S[INDEX4(k,m,0,6,numEq,numComp,8)] = tmp10;
                                        EM_S[INDEX4(k,m,2,0,numEq,numComp,8)] = tmp6 + tmp7;
                                        EM_S[INDEX4(k,m,2,2,numEq,numComp,8)] = d_1*w4 + d_2*w3 + tmp8;
                                        EM_S[INDEX4(k,m,2,4,numEq,numComp,8)] = tmp10;
                                        EM_S[INDEX4(k,m,2,6,numEq,numComp,8)] = tmp2 + tmp3;
                                        EM_S[INDEX4(k,m,4,0,numEq,numComp,8)] = tmp4 + tmp5;
                                        EM_S[INDEX4(k,m,4,2,numEq,numComp,8)] = tmp10;
                                        EM_S[INDEX4(k,m,4,4,numEq,numComp,8)] = d_1*w3 + d_2*w4 + tmp8;
                                        EM_S[INDEX4(k,m,4,6,numEq,numComp,8)] = tmp0 + tmp1;
                                        EM_S[INDEX4(k,m,6,0,numEq,numComp,8)] = tmp10;
                                        EM_S[INDEX4(k,m,6,2,numEq,numComp,8)] = tmp2 + tmp3;
                                        EM_S[INDEX4(k,m,6,4,numEq,numComp,8)] = tmp0 + tmp1;
                                        EM_S[INDEX4(k,m,6,6,numEq,numComp,8)] = d_0*w3 + d_3*w4 + tmp9;
                                    }
                                 }
                            } else { // constant data
                                for (index_t k=0; k<numEq; k++) {
                                    for (index_t m=0; m<numComp; m++) {
                                        const double wd0 = 4*d_p[INDEX2(k, m, numEq)]*w2;
                                        EM_S[INDEX4(k,m,0,0,numEq,numComp,8)] = 4*wd0;
                                        EM_S[INDEX4(k,m,0,2,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,0,4,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,0,6,numEq,numComp,8)] =   wd0;
                                        EM_S[INDEX4(k,m,2,0,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,2,2,numEq,numComp,8)] = 4*wd0;
                                        EM_S[INDEX4(k,m,2,4,numEq,numComp,8)] =   wd0;
                                        EM_S[INDEX4(k,m,2,6,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,4,0,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,4,2,numEq,numComp,8)] =   wd0;
                                        EM_S[INDEX4(k,m,4,4,numEq,numComp,8)] = 4*wd0;
                                        EM_S[INDEX4(k,m,4,6,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,6,0,numEq,numComp,8)] =   wd0;
                                        EM_S[INDEX4(k,m,6,2,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,6,4,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,6,6,numEq,numComp,8)] = 4*wd0;
                                    }
                                }
                            }
                        }
                        ///////////////
                        // process y //
                        ///////////////
                        if (add_EM_F) {
                            const double* y_p=y.getSampleDataRO(e);
                            if (y.actsExpanded()) {
                                for (index_t k=0; k<numEq; k++) {
                                    const double y_0 = y_p[INDEX2(k, 0, numEq)];
                                    const double y_1 = y_p[INDEX2(k, 1, numEq)];
                                    const double y_2 = y_p[INDEX2(k, 2, numEq)];
                                    const double y_3 = y_p[INDEX2(k, 3, numEq)];
                                    const double tmp0 = 6*w2*(y_1 + y_2);
                                    const double tmp1 = 6*w2*(y_0 + y_3);
                                    EM_F[INDEX2(k,0,numEq)] = tmp0 + 6*w0*y_3 + 6*w1*y_0;
                                    EM_F[INDEX2(k,2,numEq)] = tmp1 + 6*w0*y_2 + 6*w1*y_1;
                                    EM_F[INDEX2(k,4,numEq)] = tmp1 + 6*w0*y_1 + 6*w1*y_2;
                                    EM_F[INDEX2(k,6,numEq)] = tmp0 + 6*w0*y_0 + 6*w1*y_3;
                                }
                            } else { // constant data
                                for (index_t k=0; k<numEq; k++) {
                                    EM_F[INDEX2(k,0,numEq)] = 36*w2*y_p[k];
                                    EM_F[INDEX2(k,2,numEq)] = 36*w2*y_p[k];
                                    EM_F[INDEX2(k,4,numEq)] = 36*w2*y_p[k];
                                    EM_F[INDEX2(k,6,numEq)] = 36*w2*y_p[k];
                                }
                            }
                        }
                        const index_t firstNode=m_NN[0]*m_NN[1]*k2+m_NN[0]*k1;
                        domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F,
                                add_EM_S, add_EM_F, firstNode, numEq, numComp);
                    } // k1 loop
                } // k2 loop
            } // colouring
        } // face 0

        if (domain->m_faceOffset[1] > -1) {
            if (add_EM_S)
                fill(EM_S.begin(), EM_S.end(), 0);
            if (add_EM_F)
                fill(EM_F.begin(), EM_F.end(), 0);

            for (index_t k2_0=0; k2_0<2; k2_0++) { // colouring
#pragma omp for
                for (index_t k2=k2_0; k2<m_NE[2]; k2+=2) {
                    for (index_t k1=0; k1<m_NE[1]; ++k1) {
                        const index_t e = domain->m_faceOffset[1]+INDEX2(k1,k2,m_NE[1]);
                        ///////////////
                        // process d //
                        ///////////////
                        if (add_EM_S) {
                            const double* d_p=d.getSampleDataRO(e);
                            if (d.actsExpanded()) {
                                for (index_t k=0; k<numEq; k++) {
                                    for (index_t m=0; m<numComp; m++) {
                                        const double d_0 = d_p[INDEX3(k,m,0,numEq,numComp)];
                                        const double d_1 = d_p[INDEX3(k,m,1,numEq,numComp)];
                                        const double d_2 = d_p[INDEX3(k,m,2,numEq,numComp)];
                                        const double d_3 = d_p[INDEX3(k,m,3,numEq,numComp)];
                                        const double tmp0 = w0*(d_0 + d_2);
                                        const double tmp1 = w1*(d_1 + d_3);
                                        const double tmp2 = w0*(d_2 + d_3);
                                        const double tmp3 = w1*(d_0 + d_1);
                                        const double tmp4 = w0*(d_1 + d_3);
                                        const double tmp5 = w1*(d_0 + d_2);
                                        const double tmp6 = w2*(d_0 + d_3);
                                        const double tmp7 = w2*(d_1 + d_2);
                                        const double tmp8 = w0*(d_0 + d_1);
                                        const double tmp9 = w1*(d_2 + d_3);
                                        const double tmp10 = w2*(d_0 + d_1 + d_2 + d_3);
                                        EM_S[INDEX4(k,m,1,1,numEq,numComp,8)] = d_0*w4 + d_3*w3 + tmp7;
                                        EM_S[INDEX4(k,m,1,3,numEq,numComp,8)] = tmp2 + tmp3;
                                        EM_S[INDEX4(k,m,1,5,numEq,numComp,8)] = tmp4 + tmp5;
                                        EM_S[INDEX4(k,m,1,7,numEq,numComp,8)] = tmp10;
                                        EM_S[INDEX4(k,m,3,1,numEq,numComp,8)] = tmp2 + tmp3;
                                        EM_S[INDEX4(k,m,3,3,numEq,numComp,8)] = d_1*w4 + d_2*w3 + tmp6;
                                        EM_S[INDEX4(k,m,3,5,numEq,numComp,8)] = tmp10;
                                        EM_S[INDEX4(k,m,3,7,numEq,numComp,8)] = tmp0 + tmp1;
                                        EM_S[INDEX4(k,m,5,1,numEq,numComp,8)] = tmp4 + tmp5;
                                        EM_S[INDEX4(k,m,5,3,numEq,numComp,8)] = tmp10;
                                        EM_S[INDEX4(k,m,5,5,numEq,numComp,8)] = d_1*w3 + d_2*w4 + tmp6;
                                        EM_S[INDEX4(k,m,5,7,numEq,numComp,8)] = tmp8 + tmp9;
                                        EM_S[INDEX4(k,m,7,1,numEq,numComp,8)] = tmp10;
                                        EM_S[INDEX4(k,m,7,3,numEq,numComp,8)] = tmp0 + tmp1;
                                        EM_S[INDEX4(k,m,7,5,numEq,numComp,8)] = tmp8 + tmp9;
                                        EM_S[INDEX4(k,m,7,7,numEq,numComp,8)] = d_0*w3 + d_3*w4 + tmp7;
                                    }
                                 }
                            } else { // constant data
                                for (index_t k=0; k<numEq; k++) {
                                    for (index_t m=0; m<numComp; m++) {
                                        const double wd0 = 4*d_p[INDEX2(k, m, numEq)]*w2;
                                        EM_S[INDEX4(k,m,1,1,numEq,numComp,8)] = 4*wd0;
                                        EM_S[INDEX4(k,m,1,3,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,1,5,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,1,7,numEq,numComp,8)] =   wd0;
                                        EM_S[INDEX4(k,m,3,1,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,3,3,numEq,numComp,8)] = 4*wd0;
                                        EM_S[INDEX4(k,m,3,5,numEq,numComp,8)] =   wd0;
                                        EM_S[INDEX4(k,m,3,7,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,5,1,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,5,3,numEq,numComp,8)] =   wd0;
                                        EM_S[INDEX4(k,m,5,5,numEq,numComp,8)] = 4*wd0;
                                        EM_S[INDEX4(k,m,5,7,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,7,1,numEq,numComp,8)] =   wd0;
                                        EM_S[INDEX4(k,m,7,3,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,7,5,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,7,7,numEq,numComp,8)] = 4*wd0;
                                    }
                                }
                            }
                        }
                        ///////////////
                        // process y //
                        ///////////////
                        if (add_EM_F) {
                            const double* y_p=y.getSampleDataRO(e);
                            if (y.actsExpanded()) {
                                for (index_t k=0; k<numEq; k++) {
                                    const double y_0 = y_p[INDEX2(k, 0, numEq)];
                                    const double y_1 = y_p[INDEX2(k, 1, numEq)];
                                    const double y_2 = y_p[INDEX2(k, 2, numEq)];
                                    const double y_3 = y_p[INDEX2(k, 3, numEq)];
                                    const double tmp0 = 6*w2*(y_1 + y_2);
                                    const double tmp1 = 6*w2*(y_0 + y_3);
                                    EM_F[INDEX2(k,1,numEq)] = tmp0 + 6*w0*y_3 + 6*w1*y_0;
                                    EM_F[INDEX2(k,3,numEq)] = tmp1 + 6*w0*y_2 + 6*w1*y_1;
                                    EM_F[INDEX2(k,5,numEq)] = tmp1 + 6*w0*y_1 + 6*w1*y_2;
                                    EM_F[INDEX2(k,7,numEq)] = tmp0 + 6*w0*y_0 + 6*w1*y_3;
                                }
                            } else { // constant data
                                for (index_t k=0; k<numEq; k++) {
                                    EM_F[INDEX2(k,1,numEq)] = 36*w2*y_p[k];
                                    EM_F[INDEX2(k,3,numEq)] = 36*w2*y_p[k];
                                    EM_F[INDEX2(k,5,numEq)] = 36*w2*y_p[k];
                                    EM_F[INDEX2(k,7,numEq)] = 36*w2*y_p[k];
                                }
                            }
                        }
                        const index_t firstNode=m_NN[0]*m_NN[1]*k2+m_NN[0]*(k1+1)-2;
                        domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F,
                                add_EM_S, add_EM_F, firstNode, numEq, numComp);
                    } // k1 loop
                } // k2 loop
            } // colouring
        } // face 1

        if (domain->m_faceOffset[2] > -1) {
            if (add_EM_S)
                fill(EM_S.begin(), EM_S.end(), 0);
            if (add_EM_F)
                fill(EM_F.begin(), EM_F.end(), 0);

            for (index_t k2_0=0; k2_0<2; k2_0++) { // colouring
#pragma omp for
                for (index_t k2=k2_0; k2<m_NE[2]; k2+=2) {
                    for (index_t k0=0; k0<m_NE[0]; ++k0) {
                        const index_t e = domain->m_faceOffset[2]+INDEX2(k0,k2,m_NE[0]);
                        ///////////////
                        // process d //
                        ///////////////
                        if (add_EM_S) {
                            const double* d_p=d.getSampleDataRO(e);
                            if (d.actsExpanded()) {
                                for (index_t k=0; k<numEq; k++) {
                                    for (index_t m=0; m<numComp; m++) {
                                        const double d_0 = d_p[INDEX3(k,m,0,numEq,numComp)];
                                        const double d_1 = d_p[INDEX3(k,m,1,numEq,numComp)];
                                        const double d_2 = d_p[INDEX3(k,m,2,numEq,numComp)];
                                        const double d_3 = d_p[INDEX3(k,m,3,numEq,numComp)];
                                        const double tmp0 = w5*(d_0 + d_1);
                                        const double tmp1 = w6*(d_2 + d_3);
                                        const double tmp2 = w5*(d_0 + d_2);
                                        const double tmp3 = w6*(d_1 + d_3);
                                        const double tmp4 = w5*(d_1 + d_3);
                                        const double tmp5 = w6*(d_0 + d_2);
                                        const double tmp6 = w7*(d_0 + d_3);
                                        const double tmp7 = w7*(d_0 + d_1 + d_2 + d_3);
                                        const double tmp8 = w7*(d_1 + d_2);
                                        const double tmp9 = w5*(d_2 + d_3);
                                        const double tmp10 = w6*(d_0 + d_1);
                                        EM_S[INDEX4(k,m,0,0,numEq,numComp,8)] = d_0*w9 + d_3*w8 + tmp8;
                                        EM_S[INDEX4(k,m,0,1,numEq,numComp,8)] = tmp10 + tmp9;
                                        EM_S[INDEX4(k,m,0,4,numEq,numComp,8)] = tmp4 + tmp5;
                                        EM_S[INDEX4(k,m,0,5,numEq,numComp,8)] = tmp7;
                                        EM_S[INDEX4(k,m,1,0,numEq,numComp,8)] = tmp10 + tmp9;
                                        EM_S[INDEX4(k,m,1,1,numEq,numComp,8)] = d_1*w9 + d_2*w8 + tmp6;
                                        EM_S[INDEX4(k,m,1,4,numEq,numComp,8)] = tmp7;
                                        EM_S[INDEX4(k,m,1,5,numEq,numComp,8)] = tmp2 + tmp3;
                                        EM_S[INDEX4(k,m,4,0,numEq,numComp,8)] = tmp4 + tmp5;
                                        EM_S[INDEX4(k,m,4,1,numEq,numComp,8)] = tmp7;
                                        EM_S[INDEX4(k,m,4,4,numEq,numComp,8)] = d_1*w8 + d_2*w9 + tmp6;
                                        EM_S[INDEX4(k,m,4,5,numEq,numComp,8)] = tmp0 + tmp1;
                                        EM_S[INDEX4(k,m,5,0,numEq,numComp,8)] = tmp7;
                                        EM_S[INDEX4(k,m,5,1,numEq,numComp,8)] = tmp2 + tmp3;
                                        EM_S[INDEX4(k,m,5,4,numEq,numComp,8)] = tmp0 + tmp1;
                                        EM_S[INDEX4(k,m,5,5,numEq,numComp,8)] = d_0*w8 + d_3*w9 + tmp8;
                                    }
                                 }
                            } else { // constant data
                                for (index_t k=0; k<numEq; k++) {
                                    for (index_t m=0; m<numComp; m++) {
                                        const double wd0 = 4*d_p[INDEX2(k, m, numEq)]*w7;
                                        EM_S[INDEX4(k,m,0,0,numEq,numComp,8)] = 4*wd0;
                                        EM_S[INDEX4(k,m,0,1,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,0,4,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,0,5,numEq,numComp,8)] =   wd0;
                                        EM_S[INDEX4(k,m,1,0,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,1,1,numEq,numComp,8)] = 4*wd0;
                                        EM_S[INDEX4(k,m,1,4,numEq,numComp,8)] =   wd0;
                                        EM_S[INDEX4(k,m,1,5,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,4,0,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,4,1,numEq,numComp,8)] =   wd0;
                                        EM_S[INDEX4(k,m,4,4,numEq,numComp,8)] = 4*wd0;
                                        EM_S[INDEX4(k,m,4,5,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,5,0,numEq,numComp,8)] =   wd0;
                                        EM_S[INDEX4(k,m,5,1,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,5,4,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,5,5,numEq,numComp,8)] = 4*wd0;
                                    }
                                }
                            }
                        }
                        ///////////////
                        // process y //
                        ///////////////
                        if (add_EM_F) {
                            const double* y_p=y.getSampleDataRO(e);
                            if (y.actsExpanded()) {
                                for (index_t k=0; k<numEq; k++) {
                                    const double y_0 = y_p[INDEX2(k, 0, numEq)];
                                    const double y_1 = y_p[INDEX2(k, 1, numEq)];
                                    const double y_2 = y_p[INDEX2(k, 2, numEq)];
                                    const double y_3 = y_p[INDEX2(k, 3, numEq)];
                                    const double tmp0 = 6*w7*(y_1 + y_2);
                                    const double tmp1 = 6*w7*(y_0 + y_3);
                                    EM_F[INDEX2(k,0,numEq)] = tmp0 + 6*w5*y_3 + 6*w6*y_0;
                                    EM_F[INDEX2(k,1,numEq)] = tmp1 + 6*w5*y_2 + 6*w6*y_1;
                                    EM_F[INDEX2(k,4,numEq)] = tmp1 + 6*w5*y_1 + 6*w6*y_2;
                                    EM_F[INDEX2(k,5,numEq)] = tmp0 + 6*w5*y_0 + 6*w6*y_3;
                                }
                            } else { // constant data
                                for (index_t k=0; k<numEq; k++) {
                                    EM_F[INDEX2(k,0,numEq)] = 36*w7*y_p[k];
                                    EM_F[INDEX2(k,1,numEq)] = 36*w7*y_p[k];
                                    EM_F[INDEX2(k,4,numEq)] = 36*w7*y_p[k];
                                    EM_F[INDEX2(k,5,numEq)] = 36*w7*y_p[k];
                                }
                            }
                        }
                        const index_t firstNode=m_NN[0]*m_NN[1]*k2+k0;
                        domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F,
                                add_EM_S, add_EM_F, firstNode, numEq, numComp);
                    } // k0 loop
                } // k2 loop
            } // colouring
        } // face 2

        if (domain->m_faceOffset[3] > -1) {
            if (add_EM_S)
                fill(EM_S.begin(), EM_S.end(), 0);
            if (add_EM_F)
                fill(EM_F.begin(), EM_F.end(), 0);

            for (index_t k2_0=0; k2_0<2; k2_0++) { // colouring
#pragma omp for
                for (index_t k2=k2_0; k2<m_NE[2]; k2+=2) {
                    for (index_t k0=0; k0<m_NE[0]; ++k0) {
                        const index_t e = domain->m_faceOffset[3]+INDEX2(k0,k2,m_NE[0]);
                        ///////////////
                        // process d //
                        ///////////////
                        if (add_EM_S) {
                            const double* d_p=d.getSampleDataRO(e);
                            if (d.actsExpanded()) {
                                for (index_t k=0; k<numEq; k++) {
                                    for (index_t m=0; m<numComp; m++) {
                                        const double d_0 = d_p[INDEX3(k,m,0,numEq,numComp)];
                                        const double d_1 = d_p[INDEX3(k,m,1,numEq,numComp)];
                                        const double d_2 = d_p[INDEX3(k,m,2,numEq,numComp)];
                                        const double d_3 = d_p[INDEX3(k,m,3,numEq,numComp)];
                                        const double tmp0 = w5*(d_0 + d_2);
                                        const double tmp1 = w6*(d_1 + d_3);
                                        const double tmp2 = w5*(d_1 + d_3);
                                        const double tmp3 = w6*(d_0 + d_2);
                                        const double tmp4 = w7*(d_0 + d_1 + d_2 + d_3);
                                        const double tmp5 = w5*(d_0 + d_1);
                                        const double tmp6 = w6*(d_2 + d_3);
                                        const double tmp7 = w7*(d_0 + d_3);
                                        const double tmp8 = w7*(d_1 + d_2);
                                        const double tmp9 = w5*(d_2 + d_3);
                                        const double tmp10 = w6*(d_0 + d_1);
                                        EM_S[INDEX4(k,m,2,2,numEq,numComp,8)] = d_0*w9 + d_3*w8 + tmp8;
                                        EM_S[INDEX4(k,m,2,3,numEq,numComp,8)] = tmp10 + tmp9;
                                        EM_S[INDEX4(k,m,2,6,numEq,numComp,8)] = tmp2 + tmp3;
                                        EM_S[INDEX4(k,m,2,7,numEq,numComp,8)] = tmp4;
                                        EM_S[INDEX4(k,m,3,2,numEq,numComp,8)] = tmp10 + tmp9;
                                        EM_S[INDEX4(k,m,3,3,numEq,numComp,8)] = d_1*w9 + d_2*w8 + tmp7;
                                        EM_S[INDEX4(k,m,3,6,numEq,numComp,8)] = tmp4;
                                        EM_S[INDEX4(k,m,3,7,numEq,numComp,8)] = tmp0 + tmp1;
                                        EM_S[INDEX4(k,m,6,2,numEq,numComp,8)] = tmp2 + tmp3;
                                        EM_S[INDEX4(k,m,6,3,numEq,numComp,8)] = tmp4;
                                        EM_S[INDEX4(k,m,6,6,numEq,numComp,8)] = d_1*w8 + d_2*w9 + tmp7;
                                        EM_S[INDEX4(k,m,6,7,numEq,numComp,8)] = tmp5 + tmp6;
                                        EM_S[INDEX4(k,m,7,2,numEq,numComp,8)] = tmp4;
                                        EM_S[INDEX4(k,m,7,3,numEq,numComp,8)] = tmp0 + tmp1;
                                        EM_S[INDEX4(k,m,7,6,numEq,numComp,8)] = tmp5 + tmp6;
                                        EM_S[INDEX4(k,m,7,7,numEq,numComp,8)] = d_0*w8 + d_3*w9 + tmp8;
                                    }
                                 }
                            } else { // constant data
                                for (index_t k=0; k<numEq; k++) {
                                    for (index_t m=0; m<numComp; m++) {
                                        const double wd0 = 4*d_p[INDEX2(k, m, numEq)]*w7;
                                        EM_S[INDEX4(k,m,2,2,numEq,numComp,8)] = 4*wd0;
                                        EM_S[INDEX4(k,m,2,3,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,2,6,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,2,7,numEq,numComp,8)] =   wd0;
                                        EM_S[INDEX4(k,m,3,2,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,3,3,numEq,numComp,8)] = 4*wd0;
                                        EM_S[INDEX4(k,m,3,6,numEq,numComp,8)] =   wd0;
                                        EM_S[INDEX4(k,m,3,7,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,6,2,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,6,3,numEq,numComp,8)] =   wd0;
                                        EM_S[INDEX4(k,m,6,6,numEq,numComp,8)] = 4*wd0;
                                        EM_S[INDEX4(k,m,6,7,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,7,2,numEq,numComp,8)] =   wd0;
                                        EM_S[INDEX4(k,m,7,3,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,7,6,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,7,7,numEq,numComp,8)] = 4*wd0;
                                    }
                                }
                            }
                        }
                        ///////////////
                        // process y //
                        ///////////////
                        if (add_EM_F) {
                            const double* y_p=y.getSampleDataRO(e);
                            if (y.actsExpanded()) {
                                for (index_t k=0; k<numEq; k++) {
                                    const double y_0 = y_p[INDEX2(k, 0, numEq)];
                                    const double y_1 = y_p[INDEX2(k, 1, numEq)];
                                    const double y_2 = y_p[INDEX2(k, 2, numEq)];
                                    const double y_3 = y_p[INDEX2(k, 3, numEq)];
                                    const double tmp0 = 6*w7*(y_1 + y_2);
                                    const double tmp1 = 6*w7*(y_0 + y_3);
                                    EM_F[INDEX2(k,2,numEq)] = tmp0 + 6*w5*y_3 + 6*w6*y_0;
                                    EM_F[INDEX2(k,3,numEq)] = tmp1 + 6*w5*y_2 + 6*w6*y_1;
                                    EM_F[INDEX2(k,6,numEq)] = tmp1 + 6*w5*y_1 + 6*w6*y_2;
                                    EM_F[INDEX2(k,7,numEq)] = tmp0 + 6*w5*y_0 + 6*w6*y_3;
                                }
                            } else { // constant data
                                for (index_t k=0; k<numEq; k++) {
                                    EM_F[INDEX2(k,2,numEq)] = 36*w7*y_p[k];
                                    EM_F[INDEX2(k,3,numEq)] = 36*w7*y_p[k];
                                    EM_F[INDEX2(k,6,numEq)] = 36*w7*y_p[k];
                                    EM_F[INDEX2(k,7,numEq)] = 36*w7*y_p[k];
                                }
                            }
                        }
                        const index_t firstNode=m_NN[0]*m_NN[1]*k2+m_NN[0]*(m_NN[1]-2)+k0;
                        domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F,
                                add_EM_S, add_EM_F, firstNode, numEq, numComp);
                    } // k0 loop
                } // k2 loop
            } // colouring
        } // face 3

        if (domain->m_faceOffset[4] > -1) {
            if (add_EM_S)
                fill(EM_S.begin(), EM_S.end(), 0);
            if (add_EM_F)
                fill(EM_F.begin(), EM_F.end(), 0);

            for (index_t k1_0=0; k1_0<2; k1_0++) { // colouring
#pragma omp for
                for (index_t k1=k1_0; k1<m_NE[1]; k1+=2) {
                    for (index_t k0=0; k0<m_NE[0]; ++k0) {
                        const index_t e = domain->m_faceOffset[4]+INDEX2(k0,k1,m_NE[0]);
                        ///////////////
                        // process d //
                        ///////////////
                        if (add_EM_S) {
                            const double* d_p=d.getSampleDataRO(e);
                            if (d.actsExpanded()) {
                                for (index_t k=0; k<numEq; k++) {
                                    for (index_t m=0; m<numComp; m++) {
                                        const double d_0 = d_p[INDEX3(k,m,0,numEq,numComp)];
                                        const double d_1 = d_p[INDEX3(k,m,1,numEq,numComp)];
                                        const double d_2 = d_p[INDEX3(k,m,2,numEq,numComp)];
                                        const double d_3 = d_p[INDEX3(k,m,3,numEq,numComp)];
                                        const double tmp0 = w10*(d_0 + d_2);
                                        const double tmp1 = w11*(d_1 + d_3);
                                        const double tmp2 = w12*(d_0 + d_1 + d_2 + d_3);
                                        const double tmp3 = w12*(d_1 + d_2);
                                        const double tmp4 = w10*(d_1 + d_3);
                                        const double tmp5 = w11*(d_0 + d_2);
                                        const double tmp6 = w12*(d_0 + d_3);
                                        const double tmp7 = w10*(d_0 + d_1);
                                        const double tmp8 = w11*(d_2 + d_3);
                                        const double tmp9 = w10*(d_2 + d_3);
                                        const double tmp10 = w11*(d_0 + d_1);
                                        EM_S[INDEX4(k,m,0,0,numEq,numComp,8)] = d_0*w14 + d_3*w13 + tmp3;
                                        EM_S[INDEX4(k,m,0,1,numEq,numComp,8)] = tmp10 + tmp9;
                                        EM_S[INDEX4(k,m,0,2,numEq,numComp,8)] = tmp4 + tmp5;
                                        EM_S[INDEX4(k,m,0,3,numEq,numComp,8)] = tmp2;
                                        EM_S[INDEX4(k,m,1,0,numEq,numComp,8)] = tmp10 + tmp9;
                                        EM_S[INDEX4(k,m,1,1,numEq,numComp,8)] = d_1*w14 + d_2*w13 + tmp6;
                                        EM_S[INDEX4(k,m,1,2,numEq,numComp,8)] = tmp2;
                                        EM_S[INDEX4(k,m,1,3,numEq,numComp,8)] = tmp0 + tmp1;
                                        EM_S[INDEX4(k,m,2,0,numEq,numComp,8)] = tmp4 + tmp5;
                                        EM_S[INDEX4(k,m,2,1,numEq,numComp,8)] = tmp2;
                                        EM_S[INDEX4(k,m,2,2,numEq,numComp,8)] = d_1*w13 + d_2*w14 + tmp6;
                                        EM_S[INDEX4(k,m,2,3,numEq,numComp,8)] = tmp7 + tmp8;
                                        EM_S[INDEX4(k,m,3,0,numEq,numComp,8)] = tmp2;
                                        EM_S[INDEX4(k,m,3,1,numEq,numComp,8)] = tmp0 + tmp1;
                                        EM_S[INDEX4(k,m,3,2,numEq,numComp,8)] = tmp7 + tmp8;
                                        EM_S[INDEX4(k,m,3,3,numEq,numComp,8)] = d_0*w13 + d_3*w14 + tmp3;
                                    }
                                 }
                            } else { // constant data
                                for (index_t k=0; k<numEq; k++) {
                                    for (index_t m=0; m<numComp; m++) {
                                        const double wd0 = 4*d_p[INDEX2(k, m, numEq)]*w12;
                                        EM_S[INDEX4(k,m,0,0,numEq,numComp,8)] = 4*wd0;
                                        EM_S[INDEX4(k,m,0,1,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,0,2,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,0,3,numEq,numComp,8)] =   wd0;
                                        EM_S[INDEX4(k,m,1,0,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,1,1,numEq,numComp,8)] = 4*wd0;
                                        EM_S[INDEX4(k,m,1,2,numEq,numComp,8)] =   wd0;
                                        EM_S[INDEX4(k,m,1,3,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,2,0,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,2,1,numEq,numComp,8)] =   wd0;
                                        EM_S[INDEX4(k,m,2,2,numEq,numComp,8)] = 4*wd0;
                                        EM_S[INDEX4(k,m,2,3,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,3,0,numEq,numComp,8)] =   wd0;
                                        EM_S[INDEX4(k,m,3,1,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,3,2,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,3,3,numEq,numComp,8)] = 4*wd0;
                                    }
                                }
                            }
                        }
                        ///////////////
                        // process y //
                        ///////////////
                        if (add_EM_F) {
                            const double* y_p=y.getSampleDataRO(e);
                            if (y.actsExpanded()) {
                                for (index_t k=0; k<numEq; k++) {
                                    const double y_0 = y_p[INDEX2(k, 0, numEq)];
                                    const double y_1 = y_p[INDEX2(k, 1, numEq)];
                                    const double y_2 = y_p[INDEX2(k, 2, numEq)];
                                    const double y_3 = y_p[INDEX2(k, 3, numEq)];
                                    const double tmp0 = 6*w12*(y_1 + y_2);
                                    const double tmp1 = 6*w12*(y_0 + y_3);
                                    EM_F[INDEX2(k,0,numEq)] = tmp0 + 6*w10*y_3 + 6*w11*y_0;
                                    EM_F[INDEX2(k,1,numEq)] = tmp1 + 6*w10*y_2 + 6*w11*y_1;
                                    EM_F[INDEX2(k,2,numEq)] = tmp1 + 6*w10*y_1 + 6*w11*y_2;
                                    EM_F[INDEX2(k,3,numEq)] = tmp0 + 6*w10*y_0 + 6*w11*y_3;
                                }
                            } else { // constant data
                                for (index_t k=0; k<numEq; k++) {
                                    EM_F[INDEX2(k,0,numEq)] = 36*w12*y_p[k];
                                    EM_F[INDEX2(k,1,numEq)] = 36*w12*y_p[k];
                                    EM_F[INDEX2(k,2,numEq)] = 36*w12*y_p[k];
                                    EM_F[INDEX2(k,3,numEq)] = 36*w12*y_p[k];
                                }
                            }
                        }
                        const index_t firstNode=m_NN[0]*k1+k0;
                        domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F,
                                add_EM_S, add_EM_F, firstNode, numEq, numComp);
                    } // k0 loop
                } // k1 loop
            } // colouring
        } // face 4

        if (domain->m_faceOffset[5] > -1) {
            if (add_EM_S)
                fill(EM_S.begin(), EM_S.end(), 0);
            if (add_EM_F)
                fill(EM_F.begin(), EM_F.end(), 0);

            for (index_t k1_0=0; k1_0<2; k1_0++) { // colouring
#pragma omp for
                for (index_t k1=k1_0; k1<m_NE[1]; k1+=2) {
                    for (index_t k0=0; k0<m_NE[0]; ++k0) {
                        const index_t e = domain->m_faceOffset[5]+INDEX2(k0,k1,m_NE[0]);
                        ///////////////
                        // process d //
                        ///////////////
                        if (add_EM_S) {
                            const double* d_p=d.getSampleDataRO(e);
                            if (d.actsExpanded()) {
                                for (index_t k=0; k<numEq; k++) {
                                    for (index_t m=0; m<numComp; m++) {
                                        const double d_0 = d_p[INDEX3(k,m,0,numEq,numComp)];
                                        const double d_1 = d_p[INDEX3(k,m,1,numEq,numComp)];
                                        const double d_2 = d_p[INDEX3(k,m,2,numEq,numComp)];
                                        const double d_3 = d_p[INDEX3(k,m,3,numEq,numComp)];
                                        const double tmp0 = w12*(d_0 + d_1 + d_2 + d_3);
                                        const double tmp1 = w10*(d_1 + d_3);
                                        const double tmp2 = w11*(d_0 + d_2);
                                        const double tmp3 = w10*(d_2 + d_3);
                                        const double tmp4 = w11*(d_0 + d_1);
                                        const double tmp5 = w10*(d_0 + d_1);
                                        const double tmp6 = w11*(d_2 + d_3);
                                        const double tmp7 = w12*(d_1 + d_2);
                                        const double tmp8 = w10*(d_0 + d_2);
                                        const double tmp9 = w11*(d_1 + d_3);
                                        const double tmp10 = w12*(d_0 + d_3);
                                        EM_S[INDEX4(k,m,4,4,numEq,numComp,8)] = d_0*w14 + d_3*w13 + tmp7;
                                        EM_S[INDEX4(k,m,4,5,numEq,numComp,8)] = tmp3 + tmp4;
                                        EM_S[INDEX4(k,m,4,6,numEq,numComp,8)] = tmp1 + tmp2;
                                        EM_S[INDEX4(k,m,4,7,numEq,numComp,8)] = tmp0;
                                        EM_S[INDEX4(k,m,5,4,numEq,numComp,8)] = tmp3 + tmp4;
                                        EM_S[INDEX4(k,m,5,5,numEq,numComp,8)] = d_1*w14 + d_2*w13 + tmp10;
                                        EM_S[INDEX4(k,m,5,6,numEq,numComp,8)] = tmp0;
                                        EM_S[INDEX4(k,m,5,7,numEq,numComp,8)] = tmp8 + tmp9;
                                        EM_S[INDEX4(k,m,6,4,numEq,numComp,8)] = tmp1 + tmp2;
                                        EM_S[INDEX4(k,m,6,5,numEq,numComp,8)] = tmp0;
                                        EM_S[INDEX4(k,m,6,6,numEq,numComp,8)] = d_1*w13 + d_2*w14 + tmp10;
                                        EM_S[INDEX4(k,m,6,7,numEq,numComp,8)] = tmp5 + tmp6;
                                        EM_S[INDEX4(k,m,7,4,numEq,numComp,8)] = tmp0;
                                        EM_S[INDEX4(k,m,7,5,numEq,numComp,8)] = tmp8 + tmp9;
                                        EM_S[INDEX4(k,m,7,6,numEq,numComp,8)] = tmp5 + tmp6;
                                        EM_S[INDEX4(k,m,7,7,numEq,numComp,8)] = d_0*w13 + d_3*w14 + tmp7;
                                    }
                                 }
                            } else { // constant data
                                for (index_t k=0; k<numEq; k++) {
                                    for (index_t m=0; m<numComp; m++) {
                                        const double wd0 = 4*d_p[INDEX2(k, m, numEq)]*w12;
                                        EM_S[INDEX4(k,m,4,4,numEq,numComp,8)] = 4*wd0;
                                        EM_S[INDEX4(k,m,4,5,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,4,6,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,4,7,numEq,numComp,8)] =   wd0;
                                        EM_S[INDEX4(k,m,5,4,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,5,5,numEq,numComp,8)] = 4*wd0;
                                        EM_S[INDEX4(k,m,5,6,numEq,numComp,8)] =   wd0;
                                        EM_S[INDEX4(k,m,5,7,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,6,4,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,6,5,numEq,numComp,8)] =   wd0;
                                        EM_S[INDEX4(k,m,6,6,numEq,numComp,8)] = 4*wd0;
                                        EM_S[INDEX4(k,m,6,7,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,7,4,numEq,numComp,8)] =   wd0;
                                        EM_S[INDEX4(k,m,7,5,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,7,6,numEq,numComp,8)] = 2*wd0;
                                        EM_S[INDEX4(k,m,7,7,numEq,numComp,8)] = 4*wd0;
                                    }
                                }
                            }
                        }
                        ///////////////
                        // process y //
                        ///////////////
                        if (add_EM_F) {
                            const double* y_p=y.getSampleDataRO(e);
                            if (y.actsExpanded()) {
                                for (index_t k=0; k<numEq; k++) {
                                    const double y_0 = y_p[INDEX2(k, 0, numEq)];
                                    const double y_1 = y_p[INDEX2(k, 1, numEq)];
                                    const double y_2 = y_p[INDEX2(k, 2, numEq)];
                                    const double y_3 = y_p[INDEX2(k, 3, numEq)];
                                    const double tmp0 = 6*w12*(y_1 + y_2);
                                    const double tmp1 = 6*w12*(y_0 + y_3);
                                    EM_F[INDEX2(k,4,numEq)] = tmp0 + 6*w10*y_3 + 6*w11*y_0;
                                    EM_F[INDEX2(k,5,numEq)] = tmp1 + 6*w10*y_2 + 6*w11*y_1;
                                    EM_F[INDEX2(k,6,numEq)] = tmp1 + 6*w10*y_1 + 6*w11*y_2;
                                    EM_F[INDEX2(k,7,numEq)] = tmp0 + 6*w10*y_0 + 6*w11*y_3;
                                }
                            } else { // constant data
                                for (index_t k=0; k<numEq; k++) {
                                    EM_F[INDEX2(k,4,numEq)] = 36*w12*y_p[k];
                                    EM_F[INDEX2(k,5,numEq)] = 36*w12*y_p[k];
                                    EM_F[INDEX2(k,6,numEq)] = 36*w12*y_p[k];
                                    EM_F[INDEX2(k,7,numEq)] = 36*w12*y_p[k];
                                }
                            }
                        }
                        const index_t firstNode=m_NN[0]*m_NN[1]*(m_NN[2]-2)+m_NN[0]*k1+k0;
                        domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F,
                                add_EM_S, add_EM_F, firstNode, numEq, numComp);
                    } // k0 loop
                } // k1 loop
            } // colouring
        } // face 5
    } // end of parallel region
}

void LameAssembler3D::assemblePDESystem(AbstractSystemMatrix* mat, Data& rhs,
                                        const DataMap& coefs) const
{
    if (isNotEmpty("A", coefs))
        throw RipleyException("Coefficient A was given to LameAssembler "
                "unexpectedly. Specialised domains can't be used for general "
                "assemblage.");
    const Data& lambda = unpackData("lame_lambda", coefs);
    const Data& mu = unpackData("lame_mu", coefs);
    const Data& B = unpackData("B", coefs);
    const Data& C = unpackData("C", coefs);
    const Data& D = unpackData("D", coefs);
    const Data& X = unpackData("X", coefs);
    const Data& Y = unpackData("Y", coefs);

    dim_t numEq, numComp;
    if (!mat)
        numEq=numComp=(rhs.isEmpty() ? 1 : rhs.getDataPointSize());
    else {
        numEq=mat->getRowBlockSize();
        numComp=mat->getColumnBlockSize();
    }

    const double SQRT3 = 1.73205080756887719318;
    const double w10 = -m_dx[0]/288;
    const double w12 = w10*(-SQRT3 - 2);
    const double w6 = w10*(SQRT3 - 2);
    const double w18 = w10*(-4*SQRT3 - 7);
    const double w4 = w10*(-4*SQRT3 + 7);
    const double w11 = m_dx[1]/288;
    const double w15 = w11*(SQRT3 + 2);
    const double w5 = w11*(-SQRT3 + 2);
    const double w2 = w11*(4*SQRT3 - 7);
    const double w17 = w11*(4*SQRT3 + 7);
    const double w8 = m_dx[2]/288;
    const double w16 = w8*(SQRT3 + 2);
    const double w1 = w8*(-SQRT3 + 2);
    const double w20 = w8*(4*SQRT3 - 7);
    const double w21 = w8*(-4*SQRT3 - 7);
    const double w54 = -m_dx[0]*m_dx[1]/72;
    const double w68 = -m_dx[0]*m_dx[1]/48;
    const double w38 = w68*(-SQRT3 - 3)/36;
    const double w45 = w68*(SQRT3 - 3)/36;
    const double w35 = w68*(5*SQRT3 - 9)/36;
    const double w46 = w68*(-5*SQRT3 - 9)/36;
    const double w43 = w68*(-19*SQRT3 - 33)/36;
    const double w44 = w68*(19*SQRT3 - 33)/36;
    const double w66 = w68*(SQRT3 + 2);
    const double w70 = w68*(-SQRT3 + 2);
    const double w56 = -m_dx[0]*m_dx[2]/72;
    const double w67 = -m_dx[0]*m_dx[2]/48;
    const double w37 = w67*(-SQRT3 - 3)/36;
    const double w40 = w67*(SQRT3 - 3)/36;
    const double w34 = w67*(5*SQRT3 - 9)/36;
    const double w42 = w67*(-5*SQRT3 - 9)/36;
    const double w47 = w67*(19*SQRT3 + 33)/36;
    const double w48 = w67*(-19*SQRT3 + 33)/36;
    const double w65 = w67*(SQRT3 + 2);
    const double w71 = w67*(-SQRT3 + 2);
    const double w55 = -m_dx[1]*m_dx[2]/72;
    const double w69 = -m_dx[1]*m_dx[2]/48;
    const double w36 = w69*(SQRT3 - 3)/36;
    const double w39 = w69*(-SQRT3 - 3)/36;
    const double w33 = w69*(5*SQRT3 - 9)/36;
    const double w41 = w69*(-5*SQRT3 - 9)/36;
    const double w49 = w69*(19*SQRT3 - 33)/36;
    const double w50 = w69*(-19*SQRT3 - 33)/36;
    const double w64 = w69*(SQRT3 + 2);
    const double w72 = w69*(-SQRT3 + 2);
    const double w58 = m_dx[0]*m_dx[1]*m_dx[2]/1728;
    const double w60 = w58*(-SQRT3 + 2);
    const double w61 = w58*(SQRT3 + 2);
    const double w57 = w58*(-4*SQRT3 + 7);
    const double w59 = w58*(4*SQRT3 + 7);
    const double w62 = w58*(15*SQRT3 + 26);
    const double w63 = w58*(-15*SQRT3 + 26);
    const double w75 = w58*6*(SQRT3 + 3);
    const double w76 = w58*6*(-SQRT3 + 3);
    const double w74 = w58*6*(5*SQRT3 + 9);
    const double w77 = w58*6*(-5*SQRT3 + 9);
    const double w13 = -m_dx[0]*m_dx[1]/(288*m_dx[2]);
    const double w19 = w13*(4*SQRT3 + 7);
    const double w7 = w13*(-4*SQRT3 + 7);
    const double w23 = w13*(+SQRT3 - 2);
    const double w25 = w13*(-SQRT3 - 2);
    const double w22 = -m_dx[0]*m_dx[2]/(288*m_dx[1]);
    const double w3 = w22*(SQRT3 - 2);
    const double w9 = w22*(-SQRT3 - 2);
    const double w24 = w22*(4*SQRT3 + 7);
    const double w26 = w22*(-4*SQRT3 + 7);
    const double w27 = -m_dx[1]*m_dx[2]/(288*m_dx[0]);
    const double w0 = w27*(SQRT3 - 2);
    const double w14 = w27*(-SQRT3 - 2);
    const double w28 = w27*(-4*SQRT3 + 7);
    const double w29 = w27*(4*SQRT3 + 7);
    const bool add_EM_S = (!mu.isEmpty() || !lambda.isEmpty() || !B.isEmpty()
                                         || !C.isEmpty() || !D.isEmpty());
    const bool add_EM_F = (!X.isEmpty() || !Y.isEmpty());
    rhs.requireWrite();

#pragma omp parallel
    {
        vector<double> EM_S(8*8*numEq*numComp, 0);
        vector<double> EM_F(8*numEq, 0);

        for (index_t k2_0=0; k2_0<2; k2_0++) { // colouring
#pragma omp for
            for (index_t k2=k2_0; k2<m_NE[2]; k2+=2) {
                for (index_t k1=0; k1<m_NE[1]; ++k1) {
                    for (index_t k0=0; k0<m_NE[0]; ++k0)  {
                        const index_t e = k0 + m_NE[0]*k1 + m_NE[0]*m_NE[1]*k2;
                        if (add_EM_S)
                            fill(EM_S.begin(), EM_S.end(), 0);
                        if (add_EM_F)
                            fill(EM_F.begin(), EM_F.end(), 0);

                        ///////////////
                        // process A //
                        ///////////////
                        if (!mu.isEmpty() || !lambda.isEmpty()) {
                            if (mu.actsExpanded() || lambda.actsExpanded()) {
                                double A_0000[8] = {0};
                                double A_0011[8] = {0};
                                double A_0022[8] = {0};
                                double A_0101[8] = {0};
                                double A_0110[8] = {0};
                                double A_0202[8] = {0};
                                double A_0220[8] = {0};
                                double A_1001[8] = {0};
                                double A_1010[8] = {0};
                                double A_1100[8] = {0};
                                double A_1111[8] = {0};
                                double A_1122[8] = {0};
                                double A_1212[8] = {0};
                                double A_1221[8] = {0};
                                double A_2002[8] = {0};
                                double A_2020[8] = {0};
                                double A_2112[8] = {0};
                                double A_2121[8] = {0};
                                double A_2200[8] = {0};
                                double A_2211[8] = {0};
                                double A_2222[8] = {0};
                                if (!mu.isEmpty()) {
                                    const double *mu_p = mu.getSampleDataRO(e);
                                    A_0000[0] += 2*mu_p[0];
                                    A_0000[1] += 2*mu_p[1];
                                    A_0000[2] += 2*mu_p[2];
                                    A_0000[3] += 2*mu_p[3];
                                    A_0000[4] += 2*mu_p[4];
                                    A_0000[5] += 2*mu_p[5];
                                    A_0000[6] += 2*mu_p[6];
                                    A_0000[7] += 2*mu_p[7];
                                    A_0110[0] += mu_p[0];
                                    A_0101[0] += mu_p[0];
                                    A_0110[1] += mu_p[1];
                                    A_0101[1] += mu_p[1];
                                    A_0110[2] += mu_p[2];
                                    A_0101[2] += mu_p[2];
                                    A_0110[3] += mu_p[3];
                                    A_0101[3] += mu_p[3];
                                    A_0110[4] += mu_p[4];
                                    A_0101[4] += mu_p[4];
                                    A_0110[5] += mu_p[5];
                                    A_0101[5] += mu_p[5];
                                    A_0110[6] += mu_p[6];
                                    A_0101[6] += mu_p[6];
                                    A_0110[7] += mu_p[7];
                                    A_0101[7] += mu_p[7];
                                    A_0220[0] += mu_p[0];
                                    A_0202[0] += mu_p[0];
                                    A_0220[1] += mu_p[1];
                                    A_0202[1] += mu_p[1];
                                    A_0220[2] += mu_p[2];
                                    A_0202[2] += mu_p[2];
                                    A_0220[3] += mu_p[3];
                                    A_0202[3] += mu_p[3];
                                    A_0220[4] += mu_p[4];
                                    A_0202[4] += mu_p[4];
                                    A_0220[5] += mu_p[5];
                                    A_0202[5] += mu_p[5];
                                    A_0220[6] += mu_p[6];
                                    A_0202[6] += mu_p[6];
                                    A_0220[7] += mu_p[7];
                                    A_0202[7] += mu_p[7];
                                    A_1001[0] += mu_p[0];
                                    A_1010[0] += mu_p[0];
                                    A_1001[1] += mu_p[1];
                                    A_1010[1] += mu_p[1];
                                    A_1001[2] += mu_p[2];
                                    A_1010[2] += mu_p[2];
                                    A_1001[3] += mu_p[3];
                                    A_1010[3] += mu_p[3];
                                    A_1001[4] += mu_p[4];
                                    A_1010[4] += mu_p[4];
                                    A_1001[5] += mu_p[5];
                                    A_1010[5] += mu_p[5];
                                    A_1001[6] += mu_p[6];
                                    A_1010[6] += mu_p[6];
                                    A_1001[7] += mu_p[7];
                                    A_1010[7] += mu_p[7];
                                    A_1111[0] += 2*mu_p[0];
                                    A_1111[1] += 2*mu_p[1];
                                    A_1111[2] += 2*mu_p[2];
                                    A_1111[3] += 2*mu_p[3];
                                    A_1111[4] += 2*mu_p[4];
                                    A_1111[5] += 2*mu_p[5];
                                    A_1111[6] += 2*mu_p[6];
                                    A_1111[7] += 2*mu_p[7];
                                    A_1221[0] += mu_p[0];
                                    A_1212[0] += mu_p[0];
                                    A_1221[1] += mu_p[1];
                                    A_1212[1] += mu_p[1];
                                    A_1221[2] += mu_p[2];
                                    A_1212[2] += mu_p[2];
                                    A_1221[3] += mu_p[3];
                                    A_1212[3] += mu_p[3];
                                    A_1221[4] += mu_p[4];
                                    A_1212[4] += mu_p[4];
                                    A_1221[5] += mu_p[5];
                                    A_1212[5] += mu_p[5];
                                    A_1221[6] += mu_p[6];
                                    A_1212[6] += mu_p[6];
                                    A_1221[7] += mu_p[7];
                                    A_1212[7] += mu_p[7];
                                    A_2002[0] += mu_p[0];
                                    A_2020[0] += mu_p[0];
                                    A_2002[1] += mu_p[1];
                                    A_2020[1] += mu_p[1];
                                    A_2002[2] += mu_p[2];
                                    A_2020[2] += mu_p[2];
                                    A_2002[3] += mu_p[3];
                                    A_2020[3] += mu_p[3];
                                    A_2002[4] += mu_p[4];
                                    A_2020[4] += mu_p[4];
                                    A_2002[5] += mu_p[5];
                                    A_2020[5] += mu_p[5];
                                    A_2002[6] += mu_p[6];
                                    A_2020[6] += mu_p[6];
                                    A_2002[7] += mu_p[7];
                                    A_2020[7] += mu_p[7];
                                    A_2112[0] += mu_p[0];
                                    A_2121[0] += mu_p[0];
                                    A_2112[1] += mu_p[1];
                                    A_2121[1] += mu_p[1];
                                    A_2112[2] += mu_p[2];
                                    A_2121[2] += mu_p[2];
                                    A_2112[3] += mu_p[3];
                                    A_2121[3] += mu_p[3];
                                    A_2112[4] += mu_p[4];
                                    A_2121[4] += mu_p[4];
                                    A_2112[5] += mu_p[5];
                                    A_2121[5] += mu_p[5];
                                    A_2112[6] += mu_p[6];
                                    A_2121[6] += mu_p[6];
                                    A_2112[7] += mu_p[7];
                                    A_2121[7] += mu_p[7];
                                    A_2222[0] += 2*mu_p[0];
                                    A_2222[1] += 2*mu_p[1];
                                    A_2222[2] += 2*mu_p[2];
                                    A_2222[3] += 2*mu_p[3];
                                    A_2222[4] += 2*mu_p[4];
                                    A_2222[5] += 2*mu_p[5];
                                    A_2222[6] += 2*mu_p[6];
                                    A_2222[7] += 2*mu_p[7];
                                }
                                if (!lambda.isEmpty()) {
                                    const double *lambda_p = lambda.getSampleDataRO(e);
                                    A_0000[0] += lambda_p[0];
                                    A_0000[1] += lambda_p[1];
                                    A_0000[2] += lambda_p[2];
                                    A_0000[3] += lambda_p[3];
                                    A_0000[4] += lambda_p[4];
                                    A_0000[5] += lambda_p[5];
                                    A_0000[6] += lambda_p[6];
                                    A_0000[7] += lambda_p[7];
                                    A_0011[0] += lambda_p[0];
                                    A_0011[1] += lambda_p[1];
                                    A_0011[2] += lambda_p[2];
                                    A_0011[3] += lambda_p[3];
                                    A_0011[4] += lambda_p[4];
                                    A_0011[5] += lambda_p[5];
                                    A_0011[6] += lambda_p[6];
                                    A_0011[7] += lambda_p[7];
                                    A_0022[0] += lambda_p[0];
                                    A_0022[1] += lambda_p[1];
                                    A_0022[2] += lambda_p[2];
                                    A_0022[3] += lambda_p[3];
                                    A_0022[4] += lambda_p[4];
                                    A_0022[5] += lambda_p[5];
                                    A_0022[6] += lambda_p[6];
                                    A_0022[7] += lambda_p[7];
                                    A_1100[0] += lambda_p[0];
                                    A_1100[1] += lambda_p[1];
                                    A_1100[2] += lambda_p[2];
                                    A_1100[3] += lambda_p[3];
                                    A_1100[4] += lambda_p[4];
                                    A_1100[5] += lambda_p[5];
                                    A_1100[6] += lambda_p[6];
                                    A_1100[7] += lambda_p[7];
                                    A_1111[0] += lambda_p[0];
                                    A_1111[1] += lambda_p[1];
                                    A_1111[2] += lambda_p[2];
                                    A_1111[3] += lambda_p[3];
                                    A_1111[4] += lambda_p[4];
                                    A_1111[5] += lambda_p[5];
                                    A_1111[6] += lambda_p[6];
                                    A_1111[7] += lambda_p[7];
                                    A_1122[0] += lambda_p[0];
                                    A_1122[1] += lambda_p[1];
                                    A_1122[2] += lambda_p[2];
                                    A_1122[3] += lambda_p[3];
                                    A_1122[4] += lambda_p[4];
                                    A_1122[5] += lambda_p[5];
                                    A_1122[6] += lambda_p[6];
                                    A_1122[7] += lambda_p[7];
                                    A_2200[0] += lambda_p[0];
                                    A_2200[1] += lambda_p[1];
                                    A_2200[2] += lambda_p[2];
                                    A_2200[3] += lambda_p[3];
                                    A_2200[4] += lambda_p[4];
                                    A_2200[5] += lambda_p[5];
                                    A_2200[6] += lambda_p[6];
                                    A_2200[7] += lambda_p[7];
                                    A_2211[0] += lambda_p[0];
                                    A_2211[1] += lambda_p[1];
                                    A_2211[2] += lambda_p[2];
                                    A_2211[3] += lambda_p[3];
                                    A_2211[4] += lambda_p[4];
                                    A_2211[5] += lambda_p[5];
                                    A_2211[6] += lambda_p[6];
                                    A_2211[7] += lambda_p[7];
                                    A_2222[0] += lambda_p[0];
                                    A_2222[1] += lambda_p[1];
                                    A_2222[2] += lambda_p[2];
                                    A_2222[3] += lambda_p[3];
                                    A_2222[4] += lambda_p[4];
                                    A_2222[5] += lambda_p[5];
                                    A_2222[6] += lambda_p[6];
                                    A_2222[7] += lambda_p[7];
                                }
                                {
                                    const double tmp1 = w13*(A_0202[1] + A_0202[2] + A_0202[5] + A_0202[6]);
                                    const double tmp3 = w14*(A_0000[2] + A_0000[3] + A_0000[6] + A_0000[7]);
                                    const double tmp4 = w7*(A_0202[0] + A_0202[4]);
                                    const double tmp6 = w3*(A_0101[0] + A_0101[2] + A_0101[4] + A_0101[6]);
                                    const double tmp10 = w0*(A_0000[0] + A_0000[1] + A_0000[4] + A_0000[5]);
                                    const double tmp12 = w9*(A_0101[1] + A_0101[3] + A_0101[5] + A_0101[7]);
                                    const double tmp17 = w19*(A_0202[3] + A_0202[7]);
                                    const double tmp20 = w13*(-A_0202[0] - A_0202[1] - A_0202[2] - A_0202[3] - A_0202[4] - A_0202[5] - A_0202[6] - A_0202[7]);
                                    const double tmp22 = w14*(-A_0000[4] - A_0000[5] - A_0000[6] - A_0000[7]);
                                    const double tmp25 = w3*(-A_0101[0] - A_0101[1] - A_0101[2] - A_0101[3]);
                                    const double tmp28 = w0*(-A_0000[0] - A_0000[1] - A_0000[2] - A_0000[3]);
                                    const double tmp30 = w9*(-A_0101[4] - A_0101[5] - A_0101[6] - A_0101[7]);
                                    const double tmp39 = w14*(A_0000[0] + A_0000[1] + A_0000[2] + A_0000[3]);
                                    const double tmp40 = w26*(A_0101[4] + A_0101[6]);
                                    const double tmp41 = w0*(A_0000[4] + A_0000[5] + A_0000[6] + A_0000[7]);
                                    const double tmp43 = w22*(A_0101[0] + A_0101[2] + A_0101[5] + A_0101[7]);
                                    const double tmp45 = w25*(A_0202[1] + A_0202[3] + A_0202[5] + A_0202[7]);
                                    const double tmp52 = w24*(A_0101[1] + A_0101[3]);
                                    const double tmp55 = w23*(A_0202[0] + A_0202[2] + A_0202[4] + A_0202[6]);
                                    const double tmp57 = w14*(A_0000[4] + A_0000[5] + A_0000[6] + A_0000[7]);
                                    const double tmp58 = w26*(A_0101[1] + A_0101[3]);
                                    const double tmp61 = w25*(A_0202[0] + A_0202[2] + A_0202[4] + A_0202[6]);
                                    const double tmp64 = w0*(A_0000[0] + A_0000[1] + A_0000[2] + A_0000[3]);
                                    const double tmp66 = w24*(A_0101[4] + A_0101[6]);
                                    const double tmp71 = w23*(A_0202[1] + A_0202[3] + A_0202[5] + A_0202[7]);
                                    const double tmp73 = w14*(-A_0000[0] - A_0000[1] - A_0000[2] - A_0000[3]);
                                    const double tmp74 = w0*(-A_0000[4] - A_0000[5] - A_0000[6] - A_0000[7]);
                                    const double tmp75 = w3*(-A_0101[4] - A_0101[5] - A_0101[6] - A_0101[7]);
                                    const double tmp80 = w9*(-A_0101[0] - A_0101[1] - A_0101[2] - A_0101[3]);
                                    const double tmp88 = w3*(A_0101[0] + A_0101[1] + A_0101[2] + A_0101[3]);
                                    const double tmp89 = w23*(A_0202[2] + A_0202[3] + A_0202[6] + A_0202[7]);
                                    const double tmp91 = w25*(A_0202[0] + A_0202[1] + A_0202[4] + A_0202[5]);
                                    const double tmp95 = w28*(A_0000[2] + A_0000[3]);
                                    const double tmp97 = w29*(A_0000[4] + A_0000[5]);
                                    const double tmp100 = w9*(A_0101[4] + A_0101[5] + A_0101[6] + A_0101[7]);
                                    const double tmp101 = w27*(A_0000[0] + A_0000[1] + A_0000[6] + A_0000[7]);
                                    const double tmp104 = w13*(A_0202[0] + A_0202[1] + A_0202[2] + A_0202[3] + A_0202[4] + A_0202[5] + A_0202[6] + A_0202[7]);
                                    const double tmp106 = w22*(A_0101[0] + A_0101[1] + A_0101[2] + A_0101[3] + A_0101[4] + A_0101[5] + A_0101[6] + A_0101[7]);
                                    const double tmp113 = w27*(A_0000[0] + A_0000[1] + A_0000[2] + A_0000[3] + A_0000[4] + A_0000[5] + A_0000[6] + A_0000[7]);
                                    const double tmp123 = w13*(A_0202[0] + A_0202[3] + A_0202[4] + A_0202[7]);
                                    const double tmp125 = w7*(A_0202[1] + A_0202[5]);
                                    const double tmp127 = w3*(A_0101[1] + A_0101[3] + A_0101[5] + A_0101[7]);
                                    const double tmp131 = w9*(A_0101[0] + A_0101[2] + A_0101[4] + A_0101[6]);
                                    const double tmp132 = w19*(A_0202[2] + A_0202[6]);
                                    const double tmp141 = w14*(A_0000[0] + A_0000[1] + A_0000[4] + A_0000[5]);
                                    const double tmp142 = w7*(A_0202[2] + A_0202[6]);
                                    const double tmp146 = w0*(A_0000[2] + A_0000[3] + A_0000[6] + A_0000[7]);
                                    const double tmp149 = w19*(A_0202[1] + A_0202[5]);
                                    const double tmp175 = w14*(-A_0000[2] - A_0000[3] - A_0000[6] - A_0000[7]);
                                    const double tmp176 = w22*(-A_0101[0] - A_0101[1] - A_0101[2] - A_0101[3] - A_0101[4] - A_0101[5] - A_0101[6] - A_0101[7]);
                                    const double tmp178 = w25*(-A_0202[2] - A_0202[3] - A_0202[6] - A_0202[7]);
                                    const double tmp180 = w0*(-A_0000[0] - A_0000[1] - A_0000[4] - A_0000[5]);
                                    const double tmp187 = w23*(-A_0202[0] - A_0202[1] - A_0202[4] - A_0202[5]);
                                    const double tmp189 = w7*(A_0202[3] + A_0202[7]);
                                    const double tmp193 = w19*(A_0202[0] + A_0202[4]);
                                    const double tmp201 = w27*(A_0000[2] + A_0000[3] + A_0000[4] + A_0000[5]);
                                    const double tmp204 = w23*(A_0202[0] + A_0202[1] + A_0202[4] + A_0202[5]);
                                    const double tmp205 = w25*(A_0202[2] + A_0202[3] + A_0202[6] + A_0202[7]);
                                    const double tmp208 = w28*(A_0000[0] + A_0000[1]);
                                    const double tmp209 = w29*(A_0000[6] + A_0000[7]);
                                    const double tmp214 = w13*(-A_0202[1] - A_0202[2] - A_0202[5] - A_0202[6]);
                                    const double tmp215 = w22*(-A_0101[0] - A_0101[2] - A_0101[5] - A_0101[7]);
                                    const double tmp217 = w27*(-A_0000[0] - A_0000[1] - A_0000[6] - A_0000[7]);
                                    const double tmp221 = w26*(-A_0101[4] - A_0101[6]);
                                    const double tmp226 = w7*(-A_0202[0] - A_0202[4]);
                                    const double tmp227 = w24*(-A_0101[1] - A_0101[3]);
                                    const double tmp228 = w19*(-A_0202[3] - A_0202[7]);
                                    const double tmp231 = w28*(-A_0000[4] - A_0000[5]);
                                    const double tmp233 = w29*(-A_0000[2] - A_0000[3]);
                                    const double tmp236 = w26*(A_0101[5] + A_0101[7]);
                                    const double tmp238 = w22*(A_0101[1] + A_0101[3] + A_0101[4] + A_0101[6]);
                                    const double tmp244 = w24*(A_0101[0] + A_0101[2]);
                                    const double tmp255 = w26*(-A_0101[1] - A_0101[3]);
                                    const double tmp259 = w7*(-A_0202[3] - A_0202[7]);
                                    const double tmp261 = w24*(-A_0101[4] - A_0101[6]);
                                    const double tmp262 = w19*(-A_0202[0] - A_0202[4]);
                                    const double tmp265 = w28*(-A_0000[2] - A_0000[3]);
                                    const double tmp268 = w29*(-A_0000[4] - A_0000[5]);
                                    const double tmp288 = w13*(-A_0202[0] - A_0202[3] - A_0202[4] - A_0202[7]);
                                    const double tmp289 = w22*(-A_0101[1] - A_0101[3] - A_0101[4] - A_0101[6]);
                                    const double tmp294 = w26*(-A_0101[5] - A_0101[7]);
                                    const double tmp298 = w7*(-A_0202[1] - A_0202[5]);
                                    const double tmp299 = w24*(-A_0101[0] - A_0101[2]);
                                    const double tmp300 = w19*(-A_0202[2] - A_0202[6]);
                                    const double tmp304 = w27*(-A_0000[2] - A_0000[3] - A_0000[4] - A_0000[5]);
                                    const double tmp308 = w26*(-A_0101[0] - A_0101[2]);
                                    const double tmp313 = w24*(-A_0101[5] - A_0101[7]);
                                    const double tmp316 = w28*(-A_0000[0] - A_0000[1]);
                                    const double tmp318 = w29*(-A_0000[6] - A_0000[7]);
                                    const double tmp320 = w26*(A_0101[0] + A_0101[2]);
                                    const double tmp325 = w24*(A_0101[5] + A_0101[7]);
                                    const double tmp329 = w3*(-A_0101[0] - A_0101[2] - A_0101[4] - A_0101[6]);
                                    const double tmp332 = w25*(-A_0202[1] - A_0202[3] - A_0202[5] - A_0202[7]);
                                    const double tmp335 = w9*(-A_0101[1] - A_0101[3] - A_0101[5] - A_0101[7]);
                                    const double tmp337 = w27*(-A_0000[0] - A_0000[1] - A_0000[2] - A_0000[3] - A_0000[4] - A_0000[5] - A_0000[6] - A_0000[7]);
                                    const double tmp338 = w23*(-A_0202[0] - A_0202[2] - A_0202[4] - A_0202[6]);
                                    const double tmp339 = w14*(-A_0000[0] - A_0000[1] - A_0000[4] - A_0000[5]);
                                    const double tmp340 = w23*(-A_0202[2] - A_0202[3] - A_0202[6] - A_0202[7]);
                                    const double tmp342 = w25*(-A_0202[0] - A_0202[1] - A_0202[4] - A_0202[5]);
                                    const double tmp344 = w0*(-A_0000[2] - A_0000[3] - A_0000[6] - A_0000[7]);
                                    const double tmp358 = w7*(-A_0202[2] - A_0202[6]);
                                    const double tmp359 = w19*(-A_0202[1] - A_0202[5]);
                                    const double tmp362 = w28*(-A_0000[6] - A_0000[7]);
                                    const double tmp363 = w29*(-A_0000[0] - A_0000[1]);
                                    const double tmp371 = w3*(A_0101[4] + A_0101[5] + A_0101[6] + A_0101[7]);
                                    const double tmp374 = w9*(A_0101[0] + A_0101[1] + A_0101[2] + A_0101[3]);
                                    const double tmp375 = w29*(A_0000[2] + A_0000[3]);
                                    const double tmp377 = w28*(A_0000[4] + A_0000[5]);
                                    const double tmp422 = w3*(-A_0101[1] - A_0101[3] - A_0101[5] - A_0101[7]);
                                    const double tmp424 = w25*(-A_0202[0] - A_0202[2] - A_0202[4] - A_0202[6]);
                                    const double tmp428 = w9*(-A_0101[0] - A_0101[2] - A_0101[4] - A_0101[6]);
                                    const double tmp430 = w23*(-A_0202[1] - A_0202[3] - A_0202[5] - A_0202[7]);
                                    const double tmp455 = w29*(A_0000[0] + A_0000[1]);
                                    const double tmp456 = w28*(A_0000[6] + A_0000[7]);
                                    EM_S[INDEX4(0,0,0,0,numEq,numComp,8)]+= tmp214 + tmp259 + tmp262 + tmp289 + tmp294 + tmp299 + tmp304 + tmp362 + tmp363;
                                    EM_S[INDEX4(0,0,0,1,numEq,numComp,8)]+= tmp201 + tmp371 + tmp374 + tmp455 + tmp456 + tmp89 + tmp91;
                                    EM_S[INDEX4(0,0,0,2,numEq,numComp,8)]+= tmp236 + tmp238 + tmp244 + tmp39 + tmp41 + tmp61 + tmp71;
                                    EM_S[INDEX4(0,0,0,3,numEq,numComp,8)]+= tmp20 + tmp73 + tmp74 + tmp75 + tmp80;
                                    EM_S[INDEX4(0,0,0,4,numEq,numComp,8)]+= tmp1 + tmp127 + tmp131 + tmp141 + tmp146 + tmp189 + tmp193;
                                    EM_S[INDEX4(0,0,0,5,numEq,numComp,8)]+= tmp176 + tmp339 + tmp340 + tmp342 + tmp344;
                                    EM_S[INDEX4(0,0,0,6,numEq,numComp,8)]+= tmp337 + tmp422 + tmp424 + tmp428 + tmp430;
                                    EM_S[INDEX4(0,0,0,7,numEq,numComp,8)]+= tmp104 + tmp106 + tmp113;
                                    EM_S[INDEX4(0,0,1,0,numEq,numComp,8)]+= tmp201 + tmp371 + tmp374 + tmp455 + tmp456 + tmp89 + tmp91;
                                    EM_S[INDEX4(0,0,1,1,numEq,numComp,8)]+= tmp215 + tmp221 + tmp227 + tmp288 + tmp304 + tmp358 + tmp359 + tmp362 + tmp363;
                                    EM_S[INDEX4(0,0,1,2,numEq,numComp,8)]+= tmp20 + tmp73 + tmp74 + tmp75 + tmp80;
                                    EM_S[INDEX4(0,0,1,3,numEq,numComp,8)]+= tmp39 + tmp40 + tmp41 + tmp43 + tmp45 + tmp52 + tmp55;
                                    EM_S[INDEX4(0,0,1,4,numEq,numComp,8)]+= tmp176 + tmp339 + tmp340 + tmp342 + tmp344;
                                    EM_S[INDEX4(0,0,1,5,numEq,numComp,8)]+= tmp12 + tmp123 + tmp141 + tmp142 + tmp146 + tmp149 + tmp6;
                                    EM_S[INDEX4(0,0,1,6,numEq,numComp,8)]+= tmp104 + tmp106 + tmp113;
                                    EM_S[INDEX4(0,0,1,7,numEq,numComp,8)]+= tmp329 + tmp332 + tmp335 + tmp337 + tmp338;
                                    EM_S[INDEX4(0,0,2,0,numEq,numComp,8)]+= tmp236 + tmp238 + tmp244 + tmp39 + tmp41 + tmp61 + tmp71;
                                    EM_S[INDEX4(0,0,2,1,numEq,numComp,8)]+= tmp20 + tmp73 + tmp74 + tmp75 + tmp80;
                                    EM_S[INDEX4(0,0,2,2,numEq,numComp,8)]+= tmp217 + tmp231 + tmp233 + tmp288 + tmp289 + tmp294 + tmp298 + tmp299 + tmp300;
                                    EM_S[INDEX4(0,0,2,3,numEq,numComp,8)]+= tmp101 + tmp204 + tmp205 + tmp371 + tmp374 + tmp375 + tmp377;
                                    EM_S[INDEX4(0,0,2,4,numEq,numComp,8)]+= tmp337 + tmp422 + tmp424 + tmp428 + tmp430;
                                    EM_S[INDEX4(0,0,2,5,numEq,numComp,8)]+= tmp104 + tmp106 + tmp113;
                                    EM_S[INDEX4(0,0,2,6,numEq,numComp,8)]+= tmp10 + tmp123 + tmp125 + tmp127 + tmp131 + tmp132 + tmp3;
                                    EM_S[INDEX4(0,0,2,7,numEq,numComp,8)]+= tmp175 + tmp176 + tmp178 + tmp180 + tmp187;
                                    EM_S[INDEX4(0,0,3,0,numEq,numComp,8)]+= tmp20 + tmp73 + tmp74 + tmp75 + tmp80;
                                    EM_S[INDEX4(0,0,3,1,numEq,numComp,8)]+= tmp39 + tmp40 + tmp41 + tmp43 + tmp45 + tmp52 + tmp55;
                                    EM_S[INDEX4(0,0,3,2,numEq,numComp,8)]+= tmp101 + tmp204 + tmp205 + tmp371 + tmp374 + tmp375 + tmp377;
                                    EM_S[INDEX4(0,0,3,3,numEq,numComp,8)]+= tmp214 + tmp215 + tmp217 + tmp221 + tmp226 + tmp227 + tmp228 + tmp231 + tmp233;
                                    EM_S[INDEX4(0,0,3,4,numEq,numComp,8)]+= tmp104 + tmp106 + tmp113;
                                    EM_S[INDEX4(0,0,3,5,numEq,numComp,8)]+= tmp329 + tmp332 + tmp335 + tmp337 + tmp338;
                                    EM_S[INDEX4(0,0,3,6,numEq,numComp,8)]+= tmp175 + tmp176 + tmp178 + tmp180 + tmp187;
                                    EM_S[INDEX4(0,0,3,7,numEq,numComp,8)]+= tmp1 + tmp10 + tmp12 + tmp17 + tmp3 + tmp4 + tmp6;
                                    EM_S[INDEX4(0,0,4,0,numEq,numComp,8)]+= tmp1 + tmp127 + tmp131 + tmp141 + tmp146 + tmp189 + tmp193;
                                    EM_S[INDEX4(0,0,4,1,numEq,numComp,8)]+= tmp176 + tmp339 + tmp340 + tmp342 + tmp344;
                                    EM_S[INDEX4(0,0,4,2,numEq,numComp,8)]+= tmp337 + tmp422 + tmp424 + tmp428 + tmp430;
                                    EM_S[INDEX4(0,0,4,3,numEq,numComp,8)]+= tmp104 + tmp106 + tmp113;
                                    EM_S[INDEX4(0,0,4,4,numEq,numComp,8)]+= tmp214 + tmp215 + tmp217 + tmp255 + tmp259 + tmp261 + tmp262 + tmp265 + tmp268;
                                    EM_S[INDEX4(0,0,4,5,numEq,numComp,8)]+= tmp100 + tmp101 + tmp88 + tmp89 + tmp91 + tmp95 + tmp97;
                                    EM_S[INDEX4(0,0,4,6,numEq,numComp,8)]+= tmp43 + tmp57 + tmp58 + tmp61 + tmp64 + tmp66 + tmp71;
                                    EM_S[INDEX4(0,0,4,7,numEq,numComp,8)]+= tmp20 + tmp22 + tmp25 + tmp28 + tmp30;
                                    EM_S[INDEX4(0,0,5,0,numEq,numComp,8)]+= tmp176 + tmp339 + tmp340 + tmp342 + tmp344;
                                    EM_S[INDEX4(0,0,5,1,numEq,numComp,8)]+= tmp12 + tmp123 + tmp141 + tmp142 + tmp146 + tmp149 + tmp6;
                                    EM_S[INDEX4(0,0,5,2,numEq,numComp,8)]+= tmp104 + tmp106 + tmp113;
                                    EM_S[INDEX4(0,0,5,3,numEq,numComp,8)]+= tmp329 + tmp332 + tmp335 + tmp337 + tmp338;
                                    EM_S[INDEX4(0,0,5,4,numEq,numComp,8)]+= tmp100 + tmp101 + tmp88 + tmp89 + tmp91 + tmp95 + tmp97;
                                    EM_S[INDEX4(0,0,5,5,numEq,numComp,8)]+= tmp217 + tmp265 + tmp268 + tmp288 + tmp289 + tmp308 + tmp313 + tmp358 + tmp359;
                                    EM_S[INDEX4(0,0,5,6,numEq,numComp,8)]+= tmp20 + tmp22 + tmp25 + tmp28 + tmp30;
                                    EM_S[INDEX4(0,0,5,7,numEq,numComp,8)]+= tmp238 + tmp320 + tmp325 + tmp45 + tmp55 + tmp57 + tmp64;
                                    EM_S[INDEX4(0,0,6,0,numEq,numComp,8)]+= tmp337 + tmp422 + tmp424 + tmp428 + tmp430;
                                    EM_S[INDEX4(0,0,6,1,numEq,numComp,8)]+= tmp104 + tmp106 + tmp113;
                                    EM_S[INDEX4(0,0,6,2,numEq,numComp,8)]+= tmp10 + tmp123 + tmp125 + tmp127 + tmp131 + tmp132 + tmp3;
                                    EM_S[INDEX4(0,0,6,3,numEq,numComp,8)]+= tmp175 + tmp176 + tmp178 + tmp180 + tmp187;
                                    EM_S[INDEX4(0,0,6,4,numEq,numComp,8)]+= tmp43 + tmp57 + tmp58 + tmp61 + tmp64 + tmp66 + tmp71;
                                    EM_S[INDEX4(0,0,6,5,numEq,numComp,8)]+= tmp20 + tmp22 + tmp25 + tmp28 + tmp30;
                                    EM_S[INDEX4(0,0,6,6,numEq,numComp,8)]+= tmp215 + tmp255 + tmp261 + tmp288 + tmp298 + tmp300 + tmp304 + tmp316 + tmp318;
                                    EM_S[INDEX4(0,0,6,7,numEq,numComp,8)]+= tmp100 + tmp201 + tmp204 + tmp205 + tmp208 + tmp209 + tmp88;
                                    EM_S[INDEX4(0,0,7,0,numEq,numComp,8)]+= tmp104 + tmp106 + tmp113;
                                    EM_S[INDEX4(0,0,7,1,numEq,numComp,8)]+= tmp329 + tmp332 + tmp335 + tmp337 + tmp338;
                                    EM_S[INDEX4(0,0,7,2,numEq,numComp,8)]+= tmp175 + tmp176 + tmp178 + tmp180 + tmp187;
                                    EM_S[INDEX4(0,0,7,3,numEq,numComp,8)]+= tmp1 + tmp10 + tmp12 + tmp17 + tmp3 + tmp4 + tmp6;
                                    EM_S[INDEX4(0,0,7,4,numEq,numComp,8)]+= tmp20 + tmp22 + tmp25 + tmp28 + tmp30;
                                    EM_S[INDEX4(0,0,7,5,numEq,numComp,8)]+= tmp238 + tmp320 + tmp325 + tmp45 + tmp55 + tmp57 + tmp64;
                                    EM_S[INDEX4(0,0,7,6,numEq,numComp,8)]+= tmp100 + tmp201 + tmp204 + tmp205 + tmp208 + tmp209 + tmp88;
                                    EM_S[INDEX4(0,0,7,7,numEq,numComp,8)]+= tmp214 + tmp226 + tmp228 + tmp289 + tmp304 + tmp308 + tmp313 + tmp316 + tmp318;
                                }
                                {
                                    const double tmp7 = w1*(A_0011[0] + A_0011[4] + A_0110[0] + A_0110[4]);
                                    const double tmp11 = w16*(A_0011[3] + A_0011[7] + A_0110[3] + A_0110[7]);
                                    const double tmp15 = w8*(A_0011[1] + A_0011[2] + A_0011[5] + A_0011[6] + A_0110[1] + A_0110[2] + A_0110[5] + A_0110[6]);
                                    const double tmp26 = w1*(-A_0011[0] - A_0011[3] - A_0110[0] - A_0110[3]);
                                    const double tmp29 = w16*(-A_0011[4] - A_0011[7] - A_0110[4] - A_0110[7]);
                                    const double tmp34 = w8*(-A_0011[1] - A_0011[6] - A_0110[2] - A_0110[5]);
                                    const double tmp44 = w1*(A_0011[4] + A_0011[7] - A_0110[5] - A_0110[6]);
                                    const double tmp49 = w16*(A_0011[0] + A_0011[3] - A_0110[1] - A_0110[2]);
                                    const double tmp53 = w8*(A_0011[2] + A_0011[5] - A_0110[0] - A_0110[7]);
                                    const double tmp60 = w1*(A_0011[0] + A_0011[3] - A_0110[1] - A_0110[2]);
                                    const double tmp65 = w16*(A_0011[4] + A_0011[7] - A_0110[5] - A_0110[6]);
                                    const double tmp76 = w1*(-A_0011[4] - A_0011[7] - A_0110[4] - A_0110[7]);
                                    const double tmp79 = w16*(-A_0011[0] - A_0011[3] - A_0110[0] - A_0110[3]);
                                    const double tmp90 = w1*(-A_0011[1] - A_0011[2] + A_0110[0] + A_0110[3]);
                                    const double tmp94 = w16*(-A_0011[5] - A_0011[6] + A_0110[4] + A_0110[7]);
                                    const double tmp99 = w8*(-A_0011[0] - A_0011[7] + A_0110[1] + A_0110[6]);
                                    const double tmp107 = w1*(-A_0011[2] - A_0011[6] - A_0110[1] - A_0110[5]);
                                    const double tmp109 = w16*(-A_0011[1] - A_0011[5] - A_0110[2] - A_0110[6]);
                                    const double tmp112 = w8*(-A_0011[0] - A_0011[3] - A_0011[4] - A_0011[7] - A_0110[0] - A_0110[3] - A_0110[4] - A_0110[7]);
                                    const double tmp118 = w16*(A_0011[5] + A_0011[6] + A_0110[5] + A_0110[6]);
                                    const double tmp120 = w8*(A_0011[0] + A_0011[7] + A_0110[3] + A_0110[4]);
                                    const double tmp121 = w1*(A_0011[1] + A_0011[2] + A_0110[1] + A_0110[2]);
                                    const double tmp128 = w1*(-A_0011[1] - A_0011[5] - A_0110[1] - A_0110[5]);
                                    const double tmp130 = w16*(-A_0011[2] - A_0011[6] - A_0110[2] - A_0110[6]);
                                    const double tmp136 = w1*(A_0011[3] + A_0011[7] + A_0110[0] + A_0110[4]);
                                    const double tmp138 = w16*(A_0011[0] + A_0011[4] + A_0110[3] + A_0110[7]);
                                    const double tmp143 = w1*(-A_0011[2] - A_0011[6] - A_0110[2] - A_0110[6]);
                                    const double tmp147 = w16*(-A_0011[1] - A_0011[5] - A_0110[1] - A_0110[5]);
                                    const double tmp162 = w1*(A_0011[0] + A_0011[4] + A_0110[3] + A_0110[7]);
                                    const double tmp163 = w16*(A_0011[3] + A_0011[7] + A_0110[0] + A_0110[4]);
                                    const double tmp171 = w8*(-A_0011[2] - A_0011[5] - A_0110[1] - A_0110[6]);
                                    const double tmp177 = w1*(A_0011[1] + A_0011[5] - A_0110[0] - A_0110[4]);
                                    const double tmp181 = w16*(A_0011[2] + A_0011[6] - A_0110[3] - A_0110[7]);
                                    const double tmp184 = w8*(A_0011[0] + A_0011[3] + A_0011[4] + A_0011[7] - A_0110[1] - A_0110[2] - A_0110[5] - A_0110[6]);
                                    const double tmp190 = w1*(A_0011[3] + A_0011[7] + A_0110[3] + A_0110[7]);
                                    const double tmp192 = w16*(A_0011[0] + A_0011[4] + A_0110[0] + A_0110[4]);
                                    const double tmp198 = w16*(A_0011[1] + A_0011[2] + A_0110[1] + A_0110[2]);
                                    const double tmp199 = w8*(A_0011[3] + A_0011[4] + A_0110[0] + A_0110[7]);
                                    const double tmp200 = w1*(A_0011[5] + A_0011[6] + A_0110[5] + A_0110[6]);
                                    const double tmp210 = w8*(-A_0011[3] - A_0011[4] + A_0110[2] + A_0110[5]);
                                    const double tmp216 = w8*(A_0011[0] + A_0011[7] + A_0110[0] + A_0110[7]);
                                    const double tmp245 = w8*(A_0011[1] + A_0011[6] - A_0110[3] - A_0110[4]);
                                    const double tmp250 = w8*(A_0011[2] + A_0011[5] - A_0110[3] - A_0110[4]);
                                    const double tmp270 = w1*(-A_0011[0] - A_0011[4] + A_0110[1] + A_0110[5]);
                                    const double tmp272 = w16*(-A_0011[3] - A_0011[7] + A_0110[2] + A_0110[6]);
                                    const double tmp274 = w8*(-A_0011[1] - A_0011[2] - A_0011[5] - A_0011[6] + A_0110[0] + A_0110[3] + A_0110[4] + A_0110[7]);
                                    const double tmp290 = w8*(-A_0011[1] - A_0011[6] - A_0110[1] - A_0110[6]);
                                    const double tmp303 = w8*(A_0011[3] + A_0011[4] + A_0110[3] + A_0110[4]);
                                    const double tmp330 = w1*(A_0011[2] + A_0011[6] - A_0110[0] - A_0110[4]);
                                    const double tmp334 = w16*(A_0011[1] + A_0011[5] - A_0110[3] - A_0110[7]);
                                    const double tmp341 = w1*(A_0011[2] + A_0011[6] - A_0110[3] - A_0110[7]);
                                    const double tmp345 = w16*(A_0011[1] + A_0011[5] - A_0110[0] - A_0110[4]);
                                    const double tmp351 = w8*(-A_0011[2] - A_0011[5] - A_0110[2] - A_0110[5]);
                                    const double tmp376 = w8*(A_0011[1] + A_0011[6] - A_0110[0] - A_0110[7]);
                                    const double tmp394 = w1*(-A_0011[3] - A_0011[7] + A_0110[2] + A_0110[6]);
                                    const double tmp395 = w16*(-A_0011[0] - A_0011[4] + A_0110[1] + A_0110[5]);
                                    const double tmp399 = w1*(-A_0011[0] - A_0011[4] + A_0110[2] + A_0110[6]);
                                    const double tmp401 = w16*(-A_0011[3] - A_0011[7] + A_0110[1] + A_0110[5]);
                                    const double tmp423 = w1*(A_0011[1] + A_0011[5] - A_0110[3] - A_0110[7]);
                                    const double tmp427 = w16*(A_0011[2] + A_0011[6] - A_0110[0] - A_0110[4]);
                                    const double tmp436 = w8*(-A_0011[3] - A_0011[4] + A_0110[1] + A_0110[6]);
                                    const double tmp440 = w16*(-A_0011[1] - A_0011[2] + A_0110[0] + A_0110[3]);
                                    const double tmp441 = w1*(-A_0011[5] - A_0011[6] + A_0110[4] + A_0110[7]);
                                    const double tmp447 = w1*(-A_0011[3] - A_0011[7] + A_0110[1] + A_0110[5]);
                                    const double tmp449 = w16*(-A_0011[0] - A_0011[4] + A_0110[2] + A_0110[6]);
                                    const double tmp471 = w1*(-A_0011[1] - A_0011[5] - A_0110[2] - A_0110[6]);
                                    const double tmp473 = w16*(-A_0011[2] - A_0011[6] - A_0110[1] - A_0110[5]);
                                    const double tmp481 = w8*(-A_0011[0] - A_0011[7] + A_0110[2] + A_0110[5]);
                                    EM_S[INDEX4(0,1,0,0,numEq,numComp,8)]+= tmp198 + tmp200 + tmp303 + w20*(-A_0011[7] - A_0110[7]) + w21*(-A_0011[0] - A_0110[0]);
                                    EM_S[INDEX4(0,1,0,1,numEq,numComp,8)]+= tmp250 + tmp44 + w20*(-A_0011[6] + A_0110[7]) + w21*(-A_0011[1] + A_0110[0]) + tmp49;
                                    EM_S[INDEX4(0,1,0,2,numEq,numComp,8)]+= tmp436 + tmp440 + tmp441 + w20*(A_0011[7] - A_0110[5]) + w21*(A_0011[0] - A_0110[2]);
                                    EM_S[INDEX4(0,1,0,3,numEq,numComp,8)]+= w20*(A_0011[6] + A_0110[5]) + w21*(A_0011[1] + A_0110[2]) + tmp171 + tmp76 + tmp79;
                                    EM_S[INDEX4(0,1,0,4,numEq,numComp,8)]+= tmp15 + tmp190 + tmp192;
                                    EM_S[INDEX4(0,1,0,5,numEq,numComp,8)]+= tmp184 + tmp341 + tmp345;
                                    EM_S[INDEX4(0,1,0,6,numEq,numComp,8)]+= tmp274 + tmp447 + tmp449;
                                    EM_S[INDEX4(0,1,0,7,numEq,numComp,8)]+= tmp107 + tmp109 + tmp112;
                                    EM_S[INDEX4(0,1,1,0,numEq,numComp,8)]+= tmp210 + tmp440 + tmp441 + w20*(A_0011[7] - A_0110[6]) + w21*(A_0011[0] - A_0110[1]);
                                    EM_S[INDEX4(0,1,1,1,numEq,numComp,8)]+= tmp351 + w20*(A_0011[6] + A_0110[6]) + w21*(A_0011[1] + A_0110[1]) + tmp76 + tmp79;
                                    EM_S[INDEX4(0,1,1,2,numEq,numComp,8)]+= w20*(-A_0011[7] - A_0110[4]) + w21*(-A_0011[0] - A_0110[3]) + tmp198 + tmp199 + tmp200;
                                    EM_S[INDEX4(0,1,1,3,numEq,numComp,8)]+= w20*(-A_0011[6] + A_0110[4]) + tmp44 + w21*(-A_0011[1] + A_0110[3]) + tmp49 + tmp53;
                                    EM_S[INDEX4(0,1,1,4,numEq,numComp,8)]+= tmp274 + tmp394 + tmp395;
                                    EM_S[INDEX4(0,1,1,5,numEq,numComp,8)]+= tmp112 + tmp143 + tmp147;
                                    EM_S[INDEX4(0,1,1,6,numEq,numComp,8)]+= tmp136 + tmp138 + tmp15;
                                    EM_S[INDEX4(0,1,1,7,numEq,numComp,8)]+= tmp184 + tmp330 + tmp334;
                                    EM_S[INDEX4(0,1,2,0,numEq,numComp,8)]+= w20*(-A_0011[5] + A_0110[7]) + w21*(-A_0011[2] + A_0110[0]) + tmp245 + tmp44 + tmp49;
                                    EM_S[INDEX4(0,1,2,1,numEq,numComp,8)]+= tmp120 + tmp198 + tmp200 + w20*(-A_0011[4] - A_0110[7]) + w21*(-A_0011[3] - A_0110[0]);
                                    EM_S[INDEX4(0,1,2,2,numEq,numComp,8)]+= tmp290 + w20*(A_0011[5] + A_0110[5]) + w21*(A_0011[2] + A_0110[2]) + tmp76 + tmp79;
                                    EM_S[INDEX4(0,1,2,3,numEq,numComp,8)]+= w20*(A_0011[4] - A_0110[5]) + w21*(A_0011[3] - A_0110[2]) + tmp440 + tmp441 + tmp99;
                                    EM_S[INDEX4(0,1,2,4,numEq,numComp,8)]+= tmp184 + tmp423 + tmp427;
                                    EM_S[INDEX4(0,1,2,5,numEq,numComp,8)]+= tmp15 + tmp162 + tmp163;
                                    EM_S[INDEX4(0,1,2,6,numEq,numComp,8)]+= tmp112 + tmp128 + tmp130;
                                    EM_S[INDEX4(0,1,2,7,numEq,numComp,8)]+= tmp270 + tmp272 + tmp274;
                                    EM_S[INDEX4(0,1,3,0,numEq,numComp,8)]+= tmp34 + w20*(A_0011[5] + A_0110[6]) + tmp76 + w21*(A_0011[2] + A_0110[1]) + tmp79;
                                    EM_S[INDEX4(0,1,3,1,numEq,numComp,8)]+= tmp440 + tmp441 + tmp481 + w20*(A_0011[4] - A_0110[6]) + w21*(A_0011[3] - A_0110[1]);
                                    EM_S[INDEX4(0,1,3,2,numEq,numComp,8)]+= w20*(-A_0011[5] + A_0110[4]) + w21*(-A_0011[2] + A_0110[3]) + tmp376 + tmp44 + tmp49;
                                    EM_S[INDEX4(0,1,3,3,numEq,numComp,8)]+= tmp198 + tmp200 + tmp216 + w20*(-A_0011[4] - A_0110[4]) + w21*(-A_0011[3] - A_0110[3]);
                                    EM_S[INDEX4(0,1,3,4,numEq,numComp,8)]+= tmp112 + tmp471 + tmp473;
                                    EM_S[INDEX4(0,1,3,5,numEq,numComp,8)]+= tmp274 + tmp399 + tmp401;
                                    EM_S[INDEX4(0,1,3,6,numEq,numComp,8)]+= tmp177 + tmp181 + tmp184;
                                    EM_S[INDEX4(0,1,3,7,numEq,numComp,8)]+= tmp11 + tmp15 + tmp7;
                                    EM_S[INDEX4(0,1,4,0,numEq,numComp,8)]+= tmp15 + tmp190 + tmp192;
                                    EM_S[INDEX4(0,1,4,1,numEq,numComp,8)]+= tmp184 + tmp341 + tmp345;
                                    EM_S[INDEX4(0,1,4,2,numEq,numComp,8)]+= tmp274 + tmp447 + tmp449;
                                    EM_S[INDEX4(0,1,4,3,numEq,numComp,8)]+= tmp107 + tmp109 + tmp112;
                                    EM_S[INDEX4(0,1,4,4,numEq,numComp,8)]+= tmp118 + tmp121 + tmp216 + w20*(-A_0011[3] - A_0110[3]) + w21*(-A_0011[4] - A_0110[4]);
                                    EM_S[INDEX4(0,1,4,5,numEq,numComp,8)]+= tmp376 + w20*(-A_0011[2] + A_0110[3]) + w21*(-A_0011[5] + A_0110[4]) + tmp60 + tmp65;
                                    EM_S[INDEX4(0,1,4,6,numEq,numComp,8)]+= w20*(A_0011[3] - A_0110[1]) + w21*(A_0011[4] - A_0110[6]) + tmp481 + tmp90 + tmp94;
                                    EM_S[INDEX4(0,1,4,7,numEq,numComp,8)]+= w20*(A_0011[2] + A_0110[1]) + tmp26 + tmp29 + w21*(A_0011[5] + A_0110[6]) + tmp34;
                                    EM_S[INDEX4(0,1,5,0,numEq,numComp,8)]+= tmp274 + tmp394 + tmp395;
                                    EM_S[INDEX4(0,1,5,1,numEq,numComp,8)]+= tmp112 + tmp143 + tmp147;
                                    EM_S[INDEX4(0,1,5,2,numEq,numComp,8)]+= tmp136 + tmp138 + tmp15;
                                    EM_S[INDEX4(0,1,5,3,numEq,numComp,8)]+= tmp184 + tmp330 + tmp334;
                                    EM_S[INDEX4(0,1,5,4,numEq,numComp,8)]+= w20*(A_0011[3] - A_0110[2]) + tmp90 + w21*(A_0011[4] - A_0110[5]) + tmp94 + tmp99;
                                    EM_S[INDEX4(0,1,5,5,numEq,numComp,8)]+= tmp26 + tmp29 + tmp290 + w20*(A_0011[2] + A_0110[2]) + w21*(A_0011[5] + A_0110[5]);
                                    EM_S[INDEX4(0,1,5,6,numEq,numComp,8)]+= w21*(-A_0011[4] - A_0110[7]) + w20*(-A_0011[3] - A_0110[0]) + tmp118 + tmp120 + tmp121;
                                    EM_S[INDEX4(0,1,5,7,numEq,numComp,8)]+= tmp245 + w21*(-A_0011[5] + A_0110[7]) + w20*(-A_0011[2] + A_0110[0]) + tmp60 + tmp65;
                                    EM_S[INDEX4(0,1,6,0,numEq,numComp,8)]+= tmp184 + tmp423 + tmp427;
                                    EM_S[INDEX4(0,1,6,1,numEq,numComp,8)]+= tmp15 + tmp162 + tmp163;
                                    EM_S[INDEX4(0,1,6,2,numEq,numComp,8)]+= tmp112 + tmp128 + tmp130;
                                    EM_S[INDEX4(0,1,6,3,numEq,numComp,8)]+= tmp270 + tmp272 + tmp274;
                                    EM_S[INDEX4(0,1,6,4,numEq,numComp,8)]+= tmp53 + w20*(-A_0011[1] + A_0110[3]) + tmp60 + tmp65 + w21*(-A_0011[6] + A_0110[4]);
                                    EM_S[INDEX4(0,1,6,5,numEq,numComp,8)]+= tmp118 + tmp121 + tmp199 + w21*(-A_0011[7] - A_0110[4]) + w20*(-A_0011[0] - A_0110[3]);
                                    EM_S[INDEX4(0,1,6,6,numEq,numComp,8)]+= tmp26 + tmp29 + tmp351 + w20*(A_0011[1] + A_0110[1]) + w21*(A_0011[6] + A_0110[6]);
                                    EM_S[INDEX4(0,1,6,7,numEq,numComp,8)]+= w20*(A_0011[0] - A_0110[1]) + w21*(A_0011[7] - A_0110[6]) + tmp210 + tmp90 + tmp94;
                                    EM_S[INDEX4(0,1,7,0,numEq,numComp,8)]+= tmp112 + tmp471 + tmp473;
                                    EM_S[INDEX4(0,1,7,1,numEq,numComp,8)]+= tmp274 + tmp399 + tmp401;
                                    EM_S[INDEX4(0,1,7,2,numEq,numComp,8)]+= tmp177 + tmp181 + tmp184;
                                    EM_S[INDEX4(0,1,7,3,numEq,numComp,8)]+= tmp11 + tmp15 + tmp7;
                                    EM_S[INDEX4(0,1,7,4,numEq,numComp,8)]+= tmp171 + tmp26 + tmp29 + w20*(A_0011[1] + A_0110[2]) + w21*(A_0011[6] + A_0110[5]);
                                    EM_S[INDEX4(0,1,7,5,numEq,numComp,8)]+= w21*(A_0011[7] - A_0110[5]) + w20*(A_0011[0] - A_0110[2]) + tmp436 + tmp90 + tmp94;
                                    EM_S[INDEX4(0,1,7,6,numEq,numComp,8)]+= w20*(-A_0011[1] + A_0110[0]) + w21*(-A_0011[6] + A_0110[7]) + tmp250 + tmp60 + tmp65;
                                    EM_S[INDEX4(0,1,7,7,numEq,numComp,8)]+= tmp118 + tmp121 + tmp303 + w20*(-A_0011[0] - A_0110[0]) + w21*(-A_0011[7] - A_0110[7]);
                                }
                                {
                                    const double tmp2 = w11*(-A_0022[2] - A_0022[5] + A_0220[1] + A_0220[6]);
                                    const double tmp9 = w15*(-A_0022[3] - A_0022[6] + A_0220[2] + A_0220[7]);
                                    const double tmp14 = w5*(-A_0022[1] - A_0022[4] + A_0220[0] + A_0220[5]);
                                    const double tmp21 = w11*(-A_0022[1] - A_0022[3] - A_0022[4] - A_0022[6] + A_0220[0] + A_0220[2] + A_0220[5] + A_0220[7]);
                                    const double tmp27 = w15*(-A_0022[5] - A_0022[7] + A_0220[4] + A_0220[6]);
                                    const double tmp33 = w5*(-A_0022[0] - A_0022[2] + A_0220[1] + A_0220[3]);
                                    const double tmp38 = w11*(-A_0022[0] - A_0022[2] - A_0022[5] - A_0022[7] - A_0220[0] - A_0220[2] - A_0220[5] - A_0220[7]);
                                    const double tmp47 = w15*(-A_0022[1] - A_0022[3] - A_0220[1] - A_0220[3]);
                                    const double tmp50 = w5*(-A_0022[4] - A_0022[6] - A_0220[4] - A_0220[6]);
                                    const double tmp63 = w15*(-A_0022[4] - A_0022[6] - A_0220[4] - A_0220[6]);
                                    const double tmp69 = w5*(-A_0022[1] - A_0022[3] - A_0220[1] - A_0220[3]);
                                    const double tmp77 = w15*(-A_0022[0] - A_0022[2] + A_0220[1] + A_0220[3]);
                                    const double tmp82 = w5*(-A_0022[5] - A_0022[7] + A_0220[4] + A_0220[6]);
                                    const double tmp85 = w11*(A_0022[1] + A_0022[6] - A_0220[0] - A_0220[7]);
                                    const double tmp92 = w15*(A_0022[0] + A_0022[5] - A_0220[1] - A_0220[4]);
                                    const double tmp98 = w5*(A_0022[2] + A_0022[7] - A_0220[3] - A_0220[6]);
                                    const double tmp108 = w15*(-A_0022[1] - A_0022[3] - A_0220[4] - A_0220[6]);
                                    const double tmp111 = w5*(-A_0022[4] - A_0022[6] - A_0220[1] - A_0220[3]);
                                    const double tmp114 = w11*(A_0022[0] + A_0022[2] + A_0022[5] + A_0022[7] - A_0220[1] - A_0220[3] - A_0220[4] - A_0220[6]);
                                    const double tmp117 = w15*(A_0022[4] + A_0022[6] - A_0220[5] - A_0220[7]);
                                    const double tmp119 = w5*(A_0022[1] + A_0022[3] - A_0220[0] - A_0220[2]);
                                    const double tmp124 = w11*(-A_0022[0] - A_0022[7] + A_0220[3] + A_0220[4]);
                                    const double tmp135 = w11*(A_0022[1] + A_0022[3] + A_0022[4] + A_0022[6] + A_0220[1] + A_0220[3] + A_0220[4] + A_0220[6]);
                                    const double tmp137 = w15*(A_0022[0] + A_0022[2] + A_0220[5] + A_0220[7]);
                                    const double tmp139 = w5*(A_0022[5] + A_0022[7] + A_0220[0] + A_0220[2]);
                                    const double tmp145 = w15*(-A_0022[1] - A_0022[4] + A_0220[0] + A_0220[5]);
                                    const double tmp148 = w5*(-A_0022[3] - A_0022[6] + A_0220[2] + A_0220[7]);
                                    const double tmp153 = w11*(A_0022[1] + A_0022[6] - A_0220[2] - A_0220[5]);
                                    const double tmp156 = w15*(A_0022[2] + A_0022[7] - A_0220[3] - A_0220[6]);
                                    const double tmp157 = w5*(A_0022[0] + A_0022[5] - A_0220[1] - A_0220[4]);
                                    const double tmp167 = w15*(A_0022[1] + A_0022[3] - A_0220[0] - A_0220[2]);
                                    const double tmp170 = w5*(A_0022[4] + A_0022[6] - A_0220[5] - A_0220[7]);
                                    const double tmp174 = w11*(-A_0022[3] - A_0022[4] - A_0220[1] - A_0220[6]);
                                    const double tmp179 = w15*(-A_0022[2] - A_0022[7] - A_0220[2] - A_0220[7]);
                                    const double tmp183 = w5*(-A_0022[0] - A_0022[5] - A_0220[0] - A_0220[5]);
                                    const double tmp202 = w11*(-A_0022[2] - A_0022[5] + A_0220[3] + A_0220[4]);
                                    const double tmp220 = w11*(-A_0022[1] - A_0022[6] - A_0220[1] - A_0220[6]);
                                    const double tmp240 = w15*(A_0022[0] + A_0022[2] + A_0220[0] + A_0220[2]);
                                    const double tmp242 = w5*(A_0022[5] + A_0022[7] + A_0220[5] + A_0220[7]);
                                    const double tmp247 = w11*(A_0022[3] + A_0022[4] - A_0220[2] - A_0220[5]);
                                    const double tmp260 = w15*(-A_0022[0] - A_0022[5] - A_0220[0] - A_0220[5]);
                                    const double tmp267 = w5*(-A_0022[2] - A_0022[7] - A_0220[2] - A_0220[7]);
                                    const double tmp269 = w11*(A_0022[2] + A_0022[5] + A_0220[0] + A_0220[7]);
                                    const double tmp271 = w15*(A_0022[3] + A_0022[6] + A_0220[3] + A_0220[6]);
                                    const double tmp273 = w5*(A_0022[1] + A_0022[4] + A_0220[1] + A_0220[4]);
                                    const double tmp278 = w11*(A_0022[3] + A_0022[4] - A_0220[0] - A_0220[7]);
                                    const double tmp283 = w11*(A_0022[0] + A_0022[7] + A_0220[2] + A_0220[5]);
                                    const double tmp293 = w11*(A_0022[0] + A_0022[7] + A_0220[0] + A_0220[7]);
                                    const double tmp307 = w11*(A_0022[2] + A_0022[5] + A_0220[2] + A_0220[5]);
                                    const double tmp324 = w15*(A_0022[5] + A_0022[7] + A_0220[5] + A_0220[7]);
                                    const double tmp326 = w5*(A_0022[0] + A_0022[2] + A_0220[0] + A_0220[2]);
                                    const double tmp333 = w15*(-A_0022[5] - A_0022[7] + A_0220[1] + A_0220[3]);
                                    const double tmp336 = w5*(-A_0022[0] - A_0022[2] + A_0220[4] + A_0220[6]);
                                    const double tmp343 = w15*(A_0022[1] + A_0022[4] + A_0220[1] + A_0220[4]);
                                    const double tmp347 = w5*(A_0022[3] + A_0022[6] + A_0220[3] + A_0220[6]);
                                    const double tmp354 = w11*(-A_0022[3] - A_0022[4] - A_0220[3] - A_0220[4]);
                                    const double tmp365 = w11*(-A_0022[1] - A_0022[6] - A_0220[3] - A_0220[4]);
                                    const double tmp369 = w11*(-A_0022[0] - A_0022[7] + A_0220[1] + A_0220[6]);
                                    const double tmp426 = w15*(A_0022[4] + A_0022[6] - A_0220[0] - A_0220[2]);
                                    const double tmp429 = w5*(A_0022[1] + A_0022[3] - A_0220[5] - A_0220[7]);
                                    const double tmp464 = w15*(A_0022[1] + A_0022[3] - A_0220[5] - A_0220[7]);
                                    const double tmp465 = w5*(A_0022[4] + A_0022[6] - A_0220[0] - A_0220[2]);
                                    const double tmp472 = w15*(-A_0022[4] - A_0022[6] - A_0220[1] - A_0220[3]);
                                    const double tmp475 = w5*(-A_0022[1] - A_0022[3] - A_0220[4] - A_0220[6]);
                                    const double tmp484 = w15*(A_0022[5] + A_0022[7] + A_0220[0] + A_0220[2]);
                                    const double tmp485 = w5*(A_0022[0] + A_0022[2] + A_0220[5] + A_0220[7]);
                                    const double tmp498 = w15*(-A_0022[0] - A_0022[2] + A_0220[4] + A_0220[6]);
                                    const double tmp499 = w5*(-A_0022[5] - A_0022[7] + A_0220[1] + A_0220[3]);
                                    EM_S[INDEX4(0,2,0,0,numEq,numComp,8)]+= tmp307 + tmp343 + tmp347 + w17*(A_0022[0] + A_0220[0]) + w2*(-A_0022[7] - A_0220[7]);
                                    EM_S[INDEX4(0,2,0,1,numEq,numComp,8)]+= tmp247 + w2*(-A_0022[6] + A_0220[7]) + w17*(A_0022[1] - A_0220[0]) + tmp92 + tmp98;
                                    EM_S[INDEX4(0,2,0,2,numEq,numComp,8)]+= tmp135 + tmp240 + tmp242;
                                    EM_S[INDEX4(0,2,0,3,numEq,numComp,8)]+= tmp114 + tmp167 + tmp170;
                                    EM_S[INDEX4(0,2,0,4,numEq,numComp,8)]+= tmp145 + tmp148 + tmp2 + w17*(-A_0022[0] + A_0220[4]) + w2*(A_0022[7] - A_0220[3]);
                                    EM_S[INDEX4(0,2,0,5,numEq,numComp,8)]+= tmp174 + tmp260 + tmp267 + w2*(A_0022[6] + A_0220[3]) + w17*(-A_0022[1] - A_0220[4]);
                                    EM_S[INDEX4(0,2,0,6,numEq,numComp,8)]+= tmp21 + tmp498 + tmp499;
                                    EM_S[INDEX4(0,2,0,7,numEq,numComp,8)]+= tmp108 + tmp111 + tmp38;
                                    EM_S[INDEX4(0,2,1,0,numEq,numComp,8)]+= tmp145 + tmp148 + tmp202 + w2*(A_0022[7] - A_0220[6]) + w17*(-A_0022[0] + A_0220[1]);
                                    EM_S[INDEX4(0,2,1,1,numEq,numComp,8)]+= tmp260 + tmp267 + w17*(-A_0022[1] - A_0220[1]) + w2*(A_0022[6] + A_0220[6]) + tmp354;
                                    EM_S[INDEX4(0,2,1,2,numEq,numComp,8)]+= tmp21 + tmp77 + tmp82;
                                    EM_S[INDEX4(0,2,1,3,numEq,numComp,8)]+= tmp38 + tmp47 + tmp50;
                                    EM_S[INDEX4(0,2,1,4,numEq,numComp,8)]+= tmp269 + tmp343 + tmp347 + w17*(A_0022[0] + A_0220[5]) + w2*(-A_0022[7] - A_0220[2]);
                                    EM_S[INDEX4(0,2,1,5,numEq,numComp,8)]+= tmp278 + w17*(A_0022[1] - A_0220[5]) + w2*(-A_0022[6] + A_0220[2]) + tmp92 + tmp98;
                                    EM_S[INDEX4(0,2,1,6,numEq,numComp,8)]+= tmp135 + tmp137 + tmp139;
                                    EM_S[INDEX4(0,2,1,7,numEq,numComp,8)]+= tmp114 + tmp464 + tmp465;
                                    EM_S[INDEX4(0,2,2,0,numEq,numComp,8)]+= tmp135 + tmp240 + tmp242;
                                    EM_S[INDEX4(0,2,2,1,numEq,numComp,8)]+= tmp114 + tmp167 + tmp170;
                                    EM_S[INDEX4(0,2,2,2,numEq,numComp,8)]+= tmp271 + tmp273 + w17*(A_0022[2] + A_0220[2]) + w2*(-A_0022[5] - A_0220[5]) + tmp293;
                                    EM_S[INDEX4(0,2,2,3,numEq,numComp,8)]+= tmp156 + tmp157 + w2*(-A_0022[4] + A_0220[5]) + w17*(A_0022[3] - A_0220[2]) + tmp85;
                                    EM_S[INDEX4(0,2,2,4,numEq,numComp,8)]+= tmp21 + tmp498 + tmp499;
                                    EM_S[INDEX4(0,2,2,5,numEq,numComp,8)]+= tmp108 + tmp111 + tmp38;
                                    EM_S[INDEX4(0,2,2,6,numEq,numComp,8)]+= tmp124 + w17*(-A_0022[2] + A_0220[6]) + w2*(A_0022[5] - A_0220[1]) + tmp14 + tmp9;
                                    EM_S[INDEX4(0,2,2,7,numEq,numComp,8)]+= tmp179 + tmp183 + w2*(A_0022[4] + A_0220[1]) + tmp365 + w17*(-A_0022[3] - A_0220[6]);
                                    EM_S[INDEX4(0,2,3,0,numEq,numComp,8)]+= tmp21 + tmp77 + tmp82;
                                    EM_S[INDEX4(0,2,3,1,numEq,numComp,8)]+= tmp38 + tmp47 + tmp50;
                                    EM_S[INDEX4(0,2,3,2,numEq,numComp,8)]+= tmp14 + w2*(A_0022[5] - A_0220[4]) + tmp369 + w17*(-A_0022[2] + A_0220[3]) + tmp9;
                                    EM_S[INDEX4(0,2,3,3,numEq,numComp,8)]+= tmp179 + tmp183 + w17*(-A_0022[3] - A_0220[3]) + w2*(A_0022[4] + A_0220[4]) + tmp220;
                                    EM_S[INDEX4(0,2,3,4,numEq,numComp,8)]+= tmp135 + tmp137 + tmp139;
                                    EM_S[INDEX4(0,2,3,5,numEq,numComp,8)]+= tmp114 + tmp464 + tmp465;
                                    EM_S[INDEX4(0,2,3,6,numEq,numComp,8)]+= tmp271 + tmp273 + tmp283 + w17*(A_0022[2] + A_0220[7]) + w2*(-A_0022[5] - A_0220[0]);
                                    EM_S[INDEX4(0,2,3,7,numEq,numComp,8)]+= tmp153 + tmp156 + tmp157 + w17*(A_0022[3] - A_0220[7]) + w2*(-A_0022[4] + A_0220[0]);
                                    EM_S[INDEX4(0,2,4,0,numEq,numComp,8)]+= tmp153 + w17*(A_0022[4] - A_0220[0]) + w2*(-A_0022[3] + A_0220[7]) + tmp92 + tmp98;
                                    EM_S[INDEX4(0,2,4,1,numEq,numComp,8)]+= tmp283 + tmp343 + tmp347 + w17*(A_0022[5] + A_0220[0]) + w2*(-A_0022[2] - A_0220[7]);
                                    EM_S[INDEX4(0,2,4,2,numEq,numComp,8)]+= tmp114 + tmp426 + tmp429;
                                    EM_S[INDEX4(0,2,4,3,numEq,numComp,8)]+= tmp135 + tmp484 + tmp485;
                                    EM_S[INDEX4(0,2,4,4,numEq,numComp,8)]+= tmp220 + w17*(-A_0022[4] - A_0220[4]) + w2*(A_0022[3] + A_0220[3]) + tmp260 + tmp267;
                                    EM_S[INDEX4(0,2,4,5,numEq,numComp,8)]+= tmp145 + tmp148 + tmp369 + w17*(-A_0022[5] + A_0220[4]) + w2*(A_0022[2] - A_0220[3]);
                                    EM_S[INDEX4(0,2,4,6,numEq,numComp,8)]+= tmp38 + tmp63 + tmp69;
                                    EM_S[INDEX4(0,2,4,7,numEq,numComp,8)]+= tmp21 + tmp27 + tmp33;
                                    EM_S[INDEX4(0,2,5,0,numEq,numComp,8)]+= tmp260 + tmp267 + tmp365 + w2*(A_0022[3] + A_0220[6]) + w17*(-A_0022[4] - A_0220[1]);
                                    EM_S[INDEX4(0,2,5,1,numEq,numComp,8)]+= tmp124 + tmp145 + tmp148 + w17*(-A_0022[5] + A_0220[1]) + w2*(A_0022[2] - A_0220[6]);
                                    EM_S[INDEX4(0,2,5,2,numEq,numComp,8)]+= tmp38 + tmp472 + tmp475;
                                    EM_S[INDEX4(0,2,5,3,numEq,numComp,8)]+= tmp21 + tmp333 + tmp336;
                                    EM_S[INDEX4(0,2,5,4,numEq,numComp,8)]+= w17*(A_0022[4] - A_0220[5]) + w2*(-A_0022[3] + A_0220[2]) + tmp85 + tmp92 + tmp98;
                                    EM_S[INDEX4(0,2,5,5,numEq,numComp,8)]+= tmp293 + tmp343 + tmp347 + w17*(A_0022[5] + A_0220[5]) + w2*(-A_0022[2] - A_0220[2]);
                                    EM_S[INDEX4(0,2,5,6,numEq,numComp,8)]+= tmp114 + tmp117 + tmp119;
                                    EM_S[INDEX4(0,2,5,7,numEq,numComp,8)]+= tmp135 + tmp324 + tmp326;
                                    EM_S[INDEX4(0,2,6,0,numEq,numComp,8)]+= tmp114 + tmp426 + tmp429;
                                    EM_S[INDEX4(0,2,6,1,numEq,numComp,8)]+= tmp135 + tmp484 + tmp485;
                                    EM_S[INDEX4(0,2,6,2,numEq,numComp,8)]+= tmp156 + tmp157 + tmp278 + w17*(A_0022[6] - A_0220[2]) + w2*(-A_0022[1] + A_0220[5]);
                                    EM_S[INDEX4(0,2,6,3,numEq,numComp,8)]+= tmp269 + tmp271 + tmp273 + w17*(A_0022[7] + A_0220[2]) + w2*(-A_0022[0] - A_0220[5]);
                                    EM_S[INDEX4(0,2,6,4,numEq,numComp,8)]+= tmp38 + tmp63 + tmp69;
                                    EM_S[INDEX4(0,2,6,5,numEq,numComp,8)]+= tmp21 + tmp27 + tmp33;
                                    EM_S[INDEX4(0,2,6,6,numEq,numComp,8)]+= tmp179 + tmp183 + tmp354 + w17*(-A_0022[6] - A_0220[6]) + w2*(A_0022[1] + A_0220[1]);
                                    EM_S[INDEX4(0,2,6,7,numEq,numComp,8)]+= tmp14 + tmp202 + w17*(-A_0022[7] + A_0220[6]) + w2*(A_0022[0] - A_0220[1]) + tmp9;
                                    EM_S[INDEX4(0,2,7,0,numEq,numComp,8)]+= tmp38 + tmp472 + tmp475;
                                    EM_S[INDEX4(0,2,7,1,numEq,numComp,8)]+= tmp21 + tmp333 + tmp336;
                                    EM_S[INDEX4(0,2,7,2,numEq,numComp,8)]+= w2*(A_0022[1] + A_0220[4]) + tmp174 + tmp179 + tmp183 + w17*(-A_0022[6] - A_0220[3]);
                                    EM_S[INDEX4(0,2,7,3,numEq,numComp,8)]+= tmp14 + w17*(-A_0022[7] + A_0220[3]) + w2*(A_0022[0] - A_0220[4]) + tmp2 + tmp9;
                                    EM_S[INDEX4(0,2,7,4,numEq,numComp,8)]+= tmp114 + tmp117 + tmp119;
                                    EM_S[INDEX4(0,2,7,5,numEq,numComp,8)]+= tmp135 + tmp324 + tmp326;
                                    EM_S[INDEX4(0,2,7,6,numEq,numComp,8)]+= tmp156 + tmp157 + tmp247 + w17*(A_0022[6] - A_0220[7]) + w2*(-A_0022[1] + A_0220[0]);
                                    EM_S[INDEX4(0,2,7,7,numEq,numComp,8)]+= tmp271 + tmp273 + w17*(A_0022[7] + A_0220[7]) + w2*(-A_0022[0] - A_0220[0]) + tmp307;
                                }
                                {
                                    const double tmp7 = w1*(A_1001[0] + A_1001[4] + A_1100[0] + A_1100[4]);
                                    const double tmp11 = w16*(A_1001[3] + A_1001[7] + A_1100[3] + A_1100[7]);
                                    const double tmp15 = w8*(A_1001[1] + A_1001[2] + A_1001[5] + A_1001[6] + A_1100[1] + A_1100[2] + A_1100[5] + A_1100[6]);
                                    const double tmp26 = w1*(-A_1001[0] - A_1001[3] - A_1100[0] - A_1100[3]);
                                    const double tmp29 = w16*(-A_1001[4] - A_1001[7] - A_1100[4] - A_1100[7]);
                                    const double tmp34 = w8*(-A_1001[1] - A_1001[6] - A_1100[2] - A_1100[5]);
                                    const double tmp44 = w1*(A_1001[4] + A_1001[7] - A_1100[5] - A_1100[6]);
                                    const double tmp49 = w16*(A_1001[0] + A_1001[3] - A_1100[1] - A_1100[2]);
                                    const double tmp53 = w8*(A_1001[2] + A_1001[5] - A_1100[0] - A_1100[7]);
                                    const double tmp60 = w1*(A_1001[0] + A_1001[3] - A_1100[1] - A_1100[2]);
                                    const double tmp65 = w16*(A_1001[4] + A_1001[7] - A_1100[5] - A_1100[6]);
                                    const double tmp76 = w1*(-A_1001[4] - A_1001[7] - A_1100[4] - A_1100[7]);
                                    const double tmp79 = w16*(-A_1001[0] - A_1001[3] - A_1100[0] - A_1100[3]);
                                    const double tmp90 = w1*(-A_1001[1] - A_1001[2] + A_1100[0] + A_1100[3]);
                                    const double tmp94 = w16*(-A_1001[5] - A_1001[6] + A_1100[4] + A_1100[7]);
                                    const double tmp99 = w8*(-A_1001[0] - A_1001[7] + A_1100[1] + A_1100[6]);
                                    const double tmp107 = w1*(-A_1001[2] - A_1001[6] - A_1100[1] - A_1100[5]);
                                    const double tmp109 = w16*(-A_1001[1] - A_1001[5] - A_1100[2] - A_1100[6]);
                                    const double tmp112 = w8*(-A_1001[0] - A_1001[3] - A_1001[4] - A_1001[7] - A_1100[0] - A_1100[3] - A_1100[4] - A_1100[7]);
                                    const double tmp118 = w16*(A_1001[5] + A_1001[6] + A_1100[5] + A_1100[6]);
                                    const double tmp120 = w8*(A_1001[0] + A_1001[7] + A_1100[3] + A_1100[4]);
                                    const double tmp121 = w1*(A_1001[1] + A_1001[2] + A_1100[1] + A_1100[2]);
                                    const double tmp128 = w1*(-A_1001[1] - A_1001[5] - A_1100[1] - A_1100[5]);
                                    const double tmp130 = w16*(-A_1001[2] - A_1001[6] - A_1100[2] - A_1100[6]);
                                    const double tmp136 = w1*(A_1001[3] + A_1001[7] + A_1100[0] + A_1100[4]);
                                    const double tmp138 = w16*(A_1001[0] + A_1001[4] + A_1100[3] + A_1100[7]);
                                    const double tmp143 = w1*(-A_1001[2] - A_1001[6] - A_1100[2] - A_1100[6]);
                                    const double tmp147 = w16*(-A_1001[1] - A_1001[5] - A_1100[1] - A_1100[5]);
                                    const double tmp162 = w1*(A_1001[0] + A_1001[4] + A_1100[3] + A_1100[7]);
                                    const double tmp163 = w16*(A_1001[3] + A_1001[7] + A_1100[0] + A_1100[4]);
                                    const double tmp171 = w8*(-A_1001[2] - A_1001[5] - A_1100[1] - A_1100[6]);
                                    const double tmp177 = w1*(A_1001[1] + A_1001[5] - A_1100[0] - A_1100[4]);
                                    const double tmp181 = w16*(A_1001[2] + A_1001[6] - A_1100[3] - A_1100[7]);
                                    const double tmp184 = w8*(A_1001[0] + A_1001[3] + A_1001[4] + A_1001[7] - A_1100[1] - A_1100[2] - A_1100[5] - A_1100[6]);
                                    const double tmp190 = w1*(A_1001[3] + A_1001[7] + A_1100[3] + A_1100[7]);
                                    const double tmp192 = w16*(A_1001[0] + A_1001[4] + A_1100[0] + A_1100[4]);
                                    const double tmp198 = w16*(A_1001[1] + A_1001[2] + A_1100[1] + A_1100[2]);
                                    const double tmp199 = w8*(A_1001[3] + A_1001[4] + A_1100[0] + A_1100[7]);
                                    const double tmp200 = w1*(A_1001[5] + A_1001[6] + A_1100[5] + A_1100[6]);
                                    const double tmp210 = w8*(-A_1001[3] - A_1001[4] + A_1100[2] + A_1100[5]);
                                    const double tmp216 = w8*(A_1001[0] + A_1001[7] + A_1100[0] + A_1100[7]);
                                    const double tmp245 = w8*(A_1001[1] + A_1001[6] - A_1100[3] - A_1100[4]);
                                    const double tmp250 = w8*(A_1001[2] + A_1001[5] - A_1100[3] - A_1100[4]);
                                    const double tmp270 = w1*(-A_1001[0] - A_1001[4] + A_1100[1] + A_1100[5]);
                                    const double tmp272 = w16*(-A_1001[3] - A_1001[7] + A_1100[2] + A_1100[6]);
                                    const double tmp274 = w8*(-A_1001[1] - A_1001[2] - A_1001[5] - A_1001[6] + A_1100[0] + A_1100[3] + A_1100[4] + A_1100[7]);
                                    const double tmp290 = w8*(-A_1001[1] - A_1001[6] - A_1100[1] - A_1100[6]);
                                    const double tmp303 = w8*(A_1001[3] + A_1001[4] + A_1100[3] + A_1100[4]);
                                    const double tmp330 = w1*(A_1001[2] + A_1001[6] - A_1100[0] - A_1100[4]);
                                    const double tmp334 = w16*(A_1001[1] + A_1001[5] - A_1100[3] - A_1100[7]);
                                    const double tmp341 = w1*(A_1001[2] + A_1001[6] - A_1100[3] - A_1100[7]);
                                    const double tmp345 = w16*(A_1001[1] + A_1001[5] - A_1100[0] - A_1100[4]);
                                    const double tmp351 = w8*(-A_1001[2] - A_1001[5] - A_1100[2] - A_1100[5]);
                                    const double tmp376 = w8*(A_1001[1] + A_1001[6] - A_1100[0] - A_1100[7]);
                                    const double tmp394 = w1*(-A_1001[3] - A_1001[7] + A_1100[2] + A_1100[6]);
                                    const double tmp395 = w16*(-A_1001[0] - A_1001[4] + A_1100[1] + A_1100[5]);
                                    const double tmp399 = w1*(-A_1001[0] - A_1001[4] + A_1100[2] + A_1100[6]);
                                    const double tmp401 = w16*(-A_1001[3] - A_1001[7] + A_1100[1] + A_1100[5]);
                                    const double tmp423 = w1*(A_1001[1] + A_1001[5] - A_1100[3] - A_1100[7]);
                                    const double tmp427 = w16*(A_1001[2] + A_1001[6] - A_1100[0] - A_1100[4]);
                                    const double tmp436 = w8*(-A_1001[3] - A_1001[4] + A_1100[1] + A_1100[6]);
                                    const double tmp440 = w16*(-A_1001[1] - A_1001[2] + A_1100[0] + A_1100[3]);
                                    const double tmp441 = w1*(-A_1001[5] - A_1001[6] + A_1100[4] + A_1100[7]);
                                    const double tmp447 = w1*(-A_1001[3] - A_1001[7] + A_1100[1] + A_1100[5]);
                                    const double tmp449 = w16*(-A_1001[0] - A_1001[4] + A_1100[2] + A_1100[6]);
                                    const double tmp471 = w1*(-A_1001[1] - A_1001[5] - A_1100[2] - A_1100[6]);
                                    const double tmp473 = w16*(-A_1001[2] - A_1001[6] - A_1100[1] - A_1100[5]);
                                    const double tmp481 = w8*(-A_1001[0] - A_1001[7] + A_1100[2] + A_1100[5]);
                                    EM_S[INDEX4(1,0,0,0,numEq,numComp,8)]+= tmp198 + tmp200 + tmp303 + w20*(-A_1001[7] - A_1100[7]) + w21*(-A_1001[0] - A_1100[0]);
                                    EM_S[INDEX4(1,0,0,1,numEq,numComp,8)]+= tmp250 + tmp44 + w20*(-A_1001[6] + A_1100[7]) + w21*(-A_1001[1] + A_1100[0]) + tmp49;
                                    EM_S[INDEX4(1,0,0,2,numEq,numComp,8)]+= tmp436 + tmp440 + tmp441 + w20*(A_1001[7] - A_1100[5]) + w21*(A_1001[0] - A_1100[2]);
                                    EM_S[INDEX4(1,0,0,3,numEq,numComp,8)]+= w20*(A_1001[6] + A_1100[5]) + w21*(A_1001[1] + A_1100[2]) + tmp171 + tmp76 + tmp79;
                                    EM_S[INDEX4(1,0,0,4,numEq,numComp,8)]+= tmp15 + tmp190 + tmp192;
                                    EM_S[INDEX4(1,0,0,5,numEq,numComp,8)]+= tmp184 + tmp341 + tmp345;
                                    EM_S[INDEX4(1,0,0,6,numEq,numComp,8)]+= tmp274 + tmp447 + tmp449;
                                    EM_S[INDEX4(1,0,0,7,numEq,numComp,8)]+= tmp107 + tmp109 + tmp112;
                                    EM_S[INDEX4(1,0,1,0,numEq,numComp,8)]+= tmp210 + tmp440 + tmp441 + w20*(A_1001[7] - A_1100[6]) + w21*(A_1001[0] - A_1100[1]);
                                    EM_S[INDEX4(1,0,1,1,numEq,numComp,8)]+= tmp351 + w20*(A_1001[6] + A_1100[6]) + w21*(A_1001[1] + A_1100[1]) + tmp76 + tmp79;
                                    EM_S[INDEX4(1,0,1,2,numEq,numComp,8)]+= w20*(-A_1001[7] - A_1100[4]) + w21*(-A_1001[0] - A_1100[3]) + tmp198 + tmp199 + tmp200;
                                    EM_S[INDEX4(1,0,1,3,numEq,numComp,8)]+= w20*(-A_1001[6] + A_1100[4]) + tmp44 + w21*(-A_1001[1] + A_1100[3]) + tmp49 + tmp53;
                                    EM_S[INDEX4(1,0,1,4,numEq,numComp,8)]+= tmp274 + tmp394 + tmp395;
                                    EM_S[INDEX4(1,0,1,5,numEq,numComp,8)]+= tmp112 + tmp143 + tmp147;
                                    EM_S[INDEX4(1,0,1,6,numEq,numComp,8)]+= tmp136 + tmp138 + tmp15;
                                    EM_S[INDEX4(1,0,1,7,numEq,numComp,8)]+= tmp184 + tmp330 + tmp334;
                                    EM_S[INDEX4(1,0,2,0,numEq,numComp,8)]+= w20*(-A_1001[5] + A_1100[7]) + w21*(-A_1001[2] + A_1100[0]) + tmp245 + tmp44 + tmp49;
                                    EM_S[INDEX4(1,0,2,1,numEq,numComp,8)]+= tmp120 + tmp198 + tmp200 + w20*(-A_1001[4] - A_1100[7]) + w21*(-A_1001[3] - A_1100[0]);
                                    EM_S[INDEX4(1,0,2,2,numEq,numComp,8)]+= tmp290 + w20*(A_1001[5] + A_1100[5]) + w21*(A_1001[2] + A_1100[2]) + tmp76 + tmp79;
                                    EM_S[INDEX4(1,0,2,3,numEq,numComp,8)]+= w20*(A_1001[4] - A_1100[5]) + w21*(A_1001[3] - A_1100[2]) + tmp440 + tmp441 + tmp99;
                                    EM_S[INDEX4(1,0,2,4,numEq,numComp,8)]+= tmp184 + tmp423 + tmp427;
                                    EM_S[INDEX4(1,0,2,5,numEq,numComp,8)]+= tmp15 + tmp162 + tmp163;
                                    EM_S[INDEX4(1,0,2,6,numEq,numComp,8)]+= tmp112 + tmp128 + tmp130;
                                    EM_S[INDEX4(1,0,2,7,numEq,numComp,8)]+= tmp270 + tmp272 + tmp274;
                                    EM_S[INDEX4(1,0,3,0,numEq,numComp,8)]+= tmp34 + w20*(A_1001[5] + A_1100[6]) + tmp76 + w21*(A_1001[2] + A_1100[1]) + tmp79;
                                    EM_S[INDEX4(1,0,3,1,numEq,numComp,8)]+= tmp440 + tmp441 + tmp481 + w20*(A_1001[4] - A_1100[6]) + w21*(A_1001[3] - A_1100[1]);
                                    EM_S[INDEX4(1,0,3,2,numEq,numComp,8)]+= w20*(-A_1001[5] + A_1100[4]) + w21*(-A_1001[2] + A_1100[3]) + tmp376 + tmp44 + tmp49;
                                    EM_S[INDEX4(1,0,3,3,numEq,numComp,8)]+= tmp198 + tmp200 + tmp216 + w20*(-A_1001[4] - A_1100[4]) + w21*(-A_1001[3] - A_1100[3]);
                                    EM_S[INDEX4(1,0,3,4,numEq,numComp,8)]+= tmp112 + tmp471 + tmp473;
                                    EM_S[INDEX4(1,0,3,5,numEq,numComp,8)]+= tmp274 + tmp399 + tmp401;
                                    EM_S[INDEX4(1,0,3,6,numEq,numComp,8)]+= tmp177 + tmp181 + tmp184;
                                    EM_S[INDEX4(1,0,3,7,numEq,numComp,8)]+= tmp11 + tmp15 + tmp7;
                                    EM_S[INDEX4(1,0,4,0,numEq,numComp,8)]+= tmp15 + tmp190 + tmp192;
                                    EM_S[INDEX4(1,0,4,1,numEq,numComp,8)]+= tmp184 + tmp341 + tmp345;
                                    EM_S[INDEX4(1,0,4,2,numEq,numComp,8)]+= tmp274 + tmp447 + tmp449;
                                    EM_S[INDEX4(1,0,4,3,numEq,numComp,8)]+= tmp107 + tmp109 + tmp112;
                                    EM_S[INDEX4(1,0,4,4,numEq,numComp,8)]+= tmp118 + tmp121 + tmp216 + w20*(-A_1001[3] - A_1100[3]) + w21*(-A_1001[4] - A_1100[4]);
                                    EM_S[INDEX4(1,0,4,5,numEq,numComp,8)]+= tmp376 + w20*(-A_1001[2] + A_1100[3]) + w21*(-A_1001[5] + A_1100[4]) + tmp60 + tmp65;
                                    EM_S[INDEX4(1,0,4,6,numEq,numComp,8)]+= w20*(A_1001[3] - A_1100[1]) + w21*(A_1001[4] - A_1100[6]) + tmp481 + tmp90 + tmp94;
                                    EM_S[INDEX4(1,0,4,7,numEq,numComp,8)]+= w20*(A_1001[2] + A_1100[1]) + tmp26 + tmp29 + w21*(A_1001[5] + A_1100[6]) + tmp34;
                                    EM_S[INDEX4(1,0,5,0,numEq,numComp,8)]+= tmp274 + tmp394 + tmp395;
                                    EM_S[INDEX4(1,0,5,1,numEq,numComp,8)]+= tmp112 + tmp143 + tmp147;
                                    EM_S[INDEX4(1,0,5,2,numEq,numComp,8)]+= tmp136 + tmp138 + tmp15;
                                    EM_S[INDEX4(1,0,5,3,numEq,numComp,8)]+= tmp184 + tmp330 + tmp334;
                                    EM_S[INDEX4(1,0,5,4,numEq,numComp,8)]+= w20*(A_1001[3] - A_1100[2]) + tmp90 + w21*(A_1001[4] - A_1100[5]) + tmp94 + tmp99;
                                    EM_S[INDEX4(1,0,5,5,numEq,numComp,8)]+= tmp26 + tmp29 + tmp290 + w20*(A_1001[2] + A_1100[2]) + w21*(A_1001[5] + A_1100[5]);
                                    EM_S[INDEX4(1,0,5,6,numEq,numComp,8)]+= w21*(-A_1001[4] - A_1100[7]) + w20*(-A_1001[3] - A_1100[0]) + tmp118 + tmp120 + tmp121;
                                    EM_S[INDEX4(1,0,5,7,numEq,numComp,8)]+= tmp245 + w21*(-A_1001[5] + A_1100[7]) + w20*(-A_1001[2] + A_1100[0]) + tmp60 + tmp65;
                                    EM_S[INDEX4(1,0,6,0,numEq,numComp,8)]+= tmp184 + tmp423 + tmp427;
                                    EM_S[INDEX4(1,0,6,1,numEq,numComp,8)]+= tmp15 + tmp162 + tmp163;
                                    EM_S[INDEX4(1,0,6,2,numEq,numComp,8)]+= tmp112 + tmp128 + tmp130;
                                    EM_S[INDEX4(1,0,6,3,numEq,numComp,8)]+= tmp270 + tmp272 + tmp274;
                                    EM_S[INDEX4(1,0,6,4,numEq,numComp,8)]+= tmp53 + w20*(-A_1001[1] + A_1100[3]) + tmp60 + tmp65 + w21*(-A_1001[6] + A_1100[4]);
                                    EM_S[INDEX4(1,0,6,5,numEq,numComp,8)]+= tmp118 + tmp121 + tmp199 + w21*(-A_1001[7] - A_1100[4]) + w20*(-A_1001[0] - A_1100[3]);
                                    EM_S[INDEX4(1,0,6,6,numEq,numComp,8)]+= tmp26 + tmp29 + tmp351 + w20*(A_1001[1] + A_1100[1]) + w21*(A_1001[6] + A_1100[6]);
                                    EM_S[INDEX4(1,0,6,7,numEq,numComp,8)]+= w20*(A_1001[0] - A_1100[1]) + w21*(A_1001[7] - A_1100[6]) + tmp210 + tmp90 + tmp94;
                                    EM_S[INDEX4(1,0,7,0,numEq,numComp,8)]+= tmp112 + tmp471 + tmp473;
                                    EM_S[INDEX4(1,0,7,1,numEq,numComp,8)]+= tmp274 + tmp399 + tmp401;
                                    EM_S[INDEX4(1,0,7,2,numEq,numComp,8)]+= tmp177 + tmp181 + tmp184;
                                    EM_S[INDEX4(1,0,7,3,numEq,numComp,8)]+= tmp11 + tmp15 + tmp7;
                                    EM_S[INDEX4(1,0,7,4,numEq,numComp,8)]+= tmp171 + tmp26 + tmp29 + w20*(A_1001[1] + A_1100[2]) + w21*(A_1001[6] + A_1100[5]);
                                    EM_S[INDEX4(1,0,7,5,numEq,numComp,8)]+= w21*(A_1001[7] - A_1100[5]) + w20*(A_1001[0] - A_1100[2]) + tmp436 + tmp90 + tmp94;
                                    EM_S[INDEX4(1,0,7,6,numEq,numComp,8)]+= w20*(-A_1001[1] + A_1100[0]) + w21*(-A_1001[6] + A_1100[7]) + tmp250 + tmp60 + tmp65;
                                    EM_S[INDEX4(1,0,7,7,numEq,numComp,8)]+= tmp118 + tmp121 + tmp303 + w20*(-A_1001[0] - A_1100[0]) + w21*(-A_1001[7] - A_1100[7]);
                                }
                                {
                                    const double tmp1 = w13*(A_1212[1] + A_1212[2] + A_1212[5] + A_1212[6]);
                                    const double tmp3 = w14*(A_1010[2] + A_1010[3] + A_1010[6] + A_1010[7]);
                                    const double tmp4 = w7*(A_1212[0] + A_1212[4]);
                                    const double tmp6 = w3*(A_1111[0] + A_1111[2] + A_1111[4] + A_1111[6]);
                                    const double tmp10 = w0*(A_1010[0] + A_1010[1] + A_1010[4] + A_1010[5]);
                                    const double tmp12 = w9*(A_1111[1] + A_1111[3] + A_1111[5] + A_1111[7]);
                                    const double tmp17 = w19*(A_1212[3] + A_1212[7]);
                                    const double tmp20 = w13*(-A_1212[0] - A_1212[1] - A_1212[2] - A_1212[3] - A_1212[4] - A_1212[5] - A_1212[6] - A_1212[7]);
                                    const double tmp22 = w14*(-A_1010[4] - A_1010[5] - A_1010[6] - A_1010[7]);
                                    const double tmp25 = w3*(-A_1111[0] - A_1111[1] - A_1111[2] - A_1111[3]);
                                    const double tmp28 = w0*(-A_1010[0] - A_1010[1] - A_1010[2] - A_1010[3]);
                                    const double tmp30 = w9*(-A_1111[4] - A_1111[5] - A_1111[6] - A_1111[7]);
                                    const double tmp39 = w14*(A_1010[0] + A_1010[1] + A_1010[2] + A_1010[3]);
                                    const double tmp40 = w26*(A_1111[4] + A_1111[6]);
                                    const double tmp41 = w0*(A_1010[4] + A_1010[5] + A_1010[6] + A_1010[7]);
                                    const double tmp43 = w22*(A_1111[0] + A_1111[2] + A_1111[5] + A_1111[7]);
                                    const double tmp45 = w25*(A_1212[1] + A_1212[3] + A_1212[5] + A_1212[7]);
                                    const double tmp52 = w24*(A_1111[1] + A_1111[3]);
                                    const double tmp55 = w23*(A_1212[0] + A_1212[2] + A_1212[4] + A_1212[6]);
                                    const double tmp57 = w14*(A_1010[4] + A_1010[5] + A_1010[6] + A_1010[7]);
                                    const double tmp58 = w26*(A_1111[1] + A_1111[3]);
                                    const double tmp61 = w25*(A_1212[0] + A_1212[2] + A_1212[4] + A_1212[6]);
                                    const double tmp64 = w0*(A_1010[0] + A_1010[1] + A_1010[2] + A_1010[3]);
                                    const double tmp66 = w24*(A_1111[4] + A_1111[6]);
                                    const double tmp71 = w23*(A_1212[1] + A_1212[3] + A_1212[5] + A_1212[7]);
                                    const double tmp73 = w14*(-A_1010[0] - A_1010[1] - A_1010[2] - A_1010[3]);
                                    const double tmp74 = w0*(-A_1010[4] - A_1010[5] - A_1010[6] - A_1010[7]);
                                    const double tmp75 = w3*(-A_1111[4] - A_1111[5] - A_1111[6] - A_1111[7]);
                                    const double tmp80 = w9*(-A_1111[0] - A_1111[1] - A_1111[2] - A_1111[3]);
                                    const double tmp88 = w3*(A_1111[0] + A_1111[1] + A_1111[2] + A_1111[3]);
                                    const double tmp89 = w23*(A_1212[2] + A_1212[3] + A_1212[6] + A_1212[7]);
                                    const double tmp91 = w25*(A_1212[0] + A_1212[1] + A_1212[4] + A_1212[5]);
                                    const double tmp95 = w28*(A_1010[2] + A_1010[3]);
                                    const double tmp97 = w29*(A_1010[4] + A_1010[5]);
                                    const double tmp100 = w9*(A_1111[4] + A_1111[5] + A_1111[6] + A_1111[7]);
                                    const double tmp101 = w27*(A_1010[0] + A_1010[1] + A_1010[6] + A_1010[7]);
                                    const double tmp104 = w13*(A_1212[0] + A_1212[1] + A_1212[2] + A_1212[3] + A_1212[4] + A_1212[5] + A_1212[6] + A_1212[7]);
                                    const double tmp106 = w22*(A_1111[0] + A_1111[1] + A_1111[2] + A_1111[3] + A_1111[4] + A_1111[5] + A_1111[6] + A_1111[7]);
                                    const double tmp113 = w27*(A_1010[0] + A_1010[1] + A_1010[2] + A_1010[3] + A_1010[4] + A_1010[5] + A_1010[6] + A_1010[7]);
                                    const double tmp123 = w13*(A_1212[0] + A_1212[3] + A_1212[4] + A_1212[7]);
                                    const double tmp125 = w7*(A_1212[1] + A_1212[5]);
                                    const double tmp127 = w3*(A_1111[1] + A_1111[3] + A_1111[5] + A_1111[7]);
                                    const double tmp131 = w9*(A_1111[0] + A_1111[2] + A_1111[4] + A_1111[6]);
                                    const double tmp132 = w19*(A_1212[2] + A_1212[6]);
                                    const double tmp141 = w14*(A_1010[0] + A_1010[1] + A_1010[4] + A_1010[5]);
                                    const double tmp142 = w7*(A_1212[2] + A_1212[6]);
                                    const double tmp146 = w0*(A_1010[2] + A_1010[3] + A_1010[6] + A_1010[7]);
                                    const double tmp149 = w19*(A_1212[1] + A_1212[5]);
                                    const double tmp175 = w14*(-A_1010[2] - A_1010[3] - A_1010[6] - A_1010[7]);
                                    const double tmp176 = w22*(-A_1111[0] - A_1111[1] - A_1111[2] - A_1111[3] - A_1111[4] - A_1111[5] - A_1111[6] - A_1111[7]);
                                    const double tmp178 = w25*(-A_1212[2] - A_1212[3] - A_1212[6] - A_1212[7]);
                                    const double tmp180 = w0*(-A_1010[0] - A_1010[1] - A_1010[4] - A_1010[5]);
                                    const double tmp187 = w23*(-A_1212[0] - A_1212[1] - A_1212[4] - A_1212[5]);
                                    const double tmp189 = w7*(A_1212[3] + A_1212[7]);
                                    const double tmp193 = w19*(A_1212[0] + A_1212[4]);
                                    const double tmp201 = w27*(A_1010[2] + A_1010[3] + A_1010[4] + A_1010[5]);
                                    const double tmp204 = w23*(A_1212[0] + A_1212[1] + A_1212[4] + A_1212[5]);
                                    const double tmp205 = w25*(A_1212[2] + A_1212[3] + A_1212[6] + A_1212[7]);
                                    const double tmp208 = w28*(A_1010[0] + A_1010[1]);
                                    const double tmp209 = w29*(A_1010[6] + A_1010[7]);
                                    const double tmp214 = w13*(-A_1212[1] - A_1212[2] - A_1212[5] - A_1212[6]);
                                    const double tmp215 = w22*(-A_1111[0] - A_1111[2] - A_1111[5] - A_1111[7]);
                                    const double tmp217 = w27*(-A_1010[0] - A_1010[1] - A_1010[6] - A_1010[7]);
                                    const double tmp221 = w26*(-A_1111[4] - A_1111[6]);
                                    const double tmp226 = w7*(-A_1212[0] - A_1212[4]);
                                    const double tmp227 = w24*(-A_1111[1] - A_1111[3]);
                                    const double tmp228 = w19*(-A_1212[3] - A_1212[7]);
                                    const double tmp231 = w28*(-A_1010[4] - A_1010[5]);
                                    const double tmp233 = w29*(-A_1010[2] - A_1010[3]);
                                    const double tmp236 = w26*(A_1111[5] + A_1111[7]);
                                    const double tmp238 = w22*(A_1111[1] + A_1111[3] + A_1111[4] + A_1111[6]);
                                    const double tmp244 = w24*(A_1111[0] + A_1111[2]);
                                    const double tmp255 = w26*(-A_1111[1] - A_1111[3]);
                                    const double tmp259 = w7*(-A_1212[3] - A_1212[7]);
                                    const double tmp261 = w24*(-A_1111[4] - A_1111[6]);
                                    const double tmp262 = w19*(-A_1212[0] - A_1212[4]);
                                    const double tmp265 = w28*(-A_1010[2] - A_1010[3]);
                                    const double tmp268 = w29*(-A_1010[4] - A_1010[5]);
                                    const double tmp288 = w13*(-A_1212[0] - A_1212[3] - A_1212[4] - A_1212[7]);
                                    const double tmp289 = w22*(-A_1111[1] - A_1111[3] - A_1111[4] - A_1111[6]);
                                    const double tmp294 = w26*(-A_1111[5] - A_1111[7]);
                                    const double tmp298 = w7*(-A_1212[1] - A_1212[5]);
                                    const double tmp299 = w24*(-A_1111[0] - A_1111[2]);
                                    const double tmp300 = w19*(-A_1212[2] - A_1212[6]);
                                    const double tmp304 = w27*(-A_1010[2] - A_1010[3] - A_1010[4] - A_1010[5]);
                                    const double tmp308 = w26*(-A_1111[0] - A_1111[2]);
                                    const double tmp313 = w24*(-A_1111[5] - A_1111[7]);
                                    const double tmp316 = w28*(-A_1010[0] - A_1010[1]);
                                    const double tmp318 = w29*(-A_1010[6] - A_1010[7]);
                                    const double tmp320 = w26*(A_1111[0] + A_1111[2]);
                                    const double tmp325 = w24*(A_1111[5] + A_1111[7]);
                                    const double tmp329 = w3*(-A_1111[0] - A_1111[2] - A_1111[4] - A_1111[6]);
                                    const double tmp332 = w25*(-A_1212[1] - A_1212[3] - A_1212[5] - A_1212[7]);
                                    const double tmp335 = w9*(-A_1111[1] - A_1111[3] - A_1111[5] - A_1111[7]);
                                    const double tmp337 = w27*(-A_1010[0] - A_1010[1] - A_1010[2] - A_1010[3] - A_1010[4] - A_1010[5] - A_1010[6] - A_1010[7]);
                                    const double tmp338 = w23*(-A_1212[0] - A_1212[2] - A_1212[4] - A_1212[6]);
                                    const double tmp339 = w14*(-A_1010[0] - A_1010[1] - A_1010[4] - A_1010[5]);
                                    const double tmp340 = w23*(-A_1212[2] - A_1212[3] - A_1212[6] - A_1212[7]);
                                    const double tmp342 = w25*(-A_1212[0] - A_1212[1] - A_1212[4] - A_1212[5]);
                                    const double tmp344 = w0*(-A_1010[2] - A_1010[3] - A_1010[6] - A_1010[7]);
                                    const double tmp358 = w7*(-A_1212[2] - A_1212[6]);
                                    const double tmp359 = w19*(-A_1212[1] - A_1212[5]);
                                    const double tmp362 = w28*(-A_1010[6] - A_1010[7]);
                                    const double tmp363 = w29*(-A_1010[0] - A_1010[1]);
                                    const double tmp371 = w3*(A_1111[4] + A_1111[5] + A_1111[6] + A_1111[7]);
                                    const double tmp374 = w9*(A_1111[0] + A_1111[1] + A_1111[2] + A_1111[3]);
                                    const double tmp375 = w29*(A_1010[2] + A_1010[3]);
                                    const double tmp377 = w28*(A_1010[4] + A_1010[5]);
                                    const double tmp422 = w3*(-A_1111[1] - A_1111[3] - A_1111[5] - A_1111[7]);
                                    const double tmp424 = w25*(-A_1212[0] - A_1212[2] - A_1212[4] - A_1212[6]);
                                    const double tmp428 = w9*(-A_1111[0] - A_1111[2] - A_1111[4] - A_1111[6]);
                                    const double tmp430 = w23*(-A_1212[1] - A_1212[3] - A_1212[5] - A_1212[7]);
                                    const double tmp455 = w29*(A_1010[0] + A_1010[1]);
                                    const double tmp456 = w28*(A_1010[6] + A_1010[7]);
                                    EM_S[INDEX4(1,1,0,0,numEq,numComp,8)]+= tmp214 + tmp259 + tmp262 + tmp289 + tmp294 + tmp299 + tmp304 + tmp362 + tmp363;
                                    EM_S[INDEX4(1,1,0,1,numEq,numComp,8)]+= tmp201 + tmp371 + tmp374 + tmp455 + tmp456 + tmp89 + tmp91;
                                    EM_S[INDEX4(1,1,0,2,numEq,numComp,8)]+= tmp236 + tmp238 + tmp244 + tmp39 + tmp41 + tmp61 + tmp71;
                                    EM_S[INDEX4(1,1,0,3,numEq,numComp,8)]+= tmp20 + tmp73 + tmp74 + tmp75 + tmp80;
                                    EM_S[INDEX4(1,1,0,4,numEq,numComp,8)]+= tmp1 + tmp127 + tmp131 + tmp141 + tmp146 + tmp189 + tmp193;
                                    EM_S[INDEX4(1,1,0,5,numEq,numComp,8)]+= tmp176 + tmp339 + tmp340 + tmp342 + tmp344;
                                    EM_S[INDEX4(1,1,0,6,numEq,numComp,8)]+= tmp337 + tmp422 + tmp424 + tmp428 + tmp430;
                                    EM_S[INDEX4(1,1,0,7,numEq,numComp,8)]+= tmp104 + tmp106 + tmp113;
                                    EM_S[INDEX4(1,1,1,0,numEq,numComp,8)]+= tmp201 + tmp371 + tmp374 + tmp455 + tmp456 + tmp89 + tmp91;
                                    EM_S[INDEX4(1,1,1,1,numEq,numComp,8)]+= tmp215 + tmp221 + tmp227 + tmp288 + tmp304 + tmp358 + tmp359 + tmp362 + tmp363;
                                    EM_S[INDEX4(1,1,1,2,numEq,numComp,8)]+= tmp20 + tmp73 + tmp74 + tmp75 + tmp80;
                                    EM_S[INDEX4(1,1,1,3,numEq,numComp,8)]+= tmp39 + tmp40 + tmp41 + tmp43 + tmp45 + tmp52 + tmp55;
                                    EM_S[INDEX4(1,1,1,4,numEq,numComp,8)]+= tmp176 + tmp339 + tmp340 + tmp342 + tmp344;
                                    EM_S[INDEX4(1,1,1,5,numEq,numComp,8)]+= tmp12 + tmp123 + tmp141 + tmp142 + tmp146 + tmp149 + tmp6;
                                    EM_S[INDEX4(1,1,1,6,numEq,numComp,8)]+= tmp104 + tmp106 + tmp113;
                                    EM_S[INDEX4(1,1,1,7,numEq,numComp,8)]+= tmp329 + tmp332 + tmp335 + tmp337 + tmp338;
                                    EM_S[INDEX4(1,1,2,0,numEq,numComp,8)]+= tmp236 + tmp238 + tmp244 + tmp39 + tmp41 + tmp61 + tmp71;
                                    EM_S[INDEX4(1,1,2,1,numEq,numComp,8)]+= tmp20 + tmp73 + tmp74 + tmp75 + tmp80;
                                    EM_S[INDEX4(1,1,2,2,numEq,numComp,8)]+= tmp217 + tmp231 + tmp233 + tmp288 + tmp289 + tmp294 + tmp298 + tmp299 + tmp300;
                                    EM_S[INDEX4(1,1,2,3,numEq,numComp,8)]+= tmp101 + tmp204 + tmp205 + tmp371 + tmp374 + tmp375 + tmp377;
                                    EM_S[INDEX4(1,1,2,4,numEq,numComp,8)]+= tmp337 + tmp422 + tmp424 + tmp428 + tmp430;
                                    EM_S[INDEX4(1,1,2,5,numEq,numComp,8)]+= tmp104 + tmp106 + tmp113;
                                    EM_S[INDEX4(1,1,2,6,numEq,numComp,8)]+= tmp10 + tmp123 + tmp125 + tmp127 + tmp131 + tmp132 + tmp3;
                                    EM_S[INDEX4(1,1,2,7,numEq,numComp,8)]+= tmp175 + tmp176 + tmp178 + tmp180 + tmp187;
                                    EM_S[INDEX4(1,1,3,0,numEq,numComp,8)]+= tmp20 + tmp73 + tmp74 + tmp75 + tmp80;
                                    EM_S[INDEX4(1,1,3,1,numEq,numComp,8)]+= tmp39 + tmp40 + tmp41 + tmp43 + tmp45 + tmp52 + tmp55;
                                    EM_S[INDEX4(1,1,3,2,numEq,numComp,8)]+= tmp101 + tmp204 + tmp205 + tmp371 + tmp374 + tmp375 + tmp377;
                                    EM_S[INDEX4(1,1,3,3,numEq,numComp,8)]+= tmp214 + tmp215 + tmp217 + tmp221 + tmp226 + tmp227 + tmp228 + tmp231 + tmp233;
                                    EM_S[INDEX4(1,1,3,4,numEq,numComp,8)]+= tmp104 + tmp106 + tmp113;
                                    EM_S[INDEX4(1,1,3,5,numEq,numComp,8)]+= tmp329 + tmp332 + tmp335 + tmp337 + tmp338;
                                    EM_S[INDEX4(1,1,3,6,numEq,numComp,8)]+= tmp175 + tmp176 + tmp178 + tmp180 + tmp187;
                                    EM_S[INDEX4(1,1,3,7,numEq,numComp,8)]+= tmp1 + tmp10 + tmp12 + tmp17 + tmp3 + tmp4 + tmp6;
                                    EM_S[INDEX4(1,1,4,0,numEq,numComp,8)]+= tmp1 + tmp127 + tmp131 + tmp141 + tmp146 + tmp189 + tmp193;
                                    EM_S[INDEX4(1,1,4,1,numEq,numComp,8)]+= tmp176 + tmp339 + tmp340 + tmp342 + tmp344;
                                    EM_S[INDEX4(1,1,4,2,numEq,numComp,8)]+= tmp337 + tmp422 + tmp424 + tmp428 + tmp430;
                                    EM_S[INDEX4(1,1,4,3,numEq,numComp,8)]+= tmp104 + tmp106 + tmp113;
                                    EM_S[INDEX4(1,1,4,4,numEq,numComp,8)]+= tmp214 + tmp215 + tmp217 + tmp255 + tmp259 + tmp261 + tmp262 + tmp265 + tmp268;
                                    EM_S[INDEX4(1,1,4,5,numEq,numComp,8)]+= tmp100 + tmp101 + tmp88 + tmp89 + tmp91 + tmp95 + tmp97;
                                    EM_S[INDEX4(1,1,4,6,numEq,numComp,8)]+= tmp43 + tmp57 + tmp58 + tmp61 + tmp64 + tmp66 + tmp71;
                                    EM_S[INDEX4(1,1,4,7,numEq,numComp,8)]+= tmp20 + tmp22 + tmp25 + tmp28 + tmp30;
                                    EM_S[INDEX4(1,1,5,0,numEq,numComp,8)]+= tmp176 + tmp339 + tmp340 + tmp342 + tmp344;
                                    EM_S[INDEX4(1,1,5,1,numEq,numComp,8)]+= tmp12 + tmp123 + tmp141 + tmp142 + tmp146 + tmp149 + tmp6;
                                    EM_S[INDEX4(1,1,5,2,numEq,numComp,8)]+= tmp104 + tmp106 + tmp113;
                                    EM_S[INDEX4(1,1,5,3,numEq,numComp,8)]+= tmp329 + tmp332 + tmp335 + tmp337 + tmp338;
                                    EM_S[INDEX4(1,1,5,4,numEq,numComp,8)]+= tmp100 + tmp101 + tmp88 + tmp89 + tmp91 + tmp95 + tmp97;
                                    EM_S[INDEX4(1,1,5,5,numEq,numComp,8)]+= tmp217 + tmp265 + tmp268 + tmp288 + tmp289 + tmp308 + tmp313 + tmp358 + tmp359;
                                    EM_S[INDEX4(1,1,5,6,numEq,numComp,8)]+= tmp20 + tmp22 + tmp25 + tmp28 + tmp30;
                                    EM_S[INDEX4(1,1,5,7,numEq,numComp,8)]+= tmp238 + tmp320 + tmp325 + tmp45 + tmp55 + tmp57 + tmp64;
                                    EM_S[INDEX4(1,1,6,0,numEq,numComp,8)]+= tmp337 + tmp422 + tmp424 + tmp428 + tmp430;
                                    EM_S[INDEX4(1,1,6,1,numEq,numComp,8)]+= tmp104 + tmp106 + tmp113;
                                    EM_S[INDEX4(1,1,6,2,numEq,numComp,8)]+= tmp10 + tmp123 + tmp125 + tmp127 + tmp131 + tmp132 + tmp3;
                                    EM_S[INDEX4(1,1,6,3,numEq,numComp,8)]+= tmp175 + tmp176 + tmp178 + tmp180 + tmp187;
                                    EM_S[INDEX4(1,1,6,4,numEq,numComp,8)]+= tmp43 + tmp57 + tmp58 + tmp61 + tmp64 + tmp66 + tmp71;
                                    EM_S[INDEX4(1,1,6,5,numEq,numComp,8)]+= tmp20 + tmp22 + tmp25 + tmp28 + tmp30;
                                    EM_S[INDEX4(1,1,6,6,numEq,numComp,8)]+= tmp215 + tmp255 + tmp261 + tmp288 + tmp298 + tmp300 + tmp304 + tmp316 + tmp318;
                                    EM_S[INDEX4(1,1,6,7,numEq,numComp,8)]+= tmp100 + tmp201 + tmp204 + tmp205 + tmp208 + tmp209 + tmp88;
                                    EM_S[INDEX4(1,1,7,0,numEq,numComp,8)]+= tmp104 + tmp106 + tmp113;
                                    EM_S[INDEX4(1,1,7,1,numEq,numComp,8)]+= tmp329 + tmp332 + tmp335 + tmp337 + tmp338;
                                    EM_S[INDEX4(1,1,7,2,numEq,numComp,8)]+= tmp175 + tmp176 + tmp178 + tmp180 + tmp187;
                                    EM_S[INDEX4(1,1,7,3,numEq,numComp,8)]+= tmp1 + tmp10 + tmp12 + tmp17 + tmp3 + tmp4 + tmp6;
                                    EM_S[INDEX4(1,1,7,4,numEq,numComp,8)]+= tmp20 + tmp22 + tmp25 + tmp28 + tmp30;
                                    EM_S[INDEX4(1,1,7,5,numEq,numComp,8)]+= tmp238 + tmp320 + tmp325 + tmp45 + tmp55 + tmp57 + tmp64;
                                    EM_S[INDEX4(1,1,7,6,numEq,numComp,8)]+= tmp100 + tmp201 + tmp204 + tmp205 + tmp208 + tmp209 + tmp88;
                                    EM_S[INDEX4(1,1,7,7,numEq,numComp,8)]+= tmp214 + tmp226 + tmp228 + tmp289 + tmp304 + tmp308 + tmp313 + tmp316 + tmp318;
                                }
                                {
                                    const double tmp5 = w10*(A_1122[1] + A_1122[6] - A_1221[2] - A_1221[5]);
                                    const double tmp13 = w12*(-A_1122[3] - A_1122[5] + A_1221[1] + A_1221[7]);
                                    const double tmp16 = w6*(-A_1122[2] - A_1122[4] + A_1221[0] + A_1221[6]);
                                    const double tmp24 = w10*(A_1122[2] + A_1122[3] + A_1122[4] + A_1122[5] - A_1221[0] - A_1221[1] - A_1221[6] - A_1221[7]);
                                    const double tmp32 = w12*(-A_1122[6] - A_1122[7] + A_1221[4] + A_1221[5]);
                                    const double tmp35 = w6*(-A_1122[0] - A_1122[1] + A_1221[2] + A_1221[3]);
                                    const double tmp42 = w10*(-A_1122[2] - A_1122[5] + A_1221[0] + A_1221[7]);
                                    const double tmp51 = w12*(A_1122[1] + A_1122[7] - A_1221[3] - A_1221[5]);
                                    const double tmp54 = w6*(A_1122[0] + A_1122[6] - A_1221[2] - A_1221[4]);
                                    const double tmp68 = w12*(A_1122[0] + A_1122[6] - A_1221[2] - A_1221[4]);
                                    const double tmp70 = w6*(A_1122[1] + A_1122[7] - A_1221[3] - A_1221[5]);
                                    const double tmp81 = w12*(-A_1122[0] - A_1122[1] + A_1221[2] + A_1221[3]);
                                    const double tmp83 = w6*(-A_1122[6] - A_1122[7] + A_1221[4] + A_1221[5]);
                                    const double tmp84 = w6*(-A_1122[2] - A_1122[3] - A_1221[2] - A_1221[3]);
                                    const double tmp87 = w10*(A_1122[0] + A_1122[1] + A_1122[6] + A_1122[7] + A_1221[0] + A_1221[1] + A_1221[6] + A_1221[7]);
                                    const double tmp96 = w12*(-A_1122[4] - A_1122[5] - A_1221[4] - A_1221[5]);
                                    const double tmp105 = w6*(-A_1122[4] - A_1122[5] - A_1221[2] - A_1221[3]);
                                    const double tmp110 = w12*(-A_1122[2] - A_1122[3] - A_1221[4] - A_1221[5]);
                                    const double tmp126 = w10*(-A_1122[3] - A_1122[4] + A_1221[0] + A_1221[7]);
                                    const double tmp154 = w10*(-A_1122[2] - A_1122[5] + A_1221[1] + A_1221[6]);
                                    const double tmp160 = w6*(A_1122[6] + A_1122[7] + A_1221[0] + A_1221[1]);
                                    const double tmp161 = w10*(-A_1122[2] - A_1122[3] - A_1122[4] - A_1122[5] - A_1221[2] - A_1221[3] - A_1221[4] - A_1221[5]);
                                    const double tmp164 = w12*(A_1122[0] + A_1122[1] + A_1221[6] + A_1221[7]);
                                    const double tmp166 = w10*(-A_1122[0] - A_1122[1] - A_1122[6] - A_1122[7] + A_1221[2] + A_1221[3] + A_1221[4] + A_1221[5]);
                                    const double tmp169 = w12*(A_1122[2] + A_1122[3] - A_1221[0] - A_1221[1]);
                                    const double tmp172 = w6*(A_1122[4] + A_1122[5] - A_1221[6] - A_1221[7]);
                                    const double tmp182 = w12*(-A_1122[6] - A_1122[7] + A_1221[2] + A_1221[3]);
                                    const double tmp185 = w6*(-A_1122[0] - A_1122[1] + A_1221[4] + A_1221[5]);
                                    const double tmp207 = w12*(A_1122[6] + A_1122[7] + A_1221[6] + A_1221[7]);
                                    const double tmp211 = w6*(A_1122[0] + A_1122[1] + A_1221[0] + A_1221[1]);
                                    const double tmp222 = w10*(A_1122[2] + A_1122[5] + A_1221[2] + A_1221[5]);
                                    const double tmp225 = w6*(-A_1122[0] - A_1122[6] - A_1221[0] - A_1221[6]);
                                    const double tmp232 = w12*(-A_1122[1] - A_1122[7] - A_1221[1] - A_1221[7]);
                                    const double tmp237 = w10*(A_1122[1] + A_1122[6] - A_1221[3] - A_1221[4]);
                                    const double tmp243 = w12*(-A_1122[2] - A_1122[4] + A_1221[0] + A_1221[6]);
                                    const double tmp246 = w6*(-A_1122[3] - A_1122[5] + A_1221[1] + A_1221[7]);
                                    const double tmp258 = w6*(-A_1122[1] - A_1122[7] - A_1221[1] - A_1221[7]);
                                    const double tmp266 = w12*(-A_1122[0] - A_1122[6] - A_1221[0] - A_1221[6]);
                                    const double tmp279 = w10*(A_1122[0] + A_1122[7] - A_1221[3] - A_1221[4]);
                                    const double tmp284 = w12*(A_1122[2] + A_1122[3] - A_1221[6] - A_1221[7]);
                                    const double tmp285 = w6*(A_1122[4] + A_1122[5] - A_1221[0] - A_1221[1]);
                                    const double tmp295 = w10*(A_1122[3] + A_1122[4] + A_1221[3] + A_1221[4]);
                                    const double tmp309 = w10*(-A_1122[1] - A_1122[6] - A_1221[1] - A_1221[6]);
                                    const double tmp312 = w6*(A_1122[2] + A_1122[4] + A_1221[2] + A_1221[4]);
                                    const double tmp317 = w12*(A_1122[3] + A_1122[5] + A_1221[3] + A_1221[5]);
                                    const double tmp328 = w10*(-A_1122[1] - A_1122[6] - A_1221[0] - A_1221[7]);
                                    const double tmp346 = w12*(A_1122[4] + A_1122[5] - A_1221[0] - A_1221[1]);
                                    const double tmp348 = w6*(A_1122[2] + A_1122[3] - A_1221[6] - A_1221[7]);
                                    const double tmp355 = w10*(-A_1122[0] - A_1122[7] - A_1221[0] - A_1221[7]);
                                    const double tmp368 = w6*(-A_1122[4] - A_1122[5] - A_1221[4] - A_1221[5]);
                                    const double tmp372 = w12*(-A_1122[2] - A_1122[3] - A_1221[2] - A_1221[3]);
                                    const double tmp383 = w6*(A_1122[3] + A_1122[5] + A_1221[3] + A_1221[5]);
                                    const double tmp386 = w12*(A_1122[2] + A_1122[4] + A_1221[2] + A_1221[4]);
                                    const double tmp398 = w10*(A_1122[3] + A_1122[4] + A_1221[2] + A_1221[5]);
                                    const double tmp416 = w12*(-A_1122[0] - A_1122[1] + A_1221[4] + A_1221[5]);
                                    const double tmp417 = w6*(-A_1122[6] - A_1122[7] + A_1221[2] + A_1221[3]);
                                    const double tmp421 = w10*(A_1122[2] + A_1122[5] + A_1221[3] + A_1221[4]);
                                    const double tmp432 = w10*(-A_1122[3] - A_1122[4] + A_1221[1] + A_1221[6]);
                                    const double tmp446 = w10*(-A_1122[0] - A_1122[7] - A_1221[1] - A_1221[6]);
                                    const double tmp451 = w6*(A_1122[6] + A_1122[7] + A_1221[6] + A_1221[7]);
                                    const double tmp454 = w12*(A_1122[0] + A_1122[1] + A_1221[0] + A_1221[1]);
                                    const double tmp460 = w12*(A_1122[4] + A_1122[5] - A_1221[6] - A_1221[7]);
                                    const double tmp461 = w6*(A_1122[2] + A_1122[3] - A_1221[0] - A_1221[1]);
                                    const double tmp470 = w6*(-A_1122[2] - A_1122[3] - A_1221[4] - A_1221[5]);
                                    const double tmp474 = w12*(-A_1122[4] - A_1122[5] - A_1221[2] - A_1221[3]);
                                    const double tmp478 = w10*(A_1122[0] + A_1122[7] - A_1221[2] - A_1221[5]);
                                    const double tmp482 = w6*(A_1122[0] + A_1122[1] + A_1221[6] + A_1221[7]);
                                    const double tmp483 = w12*(A_1122[6] + A_1122[7] + A_1221[0] + A_1221[1]);
                                    EM_S[INDEX4(1,2,0,0,numEq,numComp,8)]+= tmp309 + tmp383 + w18*(A_1122[0] + A_1221[0]) + w4*(-A_1122[7] - A_1221[7]) + tmp386;
                                    EM_S[INDEX4(1,2,0,1,numEq,numComp,8)]+= tmp161 + tmp451 + tmp454;
                                    EM_S[INDEX4(1,2,0,2,numEq,numComp,8)]+= tmp432 + w18*(A_1122[2] - A_1221[0]) + w4*(-A_1122[5] + A_1221[7]) + tmp68 + tmp70;
                                    EM_S[INDEX4(1,2,0,3,numEq,numComp,8)]+= tmp166 + tmp169 + tmp172;
                                    EM_S[INDEX4(1,2,0,4,numEq,numComp,8)]+= tmp243 + tmp246 + w18*(-A_1122[0] + A_1221[4]) + w4*(A_1122[7] - A_1221[3]) + tmp5;
                                    EM_S[INDEX4(1,2,0,5,numEq,numComp,8)]+= tmp24 + tmp416 + tmp417;
                                    EM_S[INDEX4(1,2,0,6,numEq,numComp,8)]+= tmp258 + tmp266 + tmp398 + w18*(-A_1122[2] - A_1221[4]) + w4*(A_1122[5] + A_1221[3]);
                                    EM_S[INDEX4(1,2,0,7,numEq,numComp,8)]+= tmp105 + tmp110 + tmp87;
                                    EM_S[INDEX4(1,2,1,0,numEq,numComp,8)]+= tmp161 + tmp451 + tmp454;
                                    EM_S[INDEX4(1,2,1,1,numEq,numComp,8)]+= tmp312 + tmp317 + tmp355 + w18*(A_1122[1] + A_1221[1]) + w4*(-A_1122[6] - A_1221[6]);
                                    EM_S[INDEX4(1,2,1,2,numEq,numComp,8)]+= tmp166 + tmp169 + tmp172;
                                    EM_S[INDEX4(1,2,1,3,numEq,numComp,8)]+= w18*(A_1122[3] - A_1221[1]) + tmp42 + w4*(-A_1122[4] + A_1221[6]) + tmp51 + tmp54;
                                    EM_S[INDEX4(1,2,1,4,numEq,numComp,8)]+= tmp24 + tmp416 + tmp417;
                                    EM_S[INDEX4(1,2,1,5,numEq,numComp,8)]+= tmp13 + tmp16 + w18*(-A_1122[1] + A_1221[5]) + tmp279 + w4*(A_1122[6] - A_1221[2]);
                                    EM_S[INDEX4(1,2,1,6,numEq,numComp,8)]+= tmp105 + tmp110 + tmp87;
                                    EM_S[INDEX4(1,2,1,7,numEq,numComp,8)]+= tmp225 + tmp232 + tmp421 + w18*(-A_1122[3] - A_1221[5]) + w4*(A_1122[4] + A_1221[2]);
                                    EM_S[INDEX4(1,2,2,0,numEq,numComp,8)]+= w18*(-A_1122[0] + A_1221[2]) + tmp237 + w4*(A_1122[7] - A_1221[5]) + tmp243 + tmp246;
                                    EM_S[INDEX4(1,2,2,1,numEq,numComp,8)]+= tmp24 + tmp81 + tmp83;
                                    EM_S[INDEX4(1,2,2,2,numEq,numComp,8)]+= tmp258 + tmp266 + tmp295 + w18*(-A_1122[2] - A_1221[2]) + w4*(A_1122[5] + A_1221[5]);
                                    EM_S[INDEX4(1,2,2,3,numEq,numComp,8)]+= tmp368 + tmp372 + tmp87;
                                    EM_S[INDEX4(1,2,2,4,numEq,numComp,8)]+= tmp328 + tmp383 + tmp386 + w18*(A_1122[0] + A_1221[6]) + w4*(-A_1122[7] - A_1221[1]);
                                    EM_S[INDEX4(1,2,2,5,numEq,numComp,8)]+= tmp160 + tmp161 + tmp164;
                                    EM_S[INDEX4(1,2,2,6,numEq,numComp,8)]+= w18*(A_1122[2] - A_1221[6]) + tmp126 + w4*(-A_1122[5] + A_1221[1]) + tmp68 + tmp70;
                                    EM_S[INDEX4(1,2,2,7,numEq,numComp,8)]+= tmp166 + tmp284 + tmp285;
                                    EM_S[INDEX4(1,2,3,0,numEq,numComp,8)]+= tmp24 + tmp81 + tmp83;
                                    EM_S[INDEX4(1,2,3,1,numEq,numComp,8)]+= tmp13 + tmp16 + tmp478 + w18*(-A_1122[1] + A_1221[3]) + w4*(A_1122[6] - A_1221[4]);
                                    EM_S[INDEX4(1,2,3,2,numEq,numComp,8)]+= tmp368 + tmp372 + tmp87;
                                    EM_S[INDEX4(1,2,3,3,numEq,numComp,8)]+= tmp222 + tmp225 + w18*(-A_1122[3] - A_1221[3]) + w4*(A_1122[4] + A_1221[4]) + tmp232;
                                    EM_S[INDEX4(1,2,3,4,numEq,numComp,8)]+= tmp160 + tmp161 + tmp164;
                                    EM_S[INDEX4(1,2,3,5,numEq,numComp,8)]+= tmp312 + tmp317 + tmp446 + w18*(A_1122[1] + A_1221[7]) + w4*(-A_1122[6] - A_1221[0]);
                                    EM_S[INDEX4(1,2,3,6,numEq,numComp,8)]+= tmp166 + tmp284 + tmp285;
                                    EM_S[INDEX4(1,2,3,7,numEq,numComp,8)]+= w18*(A_1122[3] - A_1221[7]) + tmp154 + w4*(-A_1122[4] + A_1221[0]) + tmp51 + tmp54;
                                    EM_S[INDEX4(1,2,4,0,numEq,numComp,8)]+= tmp154 + w18*(A_1122[4] - A_1221[0]) + w4*(-A_1122[3] + A_1221[7]) + tmp68 + tmp70;
                                    EM_S[INDEX4(1,2,4,1,numEq,numComp,8)]+= tmp166 + tmp346 + tmp348;
                                    EM_S[INDEX4(1,2,4,2,numEq,numComp,8)]+= tmp383 + tmp386 + w18*(A_1122[6] + A_1221[0]) + tmp446 + w4*(-A_1122[1] - A_1221[7]);
                                    EM_S[INDEX4(1,2,4,3,numEq,numComp,8)]+= tmp161 + tmp482 + tmp483;
                                    EM_S[INDEX4(1,2,4,4,numEq,numComp,8)]+= tmp222 + tmp258 + w18*(-A_1122[4] - A_1221[4]) + w4*(A_1122[3] + A_1221[3]) + tmp266;
                                    EM_S[INDEX4(1,2,4,5,numEq,numComp,8)]+= tmp84 + tmp87 + tmp96;
                                    EM_S[INDEX4(1,2,4,6,numEq,numComp,8)]+= tmp243 + tmp246 + w18*(-A_1122[6] + A_1221[4]) + tmp478 + w4*(A_1122[1] - A_1221[3]);
                                    EM_S[INDEX4(1,2,4,7,numEq,numComp,8)]+= tmp24 + tmp32 + tmp35;
                                    EM_S[INDEX4(1,2,5,0,numEq,numComp,8)]+= tmp166 + tmp346 + tmp348;
                                    EM_S[INDEX4(1,2,5,1,numEq,numComp,8)]+= tmp126 + w18*(A_1122[5] - A_1221[1]) + w4*(-A_1122[2] + A_1221[6]) + tmp51 + tmp54;
                                    EM_S[INDEX4(1,2,5,2,numEq,numComp,8)]+= tmp161 + tmp482 + tmp483;
                                    EM_S[INDEX4(1,2,5,3,numEq,numComp,8)]+= tmp312 + tmp317 + w18*(A_1122[7] + A_1221[1]) + tmp328 + w4*(-A_1122[0] - A_1221[6]);
                                    EM_S[INDEX4(1,2,5,4,numEq,numComp,8)]+= tmp84 + tmp87 + tmp96;
                                    EM_S[INDEX4(1,2,5,5,numEq,numComp,8)]+= tmp225 + tmp232 + tmp295 + w18*(-A_1122[5] - A_1221[5]) + w4*(A_1122[2] + A_1221[2]);
                                    EM_S[INDEX4(1,2,5,6,numEq,numComp,8)]+= tmp24 + tmp32 + tmp35;
                                    EM_S[INDEX4(1,2,5,7,numEq,numComp,8)]+= tmp13 + tmp16 + tmp237 + w18*(-A_1122[7] + A_1221[5]) + w4*(A_1122[0] - A_1221[2]);
                                    EM_S[INDEX4(1,2,6,0,numEq,numComp,8)]+= tmp258 + tmp266 + w18*(-A_1122[4] - A_1221[2]) + tmp421 + w4*(A_1122[3] + A_1221[5]);
                                    EM_S[INDEX4(1,2,6,1,numEq,numComp,8)]+= tmp470 + tmp474 + tmp87;
                                    EM_S[INDEX4(1,2,6,2,numEq,numComp,8)]+= tmp243 + tmp246 + tmp279 + w18*(-A_1122[6] + A_1221[2]) + w4*(A_1122[1] - A_1221[5]);
                                    EM_S[INDEX4(1,2,6,3,numEq,numComp,8)]+= tmp182 + tmp185 + tmp24;
                                    EM_S[INDEX4(1,2,6,4,numEq,numComp,8)]+= tmp42 + w18*(A_1122[4] - A_1221[6]) + w4*(-A_1122[3] + A_1221[1]) + tmp68 + tmp70;
                                    EM_S[INDEX4(1,2,6,5,numEq,numComp,8)]+= tmp166 + tmp460 + tmp461;
                                    EM_S[INDEX4(1,2,6,6,numEq,numComp,8)]+= tmp355 + tmp383 + tmp386 + w18*(A_1122[6] + A_1221[6]) + w4*(-A_1122[1] - A_1221[1]);
                                    EM_S[INDEX4(1,2,6,7,numEq,numComp,8)]+= tmp161 + tmp207 + tmp211;
                                    EM_S[INDEX4(1,2,7,0,numEq,numComp,8)]+= tmp470 + tmp474 + tmp87;
                                    EM_S[INDEX4(1,2,7,1,numEq,numComp,8)]+= tmp225 + tmp232 + w18*(-A_1122[5] - A_1221[3]) + tmp398 + w4*(A_1122[2] + A_1221[4]);
                                    EM_S[INDEX4(1,2,7,2,numEq,numComp,8)]+= tmp182 + tmp185 + tmp24;
                                    EM_S[INDEX4(1,2,7,3,numEq,numComp,8)]+= w18*(-A_1122[7] + A_1221[3]) + tmp13 + tmp16 + tmp5 + w4*(A_1122[0] - A_1221[4]);
                                    EM_S[INDEX4(1,2,7,4,numEq,numComp,8)]+= tmp166 + tmp460 + tmp461;
                                    EM_S[INDEX4(1,2,7,5,numEq,numComp,8)]+= w18*(A_1122[5] - A_1221[7]) + tmp432 + w4*(-A_1122[2] + A_1221[0]) + tmp51 + tmp54;
                                    EM_S[INDEX4(1,2,7,6,numEq,numComp,8)]+= tmp161 + tmp207 + tmp211;
                                    EM_S[INDEX4(1,2,7,7,numEq,numComp,8)]+= tmp309 + tmp312 + w18*(A_1122[7] + A_1221[7]) + w4*(-A_1122[0] - A_1221[0]) + tmp317;
                                }
                                {
                                    const double tmp2 = w11*(-A_2002[2] - A_2002[5] + A_2200[1] + A_2200[6]);
                                    const double tmp9 = w15*(-A_2002[3] - A_2002[6] + A_2200[2] + A_2200[7]);
                                    const double tmp14 = w5*(-A_2002[1] - A_2002[4] + A_2200[0] + A_2200[5]);
                                    const double tmp21 = w11*(-A_2002[1] - A_2002[3] - A_2002[4] - A_2002[6] + A_2200[0] + A_2200[2] + A_2200[5] + A_2200[7]);
                                    const double tmp27 = w15*(-A_2002[5] - A_2002[7] + A_2200[4] + A_2200[6]);
                                    const double tmp33 = w5*(-A_2002[0] - A_2002[2] + A_2200[1] + A_2200[3]);
                                    const double tmp38 = w11*(-A_2002[0] - A_2002[2] - A_2002[5] - A_2002[7] - A_2200[0] - A_2200[2] - A_2200[5] - A_2200[7]);
                                    const double tmp47 = w15*(-A_2002[1] - A_2002[3] - A_2200[1] - A_2200[3]);
                                    const double tmp50 = w5*(-A_2002[4] - A_2002[6] - A_2200[4] - A_2200[6]);
                                    const double tmp63 = w15*(-A_2002[4] - A_2002[6] - A_2200[4] - A_2200[6]);
                                    const double tmp69 = w5*(-A_2002[1] - A_2002[3] - A_2200[1] - A_2200[3]);
                                    const double tmp77 = w15*(-A_2002[0] - A_2002[2] + A_2200[1] + A_2200[3]);
                                    const double tmp82 = w5*(-A_2002[5] - A_2002[7] + A_2200[4] + A_2200[6]);
                                    const double tmp85 = w11*(A_2002[1] + A_2002[6] - A_2200[0] - A_2200[7]);
                                    const double tmp92 = w15*(A_2002[0] + A_2002[5] - A_2200[1] - A_2200[4]);
                                    const double tmp98 = w5*(A_2002[2] + A_2002[7] - A_2200[3] - A_2200[6]);
                                    const double tmp108 = w15*(-A_2002[1] - A_2002[3] - A_2200[4] - A_2200[6]);
                                    const double tmp111 = w5*(-A_2002[4] - A_2002[6] - A_2200[1] - A_2200[3]);
                                    const double tmp114 = w11*(A_2002[0] + A_2002[2] + A_2002[5] + A_2002[7] - A_2200[1] - A_2200[3] - A_2200[4] - A_2200[6]);
                                    const double tmp117 = w15*(A_2002[4] + A_2002[6] - A_2200[5] - A_2200[7]);
                                    const double tmp119 = w5*(A_2002[1] + A_2002[3] - A_2200[0] - A_2200[2]);
                                    const double tmp124 = w11*(-A_2002[0] - A_2002[7] + A_2200[3] + A_2200[4]);
                                    const double tmp135 = w11*(A_2002[1] + A_2002[3] + A_2002[4] + A_2002[6] + A_2200[1] + A_2200[3] + A_2200[4] + A_2200[6]);
                                    const double tmp137 = w15*(A_2002[0] + A_2002[2] + A_2200[5] + A_2200[7]);
                                    const double tmp139 = w5*(A_2002[5] + A_2002[7] + A_2200[0] + A_2200[2]);
                                    const double tmp145 = w15*(-A_2002[1] - A_2002[4] + A_2200[0] + A_2200[5]);
                                    const double tmp148 = w5*(-A_2002[3] - A_2002[6] + A_2200[2] + A_2200[7]);
                                    const double tmp153 = w11*(A_2002[1] + A_2002[6] - A_2200[2] - A_2200[5]);
                                    const double tmp156 = w15*(A_2002[2] + A_2002[7] - A_2200[3] - A_2200[6]);
                                    const double tmp157 = w5*(A_2002[0] + A_2002[5] - A_2200[1] - A_2200[4]);
                                    const double tmp167 = w15*(A_2002[1] + A_2002[3] - A_2200[0] - A_2200[2]);
                                    const double tmp170 = w5*(A_2002[4] + A_2002[6] - A_2200[5] - A_2200[7]);
                                    const double tmp174 = w11*(-A_2002[3] - A_2002[4] - A_2200[1] - A_2200[6]);
                                    const double tmp179 = w15*(-A_2002[2] - A_2002[7] - A_2200[2] - A_2200[7]);
                                    const double tmp183 = w5*(-A_2002[0] - A_2002[5] - A_2200[0] - A_2200[5]);
                                    const double tmp202 = w11*(-A_2002[2] - A_2002[5] + A_2200[3] + A_2200[4]);
                                    const double tmp220 = w11*(-A_2002[1] - A_2002[6] - A_2200[1] - A_2200[6]);
                                    const double tmp240 = w15*(A_2002[0] + A_2002[2] + A_2200[0] + A_2200[2]);
                                    const double tmp242 = w5*(A_2002[5] + A_2002[7] + A_2200[5] + A_2200[7]);
                                    const double tmp247 = w11*(A_2002[3] + A_2002[4] - A_2200[2] - A_2200[5]);
                                    const double tmp260 = w15*(-A_2002[0] - A_2002[5] - A_2200[0] - A_2200[5]);
                                    const double tmp267 = w5*(-A_2002[2] - A_2002[7] - A_2200[2] - A_2200[7]);
                                    const double tmp269 = w11*(A_2002[2] + A_2002[5] + A_2200[0] + A_2200[7]);
                                    const double tmp271 = w15*(A_2002[3] + A_2002[6] + A_2200[3] + A_2200[6]);
                                    const double tmp273 = w5*(A_2002[1] + A_2002[4] + A_2200[1] + A_2200[4]);
                                    const double tmp278 = w11*(A_2002[3] + A_2002[4] - A_2200[0] - A_2200[7]);
                                    const double tmp283 = w11*(A_2002[0] + A_2002[7] + A_2200[2] + A_2200[5]);
                                    const double tmp293 = w11*(A_2002[0] + A_2002[7] + A_2200[0] + A_2200[7]);
                                    const double tmp307 = w11*(A_2002[2] + A_2002[5] + A_2200[2] + A_2200[5]);
                                    const double tmp324 = w15*(A_2002[5] + A_2002[7] + A_2200[5] + A_2200[7]);
                                    const double tmp326 = w5*(A_2002[0] + A_2002[2] + A_2200[0] + A_2200[2]);
                                    const double tmp333 = w15*(-A_2002[5] - A_2002[7] + A_2200[1] + A_2200[3]);
                                    const double tmp336 = w5*(-A_2002[0] - A_2002[2] + A_2200[4] + A_2200[6]);
                                    const double tmp343 = w15*(A_2002[1] + A_2002[4] + A_2200[1] + A_2200[4]);
                                    const double tmp347 = w5*(A_2002[3] + A_2002[6] + A_2200[3] + A_2200[6]);
                                    const double tmp354 = w11*(-A_2002[3] - A_2002[4] - A_2200[3] - A_2200[4]);
                                    const double tmp365 = w11*(-A_2002[1] - A_2002[6] - A_2200[3] - A_2200[4]);
                                    const double tmp369 = w11*(-A_2002[0] - A_2002[7] + A_2200[1] + A_2200[6]);
                                    const double tmp426 = w15*(A_2002[4] + A_2002[6] - A_2200[0] - A_2200[2]);
                                    const double tmp429 = w5*(A_2002[1] + A_2002[3] - A_2200[5] - A_2200[7]);
                                    const double tmp464 = w15*(A_2002[1] + A_2002[3] - A_2200[5] - A_2200[7]);
                                    const double tmp465 = w5*(A_2002[4] + A_2002[6] - A_2200[0] - A_2200[2]);
                                    const double tmp472 = w15*(-A_2002[4] - A_2002[6] - A_2200[1] - A_2200[3]);
                                    const double tmp475 = w5*(-A_2002[1] - A_2002[3] - A_2200[4] - A_2200[6]);
                                    const double tmp484 = w15*(A_2002[5] + A_2002[7] + A_2200[0] + A_2200[2]);
                                    const double tmp485 = w5*(A_2002[0] + A_2002[2] + A_2200[5] + A_2200[7]);
                                    const double tmp498 = w15*(-A_2002[0] - A_2002[2] + A_2200[4] + A_2200[6]);
                                    const double tmp499 = w5*(-A_2002[5] - A_2002[7] + A_2200[1] + A_2200[3]);
                                    EM_S[INDEX4(2,0,0,0,numEq,numComp,8)]+= tmp307 + tmp343 + tmp347 + w17*(A_2002[0] + A_2200[0]) + w2*(-A_2002[7] - A_2200[7]);
                                    EM_S[INDEX4(2,0,0,1,numEq,numComp,8)]+= tmp247 + w2*(-A_2002[6] + A_2200[7]) + w17*(A_2002[1] - A_2200[0]) + tmp92 + tmp98;
                                    EM_S[INDEX4(2,0,0,2,numEq,numComp,8)]+= tmp135 + tmp240 + tmp242;
                                    EM_S[INDEX4(2,0,0,3,numEq,numComp,8)]+= tmp114 + tmp167 + tmp170;
                                    EM_S[INDEX4(2,0,0,4,numEq,numComp,8)]+= tmp145 + tmp148 + tmp2 + w17*(-A_2002[0] + A_2200[4]) + w2*(A_2002[7] - A_2200[3]);
                                    EM_S[INDEX4(2,0,0,5,numEq,numComp,8)]+= tmp174 + tmp260 + tmp267 + w2*(A_2002[6] + A_2200[3]) + w17*(-A_2002[1] - A_2200[4]);
                                    EM_S[INDEX4(2,0,0,6,numEq,numComp,8)]+= tmp21 + tmp498 + tmp499;
                                    EM_S[INDEX4(2,0,0,7,numEq,numComp,8)]+= tmp108 + tmp111 + tmp38;
                                    EM_S[INDEX4(2,0,1,0,numEq,numComp,8)]+= tmp145 + tmp148 + tmp202 + w2*(A_2002[7] - A_2200[6]) + w17*(-A_2002[0] + A_2200[1]);
                                    EM_S[INDEX4(2,0,1,1,numEq,numComp,8)]+= tmp260 + tmp267 + w17*(-A_2002[1] - A_2200[1]) + w2*(A_2002[6] + A_2200[6]) + tmp354;
                                    EM_S[INDEX4(2,0,1,2,numEq,numComp,8)]+= tmp21 + tmp77 + tmp82;
                                    EM_S[INDEX4(2,0,1,3,numEq,numComp,8)]+= tmp38 + tmp47 + tmp50;
                                    EM_S[INDEX4(2,0,1,4,numEq,numComp,8)]+= tmp269 + tmp343 + tmp347 + w17*(A_2002[0] + A_2200[5]) + w2*(-A_2002[7] - A_2200[2]);
                                    EM_S[INDEX4(2,0,1,5,numEq,numComp,8)]+= tmp278 + w17*(A_2002[1] - A_2200[5]) + w2*(-A_2002[6] + A_2200[2]) + tmp92 + tmp98;
                                    EM_S[INDEX4(2,0,1,6,numEq,numComp,8)]+= tmp135 + tmp137 + tmp139;
                                    EM_S[INDEX4(2,0,1,7,numEq,numComp,8)]+= tmp114 + tmp464 + tmp465;
                                    EM_S[INDEX4(2,0,2,0,numEq,numComp,8)]+= tmp135 + tmp240 + tmp242;
                                    EM_S[INDEX4(2,0,2,1,numEq,numComp,8)]+= tmp114 + tmp167 + tmp170;
                                    EM_S[INDEX4(2,0,2,2,numEq,numComp,8)]+= tmp271 + tmp273 + w17*(A_2002[2] + A_2200[2]) + w2*(-A_2002[5] - A_2200[5]) + tmp293;
                                    EM_S[INDEX4(2,0,2,3,numEq,numComp,8)]+= tmp156 + tmp157 + w2*(-A_2002[4] + A_2200[5]) + w17*(A_2002[3] - A_2200[2]) + tmp85;
                                    EM_S[INDEX4(2,0,2,4,numEq,numComp,8)]+= tmp21 + tmp498 + tmp499;
                                    EM_S[INDEX4(2,0,2,5,numEq,numComp,8)]+= tmp108 + tmp111 + tmp38;
                                    EM_S[INDEX4(2,0,2,6,numEq,numComp,8)]+= tmp124 + w17*(-A_2002[2] + A_2200[6]) + w2*(A_2002[5] - A_2200[1]) + tmp14 + tmp9;
                                    EM_S[INDEX4(2,0,2,7,numEq,numComp,8)]+= tmp179 + tmp183 + w2*(A_2002[4] + A_2200[1]) + tmp365 + w17*(-A_2002[3] - A_2200[6]);
                                    EM_S[INDEX4(2,0,3,0,numEq,numComp,8)]+= tmp21 + tmp77 + tmp82;
                                    EM_S[INDEX4(2,0,3,1,numEq,numComp,8)]+= tmp38 + tmp47 + tmp50;
                                    EM_S[INDEX4(2,0,3,2,numEq,numComp,8)]+= tmp14 + w2*(A_2002[5] - A_2200[4]) + tmp369 + w17*(-A_2002[2] + A_2200[3]) + tmp9;
                                    EM_S[INDEX4(2,0,3,3,numEq,numComp,8)]+= tmp179 + tmp183 + w17*(-A_2002[3] - A_2200[3]) + w2*(A_2002[4] + A_2200[4]) + tmp220;
                                    EM_S[INDEX4(2,0,3,4,numEq,numComp,8)]+= tmp135 + tmp137 + tmp139;
                                    EM_S[INDEX4(2,0,3,5,numEq,numComp,8)]+= tmp114 + tmp464 + tmp465;
                                    EM_S[INDEX4(2,0,3,6,numEq,numComp,8)]+= tmp271 + tmp273 + tmp283 + w17*(A_2002[2] + A_2200[7]) + w2*(-A_2002[5] - A_2200[0]);
                                    EM_S[INDEX4(2,0,3,7,numEq,numComp,8)]+= tmp153 + tmp156 + tmp157 + w17*(A_2002[3] - A_2200[7]) + w2*(-A_2002[4] + A_2200[0]);
                                    EM_S[INDEX4(2,0,4,0,numEq,numComp,8)]+= tmp153 + w17*(A_2002[4] - A_2200[0]) + w2*(-A_2002[3] + A_2200[7]) + tmp92 + tmp98;
                                    EM_S[INDEX4(2,0,4,1,numEq,numComp,8)]+= tmp283 + tmp343 + tmp347 + w17*(A_2002[5] + A_2200[0]) + w2*(-A_2002[2] - A_2200[7]);
                                    EM_S[INDEX4(2,0,4,2,numEq,numComp,8)]+= tmp114 + tmp426 + tmp429;
                                    EM_S[INDEX4(2,0,4,3,numEq,numComp,8)]+= tmp135 + tmp484 + tmp485;
                                    EM_S[INDEX4(2,0,4,4,numEq,numComp,8)]+= tmp220 + w17*(-A_2002[4] - A_2200[4]) + w2*(A_2002[3] + A_2200[3]) + tmp260 + tmp267;
                                    EM_S[INDEX4(2,0,4,5,numEq,numComp,8)]+= tmp145 + tmp148 + tmp369 + w17*(-A_2002[5] + A_2200[4]) + w2*(A_2002[2] - A_2200[3]);
                                    EM_S[INDEX4(2,0,4,6,numEq,numComp,8)]+= tmp38 + tmp63 + tmp69;
                                    EM_S[INDEX4(2,0,4,7,numEq,numComp,8)]+= tmp21 + tmp27 + tmp33;
                                    EM_S[INDEX4(2,0,5,0,numEq,numComp,8)]+= tmp260 + tmp267 + tmp365 + w2*(A_2002[3] + A_2200[6]) + w17*(-A_2002[4] - A_2200[1]);
                                    EM_S[INDEX4(2,0,5,1,numEq,numComp,8)]+= tmp124 + tmp145 + tmp148 + w17*(-A_2002[5] + A_2200[1]) + w2*(A_2002[2] - A_2200[6]);
                                    EM_S[INDEX4(2,0,5,2,numEq,numComp,8)]+= tmp38 + tmp472 + tmp475;
                                    EM_S[INDEX4(2,0,5,3,numEq,numComp,8)]+= tmp21 + tmp333 + tmp336;
                                    EM_S[INDEX4(2,0,5,4,numEq,numComp,8)]+= w17*(A_2002[4] - A_2200[5]) + w2*(-A_2002[3] + A_2200[2]) + tmp85 + tmp92 + tmp98;
                                    EM_S[INDEX4(2,0,5,5,numEq,numComp,8)]+= tmp293 + tmp343 + tmp347 + w17*(A_2002[5] + A_2200[5]) + w2*(-A_2002[2] - A_2200[2]);
                                    EM_S[INDEX4(2,0,5,6,numEq,numComp,8)]+= tmp114 + tmp117 + tmp119;
                                    EM_S[INDEX4(2,0,5,7,numEq,numComp,8)]+= tmp135 + tmp324 + tmp326;
                                    EM_S[INDEX4(2,0,6,0,numEq,numComp,8)]+= tmp114 + tmp426 + tmp429;
                                    EM_S[INDEX4(2,0,6,1,numEq,numComp,8)]+= tmp135 + tmp484 + tmp485;
                                    EM_S[INDEX4(2,0,6,2,numEq,numComp,8)]+= tmp156 + tmp157 + tmp278 + w17*(A_2002[6] - A_2200[2]) + w2*(-A_2002[1] + A_2200[5]);
                                    EM_S[INDEX4(2,0,6,3,numEq,numComp,8)]+= tmp269 + tmp271 + tmp273 + w17*(A_2002[7] + A_2200[2]) + w2*(-A_2002[0] - A_2200[5]);
                                    EM_S[INDEX4(2,0,6,4,numEq,numComp,8)]+= tmp38 + tmp63 + tmp69;
                                    EM_S[INDEX4(2,0,6,5,numEq,numComp,8)]+= tmp21 + tmp27 + tmp33;
                                    EM_S[INDEX4(2,0,6,6,numEq,numComp,8)]+= tmp179 + tmp183 + tmp354 + w17*(-A_2002[6] - A_2200[6]) + w2*(A_2002[1] + A_2200[1]);
                                    EM_S[INDEX4(2,0,6,7,numEq,numComp,8)]+= tmp14 + tmp202 + w17*(-A_2002[7] + A_2200[6]) + w2*(A_2002[0] - A_2200[1]) + tmp9;
                                    EM_S[INDEX4(2,0,7,0,numEq,numComp,8)]+= tmp38 + tmp472 + tmp475;
                                    EM_S[INDEX4(2,0,7,1,numEq,numComp,8)]+= tmp21 + tmp333 + tmp336;
                                    EM_S[INDEX4(2,0,7,2,numEq,numComp,8)]+= w2*(A_2002[1] + A_2200[4]) + tmp174 + tmp179 + tmp183 + w17*(-A_2002[6] - A_2200[3]);
                                    EM_S[INDEX4(2,0,7,3,numEq,numComp,8)]+= tmp14 + w17*(-A_2002[7] + A_2200[3]) + w2*(A_2002[0] - A_2200[4]) + tmp2 + tmp9;
                                    EM_S[INDEX4(2,0,7,4,numEq,numComp,8)]+= tmp114 + tmp117 + tmp119;
                                    EM_S[INDEX4(2,0,7,5,numEq,numComp,8)]+= tmp135 + tmp324 + tmp326;
                                    EM_S[INDEX4(2,0,7,6,numEq,numComp,8)]+= tmp156 + tmp157 + tmp247 + w17*(A_2002[6] - A_2200[7]) + w2*(-A_2002[1] + A_2200[0]);
                                    EM_S[INDEX4(2,0,7,7,numEq,numComp,8)]+= tmp271 + tmp273 + w17*(A_2002[7] + A_2200[7]) + w2*(-A_2002[0] - A_2200[0]) + tmp307;
                                }
                                {
                                    const double tmp5 = w10*(A_2112[1] + A_2112[6] - A_2211[2] - A_2211[5]);
                                    const double tmp13 = w12*(-A_2112[3] - A_2112[5] + A_2211[1] + A_2211[7]);
                                    const double tmp16 = w6*(-A_2112[2] - A_2112[4] + A_2211[0] + A_2211[6]);
                                    const double tmp24 = w10*(A_2112[2] + A_2112[3] + A_2112[4] + A_2112[5] - A_2211[0] - A_2211[1] - A_2211[6] - A_2211[7]);
                                    const double tmp32 = w12*(-A_2112[6] - A_2112[7] + A_2211[4] + A_2211[5]);
                                    const double tmp35 = w6*(-A_2112[0] - A_2112[1] + A_2211[2] + A_2211[3]);
                                    const double tmp42 = w10*(-A_2112[2] - A_2112[5] + A_2211[0] + A_2211[7]);
                                    const double tmp51 = w12*(A_2112[1] + A_2112[7] - A_2211[3] - A_2211[5]);
                                    const double tmp54 = w6*(A_2112[0] + A_2112[6] - A_2211[2] - A_2211[4]);
                                    const double tmp68 = w12*(A_2112[0] + A_2112[6] - A_2211[2] - A_2211[4]);
                                    const double tmp70 = w6*(A_2112[1] + A_2112[7] - A_2211[3] - A_2211[5]);
                                    const double tmp81 = w12*(-A_2112[0] - A_2112[1] + A_2211[2] + A_2211[3]);
                                    const double tmp83 = w6*(-A_2112[6] - A_2112[7] + A_2211[4] + A_2211[5]);
                                    const double tmp84 = w6*(-A_2112[2] - A_2112[3] - A_2211[2] - A_2211[3]);
                                    const double tmp87 = w10*(A_2112[0] + A_2112[1] + A_2112[6] + A_2112[7] + A_2211[0] + A_2211[1] + A_2211[6] + A_2211[7]);
                                    const double tmp96 = w12*(-A_2112[4] - A_2112[5] - A_2211[4] - A_2211[5]);
                                    const double tmp105 = w6*(-A_2112[4] - A_2112[5] - A_2211[2] - A_2211[3]);
                                    const double tmp110 = w12*(-A_2112[2] - A_2112[3] - A_2211[4] - A_2211[5]);
                                    const double tmp126 = w10*(-A_2112[3] - A_2112[4] + A_2211[0] + A_2211[7]);
                                    const double tmp154 = w10*(-A_2112[2] - A_2112[5] + A_2211[1] + A_2211[6]);
                                    const double tmp160 = w6*(A_2112[6] + A_2112[7] + A_2211[0] + A_2211[1]);
                                    const double tmp161 = w10*(-A_2112[2] - A_2112[3] - A_2112[4] - A_2112[5] - A_2211[2] - A_2211[3] - A_2211[4] - A_2211[5]);
                                    const double tmp164 = w12*(A_2112[0] + A_2112[1] + A_2211[6] + A_2211[7]);
                                    const double tmp166 = w10*(-A_2112[0] - A_2112[1] - A_2112[6] - A_2112[7] + A_2211[2] + A_2211[3] + A_2211[4] + A_2211[5]);
                                    const double tmp169 = w12*(A_2112[2] + A_2112[3] - A_2211[0] - A_2211[1]);
                                    const double tmp172 = w6*(A_2112[4] + A_2112[5] - A_2211[6] - A_2211[7]);
                                    const double tmp182 = w12*(-A_2112[6] - A_2112[7] + A_2211[2] + A_2211[3]);
                                    const double tmp185 = w6*(-A_2112[0] - A_2112[1] + A_2211[4] + A_2211[5]);
                                    const double tmp207 = w12*(A_2112[6] + A_2112[7] + A_2211[6] + A_2211[7]);
                                    const double tmp211 = w6*(A_2112[0] + A_2112[1] + A_2211[0] + A_2211[1]);
                                    const double tmp222 = w10*(A_2112[2] + A_2112[5] + A_2211[2] + A_2211[5]);
                                    const double tmp225 = w6*(-A_2112[0] - A_2112[6] - A_2211[0] - A_2211[6]);
                                    const double tmp232 = w12*(-A_2112[1] - A_2112[7] - A_2211[1] - A_2211[7]);
                                    const double tmp237 = w10*(A_2112[1] + A_2112[6] - A_2211[3] - A_2211[4]);
                                    const double tmp243 = w12*(-A_2112[2] - A_2112[4] + A_2211[0] + A_2211[6]);
                                    const double tmp246 = w6*(-A_2112[3] - A_2112[5] + A_2211[1] + A_2211[7]);
                                    const double tmp258 = w6*(-A_2112[1] - A_2112[7] - A_2211[1] - A_2211[7]);
                                    const double tmp266 = w12*(-A_2112[0] - A_2112[6] - A_2211[0] - A_2211[6]);
                                    const double tmp279 = w10*(A_2112[0] + A_2112[7] - A_2211[3] - A_2211[4]);
                                    const double tmp284 = w12*(A_2112[2] + A_2112[3] - A_2211[6] - A_2211[7]);
                                    const double tmp285 = w6*(A_2112[4] + A_2112[5] - A_2211[0] - A_2211[1]);
                                    const double tmp295 = w10*(A_2112[3] + A_2112[4] + A_2211[3] + A_2211[4]);
                                    const double tmp309 = w10*(-A_2112[1] - A_2112[6] - A_2211[1] - A_2211[6]);
                                    const double tmp312 = w6*(A_2112[2] + A_2112[4] + A_2211[2] + A_2211[4]);
                                    const double tmp317 = w12*(A_2112[3] + A_2112[5] + A_2211[3] + A_2211[5]);
                                    const double tmp328 = w10*(-A_2112[1] - A_2112[6] - A_2211[0] - A_2211[7]);
                                    const double tmp346 = w12*(A_2112[4] + A_2112[5] - A_2211[0] - A_2211[1]);
                                    const double tmp348 = w6*(A_2112[2] + A_2112[3] - A_2211[6] - A_2211[7]);
                                    const double tmp355 = w10*(-A_2112[0] - A_2112[7] - A_2211[0] - A_2211[7]);
                                    const double tmp368 = w6*(-A_2112[4] - A_2112[5] - A_2211[4] - A_2211[5]);
                                    const double tmp372 = w12*(-A_2112[2] - A_2112[3] - A_2211[2] - A_2211[3]);
                                    const double tmp383 = w6*(A_2112[3] + A_2112[5] + A_2211[3] + A_2211[5]);
                                    const double tmp386 = w12*(A_2112[2] + A_2112[4] + A_2211[2] + A_2211[4]);
                                    const double tmp398 = w10*(A_2112[3] + A_2112[4] + A_2211[2] + A_2211[5]);
                                    const double tmp416 = w12*(-A_2112[0] - A_2112[1] + A_2211[4] + A_2211[5]);
                                    const double tmp417 = w6*(-A_2112[6] - A_2112[7] + A_2211[2] + A_2211[3]);
                                    const double tmp421 = w10*(A_2112[2] + A_2112[5] + A_2211[3] + A_2211[4]);
                                    const double tmp432 = w10*(-A_2112[3] - A_2112[4] + A_2211[1] + A_2211[6]);
                                    const double tmp446 = w10*(-A_2112[0] - A_2112[7] - A_2211[1] - A_2211[6]);
                                    const double tmp451 = w6*(A_2112[6] + A_2112[7] + A_2211[6] + A_2211[7]);
                                    const double tmp454 = w12*(A_2112[0] + A_2112[1] + A_2211[0] + A_2211[1]);
                                    const double tmp460 = w12*(A_2112[4] + A_2112[5] - A_2211[6] - A_2211[7]);
                                    const double tmp461 = w6*(A_2112[2] + A_2112[3] - A_2211[0] - A_2211[1]);
                                    const double tmp470 = w6*(-A_2112[2] - A_2112[3] - A_2211[4] - A_2211[5]);
                                    const double tmp474 = w12*(-A_2112[4] - A_2112[5] - A_2211[2] - A_2211[3]);
                                    const double tmp478 = w10*(A_2112[0] + A_2112[7] - A_2211[2] - A_2211[5]);
                                    const double tmp482 = w6*(A_2112[0] + A_2112[1] + A_2211[6] + A_2211[7]);
                                    const double tmp483 = w12*(A_2112[6] + A_2112[7] + A_2211[0] + A_2211[1]);
                                    EM_S[INDEX4(2,1,0,0,numEq,numComp,8)]+= tmp309 + tmp383 + w18*(A_2112[0] + A_2211[0]) + w4*(-A_2112[7] - A_2211[7]) + tmp386;
                                    EM_S[INDEX4(2,1,0,1,numEq,numComp,8)]+= tmp161 + tmp451 + tmp454;
                                    EM_S[INDEX4(2,1,0,2,numEq,numComp,8)]+= tmp432 + w18*(A_2112[2] - A_2211[0]) + w4*(-A_2112[5] + A_2211[7]) + tmp68 + tmp70;
                                    EM_S[INDEX4(2,1,0,3,numEq,numComp,8)]+= tmp166 + tmp169 + tmp172;
                                    EM_S[INDEX4(2,1,0,4,numEq,numComp,8)]+= tmp243 + tmp246 + w18*(-A_2112[0] + A_2211[4]) + w4*(A_2112[7] - A_2211[3]) + tmp5;
                                    EM_S[INDEX4(2,1,0,5,numEq,numComp,8)]+= tmp24 + tmp416 + tmp417;
                                    EM_S[INDEX4(2,1,0,6,numEq,numComp,8)]+= tmp258 + tmp266 + tmp398 + w18*(-A_2112[2] - A_2211[4]) + w4*(A_2112[5] + A_2211[3]);
                                    EM_S[INDEX4(2,1,0,7,numEq,numComp,8)]+= tmp105 + tmp110 + tmp87;
                                    EM_S[INDEX4(2,1,1,0,numEq,numComp,8)]+= tmp161 + tmp451 + tmp454;
                                    EM_S[INDEX4(2,1,1,1,numEq,numComp,8)]+= tmp312 + tmp317 + tmp355 + w18*(A_2112[1] + A_2211[1]) + w4*(-A_2112[6] - A_2211[6]);
                                    EM_S[INDEX4(2,1,1,2,numEq,numComp,8)]+= tmp166 + tmp169 + tmp172;
                                    EM_S[INDEX4(2,1,1,3,numEq,numComp,8)]+= w18*(A_2112[3] - A_2211[1]) + tmp42 + w4*(-A_2112[4] + A_2211[6]) + tmp51 + tmp54;
                                    EM_S[INDEX4(2,1,1,4,numEq,numComp,8)]+= tmp24 + tmp416 + tmp417;
                                    EM_S[INDEX4(2,1,1,5,numEq,numComp,8)]+= tmp13 + tmp16 + w18*(-A_2112[1] + A_2211[5]) + tmp279 + w4*(A_2112[6] - A_2211[2]);
                                    EM_S[INDEX4(2,1,1,6,numEq,numComp,8)]+= tmp105 + tmp110 + tmp87;
                                    EM_S[INDEX4(2,1,1,7,numEq,numComp,8)]+= tmp225 + tmp232 + tmp421 + w18*(-A_2112[3] - A_2211[5]) + w4*(A_2112[4] + A_2211[2]);
                                    EM_S[INDEX4(2,1,2,0,numEq,numComp,8)]+= w18*(-A_2112[0] + A_2211[2]) + tmp237 + w4*(A_2112[7] - A_2211[5]) + tmp243 + tmp246;
                                    EM_S[INDEX4(2,1,2,1,numEq,numComp,8)]+= tmp24 + tmp81 + tmp83;
                                    EM_S[INDEX4(2,1,2,2,numEq,numComp,8)]+= tmp258 + tmp266 + tmp295 + w18*(-A_2112[2] - A_2211[2]) + w4*(A_2112[5] + A_2211[5]);
                                    EM_S[INDEX4(2,1,2,3,numEq,numComp,8)]+= tmp368 + tmp372 + tmp87;
                                    EM_S[INDEX4(2,1,2,4,numEq,numComp,8)]+= tmp328 + tmp383 + tmp386 + w18*(A_2112[0] + A_2211[6]) + w4*(-A_2112[7] - A_2211[1]);
                                    EM_S[INDEX4(2,1,2,5,numEq,numComp,8)]+= tmp160 + tmp161 + tmp164;
                                    EM_S[INDEX4(2,1,2,6,numEq,numComp,8)]+= w18*(A_2112[2] - A_2211[6]) + tmp126 + w4*(-A_2112[5] + A_2211[1]) + tmp68 + tmp70;
                                    EM_S[INDEX4(2,1,2,7,numEq,numComp,8)]+= tmp166 + tmp284 + tmp285;
                                    EM_S[INDEX4(2,1,3,0,numEq,numComp,8)]+= tmp24 + tmp81 + tmp83;
                                    EM_S[INDEX4(2,1,3,1,numEq,numComp,8)]+= tmp13 + tmp16 + tmp478 + w18*(-A_2112[1] + A_2211[3]) + w4*(A_2112[6] - A_2211[4]);
                                    EM_S[INDEX4(2,1,3,2,numEq,numComp,8)]+= tmp368 + tmp372 + tmp87;
                                    EM_S[INDEX4(2,1,3,3,numEq,numComp,8)]+= tmp222 + tmp225 + w18*(-A_2112[3] - A_2211[3]) + w4*(A_2112[4] + A_2211[4]) + tmp232;
                                    EM_S[INDEX4(2,1,3,4,numEq,numComp,8)]+= tmp160 + tmp161 + tmp164;
                                    EM_S[INDEX4(2,1,3,5,numEq,numComp,8)]+= tmp312 + tmp317 + tmp446 + w18*(A_2112[1] + A_2211[7]) + w4*(-A_2112[6] - A_2211[0]);
                                    EM_S[INDEX4(2,1,3,6,numEq,numComp,8)]+= tmp166 + tmp284 + tmp285;
                                    EM_S[INDEX4(2,1,3,7,numEq,numComp,8)]+= w18*(A_2112[3] - A_2211[7]) + tmp154 + w4*(-A_2112[4] + A_2211[0]) + tmp51 + tmp54;
                                    EM_S[INDEX4(2,1,4,0,numEq,numComp,8)]+= tmp154 + w18*(A_2112[4] - A_2211[0]) + w4*(-A_2112[3] + A_2211[7]) + tmp68 + tmp70;
                                    EM_S[INDEX4(2,1,4,1,numEq,numComp,8)]+= tmp166 + tmp346 + tmp348;
                                    EM_S[INDEX4(2,1,4,2,numEq,numComp,8)]+= tmp383 + tmp386 + w18*(A_2112[6] + A_2211[0]) + tmp446 + w4*(-A_2112[1] - A_2211[7]);
                                    EM_S[INDEX4(2,1,4,3,numEq,numComp,8)]+= tmp161 + tmp482 + tmp483;
                                    EM_S[INDEX4(2,1,4,4,numEq,numComp,8)]+= tmp222 + tmp258 + w18*(-A_2112[4] - A_2211[4]) + w4*(A_2112[3] + A_2211[3]) + tmp266;
                                    EM_S[INDEX4(2,1,4,5,numEq,numComp,8)]+= tmp84 + tmp87 + tmp96;
                                    EM_S[INDEX4(2,1,4,6,numEq,numComp,8)]+= tmp243 + tmp246 + w18*(-A_2112[6] + A_2211[4]) + tmp478 + w4*(A_2112[1] - A_2211[3]);
                                    EM_S[INDEX4(2,1,4,7,numEq,numComp,8)]+= tmp24 + tmp32 + tmp35;
                                    EM_S[INDEX4(2,1,5,0,numEq,numComp,8)]+= tmp166 + tmp346 + tmp348;
                                    EM_S[INDEX4(2,1,5,1,numEq,numComp,8)]+= tmp126 + w18*(A_2112[5] - A_2211[1]) + w4*(-A_2112[2] + A_2211[6]) + tmp51 + tmp54;
                                    EM_S[INDEX4(2,1,5,2,numEq,numComp,8)]+= tmp161 + tmp482 + tmp483;
                                    EM_S[INDEX4(2,1,5,3,numEq,numComp,8)]+= tmp312 + tmp317 + w18*(A_2112[7] + A_2211[1]) + tmp328 + w4*(-A_2112[0] - A_2211[6]);
                                    EM_S[INDEX4(2,1,5,4,numEq,numComp,8)]+= tmp84 + tmp87 + tmp96;
                                    EM_S[INDEX4(2,1,5,5,numEq,numComp,8)]+= tmp225 + tmp232 + tmp295 + w18*(-A_2112[5] - A_2211[5]) + w4*(A_2112[2] + A_2211[2]);
                                    EM_S[INDEX4(2,1,5,6,numEq,numComp,8)]+= tmp24 + tmp32 + tmp35;
                                    EM_S[INDEX4(2,1,5,7,numEq,numComp,8)]+= tmp13 + tmp16 + tmp237 + w18*(-A_2112[7] + A_2211[5]) + w4*(A_2112[0] - A_2211[2]);
                                    EM_S[INDEX4(2,1,6,0,numEq,numComp,8)]+= tmp258 + tmp266 + w18*(-A_2112[4] - A_2211[2]) + tmp421 + w4*(A_2112[3] + A_2211[5]);
                                    EM_S[INDEX4(2,1,6,1,numEq,numComp,8)]+= tmp470 + tmp474 + tmp87;
                                    EM_S[INDEX4(2,1,6,2,numEq,numComp,8)]+= tmp243 + tmp246 + tmp279 + w18*(-A_2112[6] + A_2211[2]) + w4*(A_2112[1] - A_2211[5]);
                                    EM_S[INDEX4(2,1,6,3,numEq,numComp,8)]+= tmp182 + tmp185 + tmp24;
                                    EM_S[INDEX4(2,1,6,4,numEq,numComp,8)]+= tmp42 + w18*(A_2112[4] - A_2211[6]) + w4*(-A_2112[3] + A_2211[1]) + tmp68 + tmp70;
                                    EM_S[INDEX4(2,1,6,5,numEq,numComp,8)]+= tmp166 + tmp460 + tmp461;
                                    EM_S[INDEX4(2,1,6,6,numEq,numComp,8)]+= tmp355 + tmp383 + tmp386 + w18*(A_2112[6] + A_2211[6]) + w4*(-A_2112[1] - A_2211[1]);
                                    EM_S[INDEX4(2,1,6,7,numEq,numComp,8)]+= tmp161 + tmp207 + tmp211;
                                    EM_S[INDEX4(2,1,7,0,numEq,numComp,8)]+= tmp470 + tmp474 + tmp87;
                                    EM_S[INDEX4(2,1,7,1,numEq,numComp,8)]+= tmp225 + tmp232 + w18*(-A_2112[5] - A_2211[3]) + tmp398 + w4*(A_2112[2] + A_2211[4]);
                                    EM_S[INDEX4(2,1,7,2,numEq,numComp,8)]+= tmp182 + tmp185 + tmp24;
                                    EM_S[INDEX4(2,1,7,3,numEq,numComp,8)]+= w18*(-A_2112[7] + A_2211[3]) + tmp13 + tmp16 + tmp5 + w4*(A_2112[0] - A_2211[4]);
                                    EM_S[INDEX4(2,1,7,4,numEq,numComp,8)]+= tmp166 + tmp460 + tmp461;
                                    EM_S[INDEX4(2,1,7,5,numEq,numComp,8)]+= w18*(A_2112[5] - A_2211[7]) + tmp432 + w4*(-A_2112[2] + A_2211[0]) + tmp51 + tmp54;
                                    EM_S[INDEX4(2,1,7,6,numEq,numComp,8)]+= tmp161 + tmp207 + tmp211;
                                    EM_S[INDEX4(2,1,7,7,numEq,numComp,8)]+= tmp309 + tmp312 + w18*(A_2112[7] + A_2211[7]) + w4*(-A_2112[0] - A_2211[0]) + tmp317;
                                }
                                {
                                    const double tmp1 = w13*(A_2222[1] + A_2222[2] + A_2222[5] + A_2222[6]);
                                    const double tmp3 = w14*(A_2020[2] + A_2020[3] + A_2020[6] + A_2020[7]);
                                    const double tmp4 = w7*(A_2222[0] + A_2222[4]);
                                    const double tmp6 = w3*(A_2121[0] + A_2121[2] + A_2121[4] + A_2121[6]);
                                    const double tmp10 = w0*(A_2020[0] + A_2020[1] + A_2020[4] + A_2020[5]);
                                    const double tmp12 = w9*(A_2121[1] + A_2121[3] + A_2121[5] + A_2121[7]);
                                    const double tmp17 = w19*(A_2222[3] + A_2222[7]);
                                    const double tmp20 = w13*(-A_2222[0] - A_2222[1] - A_2222[2] - A_2222[3] - A_2222[4] - A_2222[5] - A_2222[6] - A_2222[7]);
                                    const double tmp22 = w14*(-A_2020[4] - A_2020[5] - A_2020[6] - A_2020[7]);
                                    const double tmp25 = w3*(-A_2121[0] - A_2121[1] - A_2121[2] - A_2121[3]);
                                    const double tmp28 = w0*(-A_2020[0] - A_2020[1] - A_2020[2] - A_2020[3]);
                                    const double tmp30 = w9*(-A_2121[4] - A_2121[5] - A_2121[6] - A_2121[7]);
                                    const double tmp39 = w14*(A_2020[0] + A_2020[1] + A_2020[2] + A_2020[3]);
                                    const double tmp40 = w26*(A_2121[4] + A_2121[6]);
                                    const double tmp41 = w0*(A_2020[4] + A_2020[5] + A_2020[6] + A_2020[7]);
                                    const double tmp43 = w22*(A_2121[0] + A_2121[2] + A_2121[5] + A_2121[7]);
                                    const double tmp45 = w25*(A_2222[1] + A_2222[3] + A_2222[5] + A_2222[7]);
                                    const double tmp52 = w24*(A_2121[1] + A_2121[3]);
                                    const double tmp55 = w23*(A_2222[0] + A_2222[2] + A_2222[4] + A_2222[6]);
                                    const double tmp57 = w14*(A_2020[4] + A_2020[5] + A_2020[6] + A_2020[7]);
                                    const double tmp58 = w26*(A_2121[1] + A_2121[3]);
                                    const double tmp61 = w25*(A_2222[0] + A_2222[2] + A_2222[4] + A_2222[6]);
                                    const double tmp64 = w0*(A_2020[0] + A_2020[1] + A_2020[2] + A_2020[3]);
                                    const double tmp66 = w24*(A_2121[4] + A_2121[6]);
                                    const double tmp71 = w23*(A_2222[1] + A_2222[3] + A_2222[5] + A_2222[7]);
                                    const double tmp73 = w14*(-A_2020[0] - A_2020[1] - A_2020[2] - A_2020[3]);
                                    const double tmp74 = w0*(-A_2020[4] - A_2020[5] - A_2020[6] - A_2020[7]);
                                    const double tmp75 = w3*(-A_2121[4] - A_2121[5] - A_2121[6] - A_2121[7]);
                                    const double tmp80 = w9*(-A_2121[0] - A_2121[1] - A_2121[2] - A_2121[3]);
                                    const double tmp88 = w3*(A_2121[0] + A_2121[1] + A_2121[2] + A_2121[3]);
                                    const double tmp89 = w23*(A_2222[2] + A_2222[3] + A_2222[6] + A_2222[7]);
                                    const double tmp91 = w25*(A_2222[0] + A_2222[1] + A_2222[4] + A_2222[5]);
                                    const double tmp95 = w28*(A_2020[2] + A_2020[3]);
                                    const double tmp97 = w29*(A_2020[4] + A_2020[5]);
                                    const double tmp100 = w9*(A_2121[4] + A_2121[5] + A_2121[6] + A_2121[7]);
                                    const double tmp101 = w27*(A_2020[0] + A_2020[1] + A_2020[6] + A_2020[7]);
                                    const double tmp104 = w13*(A_2222[0] + A_2222[1] + A_2222[2] + A_2222[3] + A_2222[4] + A_2222[5] + A_2222[6] + A_2222[7]);
                                    const double tmp106 = w22*(A_2121[0] + A_2121[1] + A_2121[2] + A_2121[3] + A_2121[4] + A_2121[5] + A_2121[6] + A_2121[7]);
                                    const double tmp113 = w27*(A_2020[0] + A_2020[1] + A_2020[2] + A_2020[3] + A_2020[4] + A_2020[5] + A_2020[6] + A_2020[7]);
                                    const double tmp123 = w13*(A_2222[0] + A_2222[3] + A_2222[4] + A_2222[7]);
                                    const double tmp125 = w7*(A_2222[1] + A_2222[5]);
                                    const double tmp127 = w3*(A_2121[1] + A_2121[3] + A_2121[5] + A_2121[7]);
                                    const double tmp131 = w9*(A_2121[0] + A_2121[2] + A_2121[4] + A_2121[6]);
                                    const double tmp132 = w19*(A_2222[2] + A_2222[6]);
                                    const double tmp141 = w14*(A_2020[0] + A_2020[1] + A_2020[4] + A_2020[5]);
                                    const double tmp142 = w7*(A_2222[2] + A_2222[6]);
                                    const double tmp146 = w0*(A_2020[2] + A_2020[3] + A_2020[6] + A_2020[7]);
                                    const double tmp149 = w19*(A_2222[1] + A_2222[5]);
                                    const double tmp175 = w14*(-A_2020[2] - A_2020[3] - A_2020[6] - A_2020[7]);
                                    const double tmp176 = w22*(-A_2121[0] - A_2121[1] - A_2121[2] - A_2121[3] - A_2121[4] - A_2121[5] - A_2121[6] - A_2121[7]);
                                    const double tmp178 = w25*(-A_2222[2] - A_2222[3] - A_2222[6] - A_2222[7]);
                                    const double tmp180 = w0*(-A_2020[0] - A_2020[1] - A_2020[4] - A_2020[5]);
                                    const double tmp187 = w23*(-A_2222[0] - A_2222[1] - A_2222[4] - A_2222[5]);
                                    const double tmp189 = w7*(A_2222[3] + A_2222[7]);
                                    const double tmp193 = w19*(A_2222[0] + A_2222[4]);
                                    const double tmp201 = w27*(A_2020[2] + A_2020[3] + A_2020[4] + A_2020[5]);
                                    const double tmp204 = w23*(A_2222[0] + A_2222[1] + A_2222[4] + A_2222[5]);
                                    const double tmp205 = w25*(A_2222[2] + A_2222[3] + A_2222[6] + A_2222[7]);
                                    const double tmp208 = w28*(A_2020[0] + A_2020[1]);
                                    const double tmp209 = w29*(A_2020[6] + A_2020[7]);
                                    const double tmp214 = w13*(-A_2222[1] - A_2222[2] - A_2222[5] - A_2222[6]);
                                    const double tmp215 = w22*(-A_2121[0] - A_2121[2] - A_2121[5] - A_2121[7]);
                                    const double tmp217 = w27*(-A_2020[0] - A_2020[1] - A_2020[6] - A_2020[7]);
                                    const double tmp221 = w26*(-A_2121[4] - A_2121[6]);
                                    const double tmp226 = w7*(-A_2222[0] - A_2222[4]);
                                    const double tmp227 = w24*(-A_2121[1] - A_2121[3]);
                                    const double tmp228 = w19*(-A_2222[3] - A_2222[7]);
                                    const double tmp231 = w28*(-A_2020[4] - A_2020[5]);
                                    const double tmp233 = w29*(-A_2020[2] - A_2020[3]);
                                    const double tmp236 = w26*(A_2121[5] + A_2121[7]);
                                    const double tmp238 = w22*(A_2121[1] + A_2121[3] + A_2121[4] + A_2121[6]);
                                    const double tmp244 = w24*(A_2121[0] + A_2121[2]);
                                    const double tmp255 = w26*(-A_2121[1] - A_2121[3]);
                                    const double tmp259 = w7*(-A_2222[3] - A_2222[7]);
                                    const double tmp261 = w24*(-A_2121[4] - A_2121[6]);
                                    const double tmp262 = w19*(-A_2222[0] - A_2222[4]);
                                    const double tmp265 = w28*(-A_2020[2] - A_2020[3]);
                                    const double tmp268 = w29*(-A_2020[4] - A_2020[5]);
                                    const double tmp288 = w13*(-A_2222[0] - A_2222[3] - A_2222[4] - A_2222[7]);
                                    const double tmp289 = w22*(-A_2121[1] - A_2121[3] - A_2121[4] - A_2121[6]);
                                    const double tmp294 = w26*(-A_2121[5] - A_2121[7]);
                                    const double tmp298 = w7*(-A_2222[1] - A_2222[5]);
                                    const double tmp299 = w24*(-A_2121[0] - A_2121[2]);
                                    const double tmp300 = w19*(-A_2222[2] - A_2222[6]);
                                    const double tmp304 = w27*(-A_2020[2] - A_2020[3] - A_2020[4] - A_2020[5]);
                                    const double tmp308 = w26*(-A_2121[0] - A_2121[2]);
                                    const double tmp313 = w24*(-A_2121[5] - A_2121[7]);
                                    const double tmp316 = w28*(-A_2020[0] - A_2020[1]);
                                    const double tmp318 = w29*(-A_2020[6] - A_2020[7]);
                                    const double tmp320 = w26*(A_2121[0] + A_2121[2]);
                                    const double tmp325 = w24*(A_2121[5] + A_2121[7]);
                                    const double tmp329 = w3*(-A_2121[0] - A_2121[2] - A_2121[4] - A_2121[6]);
                                    const double tmp332 = w25*(-A_2222[1] - A_2222[3] - A_2222[5] - A_2222[7]);
                                    const double tmp335 = w9*(-A_2121[1] - A_2121[3] - A_2121[5] - A_2121[7]);
                                    const double tmp337 = w27*(-A_2020[0] - A_2020[1] - A_2020[2] - A_2020[3] - A_2020[4] - A_2020[5] - A_2020[6] - A_2020[7]);
                                    const double tmp338 = w23*(-A_2222[0] - A_2222[2] - A_2222[4] - A_2222[6]);
                                    const double tmp339 = w14*(-A_2020[0] - A_2020[1] - A_2020[4] - A_2020[5]);
                                    const double tmp340 = w23*(-A_2222[2] - A_2222[3] - A_2222[6] - A_2222[7]);
                                    const double tmp342 = w25*(-A_2222[0] - A_2222[1] - A_2222[4] - A_2222[5]);
                                    const double tmp344 = w0*(-A_2020[2] - A_2020[3] - A_2020[6] - A_2020[7]);
                                    const double tmp358 = w7*(-A_2222[2] - A_2222[6]);
                                    const double tmp359 = w19*(-A_2222[1] - A_2222[5]);
                                    const double tmp362 = w28*(-A_2020[6] - A_2020[7]);
                                    const double tmp363 = w29*(-A_2020[0] - A_2020[1]);
                                    const double tmp371 = w3*(A_2121[4] + A_2121[5] + A_2121[6] + A_2121[7]);
                                    const double tmp374 = w9*(A_2121[0] + A_2121[1] + A_2121[2] + A_2121[3]);
                                    const double tmp375 = w29*(A_2020[2] + A_2020[3]);
                                    const double tmp377 = w28*(A_2020[4] + A_2020[5]);
                                    const double tmp422 = w3*(-A_2121[1] - A_2121[3] - A_2121[5] - A_2121[7]);
                                    const double tmp424 = w25*(-A_2222[0] - A_2222[2] - A_2222[4] - A_2222[6]);
                                    const double tmp428 = w9*(-A_2121[0] - A_2121[2] - A_2121[4] - A_2121[6]);
                                    const double tmp430 = w23*(-A_2222[1] - A_2222[3] - A_2222[5] - A_2222[7]);
                                    const double tmp455 = w29*(A_2020[0] + A_2020[1]);
                                    const double tmp456 = w28*(A_2020[6] + A_2020[7]);
                                    EM_S[INDEX4(2,2,0,0,numEq,numComp,8)]+= tmp214 + tmp259 + tmp262 + tmp289 + tmp294 + tmp299 + tmp304 + tmp362 + tmp363;
                                    EM_S[INDEX4(2,2,0,1,numEq,numComp,8)]+= tmp201 + tmp371 + tmp374 + tmp455 + tmp456 + tmp89 + tmp91;
                                    EM_S[INDEX4(2,2,0,2,numEq,numComp,8)]+= tmp236 + tmp238 + tmp244 + tmp39 + tmp41 + tmp61 + tmp71;
                                    EM_S[INDEX4(2,2,0,3,numEq,numComp,8)]+= tmp20 + tmp73 + tmp74 + tmp75 + tmp80;
                                    EM_S[INDEX4(2,2,0,4,numEq,numComp,8)]+= tmp1 + tmp127 + tmp131 + tmp141 + tmp146 + tmp189 + tmp193;
                                    EM_S[INDEX4(2,2,0,5,numEq,numComp,8)]+= tmp176 + tmp339 + tmp340 + tmp342 + tmp344;
                                    EM_S[INDEX4(2,2,0,6,numEq,numComp,8)]+= tmp337 + tmp422 + tmp424 + tmp428 + tmp430;
                                    EM_S[INDEX4(2,2,0,7,numEq,numComp,8)]+= tmp104 + tmp106 + tmp113;
                                    EM_S[INDEX4(2,2,1,0,numEq,numComp,8)]+= tmp201 + tmp371 + tmp374 + tmp455 + tmp456 + tmp89 + tmp91;
                                    EM_S[INDEX4(2,2,1,1,numEq,numComp,8)]+= tmp215 + tmp221 + tmp227 + tmp288 + tmp304 + tmp358 + tmp359 + tmp362 + tmp363;
                                    EM_S[INDEX4(2,2,1,2,numEq,numComp,8)]+= tmp20 + tmp73 + tmp74 + tmp75 + tmp80;
                                    EM_S[INDEX4(2,2,1,3,numEq,numComp,8)]+= tmp39 + tmp40 + tmp41 + tmp43 + tmp45 + tmp52 + tmp55;
                                    EM_S[INDEX4(2,2,1,4,numEq,numComp,8)]+= tmp176 + tmp339 + tmp340 + tmp342 + tmp344;
                                    EM_S[INDEX4(2,2,1,5,numEq,numComp,8)]+= tmp12 + tmp123 + tmp141 + tmp142 + tmp146 + tmp149 + tmp6;
                                    EM_S[INDEX4(2,2,1,6,numEq,numComp,8)]+= tmp104 + tmp106 + tmp113;
                                    EM_S[INDEX4(2,2,1,7,numEq,numComp,8)]+= tmp329 + tmp332 + tmp335 + tmp337 + tmp338;
                                    EM_S[INDEX4(2,2,2,0,numEq,numComp,8)]+= tmp236 + tmp238 + tmp244 + tmp39 + tmp41 + tmp61 + tmp71;
                                    EM_S[INDEX4(2,2,2,1,numEq,numComp,8)]+= tmp20 + tmp73 + tmp74 + tmp75 + tmp80;
                                    EM_S[INDEX4(2,2,2,2,numEq,numComp,8)]+= tmp217 + tmp231 + tmp233 + tmp288 + tmp289 + tmp294 + tmp298 + tmp299 + tmp300;
                                    EM_S[INDEX4(2,2,2,3,numEq,numComp,8)]+= tmp101 + tmp204 + tmp205 + tmp371 + tmp374 + tmp375 + tmp377;
                                    EM_S[INDEX4(2,2,2,4,numEq,numComp,8)]+= tmp337 + tmp422 + tmp424 + tmp428 + tmp430;
                                    EM_S[INDEX4(2,2,2,5,numEq,numComp,8)]+= tmp104 + tmp106 + tmp113;
                                    EM_S[INDEX4(2,2,2,6,numEq,numComp,8)]+= tmp10 + tmp123 + tmp125 + tmp127 + tmp131 + tmp132 + tmp3;
                                    EM_S[INDEX4(2,2,2,7,numEq,numComp,8)]+= tmp175 + tmp176 + tmp178 + tmp180 + tmp187;
                                    EM_S[INDEX4(2,2,3,0,numEq,numComp,8)]+= tmp20 + tmp73 + tmp74 + tmp75 + tmp80;
                                    EM_S[INDEX4(2,2,3,1,numEq,numComp,8)]+= tmp39 + tmp40 + tmp41 + tmp43 + tmp45 + tmp52 + tmp55;
                                    EM_S[INDEX4(2,2,3,2,numEq,numComp,8)]+= tmp101 + tmp204 + tmp205 + tmp371 + tmp374 + tmp375 + tmp377;
                                    EM_S[INDEX4(2,2,3,3,numEq,numComp,8)]+= tmp214 + tmp215 + tmp217 + tmp221 + tmp226 + tmp227 + tmp228 + tmp231 + tmp233;
                                    EM_S[INDEX4(2,2,3,4,numEq,numComp,8)]+= tmp104 + tmp106 + tmp113;
                                    EM_S[INDEX4(2,2,3,5,numEq,numComp,8)]+= tmp329 + tmp332 + tmp335 + tmp337 + tmp338;
                                    EM_S[INDEX4(2,2,3,6,numEq,numComp,8)]+= tmp175 + tmp176 + tmp178 + tmp180 + tmp187;
                                    EM_S[INDEX4(2,2,3,7,numEq,numComp,8)]+= tmp1 + tmp10 + tmp12 + tmp17 + tmp3 + tmp4 + tmp6;
                                    EM_S[INDEX4(2,2,4,0,numEq,numComp,8)]+= tmp1 + tmp127 + tmp131 + tmp141 + tmp146 + tmp189 + tmp193;
                                    EM_S[INDEX4(2,2,4,1,numEq,numComp,8)]+= tmp176 + tmp339 + tmp340 + tmp342 + tmp344;
                                    EM_S[INDEX4(2,2,4,2,numEq,numComp,8)]+= tmp337 + tmp422 + tmp424 + tmp428 + tmp430;
                                    EM_S[INDEX4(2,2,4,3,numEq,numComp,8)]+= tmp104 + tmp106 + tmp113;
                                    EM_S[INDEX4(2,2,4,4,numEq,numComp,8)]+= tmp214 + tmp215 + tmp217 + tmp255 + tmp259 + tmp261 + tmp262 + tmp265 + tmp268;
                                    EM_S[INDEX4(2,2,4,5,numEq,numComp,8)]+= tmp100 + tmp101 + tmp88 + tmp89 + tmp91 + tmp95 + tmp97;
                                    EM_S[INDEX4(2,2,4,6,numEq,numComp,8)]+= tmp43 + tmp57 + tmp58 + tmp61 + tmp64 + tmp66 + tmp71;
                                    EM_S[INDEX4(2,2,4,7,numEq,numComp,8)]+= tmp20 + tmp22 + tmp25 + tmp28 + tmp30;
                                    EM_S[INDEX4(2,2,5,0,numEq,numComp,8)]+= tmp176 + tmp339 + tmp340 + tmp342 + tmp344;
                                    EM_S[INDEX4(2,2,5,1,numEq,numComp,8)]+= tmp12 + tmp123 + tmp141 + tmp142 + tmp146 + tmp149 + tmp6;
                                    EM_S[INDEX4(2,2,5,2,numEq,numComp,8)]+= tmp104 + tmp106 + tmp113;
                                    EM_S[INDEX4(2,2,5,3,numEq,numComp,8)]+= tmp329 + tmp332 + tmp335 + tmp337 + tmp338;
                                    EM_S[INDEX4(2,2,5,4,numEq,numComp,8)]+= tmp100 + tmp101 + tmp88 + tmp89 + tmp91 + tmp95 + tmp97;
                                    EM_S[INDEX4(2,2,5,5,numEq,numComp,8)]+= tmp217 + tmp265 + tmp268 + tmp288 + tmp289 + tmp308 + tmp313 + tmp358 + tmp359;
                                    EM_S[INDEX4(2,2,5,6,numEq,numComp,8)]+= tmp20 + tmp22 + tmp25 + tmp28 + tmp30;
                                    EM_S[INDEX4(2,2,5,7,numEq,numComp,8)]+= tmp238 + tmp320 + tmp325 + tmp45 + tmp55 + tmp57 + tmp64;
                                    EM_S[INDEX4(2,2,6,0,numEq,numComp,8)]+= tmp337 + tmp422 + tmp424 + tmp428 + tmp430;
                                    EM_S[INDEX4(2,2,6,1,numEq,numComp,8)]+= tmp104 + tmp106 + tmp113;
                                    EM_S[INDEX4(2,2,6,2,numEq,numComp,8)]+= tmp10 + tmp123 + tmp125 + tmp127 + tmp131 + tmp132 + tmp3;
                                    EM_S[INDEX4(2,2,6,3,numEq,numComp,8)]+= tmp175 + tmp176 + tmp178 + tmp180 + tmp187;
                                    EM_S[INDEX4(2,2,6,4,numEq,numComp,8)]+= tmp43 + tmp57 + tmp58 + tmp61 + tmp64 + tmp66 + tmp71;
                                    EM_S[INDEX4(2,2,6,5,numEq,numComp,8)]+= tmp20 + tmp22 + tmp25 + tmp28 + tmp30;
                                    EM_S[INDEX4(2,2,6,6,numEq,numComp,8)]+= tmp215 + tmp255 + tmp261 + tmp288 + tmp298 + tmp300 + tmp304 + tmp316 + tmp318;
                                    EM_S[INDEX4(2,2,6,7,numEq,numComp,8)]+= tmp100 + tmp201 + tmp204 + tmp205 + tmp208 + tmp209 + tmp88;
                                    EM_S[INDEX4(2,2,7,0,numEq,numComp,8)]+= tmp104 + tmp106 + tmp113;
                                    EM_S[INDEX4(2,2,7,1,numEq,numComp,8)]+= tmp329 + tmp332 + tmp335 + tmp337 + tmp338;
                                    EM_S[INDEX4(2,2,7,2,numEq,numComp,8)]+= tmp175 + tmp176 + tmp178 + tmp180 + tmp187;
                                    EM_S[INDEX4(2,2,7,3,numEq,numComp,8)]+= tmp1 + tmp10 + tmp12 + tmp17 + tmp3 + tmp4 + tmp6;
                                    EM_S[INDEX4(2,2,7,4,numEq,numComp,8)]+= tmp20 + tmp22 + tmp25 + tmp28 + tmp30;
                                    EM_S[INDEX4(2,2,7,5,numEq,numComp,8)]+= tmp238 + tmp320 + tmp325 + tmp45 + tmp55 + tmp57 + tmp64;
                                    EM_S[INDEX4(2,2,7,6,numEq,numComp,8)]+= tmp100 + tmp201 + tmp204 + tmp205 + tmp208 + tmp209 + tmp88;
                                    EM_S[INDEX4(2,2,7,7,numEq,numComp,8)]+= tmp214 + tmp226 + tmp228 + tmp289 + tmp304 + tmp308 + tmp313 + tmp316 + tmp318;
                                }
                            } else { // constant data
                                double Aw0000 = 0;
                                double Aw0011 = 0;
                                double Aw0022 = 0;
                                double Aw0101 = 0;
                                double Aw0110 = 0;
                                double Aw0202 = 0;
                                double Aw0220 = 0;
                                double Aw1001 = 0;
                                double Aw1010 = 0;
                                double Aw1100 = 0;
                                double Aw1111 = 0;
                                double Aw1122 = 0;
                                double Aw1212 = 0;
                                double Aw1221 = 0;
                                double Aw2002 = 0;
                                double Aw2020 = 0;
                                double Aw2112 = 0;
                                double Aw2121 = 0;
                                double Aw2200 = 0;
                                double Aw2211 = 0;
                                double Aw2222 = 0;
                                if (!mu.isEmpty()) {
                                    const double *mu_p = mu.getSampleDataRO(e);
                                    Aw0000 += 2*mu_p[0];
                                    Aw0110 += mu_p[0];
                                    Aw0101 += mu_p[0];
                                    Aw0220 += mu_p[0];
                                    Aw0202 += mu_p[0];
                                    Aw1001 += mu_p[0];
                                    Aw1010 += mu_p[0];
                                    Aw1111 += 2*mu_p[0];
                                    Aw1221 += mu_p[0];
                                    Aw1212 += mu_p[0];
                                    Aw2002 += mu_p[0];
                                    Aw2020 += mu_p[0];
                                    Aw2112 += mu_p[0];
                                    Aw2121 += mu_p[0];
                                    Aw2222 += 2*mu_p[0];
                                }
                                if (!lambda.isEmpty()) {
                                    const double *lambda_p = lambda.getSampleDataRO(e);
                                    Aw0000 += lambda_p[0];
                                    Aw0011 += lambda_p[0];
                                    Aw0022 += lambda_p[0];
                                    Aw1100 += lambda_p[0];
                                    Aw1111 += lambda_p[0];
                                    Aw1122 += lambda_p[0];
                                    Aw2200 += lambda_p[0];
                                    Aw2211 += lambda_p[0];
                                    Aw2222 += lambda_p[0];
                                }
                                Aw0000 *= 8*w27;
                                Aw0101 *= 8*w22;
                                Aw0202 *= 8*w13;
                                Aw0011 *= 12*w8;
                                Aw0110 *= 12*w8;
                                Aw0022 *= 12*w11;
                                Aw0220 *= 12*w11;
                                Aw1001 *= 12*w8;
                                Aw1100 *= 12*w8;
                                Aw1010 *= 8*w27;
                                Aw1111 *= 8*w22;
                                Aw1212 *= 8*w13;
                                Aw1122 *= 12*w10;
                                Aw1221 *= 12*w10;
                                Aw2002 *= 12*w11;
                                Aw2200 *= 12*w11;
                                Aw2112 *= 12*w10;
                                Aw2211 *= 12*w10;
                                Aw2020 *= 8*w27;
                                Aw2121 *= 8*w22;
                                Aw2222 *= 8*w13;
                                {
                                    const double tmp6 = 4*Aw0000;
                                    const double tmp7 = 2*Aw0000;
                                    const double tmp8 = 4*Aw0101;
                                    const double tmp9 = 2*Aw0101;
                                    const double tmp10 = 4*Aw0202;
                                    const double tmp11 = 2*Aw0202;
                                    EM_S[INDEX4(0,0,0,0,numEq,numComp,8)]+= -tmp6 - tmp8 - tmp10;
                                    EM_S[INDEX4(0,0,0,1,numEq,numComp,8)]+= tmp6 - tmp9 - tmp11;
                                    EM_S[INDEX4(0,0,0,2,numEq,numComp,8)]+= -tmp7 + tmp8 - tmp11;
                                    EM_S[INDEX4(0,0,0,3,numEq,numComp,8)]+= tmp7 + tmp9 - Aw0202;
                                    EM_S[INDEX4(0,0,0,4,numEq,numComp,8)]+= -tmp7 - tmp9 + tmp10;
                                    EM_S[INDEX4(0,0,0,5,numEq,numComp,8)]+= tmp7 - Aw0101 + tmp11;
                                    EM_S[INDEX4(0,0,0,6,numEq,numComp,8)]+= -Aw0000 + tmp9 + tmp11;
                                    EM_S[INDEX4(0,0,0,7,numEq,numComp,8)]+= Aw0000 + Aw0101 + Aw0202;
                                    EM_S[INDEX4(0,0,1,0,numEq,numComp,8)]+= tmp6 - tmp9 - tmp11;
                                    EM_S[INDEX4(0,0,1,1,numEq,numComp,8)]+= -tmp6 - tmp8 - tmp10;
                                    EM_S[INDEX4(0,0,1,2,numEq,numComp,8)]+= tmp7 + tmp9 - Aw0202;
                                    EM_S[INDEX4(0,0,1,3,numEq,numComp,8)]+= -tmp7 + tmp8 - tmp11;
                                    EM_S[INDEX4(0,0,1,4,numEq,numComp,8)]+= tmp7 - Aw0101 + tmp11;
                                    EM_S[INDEX4(0,0,1,5,numEq,numComp,8)]+= -tmp7 - tmp9 + tmp10;
                                    EM_S[INDEX4(0,0,1,6,numEq,numComp,8)]+= Aw0000 + Aw0101 + Aw0202;
                                    EM_S[INDEX4(0,0,1,7,numEq,numComp,8)]+= -Aw0000 + tmp9 + tmp11;
                                    EM_S[INDEX4(0,0,2,0,numEq,numComp,8)]+= -tmp7 + tmp8 - tmp11;
                                    EM_S[INDEX4(0,0,2,1,numEq,numComp,8)]+= tmp7 + tmp9 - Aw0202;
                                    EM_S[INDEX4(0,0,2,2,numEq,numComp,8)]+= -tmp6 - tmp8 - tmp10;
                                    EM_S[INDEX4(0,0,2,3,numEq,numComp,8)]+= tmp6 - tmp9 - tmp11;
                                    EM_S[INDEX4(0,0,2,4,numEq,numComp,8)]+= -Aw0000 + tmp9 + tmp11;
                                    EM_S[INDEX4(0,0,2,5,numEq,numComp,8)]+= Aw0000 + Aw0101 + Aw0202;
                                    EM_S[INDEX4(0,0,2,6,numEq,numComp,8)]+= -tmp7 - tmp9 + tmp10;
                                    EM_S[INDEX4(0,0,2,7,numEq,numComp,8)]+= tmp7 - Aw0101 + tmp11;
                                    EM_S[INDEX4(0,0,3,0,numEq,numComp,8)]+= tmp7 + tmp9 - Aw0202;
                                    EM_S[INDEX4(0,0,3,1,numEq,numComp,8)]+= -tmp7 + tmp8 - tmp11;
                                    EM_S[INDEX4(0,0,3,2,numEq,numComp,8)]+= tmp6 - tmp9 - tmp11;
                                    EM_S[INDEX4(0,0,3,3,numEq,numComp,8)]+= -tmp6 - tmp8 - tmp10;
                                    EM_S[INDEX4(0,0,3,4,numEq,numComp,8)]+= Aw0000 + Aw0101 + Aw0202;
                                    EM_S[INDEX4(0,0,3,5,numEq,numComp,8)]+= -Aw0000 + tmp9 + tmp11;
                                    EM_S[INDEX4(0,0,3,6,numEq,numComp,8)]+= tmp7 - Aw0101 + tmp11;
                                    EM_S[INDEX4(0,0,3,7,numEq,numComp,8)]+= -tmp7 - tmp9 + tmp10;
                                    EM_S[INDEX4(0,0,4,0,numEq,numComp,8)]+= -tmp7 - tmp9 + tmp10;
                                    EM_S[INDEX4(0,0,4,1,numEq,numComp,8)]+= tmp7 - Aw0101 + tmp11;
                                    EM_S[INDEX4(0,0,4,2,numEq,numComp,8)]+= -Aw0000 + tmp9 + tmp11;
                                    EM_S[INDEX4(0,0,4,3,numEq,numComp,8)]+= Aw0000 + Aw0101 + Aw0202;
                                    EM_S[INDEX4(0,0,4,4,numEq,numComp,8)]+= -tmp6 - tmp8 - tmp10;
                                    EM_S[INDEX4(0,0,4,5,numEq,numComp,8)]+= tmp6 - tmp9 - tmp11;
                                    EM_S[INDEX4(0,0,4,6,numEq,numComp,8)]+= -tmp7 + tmp8 - tmp11;
                                    EM_S[INDEX4(0,0,4,7,numEq,numComp,8)]+= tmp7 + tmp9 - Aw0202;
                                    EM_S[INDEX4(0,0,5,0,numEq,numComp,8)]+= tmp7 - Aw0101 + tmp11;
                                    EM_S[INDEX4(0,0,5,1,numEq,numComp,8)]+= -tmp7 - tmp9 + tmp10;
                                    EM_S[INDEX4(0,0,5,2,numEq,numComp,8)]+= Aw0000 + Aw0101 + Aw0202;
                                    EM_S[INDEX4(0,0,5,3,numEq,numComp,8)]+= -Aw0000 + tmp9 + tmp11;
                                    EM_S[INDEX4(0,0,5,4,numEq,numComp,8)]+= tmp6 - tmp9 - tmp11;
                                    EM_S[INDEX4(0,0,5,5,numEq,numComp,8)]+= -tmp6 - tmp8 - tmp10;
                                    EM_S[INDEX4(0,0,5,6,numEq,numComp,8)]+= tmp7 + tmp9 - Aw0202;
                                    EM_S[INDEX4(0,0,5,7,numEq,numComp,8)]+= -tmp7 + tmp8 - tmp11;
                                    EM_S[INDEX4(0,0,6,0,numEq,numComp,8)]+= -Aw0000 + tmp9 + tmp11;
                                    EM_S[INDEX4(0,0,6,1,numEq,numComp,8)]+= Aw0000 + Aw0101 + Aw0202;
                                    EM_S[INDEX4(0,0,6,2,numEq,numComp,8)]+= -tmp7 - tmp9 + tmp10;
                                    EM_S[INDEX4(0,0,6,3,numEq,numComp,8)]+= tmp7 - Aw0101 + tmp11;
                                    EM_S[INDEX4(0,0,6,4,numEq,numComp,8)]+= -tmp7 + tmp8 - tmp11;
                                    EM_S[INDEX4(0,0,6,5,numEq,numComp,8)]+= tmp7 + tmp9 - Aw0202;
                                    EM_S[INDEX4(0,0,6,6,numEq,numComp,8)]+= -tmp6 - tmp8 - tmp10;
                                    EM_S[INDEX4(0,0,6,7,numEq,numComp,8)]+= tmp6 - tmp9 - tmp11;
                                    EM_S[INDEX4(0,0,7,0,numEq,numComp,8)]+= Aw0000 + Aw0101 + Aw0202;
                                    EM_S[INDEX4(0,0,7,1,numEq,numComp,8)]+= -Aw0000 + tmp9 + tmp11;
                                    EM_S[INDEX4(0,0,7,2,numEq,numComp,8)]+= tmp7 - Aw0101 + tmp11;
                                    EM_S[INDEX4(0,0,7,3,numEq,numComp,8)]+= -tmp7 - tmp9 + tmp10;
                                    EM_S[INDEX4(0,0,7,4,numEq,numComp,8)]+= tmp7 + tmp9 - Aw0202;
                                    EM_S[INDEX4(0,0,7,5,numEq,numComp,8)]+= -tmp7 + tmp8 - tmp11;
                                    EM_S[INDEX4(0,0,7,6,numEq,numComp,8)]+= tmp6 - tmp9 - tmp11;
                                    EM_S[INDEX4(0,0,7,7,numEq,numComp,8)]+= -tmp6 - tmp8 - tmp10;
                                }
                                {
                                    const double tmp0 = Aw0011 + Aw0110;
                                    const double tmp1 = Aw0011 - Aw0110;
                                    const double tmp15 = 2*tmp0;
                                    const double tmp16 = 2*tmp1;
                                    EM_S[INDEX4(0,1,0,0,numEq,numComp,8)]+= tmp15;
                                    EM_S[INDEX4(0,1,0,1,numEq,numComp,8)]+= tmp16;
                                    EM_S[INDEX4(0,1,0,2,numEq,numComp,8)]+=-tmp16;
                                    EM_S[INDEX4(0,1,0,3,numEq,numComp,8)]+=-tmp15;
                                    EM_S[INDEX4(0,1,0,4,numEq,numComp,8)]+= tmp0;
                                    EM_S[INDEX4(0,1,0,5,numEq,numComp,8)]+= tmp1;
                                    EM_S[INDEX4(0,1,0,6,numEq,numComp,8)]+=-tmp1;
                                    EM_S[INDEX4(0,1,0,7,numEq,numComp,8)]+=-tmp0;
                                    EM_S[INDEX4(0,1,1,0,numEq,numComp,8)]+=-tmp16;
                                    EM_S[INDEX4(0,1,1,1,numEq,numComp,8)]+=-tmp15;
                                    EM_S[INDEX4(0,1,1,2,numEq,numComp,8)]+= tmp15;
                                    EM_S[INDEX4(0,1,1,3,numEq,numComp,8)]+= tmp16;
                                    EM_S[INDEX4(0,1,1,4,numEq,numComp,8)]+=-tmp1;
                                    EM_S[INDEX4(0,1,1,5,numEq,numComp,8)]+=-tmp0;
                                    EM_S[INDEX4(0,1,1,6,numEq,numComp,8)]+= tmp0;
                                    EM_S[INDEX4(0,1,1,7,numEq,numComp,8)]+= tmp1;
                                    EM_S[INDEX4(0,1,2,0,numEq,numComp,8)]+= tmp16;
                                    EM_S[INDEX4(0,1,2,1,numEq,numComp,8)]+= tmp15;
                                    EM_S[INDEX4(0,1,2,2,numEq,numComp,8)]+=-tmp15;
                                    EM_S[INDEX4(0,1,2,3,numEq,numComp,8)]+=-tmp16;
                                    EM_S[INDEX4(0,1,2,4,numEq,numComp,8)]+= tmp1;
                                    EM_S[INDEX4(0,1,2,5,numEq,numComp,8)]+= tmp0;
                                    EM_S[INDEX4(0,1,2,6,numEq,numComp,8)]+=-tmp0;
                                    EM_S[INDEX4(0,1,2,7,numEq,numComp,8)]+=-tmp1;
                                    EM_S[INDEX4(0,1,3,0,numEq,numComp,8)]+=-tmp15;
                                    EM_S[INDEX4(0,1,3,1,numEq,numComp,8)]+=-tmp16;
                                    EM_S[INDEX4(0,1,3,2,numEq,numComp,8)]+= tmp16;
                                    EM_S[INDEX4(0,1,3,3,numEq,numComp,8)]+= tmp15;
                                    EM_S[INDEX4(0,1,3,4,numEq,numComp,8)]+=-tmp0;
                                    EM_S[INDEX4(0,1,3,5,numEq,numComp,8)]+=-tmp1;
                                    EM_S[INDEX4(0,1,3,6,numEq,numComp,8)]+= tmp1;
                                    EM_S[INDEX4(0,1,3,7,numEq,numComp,8)]+= tmp0;
                                    EM_S[INDEX4(0,1,4,0,numEq,numComp,8)]+= tmp0;
                                    EM_S[INDEX4(0,1,4,1,numEq,numComp,8)]+= tmp1;
                                    EM_S[INDEX4(0,1,4,2,numEq,numComp,8)]+=-tmp1;
                                    EM_S[INDEX4(0,1,4,3,numEq,numComp,8)]+=-tmp0;
                                    EM_S[INDEX4(0,1,4,4,numEq,numComp,8)]+= tmp15;
                                    EM_S[INDEX4(0,1,4,5,numEq,numComp,8)]+= tmp16;
                                    EM_S[INDEX4(0,1,4,6,numEq,numComp,8)]+=-tmp16;
                                    EM_S[INDEX4(0,1,4,7,numEq,numComp,8)]+=-tmp15;
                                    EM_S[INDEX4(0,1,5,0,numEq,numComp,8)]+=-tmp1;
                                    EM_S[INDEX4(0,1,5,1,numEq,numComp,8)]+=-tmp0;
                                    EM_S[INDEX4(0,1,5,2,numEq,numComp,8)]+= tmp0;
                                    EM_S[INDEX4(0,1,5,3,numEq,numComp,8)]+= tmp1;
                                    EM_S[INDEX4(0,1,5,4,numEq,numComp,8)]+=-tmp16;
                                    EM_S[INDEX4(0,1,5,5,numEq,numComp,8)]+=-tmp15;
                                    EM_S[INDEX4(0,1,5,6,numEq,numComp,8)]+= tmp15;
                                    EM_S[INDEX4(0,1,5,7,numEq,numComp,8)]+= tmp16;
                                    EM_S[INDEX4(0,1,6,0,numEq,numComp,8)]+= tmp1;
                                    EM_S[INDEX4(0,1,6,1,numEq,numComp,8)]+= tmp0;
                                    EM_S[INDEX4(0,1,6,2,numEq,numComp,8)]+=-tmp0;
                                    EM_S[INDEX4(0,1,6,3,numEq,numComp,8)]+=-tmp1;
                                    EM_S[INDEX4(0,1,6,4,numEq,numComp,8)]+= tmp16;
                                    EM_S[INDEX4(0,1,6,5,numEq,numComp,8)]+= tmp15;
                                    EM_S[INDEX4(0,1,6,6,numEq,numComp,8)]+=-tmp15;
                                    EM_S[INDEX4(0,1,6,7,numEq,numComp,8)]+=-tmp16;
                                    EM_S[INDEX4(0,1,7,0,numEq,numComp,8)]+=-tmp0;
                                    EM_S[INDEX4(0,1,7,1,numEq,numComp,8)]+=-tmp1;
                                    EM_S[INDEX4(0,1,7,2,numEq,numComp,8)]+= tmp1;
                                    EM_S[INDEX4(0,1,7,3,numEq,numComp,8)]+= tmp0;
                                    EM_S[INDEX4(0,1,7,4,numEq,numComp,8)]+=-tmp15;
                                    EM_S[INDEX4(0,1,7,5,numEq,numComp,8)]+=-tmp16;
                                    EM_S[INDEX4(0,1,7,6,numEq,numComp,8)]+= tmp16;
                                    EM_S[INDEX4(0,1,7,7,numEq,numComp,8)]+= tmp15;
                                }
                                {
                                    const double tmp2 = Aw0022 + Aw0220;
                                    const double tmp3 = Aw0022 - Aw0220;
                                    const double tmp17 = 2*tmp2;
                                    const double tmp18 = 2*tmp3;
                                    EM_S[INDEX4(0,2,0,0,numEq,numComp,8)]+= tmp17;
                                    EM_S[INDEX4(0,2,0,1,numEq,numComp,8)]+= tmp18;
                                    EM_S[INDEX4(0,2,0,2,numEq,numComp,8)]+= tmp2;
                                    EM_S[INDEX4(0,2,0,3,numEq,numComp,8)]+= tmp3;
                                    EM_S[INDEX4(0,2,0,4,numEq,numComp,8)]+=-tmp18;
                                    EM_S[INDEX4(0,2,0,5,numEq,numComp,8)]+=-tmp17;
                                    EM_S[INDEX4(0,2,0,6,numEq,numComp,8)]+=-tmp3;
                                    EM_S[INDEX4(0,2,0,7,numEq,numComp,8)]+=-tmp2;
                                    EM_S[INDEX4(0,2,1,0,numEq,numComp,8)]+=-tmp18;
                                    EM_S[INDEX4(0,2,1,1,numEq,numComp,8)]+=-tmp17;
                                    EM_S[INDEX4(0,2,1,2,numEq,numComp,8)]+=-tmp3;
                                    EM_S[INDEX4(0,2,1,3,numEq,numComp,8)]+=-tmp2;
                                    EM_S[INDEX4(0,2,1,4,numEq,numComp,8)]+= tmp17;
                                    EM_S[INDEX4(0,2,1,5,numEq,numComp,8)]+= tmp18;
                                    EM_S[INDEX4(0,2,1,6,numEq,numComp,8)]+= tmp2;
                                    EM_S[INDEX4(0,2,1,7,numEq,numComp,8)]+= tmp3;
                                    EM_S[INDEX4(0,2,2,0,numEq,numComp,8)]+= tmp2;
                                    EM_S[INDEX4(0,2,2,1,numEq,numComp,8)]+= tmp3;
                                    EM_S[INDEX4(0,2,2,2,numEq,numComp,8)]+= tmp17;
                                    EM_S[INDEX4(0,2,2,3,numEq,numComp,8)]+= tmp18;
                                    EM_S[INDEX4(0,2,2,4,numEq,numComp,8)]+=-tmp3;
                                    EM_S[INDEX4(0,2,2,5,numEq,numComp,8)]+=-tmp2;
                                    EM_S[INDEX4(0,2,2,6,numEq,numComp,8)]+=-tmp18;
                                    EM_S[INDEX4(0,2,2,7,numEq,numComp,8)]+=-tmp17;
                                    EM_S[INDEX4(0,2,3,0,numEq,numComp,8)]+=-tmp3;
                                    EM_S[INDEX4(0,2,3,1,numEq,numComp,8)]+=-tmp2;
                                    EM_S[INDEX4(0,2,3,2,numEq,numComp,8)]+=-tmp18;
                                    EM_S[INDEX4(0,2,3,3,numEq,numComp,8)]+=-tmp17;
                                    EM_S[INDEX4(0,2,3,4,numEq,numComp,8)]+= tmp2;
                                    EM_S[INDEX4(0,2,3,5,numEq,numComp,8)]+= tmp3;
                                    EM_S[INDEX4(0,2,3,6,numEq,numComp,8)]+= tmp17;
                                    EM_S[INDEX4(0,2,3,7,numEq,numComp,8)]+= tmp18;
                                    EM_S[INDEX4(0,2,4,0,numEq,numComp,8)]+= tmp18;
                                    EM_S[INDEX4(0,2,4,1,numEq,numComp,8)]+= tmp17;
                                    EM_S[INDEX4(0,2,4,2,numEq,numComp,8)]+= tmp3;
                                    EM_S[INDEX4(0,2,4,3,numEq,numComp,8)]+= tmp2;
                                    EM_S[INDEX4(0,2,4,4,numEq,numComp,8)]+=-tmp17;
                                    EM_S[INDEX4(0,2,4,5,numEq,numComp,8)]+=-tmp18;
                                    EM_S[INDEX4(0,2,4,6,numEq,numComp,8)]+=-tmp2;
                                    EM_S[INDEX4(0,2,4,7,numEq,numComp,8)]+=-tmp3;
                                    EM_S[INDEX4(0,2,5,0,numEq,numComp,8)]+=-tmp17;
                                    EM_S[INDEX4(0,2,5,1,numEq,numComp,8)]+=-tmp18;
                                    EM_S[INDEX4(0,2,5,2,numEq,numComp,8)]+=-tmp2;
                                    EM_S[INDEX4(0,2,5,3,numEq,numComp,8)]+=-tmp3;
                                    EM_S[INDEX4(0,2,5,4,numEq,numComp,8)]+= tmp18;
                                    EM_S[INDEX4(0,2,5,5,numEq,numComp,8)]+= tmp17;
                                    EM_S[INDEX4(0,2,5,6,numEq,numComp,8)]+= tmp3;
                                    EM_S[INDEX4(0,2,5,7,numEq,numComp,8)]+= tmp2;
                                    EM_S[INDEX4(0,2,6,0,numEq,numComp,8)]+= tmp3;
                                    EM_S[INDEX4(0,2,6,1,numEq,numComp,8)]+= tmp2;
                                    EM_S[INDEX4(0,2,6,2,numEq,numComp,8)]+= tmp18;
                                    EM_S[INDEX4(0,2,6,3,numEq,numComp,8)]+= tmp17;
                                    EM_S[INDEX4(0,2,6,4,numEq,numComp,8)]+=-tmp2;
                                    EM_S[INDEX4(0,2,6,5,numEq,numComp,8)]+=-tmp3;
                                    EM_S[INDEX4(0,2,6,6,numEq,numComp,8)]+=-tmp17;
                                    EM_S[INDEX4(0,2,6,7,numEq,numComp,8)]+=-tmp18;
                                    EM_S[INDEX4(0,2,7,0,numEq,numComp,8)]+=-tmp2;
                                    EM_S[INDEX4(0,2,7,1,numEq,numComp,8)]+=-tmp3;
                                    EM_S[INDEX4(0,2,7,2,numEq,numComp,8)]+=-tmp17;
                                    EM_S[INDEX4(0,2,7,3,numEq,numComp,8)]+=-tmp18;
                                    EM_S[INDEX4(0,2,7,4,numEq,numComp,8)]+= tmp3;
                                    EM_S[INDEX4(0,2,7,5,numEq,numComp,8)]+= tmp2;
                                    EM_S[INDEX4(0,2,7,6,numEq,numComp,8)]+= tmp18;
                                    EM_S[INDEX4(0,2,7,7,numEq,numComp,8)]+= tmp17;
                                }
                                {
                                    const double tmp0 = Aw1001 + Aw1100;
                                    const double tmp1 = Aw1001 - Aw1100;
                                    const double tmp15 = 2*tmp0;
                                    const double tmp16 = 2*tmp1;
                                    EM_S[INDEX4(1,0,0,0,numEq,numComp,8)]+= tmp15;
                                    EM_S[INDEX4(1,0,0,1,numEq,numComp,8)]+= tmp16;
                                    EM_S[INDEX4(1,0,0,2,numEq,numComp,8)]+=-tmp16;
                                    EM_S[INDEX4(1,0,0,3,numEq,numComp,8)]+=-tmp15;
                                    EM_S[INDEX4(1,0,0,4,numEq,numComp,8)]+= tmp0;
                                    EM_S[INDEX4(1,0,0,5,numEq,numComp,8)]+= tmp1;
                                    EM_S[INDEX4(1,0,0,6,numEq,numComp,8)]+=-tmp1;
                                    EM_S[INDEX4(1,0,0,7,numEq,numComp,8)]+=-tmp0;
                                    EM_S[INDEX4(1,0,1,0,numEq,numComp,8)]+=-tmp16;
                                    EM_S[INDEX4(1,0,1,1,numEq,numComp,8)]+=-tmp15;
                                    EM_S[INDEX4(1,0,1,2,numEq,numComp,8)]+= tmp15;
                                    EM_S[INDEX4(1,0,1,3,numEq,numComp,8)]+= tmp16;
                                    EM_S[INDEX4(1,0,1,4,numEq,numComp,8)]+=-tmp1;
                                    EM_S[INDEX4(1,0,1,5,numEq,numComp,8)]+=-tmp0;
                                    EM_S[INDEX4(1,0,1,6,numEq,numComp,8)]+= tmp0;
                                    EM_S[INDEX4(1,0,1,7,numEq,numComp,8)]+= tmp1;
                                    EM_S[INDEX4(1,0,2,0,numEq,numComp,8)]+= tmp16;
                                    EM_S[INDEX4(1,0,2,1,numEq,numComp,8)]+= tmp15;
                                    EM_S[INDEX4(1,0,2,2,numEq,numComp,8)]+=-tmp15;
                                    EM_S[INDEX4(1,0,2,3,numEq,numComp,8)]+=-tmp16;
                                    EM_S[INDEX4(1,0,2,4,numEq,numComp,8)]+= tmp1;
                                    EM_S[INDEX4(1,0,2,5,numEq,numComp,8)]+= tmp0;
                                    EM_S[INDEX4(1,0,2,6,numEq,numComp,8)]+=-tmp0;
                                    EM_S[INDEX4(1,0,2,7,numEq,numComp,8)]+=-tmp1;
                                    EM_S[INDEX4(1,0,3,0,numEq,numComp,8)]+=-tmp15;
                                    EM_S[INDEX4(1,0,3,1,numEq,numComp,8)]+=-tmp16;
                                    EM_S[INDEX4(1,0,3,2,numEq,numComp,8)]+= tmp16;
                                    EM_S[INDEX4(1,0,3,3,numEq,numComp,8)]+= tmp15;
                                    EM_S[INDEX4(1,0,3,4,numEq,numComp,8)]+=-tmp0;
                                    EM_S[INDEX4(1,0,3,5,numEq,numComp,8)]+=-tmp1;
                                    EM_S[INDEX4(1,0,3,6,numEq,numComp,8)]+= tmp1;
                                    EM_S[INDEX4(1,0,3,7,numEq,numComp,8)]+= tmp0;
                                    EM_S[INDEX4(1,0,4,0,numEq,numComp,8)]+= tmp0;
                                    EM_S[INDEX4(1,0,4,1,numEq,numComp,8)]+= tmp1;
                                    EM_S[INDEX4(1,0,4,2,numEq,numComp,8)]+=-tmp1;
                                    EM_S[INDEX4(1,0,4,3,numEq,numComp,8)]+=-tmp0;
                                    EM_S[INDEX4(1,0,4,4,numEq,numComp,8)]+= tmp15;
                                    EM_S[INDEX4(1,0,4,5,numEq,numComp,8)]+= tmp16;
                                    EM_S[INDEX4(1,0,4,6,numEq,numComp,8)]+=-tmp16;
                                    EM_S[INDEX4(1,0,4,7,numEq,numComp,8)]+=-tmp15;
                                    EM_S[INDEX4(1,0,5,0,numEq,numComp,8)]+=-tmp1;
                                    EM_S[INDEX4(1,0,5,1,numEq,numComp,8)]+=-tmp0;
                                    EM_S[INDEX4(1,0,5,2,numEq,numComp,8)]+= tmp0;
                                    EM_S[INDEX4(1,0,5,3,numEq,numComp,8)]+= tmp1;
                                    EM_S[INDEX4(1,0,5,4,numEq,numComp,8)]+=-tmp16;
                                    EM_S[INDEX4(1,0,5,5,numEq,numComp,8)]+=-tmp15;
                                    EM_S[INDEX4(1,0,5,6,numEq,numComp,8)]+= tmp15;
                                    EM_S[INDEX4(1,0,5,7,numEq,numComp,8)]+= tmp16;
                                    EM_S[INDEX4(1,0,6,0,numEq,numComp,8)]+= tmp1;
                                    EM_S[INDEX4(1,0,6,1,numEq,numComp,8)]+= tmp0;
                                    EM_S[INDEX4(1,0,6,2,numEq,numComp,8)]+=-tmp0;
                                    EM_S[INDEX4(1,0,6,3,numEq,numComp,8)]+=-tmp1;
                                    EM_S[INDEX4(1,0,6,4,numEq,numComp,8)]+= tmp16;
                                    EM_S[INDEX4(1,0,6,5,numEq,numComp,8)]+= tmp15;
                                    EM_S[INDEX4(1,0,6,6,numEq,numComp,8)]+=-tmp15;
                                    EM_S[INDEX4(1,0,6,7,numEq,numComp,8)]+=-tmp16;
                                    EM_S[INDEX4(1,0,7,0,numEq,numComp,8)]+=-tmp0;
                                    EM_S[INDEX4(1,0,7,1,numEq,numComp,8)]+=-tmp1;
                                    EM_S[INDEX4(1,0,7,2,numEq,numComp,8)]+= tmp1;
                                    EM_S[INDEX4(1,0,7,3,numEq,numComp,8)]+= tmp0;
                                    EM_S[INDEX4(1,0,7,4,numEq,numComp,8)]+=-tmp15;
                                    EM_S[INDEX4(1,0,7,5,numEq,numComp,8)]+=-tmp16;
                                    EM_S[INDEX4(1,0,7,6,numEq,numComp,8)]+= tmp16;
                                    EM_S[INDEX4(1,0,7,7,numEq,numComp,8)]+= tmp15;
                                }
                                {
                                    const double tmp6 = 4*Aw1010;
                                    const double tmp7 = 2*Aw1010;
                                    const double tmp8 = 4*Aw1111;
                                    const double tmp9 = 2*Aw1111;
                                    const double tmp10 = 4*Aw1212;
                                    const double tmp11 = 2*Aw1212;
                                    EM_S[INDEX4(1,1,0,0,numEq,numComp,8)]+= -tmp6 - tmp8 - tmp10;
                                    EM_S[INDEX4(1,1,0,1,numEq,numComp,8)]+= tmp6 - tmp9 - tmp11;
                                    EM_S[INDEX4(1,1,0,2,numEq,numComp,8)]+= -tmp7 + tmp8 - tmp11;
                                    EM_S[INDEX4(1,1,0,3,numEq,numComp,8)]+= tmp7 + tmp9 - Aw1212;
                                    EM_S[INDEX4(1,1,0,4,numEq,numComp,8)]+= -tmp7 - tmp9 + tmp10;
                                    EM_S[INDEX4(1,1,0,5,numEq,numComp,8)]+= tmp7 - Aw1111 + tmp11;
                                    EM_S[INDEX4(1,1,0,6,numEq,numComp,8)]+= -Aw1010 + tmp9 + tmp11;
                                    EM_S[INDEX4(1,1,0,7,numEq,numComp,8)]+= Aw1010 + Aw1111 + Aw1212;
                                    EM_S[INDEX4(1,1,1,0,numEq,numComp,8)]+= tmp6 - tmp9 - tmp11;
                                    EM_S[INDEX4(1,1,1,1,numEq,numComp,8)]+= -tmp6 - tmp8 - tmp10;
                                    EM_S[INDEX4(1,1,1,2,numEq,numComp,8)]+= tmp7 + tmp9 - Aw1212;
                                    EM_S[INDEX4(1,1,1,3,numEq,numComp,8)]+= -tmp7 + tmp8 - tmp11;
                                    EM_S[INDEX4(1,1,1,4,numEq,numComp,8)]+= tmp7 - Aw1111 + tmp11;
                                    EM_S[INDEX4(1,1,1,5,numEq,numComp,8)]+= -tmp7 - tmp9 + tmp10;
                                    EM_S[INDEX4(1,1,1,6,numEq,numComp,8)]+= Aw1010 + Aw1111 + Aw1212;
                                    EM_S[INDEX4(1,1,1,7,numEq,numComp,8)]+= -Aw1010 + tmp9 + tmp11;
                                    EM_S[INDEX4(1,1,2,0,numEq,numComp,8)]+= -tmp7 + tmp8 - tmp11;
                                    EM_S[INDEX4(1,1,2,1,numEq,numComp,8)]+= tmp7 + tmp9 - Aw1212;
                                    EM_S[INDEX4(1,1,2,2,numEq,numComp,8)]+= -tmp6 - tmp8 - tmp10;
                                    EM_S[INDEX4(1,1,2,3,numEq,numComp,8)]+= tmp6 - tmp9 - tmp11;
                                    EM_S[INDEX4(1,1,2,4,numEq,numComp,8)]+= -Aw1010 + tmp9 + tmp11;
                                    EM_S[INDEX4(1,1,2,5,numEq,numComp,8)]+= Aw1010 + Aw1111 + Aw1212;
                                    EM_S[INDEX4(1,1,2,6,numEq,numComp,8)]+= -tmp7 - tmp9 + tmp10;
                                    EM_S[INDEX4(1,1,2,7,numEq,numComp,8)]+= tmp7 - Aw1111 + tmp11;
                                    EM_S[INDEX4(1,1,3,0,numEq,numComp,8)]+= tmp7 + tmp9 - Aw1212;
                                    EM_S[INDEX4(1,1,3,1,numEq,numComp,8)]+= -tmp7 + tmp8 - tmp11;
                                    EM_S[INDEX4(1,1,3,2,numEq,numComp,8)]+= tmp6 - tmp9 - tmp11;
                                    EM_S[INDEX4(1,1,3,3,numEq,numComp,8)]+= -tmp6 - tmp8 - tmp10;
                                    EM_S[INDEX4(1,1,3,4,numEq,numComp,8)]+= Aw1010 + Aw1111 + Aw1212;
                                    EM_S[INDEX4(1,1,3,5,numEq,numComp,8)]+= -Aw1010 + tmp9 + tmp11;
                                    EM_S[INDEX4(1,1,3,6,numEq,numComp,8)]+= tmp7 - Aw1111 + tmp11;
                                    EM_S[INDEX4(1,1,3,7,numEq,numComp,8)]+= -tmp7 - tmp9 + tmp10;
                                    EM_S[INDEX4(1,1,4,0,numEq,numComp,8)]+= -tmp7 - tmp9 + tmp10;
                                    EM_S[INDEX4(1,1,4,1,numEq,numComp,8)]+= tmp7 - Aw1111 + tmp11;
                                    EM_S[INDEX4(1,1,4,2,numEq,numComp,8)]+= -Aw1010 + tmp9 + tmp11;
                                    EM_S[INDEX4(1,1,4,3,numEq,numComp,8)]+= Aw1010 + Aw1111 + Aw1212;
                                    EM_S[INDEX4(1,1,4,4,numEq,numComp,8)]+= -tmp6 - tmp8 - tmp10;
                                    EM_S[INDEX4(1,1,4,5,numEq,numComp,8)]+= tmp6 - tmp9 - tmp11;
                                    EM_S[INDEX4(1,1,4,6,numEq,numComp,8)]+= -tmp7 + tmp8 - tmp11;
                                    EM_S[INDEX4(1,1,4,7,numEq,numComp,8)]+= tmp7 + tmp9 - Aw1212;
                                    EM_S[INDEX4(1,1,5,0,numEq,numComp,8)]+= tmp7 - Aw1111 + tmp11;
                                    EM_S[INDEX4(1,1,5,1,numEq,numComp,8)]+= -tmp7 - tmp9 + tmp10;
                                    EM_S[INDEX4(1,1,5,2,numEq,numComp,8)]+= Aw1010 + Aw1111 + Aw1212;
                                    EM_S[INDEX4(1,1,5,3,numEq,numComp,8)]+= -Aw1010 + tmp9 + tmp11;
                                    EM_S[INDEX4(1,1,5,4,numEq,numComp,8)]+= tmp6 - tmp9 - tmp11;
                                    EM_S[INDEX4(1,1,5,5,numEq,numComp,8)]+= -tmp6 - tmp8 - tmp10;
                                    EM_S[INDEX4(1,1,5,6,numEq,numComp,8)]+= tmp7 + tmp9 - Aw1212;
                                    EM_S[INDEX4(1,1,5,7,numEq,numComp,8)]+= -tmp7 + tmp8 - tmp11;
                                    EM_S[INDEX4(1,1,6,0,numEq,numComp,8)]+= -Aw1010 + tmp9 + tmp11;
                                    EM_S[INDEX4(1,1,6,1,numEq,numComp,8)]+= Aw1010 + Aw1111 + Aw1212;
                                    EM_S[INDEX4(1,1,6,2,numEq,numComp,8)]+= -tmp7 - tmp9 + tmp10;
                                    EM_S[INDEX4(1,1,6,3,numEq,numComp,8)]+= tmp7 - Aw1111 + tmp11;
                                    EM_S[INDEX4(1,1,6,4,numEq,numComp,8)]+= -tmp7 + tmp8 - tmp11;
                                    EM_S[INDEX4(1,1,6,5,numEq,numComp,8)]+= tmp7 + tmp9 - Aw1212;
                                    EM_S[INDEX4(1,1,6,6,numEq,numComp,8)]+= -tmp6 - tmp8 - tmp10;
                                    EM_S[INDEX4(1,1,6,7,numEq,numComp,8)]+= tmp6 - tmp9 - tmp11;
                                    EM_S[INDEX4(1,1,7,0,numEq,numComp,8)]+= Aw1010 + Aw1111 + Aw1212;
                                    EM_S[INDEX4(1,1,7,1,numEq,numComp,8)]+= -Aw1010 + tmp9 + tmp11;
                                    EM_S[INDEX4(1,1,7,2,numEq,numComp,8)]+= tmp7 - Aw1111 + tmp11;
                                    EM_S[INDEX4(1,1,7,3,numEq,numComp,8)]+= -tmp7 - tmp9 + tmp10;
                                    EM_S[INDEX4(1,1,7,4,numEq,numComp,8)]+= tmp7 + tmp9 - Aw1212;
                                    EM_S[INDEX4(1,1,7,5,numEq,numComp,8)]+= -tmp7 + tmp8 - tmp11;
                                    EM_S[INDEX4(1,1,7,6,numEq,numComp,8)]+= tmp6 - tmp9 - tmp11;
                                    EM_S[INDEX4(1,1,7,7,numEq,numComp,8)]+= -tmp6 - tmp8 - tmp10;
                                }
                                {
                                    const double tmp4 = Aw1122 + Aw1221;
                                    const double tmp5 = Aw1122 - Aw1221;
                                    const double tmp19 = 2*tmp4;
                                    const double tmp20 = 2*tmp5;
                                    EM_S[INDEX4(1,2,0,0,numEq,numComp,8)]+=-tmp19;
                                    EM_S[INDEX4(1,2,0,1,numEq,numComp,8)]+=-tmp4;
                                    EM_S[INDEX4(1,2,0,2,numEq,numComp,8)]+=-tmp20;
                                    EM_S[INDEX4(1,2,0,3,numEq,numComp,8)]+=-tmp5;
                                    EM_S[INDEX4(1,2,0,4,numEq,numComp,8)]+= tmp20;
                                    EM_S[INDEX4(1,2,0,5,numEq,numComp,8)]+= tmp5;
                                    EM_S[INDEX4(1,2,0,6,numEq,numComp,8)]+= tmp19;
                                    EM_S[INDEX4(1,2,0,7,numEq,numComp,8)]+= tmp4;
                                    EM_S[INDEX4(1,2,1,0,numEq,numComp,8)]+=-tmp4;
                                    EM_S[INDEX4(1,2,1,1,numEq,numComp,8)]+=-tmp19;
                                    EM_S[INDEX4(1,2,1,2,numEq,numComp,8)]+=-tmp5;
                                    EM_S[INDEX4(1,2,1,3,numEq,numComp,8)]+=-tmp20;
                                    EM_S[INDEX4(1,2,1,4,numEq,numComp,8)]+= tmp5;
                                    EM_S[INDEX4(1,2,1,5,numEq,numComp,8)]+= tmp20;
                                    EM_S[INDEX4(1,2,1,6,numEq,numComp,8)]+= tmp4;
                                    EM_S[INDEX4(1,2,1,7,numEq,numComp,8)]+= tmp19;
                                    EM_S[INDEX4(1,2,2,0,numEq,numComp,8)]+= tmp20;
                                    EM_S[INDEX4(1,2,2,1,numEq,numComp,8)]+= tmp5;
                                    EM_S[INDEX4(1,2,2,2,numEq,numComp,8)]+= tmp19;
                                    EM_S[INDEX4(1,2,2,3,numEq,numComp,8)]+= tmp4;
                                    EM_S[INDEX4(1,2,2,4,numEq,numComp,8)]+=-tmp19;
                                    EM_S[INDEX4(1,2,2,5,numEq,numComp,8)]+=-tmp4;
                                    EM_S[INDEX4(1,2,2,6,numEq,numComp,8)]+=-tmp20;
                                    EM_S[INDEX4(1,2,2,7,numEq,numComp,8)]+=-tmp5;
                                    EM_S[INDEX4(1,2,3,0,numEq,numComp,8)]+= tmp5;
                                    EM_S[INDEX4(1,2,3,1,numEq,numComp,8)]+= tmp20;
                                    EM_S[INDEX4(1,2,3,2,numEq,numComp,8)]+= tmp4;
                                    EM_S[INDEX4(1,2,3,3,numEq,numComp,8)]+= tmp19;
                                    EM_S[INDEX4(1,2,3,4,numEq,numComp,8)]+=-tmp4;
                                    EM_S[INDEX4(1,2,3,5,numEq,numComp,8)]+=-tmp19;
                                    EM_S[INDEX4(1,2,3,6,numEq,numComp,8)]+=-tmp5;
                                    EM_S[INDEX4(1,2,3,7,numEq,numComp,8)]+=-tmp20;
                                    EM_S[INDEX4(1,2,4,0,numEq,numComp,8)]+=-tmp20;
                                    EM_S[INDEX4(1,2,4,1,numEq,numComp,8)]+=-tmp5;
                                    EM_S[INDEX4(1,2,4,2,numEq,numComp,8)]+=-tmp19;
                                    EM_S[INDEX4(1,2,4,3,numEq,numComp,8)]+=-tmp4;
                                    EM_S[INDEX4(1,2,4,4,numEq,numComp,8)]+= tmp19;
                                    EM_S[INDEX4(1,2,4,5,numEq,numComp,8)]+= tmp4;
                                    EM_S[INDEX4(1,2,4,6,numEq,numComp,8)]+= tmp20;
                                    EM_S[INDEX4(1,2,4,7,numEq,numComp,8)]+= tmp5;
                                    EM_S[INDEX4(1,2,5,0,numEq,numComp,8)]+=-tmp5;
                                    EM_S[INDEX4(1,2,5,1,numEq,numComp,8)]+=-tmp20;
                                    EM_S[INDEX4(1,2,5,2,numEq,numComp,8)]+=-tmp4;
                                    EM_S[INDEX4(1,2,5,3,numEq,numComp,8)]+=-tmp19;
                                    EM_S[INDEX4(1,2,5,4,numEq,numComp,8)]+= tmp4;
                                    EM_S[INDEX4(1,2,5,5,numEq,numComp,8)]+= tmp19;
                                    EM_S[INDEX4(1,2,5,6,numEq,numComp,8)]+= tmp5;
                                    EM_S[INDEX4(1,2,5,7,numEq,numComp,8)]+= tmp20;
                                    EM_S[INDEX4(1,2,6,0,numEq,numComp,8)]+= tmp19;
                                    EM_S[INDEX4(1,2,6,1,numEq,numComp,8)]+= tmp4;
                                    EM_S[INDEX4(1,2,6,2,numEq,numComp,8)]+= tmp20;
                                    EM_S[INDEX4(1,2,6,3,numEq,numComp,8)]+= tmp5;
                                    EM_S[INDEX4(1,2,6,4,numEq,numComp,8)]+=-tmp20;
                                    EM_S[INDEX4(1,2,6,5,numEq,numComp,8)]+=-tmp5;
                                    EM_S[INDEX4(1,2,6,6,numEq,numComp,8)]+=-tmp19;
                                    EM_S[INDEX4(1,2,6,7,numEq,numComp,8)]+=-tmp4;
                                    EM_S[INDEX4(1,2,7,0,numEq,numComp,8)]+= tmp4;
                                    EM_S[INDEX4(1,2,7,1,numEq,numComp,8)]+= tmp19;
                                    EM_S[INDEX4(1,2,7,2,numEq,numComp,8)]+= tmp5;
                                    EM_S[INDEX4(1,2,7,3,numEq,numComp,8)]+= tmp20;
                                    EM_S[INDEX4(1,2,7,4,numEq,numComp,8)]+=-tmp5;
                                    EM_S[INDEX4(1,2,7,5,numEq,numComp,8)]+=-tmp20;
                                    EM_S[INDEX4(1,2,7,6,numEq,numComp,8)]+=-tmp4;
                                    EM_S[INDEX4(1,2,7,7,numEq,numComp,8)]+=-tmp19;
                                }
                                {
                                    const double tmp2 = Aw2002 + Aw2200;
                                    const double tmp3 = Aw2002 - Aw2200;
                                    const double tmp17 = 2*tmp2;
                                    const double tmp18 = 2*tmp3;
                                    EM_S[INDEX4(2,0,0,0,numEq,numComp,8)]+= tmp17;
                                    EM_S[INDEX4(2,0,0,1,numEq,numComp,8)]+= tmp18;
                                    EM_S[INDEX4(2,0,0,2,numEq,numComp,8)]+= tmp2;
                                    EM_S[INDEX4(2,0,0,3,numEq,numComp,8)]+= tmp3;
                                    EM_S[INDEX4(2,0,0,4,numEq,numComp,8)]+=-tmp18;
                                    EM_S[INDEX4(2,0,0,5,numEq,numComp,8)]+=-tmp17;
                                    EM_S[INDEX4(2,0,0,6,numEq,numComp,8)]+=-tmp3;
                                    EM_S[INDEX4(2,0,0,7,numEq,numComp,8)]+=-tmp2;
                                    EM_S[INDEX4(2,0,1,0,numEq,numComp,8)]+=-tmp18;
                                    EM_S[INDEX4(2,0,1,1,numEq,numComp,8)]+=-tmp17;
                                    EM_S[INDEX4(2,0,1,2,numEq,numComp,8)]+=-tmp3;
                                    EM_S[INDEX4(2,0,1,3,numEq,numComp,8)]+=-tmp2;
                                    EM_S[INDEX4(2,0,1,4,numEq,numComp,8)]+= tmp17;
                                    EM_S[INDEX4(2,0,1,5,numEq,numComp,8)]+= tmp18;
                                    EM_S[INDEX4(2,0,1,6,numEq,numComp,8)]+= tmp2;
                                    EM_S[INDEX4(2,0,1,7,numEq,numComp,8)]+= tmp3;
                                    EM_S[INDEX4(2,0,2,0,numEq,numComp,8)]+= tmp2;
                                    EM_S[INDEX4(2,0,2,1,numEq,numComp,8)]+= tmp3;
                                    EM_S[INDEX4(2,0,2,2,numEq,numComp,8)]+= tmp17;
                                    EM_S[INDEX4(2,0,2,3,numEq,numComp,8)]+= tmp18;
                                    EM_S[INDEX4(2,0,2,4,numEq,numComp,8)]+=-tmp3;
                                    EM_S[INDEX4(2,0,2,5,numEq,numComp,8)]+=-tmp2;
                                    EM_S[INDEX4(2,0,2,6,numEq,numComp,8)]+=-tmp18;
                                    EM_S[INDEX4(2,0,2,7,numEq,numComp,8)]+=-tmp17;
                                    EM_S[INDEX4(2,0,3,0,numEq,numComp,8)]+=-tmp3;
                                    EM_S[INDEX4(2,0,3,1,numEq,numComp,8)]+=-tmp2;
                                    EM_S[INDEX4(2,0,3,2,numEq,numComp,8)]+=-tmp18;
                                    EM_S[INDEX4(2,0,3,3,numEq,numComp,8)]+=-tmp17;
                                    EM_S[INDEX4(2,0,3,4,numEq,numComp,8)]+= tmp2;
                                    EM_S[INDEX4(2,0,3,5,numEq,numComp,8)]+= tmp3;
                                    EM_S[INDEX4(2,0,3,6,numEq,numComp,8)]+= tmp17;
                                    EM_S[INDEX4(2,0,3,7,numEq,numComp,8)]+= tmp18;
                                    EM_S[INDEX4(2,0,4,0,numEq,numComp,8)]+= tmp18;
                                    EM_S[INDEX4(2,0,4,1,numEq,numComp,8)]+= tmp17;
                                    EM_S[INDEX4(2,0,4,2,numEq,numComp,8)]+= tmp3;
                                    EM_S[INDEX4(2,0,4,3,numEq,numComp,8)]+= tmp2;
                                    EM_S[INDEX4(2,0,4,4,numEq,numComp,8)]+=-tmp17;
                                    EM_S[INDEX4(2,0,4,5,numEq,numComp,8)]+=-tmp18;
                                    EM_S[INDEX4(2,0,4,6,numEq,numComp,8)]+=-tmp2;
                                    EM_S[INDEX4(2,0,4,7,numEq,numComp,8)]+=-tmp3;
                                    EM_S[INDEX4(2,0,5,0,numEq,numComp,8)]+=-tmp17;
                                    EM_S[INDEX4(2,0,5,1,numEq,numComp,8)]+=-tmp18;
                                    EM_S[INDEX4(2,0,5,2,numEq,numComp,8)]+=-tmp2;
                                    EM_S[INDEX4(2,0,5,3,numEq,numComp,8)]+=-tmp3;
                                    EM_S[INDEX4(2,0,5,4,numEq,numComp,8)]+= tmp18;
                                    EM_S[INDEX4(2,0,5,5,numEq,numComp,8)]+= tmp17;
                                    EM_S[INDEX4(2,0,5,6,numEq,numComp,8)]+= tmp3;
                                    EM_S[INDEX4(2,0,5,7,numEq,numComp,8)]+= tmp2;
                                    EM_S[INDEX4(2,0,6,0,numEq,numComp,8)]+= tmp3;
                                    EM_S[INDEX4(2,0,6,1,numEq,numComp,8)]+= tmp2;
                                    EM_S[INDEX4(2,0,6,2,numEq,numComp,8)]+= tmp18;
                                    EM_S[INDEX4(2,0,6,3,numEq,numComp,8)]+= tmp17;
                                    EM_S[INDEX4(2,0,6,4,numEq,numComp,8)]+=-tmp2;
                                    EM_S[INDEX4(2,0,6,5,numEq,numComp,8)]+=-tmp3;
                                    EM_S[INDEX4(2,0,6,6,numEq,numComp,8)]+=-tmp17;
                                    EM_S[INDEX4(2,0,6,7,numEq,numComp,8)]+=-tmp18;
                                    EM_S[INDEX4(2,0,7,0,numEq,numComp,8)]+=-tmp2;
                                    EM_S[INDEX4(2,0,7,1,numEq,numComp,8)]+=-tmp3;
                                    EM_S[INDEX4(2,0,7,2,numEq,numComp,8)]+=-tmp17;
                                    EM_S[INDEX4(2,0,7,3,numEq,numComp,8)]+=-tmp18;
                                    EM_S[INDEX4(2,0,7,4,numEq,numComp,8)]+= tmp3;
                                    EM_S[INDEX4(2,0,7,5,numEq,numComp,8)]+= tmp2;
                                    EM_S[INDEX4(2,0,7,6,numEq,numComp,8)]+= tmp18;
                                    EM_S[INDEX4(2,0,7,7,numEq,numComp,8)]+= tmp17;
                                }
                                {
                                    const double tmp4 = Aw2112 + Aw2211;
                                    const double tmp5 = Aw2112 - Aw2211;
                                    const double tmp19 = 2*tmp4;
                                    const double tmp20 = 2*tmp5;
                                    EM_S[INDEX4(2,1,0,0,numEq,numComp,8)]+=-tmp19;
                                    EM_S[INDEX4(2,1,0,1,numEq,numComp,8)]+=-tmp4;
                                    EM_S[INDEX4(2,1,0,2,numEq,numComp,8)]+=-tmp20;
                                    EM_S[INDEX4(2,1,0,3,numEq,numComp,8)]+=-tmp5;
                                    EM_S[INDEX4(2,1,0,4,numEq,numComp,8)]+= tmp20;
                                    EM_S[INDEX4(2,1,0,5,numEq,numComp,8)]+= tmp5;
                                    EM_S[INDEX4(2,1,0,6,numEq,numComp,8)]+= tmp19;
                                    EM_S[INDEX4(2,1,0,7,numEq,numComp,8)]+= tmp4;
                                    EM_S[INDEX4(2,1,1,0,numEq,numComp,8)]+=-tmp4;
                                    EM_S[INDEX4(2,1,1,1,numEq,numComp,8)]+=-tmp19;
                                    EM_S[INDEX4(2,1,1,2,numEq,numComp,8)]+=-tmp5;
                                    EM_S[INDEX4(2,1,1,3,numEq,numComp,8)]+=-tmp20;
                                    EM_S[INDEX4(2,1,1,4,numEq,numComp,8)]+= tmp5;
                                    EM_S[INDEX4(2,1,1,5,numEq,numComp,8)]+= tmp20;
                                    EM_S[INDEX4(2,1,1,6,numEq,numComp,8)]+= tmp4;
                                    EM_S[INDEX4(2,1,1,7,numEq,numComp,8)]+= tmp19;
                                    EM_S[INDEX4(2,1,2,0,numEq,numComp,8)]+= tmp20;
                                    EM_S[INDEX4(2,1,2,1,numEq,numComp,8)]+= tmp5;
                                    EM_S[INDEX4(2,1,2,2,numEq,numComp,8)]+= tmp19;
                                    EM_S[INDEX4(2,1,2,3,numEq,numComp,8)]+= tmp4;
                                    EM_S[INDEX4(2,1,2,4,numEq,numComp,8)]+=-tmp19;
                                    EM_S[INDEX4(2,1,2,5,numEq,numComp,8)]+=-tmp4;
                                    EM_S[INDEX4(2,1,2,6,numEq,numComp,8)]+=-tmp20;
                                    EM_S[INDEX4(2,1,2,7,numEq,numComp,8)]+=-tmp5;
                                    EM_S[INDEX4(2,1,3,0,numEq,numComp,8)]+= tmp5;
                                    EM_S[INDEX4(2,1,3,1,numEq,numComp,8)]+= tmp20;
                                    EM_S[INDEX4(2,1,3,2,numEq,numComp,8)]+= tmp4;
                                    EM_S[INDEX4(2,1,3,3,numEq,numComp,8)]+= tmp19;
                                    EM_S[INDEX4(2,1,3,4,numEq,numComp,8)]+=-tmp4;
                                    EM_S[INDEX4(2,1,3,5,numEq,numComp,8)]+=-tmp19;
                                    EM_S[INDEX4(2,1,3,6,numEq,numComp,8)]+=-tmp5;
                                    EM_S[INDEX4(2,1,3,7,numEq,numComp,8)]+=-tmp20;
                                    EM_S[INDEX4(2,1,4,0,numEq,numComp,8)]+=-tmp20;
                                    EM_S[INDEX4(2,1,4,1,numEq,numComp,8)]+=-tmp5;
                                    EM_S[INDEX4(2,1,4,2,numEq,numComp,8)]+=-tmp19;
                                    EM_S[INDEX4(2,1,4,3,numEq,numComp,8)]+=-tmp4;
                                    EM_S[INDEX4(2,1,4,4,numEq,numComp,8)]+= tmp19;
                                    EM_S[INDEX4(2,1,4,5,numEq,numComp,8)]+= tmp4;
                                    EM_S[INDEX4(2,1,4,6,numEq,numComp,8)]+= tmp20;
                                    EM_S[INDEX4(2,1,4,7,numEq,numComp,8)]+= tmp5;
                                    EM_S[INDEX4(2,1,5,0,numEq,numComp,8)]+=-tmp5;
                                    EM_S[INDEX4(2,1,5,1,numEq,numComp,8)]+=-tmp20;
                                    EM_S[INDEX4(2,1,5,2,numEq,numComp,8)]+=-tmp4;
                                    EM_S[INDEX4(2,1,5,3,numEq,numComp,8)]+=-tmp19;
                                    EM_S[INDEX4(2,1,5,4,numEq,numComp,8)]+= tmp4;
                                    EM_S[INDEX4(2,1,5,5,numEq,numComp,8)]+= tmp19;
                                    EM_S[INDEX4(2,1,5,6,numEq,numComp,8)]+= tmp5;
                                    EM_S[INDEX4(2,1,5,7,numEq,numComp,8)]+= tmp20;
                                    EM_S[INDEX4(2,1,6,0,numEq,numComp,8)]+= tmp19;
                                    EM_S[INDEX4(2,1,6,1,numEq,numComp,8)]+= tmp4;
                                    EM_S[INDEX4(2,1,6,2,numEq,numComp,8)]+= tmp20;
                                    EM_S[INDEX4(2,1,6,3,numEq,numComp,8)]+= tmp5;
                                    EM_S[INDEX4(2,1,6,4,numEq,numComp,8)]+=-tmp20;
                                    EM_S[INDEX4(2,1,6,5,numEq,numComp,8)]+=-tmp5;
                                    EM_S[INDEX4(2,1,6,6,numEq,numComp,8)]+=-tmp19;
                                    EM_S[INDEX4(2,1,6,7,numEq,numComp,8)]+=-tmp4;
                                    EM_S[INDEX4(2,1,7,0,numEq,numComp,8)]+= tmp4;
                                    EM_S[INDEX4(2,1,7,1,numEq,numComp,8)]+= tmp19;
                                    EM_S[INDEX4(2,1,7,2,numEq,numComp,8)]+= tmp5;
                                    EM_S[INDEX4(2,1,7,3,numEq,numComp,8)]+= tmp20;
                                    EM_S[INDEX4(2,1,7,4,numEq,numComp,8)]+=-tmp5;
                                    EM_S[INDEX4(2,1,7,5,numEq,numComp,8)]+=-tmp20;
                                    EM_S[INDEX4(2,1,7,6,numEq,numComp,8)]+=-tmp4;
                                    EM_S[INDEX4(2,1,7,7,numEq,numComp,8)]+=-tmp19;
                                }
                                {
                                    const double tmp6 = 4*Aw2020;
                                    const double tmp7 = 2*Aw2020;
                                    const double tmp8 = 4*Aw2121;
                                    const double tmp9 = 2*Aw2121;
                                    const double tmp10 = 4*Aw2222;
                                    const double tmp11 = 2*Aw2222;
                                    EM_S[INDEX4(2,2,0,0,numEq,numComp,8)]+= -tmp6 - tmp8 - tmp10;
                                    EM_S[INDEX4(2,2,0,1,numEq,numComp,8)]+= tmp6 - tmp9 - tmp11;
                                    EM_S[INDEX4(2,2,0,2,numEq,numComp,8)]+= -tmp7 + tmp8 - tmp11;
                                    EM_S[INDEX4(2,2,0,3,numEq,numComp,8)]+= tmp7 + tmp9 - Aw2222;
                                    EM_S[INDEX4(2,2,0,4,numEq,numComp,8)]+= -tmp7 - tmp9 + tmp10;
                                    EM_S[INDEX4(2,2,0,5,numEq,numComp,8)]+= tmp7 - Aw2121 + tmp11;
                                    EM_S[INDEX4(2,2,0,6,numEq,numComp,8)]+= -Aw2020 + tmp9 + tmp11;
                                    EM_S[INDEX4(2,2,0,7,numEq,numComp,8)]+= Aw2020 + Aw2121 + Aw2222;
                                    EM_S[INDEX4(2,2,1,0,numEq,numComp,8)]+= tmp6 - tmp9 - tmp11;
                                    EM_S[INDEX4(2,2,1,1,numEq,numComp,8)]+= -tmp6 - tmp8 - tmp10;
                                    EM_S[INDEX4(2,2,1,2,numEq,numComp,8)]+= tmp7 + tmp9 - Aw2222;
                                    EM_S[INDEX4(2,2,1,3,numEq,numComp,8)]+= -tmp7 + tmp8 - tmp11;
                                    EM_S[INDEX4(2,2,1,4,numEq,numComp,8)]+= tmp7 - Aw2121 + tmp11;
                                    EM_S[INDEX4(2,2,1,5,numEq,numComp,8)]+= -tmp7 - tmp9 + tmp10;
                                    EM_S[INDEX4(2,2,1,6,numEq,numComp,8)]+= Aw2020 + Aw2121 + Aw2222;
                                    EM_S[INDEX4(2,2,1,7,numEq,numComp,8)]+= -Aw2020 + tmp9 + tmp11;
                                    EM_S[INDEX4(2,2,2,0,numEq,numComp,8)]+= -tmp7 + tmp8 - tmp11;
                                    EM_S[INDEX4(2,2,2,1,numEq,numComp,8)]+= tmp7 + tmp9 - Aw2222;
                                    EM_S[INDEX4(2,2,2,2,numEq,numComp,8)]+= -tmp6 - tmp8 - tmp10;
                                    EM_S[INDEX4(2,2,2,3,numEq,numComp,8)]+= tmp6 - tmp9 - tmp11;
                                    EM_S[INDEX4(2,2,2,4,numEq,numComp,8)]+= -Aw2020 + tmp9 + tmp11;
                                    EM_S[INDEX4(2,2,2,5,numEq,numComp,8)]+= Aw2020 + Aw2121 + Aw2222;
                                    EM_S[INDEX4(2,2,2,6,numEq,numComp,8)]+= -tmp7 - tmp9 + tmp10;
                                    EM_S[INDEX4(2,2,2,7,numEq,numComp,8)]+= tmp7 - Aw2121 + tmp11;
                                    EM_S[INDEX4(2,2,3,0,numEq,numComp,8)]+= tmp7 + tmp9 - Aw2222;
                                    EM_S[INDEX4(2,2,3,1,numEq,numComp,8)]+= -tmp7 + tmp8 - tmp11;
                                    EM_S[INDEX4(2,2,3,2,numEq,numComp,8)]+= tmp6 - tmp9 - tmp11;
                                    EM_S[INDEX4(2,2,3,3,numEq,numComp,8)]+= -tmp6 - tmp8 - tmp10;
                                    EM_S[INDEX4(2,2,3,4,numEq,numComp,8)]+= Aw2020 + Aw2121 + Aw2222;
                                    EM_S[INDEX4(2,2,3,5,numEq,numComp,8)]+= -Aw2020 + tmp9 + tmp11;
                                    EM_S[INDEX4(2,2,3,6,numEq,numComp,8)]+= tmp7 - Aw2121 + tmp11;
                                    EM_S[INDEX4(2,2,3,7,numEq,numComp,8)]+= -tmp7 - tmp9 + tmp10;
                                    EM_S[INDEX4(2,2,4,0,numEq,numComp,8)]+= -tmp7 - tmp9 + tmp10;
                                    EM_S[INDEX4(2,2,4,1,numEq,numComp,8)]+= tmp7 - Aw2121 + tmp11;
                                    EM_S[INDEX4(2,2,4,2,numEq,numComp,8)]+= -Aw2020 + tmp9 + tmp11;
                                    EM_S[INDEX4(2,2,4,3,numEq,numComp,8)]+= Aw2020 + Aw2121 + Aw2222;
                                    EM_S[INDEX4(2,2,4,4,numEq,numComp,8)]+= -tmp6 - tmp8 - tmp10;
                                    EM_S[INDEX4(2,2,4,5,numEq,numComp,8)]+= tmp6 - tmp9 - tmp11;
                                    EM_S[INDEX4(2,2,4,6,numEq,numComp,8)]+= -tmp7 + tmp8 - tmp11;
                                    EM_S[INDEX4(2,2,4,7,numEq,numComp,8)]+= tmp7 + tmp9 - Aw2222;
                                    EM_S[INDEX4(2,2,5,0,numEq,numComp,8)]+= tmp7 - Aw2121 + tmp11;
                                    EM_S[INDEX4(2,2,5,1,numEq,numComp,8)]+= -tmp7 - tmp9 + tmp10;
                                    EM_S[INDEX4(2,2,5,2,numEq,numComp,8)]+= Aw2020 + Aw2121 + Aw2222;
                                    EM_S[INDEX4(2,2,5,3,numEq,numComp,8)]+= -Aw2020 + tmp9 + tmp11;
                                    EM_S[INDEX4(2,2,5,4,numEq,numComp,8)]+= tmp6 - tmp9 - tmp11;
                                    EM_S[INDEX4(2,2,5,5,numEq,numComp,8)]+= -tmp6 - tmp8 - tmp10;
                                    EM_S[INDEX4(2,2,5,6,numEq,numComp,8)]+= tmp7 + tmp9 - Aw2222;
                                    EM_S[INDEX4(2,2,5,7,numEq,numComp,8)]+= -tmp7 + tmp8 - tmp11;
                                    EM_S[INDEX4(2,2,6,0,numEq,numComp,8)]+= -Aw2020 + tmp9 + tmp11;
                                    EM_S[INDEX4(2,2,6,1,numEq,numComp,8)]+= Aw2020 + Aw2121 + Aw2222;
                                    EM_S[INDEX4(2,2,6,2,numEq,numComp,8)]+= -tmp7 - tmp9 + tmp10;
                                    EM_S[INDEX4(2,2,6,3,numEq,numComp,8)]+= tmp7 - Aw2121 + tmp11;
                                    EM_S[INDEX4(2,2,6,4,numEq,numComp,8)]+= -tmp7 + tmp8 - tmp11;
                                    EM_S[INDEX4(2,2,6,5,numEq,numComp,8)]+= tmp7 + tmp9 - Aw2222;
                                    EM_S[INDEX4(2,2,6,6,numEq,numComp,8)]+= -tmp6 - tmp8 - tmp10;
                                    EM_S[INDEX4(2,2,6,7,numEq,numComp,8)]+= tmp6 - tmp9 - tmp11;
                                    EM_S[INDEX4(2,2,7,0,numEq,numComp,8)]+= Aw2020 + Aw2121 + Aw2222;
                                    EM_S[INDEX4(2,2,7,1,numEq,numComp,8)]+= -Aw2020 + tmp9 + tmp11;
                                    EM_S[INDEX4(2,2,7,2,numEq,numComp,8)]+= tmp7 - Aw2121 + tmp11;
                                    EM_S[INDEX4(2,2,7,3,numEq,numComp,8)]+= -tmp7 - tmp9 + tmp10;
                                    EM_S[INDEX4(2,2,7,4,numEq,numComp,8)]+= tmp7 + tmp9 - Aw2222;
                                    EM_S[INDEX4(2,2,7,5,numEq,numComp,8)]+= -tmp7 + tmp8 - tmp11;
                                    EM_S[INDEX4(2,2,7,6,numEq,numComp,8)]+= tmp6 - tmp9 - tmp11;
                                    EM_S[INDEX4(2,2,7,7,numEq,numComp,8)]+= -tmp6 - tmp8 - tmp10;
                                }
                            }
                        }
                        ///////////////
                        // process B //
                        ///////////////
                        if (!B.isEmpty()) {
                            const double* B_p=B.getSampleDataRO(e);
                            if (B.actsExpanded()) {
                                for (index_t k=0; k<numEq; k++) {
                                    for (index_t m=0; m<numComp; m++) {
                                        const double B_0_0 = B_p[INDEX4(k,0,m,0, numEq,3,numComp)];
                                        const double B_1_0 = B_p[INDEX4(k,1,m,0, numEq,3,numComp)];
                                        const double B_2_0 = B_p[INDEX4(k,2,m,0, numEq,3,numComp)];
                                        const double B_0_1 = B_p[INDEX4(k,0,m,1, numEq,3,numComp)];
                                        const double B_1_1 = B_p[INDEX4(k,1,m,1, numEq,3,numComp)];
                                        const double B_2_1 = B_p[INDEX4(k,2,m,1, numEq,3,numComp)];
                                        const double B_0_2 = B_p[INDEX4(k,0,m,2, numEq,3,numComp)];
                                        const double B_1_2 = B_p[INDEX4(k,1,m,2, numEq,3,numComp)];
                                        const double B_2_2 = B_p[INDEX4(k,2,m,2, numEq,3,numComp)];
                                        const double B_0_3 = B_p[INDEX4(k,0,m,3, numEq,3,numComp)];
                                        const double B_1_3 = B_p[INDEX4(k,1,m,3, numEq,3,numComp)];
                                        const double B_2_3 = B_p[INDEX4(k,2,m,3, numEq,3,numComp)];
                                        const double B_0_4 = B_p[INDEX4(k,0,m,4, numEq,3,numComp)];
                                        const double B_1_4 = B_p[INDEX4(k,1,m,4, numEq,3,numComp)];
                                        const double B_2_4 = B_p[INDEX4(k,2,m,4, numEq,3,numComp)];
                                        const double B_0_5 = B_p[INDEX4(k,0,m,5, numEq,3,numComp)];
                                        const double B_1_5 = B_p[INDEX4(k,1,m,5, numEq,3,numComp)];
                                        const double B_2_5 = B_p[INDEX4(k,2,m,5, numEq,3,numComp)];
                                        const double B_0_6 = B_p[INDEX4(k,0,m,6, numEq,3,numComp)];
                                        const double B_1_6 = B_p[INDEX4(k,1,m,6, numEq,3,numComp)];
                                        const double B_2_6 = B_p[INDEX4(k,2,m,6, numEq,3,numComp)];
                                        const double B_0_7 = B_p[INDEX4(k,0,m,7, numEq,3,numComp)];
                                        const double B_1_7 = B_p[INDEX4(k,1,m,7, numEq,3,numComp)];
                                        const double B_2_7 = B_p[INDEX4(k,2,m,7, numEq,3,numComp)];
                                        const double tmp0 = w38*(B_2_1 + B_2_2);
                                        const double tmp1 = w42*(B_1_3 + B_1_7);
                                        const double tmp2 = w41*(B_0_3 + B_0_7);
                                        const double tmp3 = w37*(B_1_1 + B_1_5);
                                        const double tmp4 = w39*(B_0_2 + B_0_6);
                                        const double tmp5 = w45*(B_2_5 + B_2_6);
                                        const double tmp6 = w36*(B_0_1 + B_0_5);
                                        const double tmp7 = w40*(B_1_2 + B_1_6);
                                        const double tmp8 = w33*(B_0_0 + B_0_4);
                                        const double tmp9 = w34*(B_1_0 + B_1_4);
                                        const double tmp10 = w38*(B_2_4 + B_2_5 + B_2_6 + B_2_7);
                                        const double tmp11 = w42*(-B_1_6 - B_1_7);
                                        const double tmp12 = w41*(-B_0_5 - B_0_7);
                                        const double tmp13 = w37*(-B_1_4 - B_1_5);
                                        const double tmp14 = w39*(-B_0_4 - B_0_6);
                                        const double tmp15 = w45*(B_2_0 + B_2_1 + B_2_2 + B_2_3);
                                        const double tmp16 = w36*(-B_0_1 - B_0_3);
                                        const double tmp17 = w40*(-B_1_2 - B_1_3);
                                        const double tmp18 = w33*(-B_0_0 - B_0_2);
                                        const double tmp19 = w34*(-B_1_0 - B_1_1);
                                        const double tmp20 = w38*(-B_2_5 - B_2_7);
                                        const double tmp21 = w35*(-B_2_4 - B_2_6);
                                        const double tmp22 = w41*(B_0_1 + B_0_3);
                                        const double tmp23 = w37*(-B_1_2 - B_1_7);
                                        const double tmp24 = w39*(B_0_0 + B_0_2);
                                        const double tmp25 = w45*(-B_2_0 - B_2_2);
                                        const double tmp26 = w36*(B_0_5 + B_0_7);
                                        const double tmp27 = w40*(-B_1_0 - B_1_5);
                                        const double tmp28 = w33*(B_0_4 + B_0_6);
                                        const double tmp29 = w46*(-B_2_1 - B_2_3);
                                        const double tmp30 = w38*(B_2_0 + B_2_2);
                                        const double tmp31 = w35*(B_2_1 + B_2_3);
                                        const double tmp32 = w41*(-B_0_4 - B_0_6);
                                        const double tmp33 = w37*(B_1_0 + B_1_5);
                                        const double tmp34 = w39*(-B_0_5 - B_0_7);
                                        const double tmp35 = w45*(B_2_5 + B_2_7);
                                        const double tmp36 = w36*(-B_0_0 - B_0_2);
                                        const double tmp37 = w40*(B_1_2 + B_1_7);
                                        const double tmp38 = w33*(-B_0_1 - B_0_3);
                                        const double tmp39 = w46*(B_2_4 + B_2_6);
                                        const double tmp40 = w38*(-B_2_0 - B_2_1 - B_2_2 - B_2_3);
                                        const double tmp41 = w42*(B_1_0 + B_1_1);
                                        const double tmp42 = w41*(B_0_0 + B_0_2);
                                        const double tmp43 = w37*(B_1_2 + B_1_3);
                                        const double tmp44 = w39*(B_0_1 + B_0_3);
                                        const double tmp45 = w45*(-B_2_4 - B_2_5 - B_2_6 - B_2_7);
                                        const double tmp46 = w36*(B_0_4 + B_0_6);
                                        const double tmp47 = w40*(B_1_4 + B_1_5);
                                        const double tmp48 = w33*(B_0_5 + B_0_7);
                                        const double tmp49 = w34*(B_1_6 + B_1_7);
                                        const double tmp50 = w38*(B_2_0 + B_2_1);
                                        const double tmp51 = w42*(-B_1_4 - B_1_5);
                                        const double tmp52 = w35*(B_2_2 + B_2_3);
                                        const double tmp53 = w37*(-B_1_6 - B_1_7);
                                        const double tmp54 = w39*(B_0_0 + B_0_6);
                                        const double tmp55 = w45*(B_2_6 + B_2_7);
                                        const double tmp56 = w36*(B_0_1 + B_0_7);
                                        const double tmp57 = w40*(-B_1_0 - B_1_1);
                                        const double tmp58 = w46*(B_2_4 + B_2_5);
                                        const double tmp59 = w34*(-B_1_2 - B_1_3);
                                        const double tmp60 = w38*(-B_2_4 - B_2_5 - B_2_6 - B_2_7);
                                        const double tmp61 = w37*(-B_1_2 - B_1_3 - B_1_6 - B_1_7);
                                        const double tmp62 = w39*(-B_0_1 - B_0_3 - B_0_5 - B_0_7);
                                        const double tmp63 = w45*(-B_2_0 - B_2_1 - B_2_2 - B_2_3);
                                        const double tmp64 = w36*(-B_0_0 - B_0_2 - B_0_4 - B_0_6);
                                        const double tmp65 = w40*(-B_1_0 - B_1_1 - B_1_4 - B_1_5);
                                        const double tmp66 = w41*(B_0_4 + B_0_6);
                                        const double tmp67 = w39*(B_0_5 + B_0_7);
                                        const double tmp68 = w36*(B_0_0 + B_0_2);
                                        const double tmp69 = w33*(B_0_1 + B_0_3);
                                        const double tmp70 = w38*(-B_2_4 - B_2_7);
                                        const double tmp71 = w42*(B_1_2 + B_1_6);
                                        const double tmp72 = w41*(-B_0_2 - B_0_6);
                                        const double tmp73 = w37*(B_1_0 + B_1_4);
                                        const double tmp74 = w39*(-B_0_3 - B_0_7);
                                        const double tmp75 = w45*(-B_2_0 - B_2_3);
                                        const double tmp76 = w36*(-B_0_0 - B_0_4);
                                        const double tmp77 = w40*(B_1_3 + B_1_7);
                                        const double tmp78 = w33*(-B_0_1 - B_0_5);
                                        const double tmp79 = w34*(B_1_1 + B_1_5);
                                        const double tmp80 = w39*(B_0_0 + B_0_2 + B_0_4 + B_0_6);
                                        const double tmp81 = w36*(B_0_1 + B_0_3 + B_0_5 + B_0_7);
                                        const double tmp82 = w38*(B_2_0 + B_2_3);
                                        const double tmp83 = w42*(-B_1_1 - B_1_5);
                                        const double tmp84 = w41*(B_0_1 + B_0_5);
                                        const double tmp85 = w37*(-B_1_3 - B_1_7);
                                        const double tmp86 = w39*(B_0_0 + B_0_4);
                                        const double tmp87 = w45*(B_2_4 + B_2_7);
                                        const double tmp88 = w36*(B_0_3 + B_0_7);
                                        const double tmp89 = w40*(-B_1_0 - B_1_4);
                                        const double tmp90 = w33*(B_0_2 + B_0_6);
                                        const double tmp91 = w34*(-B_1_2 - B_1_6);
                                        const double tmp92 = w38*(-B_2_5 - B_2_6);
                                        const double tmp93 = w45*(-B_2_1 - B_2_2);
                                        const double tmp94 = w37*(B_1_0 + B_1_1 + B_1_4 + B_1_5);
                                        const double tmp95 = w40*(B_1_2 + B_1_3 + B_1_6 + B_1_7);
                                        const double tmp96 = w42*(-B_1_2 - B_1_3);
                                        const double tmp97 = w41*(-B_0_1 - B_0_3);
                                        const double tmp98 = w37*(-B_1_0 - B_1_1);
                                        const double tmp99 = w39*(-B_0_0 - B_0_2);
                                        const double tmp100 = w36*(-B_0_5 - B_0_7);
                                        const double tmp101 = w40*(-B_1_6 - B_1_7);
                                        const double tmp102 = w33*(-B_0_4 - B_0_6);
                                        const double tmp103 = w34*(-B_1_4 - B_1_5);
                                        const double tmp104 = w38*(B_2_6 + B_2_7);
                                        const double tmp105 = w35*(B_2_4 + B_2_5);
                                        const double tmp106 = w41*(B_0_2 + B_0_6);
                                        const double tmp107 = w37*(B_1_2 + B_1_3 + B_1_6 + B_1_7);
                                        const double tmp108 = w39*(B_0_3 + B_0_7);
                                        const double tmp109 = w45*(B_2_0 + B_2_1);
                                        const double tmp110 = w36*(B_0_0 + B_0_4);
                                        const double tmp111 = w40*(B_1_0 + B_1_1 + B_1_4 + B_1_5);
                                        const double tmp112 = w33*(B_0_1 + B_0_5);
                                        const double tmp113 = w46*(B_2_2 + B_2_3);
                                        const double tmp114 = w42*(-B_1_0 - B_1_4);
                                        const double tmp115 = w41*(-B_0_0 - B_0_4);
                                        const double tmp116 = w37*(-B_1_2 - B_1_6);
                                        const double tmp117 = w39*(-B_0_1 - B_0_5);
                                        const double tmp118 = w36*(-B_0_2 - B_0_6);
                                        const double tmp119 = w40*(-B_1_1 - B_1_5);
                                        const double tmp120 = w33*(-B_0_3 - B_0_7);
                                        const double tmp121 = w34*(-B_1_3 - B_1_7);
                                        const double tmp122 = w38*(B_2_2 + B_2_3);
                                        const double tmp123 = w42*(B_1_6 + B_1_7);
                                        const double tmp124 = w35*(B_2_0 + B_2_1);
                                        const double tmp125 = w37*(B_1_4 + B_1_5);
                                        const double tmp126 = w39*(-B_0_3 - B_0_5);
                                        const double tmp127 = w45*(B_2_4 + B_2_5);
                                        const double tmp128 = w36*(-B_0_2 - B_0_4);
                                        const double tmp129 = w40*(B_1_2 + B_1_3);
                                        const double tmp130 = w46*(B_2_6 + B_2_7);
                                        const double tmp131 = w34*(B_1_0 + B_1_1);
                                        const double tmp132 = w38*(-B_2_1 - B_2_2);
                                        const double tmp133 = w37*(B_1_2 + B_1_7);
                                        const double tmp134 = w39*(B_0_1 + B_0_7);
                                        const double tmp135 = w36*(B_0_0 + B_0_6);
                                        const double tmp136 = w40*(B_1_0 + B_1_5);
                                        const double tmp137 = w45*(-B_2_5 - B_2_6);
                                        const double tmp138 = w38*(-B_2_4 - B_2_6);
                                        const double tmp139 = w35*(-B_2_5 - B_2_7);
                                        const double tmp140 = w41*(-B_0_0 - B_0_2);
                                        const double tmp141 = w37*(B_1_1 + B_1_4);
                                        const double tmp142 = w39*(-B_0_1 - B_0_3);
                                        const double tmp143 = w45*(-B_2_1 - B_2_3);
                                        const double tmp144 = w36*(-B_0_4 - B_0_6);
                                        const double tmp145 = w40*(B_1_3 + B_1_6);
                                        const double tmp146 = w33*(-B_0_5 - B_0_7);
                                        const double tmp147 = w46*(-B_2_0 - B_2_2);
                                        const double tmp148 = w39*(B_0_2 + B_0_4);
                                        const double tmp149 = w36*(B_0_3 + B_0_5);
                                        const double tmp150 = w38*(B_2_5 + B_2_6);
                                        const double tmp151 = w37*(-B_1_0 - B_1_5);
                                        const double tmp152 = w39*(-B_0_0 - B_0_6);
                                        const double tmp153 = w45*(B_2_1 + B_2_2);
                                        const double tmp154 = w36*(-B_0_1 - B_0_7);
                                        const double tmp155 = w40*(-B_1_2 - B_1_7);
                                        const double tmp156 = w41*(-B_0_3 - B_0_7);
                                        const double tmp157 = w39*(-B_0_2 - B_0_6);
                                        const double tmp158 = w36*(-B_0_1 - B_0_5);
                                        const double tmp159 = w33*(-B_0_0 - B_0_4);
                                        const double tmp160 = w38*(-B_2_2 - B_2_3);
                                        const double tmp161 = w35*(-B_2_0 - B_2_1);
                                        const double tmp162 = w45*(-B_2_4 - B_2_5);
                                        const double tmp163 = w46*(-B_2_6 - B_2_7);
                                        const double tmp164 = w38*(-B_2_0 - B_2_3);
                                        const double tmp165 = w37*(B_1_3 + B_1_6);
                                        const double tmp166 = w40*(B_1_1 + B_1_4);
                                        const double tmp167 = w45*(-B_2_4 - B_2_7);
                                        const double tmp168 = w39*(B_0_3 + B_0_5);
                                        const double tmp169 = w36*(B_0_2 + B_0_4);
                                        const double tmp170 = w38*(B_2_1 + B_2_3);
                                        const double tmp171 = w35*(B_2_0 + B_2_2);
                                        const double tmp172 = w41*(B_0_5 + B_0_7);
                                        const double tmp173 = w37*(-B_1_3 - B_1_6);
                                        const double tmp174 = w39*(B_0_4 + B_0_6);
                                        const double tmp175 = w45*(B_2_4 + B_2_6);
                                        const double tmp176 = w36*(B_0_1 + B_0_3);
                                        const double tmp177 = w40*(-B_1_1 - B_1_4);
                                        const double tmp178 = w33*(B_0_0 + B_0_2);
                                        const double tmp179 = w46*(B_2_5 + B_2_7);
                                        const double tmp180 = w38*(B_2_5 + B_2_7);
                                        const double tmp181 = w42*(-B_1_3 - B_1_7);
                                        const double tmp182 = w35*(B_2_4 + B_2_6);
                                        const double tmp183 = w37*(-B_1_1 - B_1_5);
                                        const double tmp184 = w39*(B_0_1 + B_0_3 + B_0_5 + B_0_7);
                                        const double tmp185 = w45*(B_2_0 + B_2_2);
                                        const double tmp186 = w36*(B_0_0 + B_0_2 + B_0_4 + B_0_6);
                                        const double tmp187 = w40*(-B_1_2 - B_1_6);
                                        const double tmp188 = w46*(B_2_1 + B_2_3);
                                        const double tmp189 = w34*(-B_1_0 - B_1_4);
                                        const double tmp190 = w38*(B_2_4 + B_2_5);
                                        const double tmp191 = w35*(B_2_6 + B_2_7);
                                        const double tmp192 = w41*(-B_0_1 - B_0_5);
                                        const double tmp193 = w37*(-B_1_0 - B_1_1 - B_1_4 - B_1_5);
                                        const double tmp194 = w39*(-B_0_0 - B_0_4);
                                        const double tmp195 = w45*(B_2_2 + B_2_3);
                                        const double tmp196 = w36*(-B_0_3 - B_0_7);
                                        const double tmp197 = w40*(-B_1_2 - B_1_3 - B_1_6 - B_1_7);
                                        const double tmp198 = w33*(-B_0_2 - B_0_6);
                                        const double tmp199 = w46*(B_2_0 + B_2_1);
                                        const double tmp200 = w38*(-B_2_6 - B_2_7);
                                        const double tmp201 = w42*(B_1_2 + B_1_3);
                                        const double tmp202 = w35*(-B_2_4 - B_2_5);
                                        const double tmp203 = w37*(B_1_0 + B_1_1);
                                        const double tmp204 = w45*(-B_2_0 - B_2_1);
                                        const double tmp205 = w40*(B_1_6 + B_1_7);
                                        const double tmp206 = w46*(-B_2_2 - B_2_3);
                                        const double tmp207 = w34*(B_1_4 + B_1_5);
                                        const double tmp208 = w37*(-B_1_1 - B_1_4);
                                        const double tmp209 = w39*(-B_0_2 - B_0_4);
                                        const double tmp210 = w36*(-B_0_3 - B_0_5);
                                        const double tmp211 = w40*(-B_1_3 - B_1_6);
                                        const double tmp212 = w38*(B_2_4 + B_2_7);
                                        const double tmp213 = w45*(B_2_0 + B_2_3);
                                        const double tmp214 = w41*(B_0_0 + B_0_4);
                                        const double tmp215 = w39*(B_0_1 + B_0_5);
                                        const double tmp216 = w36*(B_0_2 + B_0_6);
                                        const double tmp217 = w33*(B_0_3 + B_0_7);
                                        const double tmp218 = w42*(B_1_1 + B_1_5);
                                        const double tmp219 = w37*(B_1_3 + B_1_7);
                                        const double tmp220 = w40*(B_1_0 + B_1_4);
                                        const double tmp221 = w34*(B_1_2 + B_1_6);
                                        const double tmp222 = w39*(-B_0_1 - B_0_7);
                                        const double tmp223 = w36*(-B_0_0 - B_0_6);
                                        const double tmp224 = w38*(-B_2_0 - B_2_1);
                                        const double tmp225 = w35*(-B_2_2 - B_2_3);
                                        const double tmp226 = w45*(-B_2_6 - B_2_7);
                                        const double tmp227 = w46*(-B_2_4 - B_2_5);
                                        const double tmp228 = w38*(B_2_4 + B_2_6);
                                        const double tmp229 = w42*(B_1_0 + B_1_4);
                                        const double tmp230 = w35*(B_2_5 + B_2_7);
                                        const double tmp231 = w37*(B_1_2 + B_1_6);
                                        const double tmp232 = w39*(-B_0_0 - B_0_2 - B_0_4 - B_0_6);
                                        const double tmp233 = w45*(B_2_1 + B_2_3);
                                        const double tmp234 = w36*(-B_0_1 - B_0_3 - B_0_5 - B_0_7);
                                        const double tmp235 = w40*(B_1_1 + B_1_5);
                                        const double tmp236 = w46*(B_2_0 + B_2_2);
                                        const double tmp237 = w34*(B_1_3 + B_1_7);
                                        const double tmp238 = w42*(-B_1_2 - B_1_6);
                                        const double tmp239 = w37*(-B_1_0 - B_1_4);
                                        const double tmp240 = w40*(-B_1_3 - B_1_7);
                                        const double tmp241 = w34*(-B_1_1 - B_1_5);
                                        const double tmp242 = w38*(-B_2_4 - B_2_5);
                                        const double tmp243 = w42*(-B_1_0 - B_1_1);
                                        const double tmp244 = w35*(-B_2_6 - B_2_7);
                                        const double tmp245 = w37*(-B_1_2 - B_1_3);
                                        const double tmp246 = w45*(-B_2_2 - B_2_3);
                                        const double tmp247 = w40*(-B_1_4 - B_1_5);
                                        const double tmp248 = w46*(-B_2_0 - B_2_1);
                                        const double tmp249 = w34*(-B_1_6 - B_1_7);
                                        const double tmp250 = w42*(B_1_4 + B_1_5);
                                        const double tmp251 = w37*(B_1_6 + B_1_7);
                                        const double tmp252 = w40*(B_1_0 + B_1_1);
                                        const double tmp253 = w34*(B_1_2 + B_1_3);
                                        const double tmp254 = w38*(-B_2_1 - B_2_3);
                                        const double tmp255 = w35*(-B_2_0 - B_2_2);
                                        const double tmp256 = w45*(-B_2_4 - B_2_6);
                                        const double tmp257 = w46*(-B_2_5 - B_2_7);
                                        const double tmp258 = w38*(B_2_0 + B_2_1 + B_2_2 + B_2_3);
                                        const double tmp259 = w45*(B_2_4 + B_2_5 + B_2_6 + B_2_7);
                                        const double tmp260 = w38*(-B_2_0 - B_2_2);
                                        const double tmp261 = w35*(-B_2_1 - B_2_3);
                                        const double tmp262 = w45*(-B_2_5 - B_2_7);
                                        const double tmp263 = w46*(-B_2_4 - B_2_6);
                                        EM_S[INDEX4(k,m,0,0,numEq,numComp,8)]+=-B_0_0*w50 - B_0_1*w41 - B_0_6*w33 - B_0_7*w49 + B_1_0*w47 - B_1_2*w42 - B_1_5*w34 + B_1_7*w48 - B_2_0*w43 - B_2_3*w35 - B_2_4*w46 - B_2_7*w44 + tmp132 + tmp137 + tmp208 + tmp209 + tmp210 + tmp211;
                                        EM_S[INDEX4(k,m,0,1,numEq,numComp,8)]+=-B_0_0*w41 - B_0_1*w50 - B_0_6*w49 - B_0_7*w33 + tmp126 + tmp128 + tmp242 + tmp243 + tmp244 + tmp245 + tmp246 + tmp247 + tmp248 + tmp249;
                                        EM_S[INDEX4(k,m,0,2,numEq,numComp,8)]+=-B_1_0*w42 + B_1_2*w47 + B_1_5*w48 - B_1_7*w34 + tmp138 + tmp139 + tmp140 + tmp142 + tmp143 + tmp144 + tmp146 + tmp147 + tmp173 + tmp177;
                                        EM_S[INDEX4(k,m,0,3,numEq,numComp,8)]+=tmp100 + tmp101 + tmp102 + tmp103 + tmp40 + tmp45 + tmp96 + tmp97 + tmp98 + tmp99;
                                        EM_S[INDEX4(k,m,0,4,numEq,numComp,8)]+=-B_2_0*w46 - B_2_3*w44 - B_2_4*w43 - B_2_7*w35 + tmp114 + tmp115 + tmp116 + tmp117 + tmp118 + tmp119 + tmp120 + tmp121 + tmp92 + tmp93;
                                        EM_S[INDEX4(k,m,0,5,numEq,numComp,8)]+=tmp192 + tmp193 + tmp194 + tmp196 + tmp197 + tmp198 + tmp224 + tmp225 + tmp226 + tmp227;
                                        EM_S[INDEX4(k,m,0,6,numEq,numComp,8)]+=tmp232 + tmp234 + tmp238 + tmp239 + tmp240 + tmp241 + tmp260 + tmp261 + tmp262 + tmp263;
                                        EM_S[INDEX4(k,m,0,7,numEq,numComp,8)]+=tmp60 + tmp61 + tmp62 + tmp63 + tmp64 + tmp65;
                                        EM_S[INDEX4(k,m,1,0,numEq,numComp,8)]+=B_0_0*w50 + B_0_1*w41 + B_0_6*w33 + B_0_7*w49 + tmp148 + tmp149 + tmp242 + tmp243 + tmp244 + tmp245 + tmp246 + tmp247 + tmp248 + tmp249;
                                        EM_S[INDEX4(k,m,1,1,numEq,numComp,8)]+=B_0_0*w41 + B_0_1*w50 + B_0_6*w49 + B_0_7*w33 + B_1_1*w47 - B_1_3*w42 - B_1_4*w34 + B_1_6*w48 - B_2_1*w43 - B_2_2*w35 - B_2_5*w46 - B_2_6*w44 + tmp151 + tmp155 + tmp164 + tmp167 + tmp168 + tmp169;
                                        EM_S[INDEX4(k,m,1,2,numEq,numComp,8)]+=tmp101 + tmp103 + tmp40 + tmp42 + tmp44 + tmp45 + tmp46 + tmp48 + tmp96 + tmp98;
                                        EM_S[INDEX4(k,m,1,3,numEq,numComp,8)]+=-B_1_1*w42 + B_1_3*w47 + B_1_4*w48 - B_1_6*w34 + tmp20 + tmp21 + tmp22 + tmp23 + tmp24 + tmp25 + tmp26 + tmp27 + tmp28 + tmp29;
                                        EM_S[INDEX4(k,m,1,4,numEq,numComp,8)]+=tmp193 + tmp197 + tmp214 + tmp215 + tmp216 + tmp217 + tmp224 + tmp225 + tmp226 + tmp227;
                                        EM_S[INDEX4(k,m,1,5,numEq,numComp,8)]+=-B_2_1*w46 - B_2_2*w44 - B_2_5*w43 - B_2_6*w35 + tmp70 + tmp75 + tmp83 + tmp84 + tmp85 + tmp86 + tmp88 + tmp89 + tmp90 + tmp91;
                                        EM_S[INDEX4(k,m,1,6,numEq,numComp,8)]+=tmp60 + tmp61 + tmp63 + tmp65 + tmp80 + tmp81;
                                        EM_S[INDEX4(k,m,1,7,numEq,numComp,8)]+=tmp181 + tmp183 + tmp184 + tmp186 + tmp187 + tmp189 + tmp254 + tmp255 + tmp256 + tmp257;
                                        EM_S[INDEX4(k,m,2,0,numEq,numComp,8)]+=-B_1_0*w47 + B_1_2*w42 + B_1_5*w34 - B_1_7*w48 + tmp138 + tmp139 + tmp140 + tmp141 + tmp142 + tmp143 + tmp144 + tmp145 + tmp146 + tmp147;
                                        EM_S[INDEX4(k,m,2,1,numEq,numComp,8)]+=tmp100 + tmp102 + tmp40 + tmp41 + tmp43 + tmp45 + tmp47 + tmp49 + tmp97 + tmp99;
                                        EM_S[INDEX4(k,m,2,2,numEq,numComp,8)]+=-B_0_2*w50 - B_0_3*w41 - B_0_4*w33 - B_0_5*w49 + B_1_0*w42 - B_1_2*w47 - B_1_5*w48 + B_1_7*w34 - B_2_1*w35 - B_2_2*w43 - B_2_5*w44 - B_2_6*w46 + tmp152 + tmp154 + tmp164 + tmp165 + tmp166 + tmp167;
                                        EM_S[INDEX4(k,m,2,3,numEq,numComp,8)]+=-B_0_2*w41 - B_0_3*w50 - B_0_4*w49 - B_0_5*w33 + tmp200 + tmp201 + tmp202 + tmp203 + tmp204 + tmp205 + tmp206 + tmp207 + tmp222 + tmp223;
                                        EM_S[INDEX4(k,m,2,4,numEq,numComp,8)]+=tmp229 + tmp231 + tmp232 + tmp234 + tmp235 + tmp237 + tmp260 + tmp261 + tmp262 + tmp263;
                                        EM_S[INDEX4(k,m,2,5,numEq,numComp,8)]+=tmp60 + tmp62 + tmp63 + tmp64 + tmp94 + tmp95;
                                        EM_S[INDEX4(k,m,2,6,numEq,numComp,8)]+=-B_2_1*w44 - B_2_2*w46 - B_2_5*w35 - B_2_6*w43 + tmp70 + tmp71 + tmp72 + tmp73 + tmp74 + tmp75 + tmp76 + tmp77 + tmp78 + tmp79;
                                        EM_S[INDEX4(k,m,2,7,numEq,numComp,8)]+=tmp107 + tmp111 + tmp156 + tmp157 + tmp158 + tmp159 + tmp160 + tmp161 + tmp162 + tmp163;
                                        EM_S[INDEX4(k,m,3,0,numEq,numComp,8)]+=tmp40 + tmp41 + tmp42 + tmp43 + tmp44 + tmp45 + tmp46 + tmp47 + tmp48 + tmp49;
                                        EM_S[INDEX4(k,m,3,1,numEq,numComp,8)]+=-B_1_1*w47 + B_1_3*w42 + B_1_4*w34 - B_1_6*w48 + tmp20 + tmp21 + tmp22 + tmp24 + tmp25 + tmp26 + tmp28 + tmp29 + tmp33 + tmp37;
                                        EM_S[INDEX4(k,m,3,2,numEq,numComp,8)]+=B_0_2*w50 + B_0_3*w41 + B_0_4*w33 + B_0_5*w49 + tmp200 + tmp201 + tmp202 + tmp203 + tmp204 + tmp205 + tmp206 + tmp207 + tmp54 + tmp56;
                                        EM_S[INDEX4(k,m,3,3,numEq,numComp,8)]+=B_0_2*w41 + B_0_3*w50 + B_0_4*w49 + B_0_5*w33 + B_1_1*w42 - B_1_3*w47 - B_1_4*w48 + B_1_6*w34 - B_2_0*w35 - B_2_3*w43 - B_2_4*w44 - B_2_7*w46 + tmp132 + tmp133 + tmp134 + tmp135 + tmp136 + tmp137;
                                        EM_S[INDEX4(k,m,3,4,numEq,numComp,8)]+=tmp60 + tmp63 + tmp80 + tmp81 + tmp94 + tmp95;
                                        EM_S[INDEX4(k,m,3,5,numEq,numComp,8)]+=tmp184 + tmp186 + tmp218 + tmp219 + tmp220 + tmp221 + tmp254 + tmp255 + tmp256 + tmp257;
                                        EM_S[INDEX4(k,m,3,6,numEq,numComp,8)]+=tmp106 + tmp107 + tmp108 + tmp110 + tmp111 + tmp112 + tmp160 + tmp161 + tmp162 + tmp163;
                                        EM_S[INDEX4(k,m,3,7,numEq,numComp,8)]+=-B_2_0*w44 - B_2_3*w46 - B_2_4*w35 - B_2_7*w43 + tmp1 + tmp2 + tmp3 + tmp4 + tmp6 + tmp7 + tmp8 + tmp9 + tmp92 + tmp93;
                                        EM_S[INDEX4(k,m,4,0,numEq,numComp,8)]+=B_2_0*w43 + B_2_3*w35 + B_2_4*w46 + B_2_7*w44 + tmp0 + tmp114 + tmp115 + tmp116 + tmp117 + tmp118 + tmp119 + tmp120 + tmp121 + tmp5;
                                        EM_S[INDEX4(k,m,4,1,numEq,numComp,8)]+=tmp190 + tmp191 + tmp192 + tmp193 + tmp194 + tmp195 + tmp196 + tmp197 + tmp198 + tmp199;
                                        EM_S[INDEX4(k,m,4,2,numEq,numComp,8)]+=tmp228 + tmp230 + tmp232 + tmp233 + tmp234 + tmp236 + tmp238 + tmp239 + tmp240 + tmp241;
                                        EM_S[INDEX4(k,m,4,3,numEq,numComp,8)]+=tmp258 + tmp259 + tmp61 + tmp62 + tmp64 + tmp65;
                                        EM_S[INDEX4(k,m,4,4,numEq,numComp,8)]+=-B_0_2*w33 - B_0_3*w49 - B_0_4*w50 - B_0_5*w41 - B_1_1*w34 + B_1_3*w48 + B_1_4*w47 - B_1_6*w42 + B_2_0*w46 + B_2_3*w44 + B_2_4*w43 + B_2_7*w35 + tmp150 + tmp151 + tmp152 + tmp153 + tmp154 + tmp155;
                                        EM_S[INDEX4(k,m,4,5,numEq,numComp,8)]+=-B_0_2*w49 - B_0_3*w33 - B_0_4*w41 - B_0_5*w50 + tmp222 + tmp223 + tmp50 + tmp51 + tmp52 + tmp53 + tmp55 + tmp57 + tmp58 + tmp59;
                                        EM_S[INDEX4(k,m,4,6,numEq,numComp,8)]+=B_1_1*w48 - B_1_3*w34 - B_1_4*w42 + B_1_6*w47 + tmp23 + tmp27 + tmp30 + tmp31 + tmp32 + tmp34 + tmp35 + tmp36 + tmp38 + tmp39;
                                        EM_S[INDEX4(k,m,4,7,numEq,numComp,8)]+=tmp10 + tmp11 + tmp12 + tmp13 + tmp14 + tmp15 + tmp16 + tmp17 + tmp18 + tmp19;
                                        EM_S[INDEX4(k,m,5,0,numEq,numComp,8)]+=tmp190 + tmp191 + tmp193 + tmp195 + tmp197 + tmp199 + tmp214 + tmp215 + tmp216 + tmp217;
                                        EM_S[INDEX4(k,m,5,1,numEq,numComp,8)]+=B_2_1*w43 + B_2_2*w35 + B_2_5*w46 + B_2_6*w44 + tmp82 + tmp83 + tmp84 + tmp85 + tmp86 + tmp87 + tmp88 + tmp89 + tmp90 + tmp91;
                                        EM_S[INDEX4(k,m,5,2,numEq,numComp,8)]+=tmp258 + tmp259 + tmp61 + tmp65 + tmp80 + tmp81;
                                        EM_S[INDEX4(k,m,5,3,numEq,numComp,8)]+=tmp180 + tmp181 + tmp182 + tmp183 + tmp184 + tmp185 + tmp186 + tmp187 + tmp188 + tmp189;
                                        EM_S[INDEX4(k,m,5,4,numEq,numComp,8)]+=B_0_2*w33 + B_0_3*w49 + B_0_4*w50 + B_0_5*w41 + tmp50 + tmp51 + tmp52 + tmp53 + tmp54 + tmp55 + tmp56 + tmp57 + tmp58 + tmp59;
                                        EM_S[INDEX4(k,m,5,5,numEq,numComp,8)]+=B_0_2*w49 + B_0_3*w33 + B_0_4*w41 + B_0_5*w50 - B_1_0*w34 + B_1_2*w48 + B_1_5*w47 - B_1_7*w42 + B_2_1*w46 + B_2_2*w44 + B_2_5*w43 + B_2_6*w35 + tmp134 + tmp135 + tmp208 + tmp211 + tmp212 + tmp213;
                                        EM_S[INDEX4(k,m,5,6,numEq,numComp,8)]+=tmp10 + tmp11 + tmp13 + tmp15 + tmp17 + tmp19 + tmp66 + tmp67 + tmp68 + tmp69;
                                        EM_S[INDEX4(k,m,5,7,numEq,numComp,8)]+=B_1_0*w48 - B_1_2*w34 - B_1_5*w42 + B_1_7*w47 + tmp170 + tmp171 + tmp172 + tmp173 + tmp174 + tmp175 + tmp176 + tmp177 + tmp178 + tmp179;
                                        EM_S[INDEX4(k,m,6,0,numEq,numComp,8)]+=tmp228 + tmp229 + tmp230 + tmp231 + tmp232 + tmp233 + tmp234 + tmp235 + tmp236 + tmp237;
                                        EM_S[INDEX4(k,m,6,1,numEq,numComp,8)]+=tmp258 + tmp259 + tmp62 + tmp64 + tmp94 + tmp95;
                                        EM_S[INDEX4(k,m,6,2,numEq,numComp,8)]+=B_2_1*w35 + B_2_2*w43 + B_2_5*w44 + B_2_6*w46 + tmp71 + tmp72 + tmp73 + tmp74 + tmp76 + tmp77 + tmp78 + tmp79 + tmp82 + tmp87;
                                        EM_S[INDEX4(k,m,6,3,numEq,numComp,8)]+=tmp104 + tmp105 + tmp107 + tmp109 + tmp111 + tmp113 + tmp156 + tmp157 + tmp158 + tmp159;
                                        EM_S[INDEX4(k,m,6,4,numEq,numComp,8)]+=B_1_1*w34 - B_1_3*w48 - B_1_4*w47 + B_1_6*w42 + tmp30 + tmp31 + tmp32 + tmp33 + tmp34 + tmp35 + tmp36 + tmp37 + tmp38 + tmp39;
                                        EM_S[INDEX4(k,m,6,5,numEq,numComp,8)]+=tmp10 + tmp12 + tmp14 + tmp15 + tmp16 + tmp18 + tmp250 + tmp251 + tmp252 + tmp253;
                                        EM_S[INDEX4(k,m,6,6,numEq,numComp,8)]+=-B_0_0*w33 - B_0_1*w49 - B_0_6*w50 - B_0_7*w41 - B_1_1*w48 + B_1_3*w34 + B_1_4*w42 - B_1_6*w47 + B_2_1*w44 + B_2_2*w46 + B_2_5*w35 + B_2_6*w43 + tmp133 + tmp136 + tmp209 + tmp210 + tmp212 + tmp213;
                                        EM_S[INDEX4(k,m,6,7,numEq,numComp,8)]+=-B_0_0*w49 - B_0_1*w33 - B_0_6*w41 - B_0_7*w50 + tmp122 + tmp123 + tmp124 + tmp125 + tmp126 + tmp127 + tmp128 + tmp129 + tmp130 + tmp131;
                                        EM_S[INDEX4(k,m,7,0,numEq,numComp,8)]+=tmp258 + tmp259 + tmp80 + tmp81 + tmp94 + tmp95;
                                        EM_S[INDEX4(k,m,7,1,numEq,numComp,8)]+=tmp180 + tmp182 + tmp184 + tmp185 + tmp186 + tmp188 + tmp218 + tmp219 + tmp220 + tmp221;
                                        EM_S[INDEX4(k,m,7,2,numEq,numComp,8)]+=tmp104 + tmp105 + tmp106 + tmp107 + tmp108 + tmp109 + tmp110 + tmp111 + tmp112 + tmp113;
                                        EM_S[INDEX4(k,m,7,3,numEq,numComp,8)]+=B_2_0*w35 + B_2_3*w43 + B_2_4*w44 + B_2_7*w46 + tmp0 + tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6 + tmp7 + tmp8 + tmp9;
                                        EM_S[INDEX4(k,m,7,4,numEq,numComp,8)]+=tmp10 + tmp15 + tmp250 + tmp251 + tmp252 + tmp253 + tmp66 + tmp67 + tmp68 + tmp69;
                                        EM_S[INDEX4(k,m,7,5,numEq,numComp,8)]+=B_1_0*w34 - B_1_2*w48 - B_1_5*w47 + B_1_7*w42 + tmp141 + tmp145 + tmp170 + tmp171 + tmp172 + tmp174 + tmp175 + tmp176 + tmp178 + tmp179;
                                        EM_S[INDEX4(k,m,7,6,numEq,numComp,8)]+=B_0_0*w33 + B_0_1*w49 + B_0_6*w50 + B_0_7*w41 + tmp122 + tmp123 + tmp124 + tmp125 + tmp127 + tmp129 + tmp130 + tmp131 + tmp148 + tmp149;
                                        EM_S[INDEX4(k,m,7,7,numEq,numComp,8)]+=B_0_0*w49 + B_0_1*w33 + B_0_6*w41 + B_0_7*w50 - B_1_0*w48 + B_1_2*w34 + B_1_5*w42 - B_1_7*w47 + B_2_0*w44 + B_2_3*w46 + B_2_4*w35 + B_2_7*w43 + tmp150 + tmp153 + tmp165 + tmp166 + tmp168 + tmp169;
                                    }
                                }
                            } else { // constant data
                                for (index_t k=0; k<numEq; k++) {
                                    for (index_t m=0; m<numComp; m++) {
                                        const double wB0 = B_p[INDEX3(k,0,m,numEq,3)]*w55;
                                        const double wB1 = B_p[INDEX3(k,1,m,numEq,3)]*w56;
                                        const double wB2 = B_p[INDEX3(k,2,m,numEq,3)]*w54;
                                        EM_S[INDEX4(k,m,0,0,numEq,numComp,8)]+= 4*wB0 + 4*wB1 + 4*wB2;
                                        EM_S[INDEX4(k,m,0,1,numEq,numComp,8)]+= 4*wB0 + 2*wB1 + 2*wB2;
                                        EM_S[INDEX4(k,m,0,2,numEq,numComp,8)]+= 2*wB0 + 4*wB1 + 2*wB2;
                                        EM_S[INDEX4(k,m,0,3,numEq,numComp,8)]+= 2*wB0 + 2*wB1 +   wB2;
                                        EM_S[INDEX4(k,m,0,4,numEq,numComp,8)]+= 2*wB0 + 2*wB1 + 4*wB2;
                                        EM_S[INDEX4(k,m,0,5,numEq,numComp,8)]+= 2*wB0 +   wB1 + 2*wB2;
                                        EM_S[INDEX4(k,m,0,6,numEq,numComp,8)]+=   wB0 + 2*wB1 + 2*wB2;
                                        EM_S[INDEX4(k,m,0,7,numEq,numComp,8)]+=   wB0 +   wB1 +   wB2;
                                        EM_S[INDEX4(k,m,1,0,numEq,numComp,8)]+=-4*wB0 + 2*wB1 + 2*wB2;
                                        EM_S[INDEX4(k,m,1,1,numEq,numComp,8)]+=-4*wB0 + 4*wB1 + 4*wB2;
                                        EM_S[INDEX4(k,m,1,2,numEq,numComp,8)]+=-2*wB0 + 2*wB1 +   wB2;
                                        EM_S[INDEX4(k,m,1,3,numEq,numComp,8)]+=-2*wB0 + 4*wB1 + 2*wB2;
                                        EM_S[INDEX4(k,m,1,4,numEq,numComp,8)]+=-2*wB0 +   wB1 + 2*wB2;
                                        EM_S[INDEX4(k,m,1,5,numEq,numComp,8)]+=-2*wB0 + 2*wB1 + 4*wB2;
                                        EM_S[INDEX4(k,m,1,6,numEq,numComp,8)]+=  -wB0 +   wB1 +   wB2;
                                        EM_S[INDEX4(k,m,1,7,numEq,numComp,8)]+=  -wB0 + 2*wB1 + 2*wB2;
                                        EM_S[INDEX4(k,m,2,0,numEq,numComp,8)]+= 2*wB0 - 4*wB1 + 2*wB2;
                                        EM_S[INDEX4(k,m,2,1,numEq,numComp,8)]+= 2*wB0 - 2*wB1 +   wB2;
                                        EM_S[INDEX4(k,m,2,2,numEq,numComp,8)]+= 4*wB0 - 4*wB1 + 4*wB2;
                                        EM_S[INDEX4(k,m,2,3,numEq,numComp,8)]+= 4*wB0 - 2*wB1 + 2*wB2;
                                        EM_S[INDEX4(k,m,2,4,numEq,numComp,8)]+=   wB0 - 2*wB1 + 2*wB2;
                                        EM_S[INDEX4(k,m,2,5,numEq,numComp,8)]+=   wB0 -   wB1 +   wB2;
                                        EM_S[INDEX4(k,m,2,6,numEq,numComp,8)]+= 2*wB0 - 2*wB1 + 4*wB2;
                                        EM_S[INDEX4(k,m,2,7,numEq,numComp,8)]+= 2*wB0 -   wB1 + 2*wB2;
                                        EM_S[INDEX4(k,m,3,0,numEq,numComp,8)]+=-2*wB0 - 2*wB1 +   wB2;
                                        EM_S[INDEX4(k,m,3,1,numEq,numComp,8)]+=-2*wB0 - 4*wB1 + 2*wB2;
                                        EM_S[INDEX4(k,m,3,2,numEq,numComp,8)]+=-4*wB0 - 2*wB1 + 2*wB2;
                                        EM_S[INDEX4(k,m,3,3,numEq,numComp,8)]+=-4*wB0 - 4*wB1 + 4*wB2;
                                        EM_S[INDEX4(k,m,3,4,numEq,numComp,8)]+=  -wB0 -   wB1 +   wB2;
                                        EM_S[INDEX4(k,m,3,5,numEq,numComp,8)]+=  -wB0 - 2*wB1 + 2*wB2;
                                        EM_S[INDEX4(k,m,3,6,numEq,numComp,8)]+=-2*wB0 -   wB1 + 2*wB2;
                                        EM_S[INDEX4(k,m,3,7,numEq,numComp,8)]+=-2*wB0 - 2*wB1 + 4*wB2;
                                        EM_S[INDEX4(k,m,4,0,numEq,numComp,8)]+= 2*wB0 + 2*wB1 - 4*wB2;
                                        EM_S[INDEX4(k,m,4,1,numEq,numComp,8)]+= 2*wB0 +   wB1 - 2*wB2;
                                        EM_S[INDEX4(k,m,4,2,numEq,numComp,8)]+=   wB0 + 2*wB1 - 2*wB2;
                                        EM_S[INDEX4(k,m,4,3,numEq,numComp,8)]+=   wB0 +   wB1 -   wB2;
                                        EM_S[INDEX4(k,m,4,4,numEq,numComp,8)]+= 4*wB0 + 4*wB1 - 4*wB2;
                                        EM_S[INDEX4(k,m,4,5,numEq,numComp,8)]+= 4*wB0 + 2*wB1 - 2*wB2;
                                        EM_S[INDEX4(k,m,4,6,numEq,numComp,8)]+= 2*wB0 + 4*wB1 - 2*wB2;
                                        EM_S[INDEX4(k,m,4,7,numEq,numComp,8)]+= 2*wB0 + 2*wB1 -   wB2;
                                        EM_S[INDEX4(k,m,5,0,numEq,numComp,8)]+=-2*wB0 +   wB1 - 2*wB2;
                                        EM_S[INDEX4(k,m,5,1,numEq,numComp,8)]+=-2*wB0 + 2*wB1 - 4*wB2;
                                        EM_S[INDEX4(k,m,5,2,numEq,numComp,8)]+=  -wB0 +   wB1 -   wB2;
                                        EM_S[INDEX4(k,m,5,3,numEq,numComp,8)]+=  -wB0 + 2*wB1 - 2*wB2;
                                        EM_S[INDEX4(k,m,5,4,numEq,numComp,8)]+=-4*wB0 + 2*wB1 - 2*wB2;
                                        EM_S[INDEX4(k,m,5,5,numEq,numComp,8)]+=-4*wB0 + 4*wB1 - 4*wB2;
                                        EM_S[INDEX4(k,m,5,6,numEq,numComp,8)]+=-2*wB0 + 2*wB1 -   wB2;
                                        EM_S[INDEX4(k,m,5,7,numEq,numComp,8)]+=-2*wB0 + 4*wB1 - 2*wB2;
                                        EM_S[INDEX4(k,m,6,0,numEq,numComp,8)]+=   wB0 - 2*wB1 - 2*wB2;
                                        EM_S[INDEX4(k,m,6,1,numEq,numComp,8)]+=   wB0 -   wB1 -   wB2;
                                        EM_S[INDEX4(k,m,6,2,numEq,numComp,8)]+= 2*wB0 - 2*wB1 - 4*wB2;
                                        EM_S[INDEX4(k,m,6,3,numEq,numComp,8)]+= 2*wB0 -   wB1 - 2*wB2;
                                        EM_S[INDEX4(k,m,6,4,numEq,numComp,8)]+= 2*wB0 - 4*wB1 - 2*wB2;
                                        EM_S[INDEX4(k,m,6,5,numEq,numComp,8)]+= 2*wB0 - 2*wB1 -   wB2;
                                        EM_S[INDEX4(k,m,6,6,numEq,numComp,8)]+= 4*wB0 - 4*wB1 - 4*wB2;
                                        EM_S[INDEX4(k,m,6,7,numEq,numComp,8)]+= 4*wB0 - 2*wB1 - 2*wB2;
                                        EM_S[INDEX4(k,m,7,0,numEq,numComp,8)]+=  -wB0 -   wB1 -   wB2;
                                        EM_S[INDEX4(k,m,7,1,numEq,numComp,8)]+=  -wB0 - 2*wB1 - 2*wB2;
                                        EM_S[INDEX4(k,m,7,2,numEq,numComp,8)]+=-2*wB0 -   wB1 - 2*wB2;
                                        EM_S[INDEX4(k,m,7,3,numEq,numComp,8)]+=-2*wB0 - 2*wB1 - 4*wB2;
                                        EM_S[INDEX4(k,m,7,4,numEq,numComp,8)]+=-2*wB0 - 2*wB1 -   wB2;
                                        EM_S[INDEX4(k,m,7,5,numEq,numComp,8)]+=-2*wB0 - 4*wB1 - 2*wB2;
                                        EM_S[INDEX4(k,m,7,6,numEq,numComp,8)]+=-4*wB0 - 2*wB1 - 2*wB2;
                                        EM_S[INDEX4(k,m,7,7,numEq,numComp,8)]+=-4*wB0 - 4*wB1 - 4*wB2;
                                    }
                                }
                            }
                        }
                        ///////////////
                        // process C //
                        ///////////////
                        if (!C.isEmpty()) {
                            const double* C_p=C.getSampleDataRO(e);
                            if (C.actsExpanded()) {
                                for (index_t k=0; k<numEq; k++) {
                                    for (index_t m=0; m<numComp; m++) {
                                        const double C_0_0 = C_p[INDEX4(k,m,0, 0, numEq,numComp,3)];
                                        const double C_1_0 = C_p[INDEX4(k,m,1, 0, numEq,numComp,3)];
                                        const double C_2_0 = C_p[INDEX4(k,m,2, 0, numEq,numComp,3)];
                                        const double C_0_1 = C_p[INDEX4(k,m,0, 1, numEq,numComp,3)];
                                        const double C_1_1 = C_p[INDEX4(k,m,1, 1, numEq,numComp,3)];
                                        const double C_2_1 = C_p[INDEX4(k,m,2, 1, numEq,numComp,3)];
                                        const double C_0_2 = C_p[INDEX4(k,m,0, 2, numEq,numComp,3)];
                                        const double C_1_2 = C_p[INDEX4(k,m,1, 2, numEq,numComp,3)];
                                        const double C_2_2 = C_p[INDEX4(k,m,2, 2, numEq,numComp,3)];
                                        const double C_0_3 = C_p[INDEX4(k,m,0, 3, numEq,numComp,3)];
                                        const double C_1_3 = C_p[INDEX4(k,m,1, 3, numEq,numComp,3)];
                                        const double C_2_3 = C_p[INDEX4(k,m,2, 3, numEq,numComp,3)];
                                        const double C_0_4 = C_p[INDEX4(k,m,0, 4, numEq,numComp,3)];
                                        const double C_1_4 = C_p[INDEX4(k,m,1, 4, numEq,numComp,3)];
                                        const double C_2_4 = C_p[INDEX4(k,m,2, 4, numEq,numComp,3)];
                                        const double C_0_5 = C_p[INDEX4(k,m,0, 5, numEq,numComp,3)];
                                        const double C_1_5 = C_p[INDEX4(k,m,1, 5, numEq,numComp,3)];
                                        const double C_2_5 = C_p[INDEX4(k,m,2, 5, numEq,numComp,3)];
                                        const double C_0_6 = C_p[INDEX4(k,m,0, 6, numEq,numComp,3)];
                                        const double C_1_6 = C_p[INDEX4(k,m,1, 6, numEq,numComp,3)];
                                        const double C_2_6 = C_p[INDEX4(k,m,2, 6, numEq,numComp,3)];
                                        const double C_0_7 = C_p[INDEX4(k,m,0, 7, numEq,numComp,3)];
                                        const double C_1_7 = C_p[INDEX4(k,m,1, 7, numEq,numComp,3)];
                                        const double C_2_7 = C_p[INDEX4(k,m,2, 7, numEq,numComp,3)];
                                        const double tmp0 = w38*(-C_2_5 - C_2_6);
                                        const double tmp1 = w42*(C_1_3 + C_1_7);
                                        const double tmp2 = w41*(C_0_3 + C_0_7);
                                        const double tmp3 = w37*(C_1_1 + C_1_5);
                                        const double tmp4 = w39*(C_0_2 + C_0_6);
                                        const double tmp5 = w45*(-C_2_1 - C_2_2);
                                        const double tmp6 = w36*(C_0_1 + C_0_5);
                                        const double tmp7 = w40*(C_1_2 + C_1_6);
                                        const double tmp8 = w33*(C_0_0 + C_0_4);
                                        const double tmp9 = w34*(C_1_0 + C_1_4);
                                        const double tmp10 = w38*(C_2_4 + C_2_5 + C_2_6 + C_2_7);
                                        const double tmp11 = w42*(C_1_4 + C_1_5);
                                        const double tmp12 = w41*(C_0_4 + C_0_6);
                                        const double tmp13 = w37*(C_1_6 + C_1_7);
                                        const double tmp14 = w39*(C_0_5 + C_0_7);
                                        const double tmp15 = w45*(C_2_0 + C_2_1 + C_2_2 + C_2_3);
                                        const double tmp16 = w36*(C_0_0 + C_0_2);
                                        const double tmp17 = w40*(C_1_0 + C_1_1);
                                        const double tmp18 = w33*(C_0_1 + C_0_3);
                                        const double tmp19 = w34*(C_1_2 + C_1_3);
                                        const double tmp20 = w38*(-C_2_5 - C_2_7);
                                        const double tmp21 = w35*(-C_2_4 - C_2_6);
                                        const double tmp22 = w41*(C_0_1 + C_0_3);
                                        const double tmp23 = w37*(C_1_0 + C_1_5);
                                        const double tmp24 = w39*(C_0_0 + C_0_2);
                                        const double tmp25 = w45*(-C_2_0 - C_2_2);
                                        const double tmp26 = w36*(C_0_5 + C_0_7);
                                        const double tmp27 = w40*(C_1_2 + C_1_7);
                                        const double tmp28 = w33*(C_0_4 + C_0_6);
                                        const double tmp29 = w46*(-C_2_1 - C_2_3);
                                        const double tmp30 = w38*(C_2_0 + C_2_2);
                                        const double tmp31 = w35*(C_2_1 + C_2_3);
                                        const double tmp32 = w41*(-C_0_4 - C_0_6);
                                        const double tmp33 = w37*(-C_1_2 - C_1_7);
                                        const double tmp34 = w39*(-C_0_5 - C_0_7);
                                        const double tmp35 = w45*(C_2_5 + C_2_7);
                                        const double tmp36 = w36*(-C_0_0 - C_0_2);
                                        const double tmp37 = w40*(-C_1_0 - C_1_5);
                                        const double tmp38 = w33*(-C_0_1 - C_0_3);
                                        const double tmp39 = w46*(C_2_4 + C_2_6);
                                        const double tmp40 = w38*(-C_2_0 - C_2_1 - C_2_2 - C_2_3);
                                        const double tmp41 = w42*(-C_1_2 - C_1_3);
                                        const double tmp42 = w41*(-C_0_1 - C_0_3);
                                        const double tmp43 = w37*(-C_1_0 - C_1_1);
                                        const double tmp44 = w39*(-C_0_0 - C_0_2);
                                        const double tmp45 = w45*(-C_2_4 - C_2_5 - C_2_6 - C_2_7);
                                        const double tmp46 = w36*(-C_0_5 - C_0_7);
                                        const double tmp47 = w40*(-C_1_6 - C_1_7);
                                        const double tmp48 = w33*(-C_0_4 - C_0_6);
                                        const double tmp49 = w34*(-C_1_4 - C_1_5);
                                        const double tmp50 = w38*(C_2_0 + C_2_1);
                                        const double tmp51 = w42*(-C_1_4 - C_1_5);
                                        const double tmp52 = w35*(C_2_2 + C_2_3);
                                        const double tmp53 = w37*(-C_1_6 - C_1_7);
                                        const double tmp54 = w39*(-C_0_1 - C_0_7);
                                        const double tmp55 = w45*(C_2_6 + C_2_7);
                                        const double tmp56 = w36*(-C_0_0 - C_0_6);
                                        const double tmp57 = w40*(-C_1_0 - C_1_1);
                                        const double tmp58 = w46*(C_2_4 + C_2_5);
                                        const double tmp59 = w34*(-C_1_2 - C_1_3);
                                        const double tmp60 = w38*(C_2_0 + C_2_1 + C_2_2 + C_2_3);
                                        const double tmp61 = w37*(C_1_0 + C_1_1 + C_1_4 + C_1_5);
                                        const double tmp62 = w39*(C_0_0 + C_0_2 + C_0_4 + C_0_6);
                                        const double tmp63 = w45*(C_2_4 + C_2_5 + C_2_6 + C_2_7);
                                        const double tmp64 = w36*(C_0_1 + C_0_3 + C_0_5 + C_0_7);
                                        const double tmp65 = w40*(C_1_2 + C_1_3 + C_1_6 + C_1_7);
                                        const double tmp66 = w41*(-C_0_5 - C_0_7);
                                        const double tmp67 = w39*(-C_0_4 - C_0_6);
                                        const double tmp68 = w36*(-C_0_1 - C_0_3);
                                        const double tmp69 = w33*(-C_0_0 - C_0_2);
                                        const double tmp70 = w38*(C_2_0 + C_2_3);
                                        const double tmp71 = w42*(C_1_2 + C_1_6);
                                        const double tmp72 = w41*(-C_0_2 - C_0_6);
                                        const double tmp73 = w37*(C_1_0 + C_1_4);
                                        const double tmp74 = w39*(-C_0_3 - C_0_7);
                                        const double tmp75 = w45*(C_2_4 + C_2_7);
                                        const double tmp76 = w36*(-C_0_0 - C_0_4);
                                        const double tmp77 = w40*(C_1_3 + C_1_7);
                                        const double tmp78 = w33*(-C_0_1 - C_0_5);
                                        const double tmp79 = w34*(C_1_1 + C_1_5);
                                        const double tmp80 = w39*(-C_0_1 - C_0_3 - C_0_5 - C_0_7);
                                        const double tmp81 = w36*(-C_0_0 - C_0_2 - C_0_4 - C_0_6);
                                        const double tmp82 = w38*(-C_2_4 - C_2_7);
                                        const double tmp83 = w42*(-C_1_1 - C_1_5);
                                        const double tmp84 = w41*(C_0_1 + C_0_5);
                                        const double tmp85 = w37*(-C_1_3 - C_1_7);
                                        const double tmp86 = w39*(C_0_0 + C_0_4);
                                        const double tmp87 = w45*(-C_2_0 - C_2_3);
                                        const double tmp88 = w36*(C_0_3 + C_0_7);
                                        const double tmp89 = w40*(-C_1_0 - C_1_4);
                                        const double tmp90 = w33*(C_0_2 + C_0_6);
                                        const double tmp91 = w34*(-C_1_2 - C_1_6);
                                        const double tmp92 = w38*(C_2_1 + C_2_2);
                                        const double tmp93 = w45*(C_2_5 + C_2_6);
                                        const double tmp94 = w37*(-C_1_2 - C_1_3 - C_1_6 - C_1_7);
                                        const double tmp95 = w40*(-C_1_0 - C_1_1 - C_1_4 - C_1_5);
                                        const double tmp96 = w42*(C_1_0 + C_1_1);
                                        const double tmp97 = w41*(C_0_0 + C_0_2);
                                        const double tmp98 = w37*(C_1_2 + C_1_3);
                                        const double tmp99 = w39*(C_0_1 + C_0_3);
                                        const double tmp100 = w36*(C_0_4 + C_0_6);
                                        const double tmp101 = w40*(C_1_4 + C_1_5);
                                        const double tmp102 = w33*(C_0_5 + C_0_7);
                                        const double tmp103 = w34*(C_1_6 + C_1_7);
                                        const double tmp104 = w38*(-C_2_2 - C_2_3);
                                        const double tmp105 = w35*(-C_2_0 - C_2_1);
                                        const double tmp106 = w41*(-C_0_3 - C_0_7);
                                        const double tmp107 = w37*(C_1_2 + C_1_3 + C_1_6 + C_1_7);
                                        const double tmp108 = w39*(-C_0_2 - C_0_6);
                                        const double tmp109 = w45*(-C_2_4 - C_2_5);
                                        const double tmp110 = w36*(-C_0_1 - C_0_5);
                                        const double tmp111 = w40*(C_1_0 + C_1_1 + C_1_4 + C_1_5);
                                        const double tmp112 = w33*(-C_0_0 - C_0_4);
                                        const double tmp113 = w46*(-C_2_6 - C_2_7);
                                        const double tmp114 = w42*(-C_1_0 - C_1_4);
                                        const double tmp115 = w41*(-C_0_0 - C_0_4);
                                        const double tmp116 = w37*(-C_1_2 - C_1_6);
                                        const double tmp117 = w39*(-C_0_1 - C_0_5);
                                        const double tmp118 = w36*(-C_0_2 - C_0_6);
                                        const double tmp119 = w40*(-C_1_1 - C_1_5);
                                        const double tmp120 = w33*(-C_0_3 - C_0_7);
                                        const double tmp121 = w34*(-C_1_3 - C_1_7);
                                        const double tmp122 = w38*(C_2_2 + C_2_3);
                                        const double tmp123 = w42*(C_1_6 + C_1_7);
                                        const double tmp124 = w35*(C_2_0 + C_2_1);
                                        const double tmp125 = w37*(C_1_4 + C_1_5);
                                        const double tmp126 = w39*(C_0_2 + C_0_4);
                                        const double tmp127 = w45*(C_2_4 + C_2_5);
                                        const double tmp128 = w36*(C_0_3 + C_0_5);
                                        const double tmp129 = w40*(C_1_2 + C_1_3);
                                        const double tmp130 = w46*(C_2_6 + C_2_7);
                                        const double tmp131 = w34*(C_1_0 + C_1_1);
                                        const double tmp132 = w38*(-C_2_1 - C_2_2);
                                        const double tmp133 = w37*(C_1_2 + C_1_7);
                                        const double tmp134 = w39*(C_0_1 + C_0_7);
                                        const double tmp135 = w36*(C_0_0 + C_0_6);
                                        const double tmp136 = w40*(C_1_0 + C_1_5);
                                        const double tmp137 = w45*(-C_2_5 - C_2_6);
                                        const double tmp138 = w38*(-C_2_4 - C_2_6);
                                        const double tmp139 = w35*(-C_2_5 - C_2_7);
                                        const double tmp140 = w41*(-C_0_0 - C_0_2);
                                        const double tmp141 = w37*(-C_1_3 - C_1_6);
                                        const double tmp142 = w39*(-C_0_1 - C_0_3);
                                        const double tmp143 = w45*(-C_2_1 - C_2_3);
                                        const double tmp144 = w36*(-C_0_4 - C_0_6);
                                        const double tmp145 = w40*(-C_1_1 - C_1_4);
                                        const double tmp146 = w33*(-C_0_5 - C_0_7);
                                        const double tmp147 = w46*(-C_2_0 - C_2_2);
                                        const double tmp148 = w39*(-C_0_3 - C_0_5);
                                        const double tmp149 = w36*(-C_0_2 - C_0_4);
                                        const double tmp150 = w38*(C_2_5 + C_2_6);
                                        const double tmp151 = w37*(-C_1_0 - C_1_5);
                                        const double tmp152 = w39*(-C_0_0 - C_0_6);
                                        const double tmp153 = w45*(C_2_1 + C_2_2);
                                        const double tmp154 = w36*(-C_0_1 - C_0_7);
                                        const double tmp155 = w40*(-C_1_2 - C_1_7);
                                        const double tmp156 = w41*(C_0_2 + C_0_6);
                                        const double tmp157 = w39*(C_0_3 + C_0_7);
                                        const double tmp158 = w36*(C_0_0 + C_0_4);
                                        const double tmp159 = w33*(C_0_1 + C_0_5);
                                        const double tmp160 = w38*(C_2_6 + C_2_7);
                                        const double tmp161 = w35*(C_2_4 + C_2_5);
                                        const double tmp162 = w45*(C_2_0 + C_2_1);
                                        const double tmp163 = w46*(C_2_2 + C_2_3);
                                        const double tmp164 = w38*(-C_2_0 - C_2_3);
                                        const double tmp165 = w37*(C_1_3 + C_1_6);
                                        const double tmp166 = w40*(C_1_1 + C_1_4);
                                        const double tmp167 = w45*(-C_2_4 - C_2_7);
                                        const double tmp168 = w39*(C_0_3 + C_0_5);
                                        const double tmp169 = w36*(C_0_2 + C_0_4);
                                        const double tmp170 = w38*(C_2_1 + C_2_3);
                                        const double tmp171 = w35*(C_2_0 + C_2_2);
                                        const double tmp172 = w41*(C_0_5 + C_0_7);
                                        const double tmp173 = w37*(C_1_1 + C_1_4);
                                        const double tmp174 = w39*(C_0_4 + C_0_6);
                                        const double tmp175 = w45*(C_2_4 + C_2_6);
                                        const double tmp176 = w36*(C_0_1 + C_0_3);
                                        const double tmp177 = w40*(C_1_3 + C_1_6);
                                        const double tmp178 = w33*(C_0_0 + C_0_2);
                                        const double tmp179 = w46*(C_2_5 + C_2_7);
                                        const double tmp180 = w38*(-C_2_1 - C_2_3);
                                        const double tmp181 = w42*(C_1_1 + C_1_5);
                                        const double tmp182 = w35*(-C_2_0 - C_2_2);
                                        const double tmp183 = w37*(C_1_3 + C_1_7);
                                        const double tmp184 = w39*(C_0_1 + C_0_3 + C_0_5 + C_0_7);
                                        const double tmp185 = w45*(-C_2_4 - C_2_6);
                                        const double tmp186 = w36*(C_0_0 + C_0_2 + C_0_4 + C_0_6);
                                        const double tmp187 = w40*(C_1_0 + C_1_4);
                                        const double tmp188 = w46*(-C_2_5 - C_2_7);
                                        const double tmp189 = w34*(C_1_2 + C_1_6);
                                        const double tmp190 = w38*(-C_2_0 - C_2_1);
                                        const double tmp191 = w35*(-C_2_2 - C_2_3);
                                        const double tmp192 = w41*(C_0_0 + C_0_4);
                                        const double tmp193 = w37*(-C_1_0 - C_1_1 - C_1_4 - C_1_5);
                                        const double tmp194 = w39*(C_0_1 + C_0_5);
                                        const double tmp195 = w45*(-C_2_6 - C_2_7);
                                        const double tmp196 = w36*(C_0_2 + C_0_6);
                                        const double tmp197 = w40*(-C_1_2 - C_1_3 - C_1_6 - C_1_7);
                                        const double tmp198 = w33*(C_0_3 + C_0_7);
                                        const double tmp199 = w46*(-C_2_4 - C_2_5);
                                        const double tmp200 = w38*(-C_2_6 - C_2_7);
                                        const double tmp201 = w42*(C_1_2 + C_1_3);
                                        const double tmp202 = w35*(-C_2_4 - C_2_5);
                                        const double tmp203 = w37*(C_1_0 + C_1_1);
                                        const double tmp204 = w45*(-C_2_0 - C_2_1);
                                        const double tmp205 = w40*(C_1_6 + C_1_7);
                                        const double tmp206 = w46*(-C_2_2 - C_2_3);
                                        const double tmp207 = w34*(C_1_4 + C_1_5);
                                        const double tmp208 = w37*(-C_1_1 - C_1_4);
                                        const double tmp209 = w39*(-C_0_2 - C_0_4);
                                        const double tmp210 = w36*(-C_0_3 - C_0_5);
                                        const double tmp211 = w40*(-C_1_3 - C_1_6);
                                        const double tmp212 = w38*(C_2_4 + C_2_7);
                                        const double tmp213 = w45*(C_2_0 + C_2_3);
                                        const double tmp214 = w41*(-C_0_1 - C_0_5);
                                        const double tmp215 = w39*(-C_0_0 - C_0_4);
                                        const double tmp216 = w36*(-C_0_3 - C_0_7);
                                        const double tmp217 = w33*(-C_0_2 - C_0_6);
                                        const double tmp218 = w42*(-C_1_3 - C_1_7);
                                        const double tmp219 = w37*(-C_1_1 - C_1_5);
                                        const double tmp220 = w40*(-C_1_2 - C_1_6);
                                        const double tmp221 = w34*(-C_1_0 - C_1_4);
                                        const double tmp222 = w39*(C_0_0 + C_0_6);
                                        const double tmp223 = w36*(C_0_1 + C_0_7);
                                        const double tmp224 = w38*(C_2_4 + C_2_5);
                                        const double tmp225 = w35*(C_2_6 + C_2_7);
                                        const double tmp226 = w45*(C_2_2 + C_2_3);
                                        const double tmp227 = w46*(C_2_0 + C_2_1);
                                        const double tmp228 = w38*(-C_2_0 - C_2_2);
                                        const double tmp229 = w42*(-C_1_2 - C_1_6);
                                        const double tmp230 = w35*(-C_2_1 - C_2_3);
                                        const double tmp231 = w37*(-C_1_0 - C_1_4);
                                        const double tmp232 = w39*(-C_0_0 - C_0_2 - C_0_4 - C_0_6);
                                        const double tmp233 = w45*(-C_2_5 - C_2_7);
                                        const double tmp234 = w36*(-C_0_1 - C_0_3 - C_0_5 - C_0_7);
                                        const double tmp235 = w40*(-C_1_3 - C_1_7);
                                        const double tmp236 = w46*(-C_2_4 - C_2_6);
                                        const double tmp237 = w34*(-C_1_1 - C_1_5);
                                        const double tmp238 = w42*(C_1_0 + C_1_4);
                                        const double tmp239 = w37*(C_1_2 + C_1_6);
                                        const double tmp240 = w40*(C_1_1 + C_1_5);
                                        const double tmp241 = w34*(C_1_3 + C_1_7);
                                        const double tmp242 = w38*(-C_2_4 - C_2_5);
                                        const double tmp243 = w42*(-C_1_0 - C_1_1);
                                        const double tmp244 = w35*(-C_2_6 - C_2_7);
                                        const double tmp245 = w37*(-C_1_2 - C_1_3);
                                        const double tmp246 = w45*(-C_2_2 - C_2_3);
                                        const double tmp247 = w40*(-C_1_4 - C_1_5);
                                        const double tmp248 = w46*(-C_2_0 - C_2_1);
                                        const double tmp249 = w34*(-C_1_6 - C_1_7);
                                        const double tmp250 = w42*(-C_1_6 - C_1_7);
                                        const double tmp251 = w37*(-C_1_4 - C_1_5);
                                        const double tmp252 = w40*(-C_1_2 - C_1_3);
                                        const double tmp253 = w34*(-C_1_0 - C_1_1);
                                        const double tmp254 = w38*(C_2_5 + C_2_7);
                                        const double tmp255 = w35*(C_2_4 + C_2_6);
                                        const double tmp256 = w45*(C_2_0 + C_2_2);
                                        const double tmp257 = w46*(C_2_1 + C_2_3);
                                        const double tmp258 = w38*(-C_2_4 - C_2_5 - C_2_6 - C_2_7);
                                        const double tmp259 = w45*(-C_2_0 - C_2_1 - C_2_2 - C_2_3);
                                        const double tmp260 = w38*(C_2_4 + C_2_6);
                                        const double tmp261 = w35*(C_2_5 + C_2_7);
                                        const double tmp262 = w45*(C_2_1 + C_2_3);
                                        const double tmp263 = w46*(C_2_0 + C_2_2);
                                        EM_S[INDEX4(k,m,0,0,numEq,numComp,8)]+=-C_0_0*w50 - C_0_1*w41 - C_0_6*w33 - C_0_7*w49 + C_1_0*w47 - C_1_2*w42 - C_1_5*w34 + C_1_7*w48 - C_2_0*w43 - C_2_3*w35 - C_2_4*w46 - C_2_7*w44 + tmp132 + tmp137 + tmp208 + tmp209 + tmp210 + tmp211;
                                        EM_S[INDEX4(k,m,0,1,numEq,numComp,8)]+=C_0_0*w50 + C_0_1*w41 + C_0_6*w33 + C_0_7*w49 + tmp126 + tmp128 + tmp242 + tmp243 + tmp244 + tmp245 + tmp246 + tmp247 + tmp248 + tmp249;
                                        EM_S[INDEX4(k,m,0,2,numEq,numComp,8)]+=-C_1_0*w47 + C_1_2*w42 + C_1_5*w34 - C_1_7*w48 + tmp138 + tmp139 + tmp140 + tmp142 + tmp143 + tmp144 + tmp146 + tmp147 + tmp173 + tmp177;
                                        EM_S[INDEX4(k,m,0,3,numEq,numComp,8)]+=tmp100 + tmp101 + tmp102 + tmp103 + tmp40 + tmp45 + tmp96 + tmp97 + tmp98 + tmp99;
                                        EM_S[INDEX4(k,m,0,4,numEq,numComp,8)]+=C_2_0*w43 + C_2_3*w35 + C_2_4*w46 + C_2_7*w44 + tmp114 + tmp115 + tmp116 + tmp117 + tmp118 + tmp119 + tmp120 + tmp121 + tmp92 + tmp93;
                                        EM_S[INDEX4(k,m,0,5,numEq,numComp,8)]+=tmp192 + tmp193 + tmp194 + tmp196 + tmp197 + tmp198 + tmp224 + tmp225 + tmp226 + tmp227;
                                        EM_S[INDEX4(k,m,0,6,numEq,numComp,8)]+=tmp232 + tmp234 + tmp238 + tmp239 + tmp240 + tmp241 + tmp260 + tmp261 + tmp262 + tmp263;
                                        EM_S[INDEX4(k,m,0,7,numEq,numComp,8)]+=tmp60 + tmp61 + tmp62 + tmp63 + tmp64 + tmp65;
                                        EM_S[INDEX4(k,m,1,0,numEq,numComp,8)]+=-C_0_0*w41 - C_0_1*w50 - C_0_6*w49 - C_0_7*w33 + tmp148 + tmp149 + tmp242 + tmp243 + tmp244 + tmp245 + tmp246 + tmp247 + tmp248 + tmp249;
                                        EM_S[INDEX4(k,m,1,1,numEq,numComp,8)]+=C_0_0*w41 + C_0_1*w50 + C_0_6*w49 + C_0_7*w33 + C_1_1*w47 - C_1_3*w42 - C_1_4*w34 + C_1_6*w48 - C_2_1*w43 - C_2_2*w35 - C_2_5*w46 - C_2_6*w44 + tmp151 + tmp155 + tmp164 + tmp167 + tmp168 + tmp169;
                                        EM_S[INDEX4(k,m,1,2,numEq,numComp,8)]+=tmp101 + tmp103 + tmp40 + tmp42 + tmp44 + tmp45 + tmp46 + tmp48 + tmp96 + tmp98;
                                        EM_S[INDEX4(k,m,1,3,numEq,numComp,8)]+=-C_1_1*w47 + C_1_3*w42 + C_1_4*w34 - C_1_6*w48 + tmp20 + tmp21 + tmp22 + tmp23 + tmp24 + tmp25 + tmp26 + tmp27 + tmp28 + tmp29;
                                        EM_S[INDEX4(k,m,1,4,numEq,numComp,8)]+=tmp193 + tmp197 + tmp214 + tmp215 + tmp216 + tmp217 + tmp224 + tmp225 + tmp226 + tmp227;
                                        EM_S[INDEX4(k,m,1,5,numEq,numComp,8)]+=C_2_1*w43 + C_2_2*w35 + C_2_5*w46 + C_2_6*w44 + tmp70 + tmp75 + tmp83 + tmp84 + tmp85 + tmp86 + tmp88 + tmp89 + tmp90 + tmp91;
                                        EM_S[INDEX4(k,m,1,6,numEq,numComp,8)]+=tmp60 + tmp61 + tmp63 + tmp65 + tmp80 + tmp81;
                                        EM_S[INDEX4(k,m,1,7,numEq,numComp,8)]+=tmp181 + tmp183 + tmp184 + tmp186 + tmp187 + tmp189 + tmp254 + tmp255 + tmp256 + tmp257;
                                        EM_S[INDEX4(k,m,2,0,numEq,numComp,8)]+=-C_1_0*w42 + C_1_2*w47 + C_1_5*w48 - C_1_7*w34 + tmp138 + tmp139 + tmp140 + tmp141 + tmp142 + tmp143 + tmp144 + tmp145 + tmp146 + tmp147;
                                        EM_S[INDEX4(k,m,2,1,numEq,numComp,8)]+=tmp100 + tmp102 + tmp40 + tmp41 + tmp43 + tmp45 + tmp47 + tmp49 + tmp97 + tmp99;
                                        EM_S[INDEX4(k,m,2,2,numEq,numComp,8)]+=-C_0_2*w50 - C_0_3*w41 - C_0_4*w33 - C_0_5*w49 + C_1_0*w42 - C_1_2*w47 - C_1_5*w48 + C_1_7*w34 - C_2_1*w35 - C_2_2*w43 - C_2_5*w44 - C_2_6*w46 + tmp152 + tmp154 + tmp164 + tmp165 + tmp166 + tmp167;
                                        EM_S[INDEX4(k,m,2,3,numEq,numComp,8)]+=C_0_2*w50 + C_0_3*w41 + C_0_4*w33 + C_0_5*w49 + tmp200 + tmp201 + tmp202 + tmp203 + tmp204 + tmp205 + tmp206 + tmp207 + tmp222 + tmp223;
                                        EM_S[INDEX4(k,m,2,4,numEq,numComp,8)]+=tmp229 + tmp231 + tmp232 + tmp234 + tmp235 + tmp237 + tmp260 + tmp261 + tmp262 + tmp263;
                                        EM_S[INDEX4(k,m,2,5,numEq,numComp,8)]+=tmp60 + tmp62 + tmp63 + tmp64 + tmp94 + tmp95;
                                        EM_S[INDEX4(k,m,2,6,numEq,numComp,8)]+=C_2_1*w35 + C_2_2*w43 + C_2_5*w44 + C_2_6*w46 + tmp70 + tmp71 + tmp72 + tmp73 + tmp74 + tmp75 + tmp76 + tmp77 + tmp78 + tmp79;
                                        EM_S[INDEX4(k,m,2,7,numEq,numComp,8)]+=tmp107 + tmp111 + tmp156 + tmp157 + tmp158 + tmp159 + tmp160 + tmp161 + tmp162 + tmp163;
                                        EM_S[INDEX4(k,m,3,0,numEq,numComp,8)]+=tmp40 + tmp41 + tmp42 + tmp43 + tmp44 + tmp45 + tmp46 + tmp47 + tmp48 + tmp49;
                                        EM_S[INDEX4(k,m,3,1,numEq,numComp,8)]+=-C_1_1*w42 + C_1_3*w47 + C_1_4*w48 - C_1_6*w34 + tmp20 + tmp21 + tmp22 + tmp24 + tmp25 + tmp26 + tmp28 + tmp29 + tmp33 + tmp37;
                                        EM_S[INDEX4(k,m,3,2,numEq,numComp,8)]+=-C_0_2*w41 - C_0_3*w50 - C_0_4*w49 - C_0_5*w33 + tmp200 + tmp201 + tmp202 + tmp203 + tmp204 + tmp205 + tmp206 + tmp207 + tmp54 + tmp56;
                                        EM_S[INDEX4(k,m,3,3,numEq,numComp,8)]+=C_0_2*w41 + C_0_3*w50 + C_0_4*w49 + C_0_5*w33 + C_1_1*w42 - C_1_3*w47 - C_1_4*w48 + C_1_6*w34 - C_2_0*w35 - C_2_3*w43 - C_2_4*w44 - C_2_7*w46 + tmp132 + tmp133 + tmp134 + tmp135 + tmp136 + tmp137;
                                        EM_S[INDEX4(k,m,3,4,numEq,numComp,8)]+=tmp60 + tmp63 + tmp80 + tmp81 + tmp94 + tmp95;
                                        EM_S[INDEX4(k,m,3,5,numEq,numComp,8)]+=tmp184 + tmp186 + tmp218 + tmp219 + tmp220 + tmp221 + tmp254 + tmp255 + tmp256 + tmp257;
                                        EM_S[INDEX4(k,m,3,6,numEq,numComp,8)]+=tmp106 + tmp107 + tmp108 + tmp110 + tmp111 + tmp112 + tmp160 + tmp161 + tmp162 + tmp163;
                                        EM_S[INDEX4(k,m,3,7,numEq,numComp,8)]+=C_2_0*w35 + C_2_3*w43 + C_2_4*w44 + C_2_7*w46 + tmp1 + tmp2 + tmp3 + tmp4 + tmp6 + tmp7 + tmp8 + tmp9 + tmp92 + tmp93;
                                        EM_S[INDEX4(k,m,4,0,numEq,numComp,8)]+=-C_2_0*w46 - C_2_3*w44 - C_2_4*w43 - C_2_7*w35 + tmp0 + tmp114 + tmp115 + tmp116 + tmp117 + tmp118 + tmp119 + tmp120 + tmp121 + tmp5;
                                        EM_S[INDEX4(k,m,4,1,numEq,numComp,8)]+=tmp190 + tmp191 + tmp192 + tmp193 + tmp194 + tmp195 + tmp196 + tmp197 + tmp198 + tmp199;
                                        EM_S[INDEX4(k,m,4,2,numEq,numComp,8)]+=tmp228 + tmp230 + tmp232 + tmp233 + tmp234 + tmp236 + tmp238 + tmp239 + tmp240 + tmp241;
                                        EM_S[INDEX4(k,m,4,3,numEq,numComp,8)]+=tmp258 + tmp259 + tmp61 + tmp62 + tmp64 + tmp65;
                                        EM_S[INDEX4(k,m,4,4,numEq,numComp,8)]+=-C_0_2*w33 - C_0_3*w49 - C_0_4*w50 - C_0_5*w41 - C_1_1*w34 + C_1_3*w48 + C_1_4*w47 - C_1_6*w42 + C_2_0*w46 + C_2_3*w44 + C_2_4*w43 + C_2_7*w35 + tmp150 + tmp151 + tmp152 + tmp153 + tmp154 + tmp155;
                                        EM_S[INDEX4(k,m,4,5,numEq,numComp,8)]+=C_0_2*w33 + C_0_3*w49 + C_0_4*w50 + C_0_5*w41 + tmp222 + tmp223 + tmp50 + tmp51 + tmp52 + tmp53 + tmp55 + tmp57 + tmp58 + tmp59;
                                        EM_S[INDEX4(k,m,4,6,numEq,numComp,8)]+=C_1_1*w34 - C_1_3*w48 - C_1_4*w47 + C_1_6*w42 + tmp23 + tmp27 + tmp30 + tmp31 + tmp32 + tmp34 + tmp35 + tmp36 + tmp38 + tmp39;
                                        EM_S[INDEX4(k,m,4,7,numEq,numComp,8)]+=tmp10 + tmp11 + tmp12 + tmp13 + tmp14 + tmp15 + tmp16 + tmp17 + tmp18 + tmp19;
                                        EM_S[INDEX4(k,m,5,0,numEq,numComp,8)]+=tmp190 + tmp191 + tmp193 + tmp195 + tmp197 + tmp199 + tmp214 + tmp215 + tmp216 + tmp217;
                                        EM_S[INDEX4(k,m,5,1,numEq,numComp,8)]+=-C_2_1*w46 - C_2_2*w44 - C_2_5*w43 - C_2_6*w35 + tmp82 + tmp83 + tmp84 + tmp85 + tmp86 + tmp87 + tmp88 + tmp89 + tmp90 + tmp91;
                                        EM_S[INDEX4(k,m,5,2,numEq,numComp,8)]+=tmp258 + tmp259 + tmp61 + tmp65 + tmp80 + tmp81;
                                        EM_S[INDEX4(k,m,5,3,numEq,numComp,8)]+=tmp180 + tmp181 + tmp182 + tmp183 + tmp184 + tmp185 + tmp186 + tmp187 + tmp188 + tmp189;
                                        EM_S[INDEX4(k,m,5,4,numEq,numComp,8)]+=-C_0_2*w49 - C_0_3*w33 - C_0_4*w41 - C_0_5*w50 + tmp50 + tmp51 + tmp52 + tmp53 + tmp54 + tmp55 + tmp56 + tmp57 + tmp58 + tmp59;
                                        EM_S[INDEX4(k,m,5,5,numEq,numComp,8)]+=C_0_2*w49 + C_0_3*w33 + C_0_4*w41 + C_0_5*w50 - C_1_0*w34 + C_1_2*w48 + C_1_5*w47 - C_1_7*w42 + C_2_1*w46 + C_2_2*w44 + C_2_5*w43 + C_2_6*w35 + tmp134 + tmp135 + tmp208 + tmp211 + tmp212 + tmp213;
                                        EM_S[INDEX4(k,m,5,6,numEq,numComp,8)]+=tmp10 + tmp11 + tmp13 + tmp15 + tmp17 + tmp19 + tmp66 + tmp67 + tmp68 + tmp69;
                                        EM_S[INDEX4(k,m,5,7,numEq,numComp,8)]+=C_1_0*w34 - C_1_2*w48 - C_1_5*w47 + C_1_7*w42 + tmp170 + tmp171 + tmp172 + tmp173 + tmp174 + tmp175 + tmp176 + tmp177 + tmp178 + tmp179;
                                        EM_S[INDEX4(k,m,6,0,numEq,numComp,8)]+=tmp228 + tmp229 + tmp230 + tmp231 + tmp232 + tmp233 + tmp234 + tmp235 + tmp236 + tmp237;
                                        EM_S[INDEX4(k,m,6,1,numEq,numComp,8)]+=tmp258 + tmp259 + tmp62 + tmp64 + tmp94 + tmp95;
                                        EM_S[INDEX4(k,m,6,2,numEq,numComp,8)]+=-C_2_1*w44 - C_2_2*w46 - C_2_5*w35 - C_2_6*w43 + tmp71 + tmp72 + tmp73 + tmp74 + tmp76 + tmp77 + tmp78 + tmp79 + tmp82 + tmp87;
                                        EM_S[INDEX4(k,m,6,3,numEq,numComp,8)]+=tmp104 + tmp105 + tmp107 + tmp109 + tmp111 + tmp113 + tmp156 + tmp157 + tmp158 + tmp159;
                                        EM_S[INDEX4(k,m,6,4,numEq,numComp,8)]+=C_1_1*w48 - C_1_3*w34 - C_1_4*w42 + C_1_6*w47 + tmp30 + tmp31 + tmp32 + tmp33 + tmp34 + tmp35 + tmp36 + tmp37 + tmp38 + tmp39;
                                        EM_S[INDEX4(k,m,6,5,numEq,numComp,8)]+=tmp10 + tmp12 + tmp14 + tmp15 + tmp16 + tmp18 + tmp250 + tmp251 + tmp252 + tmp253;
                                        EM_S[INDEX4(k,m,6,6,numEq,numComp,8)]+=-C_0_0*w33 - C_0_1*w49 - C_0_6*w50 - C_0_7*w41 - C_1_1*w48 + C_1_3*w34 + C_1_4*w42 - C_1_6*w47 + C_2_1*w44 + C_2_2*w46 + C_2_5*w35 + C_2_6*w43 + tmp133 + tmp136 + tmp209 + tmp210 + tmp212 + tmp213;
                                        EM_S[INDEX4(k,m,6,7,numEq,numComp,8)]+=C_0_0*w33 + C_0_1*w49 + C_0_6*w50 + C_0_7*w41 + tmp122 + tmp123 + tmp124 + tmp125 + tmp126 + tmp127 + tmp128 + tmp129 + tmp130 + tmp131;
                                        EM_S[INDEX4(k,m,7,0,numEq,numComp,8)]+=tmp258 + tmp259 + tmp80 + tmp81 + tmp94 + tmp95;
                                        EM_S[INDEX4(k,m,7,1,numEq,numComp,8)]+=tmp180 + tmp182 + tmp184 + tmp185 + tmp186 + tmp188 + tmp218 + tmp219 + tmp220 + tmp221;
                                        EM_S[INDEX4(k,m,7,2,numEq,numComp,8)]+=tmp104 + tmp105 + tmp106 + tmp107 + tmp108 + tmp109 + tmp110 + tmp111 + tmp112 + tmp113;
                                        EM_S[INDEX4(k,m,7,3,numEq,numComp,8)]+=-C_2_0*w44 - C_2_3*w46 - C_2_4*w35 - C_2_7*w43 + tmp0 + tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6 + tmp7 + tmp8 + tmp9;
                                        EM_S[INDEX4(k,m,7,4,numEq,numComp,8)]+=tmp10 + tmp15 + tmp250 + tmp251 + tmp252 + tmp253 + tmp66 + tmp67 + tmp68 + tmp69;
                                        EM_S[INDEX4(k,m,7,5,numEq,numComp,8)]+=C_1_0*w48 - C_1_2*w34 - C_1_5*w42 + C_1_7*w47 + tmp141 + tmp145 + tmp170 + tmp171 + tmp172 + tmp174 + tmp175 + tmp176 + tmp178 + tmp179;
                                        EM_S[INDEX4(k,m,7,6,numEq,numComp,8)]+=-C_0_0*w49 - C_0_1*w33 - C_0_6*w41 - C_0_7*w50 + tmp122 + tmp123 + tmp124 + tmp125 + tmp127 + tmp129 + tmp130 + tmp131 + tmp148 + tmp149;
                                        EM_S[INDEX4(k,m,7,7,numEq,numComp,8)]+=C_0_0*w49 + C_0_1*w33 + C_0_6*w41 + C_0_7*w50 - C_1_0*w48 + C_1_2*w34 + C_1_5*w42 - C_1_7*w47 + C_2_0*w44 + C_2_3*w46 + C_2_4*w35 + C_2_7*w43 + tmp150 + tmp153 + tmp165 + tmp166 + tmp168 + tmp169;
                                    }
                                }
                            } else { // constant data
                                for (index_t k=0; k<numEq; k++) {
                                    for (index_t m=0; m<numComp; m++) {
                                        const double wC0 = C_p[INDEX3(k,m,0,numEq,numComp)]*w55;
                                        const double wC1 = C_p[INDEX3(k,m,1,numEq,numComp)]*w56;
                                        const double wC2 = C_p[INDEX3(k,m,2,numEq,numComp)]*w54;
                                        EM_S[INDEX4(k,m,0,0,numEq,numComp,8)]+= 4*wC0 + 4*wC1 + 4*wC2;
                                        EM_S[INDEX4(k,m,0,1,numEq,numComp,8)]+=-4*wC0 + 2*wC1 + 2*wC2;
                                        EM_S[INDEX4(k,m,0,2,numEq,numComp,8)]+= 2*wC0 - 4*wC1 + 2*wC2;
                                        EM_S[INDEX4(k,m,0,3,numEq,numComp,8)]+=-2*wC0 - 2*wC1 +   wC2;
                                        EM_S[INDEX4(k,m,0,4,numEq,numComp,8)]+= 2*wC0 + 2*wC1 - 4*wC2;
                                        EM_S[INDEX4(k,m,0,5,numEq,numComp,8)]+=-2*wC0 +   wC1 - 2*wC2;
                                        EM_S[INDEX4(k,m,0,6,numEq,numComp,8)]+=   wC0 - 2*wC1 - 2*wC2;
                                        EM_S[INDEX4(k,m,0,7,numEq,numComp,8)]+=  -wC0 -   wC1 -   wC2;
                                        EM_S[INDEX4(k,m,1,0,numEq,numComp,8)]+= 4*wC0 + 2*wC1 + 2*wC2;
                                        EM_S[INDEX4(k,m,1,1,numEq,numComp,8)]+=-4*wC0 + 4*wC1 + 4*wC2;
                                        EM_S[INDEX4(k,m,1,2,numEq,numComp,8)]+= 2*wC0 - 2*wC1 +   wC2;
                                        EM_S[INDEX4(k,m,1,3,numEq,numComp,8)]+=-2*wC0 - 4*wC1 + 2*wC2;
                                        EM_S[INDEX4(k,m,1,4,numEq,numComp,8)]+= 2*wC0 +   wC1 - 2*wC2;
                                        EM_S[INDEX4(k,m,1,5,numEq,numComp,8)]+=-2*wC0 + 2*wC1 - 4*wC2;
                                        EM_S[INDEX4(k,m,1,6,numEq,numComp,8)]+=   wC0 -   wC1 -   wC2;
                                        EM_S[INDEX4(k,m,1,7,numEq,numComp,8)]+=  -wC0 - 2*wC1 - 2*wC2;
                                        EM_S[INDEX4(k,m,2,0,numEq,numComp,8)]+= 2*wC0 + 4*wC1 + 2*wC2;
                                        EM_S[INDEX4(k,m,2,1,numEq,numComp,8)]+=-2*wC0 + 2*wC1 +   wC2;
                                        EM_S[INDEX4(k,m,2,2,numEq,numComp,8)]+= 4*wC0 - 4*wC1 + 4*wC2;
                                        EM_S[INDEX4(k,m,2,3,numEq,numComp,8)]+=-4*wC0 - 2*wC1 + 2*wC2;
                                        EM_S[INDEX4(k,m,2,4,numEq,numComp,8)]+=   wC0 + 2*wC1 - 2*wC2;
                                        EM_S[INDEX4(k,m,2,5,numEq,numComp,8)]+=  -wC0 +   wC1 -   wC2;
                                        EM_S[INDEX4(k,m,2,6,numEq,numComp,8)]+= 2*wC0 - 2*wC1 - 4*wC2;
                                        EM_S[INDEX4(k,m,2,7,numEq,numComp,8)]+=-2*wC0 -   wC1 - 2*wC2;
                                        EM_S[INDEX4(k,m,3,0,numEq,numComp,8)]+= 2*wC0 + 2*wC1 +   wC2;
                                        EM_S[INDEX4(k,m,3,1,numEq,numComp,8)]+=-2*wC0 + 4*wC1 + 2*wC2;
                                        EM_S[INDEX4(k,m,3,2,numEq,numComp,8)]+= 4*wC0 - 2*wC1 + 2*wC2;
                                        EM_S[INDEX4(k,m,3,3,numEq,numComp,8)]+=-4*wC0 - 4*wC1 + 4*wC2;
                                        EM_S[INDEX4(k,m,3,4,numEq,numComp,8)]+=   wC0 +   wC1 -   wC2;
                                        EM_S[INDEX4(k,m,3,5,numEq,numComp,8)]+=  -wC0 + 2*wC1 - 2*wC2;
                                        EM_S[INDEX4(k,m,3,6,numEq,numComp,8)]+= 2*wC0 -   wC1 - 2*wC2;
                                        EM_S[INDEX4(k,m,3,7,numEq,numComp,8)]+=-2*wC0 - 2*wC1 - 4*wC2;
                                        EM_S[INDEX4(k,m,4,0,numEq,numComp,8)]+= 2*wC0 + 2*wC1 + 4*wC2;
                                        EM_S[INDEX4(k,m,4,1,numEq,numComp,8)]+=-2*wC0 +   wC1 + 2*wC2;
                                        EM_S[INDEX4(k,m,4,2,numEq,numComp,8)]+=   wC0 - 2*wC1 + 2*wC2;
                                        EM_S[INDEX4(k,m,4,3,numEq,numComp,8)]+=  -wC0 -   wC1 +   wC2;
                                        EM_S[INDEX4(k,m,4,4,numEq,numComp,8)]+= 4*wC0 + 4*wC1 - 4*wC2;
                                        EM_S[INDEX4(k,m,4,5,numEq,numComp,8)]+=-4*wC0 + 2*wC1 - 2*wC2;
                                        EM_S[INDEX4(k,m,4,6,numEq,numComp,8)]+= 2*wC0 - 4*wC1 - 2*wC2;
                                        EM_S[INDEX4(k,m,4,7,numEq,numComp,8)]+=-2*wC0 - 2*wC1 -   wC2;
                                        EM_S[INDEX4(k,m,5,0,numEq,numComp,8)]+= 2*wC0 +   wC1 + 2*wC2;
                                        EM_S[INDEX4(k,m,5,1,numEq,numComp,8)]+=-2*wC0 + 2*wC1 + 4*wC2;
                                        EM_S[INDEX4(k,m,5,2,numEq,numComp,8)]+=   wC0 -   wC1 +   wC2;
                                        EM_S[INDEX4(k,m,5,3,numEq,numComp,8)]+=  -wC0 - 2*wC1 + 2*wC2;
                                        EM_S[INDEX4(k,m,5,4,numEq,numComp,8)]+= 4*wC0 + 2*wC1 - 2*wC2;
                                        EM_S[INDEX4(k,m,5,5,numEq,numComp,8)]+=-4*wC0 + 4*wC1 - 4*wC2;
                                        EM_S[INDEX4(k,m,5,6,numEq,numComp,8)]+= 2*wC0 - 2*wC1 -   wC2;
                                        EM_S[INDEX4(k,m,5,7,numEq,numComp,8)]+=-2*wC0 - 4*wC1 - 2*wC2;
                                        EM_S[INDEX4(k,m,6,0,numEq,numComp,8)]+=   wC0 + 2*wC1 + 2*wC2;
                                        EM_S[INDEX4(k,m,6,1,numEq,numComp,8)]+=  -wC0 +   wC1 +   wC2;
                                        EM_S[INDEX4(k,m,6,2,numEq,numComp,8)]+= 2*wC0 - 2*wC1 + 4*wC2;
                                        EM_S[INDEX4(k,m,6,3,numEq,numComp,8)]+=-2*wC0 -   wC1 + 2*wC2;
                                        EM_S[INDEX4(k,m,6,4,numEq,numComp,8)]+= 2*wC0 + 4*wC1 - 2*wC2;
                                        EM_S[INDEX4(k,m,6,5,numEq,numComp,8)]+=-2*wC0 + 2*wC1 -   wC2;
                                        EM_S[INDEX4(k,m,6,6,numEq,numComp,8)]+= 4*wC0 - 4*wC1 - 4*wC2;
                                        EM_S[INDEX4(k,m,6,7,numEq,numComp,8)]+=-4*wC0 - 2*wC1 - 2*wC2;
                                        EM_S[INDEX4(k,m,7,0,numEq,numComp,8)]+=   wC0 +   wC1 +   wC2;
                                        EM_S[INDEX4(k,m,7,1,numEq,numComp,8)]+=  -wC0 + 2*wC1 + 2*wC2;
                                        EM_S[INDEX4(k,m,7,2,numEq,numComp,8)]+= 2*wC0 -   wC1 + 2*wC2;
                                        EM_S[INDEX4(k,m,7,3,numEq,numComp,8)]+=-2*wC0 - 2*wC1 + 4*wC2;
                                        EM_S[INDEX4(k,m,7,4,numEq,numComp,8)]+= 2*wC0 + 2*wC1 -   wC2;
                                        EM_S[INDEX4(k,m,7,5,numEq,numComp,8)]+=-2*wC0 + 4*wC1 - 2*wC2;
                                        EM_S[INDEX4(k,m,7,6,numEq,numComp,8)]+= 4*wC0 - 2*wC1 - 2*wC2;
                                        EM_S[INDEX4(k,m,7,7,numEq,numComp,8)]+=-4*wC0 - 4*wC1 - 4*wC2;
                                    }
                                }
                            }
                        }
                        ///////////////
                        // process D //
                        ///////////////
                        if (!D.isEmpty()) {
                            const double* D_p=D.getSampleDataRO(e);
                            if (D.actsExpanded()) {
                                for (index_t k=0; k<numEq; k++) {
                                    for (index_t m=0; m<numComp; m++) {
                                        const double D_0 = D_p[INDEX3(k,m,0,numEq,numComp)];
                                        const double D_1 = D_p[INDEX3(k,m,1,numEq,numComp)];
                                        const double D_2 = D_p[INDEX3(k,m,2,numEq,numComp)];
                                        const double D_3 = D_p[INDEX3(k,m,3,numEq,numComp)];
                                        const double D_4 = D_p[INDEX3(k,m,4,numEq,numComp)];
                                        const double D_5 = D_p[INDEX3(k,m,5,numEq,numComp)];
                                        const double D_6 = D_p[INDEX3(k,m,6,numEq,numComp)];
                                        const double D_7 = D_p[INDEX3(k,m,7,numEq,numComp)];
                                        const double tmp0 = w59*(D_3 + D_7);
                                        const double tmp1 = w57*(D_0 + D_4);
                                        const double tmp2 = w58*(D_1 + D_2 + D_5 + D_6);
                                        const double tmp3 = w60*(D_0 + D_1 + D_2 + D_3);
                                        const double tmp4 = w61*(D_4 + D_5 + D_6 + D_7);
                                        const double tmp5 = w59*(D_1 + D_3);
                                        const double tmp6 = w57*(D_4 + D_6);
                                        const double tmp7 = w58*(D_0 + D_2 + D_5 + D_7);
                                        const double tmp8 = w59*(D_4 + D_6);
                                        const double tmp9 = w57*(D_1 + D_3);
                                        const double tmp10 = w60*(D_4 + D_5 + D_6 + D_7);
                                        const double tmp11 = w61*(D_0 + D_1 + D_2 + D_3);
                                        const double tmp12 = w59*(D_4 + D_5);
                                        const double tmp13 = w57*(D_2 + D_3);
                                        const double tmp14 = w58*(D_0 + D_1 + D_6 + D_7);
                                        const double tmp15 = w58*(D_0 + D_1 + D_2 + D_3 + D_4 + D_5 + D_6 + D_7);
                                        const double tmp16 = w59*(D_2 + D_6);
                                        const double tmp17 = w57*(D_1 + D_5);
                                        const double tmp18 = w58*(D_0 + D_3 + D_4 + D_7);
                                        const double tmp19 = w59*(D_1 + D_5);
                                        const double tmp20 = w57*(D_2 + D_6);
                                        const double tmp21 = w60*(D_0 + D_1 + D_4 + D_5);
                                        const double tmp22 = w61*(D_2 + D_3 + D_6 + D_7);
                                        const double tmp23 = w59*(D_0 + D_4);
                                        const double tmp24 = w57*(D_3 + D_7);
                                        const double tmp25 = w59*(D_6 + D_7);
                                        const double tmp26 = w57*(D_0 + D_1);
                                        const double tmp27 = w58*(D_2 + D_3 + D_4 + D_5);
                                        const double tmp28 = w60*(D_0 + D_5 + D_6);
                                        const double tmp29 = w61*(D_1 + D_2 + D_7);
                                        const double tmp30 = w59*(D_0 + D_2);
                                        const double tmp31 = w57*(D_5 + D_7);
                                        const double tmp32 = w58*(D_1 + D_3 + D_4 + D_6);
                                        const double tmp33 = w60*(D_1 + D_2 + D_7);
                                        const double tmp34 = w61*(D_0 + D_5 + D_6);
                                        const double tmp35 = w60*(D_1 + D_4 + D_7);
                                        const double tmp36 = w61*(D_0 + D_3 + D_6);
                                        const double tmp37 = w60*(D_1 + D_2 + D_4);
                                        const double tmp38 = w61*(D_3 + D_5 + D_6);
                                        const double tmp39 = w59*(D_5 + D_7);
                                        const double tmp40 = w57*(D_0 + D_2);
                                        const double tmp41 = w60*(D_0 + D_2 + D_4 + D_6);
                                        const double tmp42 = w61*(D_1 + D_3 + D_5 + D_7);
                                        const double tmp43 = w60*(D_2 + D_3 + D_6 + D_7);
                                        const double tmp44 = w61*(D_0 + D_1 + D_4 + D_5);
                                        const double tmp45 = w60*(D_2 + D_4 + D_7);
                                        const double tmp46 = w61*(D_0 + D_3 + D_5);
                                        const double tmp47 = w59*(D_2 + D_3);
                                        const double tmp48 = w57*(D_4 + D_5);
                                        const double tmp49 = w60*(D_3 + D_5 + D_6);
                                        const double tmp50 = w61*(D_1 + D_2 + D_4);
                                        const double tmp51 = w60*(D_0 + D_3 + D_5);
                                        const double tmp52 = w61*(D_2 + D_4 + D_7);
                                        const double tmp53 = w60*(D_0 + D_3 + D_6);
                                        const double tmp54 = w61*(D_1 + D_4 + D_7);
                                        const double tmp55 = w60*(D_1 + D_3 + D_5 + D_7);
                                        const double tmp56 = w61*(D_0 + D_2 + D_4 + D_6);
                                        const double tmp57 = w59*(D_0 + D_1);
                                        const double tmp58 = w57*(D_6 + D_7);
                                        EM_S[INDEX4(k,m,0,0,numEq,numComp,8)]+=D_0*w62 + D_7*w63 + tmp49 + tmp50;
                                        EM_S[INDEX4(k,m,0,1,numEq,numComp,8)]+=tmp27 + tmp57 + tmp58;
                                        EM_S[INDEX4(k,m,0,2,numEq,numComp,8)]+=tmp30 + tmp31 + tmp32;
                                        EM_S[INDEX4(k,m,0,3,numEq,numComp,8)]+=tmp10 + tmp11;
                                        EM_S[INDEX4(k,m,0,4,numEq,numComp,8)]+=tmp2 + tmp23 + tmp24;
                                        EM_S[INDEX4(k,m,0,5,numEq,numComp,8)]+=tmp43 + tmp44;
                                        EM_S[INDEX4(k,m,0,6,numEq,numComp,8)]+=tmp55 + tmp56;
                                        EM_S[INDEX4(k,m,0,7,numEq,numComp,8)]+=tmp15;
                                        EM_S[INDEX4(k,m,1,0,numEq,numComp,8)]+=tmp27 + tmp57 + tmp58;
                                        EM_S[INDEX4(k,m,1,1,numEq,numComp,8)]+=D_1*w62 + D_6*w63 + tmp45 + tmp46;
                                        EM_S[INDEX4(k,m,1,2,numEq,numComp,8)]+=tmp10 + tmp11;
                                        EM_S[INDEX4(k,m,1,3,numEq,numComp,8)]+=tmp5 + tmp6 + tmp7;
                                        EM_S[INDEX4(k,m,1,4,numEq,numComp,8)]+=tmp43 + tmp44;
                                        EM_S[INDEX4(k,m,1,5,numEq,numComp,8)]+=tmp18 + tmp19 + tmp20;
                                        EM_S[INDEX4(k,m,1,6,numEq,numComp,8)]+=tmp15;
                                        EM_S[INDEX4(k,m,1,7,numEq,numComp,8)]+=tmp41 + tmp42;
                                        EM_S[INDEX4(k,m,2,0,numEq,numComp,8)]+=tmp30 + tmp31 + tmp32;
                                        EM_S[INDEX4(k,m,2,1,numEq,numComp,8)]+=tmp10 + tmp11;
                                        EM_S[INDEX4(k,m,2,2,numEq,numComp,8)]+=D_2*w62 + D_5*w63 + tmp35 + tmp36;
                                        EM_S[INDEX4(k,m,2,3,numEq,numComp,8)]+=tmp14 + tmp47 + tmp48;
                                        EM_S[INDEX4(k,m,2,4,numEq,numComp,8)]+=tmp55 + tmp56;
                                        EM_S[INDEX4(k,m,2,5,numEq,numComp,8)]+=tmp15;
                                        EM_S[INDEX4(k,m,2,6,numEq,numComp,8)]+=tmp16 + tmp17 + tmp18;
                                        EM_S[INDEX4(k,m,2,7,numEq,numComp,8)]+=tmp21 + tmp22;
                                        EM_S[INDEX4(k,m,3,0,numEq,numComp,8)]+=tmp10 + tmp11;
                                        EM_S[INDEX4(k,m,3,1,numEq,numComp,8)]+=tmp5 + tmp6 + tmp7;
                                        EM_S[INDEX4(k,m,3,2,numEq,numComp,8)]+=tmp14 + tmp47 + tmp48;
                                        EM_S[INDEX4(k,m,3,3,numEq,numComp,8)]+=D_3*w62 + D_4*w63 + tmp28 + tmp29;
                                        EM_S[INDEX4(k,m,3,4,numEq,numComp,8)]+=tmp15;
                                        EM_S[INDEX4(k,m,3,5,numEq,numComp,8)]+=tmp41 + tmp42;
                                        EM_S[INDEX4(k,m,3,6,numEq,numComp,8)]+=tmp21 + tmp22;
                                        EM_S[INDEX4(k,m,3,7,numEq,numComp,8)]+=tmp0 + tmp1 + tmp2;
                                        EM_S[INDEX4(k,m,4,0,numEq,numComp,8)]+=tmp2 + tmp23 + tmp24;
                                        EM_S[INDEX4(k,m,4,1,numEq,numComp,8)]+=tmp43 + tmp44;
                                        EM_S[INDEX4(k,m,4,2,numEq,numComp,8)]+=tmp55 + tmp56;
                                        EM_S[INDEX4(k,m,4,3,numEq,numComp,8)]+=tmp15;
                                        EM_S[INDEX4(k,m,4,4,numEq,numComp,8)]+=D_3*w63 + D_4*w62 + tmp33 + tmp34;
                                        EM_S[INDEX4(k,m,4,5,numEq,numComp,8)]+=tmp12 + tmp13 + tmp14;
                                        EM_S[INDEX4(k,m,4,6,numEq,numComp,8)]+=tmp7 + tmp8 + tmp9;
                                        EM_S[INDEX4(k,m,4,7,numEq,numComp,8)]+=tmp3 + tmp4;
                                        EM_S[INDEX4(k,m,5,0,numEq,numComp,8)]+=tmp43 + tmp44;
                                        EM_S[INDEX4(k,m,5,1,numEq,numComp,8)]+=tmp18 + tmp19 + tmp20;
                                        EM_S[INDEX4(k,m,5,2,numEq,numComp,8)]+=tmp15;
                                        EM_S[INDEX4(k,m,5,3,numEq,numComp,8)]+=tmp41 + tmp42;
                                        EM_S[INDEX4(k,m,5,4,numEq,numComp,8)]+=tmp12 + tmp13 + tmp14;
                                        EM_S[INDEX4(k,m,5,5,numEq,numComp,8)]+=D_2*w63 + D_5*w62 + tmp53 + tmp54;
                                        EM_S[INDEX4(k,m,5,6,numEq,numComp,8)]+=tmp3 + tmp4;
                                        EM_S[INDEX4(k,m,5,7,numEq,numComp,8)]+=tmp32 + tmp39 + tmp40;
                                        EM_S[INDEX4(k,m,6,0,numEq,numComp,8)]+=tmp55 + tmp56;
                                        EM_S[INDEX4(k,m,6,1,numEq,numComp,8)]+=tmp15;
                                        EM_S[INDEX4(k,m,6,2,numEq,numComp,8)]+=tmp16 + tmp17 + tmp18;
                                        EM_S[INDEX4(k,m,6,3,numEq,numComp,8)]+=tmp21 + tmp22;
                                        EM_S[INDEX4(k,m,6,4,numEq,numComp,8)]+=tmp7 + tmp8 + tmp9;
                                        EM_S[INDEX4(k,m,6,5,numEq,numComp,8)]+=tmp3 + tmp4;
                                        EM_S[INDEX4(k,m,6,6,numEq,numComp,8)]+=D_1*w63 + D_6*w62 + tmp51 + tmp52;
                                        EM_S[INDEX4(k,m,6,7,numEq,numComp,8)]+=tmp25 + tmp26 + tmp27;
                                        EM_S[INDEX4(k,m,7,0,numEq,numComp,8)]+=tmp15;
                                        EM_S[INDEX4(k,m,7,1,numEq,numComp,8)]+=tmp41 + tmp42;
                                        EM_S[INDEX4(k,m,7,2,numEq,numComp,8)]+=tmp21 + tmp22;
                                        EM_S[INDEX4(k,m,7,3,numEq,numComp,8)]+=tmp0 + tmp1 + tmp2;
                                        EM_S[INDEX4(k,m,7,4,numEq,numComp,8)]+=tmp3 + tmp4;
                                        EM_S[INDEX4(k,m,7,5,numEq,numComp,8)]+=tmp32 + tmp39 + tmp40;
                                        EM_S[INDEX4(k,m,7,6,numEq,numComp,8)]+=tmp25 + tmp26 + tmp27;
                                        EM_S[INDEX4(k,m,7,7,numEq,numComp,8)]+=D_0*w63 + D_7*w62 + tmp37 + tmp38;
                                    }
                                 }
                            } else { // constant data
                                for (index_t k=0; k<numEq; k++) {
                                    for (index_t m=0; m<numComp; m++) {
                                        const double wD0 = 8*D_p[INDEX2(k, m, numEq)]*w58;
                                        EM_S[INDEX4(k,m,0,0,numEq,numComp,8)]+=8*wD0;
                                        EM_S[INDEX4(k,m,0,1,numEq,numComp,8)]+=4*wD0;
                                        EM_S[INDEX4(k,m,0,2,numEq,numComp,8)]+=4*wD0;
                                        EM_S[INDEX4(k,m,0,3,numEq,numComp,8)]+=2*wD0;
                                        EM_S[INDEX4(k,m,0,4,numEq,numComp,8)]+=4*wD0;
                                        EM_S[INDEX4(k,m,0,5,numEq,numComp,8)]+=2*wD0;
                                        EM_S[INDEX4(k,m,0,6,numEq,numComp,8)]+=2*wD0;
                                        EM_S[INDEX4(k,m,0,7,numEq,numComp,8)]+=  wD0;
                                        EM_S[INDEX4(k,m,1,0,numEq,numComp,8)]+=4*wD0;
                                        EM_S[INDEX4(k,m,1,1,numEq,numComp,8)]+=8*wD0;
                                        EM_S[INDEX4(k,m,1,2,numEq,numComp,8)]+=2*wD0;
                                        EM_S[INDEX4(k,m,1,3,numEq,numComp,8)]+=4*wD0;
                                        EM_S[INDEX4(k,m,1,4,numEq,numComp,8)]+=2*wD0;
                                        EM_S[INDEX4(k,m,1,5,numEq,numComp,8)]+=4*wD0;
                                        EM_S[INDEX4(k,m,1,6,numEq,numComp,8)]+=  wD0;
                                        EM_S[INDEX4(k,m,1,7,numEq,numComp,8)]+=2*wD0;
                                        EM_S[INDEX4(k,m,2,0,numEq,numComp,8)]+=4*wD0;
                                        EM_S[INDEX4(k,m,2,1,numEq,numComp,8)]+=2*wD0;
                                        EM_S[INDEX4(k,m,2,2,numEq,numComp,8)]+=8*wD0;
                                        EM_S[INDEX4(k,m,2,3,numEq,numComp,8)]+=4*wD0;
                                        EM_S[INDEX4(k,m,2,4,numEq,numComp,8)]+=2*wD0;
                                        EM_S[INDEX4(k,m,2,5,numEq,numComp,8)]+=  wD0;
                                        EM_S[INDEX4(k,m,2,6,numEq,numComp,8)]+=4*wD0;
                                        EM_S[INDEX4(k,m,2,7,numEq,numComp,8)]+=2*wD0;
                                        EM_S[INDEX4(k,m,3,0,numEq,numComp,8)]+=2*wD0;
                                        EM_S[INDEX4(k,m,3,1,numEq,numComp,8)]+=4*wD0;
                                        EM_S[INDEX4(k,m,3,2,numEq,numComp,8)]+=4*wD0;
                                        EM_S[INDEX4(k,m,3,3,numEq,numComp,8)]+=8*wD0;
                                        EM_S[INDEX4(k,m,3,4,numEq,numComp,8)]+=  wD0;
                                        EM_S[INDEX4(k,m,3,5,numEq,numComp,8)]+=2*wD0;
                                        EM_S[INDEX4(k,m,3,6,numEq,numComp,8)]+=2*wD0;
                                        EM_S[INDEX4(k,m,3,7,numEq,numComp,8)]+=4*wD0;
                                        EM_S[INDEX4(k,m,4,0,numEq,numComp,8)]+=4*wD0;
                                        EM_S[INDEX4(k,m,4,1,numEq,numComp,8)]+=2*wD0;
                                        EM_S[INDEX4(k,m,4,2,numEq,numComp,8)]+=2*wD0;
                                        EM_S[INDEX4(k,m,4,3,numEq,numComp,8)]+=  wD0;
                                        EM_S[INDEX4(k,m,4,4,numEq,numComp,8)]+=8*wD0;
                                        EM_S[INDEX4(k,m,4,5,numEq,numComp,8)]+=4*wD0;
                                        EM_S[INDEX4(k,m,4,6,numEq,numComp,8)]+=4*wD0;
                                        EM_S[INDEX4(k,m,4,7,numEq,numComp,8)]+=2*wD0;
                                        EM_S[INDEX4(k,m,5,0,numEq,numComp,8)]+=2*wD0;
                                        EM_S[INDEX4(k,m,5,1,numEq,numComp,8)]+=4*wD0;
                                        EM_S[INDEX4(k,m,5,2,numEq,numComp,8)]+=  wD0;
                                        EM_S[INDEX4(k,m,5,3,numEq,numComp,8)]+=2*wD0;
                                        EM_S[INDEX4(k,m,5,4,numEq,numComp,8)]+=4*wD0;
                                        EM_S[INDEX4(k,m,5,5,numEq,numComp,8)]+=8*wD0;
                                        EM_S[INDEX4(k,m,5,6,numEq,numComp,8)]+=2*wD0;
                                        EM_S[INDEX4(k,m,5,7,numEq,numComp,8)]+=4*wD0;
                                        EM_S[INDEX4(k,m,6,0,numEq,numComp,8)]+=2*wD0;
                                        EM_S[INDEX4(k,m,6,1,numEq,numComp,8)]+=  wD0;
                                        EM_S[INDEX4(k,m,6,2,numEq,numComp,8)]+=4*wD0;
                                        EM_S[INDEX4(k,m,6,3,numEq,numComp,8)]+=2*wD0;
                                        EM_S[INDEX4(k,m,6,4,numEq,numComp,8)]+=4*wD0;
                                        EM_S[INDEX4(k,m,6,5,numEq,numComp,8)]+=2*wD0;
                                        EM_S[INDEX4(k,m,6,6,numEq,numComp,8)]+=8*wD0;
                                        EM_S[INDEX4(k,m,6,7,numEq,numComp,8)]+=4*wD0;
                                        EM_S[INDEX4(k,m,7,0,numEq,numComp,8)]+=  wD0;
                                        EM_S[INDEX4(k,m,7,1,numEq,numComp,8)]+=2*wD0;
                                        EM_S[INDEX4(k,m,7,2,numEq,numComp,8)]+=2*wD0;
                                        EM_S[INDEX4(k,m,7,3,numEq,numComp,8)]+=4*wD0;
                                        EM_S[INDEX4(k,m,7,4,numEq,numComp,8)]+=2*wD0;
                                        EM_S[INDEX4(k,m,7,5,numEq,numComp,8)]+=4*wD0;
                                        EM_S[INDEX4(k,m,7,6,numEq,numComp,8)]+=4*wD0;
                                        EM_S[INDEX4(k,m,7,7,numEq,numComp,8)]+=8*wD0;
                                    }
                                }
                            }
                        }
                        ///////////////
                        // process X //
                        ///////////////
                        if (!X.isEmpty()) {
                            const double* X_p=X.getSampleDataRO(e);
                            if (X.actsExpanded()) {
                                for (index_t k=0; k<numEq; k++) {
                                    const double X_0_0 = X_p[INDEX3(k,0,0,numEq,3)];
                                    const double X_1_0 = X_p[INDEX3(k,1,0,numEq,3)];
                                    const double X_2_0 = X_p[INDEX3(k,2,0,numEq,3)];
                                    const double X_0_1 = X_p[INDEX3(k,0,1,numEq,3)];
                                    const double X_1_1 = X_p[INDEX3(k,1,1,numEq,3)];
                                    const double X_2_1 = X_p[INDEX3(k,2,1,numEq,3)];
                                    const double X_0_2 = X_p[INDEX3(k,0,2,numEq,3)];
                                    const double X_1_2 = X_p[INDEX3(k,1,2,numEq,3)];
                                    const double X_2_2 = X_p[INDEX3(k,2,2,numEq,3)];
                                    const double X_0_3 = X_p[INDEX3(k,0,3,numEq,3)];
                                    const double X_1_3 = X_p[INDEX3(k,1,3,numEq,3)];
                                    const double X_2_3 = X_p[INDEX3(k,2,3,numEq,3)];
                                    const double X_0_4 = X_p[INDEX3(k,0,4,numEq,3)];
                                    const double X_1_4 = X_p[INDEX3(k,1,4,numEq,3)];
                                    const double X_2_4 = X_p[INDEX3(k,2,4,numEq,3)];
                                    const double X_0_5 = X_p[INDEX3(k,0,5,numEq,3)];
                                    const double X_1_5 = X_p[INDEX3(k,1,5,numEq,3)];
                                    const double X_2_5 = X_p[INDEX3(k,2,5,numEq,3)];
                                    const double X_0_6 = X_p[INDEX3(k,0,6,numEq,3)];
                                    const double X_1_6 = X_p[INDEX3(k,1,6,numEq,3)];
                                    const double X_2_6 = X_p[INDEX3(k,2,6,numEq,3)];
                                    const double X_0_7 = X_p[INDEX3(k,0,7,numEq,3)];
                                    const double X_1_7 = X_p[INDEX3(k,1,7,numEq,3)];
                                    const double X_2_7 = X_p[INDEX3(k,2,7,numEq,3)];
                                    const double tmp0 = w72*(X_0_6 + X_0_7);
                                    const double tmp1 = w66*(X_2_0 + X_2_4);
                                    const double tmp2 = w64*(X_0_0 + X_0_1);
                                    const double tmp3 = w68*(X_2_1 + X_2_2 + X_2_5 + X_2_6);
                                    const double tmp4 = w65*(X_1_0 + X_1_2);
                                    const double tmp5 = w70*(X_2_3 + X_2_7);
                                    const double tmp6 = w67*(X_1_1 + X_1_3 + X_1_4 + X_1_6);
                                    const double tmp7 = w71*(X_1_5 + X_1_7);
                                    const double tmp8 = w69*(X_0_2 + X_0_3 + X_0_4 + X_0_5);
                                    const double tmp9 = w72*(-X_0_6 - X_0_7);
                                    const double tmp10 = w66*(X_2_1 + X_2_5);
                                    const double tmp11 = w64*(-X_0_0 - X_0_1);
                                    const double tmp12 = w68*(X_2_0 + X_2_3 + X_2_4 + X_2_7);
                                    const double tmp13 = w65*(X_1_1 + X_1_3);
                                    const double tmp14 = w70*(X_2_2 + X_2_6);
                                    const double tmp15 = w67*(X_1_0 + X_1_2 + X_1_5 + X_1_7);
                                    const double tmp16 = w71*(X_1_4 + X_1_6);
                                    const double tmp17 = w69*(-X_0_2 - X_0_3 - X_0_4 - X_0_5);
                                    const double tmp18 = w72*(X_0_4 + X_0_5);
                                    const double tmp19 = w66*(X_2_2 + X_2_6);
                                    const double tmp20 = w64*(X_0_2 + X_0_3);
                                    const double tmp21 = w65*(-X_1_0 - X_1_2);
                                    const double tmp22 = w70*(X_2_1 + X_2_5);
                                    const double tmp23 = w67*(-X_1_1 - X_1_3 - X_1_4 - X_1_6);
                                    const double tmp24 = w71*(-X_1_5 - X_1_7);
                                    const double tmp25 = w69*(X_0_0 + X_0_1 + X_0_6 + X_0_7);
                                    const double tmp26 = w72*(-X_0_4 - X_0_5);
                                    const double tmp27 = w66*(X_2_3 + X_2_7);
                                    const double tmp28 = w64*(-X_0_2 - X_0_3);
                                    const double tmp29 = w65*(-X_1_1 - X_1_3);
                                    const double tmp30 = w70*(X_2_0 + X_2_4);
                                    const double tmp31 = w67*(-X_1_0 - X_1_2 - X_1_5 - X_1_7);
                                    const double tmp32 = w71*(-X_1_4 - X_1_6);
                                    const double tmp33 = w69*(-X_0_0 - X_0_1 - X_0_6 - X_0_7);
                                    const double tmp34 = w72*(X_0_2 + X_0_3);
                                    const double tmp35 = w66*(-X_2_0 - X_2_4);
                                    const double tmp36 = w64*(X_0_4 + X_0_5);
                                    const double tmp37 = w68*(-X_2_1 - X_2_2 - X_2_5 - X_2_6);
                                    const double tmp38 = w65*(X_1_4 + X_1_6);
                                    const double tmp39 = w70*(-X_2_3 - X_2_7);
                                    const double tmp40 = w71*(X_1_1 + X_1_3);
                                    const double tmp41 = w72*(-X_0_2 - X_0_3);
                                    const double tmp42 = w66*(-X_2_1 - X_2_5);
                                    const double tmp43 = w64*(-X_0_4 - X_0_5);
                                    const double tmp44 = w68*(-X_2_0 - X_2_3 - X_2_4 - X_2_7);
                                    const double tmp45 = w65*(X_1_5 + X_1_7);
                                    const double tmp46 = w70*(-X_2_2 - X_2_6);
                                    const double tmp47 = w71*(X_1_0 + X_1_2);
                                    const double tmp48 = w72*(X_0_0 + X_0_1);
                                    const double tmp49 = w66*(-X_2_2 - X_2_6);
                                    const double tmp50 = w64*(X_0_6 + X_0_7);
                                    const double tmp51 = w65*(-X_1_4 - X_1_6);
                                    const double tmp52 = w70*(-X_2_1 - X_2_5);
                                    const double tmp53 = w71*(-X_1_1 - X_1_3);
                                    const double tmp54 = w72*(-X_0_0 - X_0_1);
                                    const double tmp55 = w66*(-X_2_3 - X_2_7);
                                    const double tmp56 = w64*(-X_0_6 - X_0_7);
                                    const double tmp57 = w65*(-X_1_5 - X_1_7);
                                    const double tmp58 = w70*(-X_2_0 - X_2_4);
                                    const double tmp59 = w71*(-X_1_0 - X_1_2);
                                    EM_F[INDEX2(k,0,numEq)]+=tmp0 + tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6 + tmp7 + tmp8;
                                    EM_F[INDEX2(k,1,numEq)]+=tmp10 + tmp11 + tmp12 + tmp13 + tmp14 + tmp15 + tmp16 + tmp17 + tmp9;
                                    EM_F[INDEX2(k,2,numEq)]+=tmp12 + tmp18 + tmp19 + tmp20 + tmp21 + tmp22 + tmp23 + tmp24 + tmp25;
                                    EM_F[INDEX2(k,3,numEq)]+=tmp26 + tmp27 + tmp28 + tmp29 + tmp3 + tmp30 + tmp31 + tmp32 + tmp33;
                                    EM_F[INDEX2(k,4,numEq)]+=tmp15 + tmp25 + tmp34 + tmp35 + tmp36 + tmp37 + tmp38 + tmp39 + tmp40;
                                    EM_F[INDEX2(k,5,numEq)]+=tmp33 + tmp41 + tmp42 + tmp43 + tmp44 + tmp45 + tmp46 + tmp47 + tmp6;
                                    EM_F[INDEX2(k,6,numEq)]+=tmp31 + tmp44 + tmp48 + tmp49 + tmp50 + tmp51 + tmp52 + tmp53 + tmp8;
                                    EM_F[INDEX2(k,7,numEq)]+=tmp17 + tmp23 + tmp37 + tmp54 + tmp55 + tmp56 + tmp57 + tmp58 + tmp59;
                                }
                            } else { // constant data
                                for (index_t k=0; k<numEq; k++) {
                                    const double wX0 = 18*X_p[INDEX2(k, 0, numEq)]*w55;
                                    const double wX1 = 18*X_p[INDEX2(k, 1, numEq)]*w56;
                                    const double wX2 = 18*X_p[INDEX2(k, 2, numEq)]*w54;
                                    EM_F[INDEX2(k,0,numEq)]+= wX0 + wX1 + wX2;
                                    EM_F[INDEX2(k,1,numEq)]+=-wX0 + wX1 + wX2;
                                    EM_F[INDEX2(k,2,numEq)]+= wX0 - wX1 + wX2;
                                    EM_F[INDEX2(k,3,numEq)]+=-wX0 - wX1 + wX2;
                                    EM_F[INDEX2(k,4,numEq)]+= wX0 + wX1 - wX2;
                                    EM_F[INDEX2(k,5,numEq)]+=-wX0 + wX1 - wX2;
                                    EM_F[INDEX2(k,6,numEq)]+= wX0 - wX1 - wX2;
                                    EM_F[INDEX2(k,7,numEq)]+=-wX0 - wX1 - wX2;
                                }
                            }
                        }
                        ///////////////
                        // process Y //
                        ///////////////
                        if (!Y.isEmpty()) {
                            const double* Y_p=Y.getSampleDataRO(e);
                            if (Y.actsExpanded()) {
                                for (index_t k=0; k<numEq; k++) {
                                    const double Y_0 = Y_p[INDEX2(k, 0, numEq)];
                                    const double Y_1 = Y_p[INDEX2(k, 1, numEq)];
                                    const double Y_2 = Y_p[INDEX2(k, 2, numEq)];
                                    const double Y_3 = Y_p[INDEX2(k, 3, numEq)];
                                    const double Y_4 = Y_p[INDEX2(k, 4, numEq)];
                                    const double Y_5 = Y_p[INDEX2(k, 5, numEq)];
                                    const double Y_6 = Y_p[INDEX2(k, 6, numEq)];
                                    const double Y_7 = Y_p[INDEX2(k, 7, numEq)];
                                    const double tmp0 = w76*(Y_3 + Y_5 + Y_6);
                                    const double tmp1 = w75*(Y_1 + Y_2 + Y_4);
                                    const double tmp2 = w76*(Y_2 + Y_4 + Y_7);
                                    const double tmp3 = w75*(Y_0 + Y_3 + Y_5);
                                    const double tmp4 = w76*(Y_1 + Y_4 + Y_7);
                                    const double tmp5 = w75*(Y_0 + Y_3 + Y_6);
                                    const double tmp6 = w76*(Y_0 + Y_5 + Y_6);
                                    const double tmp7 = w75*(Y_1 + Y_2 + Y_7);
                                    const double tmp8 = w76*(Y_1 + Y_2 + Y_7);
                                    const double tmp9 = w75*(Y_0 + Y_5 + Y_6);
                                    const double tmp10 = w76*(Y_0 + Y_3 + Y_6);
                                    const double tmp11 = w75*(Y_1 + Y_4 + Y_7);
                                    const double tmp12 = w76*(Y_0 + Y_3 + Y_5);
                                    const double tmp13 = w75*(Y_2 + Y_4 + Y_7);
                                    const double tmp14 = w76*(Y_1 + Y_2 + Y_4);
                                    const double tmp15 = w75*(Y_3 + Y_5 + Y_6);
                                    EM_F[INDEX2(k,0,numEq)]+=Y_0*w74 + Y_7*w77 + tmp0 + tmp1;
                                    EM_F[INDEX2(k,1,numEq)]+=Y_1*w74 + Y_6*w77 + tmp2 + tmp3;
                                    EM_F[INDEX2(k,2,numEq)]+=Y_2*w74 + Y_5*w77 + tmp4 + tmp5;
                                    EM_F[INDEX2(k,3,numEq)]+=Y_3*w74 + Y_4*w77 + tmp6 + tmp7;
                                    EM_F[INDEX2(k,4,numEq)]+=Y_3*w77 + Y_4*w74 + tmp8 + tmp9;
                                    EM_F[INDEX2(k,5,numEq)]+=Y_2*w77 + Y_5*w74 + tmp10 + tmp11;
                                    EM_F[INDEX2(k,6,numEq)]+=Y_1*w77 + Y_6*w74 + tmp12 + tmp13;
                                    EM_F[INDEX2(k,7,numEq)]+=Y_0*w77 + Y_7*w74 + tmp14 + tmp15;
                                }
                            } else { // constant data
                                for (index_t k=0; k<numEq; k++) {
                                    EM_F[INDEX2(k,0,numEq)]+=216*Y_p[k]*w58;
                                    EM_F[INDEX2(k,1,numEq)]+=216*Y_p[k]*w58;
                                    EM_F[INDEX2(k,2,numEq)]+=216*Y_p[k]*w58;
                                    EM_F[INDEX2(k,3,numEq)]+=216*Y_p[k]*w58;
                                    EM_F[INDEX2(k,4,numEq)]+=216*Y_p[k]*w58;
                                    EM_F[INDEX2(k,5,numEq)]+=216*Y_p[k]*w58;
                                    EM_F[INDEX2(k,6,numEq)]+=216*Y_p[k]*w58;
                                    EM_F[INDEX2(k,7,numEq)]+=216*Y_p[k]*w58;
                                }
                            }
                        }

                        // add to matrix (if add_EM_S) and RHS (if add_EM_F)
                        const index_t firstNode=m_NN[0]*m_NN[1]*k2+m_NN[0]*k1+k0;
                        domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F,
                                add_EM_S, add_EM_F, firstNode, numEq, numComp);
                    } // end k0 loop
                } // end k1 loop
            } // end k2 loop
        } // end of colouring
    } // end of parallel region
}

} // namespace ripley

