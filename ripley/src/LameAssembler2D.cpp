
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

#include <ripley/LameAssembler2D.h>
#include <ripley/domainhelpers.h>

#include <escript/index.h>

using namespace std;
using escript::AbstractSystemMatrix;
using escript::Data;

namespace ripley {

void LameAssembler2D::collateFunctionSpaceTypes(vector<int>& fsTypes,
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

void LameAssembler2D::assemblePDESingle(AbstractSystemMatrix* mat, Data& rhs,
                                        const DataMap& coefs) const
{
    throw RipleyException("assemblePDESingle not implemented in LameAssembler2D");
}

void LameAssembler2D::assemblePDEBoundarySingle(AbstractSystemMatrix* mat,
                                        Data& rhs, const DataMap& coefs) const 
{
    throw RipleyException("assemblePDEBoundarySingle not implemented in LameAssembler2D");
}

void LameAssembler2D::assemblePDESingleReduced(AbstractSystemMatrix* mat,
                                        Data& rhs, const DataMap& coefs) const
{
    throw RipleyException("assemblePDESingleReduced not implemented in LameAssembler2D");
}

void LameAssembler2D::assemblePDEBoundarySingleReduced(
                                         AbstractSystemMatrix* mat, Data& rhs,
                                         const DataMap& coefs) const
{
    throw RipleyException("assemblePDEBoundarySingleReduced not implemented in LameAssembler2D");
}

void LameAssembler2D::assemblePDESystemReduced(AbstractSystemMatrix* mat,
                                        Data& rhs, const DataMap& coefs) const
{
    throw RipleyException("assemblePDEBoundarySystem not implemented in LameAssembler2D");
}

void LameAssembler2D::assemblePDEBoundarySystemReduced(
                                        AbstractSystemMatrix* mat,
                                        Data& rhs, const DataMap& coefs) const
{
    throw RipleyException("assemblePDEBoundarySystemReduced not implemented in LameAssembler2D");
}

void LameAssembler2D::assemblePDEBoundarySystem(AbstractSystemMatrix* mat,
                                        Data& rhs, const DataMap& coefs) const
{
    const Data& d = unpackData("d", coefs);
    const Data& y = unpackData("y", coefs);
    dim_t numEq, numComp;
    if (!mat) {
        numEq=numComp=(rhs.isEmpty() ? 1 : rhs.getDataPointSize());
    } else {
        numEq=mat->getRowBlockSize();
        numComp=mat->getColumnBlockSize();
    }
    const double SQRT3 = 1.73205080756887719318;
    const double w5 = m_dx[0]/12;
    const double w6 = w5*(SQRT3 + 2);
    const double w7 = w5*(-SQRT3 + 2);
    const double w8 = w5*(SQRT3 + 3);
    const double w9 = w5*(-SQRT3 + 3);
    const double w2 = m_dx[1]/12;
    const double w0 = w2*(SQRT3 + 2);
    const double w1 = w2*(-SQRT3 + 2);
    const double w3 = w2*(SQRT3 + 3);
    const double w4 = w2*(-SQRT3 + 3);
    const bool addEM_S = !d.isEmpty();
    const bool addEM_F = !y.isEmpty();
    rhs.requireWrite();

#pragma omp parallel
    {
        vector<double> EM_S(4*4*numEq*numComp, 0);
        vector<double> EM_F(4*numEq, 0);

        if (domain->m_faceOffset[0] > -1) {
            if (addEM_S)
                fill(EM_S.begin(), EM_S.end(), 0);
            if (addEM_F)
                fill(EM_F.begin(), EM_F.end(), 0);

            for (index_t k1_0=0; k1_0<2; k1_0++) { // colouring
#pragma omp for
                for (index_t k1=k1_0; k1<m_NE[1]; k1+=2) {
                    const index_t e = k1;
                    ///////////////
                    // process d //
                    ///////////////
                    if (addEM_S) {
                        const double* d_p=d.getSampleDataRO(e);
                        if (d.actsExpanded()) {
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const double d_0 = d_p[INDEX3(k,m,0,numEq,numComp)];
                                    const double d_1 = d_p[INDEX3(k,m,1,numEq,numComp)];
                                    const double tmp0 = w2*(d_0 + d_1);
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,4)] = d_0*w0 + d_1*w1;
                                    EM_S[INDEX4(k,m,2,0,numEq,numComp,4)] = tmp0;
                                    EM_S[INDEX4(k,m,0,2,numEq,numComp,4)] = tmp0;
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,4)] = d_0*w1 + d_1*w0;
                                }
                             }
                        } else { // constant data
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const double d_0 = d_p[INDEX2(k, m, numEq)];
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,4)] = 4*d_0*w2;
                                    EM_S[INDEX4(k,m,2,0,numEq,numComp,4)] = 2*d_0*w2;
                                    EM_S[INDEX4(k,m,0,2,numEq,numComp,4)] = 2*d_0*w2;
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,4)] = 4*d_0*w2;
                                }
                            }
                        }
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (addEM_F) {
                        const double* y_p=y.getSampleDataRO(e);
                        if (y.actsExpanded()) {
                            for (index_t k=0; k<numEq; k++) {
                                const double y_0 = y_p[INDEX2(k, 0, numEq)];
                                const double y_1 = y_p[INDEX2(k, 1, numEq)];
                                EM_F[INDEX2(k,0,numEq)] = w3*y_0 + w4*y_1;
                                EM_F[INDEX2(k,2,numEq)] = w3*y_1 + w4*y_0;
                            }
                        } else { // constant data
                            for (index_t k=0; k<numEq; k++) {
                                EM_F[INDEX2(k,0,numEq)] = 6*w2*y_p[k];
                                EM_F[INDEX2(k,2,numEq)] = 6*w2*y_p[k];
                            }
                        }
                    }
                    const index_t firstNode=m_NN[0]*k1;
                    domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, addEM_S,
                                          addEM_F, firstNode, numEq, numComp);
                }
            } // end colouring
        }

        if (domain->m_faceOffset[1] > -1) {
            if (addEM_S)
                fill(EM_S.begin(), EM_S.end(), 0);
            if (addEM_F)
                fill(EM_F.begin(), EM_F.end(), 0);

            for (index_t k1_0=0; k1_0<2; k1_0++) { // colouring            
#pragma omp for
                for (index_t k1=k1_0; k1<m_NE[1]; k1+=2) {
                    const index_t e = domain->m_faceOffset[1]+k1;
                    ///////////////
                    // process d //
                    ///////////////
                    if (addEM_S) {
                        const double* d_p=d.getSampleDataRO(e);
                        if (d.actsExpanded()) {
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const double d_0 = d_p[INDEX3(k,m,0,numEq,numComp)];
                                    const double d_1 = d_p[INDEX3(k,m,1,numEq,numComp)];
                                    const double tmp0 = w2*(d_0 + d_1);
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,4)] = d_0*w0 + d_1*w1;
                                    EM_S[INDEX4(k,m,3,1,numEq,numComp,4)] = tmp0;
                                    EM_S[INDEX4(k,m,1,3,numEq,numComp,4)] = tmp0;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,4)] = d_0*w1 + d_1*w0;
                                }
                             }
                        } else { // constant data
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const double d_0 = d_p[INDEX2(k, m, numEq)];
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,4)] = 4*d_0*w2;
                                    EM_S[INDEX4(k,m,3,1,numEq,numComp,4)] = 2*d_0*w2;
                                    EM_S[INDEX4(k,m,1,3,numEq,numComp,4)] = 2*d_0*w2;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,4)] = 4*d_0*w2;
                                }
                            }
                        }
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (addEM_F) {
                        const double* y_p=y.getSampleDataRO(e);
                        if (y.actsExpanded()) {
                            for (index_t k=0; k<numEq; k++) {
                                const double y_0 = y_p[INDEX2(k, 0, numEq)];
                                const double y_1 = y_p[INDEX2(k, 1, numEq)];
                                EM_F[INDEX2(k,1,numEq)] = w3*y_0 + w4*y_1;
                                EM_F[INDEX2(k,3,numEq)] = w3*y_1 + w4*y_0;
                            }
                        } else { // constant data
                            for (index_t k=0; k<numEq; k++) {
                                EM_F[INDEX2(k,1,numEq)] = 6*w2*y_p[k];
                                EM_F[INDEX2(k,3,numEq)] = 6*w2*y_p[k];
                            }
                        }
                    }
                    const index_t firstNode=m_NN[0]*(k1+1)-2;
                    domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, addEM_S,
                                          addEM_F, firstNode, numEq, numComp);
                }
            } // end colouring
        }

        if (domain->m_faceOffset[2] > -1) {
            if (addEM_S)
                fill(EM_S.begin(), EM_S.end(), 0);
            if (addEM_F)
                fill(EM_F.begin(), EM_F.end(), 0);

            for (index_t k0_0=0; k0_0<2; k0_0++) { // colouring
#pragma omp for
                for (index_t k0 = k0_0; k0 < m_NE[0]; k0+=2) {
                    const index_t e = domain->m_faceOffset[2]+k0;
                    ///////////////
                    // process d //
                    ///////////////
                    if (addEM_S) {
                        const double* d_p=d.getSampleDataRO(e);
                        if (d.actsExpanded()) {
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const double d_0 = d_p[INDEX3(k,m,0,numEq,numComp)];
                                    const double d_1 = d_p[INDEX3(k,m,1,numEq,numComp)];
                                    const double tmp0 = w5*(d_0 + d_1);
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,4)] = d_0*w6 + d_1*w7;
                                    EM_S[INDEX4(k,m,1,0,numEq,numComp,4)] = tmp0;
                                    EM_S[INDEX4(k,m,0,1,numEq,numComp,4)] = tmp0;
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,4)] = d_0*w7 + d_1*w6;
                                }
                             }
                        } else { // constant data
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const double d_0 = d_p[INDEX2(k, m, numEq)];
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,4)] = 4*d_0*w5;
                                    EM_S[INDEX4(k,m,1,0,numEq,numComp,4)] = 2*d_0*w5;
                                    EM_S[INDEX4(k,m,0,1,numEq,numComp,4)] = 2*d_0*w5;
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,4)] = 4*d_0*w5;
                                }
                            }
                        }
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (addEM_F) {
                        const double* y_p=y.getSampleDataRO(e);
                        if (y.actsExpanded()) {
                            for (index_t k=0; k<numEq; k++) {
                                const double y_0 = y_p[INDEX2(k, 0, numEq)];
                                const double y_1 = y_p[INDEX2(k, 1, numEq)];
                                EM_F[INDEX2(k,0,numEq)] = w8*y_0 + w9*y_1;
                                EM_F[INDEX2(k,1,numEq)] = w8*y_1 + w9*y_0;
                            }
                        } else { // constant data
                            for (index_t k=0; k<numEq; k++) {
                                EM_F[INDEX2(k,0,numEq)] = 6*w5*y_p[k];
                                EM_F[INDEX2(k,1,numEq)] = 6*w5*y_p[k];
                            }
                        }
                    }
                    const index_t firstNode=k0;
                    domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, addEM_S,
                                          addEM_F, firstNode, numEq, numComp);
                }
            } // end colouring
        }

        if (domain->m_faceOffset[3] > -1) {
            if (addEM_S)
                fill(EM_S.begin(), EM_S.end(), 0);
            if (addEM_F)
                fill(EM_F.begin(), EM_F.end(), 0);

            for (index_t k0_0=0; k0_0<2; k0_0++) { // colouring
#pragma omp for
                for (index_t k0 = k0_0; k0 < m_NE[0]; k0+=2) {
                    const index_t e = domain->m_faceOffset[3]+k0;
                    ///////////////
                    // process d //
                    ///////////////
                    if (addEM_S) {
                        const double* d_p=d.getSampleDataRO(e);
                        if (d.actsExpanded()) {
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const double d_0 = d_p[INDEX3(k,m,0,numEq,numComp)];
                                    const double d_1 = d_p[INDEX3(k,m,1,numEq,numComp)];
                                    const double tmp0 = w5*(d_0 + d_1);
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,4)] = d_0*w6 + d_1*w7;
                                    EM_S[INDEX4(k,m,3,2,numEq,numComp,4)] = tmp0;
                                    EM_S[INDEX4(k,m,2,3,numEq,numComp,4)] = tmp0;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,4)] = d_0*w7 + d_1*w6;
                                }
                             }
                        } else { // constant data
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const double d_0 = d_p[INDEX2(k, m, numEq)];
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,4)] = 4*d_0*w5;
                                    EM_S[INDEX4(k,m,3,2,numEq,numComp,4)] = 2*d_0*w5;
                                    EM_S[INDEX4(k,m,2,3,numEq,numComp,4)] = 2*d_0*w5;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,4)] = 4*d_0*w5;
                                }
                            }
                        }
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (addEM_F) {
                        const double* y_p=y.getSampleDataRO(e);
                        if (y.actsExpanded()) {
                            for (index_t k=0; k<numEq; k++) {
                                const double y_0 = y_p[INDEX2(k, 0, numEq)];
                                const double y_1 = y_p[INDEX2(k, 1, numEq)];
                                EM_F[INDEX2(k,2,numEq)] = w8*y_0 + w9*y_1;
                                EM_F[INDEX2(k,3,numEq)] = w8*y_1 + w9*y_0;
                            }
                        } else { // constant data
                            for (index_t k=0; k<numEq; k++) {
                                EM_F[INDEX2(k,2,numEq)] = 6*w5*y_p[k];
                                EM_F[INDEX2(k,3,numEq)] = 6*w5*y_p[k];
                            }
                        }
                    }
                    const index_t firstNode=m_NN[0]*(m_NN[1]-2)+k0;
                    domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, addEM_S,
                                          addEM_F, firstNode, numEq, numComp);
                }
            } // end colouring
        }
    } // end of parallel section
}

void LameAssembler2D::assemblePDESystem(AbstractSystemMatrix* mat,
            Data& rhs, const DataMap& coefs) const
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
    const double w1 = 1.0/24;
    const double w5 = -SQRT3/24 + 1.0/12;
    const double w2 = -SQRT3/24 - 1.0/12;
    const double w19 = -m_dx[0]/12;
    const double w11 = w19*(SQRT3 + 3)/12;
    const double w14 = w19*(-SQRT3 + 3)/12;
    const double w16 = w19*(5*SQRT3 + 9)/12;
    const double w17 = w19*(-5*SQRT3 + 9)/12;
    const double w27 = w19*(-SQRT3 - 3)/2;
    const double w28 = w19*(SQRT3 - 3)/2;
    const double w18 = -m_dx[1]/12;
    const double w10 = w18*(SQRT3 + 3)/12;
    const double w15 = w18*(-SQRT3 + 3)/12;
    const double w12 = w18*(5*SQRT3 + 9)/12;
    const double w13 = w18*(-5*SQRT3 + 9)/12;
    const double w25 = w18*(-SQRT3 - 3)/2;
    const double w26 = w18*(SQRT3 - 3)/2;
    const double w22 = m_dx[0]*m_dx[1]/144;
    const double w20 = w22*(SQRT3 + 2);
    const double w21 = w22*(-SQRT3 + 2);
    const double w23 = w22*(4*SQRT3 + 7);
    const double w24 = w22*(-4*SQRT3 + 7);
    const double w3 = m_dx[0]/(24*m_dx[1]);
    const double w7 = w3*(SQRT3 + 2);
    const double w8 = w3*(-SQRT3 + 2);
    const double w6 = -m_dx[1]/(24*m_dx[0]);
    const double w0 = w6*(SQRT3 + 2);
    const double w4 = w6*(-SQRT3 + 2);
    const bool addEM_S = (!mu.isEmpty() || !lambda.isEmpty() || !B.isEmpty()
                                        || !C.isEmpty() || !D.isEmpty());
    const bool addEM_F = (!X.isEmpty() || !Y.isEmpty());
    rhs.requireWrite();

#pragma omp parallel
    {
        vector<double> EM_S(4*4*numEq*numComp, 0);
        vector<double> EM_F(4*numEq, 0);

        for (index_t k1_0=0; k1_0<2; k1_0++) { // colouring
#pragma omp for
            for (index_t k1=k1_0; k1 < m_NE[1]; k1+=2) {
                for (index_t k0=0; k0 < m_NE[0]; ++k0)  {
                    const index_t e = k0 + m_NE[0]*k1;
                    if (addEM_S)
                        fill(EM_S.begin(), EM_S.end(), 0);
                    if (addEM_F)
                        fill(EM_F.begin(), EM_F.end(), 0);

                    ///////////////
                    // process A //
                    ///////////////
                    if (!mu.isEmpty() || !lambda.isEmpty()) {
                        if (mu.actsExpanded() || lambda.actsExpanded()) {
                            double A_0000[4] = {0};
                            double A_0011[4] = {0};
                            double A_0101[4] = {0};
                            double A_0110[4] = {0};
                            double A_1001[4] = {0};
                            double A_1010[4] = {0};
                            double A_1100[4] = {0};
                            double A_1111[4] = {0};
                            if (!mu.isEmpty()) {
                                const double *mu_p = mu.getSampleDataRO(e);
                                A_0000[0] += 2*mu_p[0];
                                A_0000[1] += 2*mu_p[1];
                                A_0000[2] += 2*mu_p[2];
                                A_0000[3] += 2*mu_p[3];
                                A_0110[0] += mu_p[0];
                                A_0101[0] += mu_p[0];
                                A_0110[1] += mu_p[1];
                                A_0101[1] += mu_p[1];
                                A_0110[2] += mu_p[2];
                                A_0101[2] += mu_p[2];
                                A_0110[3] += mu_p[3];
                                A_0101[3] += mu_p[3];
                                A_1001[0] += mu_p[0];
                                A_1010[0] += mu_p[0];
                                A_1001[1] += mu_p[1];
                                A_1010[1] += mu_p[1];
                                A_1001[2] += mu_p[2];
                                A_1010[2] += mu_p[2];
                                A_1001[3] += mu_p[3];
                                A_1010[3] += mu_p[3];
                                A_1111[0] += 2*mu_p[0];
                                A_1111[1] += 2*mu_p[1];
                                A_1111[2] += 2*mu_p[2];
                                A_1111[3] += 2*mu_p[3];
                            }
                            if (!lambda.isEmpty()) {
                                const double *lambda_p = lambda.getSampleDataRO(e);
                                A_0000[0] += lambda_p[0];
                                A_0000[1] += lambda_p[1];
                                A_0000[2] += lambda_p[2];
                                A_0000[3] += lambda_p[3];
                                A_0011[0] += lambda_p[0];
                                A_0011[1] += lambda_p[1];
                                A_0011[2] += lambda_p[2];
                                A_0011[3] += lambda_p[3];
                                A_1100[0] += lambda_p[0];
                                A_1100[1] += lambda_p[1];
                                A_1100[2] += lambda_p[2];
                                A_1100[3] += lambda_p[3];
                                A_1111[0] += lambda_p[0];
                                A_1111[1] += lambda_p[1];
                                A_1111[2] += lambda_p[2];
                                A_1111[3] += lambda_p[3];
                            }
                            {
                                const double tmp0 = w3*(A_0101[0] + A_0101[1] + A_0101[2] + A_0101[3]);
                                const double tmp2 = w4*(A_0000[2] + A_0000[3]);
                                const double tmp3 = w0*(A_0000[0] + A_0000[1]);
                                const double tmp7 = w3*(-A_0101[0] - A_0101[1] - A_0101[2] - A_0101[3]);
                                const double tmp8 = w6*(A_0000[0] + A_0000[1] + A_0000[2] + A_0000[3]);
                                const double tmp11 = w4*(A_0000[0] + A_0000[1]);
                                const double tmp12 = w0*(A_0000[2] + A_0000[3]);
                                const double tmp15 = w7*(A_0101[0] + A_0101[2]);
                                const double tmp16 = w4*(-A_0000[2] - A_0000[3]);
                                const double tmp17 = w0*(-A_0000[0] - A_0000[1]);
                                const double tmp19 = w8*(A_0101[1] + A_0101[3]);
                                const double tmp21 = w7*(A_0101[1] + A_0101[3]);
                                const double tmp22 = w4*(-A_0000[0] - A_0000[1]);
                                const double tmp23 = w0*(-A_0000[2] - A_0000[3]);
                                const double tmp25 = w8*(A_0101[0] + A_0101[2]);
                                const double tmp30 = w7*(-A_0101[1] - A_0101[3]);
                                const double tmp33 = w8*(-A_0101[0] - A_0101[2]);
                                const double tmp34 = w6*(-A_0000[0] - A_0000[1] - A_0000[2] - A_0000[3]);
                                const double tmp38 = w7*(-A_0101[0] - A_0101[2]);
                                const double tmp40 = w8*(-A_0101[1] - A_0101[3]);
                                EM_S[INDEX4(0,0,0,0,numEq,numComp,4)]+= tmp15 + tmp16 + tmp17 + tmp19;
                                EM_S[INDEX4(0,0,0,1,numEq,numComp,4)]+= tmp0 + tmp2 + tmp3;
                                EM_S[INDEX4(0,0,0,2,numEq,numComp,4)]+= tmp34 + tmp38 + tmp40;
                                EM_S[INDEX4(0,0,0,3,numEq,numComp,4)]+= tmp7 + tmp8;
                                EM_S[INDEX4(0,0,1,0,numEq,numComp,4)]+= tmp0 + tmp2 + tmp3;
                                EM_S[INDEX4(0,0,1,1,numEq,numComp,4)]+= tmp16 + tmp17 + tmp21 + tmp25;
                                EM_S[INDEX4(0,0,1,2,numEq,numComp,4)]+= tmp7 + tmp8;
                                EM_S[INDEX4(0,0,1,3,numEq,numComp,4)]+= tmp30 + tmp33 + tmp34;
                                EM_S[INDEX4(0,0,2,0,numEq,numComp,4)]+= tmp34 + tmp38 + tmp40;
                                EM_S[INDEX4(0,0,2,1,numEq,numComp,4)]+= tmp7 + tmp8;
                                EM_S[INDEX4(0,0,2,2,numEq,numComp,4)]+= tmp15 + tmp19 + tmp22 + tmp23;
                                EM_S[INDEX4(0,0,2,3,numEq,numComp,4)]+= tmp0 + tmp11 + tmp12;
                                EM_S[INDEX4(0,0,3,0,numEq,numComp,4)]+= tmp7 + tmp8;
                                EM_S[INDEX4(0,0,3,1,numEq,numComp,4)]+= tmp30 + tmp33 + tmp34;
                                EM_S[INDEX4(0,0,3,2,numEq,numComp,4)]+= tmp0 + tmp11 + tmp12;
                                EM_S[INDEX4(0,0,3,3,numEq,numComp,4)]+= tmp21 + tmp22 + tmp23 + tmp25;
                            }
                            {
                                const double tmp1 = w1*(A_0011[0] + A_0011[3] - A_0110[1] - A_0110[2]);
                                const double tmp9 = w1*(A_0011[1] + A_0011[2] + A_0110[1] + A_0110[2]);
                                const double tmp28 = w1*(-A_0011[0] - A_0011[3] - A_0110[0] - A_0110[3]);
                                const double tmp31 = w1*(-A_0011[1] - A_0011[2] + A_0110[0] + A_0110[3]);
                                EM_S[INDEX4(0,1,0,0,numEq,numComp,4)]+= w5*(A_0011[3] + A_0110[3]) + w2*(-A_0011[0] - A_0110[0]) + tmp9;
                                EM_S[INDEX4(0,1,0,1,numEq,numComp,4)]+= tmp1 + w5*(A_0011[2] - A_0110[3]) + w2*(-A_0011[1] + A_0110[0]);
                                EM_S[INDEX4(0,1,0,2,numEq,numComp,4)]+= tmp31 + w5*(-A_0011[3] + A_0110[1]) + w2*(A_0011[0] - A_0110[2]);
                                EM_S[INDEX4(0,1,0,3,numEq,numComp,4)]+= tmp28 + w5*(-A_0011[2] - A_0110[1]) + w2*(A_0011[1] + A_0110[2]);
                                EM_S[INDEX4(0,1,1,0,numEq,numComp,4)]+= tmp31 + w5*(-A_0011[3] + A_0110[2]) + w2*(A_0011[0] - A_0110[1]);
                                EM_S[INDEX4(0,1,1,1,numEq,numComp,4)]+= tmp28 + w5*(-A_0011[2] - A_0110[2]) + w2*(A_0011[1] + A_0110[1]);
                                EM_S[INDEX4(0,1,1,2,numEq,numComp,4)]+= w2*(-A_0011[0] - A_0110[3]) + w5*(A_0011[3] + A_0110[0]) + tmp9;
                                EM_S[INDEX4(0,1,1,3,numEq,numComp,4)]+= tmp1 + w5*(A_0011[2] - A_0110[0]) + w2*(-A_0011[1] + A_0110[3]);
                                EM_S[INDEX4(0,1,2,0,numEq,numComp,4)]+= tmp1 + w5*(A_0011[1] - A_0110[3]) + w2*(-A_0011[2] + A_0110[0]);
                                EM_S[INDEX4(0,1,2,1,numEq,numComp,4)]+= w5*(A_0011[0] + A_0110[3]) + w2*(-A_0011[3] - A_0110[0]) + tmp9;
                                EM_S[INDEX4(0,1,2,2,numEq,numComp,4)]+= tmp28 + w5*(-A_0011[1] - A_0110[1]) + w2*(A_0011[2] + A_0110[2]);
                                EM_S[INDEX4(0,1,2,3,numEq,numComp,4)]+= tmp31 + w5*(-A_0011[0] + A_0110[1]) + w2*(A_0011[3] - A_0110[2]);
                                EM_S[INDEX4(0,1,3,0,numEq,numComp,4)]+= w5*(-A_0011[1] - A_0110[2]) + tmp28 + w2*(A_0011[2] + A_0110[1]);
                                EM_S[INDEX4(0,1,3,1,numEq,numComp,4)]+= tmp31 + w5*(-A_0011[0] + A_0110[2]) + w2*(A_0011[3] - A_0110[1]);
                                EM_S[INDEX4(0,1,3,2,numEq,numComp,4)]+= tmp1 + w5*(A_0011[1] - A_0110[0]) + w2*(-A_0011[2] + A_0110[3]);
                                EM_S[INDEX4(0,1,3,3,numEq,numComp,4)]+= w5*(A_0011[0] + A_0110[0]) + w2*(-A_0011[3] - A_0110[3]) + tmp9;
                            }
                            {
                                const double tmp1 = w1*(A_1001[0] + A_1001[3] - A_1100[1] - A_1100[2]);
                                const double tmp9 = w1*(A_1001[1] + A_1001[2] + A_1100[1] + A_1100[2]);
                                const double tmp28 = w1*(-A_1001[0] - A_1001[3] - A_1100[0] - A_1100[3]);
                                const double tmp31 = w1*(-A_1001[1] - A_1001[2] + A_1100[0] + A_1100[3]);
                                EM_S[INDEX4(1,0,0,0,numEq,numComp,4)]+= w5*(A_1001[3] + A_1100[3]) + w2*(-A_1001[0] - A_1100[0]) + tmp9;
                                EM_S[INDEX4(1,0,0,1,numEq,numComp,4)]+= tmp1 + w5*(A_1001[2] - A_1100[3]) + w2*(-A_1001[1] + A_1100[0]);
                                EM_S[INDEX4(1,0,0,2,numEq,numComp,4)]+= tmp31 + w5*(-A_1001[3] + A_1100[1]) + w2*(A_1001[0] - A_1100[2]);
                                EM_S[INDEX4(1,0,0,3,numEq,numComp,4)]+= tmp28 + w5*(-A_1001[2] - A_1100[1]) + w2*(A_1001[1] + A_1100[2]);
                                EM_S[INDEX4(1,0,1,0,numEq,numComp,4)]+= tmp31 + w5*(-A_1001[3] + A_1100[2]) + w2*(A_1001[0] - A_1100[1]);
                                EM_S[INDEX4(1,0,1,1,numEq,numComp,4)]+= tmp28 + w5*(-A_1001[2] - A_1100[2]) + w2*(A_1001[1] + A_1100[1]);
                                EM_S[INDEX4(1,0,1,2,numEq,numComp,4)]+= w2*(-A_1001[0] - A_1100[3]) + w5*(A_1001[3] + A_1100[0]) + tmp9;
                                EM_S[INDEX4(1,0,1,3,numEq,numComp,4)]+= tmp1 + w5*(A_1001[2] - A_1100[0]) + w2*(-A_1001[1] + A_1100[3]);
                                EM_S[INDEX4(1,0,2,0,numEq,numComp,4)]+= tmp1 + w5*(A_1001[1] - A_1100[3]) + w2*(-A_1001[2] + A_1100[0]);
                                EM_S[INDEX4(1,0,2,1,numEq,numComp,4)]+= w5*(A_1001[0] + A_1100[3]) + w2*(-A_1001[3] - A_1100[0]) + tmp9;
                                EM_S[INDEX4(1,0,2,2,numEq,numComp,4)]+= tmp28 + w5*(-A_1001[1] - A_1100[1]) + w2*(A_1001[2] + A_1100[2]);
                                EM_S[INDEX4(1,0,2,3,numEq,numComp,4)]+= tmp31 + w5*(-A_1001[0] + A_1100[1]) + w2*(A_1001[3] - A_1100[2]);
                                EM_S[INDEX4(1,0,3,0,numEq,numComp,4)]+= w5*(-A_1001[1] - A_1100[2]) + tmp28 + w2*(A_1001[2] + A_1100[1]);
                                EM_S[INDEX4(1,0,3,1,numEq,numComp,4)]+= tmp31 + w5*(-A_1001[0] + A_1100[2]) + w2*(A_1001[3] - A_1100[1]);
                                EM_S[INDEX4(1,0,3,2,numEq,numComp,4)]+= tmp1 + w5*(A_1001[1] - A_1100[0]) + w2*(-A_1001[2] + A_1100[3]);
                                EM_S[INDEX4(1,0,3,3,numEq,numComp,4)]+= w5*(A_1001[0] + A_1100[0]) + w2*(-A_1001[3] - A_1100[3]) + tmp9;
                            }
                            {
                                const double tmp0 = w3*(A_1111[0] + A_1111[1] + A_1111[2] + A_1111[3]);
                                const double tmp2 = w4*(A_1010[2] + A_1010[3]);
                                const double tmp3 = w0*(A_1010[0] + A_1010[1]);
                                const double tmp7 = w3*(-A_1111[0] - A_1111[1] - A_1111[2] - A_1111[3]);
                                const double tmp8 = w6*(A_1010[0] + A_1010[1] + A_1010[2] + A_1010[3]);
                                const double tmp11 = w4*(A_1010[0] + A_1010[1]);
                                const double tmp12 = w0*(A_1010[2] + A_1010[3]);
                                const double tmp15 = w7*(A_1111[0] + A_1111[2]);
                                const double tmp16 = w4*(-A_1010[2] - A_1010[3]);
                                const double tmp17 = w0*(-A_1010[0] - A_1010[1]);
                                const double tmp19 = w8*(A_1111[1] + A_1111[3]);
                                const double tmp21 = w7*(A_1111[1] + A_1111[3]);
                                const double tmp22 = w4*(-A_1010[0] - A_1010[1]);
                                const double tmp23 = w0*(-A_1010[2] - A_1010[3]);
                                const double tmp25 = w8*(A_1111[0] + A_1111[2]);
                                const double tmp30 = w7*(-A_1111[1] - A_1111[3]);
                                const double tmp33 = w8*(-A_1111[0] - A_1111[2]);
                                const double tmp34 = w6*(-A_1010[0] - A_1010[1] - A_1010[2] - A_1010[3]);
                                const double tmp38 = w7*(-A_1111[0] - A_1111[2]);
                                const double tmp40 = w8*(-A_1111[1] - A_1111[3]);
                                EM_S[INDEX4(1,1,0,0,numEq,numComp,4)]+= tmp15 + tmp16 + tmp17 + tmp19;
                                EM_S[INDEX4(1,1,0,1,numEq,numComp,4)]+= tmp0 + tmp2 + tmp3;
                                EM_S[INDEX4(1,1,0,2,numEq,numComp,4)]+= tmp34 + tmp38 + tmp40;
                                EM_S[INDEX4(1,1,0,3,numEq,numComp,4)]+= tmp7 + tmp8;
                                EM_S[INDEX4(1,1,1,0,numEq,numComp,4)]+= tmp0 + tmp2 + tmp3;
                                EM_S[INDEX4(1,1,1,1,numEq,numComp,4)]+= tmp16 + tmp17 + tmp21 + tmp25;
                                EM_S[INDEX4(1,1,1,2,numEq,numComp,4)]+= tmp7 + tmp8;
                                EM_S[INDEX4(1,1,1,3,numEq,numComp,4)]+= tmp30 + tmp33 + tmp34;
                                EM_S[INDEX4(1,1,2,0,numEq,numComp,4)]+= tmp34 + tmp38 + tmp40;
                                EM_S[INDEX4(1,1,2,1,numEq,numComp,4)]+= tmp7 + tmp8;
                                EM_S[INDEX4(1,1,2,2,numEq,numComp,4)]+= tmp15 + tmp19 + tmp22 + tmp23;
                                EM_S[INDEX4(1,1,2,3,numEq,numComp,4)]+= tmp0 + tmp11 + tmp12;
                                EM_S[INDEX4(1,1,3,0,numEq,numComp,4)]+= tmp7 + tmp8;
                                EM_S[INDEX4(1,1,3,1,numEq,numComp,4)]+= tmp30 + tmp33 + tmp34;
                                EM_S[INDEX4(1,1,3,2,numEq,numComp,4)]+= tmp0 + tmp11 + tmp12;
                                EM_S[INDEX4(1,1,3,3,numEq,numComp,4)]+= tmp21 + tmp22 + tmp23 + tmp25;
                            }
                        } else { // constant data
                            double A_0000 = 0;
                            double A_0011 = 0;
                            double A_0101 = 0;
                            double A_0110 = 0;
                            double A_1001 = 0;
                            double A_1010 = 0;
                            double A_1100 = 0;
                            double A_1111 = 0;
                            if (!mu.isEmpty()) {
                                const double *mu_p = mu.getSampleDataRO(e);
                                A_0000 += 2*mu_p[0];
                                A_0110 += mu_p[0];
                                A_0101 += mu_p[0];
                                A_1001 += mu_p[0];
                                A_1010 += mu_p[0];
                                A_1111 += 2*mu_p[0];
                            }
                            if (!lambda.isEmpty()) {
                                const double *lambda_p = lambda.getSampleDataRO(e);
                                A_0000 += lambda_p[0];
                                A_0011 += lambda_p[0];
                                A_1100 += lambda_p[0];
                                A_1111 += lambda_p[0];
                            }

                            const double tmp01_0 = 6*w1*(A_0011 - A_0110);
                            const double tmp01_1 = 6*w1*(A_0011 + A_0110);
                            const double tmp01_2 = 6*w1*(-A_0011 - A_0110);
                            const double tmp01_3 = 6*w1*(-A_0011 + A_0110);
                            const double tmp10_0 = 6*w1*(A_1001 - A_1100);
                            const double tmp10_1 = 6*w1*(A_1001 + A_1100);
                            const double tmp10_2 = 6*w1*(-A_1001 - A_1100);
                            const double tmp10_3 = 6*w1*(-A_1001 + A_1100);
                            EM_S[INDEX4(0,0,0,0,numEq,numComp,4)]+=-8*A_0000*w6 + 8*A_0101*w3;
                            EM_S[INDEX4(0,0,0,1,numEq,numComp,4)]+= 8*A_0000*w6 + 4*A_0101*w3;
                            EM_S[INDEX4(0,0,0,2,numEq,numComp,4)]+=-4*A_0000*w6 - 8*A_0101*w3;
                            EM_S[INDEX4(0,0,0,3,numEq,numComp,4)]+= 4*A_0000*w6 - 4*A_0101*w3;
                            EM_S[INDEX4(0,0,1,0,numEq,numComp,4)]+= 8*A_0000*w6 + 4*A_0101*w3;
                            EM_S[INDEX4(0,0,1,1,numEq,numComp,4)]+=-8*A_0000*w6 + 8*A_0101*w3;
                            EM_S[INDEX4(0,0,1,2,numEq,numComp,4)]+= 4*A_0000*w6 - 4*A_0101*w3;
                            EM_S[INDEX4(0,0,1,3,numEq,numComp,4)]+=-4*A_0000*w6 - 8*A_0101*w3;
                            EM_S[INDEX4(0,0,2,0,numEq,numComp,4)]+=-4*A_0000*w6 - 8*A_0101*w3;
                            EM_S[INDEX4(0,0,2,1,numEq,numComp,4)]+= 4*A_0000*w6 - 4*A_0101*w3;
                            EM_S[INDEX4(0,0,2,2,numEq,numComp,4)]+=-8*A_0000*w6 + 8*A_0101*w3;
                            EM_S[INDEX4(0,0,2,3,numEq,numComp,4)]+= 8*A_0000*w6 + 4*A_0101*w3;
                            EM_S[INDEX4(0,0,3,0,numEq,numComp,4)]+= 4*A_0000*w6 - 4*A_0101*w3;
                            EM_S[INDEX4(0,0,3,1,numEq,numComp,4)]+=-4*A_0000*w6 - 8*A_0101*w3;
                            EM_S[INDEX4(0,0,3,2,numEq,numComp,4)]+= 8*A_0000*w6 + 4*A_0101*w3;
                            EM_S[INDEX4(0,0,3,3,numEq,numComp,4)]+=-8*A_0000*w6 + 8*A_0101*w3;
                            EM_S[INDEX4(0,1,0,0,numEq,numComp,4)]+= tmp01_1;
                            EM_S[INDEX4(0,1,0,1,numEq,numComp,4)]+= tmp01_0;
                            EM_S[INDEX4(0,1,0,2,numEq,numComp,4)]+= tmp01_3;
                            EM_S[INDEX4(0,1,0,3,numEq,numComp,4)]+= tmp01_2;
                            EM_S[INDEX4(0,1,1,0,numEq,numComp,4)]+= tmp01_3;
                            EM_S[INDEX4(0,1,1,1,numEq,numComp,4)]+= tmp01_2;
                            EM_S[INDEX4(0,1,1,2,numEq,numComp,4)]+= tmp01_1;
                            EM_S[INDEX4(0,1,1,3,numEq,numComp,4)]+= tmp01_0;
                            EM_S[INDEX4(0,1,2,0,numEq,numComp,4)]+= tmp01_0;
                            EM_S[INDEX4(0,1,2,1,numEq,numComp,4)]+= tmp01_1;
                            EM_S[INDEX4(0,1,2,2,numEq,numComp,4)]+= tmp01_2;
                            EM_S[INDEX4(0,1,2,3,numEq,numComp,4)]+= tmp01_3;
                            EM_S[INDEX4(0,1,3,0,numEq,numComp,4)]+= tmp01_2;
                            EM_S[INDEX4(0,1,3,1,numEq,numComp,4)]+= tmp01_3;
                            EM_S[INDEX4(0,1,3,2,numEq,numComp,4)]+= tmp01_0;
                            EM_S[INDEX4(0,1,3,3,numEq,numComp,4)]+= tmp01_1;
                            EM_S[INDEX4(1,0,0,0,numEq,numComp,4)]+= tmp10_1;
                            EM_S[INDEX4(1,0,0,1,numEq,numComp,4)]+= tmp10_0;
                            EM_S[INDEX4(1,0,0,2,numEq,numComp,4)]+= tmp10_3;
                            EM_S[INDEX4(1,0,0,3,numEq,numComp,4)]+= tmp10_2;
                            EM_S[INDEX4(1,0,1,0,numEq,numComp,4)]+= tmp10_3;
                            EM_S[INDEX4(1,0,1,1,numEq,numComp,4)]+= tmp10_2;
                            EM_S[INDEX4(1,0,1,2,numEq,numComp,4)]+= tmp10_1;
                            EM_S[INDEX4(1,0,1,3,numEq,numComp,4)]+= tmp10_0;
                            EM_S[INDEX4(1,0,2,0,numEq,numComp,4)]+= tmp10_0;
                            EM_S[INDEX4(1,0,2,1,numEq,numComp,4)]+= tmp10_1;
                            EM_S[INDEX4(1,0,2,2,numEq,numComp,4)]+= tmp10_2;
                            EM_S[INDEX4(1,0,2,3,numEq,numComp,4)]+= tmp10_3;
                            EM_S[INDEX4(1,0,3,0,numEq,numComp,4)]+= tmp10_2;
                            EM_S[INDEX4(1,0,3,1,numEq,numComp,4)]+= tmp10_3;
                            EM_S[INDEX4(1,0,3,2,numEq,numComp,4)]+= tmp10_0;
                            EM_S[INDEX4(1,0,3,3,numEq,numComp,4)]+= tmp10_1;
                            EM_S[INDEX4(1,1,0,0,numEq,numComp,4)]+=-8*A_1010*w6 + 8*A_1111*w3;
                            EM_S[INDEX4(1,1,0,1,numEq,numComp,4)]+= 8*A_1010*w6 + 4*A_1111*w3;
                            EM_S[INDEX4(1,1,0,2,numEq,numComp,4)]+=-4*A_1010*w6 - 8*A_1111*w3;
                            EM_S[INDEX4(1,1,0,3,numEq,numComp,4)]+= 4*A_1010*w6 - 4*A_1111*w3;
                            EM_S[INDEX4(1,1,1,0,numEq,numComp,4)]+= 8*A_1010*w6 + 4*A_1111*w3;
                            EM_S[INDEX4(1,1,1,1,numEq,numComp,4)]+=-8*A_1010*w6 + 8*A_1111*w3;
                            EM_S[INDEX4(1,1,1,2,numEq,numComp,4)]+= 4*A_1010*w6 - 4*A_1111*w3;
                            EM_S[INDEX4(1,1,1,3,numEq,numComp,4)]+=-4*A_1010*w6 - 8*A_1111*w3;
                            EM_S[INDEX4(1,1,2,0,numEq,numComp,4)]+=-4*A_1010*w6 - 8*A_1111*w3;
                            EM_S[INDEX4(1,1,2,1,numEq,numComp,4)]+= 4*A_1010*w6 - 4*A_1111*w3;
                            EM_S[INDEX4(1,1,2,2,numEq,numComp,4)]+=-8*A_1010*w6 + 8*A_1111*w3;
                            EM_S[INDEX4(1,1,2,3,numEq,numComp,4)]+= 8*A_1010*w6 + 4*A_1111*w3;
                            EM_S[INDEX4(1,1,3,0,numEq,numComp,4)]+= 4*A_1010*w6 - 4*A_1111*w3;
                            EM_S[INDEX4(1,1,3,1,numEq,numComp,4)]+=-4*A_1010*w6 - 8*A_1111*w3;
                            EM_S[INDEX4(1,1,3,2,numEq,numComp,4)]+= 8*A_1010*w6 + 4*A_1111*w3;
                            EM_S[INDEX4(1,1,3,3,numEq,numComp,4)]+=-8*A_1010*w6 + 8*A_1111*w3;
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
                                    const double B_0_0 = B_p[INDEX4(k,0,m,0, numEq,2,numComp)];
                                    const double B_1_0 = B_p[INDEX4(k,1,m,0, numEq,2,numComp)];
                                    const double B_0_1 = B_p[INDEX4(k,0,m,1, numEq,2,numComp)];
                                    const double B_1_1 = B_p[INDEX4(k,1,m,1, numEq,2,numComp)];
                                    const double B_0_2 = B_p[INDEX4(k,0,m,2, numEq,2,numComp)];
                                    const double B_1_2 = B_p[INDEX4(k,1,m,2, numEq,2,numComp)];
                                    const double B_0_3 = B_p[INDEX4(k,0,m,3, numEq,2,numComp)];
                                    const double B_1_3 = B_p[INDEX4(k,1,m,3, numEq,2,numComp)];
                                    const double tmp0 = w11*(B_1_0 + B_1_1);
                                    const double tmp1 = w14*(B_1_2 + B_1_3);
                                    const double tmp2 = w15*(-B_0_1 - B_0_3);
                                    const double tmp3 = w10*(-B_0_0 - B_0_2);
                                    const double tmp4 = w11*(B_1_2 + B_1_3);
                                    const double tmp5 = w14*(B_1_0 + B_1_1);
                                    const double tmp6 = w11*(-B_1_2 - B_1_3);
                                    const double tmp7 = w14*(-B_1_0 - B_1_1);
                                    const double tmp8 = w11*(-B_1_0 - B_1_1);
                                    const double tmp9 = w14*(-B_1_2 - B_1_3);
                                    const double tmp10 = w10*(-B_0_1 - B_0_3);
                                    const double tmp11 = w15*(-B_0_0 - B_0_2);
                                    const double tmp12 = w15*(B_0_0 + B_0_2);
                                    const double tmp13 = w10*(B_0_1 + B_0_3);
                                    const double tmp14 = w10*(B_0_0 + B_0_2);
                                    const double tmp15 = w15*(B_0_1 + B_0_3);
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=B_0_0*w12 + B_0_1*w10 + B_0_2*w15 + B_0_3*w13 + B_1_0*w16 + B_1_1*w14 + B_1_2*w11 + B_1_3*w17;
                                    EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+=B_0_0*w10 + B_0_1*w12 + B_0_2*w13 + B_0_3*w15 + tmp0 + tmp1;
                                    EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=B_1_0*w11 + B_1_1*w17 + B_1_2*w16 + B_1_3*w14 + tmp14 + tmp15;
                                    EM_S[INDEX4(k,m,0,3,numEq,numComp,4)]+=tmp12 + tmp13 + tmp4 + tmp5;
                                    EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+=-B_0_0*w12 - B_0_1*w10 - B_0_2*w15 - B_0_3*w13 + tmp0 + tmp1;
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=-B_0_0*w10 - B_0_1*w12 - B_0_2*w13 - B_0_3*w15 + B_1_0*w14 + B_1_1*w16 + B_1_2*w17 + B_1_3*w11;
                                    EM_S[INDEX4(k,m,1,2,numEq,numComp,4)]+=tmp2 + tmp3 + tmp4 + tmp5;
                                    EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+=B_1_0*w17 + B_1_1*w11 + B_1_2*w14 + B_1_3*w16 + tmp10 + tmp11;
                                    EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=-B_1_0*w16 - B_1_1*w14 - B_1_2*w11 - B_1_3*w17 + tmp14 + tmp15;
                                    EM_S[INDEX4(k,m,2,1,numEq,numComp,4)]+=tmp12 + tmp13 + tmp8 + tmp9;
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=B_0_0*w15 + B_0_1*w13 + B_0_2*w12 + B_0_3*w10 - B_1_0*w11 - B_1_1*w17 - B_1_2*w16 - B_1_3*w14;
                                    EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+=B_0_0*w13 + B_0_1*w15 + B_0_2*w10 + B_0_3*w12 + tmp6 + tmp7;
                                    EM_S[INDEX4(k,m,3,0,numEq,numComp,4)]+=tmp2 + tmp3 + tmp8 + tmp9;
                                    EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+=-B_1_0*w14 - B_1_1*w16 - B_1_2*w17 - B_1_3*w11 + tmp10 + tmp11;
                                    EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+=-B_0_0*w15 - B_0_1*w13 - B_0_2*w12 - B_0_3*w10 + tmp6 + tmp7;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=-B_0_0*w13 - B_0_1*w15 - B_0_2*w10 - B_0_3*w12 - B_1_0*w17 - B_1_1*w11 - B_1_2*w14 - B_1_3*w16;
                                }
                            }
                        } else { // constant data
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const double wB0 = B_p[INDEX3(k,0,m,numEq,2)]*w18;
                                    const double wB1 = B_p[INDEX3(k,1,m,numEq,2)]*w19;
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+= 2*wB0 + 2*wB1;
                                    EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+= 2*wB0 +   wB1;
                                    EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=   wB0 + 2*wB1;
                                    EM_S[INDEX4(k,m,0,3,numEq,numComp,4)]+=   wB0 +   wB1;
                                    EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+=-2*wB0 +   wB1;
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=-2*wB0 + 2*wB1;
                                    EM_S[INDEX4(k,m,1,2,numEq,numComp,4)]+=  -wB0 +   wB1;
                                    EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+=  -wB0 + 2*wB1;
                                    EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=   wB0 - 2*wB1;
                                    EM_S[INDEX4(k,m,2,1,numEq,numComp,4)]+=   wB0 -   wB1;
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+= 2*wB0 - 2*wB1;
                                    EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+= 2*wB0 -   wB1;
                                    EM_S[INDEX4(k,m,3,0,numEq,numComp,4)]+=  -wB0 -   wB1;
                                    EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+=  -wB0 - 2*wB1;
                                    EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+=-2*wB0 -   wB1;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=-2*wB0 - 2*wB1;
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
                                    const double C_0_0 = C_p[INDEX4(k,m,0, 0, numEq,numComp,2)];
                                    const double C_1_0 = C_p[INDEX4(k,m,1, 0, numEq,numComp,2)];
                                    const double C_0_1 = C_p[INDEX4(k,m,0, 1, numEq,numComp,2)];
                                    const double C_1_1 = C_p[INDEX4(k,m,1, 1, numEq,numComp,2)];
                                    const double C_0_2 = C_p[INDEX4(k,m,0, 2, numEq,numComp,2)];
                                    const double C_1_2 = C_p[INDEX4(k,m,1, 2, numEq,numComp,2)];
                                    const double C_0_3 = C_p[INDEX4(k,m,0, 3, numEq,numComp,2)];
                                    const double C_1_3 = C_p[INDEX4(k,m,1, 3, numEq,numComp,2)];
                                    const double tmp0 = w11*(C_1_0 + C_1_1);
                                    const double tmp1 = w14*(C_1_2 + C_1_3);
                                    const double tmp2 = w15*(C_0_0 + C_0_2);
                                    const double tmp3 = w10*(C_0_1 + C_0_3);
                                    const double tmp4 = w11*(-C_1_0 - C_1_1);
                                    const double tmp5 = w14*(-C_1_2 - C_1_3);
                                    const double tmp6 = w11*(-C_1_2 - C_1_3);
                                    const double tmp7 = w14*(-C_1_0 - C_1_1);
                                    const double tmp8 = w11*(C_1_2 + C_1_3);
                                    const double tmp9 = w14*(C_1_0 + C_1_1);
                                    const double tmp10 = w10*(-C_0_1 - C_0_3);
                                    const double tmp11 = w15*(-C_0_0 - C_0_2);
                                    const double tmp12 = w15*(-C_0_1 - C_0_3);
                                    const double tmp13 = w10*(-C_0_0 - C_0_2);
                                    const double tmp14 = w10*(C_0_0 + C_0_2);
                                    const double tmp15 = w15*(C_0_1 + C_0_3);
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=C_0_0*w12 + C_0_1*w10 + C_0_2*w15 + C_0_3*w13 + C_1_0*w16 + C_1_1*w14 + C_1_2*w11 + C_1_3*w17;
                                    EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+=-C_0_0*w12 - C_0_1*w10 - C_0_2*w15 - C_0_3*w13 + tmp0 + tmp1;
                                    EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=-C_1_0*w16 - C_1_1*w14 - C_1_2*w11 - C_1_3*w17 + tmp14 + tmp15;
                                    EM_S[INDEX4(k,m,0,3,numEq,numComp,4)]+=tmp12 + tmp13 + tmp4 + tmp5;
                                    EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+=C_0_0*w10 + C_0_1*w12 + C_0_2*w13 + C_0_3*w15 + tmp0 + tmp1;
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=-C_0_0*w10 - C_0_1*w12 - C_0_2*w13 - C_0_3*w15 + C_1_0*w14 + C_1_1*w16 + C_1_2*w17 + C_1_3*w11;
                                    EM_S[INDEX4(k,m,1,2,numEq,numComp,4)]+=tmp2 + tmp3 + tmp4 + tmp5;
                                    EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+=-C_1_0*w14 - C_1_1*w16 - C_1_2*w17 - C_1_3*w11 + tmp10 + tmp11;
                                    EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=C_1_0*w11 + C_1_1*w17 + C_1_2*w16 + C_1_3*w14 + tmp14 + tmp15;
                                    EM_S[INDEX4(k,m,2,1,numEq,numComp,4)]+=tmp12 + tmp13 + tmp8 + tmp9;
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=C_0_0*w15 + C_0_1*w13 + C_0_2*w12 + C_0_3*w10 - C_1_0*w11 - C_1_1*w17 - C_1_2*w16 - C_1_3*w14;
                                    EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+=-C_0_0*w15 - C_0_1*w13 - C_0_2*w12 - C_0_3*w10 + tmp6 + tmp7;
                                    EM_S[INDEX4(k,m,3,0,numEq,numComp,4)]+=tmp2 + tmp3 + tmp8 + tmp9;
                                    EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+=C_1_0*w17 + C_1_1*w11 + C_1_2*w14 + C_1_3*w16 + tmp10 + tmp11;
                                    EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+=C_0_0*w13 + C_0_1*w15 + C_0_2*w10 + C_0_3*w12 + tmp6 + tmp7;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=-C_0_0*w13 - C_0_1*w15 - C_0_2*w10 - C_0_3*w12 - C_1_0*w17 - C_1_1*w11 - C_1_2*w14 - C_1_3*w16;
                                }
                            }
                        } else { // constant data
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const double wC0 = C_p[INDEX3(k,m,0,numEq,numComp)]*w18;
                                    const double wC1 = C_p[INDEX3(k,m,1,numEq,numComp)]*w19;
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+= 2*wC0 + 2*wC1;
                                    EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+=-2*wC0 +   wC1;
                                    EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=   wC0 - 2*wC1;
                                    EM_S[INDEX4(k,m,0,3,numEq,numComp,4)]+=  -wC0 -   wC1;
                                    EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+= 2*wC0 +   wC1;
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=-2*wC0 + 2*wC1;
                                    EM_S[INDEX4(k,m,1,2,numEq,numComp,4)]+=   wC0 -   wC1;
                                    EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+=  -wC0 - 2*wC1;
                                    EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=   wC0 + 2*wC1;
                                    EM_S[INDEX4(k,m,2,1,numEq,numComp,4)]+=  -wC0 +   wC1;
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+= 2*wC0 - 2*wC1;
                                    EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+=-2*wC0 -   wC1;
                                    EM_S[INDEX4(k,m,3,0,numEq,numComp,4)]+=   wC0 +   wC1;
                                    EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+=  -wC0 + 2*wC1;
                                    EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+= 2*wC0 -   wC1;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=-2*wC0 - 2*wC1;
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
                                    const double tmp0 = w21*(D_2 + D_3);
                                    const double tmp1 = w20*(D_0 + D_1);
                                    const double tmp2 = w22*(D_0 + D_1 + D_2 + D_3);
                                    const double tmp3 = w21*(D_0 + D_1);
                                    const double tmp4 = w20*(D_2 + D_3);
                                    const double tmp5 = w22*(D_1 + D_2);
                                    const double tmp6 = w21*(D_0 + D_2);
                                    const double tmp7 = w20*(D_1 + D_3);
                                    const double tmp8 = w21*(D_1 + D_3);
                                    const double tmp9 = w20*(D_0 + D_2);
                                    const double tmp10 = w22*(D_0 + D_3);
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=D_0*w23 + D_3*w24 + tmp5;
                                    EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+=tmp0 + tmp1;
                                    EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=tmp8 + tmp9;
                                    EM_S[INDEX4(k,m,0,3,numEq,numComp,4)]+=tmp2;
                                    EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+=tmp0 + tmp1;
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=D_1*w23 + D_2*w24 + tmp10;
                                    EM_S[INDEX4(k,m,1,2,numEq,numComp,4)]+=tmp2;
                                    EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+=tmp6 + tmp7;
                                    EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=tmp8 + tmp9;
                                    EM_S[INDEX4(k,m,2,1,numEq,numComp,4)]+=tmp2;
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=D_1*w24 + D_2*w23 + tmp10;
                                    EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+=tmp3 + tmp4;
                                    EM_S[INDEX4(k,m,3,0,numEq,numComp,4)]+=tmp2;
                                    EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+=tmp6 + tmp7;
                                    EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+=tmp3 + tmp4;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=D_0*w24 + D_3*w23 + tmp5;
                                }
                             }
                        } else { // constant data
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const double D_0 = D_p[INDEX2(k, m, numEq)];
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=16*D_0*w22;
                                    EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+=8*D_0*w22;
                                    EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=8*D_0*w22;
                                    EM_S[INDEX4(k,m,0,3,numEq,numComp,4)]+=4*D_0*w22;
                                    EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+=8*D_0*w22;
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=16*D_0*w22;
                                    EM_S[INDEX4(k,m,1,2,numEq,numComp,4)]+=4*D_0*w22;
                                    EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+=8*D_0*w22;
                                    EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=8*D_0*w22;
                                    EM_S[INDEX4(k,m,2,1,numEq,numComp,4)]+=4*D_0*w22;
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=16*D_0*w22;
                                    EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+=8*D_0*w22;
                                    EM_S[INDEX4(k,m,3,0,numEq,numComp,4)]+=4*D_0*w22;
                                    EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+=8*D_0*w22;
                                    EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+=8*D_0*w22;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=16*D_0*w22;
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
                                const double X_0_0 = X_p[INDEX3(k,0,0,numEq,2)];
                                const double X_1_0 = X_p[INDEX3(k,1,0,numEq,2)];
                                const double X_0_1 = X_p[INDEX3(k,0,1,numEq,2)];
                                const double X_1_1 = X_p[INDEX3(k,1,1,numEq,2)];
                                const double X_0_2 = X_p[INDEX3(k,0,2,numEq,2)];
                                const double X_1_2 = X_p[INDEX3(k,1,2,numEq,2)];
                                const double X_0_3 = X_p[INDEX3(k,0,3,numEq,2)];
                                const double X_1_3 = X_p[INDEX3(k,1,3,numEq,2)];
                                const double tmp0 = 6*w15*(X_0_2 + X_0_3);
                                const double tmp1 = 6*w10*(X_0_0 + X_0_1);
                                const double tmp2 = 6*w11*(X_1_0 + X_1_2);
                                const double tmp3 = 6*w14*(X_1_1 + X_1_3);
                                const double tmp4 = 6*w11*(X_1_1 + X_1_3);
                                const double tmp5 = w25*(X_0_0 + X_0_1);
                                const double tmp6 = w26*(X_0_2 + X_0_3);
                                const double tmp7 = 6*w14*(X_1_0 + X_1_2);
                                const double tmp8 = w27*(X_1_0 + X_1_2);
                                const double tmp9 = w28*(X_1_1 + X_1_3);
                                const double tmp10 = w25*(-X_0_2 - X_0_3);
                                const double tmp11 = w26*(-X_0_0 - X_0_1);
                                const double tmp12 = w27*(X_1_1 + X_1_3);
                                const double tmp13 = w28*(X_1_0 + X_1_2);
                                const double tmp14 = w25*(X_0_2 + X_0_3);
                                const double tmp15 = w26*(X_0_0 + X_0_1);
                                EM_F[INDEX2(k,0,numEq)]+=tmp0 + tmp1 + tmp2 + tmp3;
                                EM_F[INDEX2(k,1,numEq)]+=tmp4 + tmp5 + tmp6 + tmp7;
                                EM_F[INDEX2(k,2,numEq)]+=tmp10 + tmp11 + tmp8 + tmp9;
                                EM_F[INDEX2(k,3,numEq)]+=tmp12 + tmp13 + tmp14 + tmp15;
                            }
                        } else { // constant data
                            for (index_t k=0; k<numEq; k++) {
                                const double wX0 = X_p[INDEX2(k, 0, numEq)]*w18;
                                const double wX1 = X_p[INDEX2(k, 1, numEq)]*w19;
                                EM_F[INDEX2(k,0,numEq)]+= 6*wX0 + 6*wX1;
                                EM_F[INDEX2(k,1,numEq)]+=-6*wX0 + 6*wX1;
                                EM_F[INDEX2(k,2,numEq)]+= 6*wX0 - 6*wX1;
                                EM_F[INDEX2(k,3,numEq)]+=-6*wX0 - 6*wX1;
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
                                const double tmp0 = 6*w22*(Y_1 + Y_2);
                                const double tmp1 = 6*w22*(Y_0 + Y_3);
                                EM_F[INDEX2(k,0,numEq)]+=6*Y_0*w20 + 6*Y_3*w21 + tmp0;
                                EM_F[INDEX2(k,1,numEq)]+=6*Y_1*w20 + 6*Y_2*w21 + tmp1;
                                EM_F[INDEX2(k,2,numEq)]+=6*Y_1*w21 + 6*Y_2*w20 + tmp1;
                                EM_F[INDEX2(k,3,numEq)]+=6*Y_0*w21 + 6*Y_3*w20 + tmp0;
                            }
                        } else { // constant data
                            for (index_t k=0; k<numEq; k++) {
                                EM_F[INDEX2(k,0,numEq)]+=36*Y_p[k]*w22;
                                EM_F[INDEX2(k,1,numEq)]+=36*Y_p[k]*w22;
                                EM_F[INDEX2(k,2,numEq)]+=36*Y_p[k]*w22;
                                EM_F[INDEX2(k,3,numEq)]+=36*Y_p[k]*w22;
                            }
                        }
                    }

                    // add to matrix (if addEM_S) and RHS (if addEM_F)
                    const index_t firstNode=m_NN[0]*k1+k0;
                    domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, addEM_S,
                                          addEM_F, firstNode, numEq, numComp);
                } // end k0 loop
            } // end k1 loop
        } // end of colouring
    } // end of parallel region
}

} // namespace ripley

