
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
**
*****************************************************************************/

#include <ripley/WaveAssembler2D.h>
#include <ripley/domainhelpers.h>

#include <escript/index.h>

using escript::Data;

namespace ripley {

WaveAssembler2D::WaveAssembler2D(escript::const_Domain_ptr dom,
                                 const double *dx, const dim_t *NE,
                                 const dim_t *NN, const DataMap& c)
    : AbstractAssembler(),
    m_dx(dx),
    m_NE(NE),
    m_NN(NN)
{
    domain = REFCOUNTNS::static_pointer_cast<const Rectangle>(dom);
    isHTI = isVTI = false;
    DataMap::const_iterator a = c.find("c12"), b = c.find("c23");
    if (c.find("c11") == c.end()
                || c.find("c13") == c.end() || c.find("c33") == c.end()
                || c.find("c44") == c.end() || c.find("c66") == c.end()
                || (a == c.end() && b == c.end()))
        throw escript::ValueError("required constants missing for WaveAssembler");

    if (a != c.end() && b != c.end()) {
        throw escript::NotImplementedError("WaveAssembler2D() doesn't support general "
                              "form waves (yet)");
    } else if (a == c.end()) {
        c23 = b->second;
        isHTI = true;
    } else if (b == c.end()) {
        c12 = a->second;
        isVTI = true;
    } // final else case taken care of with the missing constants above
    c11 = c.find("c11")->second;
    c13 = c.find("c13")->second;
    c33 = c.find("c33")->second;
    c44 = c.find("c44")->second;
    c66 = c.find("c66")->second;

    int fs = c11.getFunctionSpace().getTypeCode();

    if (fs != c13.getFunctionSpace().getTypeCode()
            || fs != c33.getFunctionSpace().getTypeCode()
            || fs != c44.getFunctionSpace().getTypeCode()
            || fs != c66.getFunctionSpace().getTypeCode()) {
        throw escript::ValueError("C tensor elements are in mismatching function spaces");
    }
}

void WaveAssembler2D::collateFunctionSpaceTypes(std::vector<int>& fsTypes,
                                                const DataMap& coefs) const
{
    if (isNotEmpty("A", coefs))
        fsTypes.push_back(coefs.find("A")->second.getFunctionSpace().getTypeCode());
    if (isNotEmpty("B", coefs))
        fsTypes.push_back(coefs.find("B")->second.getFunctionSpace().getTypeCode());
    if (isNotEmpty("C", coefs))
        fsTypes.push_back(coefs.find("C")->second.getFunctionSpace().getTypeCode());
    if (isNotEmpty("D", coefs))
        fsTypes.push_back(coefs.find("D")->second.getFunctionSpace().getTypeCode());
    if (isNotEmpty("du", coefs))
        fsTypes.push_back(coefs.find("du")->second.getFunctionSpace().getTypeCode());
    if (isNotEmpty("Y", coefs))
        fsTypes.push_back(coefs.find("Y")->second.getFunctionSpace().getTypeCode());
}

void WaveAssembler2D::assemblePDESystem(escript::AbstractSystemMatrix* mat,
                                        Data& rhs, const DataMap& coefs) const
{
    if (isNotEmpty("X", coefs))
        throw escript::ValueError("Coefficient X was given to WaveAssembler "
                "unexpectedly. Specialised domains can't be used for general "
                "assemblage.");
    const Data& A = unpackData("A", coefs);
    const Data& B = unpackData("B", coefs);
    const Data& C = unpackData("C", coefs);
    const Data& D = unpackData("D", coefs);
    const Data& Y = unpackData("Y", coefs);
    const Data& du = unpackData("du", coefs);

    if ((!du.isEmpty()) && du.getFunctionSpace().getTypeCode() != c11.getFunctionSpace().getTypeCode()) {
        throw escript::ValueError("WaveAssembler2D: du and C tensor in mismatching function spaces");
    }

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
    const int NE0 = m_NE[0];
    const int NE1 = m_NE[1];
    const bool addEM_S = (!A.isEmpty() || !B.isEmpty() || !C.isEmpty() || !D.isEmpty());
    const bool addEM_F = (!du.isEmpty() || !Y.isEmpty());
    rhs.requireWrite();

#pragma omp parallel
    {
        std::vector<double> EM_S(4*4*numEq*numComp, 0);
        std::vector<double> EM_F(4*numEq, 0);

        for (index_t k1_0=0; k1_0<2; k1_0++) { // colouring
#pragma omp for
            for (index_t k1=k1_0; k1 < NE1; k1+=2) {
                for (index_t k0=0; k0 < NE0; ++k0)  {
                    const index_t e = k0 + NE0*k1;
                    if (addEM_S)
                        fill(EM_S.begin(), EM_S.end(), 0);
                    if (addEM_F)
                        fill(EM_F.begin(), EM_F.end(), 0);

                    ///////////////
                    // process A //
                    ///////////////
                    if (!A.isEmpty()) {
                        const double* A_p = A.getSampleDataRO(e);
                        if (A.actsExpanded()) {
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const double A_00_0 = A_p[INDEX5(k,0,m,0,0,numEq,2,numComp,2)];
                                    const double A_01_0 = A_p[INDEX5(k,0,m,1,0,numEq,2,numComp,2)];
                                    const double A_10_0 = A_p[INDEX5(k,1,m,0,0,numEq,2,numComp,2)];
                                    const double A_11_0 = A_p[INDEX5(k,1,m,1,0,numEq,2,numComp,2)];
                                    const double A_00_1 = A_p[INDEX5(k,0,m,0,1,numEq,2,numComp,2)];
                                    const double A_01_1 = A_p[INDEX5(k,0,m,1,1,numEq,2,numComp,2)];
                                    const double A_10_1 = A_p[INDEX5(k,1,m,0,1,numEq,2,numComp,2)];
                                    const double A_11_1 = A_p[INDEX5(k,1,m,1,1,numEq,2,numComp,2)];
                                    const double A_00_2 = A_p[INDEX5(k,0,m,0,2,numEq,2,numComp,2)];
                                    const double A_01_2 = A_p[INDEX5(k,0,m,1,2,numEq,2,numComp,2)];
                                    const double A_10_2 = A_p[INDEX5(k,1,m,0,2,numEq,2,numComp,2)];
                                    const double A_11_2 = A_p[INDEX5(k,1,m,1,2,numEq,2,numComp,2)];
                                    const double A_00_3 = A_p[INDEX5(k,0,m,0,3,numEq,2,numComp,2)];
                                    const double A_01_3 = A_p[INDEX5(k,0,m,1,3,numEq,2,numComp,2)];
                                    const double A_10_3 = A_p[INDEX5(k,1,m,0,3,numEq,2,numComp,2)];
                                    const double A_11_3 = A_p[INDEX5(k,1,m,1,3,numEq,2,numComp,2)];
                                    const double tmp0 = w3*(A_11_0 + A_11_1 + A_11_2 + A_11_3);
                                    const double tmp1 = w1*(A_01_0 + A_01_3 - A_10_1 - A_10_2);
                                    const double tmp2 = w4*(A_00_2 + A_00_3);
                                    const double tmp3 = w0*(A_00_0 + A_00_1);
                                    const double tmp4 = w5*(A_01_2 - A_10_3);
                                    const double tmp5 = w2*(-A_01_1 + A_10_0);
                                    const double tmp6 = w5*(A_01_3 + A_10_0);
                                    const double tmp7 = w3*(-A_11_0 - A_11_1 - A_11_2 - A_11_3);
                                    const double tmp8 = w6*(A_00_0 + A_00_1 + A_00_2 + A_00_3);
                                    const double tmp9 = w1*(A_01_1 + A_01_2 + A_10_1 + A_10_2);
                                    const double tmp10 = w2*(-A_01_0 - A_10_3);
                                    const double tmp11 = w4*(A_00_0 + A_00_1);
                                    const double tmp12 = w0*(A_00_2 + A_00_3);
                                    const double tmp13 = w5*(A_01_1 - A_10_0);
                                    const double tmp14 = w2*(-A_01_2 + A_10_3);
                                    const double tmp15 = w7*(A_11_0 + A_11_2);
                                    const double tmp16 = w4*(-A_00_2 - A_00_3);
                                    const double tmp17 = w0*(-A_00_0 - A_00_1);
                                    const double tmp18 = w5*(A_01_3 + A_10_3);
                                    const double tmp19 = w8*(A_11_1 + A_11_3);
                                    const double tmp20 = w2*(-A_01_0 - A_10_0);
                                    const double tmp21 = w7*(A_11_1 + A_11_3);
                                    const double tmp22 = w4*(-A_00_0 - A_00_1);
                                    const double tmp23 = w0*(-A_00_2 - A_00_3);
                                    const double tmp24 = w5*(A_01_0 + A_10_0);
                                    const double tmp25 = w8*(A_11_0 + A_11_2);
                                    const double tmp26 = w2*(-A_01_3 - A_10_3);
                                    const double tmp27 = w5*(-A_01_1 - A_10_2);
                                    const double tmp28 = w1*(-A_01_0 - A_01_3 - A_10_0 - A_10_3);
                                    const double tmp29 = w2*(A_01_2 + A_10_1);
                                    const double tmp30 = w7*(-A_11_1 - A_11_3);
                                    const double tmp31 = w1*(-A_01_1 - A_01_2 + A_10_0 + A_10_3);
                                    const double tmp32 = w5*(-A_01_0 + A_10_2);
                                    const double tmp33 = w8*(-A_11_0 - A_11_2);
                                    const double tmp34 = w6*(-A_00_0 - A_00_1 - A_00_2 - A_00_3);
                                    const double tmp35 = w2*(A_01_3 - A_10_1);
                                    const double tmp36 = w5*(A_01_0 + A_10_3);
                                    const double tmp37 = w2*(-A_01_3 - A_10_0);
                                    const double tmp38 = w7*(-A_11_0 - A_11_2);
                                    const double tmp39 = w5*(-A_01_3 + A_10_1);
                                    const double tmp40 = w8*(-A_11_1 - A_11_3);
                                    const double tmp41 = w2*(A_01_0 - A_10_2);
                                    const double tmp42 = w5*(A_01_1 - A_10_3);
                                    const double tmp43 = w2*(-A_01_2 + A_10_0);
                                    const double tmp44 = w5*(A_01_2 - A_10_0);
                                    const double tmp45 = w2*(-A_01_1 + A_10_3);
                                    const double tmp46 = w5*(-A_01_0 + A_10_1);
                                    const double tmp47 = w2*(A_01_3 - A_10_2);
                                    const double tmp48 = w5*(-A_01_1 - A_10_1);
                                    const double tmp49 = w2*(A_01_2 + A_10_2);
                                    const double tmp50 = w5*(-A_01_3 + A_10_2);
                                    const double tmp51 = w2*(A_01_0 - A_10_1);
                                    const double tmp52 = w5*(-A_01_2 - A_10_1);
                                    const double tmp53 = w2*(A_01_1 + A_10_2);
                                    const double tmp54 = w5*(-A_01_2 - A_10_2);
                                    const double tmp55 = w2*(A_01_1 + A_10_1);
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=tmp15 + tmp16 + tmp17 + tmp18 + tmp19 + tmp20 + tmp9;
                                    EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+=tmp0 + tmp1 + tmp2 + tmp3 + tmp4 + tmp5;
                                    EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=tmp31 + tmp34 + tmp38 + tmp39 + tmp40 + tmp41;
                                    EM_S[INDEX4(k,m,0,3,numEq,numComp,4)]+=tmp28 + tmp52 + tmp53 + tmp7 + tmp8;
                                    EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+=tmp0 + tmp2 + tmp3 + tmp31 + tmp50 + tmp51;
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=tmp16 + tmp17 + tmp21 + tmp25 + tmp28 + tmp54 + tmp55;
                                    EM_S[INDEX4(k,m,1,2,numEq,numComp,4)]+=tmp10 + tmp6 + tmp7 + tmp8 + tmp9;
                                    EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+=tmp1 + tmp30 + tmp33 + tmp34 + tmp44 + tmp45;
                                    EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=tmp1 + tmp34 + tmp38 + tmp40 + tmp42 + tmp43;
                                    EM_S[INDEX4(k,m,2,1,numEq,numComp,4)]+=tmp36 + tmp37 + tmp7 + tmp8 + tmp9;
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=tmp15 + tmp19 + tmp22 + tmp23 + tmp28 + tmp48 + tmp49;
                                    EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+=tmp0 + tmp11 + tmp12 + tmp31 + tmp46 + tmp47;
                                    EM_S[INDEX4(k,m,3,0,numEq,numComp,4)]+=tmp27 + tmp28 + tmp29 + tmp7 + tmp8;
                                    EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+=tmp30 + tmp31 + tmp32 + tmp33 + tmp34 + tmp35;
                                    EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+=tmp0 + tmp1 + tmp11 + tmp12 + tmp13 + tmp14;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=tmp21 + tmp22 + tmp23 + tmp24 + tmp25 + tmp26 + tmp9;
                                }
                            }
                        } else { // constant data
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const double A_00 = A_p[INDEX4(k,0,m,0, numEq,2, numComp)];
                                    const double A_01 = A_p[INDEX4(k,0,m,1, numEq,2, numComp)];
                                    const double A_10 = A_p[INDEX4(k,1,m,0, numEq,2, numComp)];
                                    const double A_11 = A_p[INDEX4(k,1,m,1, numEq,2, numComp)];
                                    const double tmp0 = 6*w1*(A_01 - A_10);
                                    const double tmp1 = 6*w1*(A_01 + A_10);
                                    const double tmp2 = 6*w1*(-A_01 - A_10);
                                    const double tmp3 = 6*w1*(-A_01 + A_10);
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=-8*A_00*w6 + 8*A_11*w3 + tmp1;
                                    EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+=8*A_00*w6 + 4*A_11*w3 + tmp0;
                                    EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=-4*A_00*w6 - 8*A_11*w3 + tmp3;
                                    EM_S[INDEX4(k,m,0,3,numEq,numComp,4)]+=4*A_00*w6 - 4*A_11*w3 + tmp2;
                                    EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+=8*A_00*w6 + 4*A_11*w3 + tmp3;
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=-8*A_00*w6 + 8*A_11*w3 + tmp2;
                                    EM_S[INDEX4(k,m,1,2,numEq,numComp,4)]+=4*A_00*w6 - 4*A_11*w3 + tmp1;
                                    EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+=-4*A_00*w6 - 8*A_11*w3 + tmp0;
                                    EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=-4*A_00*w6 - 8*A_11*w3 + tmp0;
                                    EM_S[INDEX4(k,m,2,1,numEq,numComp,4)]+=4*A_00*w6 - 4*A_11*w3 + tmp1;
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=-8*A_00*w6 + 8*A_11*w3 + tmp2;
                                    EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+=8*A_00*w6 + 4*A_11*w3 + tmp3;
                                    EM_S[INDEX4(k,m,3,0,numEq,numComp,4)]+=4*A_00*w6 - 4*A_11*w3 + tmp2;
                                    EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+=-4*A_00*w6 - 8*A_11*w3 + tmp3;
                                    EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+=8*A_00*w6 + 4*A_11*w3 + tmp0;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=-8*A_00*w6 + 8*A_11*w3 + tmp1;
                                }
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
                    ////////////////
                    // process du //
                    ////////////////
                    if (!du.isEmpty()) {
                        const double *du_p = du.getSampleDataRO(e);
                        const double *c11_p = c11.getSampleDataRO(e);
                        const double *c13_p = c13.getSampleDataRO(e);
                        const double *c33_p = c33.getSampleDataRO(e);
                        if (du.actsExpanded()) {
                            double X_00_0, X_00_1, X_00_2, X_00_3;
                            double X_11_0, X_11_1, X_11_2, X_11_3;
                            double X_01_0 = -(du_p[INDEX3(1,0,0,numEq,2)] + du_p[INDEX3(0,1,0,numEq,2)]);
                            double X_01_1 = -(du_p[INDEX3(1,0,1,numEq,2)] + du_p[INDEX3(0,1,1,numEq,2)]);
                            double X_01_2 = -(du_p[INDEX3(1,0,2,numEq,2)] + du_p[INDEX3(0,1,2,numEq,2)]);
                            double X_01_3 = -(du_p[INDEX3(1,0,3,numEq,2)] + du_p[INDEX3(0,1,3,numEq,2)]);

                            if (isVTI) {
                                const double *c44_p = c44.getSampleDataRO(e);
                                X_00_0 = -(du_p[INDEX3(0,0,0,numEq,2)] * c11_p[0]
                                        + du_p[INDEX3(1,1,0,numEq,2)] * c13_p[0]);
                                X_00_1 = -(du_p[INDEX3(0,0,1,numEq,2)] * c11_p[1]
                                        + du_p[INDEX3(1,1,1,numEq,2)] * c13_p[1]);
                                X_00_2 = -(du_p[INDEX3(0,0,2,numEq,2)] * c11_p[2]
                                        + du_p[INDEX3(1,1,2,numEq,2)] * c13_p[2]);
                                X_00_3 = -(du_p[INDEX3(0,0,3,numEq,2)] * c11_p[3]
                                        + du_p[INDEX3(1,1,3,numEq,2)] * c13_p[3]);
                                X_11_0 = -(du_p[INDEX3(0,0,0,numEq,2)] * c13_p[0]
                                        + du_p[INDEX3(1,1,0,numEq,2)] * c33_p[0]);
                                X_11_1 = -(du_p[INDEX3(0,0,1,numEq,2)] * c13_p[1]
                                        + du_p[INDEX3(1,1,1,numEq,2)] * c33_p[1]);
                                X_11_2 = -(du_p[INDEX3(0,0,2,numEq,2)] * c13_p[2]
                                        + du_p[INDEX3(1,1,2,numEq,2)] * c33_p[2]);
                                X_11_3 = -(du_p[INDEX3(0,0,3,numEq,2)] * c13_p[3]
                                        + du_p[INDEX3(1,1,3,numEq,2)] * c33_p[3]);
                                X_01_0 *= c44_p[0];
                                X_01_1 *= c44_p[1];
                                X_01_2 *= c44_p[2];
                                X_01_3 *= c44_p[3];
                            } else { // isHTI
                                const double *c66_p = c66.getSampleDataRO(e);
                                X_00_0 = -(du_p[INDEX3(0,0,0,numEq,2)] * c11_p[0]
                                        + du_p[INDEX3(1,1,0,numEq,2)] * c13_p[0]);
                                X_00_1 = -(du_p[INDEX3(0,0,1,numEq,2)] * c11_p[1]
                                        + du_p[INDEX3(1,1,1,numEq,2)] * c13_p[1]);
                                X_00_2 = -(du_p[INDEX3(0,0,2,numEq,2)] * c11_p[2]
                                        + du_p[INDEX3(1,1,2,numEq,2)] * c13_p[2]);
                                X_00_3 = -(du_p[INDEX3(0,0,3,numEq,2)] * c11_p[3]
                                        + du_p[INDEX3(1,1,3,numEq,2)] * c13_p[3]);
                                X_11_0 = -(du_p[INDEX3(0,0,0,numEq,2)] * c13_p[0]
                                        + du_p[INDEX3(1,1,0,numEq,2)] * c33_p[0]);
                                X_11_1 = -(du_p[INDEX3(0,0,1,numEq,2)] * c13_p[1]
                                        + du_p[INDEX3(1,1,1,numEq,2)] * c33_p[1]);
                                X_11_2 = -(du_p[INDEX3(0,0,2,numEq,2)] * c13_p[2]
                                        + du_p[INDEX3(1,1,2,numEq,2)] * c33_p[2]);
                                X_11_3 = -(du_p[INDEX3(0,0,3,numEq,2)] * c13_p[3]
                                        + du_p[INDEX3(1,1,3,numEq,2)] * c33_p[3]);
                                X_01_0 *= c66_p[0];
                                X_01_1 *= c66_p[1];
                                X_01_2 *= c66_p[2];
                                X_01_3 *= c66_p[3];
                            }
                            const double X_10_0 = X_01_0;
                            const double X_10_1 = X_01_1;
                            const double X_10_2 = X_01_2;
                            const double X_10_3 = X_01_3;

                            const double Atmp0 = 6*w15*(X_00_2 + X_00_3);
                            const double Atmp1 = 6*w10*(X_00_0 + X_00_1);
                            const double Atmp2 = 6*w11*(X_01_0 + X_01_2);
                            const double Atmp3 = 6*w14*(X_01_1 + X_01_3);
                            const double Atmp4 = 6*w11*(X_01_1 + X_01_3);
                            const double Atmp5 = w25*(X_00_0 + X_00_1);
                            const double Atmp6 = w26*(X_00_2 + X_00_3);
                            const double Atmp7 = 6*w14*(X_01_0 + X_01_2);
                            const double Atmp8 = w27*(X_01_0 + X_01_2);
                            const double Atmp9 = w28*(X_01_1 + X_01_3);
                            const double Atmp10 = w25*(-X_00_2 - X_00_3);
                            const double Atmp11 = w26*(-X_00_0 - X_00_1);
                            const double Atmp12 = w27*(X_01_1 + X_01_3);
                            const double Atmp13 = w28*(X_01_0 + X_01_2);
                            const double Atmp14 = w25*(X_00_2 + X_00_3);
                            const double Atmp15 = w26*(X_00_0 + X_00_1);
                            const double Btmp0 = 6*w15*(X_10_2 + X_10_3);
                            const double Btmp1 = 6*w10*(X_10_0 + X_10_1);
                            const double Btmp2 = 6*w11*(X_11_0 + X_11_2);
                            const double Btmp3 = 6*w14*(X_11_1 + X_11_3);
                            const double Btmp4 = 6*w11*(X_11_1 + X_11_3);
                            const double Btmp5 = w25*(X_10_0 + X_10_1);
                            const double Btmp6 = w26*(X_10_2 + X_10_3);
                            const double Btmp7 = 6*w14*(X_11_0 + X_11_2);
                            const double Btmp8 = w27*(X_11_0 + X_11_2);
                            const double Btmp9 = w28*(X_11_1 + X_11_3);
                            const double Btmp10 = w25*(-X_10_2 - X_10_3);
                            const double Btmp11 = w26*(-X_10_0 - X_10_1);
                            const double Btmp12 = w27*(X_11_1 + X_11_3);
                            const double Btmp13 = w28*(X_11_0 + X_11_2);
                            const double Btmp14 = w25*(X_10_2 + X_10_3);
                            const double Btmp15 = w26*(X_10_0 + X_10_1);
                            EM_F[INDEX2(0,0,numEq)]+=Atmp0 + Atmp1 + Atmp2 + Atmp3;
                            EM_F[INDEX2(1,0,numEq)]+=Btmp0 + Btmp1 + Btmp2 + Btmp3;
                            EM_F[INDEX2(0,1,numEq)]+=Atmp4 + Atmp5 + Atmp6 + Atmp7;
                            EM_F[INDEX2(1,1,numEq)]+=Btmp4 + Btmp5 + Btmp6 + Btmp7;
                            EM_F[INDEX2(0,2,numEq)]+=Atmp10 + Atmp11 + Atmp8 + Atmp9;
                            EM_F[INDEX2(1,2,numEq)]+=Btmp10 + Btmp11 + Btmp8 + Btmp9;
                            EM_F[INDEX2(0,3,numEq)]+=Atmp12 + Atmp13 + Atmp14 + Atmp15;
                            EM_F[INDEX2(1,3,numEq)]+=Btmp12 + Btmp13 + Btmp14 + Btmp15;

                        } else { // constant data
                            double wX_00, wX_01, wX_10, wX_11;
                            if (isVTI) {
                                const double *c44_p = c44.getSampleDataRO(e);
                                wX_00 = -(du_p[INDEX2(0,0,numEq)] * c11_p[0]
                                                    + du_p[INDEX2(1,1,numEq)] * c13_p[0])*w18;
                                wX_01 = -(c44_p[0] *
                                                    (du_p[INDEX2(1,0,numEq)] + du_p[INDEX2(0,1,numEq)]))*w19;
                                wX_10 = -(c44_p[0] *
                                                    (du_p[INDEX2(1,0,numEq)] + du_p[INDEX2(0,1,numEq)]))*w18;
                                wX_11 = -(du_p[INDEX2(0,0,numEq)] * c13_p[0]
                                                    + du_p[INDEX2(1,1,numEq)] * c33_p[0])*w19;
                            } else { // isHTI
                                const double *c66_p = c66.getSampleDataRO(e);
                                wX_00 = -(du_p[INDEX2(0,0,numEq)] * c11_p[0]
                                        + du_p[INDEX2(1,1,numEq)] * c13_p[0])*w18;
                                wX_01 = -(c66_p[0] *
                                                    (du_p[INDEX2(1,0,numEq)] + du_p[INDEX2(0,1,numEq)]))*w19;
                                wX_10 = -(c66_p[0] *
                                                    (du_p[INDEX2(1,0,numEq)] + du_p[INDEX2(0,1,numEq)]))*w18;
                                wX_11 = -(du_p[INDEX2(0,0,numEq)] * c13_p[0]
                                        + du_p[INDEX2(1,1,numEq)] * c33_p[0])*w19;
                            }
                            EM_F[INDEX2(0,0,numEq)]+= wX_00 + wX_01;
                            EM_F[INDEX2(1,0,numEq)]+= wX_10 + wX_11;
                            EM_F[INDEX2(0,1,numEq)]+=-wX_00 + wX_01;
                            EM_F[INDEX2(1,1,numEq)]+=-wX_10 + wX_11;
                            EM_F[INDEX2(0,2,numEq)]+= wX_00 - wX_01;
                            EM_F[INDEX2(1,2,numEq)]+= wX_10 - wX_11;
                            EM_F[INDEX2(0,3,numEq)]+=-wX_00 - wX_01;
                            EM_F[INDEX2(1,3,numEq)]+=-wX_10 - wX_11;
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

