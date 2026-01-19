
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

#include <ripley/DefaultAssembler2D.h>
#include <ripley/domainhelpers.h>

#include <escript/DataTypes.h>
#include <escript/index.h>

using namespace std;
using escript::AbstractSystemMatrix;
using escript::Data;

namespace ripley {

template<class Scalar>
void DefaultAssembler2D<Scalar>::collateFunctionSpaceTypes(
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
void DefaultAssembler2D<Scalar>::assemblePDESingle(AbstractSystemMatrix* mat,
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
void DefaultAssembler2D<Scalar>::assemblePDEBoundarySingle(
                                        AbstractSystemMatrix* mat,
                                        Data& rhs, const DataMap& coefs) const 
{
    const Data& d = unpackData("d", coefs);
    const Data& y = unpackData("y", coefs);
    assemblePDEBoundarySingle(mat, rhs, d, y);
}

template<class Scalar>
void DefaultAssembler2D<Scalar>::assemblePDESingleReduced(
                                        AbstractSystemMatrix* mat,
                                        Data& rhs, const DataMap& coefs) const
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
void DefaultAssembler2D<Scalar>::assemblePDEBoundarySingleReduced(
                                        AbstractSystemMatrix* mat, Data& rhs,
                                        const DataMap& coefs) const
{
    const Data& d = unpackData("d", coefs);
    const Data& y = unpackData("y", coefs);
    assemblePDEBoundarySingleReduced(mat, rhs, d, y);
}

template<class Scalar>
void DefaultAssembler2D<Scalar>::assemblePDESystem(AbstractSystemMatrix* mat,
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
void DefaultAssembler2D<Scalar>::assemblePDEBoundarySystem(
                                        AbstractSystemMatrix* mat,
                                        Data& rhs, const DataMap& coefs) const
{
    const Data& d = unpackData("d", coefs);
    const Data& y = unpackData("y", coefs);
    assemblePDEBoundarySystem(mat, rhs, d, y);
}

template<class Scalar>
void DefaultAssembler2D<Scalar>::assemblePDESystemReduced(
                                        AbstractSystemMatrix* mat,
                                        Data& rhs, const DataMap& coefs) const
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
void DefaultAssembler2D<Scalar>::assemblePDEBoundarySystemReduced(
                                        AbstractSystemMatrix* mat, Data& rhs,
                                        const DataMap& coefs) const
{
    const Data& d = unpackData("d", coefs);
    const Data& y = unpackData("y", coefs);
    assemblePDEBoundarySystemReduced(mat, rhs, d, y);
}

/****************************************************************************/
// PDE SINGLE
/****************************************************************************/

template<class Scalar>
void DefaultAssembler2D<Scalar>::assemblePDESingle(AbstractSystemMatrix* mat,
                                      Data& rhs, const Data& A, const Data& B,
                                      const Data& C, const Data& D,
                                      const Data& X, const Data& Y) const
{
    const double SQRT3 = 1.73205080756887719318;
    const double w1 = 1.0/24.0;
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
    const double w12 = w18*(5*SQRT3 + 9)/12;
    const double w13 = w18*(-5*SQRT3 + 9)/12;
    const double w10 = w18*(SQRT3 + 3)/12;
    const double w15 = w18*(-SQRT3 + 3)/12;
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
    const bool addEM_F = (!X.isEmpty() || !Y.isEmpty());
    const Scalar zero = static_cast<Scalar>(0);
    rhs.requireWrite();

#pragma omp parallel
    {
        vector<Scalar> EM_S(4*4, zero);
        vector<Scalar> EM_F(4, zero);

        for (index_t k1_0 = 0; k1_0 < 2; k1_0++) { // colouring
#pragma omp for
            for (index_t k1 = k1_0; k1 < NE1; k1+=2) {
                for (index_t k0 = 0; k0 < NE0; ++k0)  {
                    const index_t e = k0 + NE0*k1;
                    if (addEM_S)
                        fill(EM_S.begin(), EM_S.end(), zero);
                    if (addEM_F)
                        fill(EM_F.begin(), EM_F.end(), zero);
                    ///////////////
                    // process A //
                    ///////////////
                    if (!A.isEmpty()) {
                        const Scalar* A_p = A.getSampleDataRO(e, zero);
                        if (A.actsExpanded()) {
                            const Scalar A_00_0 = A_p[INDEX3(0,0,0,2,2)];
                            const Scalar A_01_0 = A_p[INDEX3(0,1,0,2,2)];
                            const Scalar A_10_0 = A_p[INDEX3(1,0,0,2,2)];
                            const Scalar A_11_0 = A_p[INDEX3(1,1,0,2,2)];
                            const Scalar A_00_1 = A_p[INDEX3(0,0,1,2,2)];
                            const Scalar A_01_1 = A_p[INDEX3(0,1,1,2,2)];
                            const Scalar A_10_1 = A_p[INDEX3(1,0,1,2,2)];
                            const Scalar A_11_1 = A_p[INDEX3(1,1,1,2,2)];
                            const Scalar A_00_2 = A_p[INDEX3(0,0,2,2,2)];
                            const Scalar A_01_2 = A_p[INDEX3(0,1,2,2,2)];
                            const Scalar A_10_2 = A_p[INDEX3(1,0,2,2,2)];
                            const Scalar A_11_2 = A_p[INDEX3(1,1,2,2,2)];
                            const Scalar A_00_3 = A_p[INDEX3(0,0,3,2,2)];
                            const Scalar A_01_3 = A_p[INDEX3(0,1,3,2,2)];
                            const Scalar A_10_3 = A_p[INDEX3(1,0,3,2,2)];
                            const Scalar A_11_3 = A_p[INDEX3(1,1,3,2,2)];
                            const Scalar tmp0 = w3*(A_11_0 + A_11_1 + A_11_2 + A_11_3);
                            const Scalar tmp1 = w1*(A_01_0 + A_01_3 - A_10_1 - A_10_2);
                            const Scalar tmp2 = w4*(A_00_2 + A_00_3);
                            const Scalar tmp3 = w0*(A_00_0 + A_00_1);
                            const Scalar tmp4 = w5*(A_01_2 - A_10_3);
                            const Scalar tmp5 = w2*(-A_01_1 + A_10_0);
                            const Scalar tmp6 = w5*(A_01_3 + A_10_0);
                            const Scalar tmp7 = w3*(-A_11_0 - A_11_1 - A_11_2 - A_11_3);
                            const Scalar tmp8 = w6*(A_00_0 + A_00_1 + A_00_2 + A_00_3);
                            const Scalar tmp9 = w1*(A_01_1 + A_01_2 + A_10_1 + A_10_2);
                            const Scalar tmp10 = w2*(-A_01_0 - A_10_3);
                            const Scalar tmp11 = w4*(A_00_0 + A_00_1);
                            const Scalar tmp12 = w0*(A_00_2 + A_00_3);
                            const Scalar tmp13 = w5*(A_01_1 - A_10_0);
                            const Scalar tmp14 = w2*(-A_01_2 + A_10_3);
                            const Scalar tmp15 = w7*(A_11_0 + A_11_2);
                            const Scalar tmp16 = w4*(-A_00_2 - A_00_3);
                            const Scalar tmp17 = w0*(-A_00_0 - A_00_1);
                            const Scalar tmp18 = w5*(A_01_3 + A_10_3);
                            const Scalar tmp19 = w8*(A_11_1 + A_11_3);
                            const Scalar tmp20 = w2*(-A_01_0 - A_10_0);
                            const Scalar tmp21 = w7*(A_11_1 + A_11_3);
                            const Scalar tmp22 = w4*(-A_00_0 - A_00_1);
                            const Scalar tmp23 = w0*(-A_00_2 - A_00_3);
                            const Scalar tmp24 = w5*(A_01_0 + A_10_0);
                            const Scalar tmp25 = w8*(A_11_0 + A_11_2);
                            const Scalar tmp26 = w2*(-A_01_3 - A_10_3);
                            const Scalar tmp27 = w5*(-A_01_1 - A_10_2);
                            const Scalar tmp28 = w1*(-A_01_0 - A_01_3 - A_10_0 - A_10_3);
                            const Scalar tmp29 = w2*(A_01_2 + A_10_1);
                            const Scalar tmp30 = w7*(-A_11_1 - A_11_3);
                            const Scalar tmp31 = w1*(-A_01_1 - A_01_2 + A_10_0 + A_10_3);
                            const Scalar tmp32 = w5*(-A_01_0 + A_10_2);
                            const Scalar tmp33 = w8*(-A_11_0 - A_11_2);
                            const Scalar tmp34 = w6*(-A_00_0 - A_00_1 - A_00_2 - A_00_3);
                            const Scalar tmp35 = w2*(A_01_3 - A_10_1);
                            const Scalar tmp36 = w5*(A_01_0 + A_10_3);
                            const Scalar tmp37 = w2*(-A_01_3 - A_10_0);
                            const Scalar tmp38 = w7*(-A_11_0 - A_11_2);
                            const Scalar tmp39 = w5*(-A_01_3 + A_10_1);
                            const Scalar tmp40 = w8*(-A_11_1 - A_11_3);
                            const Scalar tmp41 = w2*(A_01_0 - A_10_2);
                            const Scalar tmp42 = w5*(A_01_1 - A_10_3);
                            const Scalar tmp43 = w2*(-A_01_2 + A_10_0);
                            const Scalar tmp44 = w5*(A_01_2 - A_10_0);
                            const Scalar tmp45 = w2*(-A_01_1 + A_10_3);
                            const Scalar tmp46 = w5*(-A_01_0 + A_10_1);
                            const Scalar tmp47 = w2*(A_01_3 - A_10_2);
                            const Scalar tmp48 = w5*(-A_01_1 - A_10_1);
                            const Scalar tmp49 = w2*(A_01_2 + A_10_2);
                            const Scalar tmp50 = w5*(-A_01_3 + A_10_2);
                            const Scalar tmp51 = w2*(A_01_0 - A_10_1);
                            const Scalar tmp52 = w5*(-A_01_2 - A_10_1);
                            const Scalar tmp53 = w2*(A_01_1 + A_10_2);
                            const Scalar tmp54 = w5*(-A_01_2 - A_10_2);
                            const Scalar tmp55 = w2*(A_01_1 + A_10_1);
                            EM_S[INDEX2(0,0,4)]+=tmp15 + tmp16 + tmp17 + tmp18 + tmp19 + tmp20 + tmp9;
                            EM_S[INDEX2(0,1,4)]+=tmp0 + tmp1 + tmp2 + tmp3 + tmp4 + tmp5;
                            EM_S[INDEX2(0,2,4)]+=tmp31 + tmp34 + tmp38 + tmp39 + tmp40 + tmp41;
                            EM_S[INDEX2(0,3,4)]+=tmp28 + tmp52 + tmp53 + tmp7 + tmp8;
                            EM_S[INDEX2(1,0,4)]+=tmp0 + tmp2 + tmp3 + tmp31 + tmp50 + tmp51;
                            EM_S[INDEX2(1,1,4)]+=tmp16 + tmp17 + tmp21 + tmp25 + tmp28 + tmp54 + tmp55;
                            EM_S[INDEX2(1,2,4)]+=tmp10 + tmp6 + tmp7 + tmp8 + tmp9;
                            EM_S[INDEX2(1,3,4)]+=tmp1 + tmp30 + tmp33 + tmp34 + tmp44 + tmp45;
                            EM_S[INDEX2(2,0,4)]+=tmp1 + tmp34 + tmp38 + tmp40 + tmp42 + tmp43;
                            EM_S[INDEX2(2,1,4)]+=tmp36 + tmp37 + tmp7 + tmp8 + tmp9;
                            EM_S[INDEX2(2,2,4)]+=tmp15 + tmp19 + tmp22 + tmp23 + tmp28 + tmp48 + tmp49;
                            EM_S[INDEX2(2,3,4)]+=tmp0 + tmp11 + tmp12 + tmp31 + tmp46 + tmp47;
                            EM_S[INDEX2(3,0,4)]+=tmp27 + tmp28 + tmp29 + tmp7 + tmp8;
                            EM_S[INDEX2(3,1,4)]+=tmp30 + tmp31 + tmp32 + tmp33 + tmp34 + tmp35;
                            EM_S[INDEX2(3,2,4)]+=tmp0 + tmp1 + tmp11 + tmp12 + tmp13 + tmp14;
                            EM_S[INDEX2(3,3,4)]+=tmp21 + tmp22 + tmp23 + tmp24 + tmp25 + tmp26 + tmp9;
                        } else { // constant data
                            const Scalar A_00 = A_p[INDEX2(0,0,2)];
                            const Scalar A_01 = A_p[INDEX2(0,1,2)];
                            const Scalar A_10 = A_p[INDEX2(1,0,2)];
                            const Scalar A_11 = A_p[INDEX2(1,1,2)];
                            const Scalar tmp0 = 6.*w1*(A_01 - A_10);
                            const Scalar tmp1 = 6.*w1*(A_01 + A_10);
                            const Scalar tmp2 = 6.*w1*(-A_01 - A_10);
                            const Scalar tmp3 = 6.*w1*(-A_01 + A_10);
                            EM_S[INDEX2(0,0,4)]+=-8.*A_00*w6 + 8.*A_11*w3 + tmp1;
                            EM_S[INDEX2(0,1,4)]+=8.*A_00*w6 + 4.*A_11*w3 + tmp0;
                            EM_S[INDEX2(0,2,4)]+=-4.*A_00*w6 - 8.*A_11*w3 + tmp3;
                            EM_S[INDEX2(0,3,4)]+=4.*A_00*w6 - 4.*A_11*w3 + tmp2;
                            EM_S[INDEX2(1,0,4)]+=8.*A_00*w6 + 4.*A_11*w3 + tmp3;
                            EM_S[INDEX2(1,1,4)]+=-8.*A_00*w6 + 8.*A_11*w3 + tmp2;
                            EM_S[INDEX2(1,2,4)]+=4.*A_00*w6 - 4.*A_11*w3 + tmp1;
                            EM_S[INDEX2(1,3,4)]+=-4.*A_00*w6 - 8.*A_11*w3 + tmp0;
                            EM_S[INDEX2(2,0,4)]+=-4.*A_00*w6 - 8.*A_11*w3 + tmp0;
                            EM_S[INDEX2(2,1,4)]+=4.*A_00*w6 - 4.*A_11*w3 + tmp1;
                            EM_S[INDEX2(2,2,4)]+=-8.*A_00*w6 + 8.*A_11*w3 + tmp2;
                            EM_S[INDEX2(2,3,4)]+=8.*A_00*w6 + 4.*A_11*w3 + tmp3;
                            EM_S[INDEX2(3,0,4)]+=4.*A_00*w6 - 4.*A_11*w3 + tmp2;
                            EM_S[INDEX2(3,1,4)]+=-4.*A_00*w6 - 8.*A_11*w3 + tmp3;
                            EM_S[INDEX2(3,2,4)]+=8.*A_00*w6 + 4.*A_11*w3 + tmp0;
                            EM_S[INDEX2(3,3,4)]+=-8.*A_00*w6 + 8.*A_11*w3 + tmp1;
                        }
                    }
                    ///////////////
                    // process B //
                    ///////////////
                    if (!B.isEmpty()) {
                        const Scalar* B_p = B.getSampleDataRO(e, zero);
                        if (B.actsExpanded()) {
                            const Scalar B_0_0 = B_p[INDEX2(0,0,2)];
                            const Scalar B_1_0 = B_p[INDEX2(1,0,2)];
                            const Scalar B_0_1 = B_p[INDEX2(0,1,2)];
                            const Scalar B_1_1 = B_p[INDEX2(1,1,2)];
                            const Scalar B_0_2 = B_p[INDEX2(0,2,2)];
                            const Scalar B_1_2 = B_p[INDEX2(1,2,2)];
                            const Scalar B_0_3 = B_p[INDEX2(0,3,2)];
                            const Scalar B_1_3 = B_p[INDEX2(1,3,2)];
                            const Scalar tmp0 = w11*(B_1_0 + B_1_1);
                            const Scalar tmp1 = w14*(B_1_2 + B_1_3);
                            const Scalar tmp2 = w15*(-B_0_1 - B_0_3);
                            const Scalar tmp3 = w10*(-B_0_0 - B_0_2);
                            const Scalar tmp4 = w11*(B_1_2 + B_1_3);
                            const Scalar tmp5 = w14*(B_1_0 + B_1_1);
                            const Scalar tmp6 = w11*(-B_1_2 - B_1_3);
                            const Scalar tmp7 = w14*(-B_1_0 - B_1_1);
                            const Scalar tmp8 = w11*(-B_1_0 - B_1_1);
                            const Scalar tmp9 = w14*(-B_1_2 - B_1_3);
                            const Scalar tmp10 = w10*(-B_0_1 - B_0_3);
                            const Scalar tmp11 = w15*(-B_0_0 - B_0_2);
                            const Scalar tmp12 = w15*(B_0_0 + B_0_2);
                            const Scalar tmp13 = w10*(B_0_1 + B_0_3);
                            const Scalar tmp14 = w10*(B_0_0 + B_0_2);
                            const Scalar tmp15 = w15*(B_0_1 + B_0_3);
                            EM_S[INDEX2(0,0,4)]+=B_0_0*w12 + B_0_1*w10 + B_0_2*w15 + B_0_3*w13 + B_1_0*w16 + B_1_1*w14 + B_1_2*w11 + B_1_3*w17;
                            EM_S[INDEX2(0,1,4)]+=B_0_0*w10 + B_0_1*w12 + B_0_2*w13 + B_0_3*w15 + tmp0 + tmp1;
                            EM_S[INDEX2(0,2,4)]+=B_1_0*w11 + B_1_1*w17 + B_1_2*w16 + B_1_3*w14 + tmp14 + tmp15;
                            EM_S[INDEX2(0,3,4)]+=tmp12 + tmp13 + tmp4 + tmp5;
                            EM_S[INDEX2(1,0,4)]+=-B_0_0*w12 - B_0_1*w10 - B_0_2*w15 - B_0_3*w13 + tmp0 + tmp1;
                            EM_S[INDEX2(1,1,4)]+=-B_0_0*w10 - B_0_1*w12 - B_0_2*w13 - B_0_3*w15 + B_1_0*w14 + B_1_1*w16 + B_1_2*w17 + B_1_3*w11;
                            EM_S[INDEX2(1,2,4)]+=tmp2 + tmp3 + tmp4 + tmp5;
                            EM_S[INDEX2(1,3,4)]+=B_1_0*w17 + B_1_1*w11 + B_1_2*w14 + B_1_3*w16 + tmp10 + tmp11;
                            EM_S[INDEX2(2,0,4)]+=-B_1_0*w16 - B_1_1*w14 - B_1_2*w11 - B_1_3*w17 + tmp14 + tmp15;
                            EM_S[INDEX2(2,1,4)]+=tmp12 + tmp13 + tmp8 + tmp9;
                            EM_S[INDEX2(2,2,4)]+=B_0_0*w15 + B_0_1*w13 + B_0_2*w12 + B_0_3*w10 - B_1_0*w11 - B_1_1*w17 - B_1_2*w16 - B_1_3*w14;
                            EM_S[INDEX2(2,3,4)]+=B_0_0*w13 + B_0_1*w15 + B_0_2*w10 + B_0_3*w12 + tmp6 + tmp7;
                            EM_S[INDEX2(3,0,4)]+=tmp2 + tmp3 + tmp8 + tmp9;
                            EM_S[INDEX2(3,1,4)]+=-B_1_0*w14 - B_1_1*w16 - B_1_2*w17 - B_1_3*w11 + tmp10 + tmp11;
                            EM_S[INDEX2(3,2,4)]+=-B_0_0*w15 - B_0_1*w13 - B_0_2*w12 - B_0_3*w10 + tmp6 + tmp7;
                            EM_S[INDEX2(3,3,4)]+=-B_0_0*w13 - B_0_1*w15 - B_0_2*w10 - B_0_3*w12 - B_1_0*w17 - B_1_1*w11 - B_1_2*w14 - B_1_3*w16;
                        } else { // constant data
                            const Scalar B_0 = B_p[0];
                            const Scalar B_1 = B_p[1];
                            EM_S[INDEX2(0,0,4)]+= 2.*B_0*w18 + 2.*B_1*w19;
                            EM_S[INDEX2(0,1,4)]+= 2.*B_0*w18 +    B_1*w19;
                            EM_S[INDEX2(0,2,4)]+=    B_0*w18 + 2.*B_1*w19;
                            EM_S[INDEX2(0,3,4)]+=    B_0*w18 +    B_1*w19;
                            EM_S[INDEX2(1,0,4)]+=-2.*B_0*w18 +    B_1*w19;
                            EM_S[INDEX2(1,1,4)]+=-2.*B_0*w18 + 2.*B_1*w19;
                            EM_S[INDEX2(1,2,4)]+=   -B_0*w18 +    B_1*w19;
                            EM_S[INDEX2(1,3,4)]+=   -B_0*w18 + 2.*B_1*w19;
                            EM_S[INDEX2(2,0,4)]+=    B_0*w18 - 2.*B_1*w19;
                            EM_S[INDEX2(2,1,4)]+=    B_0*w18 -    B_1*w19;
                            EM_S[INDEX2(2,2,4)]+= 2.*B_0*w18 - 2.*B_1*w19;
                            EM_S[INDEX2(2,3,4)]+= 2.*B_0*w18 -    B_1*w19;
                            EM_S[INDEX2(3,0,4)]+=   -B_0*w18 -    B_1*w19;
                            EM_S[INDEX2(3,1,4)]+=   -B_0*w18 - 2.*B_1*w19;
                            EM_S[INDEX2(3,2,4)]+=-2.*B_0*w18 -    B_1*w19;
                            EM_S[INDEX2(3,3,4)]+=-2.*B_0*w18 - 2.*B_1*w19;
                        }
                    }
                    ///////////////
                    // process C //
                    ///////////////
                    if (!C.isEmpty()) {
                        const Scalar* C_p = C.getSampleDataRO(e, zero);
                        if (C.actsExpanded()) {
                            const Scalar C_0_0 = C_p[INDEX2(0,0,2)];
                            const Scalar C_1_0 = C_p[INDEX2(1,0,2)];
                            const Scalar C_0_1 = C_p[INDEX2(0,1,2)];
                            const Scalar C_1_1 = C_p[INDEX2(1,1,2)];
                            const Scalar C_0_2 = C_p[INDEX2(0,2,2)];
                            const Scalar C_1_2 = C_p[INDEX2(1,2,2)];
                            const Scalar C_0_3 = C_p[INDEX2(0,3,2)];
                            const Scalar C_1_3 = C_p[INDEX2(1,3,2)];
                            const Scalar tmp0 = w11*(C_1_0 + C_1_1);
                            const Scalar tmp1 = w14*(C_1_2 + C_1_3);
                            const Scalar tmp2 = w15*(C_0_0 + C_0_2);
                            const Scalar tmp3 = w10*(C_0_1 + C_0_3);
                            const Scalar tmp4 = w11*(-C_1_0 - C_1_1);
                            const Scalar tmp5 = w14*(-C_1_2 - C_1_3);
                            const Scalar tmp6 = w11*(-C_1_2 - C_1_3);
                            const Scalar tmp7 = w14*(-C_1_0 - C_1_1);
                            const Scalar tmp8 = w11*(C_1_2 + C_1_3);
                            const Scalar tmp9 = w14*(C_1_0 + C_1_1);
                            const Scalar tmp10 = w10*(-C_0_1 - C_0_3);
                            const Scalar tmp11 = w15*(-C_0_0 - C_0_2);
                            const Scalar tmp12 = w15*(-C_0_1 - C_0_3);
                            const Scalar tmp13 = w10*(-C_0_0 - C_0_2);
                            const Scalar tmp14 = w10*(C_0_0 + C_0_2);
                            const Scalar tmp15 = w15*(C_0_1 + C_0_3);
                            EM_S[INDEX2(0,0,4)]+=C_0_0*w12 + C_0_1*w10 + C_0_2*w15 + C_0_3*w13 + C_1_0*w16 + C_1_1*w14 + C_1_2*w11 + C_1_3*w17;
                            EM_S[INDEX2(0,1,4)]+=-C_0_0*w12 - C_0_1*w10 - C_0_2*w15 - C_0_3*w13 + tmp0 + tmp1;
                            EM_S[INDEX2(0,2,4)]+=-C_1_0*w16 - C_1_1*w14 - C_1_2*w11 - C_1_3*w17 + tmp14 + tmp15;
                            EM_S[INDEX2(0,3,4)]+=tmp12 + tmp13 + tmp4 + tmp5;
                            EM_S[INDEX2(1,0,4)]+=C_0_0*w10 + C_0_1*w12 + C_0_2*w13 + C_0_3*w15 + tmp0 + tmp1;
                            EM_S[INDEX2(1,1,4)]+=-C_0_0*w10 - C_0_1*w12 - C_0_2*w13 - C_0_3*w15 + C_1_0*w14 + C_1_1*w16 + C_1_2*w17 + C_1_3*w11;
                            EM_S[INDEX2(1,2,4)]+=tmp2 + tmp3 + tmp4 + tmp5;
                            EM_S[INDEX2(1,3,4)]+=-C_1_0*w14 - C_1_1*w16 - C_1_2*w17 - C_1_3*w11 + tmp10 + tmp11;
                            EM_S[INDEX2(2,0,4)]+=C_1_0*w11 + C_1_1*w17 + C_1_2*w16 + C_1_3*w14 + tmp14 + tmp15;
                            EM_S[INDEX2(2,1,4)]+=tmp12 + tmp13 + tmp8 + tmp9;
                            EM_S[INDEX2(2,2,4)]+=C_0_0*w15 + C_0_1*w13 + C_0_2*w12 + C_0_3*w10 - C_1_0*w11 - C_1_1*w17 - C_1_2*w16 - C_1_3*w14;
                            EM_S[INDEX2(2,3,4)]+=-C_0_0*w15 - C_0_1*w13 - C_0_2*w12 - C_0_3*w10 + tmp6 + tmp7;
                            EM_S[INDEX2(3,0,4)]+=tmp2 + tmp3 + tmp8 + tmp9;
                            EM_S[INDEX2(3,1,4)]+=C_1_0*w17 + C_1_1*w11 + C_1_2*w14 + C_1_3*w16 + tmp10 + tmp11;
                            EM_S[INDEX2(3,2,4)]+=C_0_0*w13 + C_0_1*w15 + C_0_2*w10 + C_0_3*w12 + tmp6 + tmp7;
                            EM_S[INDEX2(3,3,4)]+=-C_0_0*w13 - C_0_1*w15 - C_0_2*w10 - C_0_3*w12 - C_1_0*w17 - C_1_1*w11 - C_1_2*w14 - C_1_3*w16;
                        } else { // constant data
                            const Scalar C_0 = C_p[0];
                            const Scalar C_1 = C_p[1];
                            EM_S[INDEX2(0,0,4)]+= 2.*C_0*w18 + 2.*C_1*w19;
                            EM_S[INDEX2(0,1,4)]+=-2.*C_0*w18 +    C_1*w19;
                            EM_S[INDEX2(0,2,4)]+=    C_0*w18 - 2.*C_1*w19;
                            EM_S[INDEX2(0,3,4)]+=   -C_0*w18 -    C_1*w19;
                            EM_S[INDEX2(1,0,4)]+= 2.*C_0*w18 +    C_1*w19;
                            EM_S[INDEX2(1,1,4)]+=-2.*C_0*w18 + 2.*C_1*w19;
                            EM_S[INDEX2(1,2,4)]+=    C_0*w18 -    C_1*w19;
                            EM_S[INDEX2(1,3,4)]+=   -C_0*w18 - 2.*C_1*w19;
                            EM_S[INDEX2(2,0,4)]+=    C_0*w18 + 2.*C_1*w19;
                            EM_S[INDEX2(2,1,4)]+=   -C_0*w18 +    C_1*w19;
                            EM_S[INDEX2(2,2,4)]+= 2.*C_0*w18 - 2.*C_1*w19;
                            EM_S[INDEX2(2,3,4)]+=-2.*C_0*w18 -    C_1*w19;
                            EM_S[INDEX2(3,0,4)]+=    C_0*w18 +    C_1*w19;
                            EM_S[INDEX2(3,1,4)]+=   -C_0*w18 + 2.*C_1*w19;
                            EM_S[INDEX2(3,2,4)]+= 2.*C_0*w18 -    C_1*w19;
                            EM_S[INDEX2(3,3,4)]+=-2.*C_0*w18 - 2.*C_1*w19;
                        }
                    }
                    ///////////////
                    // process D //
                    ///////////////
                    if (!D.isEmpty()) {
                        const Scalar* D_p = D.getSampleDataRO(e, zero);
                        if (D.actsExpanded()) {
                            const Scalar D_0 = D_p[0];
                            const Scalar D_1 = D_p[1];
                            const Scalar D_2 = D_p[2];
                            const Scalar D_3 = D_p[3];
                            const Scalar tmp0 = w21*(D_2 + D_3);
                            const Scalar tmp1 = w20*(D_0 + D_1);
                            const Scalar tmp2 = w22*(D_0 + D_1 + D_2 + D_3);
                            const Scalar tmp3 = w21*(D_0 + D_1);
                            const Scalar tmp4 = w20*(D_2 + D_3);
                            const Scalar tmp5 = w22*(D_1 + D_2);
                            const Scalar tmp6 = w21*(D_0 + D_2);
                            const Scalar tmp7 = w20*(D_1 + D_3);
                            const Scalar tmp8 = w21*(D_1 + D_3);
                            const Scalar tmp9 = w20*(D_0 + D_2);
                            const Scalar tmp10 = w22*(D_0 + D_3);
                            EM_S[INDEX2(0,0,4)]+=D_0*w23 + D_3*w24 + tmp5;
                            EM_S[INDEX2(0,1,4)]+=tmp0 + tmp1;
                            EM_S[INDEX2(0,2,4)]+=tmp8 + tmp9;
                            EM_S[INDEX2(0,3,4)]+=tmp2;
                            EM_S[INDEX2(1,0,4)]+=tmp0 + tmp1;
                            EM_S[INDEX2(1,1,4)]+=D_1*w23 + D_2*w24 + tmp10;
                            EM_S[INDEX2(1,2,4)]+=tmp2;
                            EM_S[INDEX2(1,3,4)]+=tmp6 + tmp7;
                            EM_S[INDEX2(2,0,4)]+=tmp8 + tmp9;
                            EM_S[INDEX2(2,1,4)]+=tmp2;
                            EM_S[INDEX2(2,2,4)]+=D_1*w24 + D_2*w23 + tmp10;
                            EM_S[INDEX2(2,3,4)]+=tmp3 + tmp4;
                            EM_S[INDEX2(3,0,4)]+=tmp2;
                            EM_S[INDEX2(3,1,4)]+=tmp6 + tmp7;
                            EM_S[INDEX2(3,2,4)]+=tmp3 + tmp4;
                            EM_S[INDEX2(3,3,4)]+=D_0*w24 + D_3*w23 + tmp5;
                        } else { // constant data
                            const Scalar D_0 = D_p[0];
                            EM_S[INDEX2(0,0,4)]+=16.*D_0*w22;
                            EM_S[INDEX2(0,1,4)]+= 8.*D_0*w22;
                            EM_S[INDEX2(0,2,4)]+= 8.*D_0*w22;
                            EM_S[INDEX2(0,3,4)]+= 4.*D_0*w22;
                            EM_S[INDEX2(1,0,4)]+= 8.*D_0*w22;
                            EM_S[INDEX2(1,1,4)]+=16.*D_0*w22;
                            EM_S[INDEX2(1,2,4)]+= 4.*D_0*w22;
                            EM_S[INDEX2(1,3,4)]+= 8.*D_0*w22;
                            EM_S[INDEX2(2,0,4)]+= 8.*D_0*w22;
                            EM_S[INDEX2(2,1,4)]+= 4.*D_0*w22;
                            EM_S[INDEX2(2,2,4)]+=16.*D_0*w22;
                            EM_S[INDEX2(2,3,4)]+= 8.*D_0*w22;
                            EM_S[INDEX2(3,0,4)]+= 4.*D_0*w22;
                            EM_S[INDEX2(3,1,4)]+= 8.*D_0*w22;
                            EM_S[INDEX2(3,2,4)]+= 8.*D_0*w22;
                            EM_S[INDEX2(3,3,4)]+=16.*D_0*w22;
                        }
                    }
                    ///////////////
                    // process X //
                    ///////////////
                    if (!X.isEmpty()) {
                        const Scalar* X_p = X.getSampleDataRO(e, zero);
                        if (X.actsExpanded()) {
                            const Scalar X_0_0 = X_p[INDEX2(0,0,2)];
                            const Scalar X_1_0 = X_p[INDEX2(1,0,2)];
                            const Scalar X_0_1 = X_p[INDEX2(0,1,2)];
                            const Scalar X_1_1 = X_p[INDEX2(1,1,2)];
                            const Scalar X_0_2 = X_p[INDEX2(0,2,2)];
                            const Scalar X_1_2 = X_p[INDEX2(1,2,2)];
                            const Scalar X_0_3 = X_p[INDEX2(0,3,2)];
                            const Scalar X_1_3 = X_p[INDEX2(1,3,2)];
                            const Scalar tmp0 = 6.*w15*(X_0_2 + X_0_3);
                            const Scalar tmp1 = 6.*w10*(X_0_0 + X_0_1);
                            const Scalar tmp2 = 6.*w11*(X_1_0 + X_1_2);
                            const Scalar tmp3 = 6.*w14*(X_1_1 + X_1_3);
                            const Scalar tmp4 = 6.*w11*(X_1_1 + X_1_3);
                            const Scalar tmp5 = w25*(X_0_0 + X_0_1);
                            const Scalar tmp6 = w26*(X_0_2 + X_0_3);
                            const Scalar tmp7 = 6.*w14*(X_1_0 + X_1_2);
                            const Scalar tmp8 = w27*(X_1_0 + X_1_2);
                            const Scalar tmp9 = w28*(X_1_1 + X_1_3);
                            const Scalar tmp10 = w25*(-X_0_2 - X_0_3);
                            const Scalar tmp11 = w26*(-X_0_0 - X_0_1);
                            const Scalar tmp12 = w27*(X_1_1 + X_1_3);
                            const Scalar tmp13 = w28*(X_1_0 + X_1_2);
                            const Scalar tmp14 = w25*(X_0_2 + X_0_3);
                            const Scalar tmp15 = w26*(X_0_0 + X_0_1);
                            EM_F[0]+=tmp0 + tmp1 + tmp2 + tmp3;
                            EM_F[1]+=tmp4 + tmp5 + tmp6 + tmp7;
                            EM_F[2]+=tmp10 + tmp11 + tmp8 + tmp9;
                            EM_F[3]+=tmp12 + tmp13 + tmp14 + tmp15;
                        } else { // constant data
                            const Scalar X_0 = X_p[0];
                            const Scalar X_1 = X_p[1];
                            EM_F[0]+= 6.*X_0*w18 + 6.*X_1*w19;
                            EM_F[1]+=-6.*X_0*w18 + 6.*X_1*w19;
                            EM_F[2]+= 6.*X_0*w18 - 6.*X_1*w19;
                            EM_F[3]+=-6.*X_0*w18 - 6.*X_1*w19;
                        }
                    }
                    ///////////////
                    // process Y //
                    ///////////////
                    if (!Y.isEmpty()) {
                        const Scalar* Y_p = Y.getSampleDataRO(e, zero);
                        if (Y.actsExpanded()) {
                            const Scalar Y_0 = Y_p[0];
                            const Scalar Y_1 = Y_p[1];
                            const Scalar Y_2 = Y_p[2];
                            const Scalar Y_3 = Y_p[3];
                            const Scalar tmp0 = 6.*w22*(Y_1 + Y_2);
                            const Scalar tmp1 = 6.*w22*(Y_0 + Y_3);
                            EM_F[0]+=6.*Y_0*w20 + 6.*Y_3*w21 + tmp0;
                            EM_F[1]+=6.*Y_1*w20 + 6.*Y_2*w21 + tmp1;
                            EM_F[2]+=6.*Y_1*w21 + 6.*Y_2*w20 + tmp1;
                            EM_F[3]+=6.*Y_0*w21 + 6.*Y_3*w20 + tmp0;
                        } else { // constant data
                            EM_F[0]+=36.*Y_p[0]*w22;
                            EM_F[1]+=36.*Y_p[0]*w22;
                            EM_F[2]+=36.*Y_p[0]*w22;
                            EM_F[3]+=36.*Y_p[0]*w22;
                        }
                    }

                    // add to matrix (if addEM_S) and RHS (if addEM_F)
                    const index_t firstNode = m_NN[0]*k1+k0;
                    domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, addEM_S, addEM_F, firstNode);
                } // end k0 loop
            } // end k1 loop
        } // end of colouring
    } // end of parallel region
}

/****************************************************************************/
// PDE SINGLE BOUNDARY
/****************************************************************************/

template<class Scalar>
void DefaultAssembler2D<Scalar>::assemblePDEBoundarySingle(
                                        AbstractSystemMatrix* mat, Data& rhs,
                                        const Data& d, const Data& y) const
{
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
    const int NE0 = m_NE[0];
    const int NE1 = m_NE[1];
    const bool addEM_S = !d.isEmpty();
    const bool addEM_F = !y.isEmpty();
    const Scalar zero = static_cast<Scalar>(0);
    rhs.requireWrite();

#pragma omp parallel
    {
        vector<Scalar> EM_S(4*4);
        vector<Scalar> EM_F(4);

        if (domain->m_faceOffset[0] > -1) {
            if (addEM_S)
                fill(EM_S.begin(), EM_S.end(), zero);
            if (addEM_F) {
                EM_F[1] = zero;
                EM_F[3] = zero;
            }

            for (index_t k1_0=0; k1_0<2; k1_0++) { // colouring
#pragma omp for
                for (index_t k1=k1_0; k1<NE1; k1+=2) {
                    ///////////////
                    // process d //
                    ///////////////
                    if (addEM_S) {
                        const Scalar* d_p = d.getSampleDataRO(k1, zero);
                        if (d.actsExpanded()) {
                            const Scalar d_0 = d_p[0];
                            const Scalar d_1 = d_p[1];
                            const Scalar tmp0 = w2*(d_0 + d_1);
                            EM_S[INDEX2(0,0,4)] = d_0*w0 + d_1*w1;
                            EM_S[INDEX2(2,0,4)] = tmp0;
                            EM_S[INDEX2(0,2,4)] = tmp0;
                            EM_S[INDEX2(2,2,4)] = d_0*w1 + d_1*w0;
                        } else { // constant data
                            EM_S[INDEX2(0,0,4)] = 4.*d_p[0]*w2;
                            EM_S[INDEX2(2,0,4)] = 2.*d_p[0]*w2;
                            EM_S[INDEX2(0,2,4)] = 2.*d_p[0]*w2;
                            EM_S[INDEX2(2,2,4)] = 4.*d_p[0]*w2;
                        }
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (addEM_F) {
                        const Scalar* y_p = y.getSampleDataRO(k1, zero);
                        if (y.actsExpanded()) {
                            EM_F[0] = w3*y_p[0] + w4*y_p[1];
                            EM_F[2] = w3*y_p[1] + w4*y_p[0];
                        } else { // constant data
                            EM_F[0] = 6.*w2*y_p[0];
                            EM_F[2] = 6.*w2*y_p[0];
                        }
                    }
                    const index_t firstNode=m_NN[0]*k1;
                    domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, addEM_S,
                                              addEM_F, firstNode);
                }
            } // end colouring
        }

        if (domain->m_faceOffset[1] > -1) {
            if (addEM_S)
                fill(EM_S.begin(), EM_S.end(), zero);
            if (addEM_F) {
                EM_F[0] = zero;
                EM_F[2] = zero;
            }

            for (index_t k1_0=0; k1_0<2; k1_0++) { // colouring            
#pragma omp for
                for (index_t k1=k1_0; k1<NE1; k1+=2) {
                    const index_t e = domain->m_faceOffset[1]+k1;
                    ///////////////
                    // process d //
                    ///////////////
                    if (addEM_S) {
                        const Scalar* d_p = d.getSampleDataRO(e, zero);
                        if (d.actsExpanded()) {
                            const Scalar d_0 = d_p[0];
                            const Scalar d_1 = d_p[1];
                            const Scalar tmp0 = w2*(d_0 + d_1);
                            EM_S[INDEX2(1,1,4)] = d_0*w0 + d_1*w1;
                            EM_S[INDEX2(3,1,4)] = tmp0;
                            EM_S[INDEX2(1,3,4)] = tmp0;
                            EM_S[INDEX2(3,3,4)] = d_0*w1 + d_1*w0;
                        } else { // constant data
                            EM_S[INDEX2(1,1,4)] = 4.*d_p[0]*w2;
                            EM_S[INDEX2(3,1,4)] = 2.*d_p[0]*w2;
                            EM_S[INDEX2(1,3,4)] = 2.*d_p[0]*w2;
                            EM_S[INDEX2(3,3,4)] = 4.*d_p[0]*w2;
                        }
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (addEM_F) {
                        const Scalar* y_p = y.getSampleDataRO(e, zero);
                        if (y.actsExpanded()) {
                            EM_F[1] = w3*y_p[0] + w4*y_p[1];
                            EM_F[3] = w3*y_p[1] + w4*y_p[0];
                        } else { // constant data
                            EM_F[1] = 6.*w2*y_p[0];
                            EM_F[3] = 6.*w2*y_p[0];
                        }
                    }
                    const index_t firstNode = m_NN[0]*(k1+1)-2;
                    domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, addEM_S,
                                              addEM_F, firstNode);
                }
            } // end colouring
        }

        if (domain->m_faceOffset[2] > -1) {
            if (addEM_S)
                fill(EM_S.begin(), EM_S.end(), zero);
            if (addEM_F) {
                EM_F[2] = zero;
                EM_F[3] = zero;
            }

            for (index_t k0_0=0; k0_0<2; k0_0++) { // colouring
#pragma omp for
                for (index_t k0 = k0_0; k0 < NE0; k0+=2) {
                    const index_t e = domain->m_faceOffset[2]+k0;
                    ///////////////
                    // process d //
                    ///////////////
                    if (addEM_S) {
                        const Scalar* d_p = d.getSampleDataRO(e, zero);
                        if (d.actsExpanded()) {
                            const Scalar d_0 = d_p[0];
                            const Scalar d_1 = d_p[1];
                            const Scalar tmp0 = w5*(d_0 + d_1);
                            EM_S[INDEX2(0,0,4)] = d_0*w6 + d_1*w7;
                            EM_S[INDEX2(1,0,4)] = tmp0;
                            EM_S[INDEX2(0,1,4)] = tmp0;
                            EM_S[INDEX2(1,1,4)] = d_0*w7 + d_1*w6;
                        } else { // constant data
                            EM_S[INDEX2(0,0,4)] = 4.*d_p[0]*w5;
                            EM_S[INDEX2(1,0,4)] = 2.*d_p[0]*w5;
                            EM_S[INDEX2(0,1,4)] = 2.*d_p[0]*w5;
                            EM_S[INDEX2(1,1,4)] = 4.*d_p[0]*w5;
                        }
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (addEM_F) {
                        const Scalar* y_p = y.getSampleDataRO(e, zero);
                        if (y.actsExpanded()) {
                            EM_F[0] = w8*y_p[0] + w9*y_p[1];
                            EM_F[1] = w8*y_p[1] + w9*y_p[0];
                        } else { // constant data
                            EM_F[0] = 6.*w5*y_p[0];
                            EM_F[1] = 6.*w5*y_p[0];
                        }
                    }
                    const index_t firstNode = k0;
                    domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, addEM_S,
                                              addEM_F, firstNode);
                }
            } // end colouring
        }

        if (domain->m_faceOffset[3] > -1) {
            if (addEM_S)
                fill(EM_S.begin(), EM_S.end(), zero);
            if (addEM_F) {
                EM_F[0] = zero;
                EM_F[1] = zero;
            }

            for (index_t k0_0=0; k0_0<2; k0_0++) { // colouring
#pragma omp for
                for (index_t k0 = k0_0; k0 < NE0; k0+=2) {
                    const index_t e = domain->m_faceOffset[3]+k0;
                    ///////////////
                    // process d //
                    ///////////////
                    if (addEM_S) {
                        const Scalar* d_p = d.getSampleDataRO(e, zero);
                        if (d.actsExpanded()) {
                            const Scalar d_0 = d_p[0];
                            const Scalar d_1 = d_p[1];
                            const Scalar tmp0 = w5*(d_0 + d_1);
                            EM_S[INDEX2(2,2,4)] = d_0*w6 + d_1*w7;
                            EM_S[INDEX2(3,2,4)] = tmp0;
                            EM_S[INDEX2(2,3,4)] = tmp0;
                            EM_S[INDEX2(3,3,4)] = d_0*w7 + d_1*w6;
                        } else { // constant data
                            EM_S[INDEX2(2,2,4)] = 4.*d_p[0]*w5;
                            EM_S[INDEX2(3,2,4)] = 2.*d_p[0]*w5;
                            EM_S[INDEX2(2,3,4)] = 2.*d_p[0]*w5;
                            EM_S[INDEX2(3,3,4)] = 4.*d_p[0]*w5;
                        }
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (addEM_F) {
                        const Scalar* y_p = y.getSampleDataRO(e, zero);
                        if (y.actsExpanded()) {
                            EM_F[2] = w8*y_p[0] + w9*y_p[1];
                            EM_F[3] = w8*y_p[1] + w9*y_p[0];
                        } else { // constant data
                            EM_F[2] = 6.*w5*y_p[0];
                            EM_F[3] = 6.*w5*y_p[0];
                        }
                    }
                    const index_t firstNode=m_NN[0]*(m_NN[1]-2)+k0;
                    domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, addEM_S,
                                              addEM_F, firstNode);
                }
            } // end colouring
        }
    } // end of parallel section
}

/****************************************************************************/
// PDE SINGLE REDUCED
/****************************************************************************/

template<class Scalar>
void DefaultAssembler2D<Scalar>::assemblePDESingleReduced(
                                    AbstractSystemMatrix* mat,
                                    Data& rhs, const Data& A, const Data& B,
                                    const Data& C, const Data& D,
                                    const Data& X, const Data& Y) const
{
    const double w0 = 1./4;
    const double w1 = m_dx[0]/8;
    const double w2 = m_dx[1]/8;
    const double w3 = m_dx[0]*m_dx[1]/16;
    const double w4 = m_dx[0]/(4*m_dx[1]);
    const double w5 = m_dx[1]/(4*m_dx[0]);
    const int NE0 = m_NE[0];
    const int NE1 = m_NE[1];
    const bool addEM_S = (!A.isEmpty() || !B.isEmpty() || !C.isEmpty() || !D.isEmpty());
    const bool addEM_F = (!X.isEmpty() || !Y.isEmpty());
    const Scalar zero = static_cast<Scalar>(0);
    rhs.requireWrite();

#pragma omp parallel
    {
        vector<Scalar> EM_S(4*4, zero);
        vector<Scalar> EM_F(4, zero);

        for (index_t k1_0=0; k1_0<2; k1_0++) { // colouring
#pragma omp for
            for (index_t k1=k1_0; k1<NE1; k1+=2) {
                for (index_t k0=0; k0<NE0; ++k0)  {
                    const index_t e = k0 + NE0*k1;
                    if (addEM_S)
                        fill(EM_S.begin(), EM_S.end(), zero);
                    if (addEM_F)
                        fill(EM_F.begin(), EM_F.end(), zero);
                    ///////////////
                    // process A //
                    ///////////////
                    if (!A.isEmpty()) {
                        const Scalar* A_p = A.getSampleDataRO(e, zero);
                        const Scalar A_00 = A_p[INDEX2(0,0,2)];
                        const Scalar A_10 = A_p[INDEX2(1,0,2)];
                        const Scalar A_01 = A_p[INDEX2(0,1,2)];
                        const Scalar A_11 = A_p[INDEX2(1,1,2)];
                        const Scalar tmp0 = (A_01 + A_10)*w0;
                        const Scalar tmp1 = A_00*w5;
                        const Scalar tmp2 = A_01*w0;
                        const Scalar tmp3 = A_10*w0;
                        const Scalar tmp4 = A_11*w4;
                        EM_S[INDEX2(0,0,4)]+=tmp4 + tmp0 + tmp1;
                        EM_S[INDEX2(1,0,4)]+=tmp4 - tmp1 + tmp3 - tmp2;
                        EM_S[INDEX2(2,0,4)]+=tmp2 - tmp3 - tmp4 + tmp1;
                        EM_S[INDEX2(3,0,4)]+=-tmp1 - tmp4 - tmp0;
                        EM_S[INDEX2(0,1,4)]+=tmp4 - tmp1 + tmp2 - tmp3;
                        EM_S[INDEX2(1,1,4)]+=tmp4 + tmp1 - tmp0;
                        EM_S[INDEX2(2,1,4)]+=-tmp1 + tmp0 - tmp4;
                        EM_S[INDEX2(3,1,4)]+=-tmp4 + tmp1 + tmp3 - tmp2;
                        EM_S[INDEX2(0,2,4)]+=-tmp4 + tmp1 + tmp3 - tmp2;
                        EM_S[INDEX2(1,2,4)]+=-tmp1 + tmp0 - tmp4;
                        EM_S[INDEX2(2,2,4)]+=tmp4 + tmp1 - tmp0;
                        EM_S[INDEX2(3,2,4)]+=tmp4 - tmp1 + tmp2 - tmp3;
                        EM_S[INDEX2(0,3,4)]+=-tmp1 - tmp4 - tmp0;
                        EM_S[INDEX2(1,3,4)]+=tmp2 - tmp3 - tmp4 + tmp1;
                        EM_S[INDEX2(2,3,4)]+=tmp4 - tmp1 + tmp3 - tmp2;
                        EM_S[INDEX2(3,3,4)]+=tmp4 + tmp0 + tmp1;
                    }
                    ///////////////
                    // process B //
                    ///////////////
                    if (!B.isEmpty()) {
                        const Scalar* B_p = B.getSampleDataRO(e, zero);
                        const Scalar tmp0 = B_p[0]*w2;
                        const Scalar tmp1 = B_p[1]*w1;
                        EM_S[INDEX2(0,0,4)]+=-tmp0 - tmp1;
                        EM_S[INDEX2(1,0,4)]+= tmp0 - tmp1;
                        EM_S[INDEX2(2,0,4)]+= tmp1 - tmp0;
                        EM_S[INDEX2(3,0,4)]+= tmp0 + tmp1;
                        EM_S[INDEX2(0,1,4)]+=-tmp0 - tmp1;
                        EM_S[INDEX2(1,1,4)]+= tmp0 - tmp1;
                        EM_S[INDEX2(2,1,4)]+= tmp1 - tmp0;
                        EM_S[INDEX2(3,1,4)]+= tmp0 + tmp1;
                        EM_S[INDEX2(0,2,4)]+=-tmp0 - tmp1;
                        EM_S[INDEX2(1,2,4)]+= tmp0 - tmp1;
                        EM_S[INDEX2(2,2,4)]+= tmp1 - tmp0;
                        EM_S[INDEX2(3,2,4)]+= tmp0 + tmp1;
                        EM_S[INDEX2(0,3,4)]+=-tmp0 - tmp1;
                        EM_S[INDEX2(1,3,4)]+= tmp0 - tmp1;
                        EM_S[INDEX2(2,3,4)]+= tmp1 - tmp0;
                        EM_S[INDEX2(3,3,4)]+= tmp0 + tmp1;
                    }
                    ///////////////
                    // process C //
                    ///////////////
                    if (!C.isEmpty()) {
                        const Scalar* C_p = C.getSampleDataRO(e, zero);
                        const Scalar tmp0 = C_p[0]*w2;
                        const Scalar tmp1 = C_p[1]*w1;
                        EM_S[INDEX2(0,0,4)]+=-tmp1 - tmp0;
                        EM_S[INDEX2(1,0,4)]+=-tmp1 - tmp0;
                        EM_S[INDEX2(2,0,4)]+=-tmp1 - tmp0;
                        EM_S[INDEX2(3,0,4)]+=-tmp1 - tmp0;
                        EM_S[INDEX2(0,1,4)]+= tmp0 - tmp1;
                        EM_S[INDEX2(1,1,4)]+= tmp0 - tmp1;
                        EM_S[INDEX2(2,1,4)]+= tmp0 - tmp1;
                        EM_S[INDEX2(3,1,4)]+= tmp0 - tmp1;
                        EM_S[INDEX2(0,2,4)]+= tmp1 - tmp0;
                        EM_S[INDEX2(1,2,4)]+= tmp1 - tmp0;
                        EM_S[INDEX2(2,2,4)]+= tmp1 - tmp0;
                        EM_S[INDEX2(3,2,4)]+= tmp1 - tmp0;
                        EM_S[INDEX2(0,3,4)]+= tmp0 + tmp1;
                        EM_S[INDEX2(1,3,4)]+= tmp0 + tmp1;
                        EM_S[INDEX2(2,3,4)]+= tmp0 + tmp1;
                        EM_S[INDEX2(3,3,4)]+= tmp0 + tmp1;
                    }
                    ///////////////
                    // process D //
                    ///////////////
                    if (!D.isEmpty()) {
                        const Scalar* D_p = D.getSampleDataRO(e, zero);
                        EM_S[INDEX2(0,0,4)]+=D_p[0]*w3;
                        EM_S[INDEX2(1,0,4)]+=D_p[0]*w3;
                        EM_S[INDEX2(2,0,4)]+=D_p[0]*w3;
                        EM_S[INDEX2(3,0,4)]+=D_p[0]*w3;
                        EM_S[INDEX2(0,1,4)]+=D_p[0]*w3;
                        EM_S[INDEX2(1,1,4)]+=D_p[0]*w3;
                        EM_S[INDEX2(2,1,4)]+=D_p[0]*w3;
                        EM_S[INDEX2(3,1,4)]+=D_p[0]*w3;
                        EM_S[INDEX2(0,2,4)]+=D_p[0]*w3;
                        EM_S[INDEX2(1,2,4)]+=D_p[0]*w3;
                        EM_S[INDEX2(2,2,4)]+=D_p[0]*w3;
                        EM_S[INDEX2(3,2,4)]+=D_p[0]*w3;
                        EM_S[INDEX2(0,3,4)]+=D_p[0]*w3;
                        EM_S[INDEX2(1,3,4)]+=D_p[0]*w3;
                        EM_S[INDEX2(2,3,4)]+=D_p[0]*w3;
                        EM_S[INDEX2(3,3,4)]+=D_p[0]*w3;
                    }
                    ///////////////
                    // process X //
                    ///////////////
                    if (!X.isEmpty()) {
                        const Scalar* X_p = X.getSampleDataRO(e, zero);
                        const Scalar wX0 = 4.*X_p[0]*w2;
                        const Scalar wX1 = 4.*X_p[1]*w1;
                        EM_F[0]+=-wX0 - wX1;
                        EM_F[1]+=-wX1 + wX0;
                        EM_F[2]+=-wX0 + wX1;
                        EM_F[3]+= wX0 + wX1;
                    }
                    ///////////////
                    // process Y //
                    ///////////////
                    if (!Y.isEmpty()) {
                        const Scalar* Y_p = Y.getSampleDataRO(e, zero);
                        EM_F[0]+=4.*Y_p[0]*w3;
                        EM_F[1]+=4.*Y_p[0]*w3;
                        EM_F[2]+=4.*Y_p[0]*w3;
                        EM_F[3]+=4.*Y_p[0]*w3;
                    }

                    // add to matrix (if addEM_S) and RHS (if addEM_F)
                    const index_t firstNode=m_NN[0]*k1+k0;
                    domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, addEM_S,
                                              addEM_F, firstNode);
                } // end k0 loop
            } // end k1 loop
        } // end of colouring
    } // end of parallel region
}

/****************************************************************************/
// PDE SINGLE REDUCED BOUNDARY
/****************************************************************************/

template<class Scalar>
void DefaultAssembler2D<Scalar>::assemblePDEBoundarySingleReduced(
                                        AbstractSystemMatrix* mat, Data& rhs,
                                        const Data& d, const Data& y) const
{
    const double w0 = m_dx[0]/4;
    const double w1 = m_dx[1]/4;
    const int NE0 = m_NE[0];
    const int NE1 = m_NE[1];
    const bool addEM_S = !d.isEmpty();
    const bool addEM_F = !y.isEmpty();
    const Scalar zero = static_cast<Scalar>(0);
    rhs.requireWrite();

#pragma omp parallel
    {
        vector<Scalar> EM_S(4*4, zero);
        vector<Scalar> EM_F(4, zero);

        if (domain->m_faceOffset[0] > -1) {
            if (addEM_S)
                fill(EM_S.begin(), EM_S.end(), zero);
            if (addEM_F) {
                EM_F[1] = zero;
                EM_F[3] = zero;
            }

            for (index_t k1_0=0; k1_0<2; k1_0++) { // colouring
#pragma omp for
                for (index_t k1=k1_0; k1<NE1; k1+=2) {
                    ///////////////
                    // process d //
                    ///////////////
                    if (addEM_S) {
                        const Scalar* d_p = d.getSampleDataRO(k1, zero);
                        EM_S[INDEX2(0,0,4)] = d_p[0]*w1;
                        EM_S[INDEX2(2,0,4)] = d_p[0]*w1;
                        EM_S[INDEX2(0,2,4)] = d_p[0]*w1;
                        EM_S[INDEX2(2,2,4)] = d_p[0]*w1;
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (addEM_F) {
                        const Scalar* y_p = y.getSampleDataRO(k1, zero);
                        EM_F[0] = 2.*w1*y_p[0];
                        EM_F[2] = 2.*w1*y_p[0];
                    }
                    const index_t firstNode=m_NN[0]*k1;
                    domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, addEM_S,
                                              addEM_F, firstNode);
                }
            } // end colouring
        }

        if (domain->m_faceOffset[1] > -1) {
            if (addEM_S)
                fill(EM_S.begin(), EM_S.end(), zero);
            if (addEM_F) {
                EM_F[0] = zero;
                EM_F[2] = zero;
            }

            for (index_t k1_0=0; k1_0<2; k1_0++) { // colouring            
#pragma omp for
                for (index_t k1=k1_0; k1<NE1; k1+=2) {
                    const index_t e = domain->m_faceOffset[1]+k1;
                    ///////////////
                    // process d //
                    ///////////////
                    if (addEM_S) {
                        const Scalar* d_p = d.getSampleDataRO(e, zero);
                        EM_S[INDEX2(1,1,4)] = d_p[0]*w1;
                        EM_S[INDEX2(3,1,4)] = d_p[0]*w1;
                        EM_S[INDEX2(1,3,4)] = d_p[0]*w1;
                        EM_S[INDEX2(3,3,4)] = d_p[0]*w1;
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (addEM_F) {
                        const Scalar* y_p = y.getSampleDataRO(e, zero);
                        EM_F[1] = 2.*w1*y_p[0];
                        EM_F[3] = 2.*w1*y_p[0];
                    }
                    const index_t firstNode=m_NN[0]*(k1+1)-2;
                    domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, addEM_S,
                                              addEM_F, firstNode);
                }
            } // end colouring
        }

        if (domain->m_faceOffset[2] > -1) {
            if (addEM_S)
                fill(EM_S.begin(), EM_S.end(), zero);
            if (addEM_F) {
                EM_F[2] = zero;
                EM_F[3] = zero;
            }

            for (index_t k0_0=0; k0_0<2; k0_0++) { // colouring
#pragma omp for
                for (index_t k0 = k0_0; k0 < NE0; k0+=2) {
                    const index_t e = domain->m_faceOffset[2]+k0;
                    ///////////////
                    // process d //
                    ///////////////
                    if (addEM_S) {
                        const Scalar* d_p = d.getSampleDataRO(e, zero);
                        EM_S[INDEX2(0,0,4)] = d_p[0]*w0;
                        EM_S[INDEX2(1,0,4)] = d_p[0]*w0;
                        EM_S[INDEX2(0,1,4)] = d_p[0]*w0;
                        EM_S[INDEX2(1,1,4)] = d_p[0]*w0;
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (addEM_F) {
                        const Scalar* y_p = y.getSampleDataRO(e, zero);
                        EM_F[0] = 2.*w0*y_p[0];
                        EM_F[1] = 2.*w0*y_p[0];
                    }
                    const index_t firstNode = k0;
                    domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, addEM_S,
                                              addEM_F, firstNode);
                }
            } // end colouring
        }

        if (domain->m_faceOffset[3] > -1) {
            if (addEM_S)
                fill(EM_S.begin(), EM_S.end(), zero);
            if (addEM_F) {
                EM_F[0] = zero;
                EM_F[1] = zero;
            }

            for (index_t k0_0=0; k0_0<2; k0_0++) { // colouring
#pragma omp for
                for (index_t k0 = k0_0; k0 < NE0; k0+=2) {
                    const index_t e = domain->m_faceOffset[3]+k0;
                    ///////////////
                    // process d //
                    ///////////////
                    if (addEM_S) {
                        const Scalar* d_p = d.getSampleDataRO(e, zero);
                        EM_S[INDEX2(2,2,4)] = d_p[0]*w0;
                        EM_S[INDEX2(3,2,4)] = d_p[0]*w0;
                        EM_S[INDEX2(2,3,4)] = d_p[0]*w0;
                        EM_S[INDEX2(3,3,4)] = d_p[0]*w0;
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (addEM_F) {
                        const Scalar* y_p = y.getSampleDataRO(e, zero);
                        EM_F[2] = 2.*w0*y_p[0];
                        EM_F[3] = 2.*w0*y_p[0];
                    }
                    const index_t firstNode=m_NN[0]*(m_NN[1]-2)+k0;
                    domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, addEM_S,
                                              addEM_F, firstNode);
                }
            } // end colouring
        }
    } // end of parallel section
}

/****************************************************************************/
// PDE SYSTEM
/****************************************************************************/

template<class Scalar>
void DefaultAssembler2D<Scalar>::assemblePDESystem(AbstractSystemMatrix* mat,
                                     Data& rhs, const Data& A, const Data& B,
                                     const Data& C, const Data& D,
                                     const Data& X, const Data& Y) const
{
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
    const bool addEM_F = (!X.isEmpty() || !Y.isEmpty());
    const Scalar zero = static_cast<Scalar>(0);
    rhs.requireWrite();

#pragma omp parallel
    {
        vector<Scalar> EM_S(4*4*numEq*numComp, zero);
        vector<Scalar> EM_F(4*numEq, zero);

        for (index_t k1_0=0; k1_0<2; k1_0++) { // colouring
#pragma omp for
            for (index_t k1=k1_0; k1 < NE1; k1+=2) {
                for (index_t k0=0; k0 < NE0; ++k0)  {
                    const index_t e = k0 + NE0*k1;
                    if (addEM_S)
                        fill(EM_S.begin(), EM_S.end(), zero);
                    if (addEM_F)
                        fill(EM_F.begin(), EM_F.end(), zero);
                    ///////////////
                    // process A //
                    ///////////////
                    if (!A.isEmpty()) {
                        const Scalar* A_p = A.getSampleDataRO(e, zero);
                        if (A.actsExpanded()) {
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const Scalar A_00_0 = A_p[INDEX5(k,0,m,0,0,numEq,2,numComp,2)];
                                    const Scalar A_01_0 = A_p[INDEX5(k,0,m,1,0,numEq,2,numComp,2)];
                                    const Scalar A_10_0 = A_p[INDEX5(k,1,m,0,0,numEq,2,numComp,2)];
                                    const Scalar A_11_0 = A_p[INDEX5(k,1,m,1,0,numEq,2,numComp,2)];
                                    const Scalar A_00_1 = A_p[INDEX5(k,0,m,0,1,numEq,2,numComp,2)];
                                    const Scalar A_01_1 = A_p[INDEX5(k,0,m,1,1,numEq,2,numComp,2)];
                                    const Scalar A_10_1 = A_p[INDEX5(k,1,m,0,1,numEq,2,numComp,2)];
                                    const Scalar A_11_1 = A_p[INDEX5(k,1,m,1,1,numEq,2,numComp,2)];
                                    const Scalar A_00_2 = A_p[INDEX5(k,0,m,0,2,numEq,2,numComp,2)];
                                    const Scalar A_01_2 = A_p[INDEX5(k,0,m,1,2,numEq,2,numComp,2)];
                                    const Scalar A_10_2 = A_p[INDEX5(k,1,m,0,2,numEq,2,numComp,2)];
                                    const Scalar A_11_2 = A_p[INDEX5(k,1,m,1,2,numEq,2,numComp,2)];
                                    const Scalar A_00_3 = A_p[INDEX5(k,0,m,0,3,numEq,2,numComp,2)];
                                    const Scalar A_01_3 = A_p[INDEX5(k,0,m,1,3,numEq,2,numComp,2)];
                                    const Scalar A_10_3 = A_p[INDEX5(k,1,m,0,3,numEq,2,numComp,2)];
                                    const Scalar A_11_3 = A_p[INDEX5(k,1,m,1,3,numEq,2,numComp,2)];
                                    const Scalar tmp0 = w3*(A_11_0 + A_11_1 + A_11_2 + A_11_3);
                                    const Scalar tmp1 = w1*(A_01_0 + A_01_3 - A_10_1 - A_10_2);
                                    const Scalar tmp2 = w4*(A_00_2 + A_00_3);
                                    const Scalar tmp3 = w0*(A_00_0 + A_00_1);
                                    const Scalar tmp4 = w5*(A_01_2 - A_10_3);
                                    const Scalar tmp5 = w2*(-A_01_1 + A_10_0);
                                    const Scalar tmp6 = w5*(A_01_3 + A_10_0);
                                    const Scalar tmp7 = w3*(-A_11_0 - A_11_1 - A_11_2 - A_11_3);
                                    const Scalar tmp8 = w6*(A_00_0 + A_00_1 + A_00_2 + A_00_3);
                                    const Scalar tmp9 = w1*(A_01_1 + A_01_2 + A_10_1 + A_10_2);
                                    const Scalar tmp10 = w2*(-A_01_0 - A_10_3);
                                    const Scalar tmp11 = w4*(A_00_0 + A_00_1);
                                    const Scalar tmp12 = w0*(A_00_2 + A_00_3);
                                    const Scalar tmp13 = w5*(A_01_1 - A_10_0);
                                    const Scalar tmp14 = w2*(-A_01_2 + A_10_3);
                                    const Scalar tmp15 = w7*(A_11_0 + A_11_2);
                                    const Scalar tmp16 = w4*(-A_00_2 - A_00_3);
                                    const Scalar tmp17 = w0*(-A_00_0 - A_00_1);
                                    const Scalar tmp18 = w5*(A_01_3 + A_10_3);
                                    const Scalar tmp19 = w8*(A_11_1 + A_11_3);
                                    const Scalar tmp20 = w2*(-A_01_0 - A_10_0);
                                    const Scalar tmp21 = w7*(A_11_1 + A_11_3);
                                    const Scalar tmp22 = w4*(-A_00_0 - A_00_1);
                                    const Scalar tmp23 = w0*(-A_00_2 - A_00_3);
                                    const Scalar tmp24 = w5*(A_01_0 + A_10_0);
                                    const Scalar tmp25 = w8*(A_11_0 + A_11_2);
                                    const Scalar tmp26 = w2*(-A_01_3 - A_10_3);
                                    const Scalar tmp27 = w5*(-A_01_1 - A_10_2);
                                    const Scalar tmp28 = w1*(-A_01_0 - A_01_3 - A_10_0 - A_10_3);
                                    const Scalar tmp29 = w2*(A_01_2 + A_10_1);
                                    const Scalar tmp30 = w7*(-A_11_1 - A_11_3);
                                    const Scalar tmp31 = w1*(-A_01_1 - A_01_2 + A_10_0 + A_10_3);
                                    const Scalar tmp32 = w5*(-A_01_0 + A_10_2);
                                    const Scalar tmp33 = w8*(-A_11_0 - A_11_2);
                                    const Scalar tmp34 = w6*(-A_00_0 - A_00_1 - A_00_2 - A_00_3);
                                    const Scalar tmp35 = w2*(A_01_3 - A_10_1);
                                    const Scalar tmp36 = w5*(A_01_0 + A_10_3);
                                    const Scalar tmp37 = w2*(-A_01_3 - A_10_0);
                                    const Scalar tmp38 = w7*(-A_11_0 - A_11_2);
                                    const Scalar tmp39 = w5*(-A_01_3 + A_10_1);
                                    const Scalar tmp40 = w8*(-A_11_1 - A_11_3);
                                    const Scalar tmp41 = w2*(A_01_0 - A_10_2);
                                    const Scalar tmp42 = w5*(A_01_1 - A_10_3);
                                    const Scalar tmp43 = w2*(-A_01_2 + A_10_0);
                                    const Scalar tmp44 = w5*(A_01_2 - A_10_0);
                                    const Scalar tmp45 = w2*(-A_01_1 + A_10_3);
                                    const Scalar tmp46 = w5*(-A_01_0 + A_10_1);
                                    const Scalar tmp47 = w2*(A_01_3 - A_10_2);
                                    const Scalar tmp48 = w5*(-A_01_1 - A_10_1);
                                    const Scalar tmp49 = w2*(A_01_2 + A_10_2);
                                    const Scalar tmp50 = w5*(-A_01_3 + A_10_2);
                                    const Scalar tmp51 = w2*(A_01_0 - A_10_1);
                                    const Scalar tmp52 = w5*(-A_01_2 - A_10_1);
                                    const Scalar tmp53 = w2*(A_01_1 + A_10_2);
                                    const Scalar tmp54 = w5*(-A_01_2 - A_10_2);
                                    const Scalar tmp55 = w2*(A_01_1 + A_10_1);
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
                                    const Scalar A_00 = A_p[INDEX4(k,0,m,0, numEq,2, numComp)];
                                    const Scalar A_01 = A_p[INDEX4(k,0,m,1, numEq,2, numComp)];
                                    const Scalar A_10 = A_p[INDEX4(k,1,m,0, numEq,2, numComp)];
                                    const Scalar A_11 = A_p[INDEX4(k,1,m,1, numEq,2, numComp)];
                                    const Scalar tmp0 = 6.*w1*(A_01 - A_10);
                                    const Scalar tmp1 = 6.*w1*(A_01 + A_10);
                                    const Scalar tmp2 = 6.*w1*(-A_01 - A_10);
                                    const Scalar tmp3 = 6.*w1*(-A_01 + A_10);
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=-8.*A_00*w6 + 8.*A_11*w3 + tmp1;
                                    EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+= 8.*A_00*w6 + 4.*A_11*w3 + tmp0;
                                    EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=-4.*A_00*w6 - 8.*A_11*w3 + tmp3;
                                    EM_S[INDEX4(k,m,0,3,numEq,numComp,4)]+= 4.*A_00*w6 - 4.*A_11*w3 + tmp2;
                                    EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+= 8.*A_00*w6 + 4.*A_11*w3 + tmp3;
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=-8.*A_00*w6 + 8.*A_11*w3 + tmp2;
                                    EM_S[INDEX4(k,m,1,2,numEq,numComp,4)]+= 4.*A_00*w6 - 4.*A_11*w3 + tmp1;
                                    EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+=-4.*A_00*w6 - 8.*A_11*w3 + tmp0;
                                    EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=-4.*A_00*w6 - 8.*A_11*w3 + tmp0;
                                    EM_S[INDEX4(k,m,2,1,numEq,numComp,4)]+= 4.*A_00*w6 - 4.*A_11*w3 + tmp1;
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=-8.*A_00*w6 + 8.*A_11*w3 + tmp2;
                                    EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+= 8.*A_00*w6 + 4.*A_11*w3 + tmp3;
                                    EM_S[INDEX4(k,m,3,0,numEq,numComp,4)]+= 4.*A_00*w6 - 4.*A_11*w3 + tmp2;
                                    EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+=-4.*A_00*w6 - 8.*A_11*w3 + tmp3;
                                    EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+= 8.*A_00*w6 + 4.*A_11*w3 + tmp0;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=-8.*A_00*w6 + 8.*A_11*w3 + tmp1;
                                }
                            }
                        }
                    }
                    ///////////////
                    // process B //
                    ///////////////
                    if (!B.isEmpty()) {
                        const Scalar* B_p = B.getSampleDataRO(e, zero);
                        if (B.actsExpanded()) {
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const Scalar B_0_0 = B_p[INDEX4(k,0,m,0, numEq,2,numComp)];
                                    const Scalar B_1_0 = B_p[INDEX4(k,1,m,0, numEq,2,numComp)];
                                    const Scalar B_0_1 = B_p[INDEX4(k,0,m,1, numEq,2,numComp)];
                                    const Scalar B_1_1 = B_p[INDEX4(k,1,m,1, numEq,2,numComp)];
                                    const Scalar B_0_2 = B_p[INDEX4(k,0,m,2, numEq,2,numComp)];
                                    const Scalar B_1_2 = B_p[INDEX4(k,1,m,2, numEq,2,numComp)];
                                    const Scalar B_0_3 = B_p[INDEX4(k,0,m,3, numEq,2,numComp)];
                                    const Scalar B_1_3 = B_p[INDEX4(k,1,m,3, numEq,2,numComp)];
                                    const Scalar tmp0 = w11*(B_1_0 + B_1_1);
                                    const Scalar tmp1 = w14*(B_1_2 + B_1_3);
                                    const Scalar tmp2 = w15*(-B_0_1 - B_0_3);
                                    const Scalar tmp3 = w10*(-B_0_0 - B_0_2);
                                    const Scalar tmp4 = w11*(B_1_2 + B_1_3);
                                    const Scalar tmp5 = w14*(B_1_0 + B_1_1);
                                    const Scalar tmp6 = w11*(-B_1_2 - B_1_3);
                                    const Scalar tmp7 = w14*(-B_1_0 - B_1_1);
                                    const Scalar tmp8 = w11*(-B_1_0 - B_1_1);
                                    const Scalar tmp9 = w14*(-B_1_2 - B_1_3);
                                    const Scalar tmp10 = w10*(-B_0_1 - B_0_3);
                                    const Scalar tmp11 = w15*(-B_0_0 - B_0_2);
                                    const Scalar tmp12 = w15*(B_0_0 + B_0_2);
                                    const Scalar tmp13 = w10*(B_0_1 + B_0_3);
                                    const Scalar tmp14 = w10*(B_0_0 + B_0_2);
                                    const Scalar tmp15 = w15*(B_0_1 + B_0_3);
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
                                    const Scalar wB0 = B_p[INDEX3(k,0,m,numEq,2)]*w18;
                                    const Scalar wB1 = B_p[INDEX3(k,1,m,numEq,2)]*w19;
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+= 2.*wB0 + 2.*wB1;
                                    EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+= 2.*wB0 +    wB1;
                                    EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=    wB0 + 2.*wB1;
                                    EM_S[INDEX4(k,m,0,3,numEq,numComp,4)]+=    wB0 +    wB1;
                                    EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+=-2.*wB0 +    wB1;
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=-2.*wB0 + 2.*wB1;
                                    EM_S[INDEX4(k,m,1,2,numEq,numComp,4)]+=   -wB0 +    wB1;
                                    EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+=   -wB0 + 2.*wB1;
                                    EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=    wB0 - 2.*wB1;
                                    EM_S[INDEX4(k,m,2,1,numEq,numComp,4)]+=    wB0 -    wB1;
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+= 2.*wB0 - 2.*wB1;
                                    EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+= 2.*wB0 -    wB1;
                                    EM_S[INDEX4(k,m,3,0,numEq,numComp,4)]+=   -wB0 -    wB1;
                                    EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+=   -wB0 - 2.*wB1;
                                    EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+=-2.*wB0 -    wB1;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=-2.*wB0 - 2.*wB1;
                                }
                            }
                        }
                    }
                    ///////////////
                    // process C //
                    ///////////////
                    if (!C.isEmpty()) {
                        const Scalar* C_p = C.getSampleDataRO(e, zero);
                        if (C.actsExpanded()) {
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const Scalar C_0_0 = C_p[INDEX4(k,m,0, 0, numEq,numComp,2)];
                                    const Scalar C_1_0 = C_p[INDEX4(k,m,1, 0, numEq,numComp,2)];
                                    const Scalar C_0_1 = C_p[INDEX4(k,m,0, 1, numEq,numComp,2)];
                                    const Scalar C_1_1 = C_p[INDEX4(k,m,1, 1, numEq,numComp,2)];
                                    const Scalar C_0_2 = C_p[INDEX4(k,m,0, 2, numEq,numComp,2)];
                                    const Scalar C_1_2 = C_p[INDEX4(k,m,1, 2, numEq,numComp,2)];
                                    const Scalar C_0_3 = C_p[INDEX4(k,m,0, 3, numEq,numComp,2)];
                                    const Scalar C_1_3 = C_p[INDEX4(k,m,1, 3, numEq,numComp,2)];
                                    const Scalar tmp0 = w11*(C_1_0 + C_1_1);
                                    const Scalar tmp1 = w14*(C_1_2 + C_1_3);
                                    const Scalar tmp2 = w15*(C_0_0 + C_0_2);
                                    const Scalar tmp3 = w10*(C_0_1 + C_0_3);
                                    const Scalar tmp4 = w11*(-C_1_0 - C_1_1);
                                    const Scalar tmp5 = w14*(-C_1_2 - C_1_3);
                                    const Scalar tmp6 = w11*(-C_1_2 - C_1_3);
                                    const Scalar tmp7 = w14*(-C_1_0 - C_1_1);
                                    const Scalar tmp8 = w11*(C_1_2 + C_1_3);
                                    const Scalar tmp9 = w14*(C_1_0 + C_1_1);
                                    const Scalar tmp10 = w10*(-C_0_1 - C_0_3);
                                    const Scalar tmp11 = w15*(-C_0_0 - C_0_2);
                                    const Scalar tmp12 = w15*(-C_0_1 - C_0_3);
                                    const Scalar tmp13 = w10*(-C_0_0 - C_0_2);
                                    const Scalar tmp14 = w10*(C_0_0 + C_0_2);
                                    const Scalar tmp15 = w15*(C_0_1 + C_0_3);
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
                                    const Scalar wC0 = C_p[INDEX3(k,m,0,numEq,numComp)]*w18;
                                    const Scalar wC1 = C_p[INDEX3(k,m,1,numEq,numComp)]*w19;
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+= 2.*wC0 + 2.*wC1;
                                    EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+=-2.*wC0 +    wC1;
                                    EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=    wC0 - 2.*wC1;
                                    EM_S[INDEX4(k,m,0,3,numEq,numComp,4)]+=   -wC0 -    wC1;
                                    EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+= 2.*wC0 +    wC1;
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=-2.*wC0 + 2.*wC1;
                                    EM_S[INDEX4(k,m,1,2,numEq,numComp,4)]+=    wC0 -    wC1;
                                    EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+=   -wC0 - 2.*wC1;
                                    EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=    wC0 + 2.*wC1;
                                    EM_S[INDEX4(k,m,2,1,numEq,numComp,4)]+=   -wC0 +    wC1;
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+= 2.*wC0 - 2.*wC1;
                                    EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+=-2.*wC0 -    wC1;
                                    EM_S[INDEX4(k,m,3,0,numEq,numComp,4)]+=    wC0 +    wC1;
                                    EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+=   -wC0 + 2.*wC1;
                                    EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+= 2.*wC0 -    wC1;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=-2.*wC0 - 2.*wC1;
                                }
                            }
                        }
                    }
                    ///////////////
                    // process D //
                    ///////////////
                    if (!D.isEmpty()) {
                        const Scalar* D_p = D.getSampleDataRO(e, zero);
                        if (D.actsExpanded()) {
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const Scalar D_0 = D_p[INDEX3(k,m,0,numEq,numComp)];
                                    const Scalar D_1 = D_p[INDEX3(k,m,1,numEq,numComp)];
                                    const Scalar D_2 = D_p[INDEX3(k,m,2,numEq,numComp)];
                                    const Scalar D_3 = D_p[INDEX3(k,m,3,numEq,numComp)];
                                    const Scalar tmp0 = w21*(D_2 + D_3);
                                    const Scalar tmp1 = w20*(D_0 + D_1);
                                    const Scalar tmp2 = w22*(D_0 + D_1 + D_2 + D_3);
                                    const Scalar tmp3 = w21*(D_0 + D_1);
                                    const Scalar tmp4 = w20*(D_2 + D_3);
                                    const Scalar tmp5 = w22*(D_1 + D_2);
                                    const Scalar tmp6 = w21*(D_0 + D_2);
                                    const Scalar tmp7 = w20*(D_1 + D_3);
                                    const Scalar tmp8 = w21*(D_1 + D_3);
                                    const Scalar tmp9 = w20*(D_0 + D_2);
                                    const Scalar tmp10 = w22*(D_0 + D_3);
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
                                    const Scalar D_0 = D_p[INDEX2(k, m, numEq)];
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=16.*D_0*w22;
                                    EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+= 8.*D_0*w22;
                                    EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+= 8.*D_0*w22;
                                    EM_S[INDEX4(k,m,0,3,numEq,numComp,4)]+= 4.*D_0*w22;
                                    EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+= 8.*D_0*w22;
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=16.*D_0*w22;
                                    EM_S[INDEX4(k,m,1,2,numEq,numComp,4)]+= 4.*D_0*w22;
                                    EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+= 8.*D_0*w22;
                                    EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+= 8.*D_0*w22;
                                    EM_S[INDEX4(k,m,2,1,numEq,numComp,4)]+= 4.*D_0*w22;
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=16.*D_0*w22;
                                    EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+= 8.*D_0*w22;
                                    EM_S[INDEX4(k,m,3,0,numEq,numComp,4)]+= 4.*D_0*w22;
                                    EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+= 8.*D_0*w22;
                                    EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+= 8.*D_0*w22;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=16.*D_0*w22;
                                }
                            }
                        }
                    }
                    ///////////////
                    // process X //
                    ///////////////
                    if (!X.isEmpty()) {
                        const Scalar* X_p = X.getSampleDataRO(e, zero);
                        if (X.actsExpanded()) {
                            for (index_t k=0; k<numEq; k++) {
                                const Scalar X_0_0 = X_p[INDEX3(k,0,0,numEq,2)];
                                const Scalar X_1_0 = X_p[INDEX3(k,1,0,numEq,2)];
                                const Scalar X_0_1 = X_p[INDEX3(k,0,1,numEq,2)];
                                const Scalar X_1_1 = X_p[INDEX3(k,1,1,numEq,2)];
                                const Scalar X_0_2 = X_p[INDEX3(k,0,2,numEq,2)];
                                const Scalar X_1_2 = X_p[INDEX3(k,1,2,numEq,2)];
                                const Scalar X_0_3 = X_p[INDEX3(k,0,3,numEq,2)];
                                const Scalar X_1_3 = X_p[INDEX3(k,1,3,numEq,2)];
                                const Scalar tmp0 = 6.*w15*(X_0_2 + X_0_3);
                                const Scalar tmp1 = 6.*w10*(X_0_0 + X_0_1);
                                const Scalar tmp2 = 6.*w11*(X_1_0 + X_1_2);
                                const Scalar tmp3 = 6.*w14*(X_1_1 + X_1_3);
                                const Scalar tmp4 = 6.*w11*(X_1_1 + X_1_3);
                                const Scalar tmp5 =    w25*(X_0_0 + X_0_1);
                                const Scalar tmp6 =    w26*(X_0_2 + X_0_3);
                                const Scalar tmp7 = 6.*w14*(X_1_0 + X_1_2);
                                const Scalar tmp8 =    w27*(X_1_0 + X_1_2);
                                const Scalar tmp9 =    w28*(X_1_1 + X_1_3);
                                const Scalar tmp10 =   w25*(-X_0_2- X_0_3);
                                const Scalar tmp11 =   w26*(-X_0_0- X_0_1);
                                const Scalar tmp12 =   w27*(X_1_1 + X_1_3);
                                const Scalar tmp13 =   w28*(X_1_0 + X_1_2);
                                const Scalar tmp14 =   w25*(X_0_2 + X_0_3);
                                const Scalar tmp15 =   w26*(X_0_0 + X_0_1);
                                EM_F[INDEX2(k,0,numEq)]+=tmp0 + tmp1 + tmp2 + tmp3;
                                EM_F[INDEX2(k,1,numEq)]+=tmp4 + tmp5 + tmp6 + tmp7;
                                EM_F[INDEX2(k,2,numEq)]+=tmp10 + tmp11 + tmp8 + tmp9;
                                EM_F[INDEX2(k,3,numEq)]+=tmp12 + tmp13 + tmp14 + tmp15;
                            }
                        } else { // constant data
                            for (index_t k=0; k<numEq; k++) {
                                const Scalar wX0 = X_p[INDEX2(k, 0, numEq)]*w18;
                                const Scalar wX1 = X_p[INDEX2(k, 1, numEq)]*w19;
                                EM_F[INDEX2(k,0,numEq)]+= 6.*wX0 + 6.*wX1;
                                EM_F[INDEX2(k,1,numEq)]+=-6.*wX0 + 6.*wX1;
                                EM_F[INDEX2(k,2,numEq)]+= 6.*wX0 - 6.*wX1;
                                EM_F[INDEX2(k,3,numEq)]+=-6.*wX0 - 6.*wX1;
                            }
                        }
                    }
                    ///////////////
                    // process Y //
                    ///////////////
                    if (!Y.isEmpty()) {
                        const Scalar* Y_p = Y.getSampleDataRO(e, zero);
                        if (Y.actsExpanded()) {
                            for (index_t k=0; k<numEq; k++) {
                                const Scalar Y_0 = Y_p[INDEX2(k, 0, numEq)];
                                const Scalar Y_1 = Y_p[INDEX2(k, 1, numEq)];
                                const Scalar Y_2 = Y_p[INDEX2(k, 2, numEq)];
                                const Scalar Y_3 = Y_p[INDEX2(k, 3, numEq)];
                                const Scalar tmp0 = 6.*w22*(Y_1 + Y_2);
                                const Scalar tmp1 = 6.*w22*(Y_0 + Y_3);
                                EM_F[INDEX2(k,0,numEq)]+=6.*Y_0*w20 + 6.*Y_3*w21 + tmp0;
                                EM_F[INDEX2(k,1,numEq)]+=6.*Y_1*w20 + 6.*Y_2*w21 + tmp1;
                                EM_F[INDEX2(k,2,numEq)]+=6.*Y_1*w21 + 6.*Y_2*w20 + tmp1;
                                EM_F[INDEX2(k,3,numEq)]+=6.*Y_0*w21 + 6.*Y_3*w20 + tmp0;
                            }
                        } else { // constant data
                            for (index_t k=0; k<numEq; k++) {
                                EM_F[INDEX2(k,0,numEq)]+=36.*Y_p[k]*w22;
                                EM_F[INDEX2(k,1,numEq)]+=36.*Y_p[k]*w22;
                                EM_F[INDEX2(k,2,numEq)]+=36.*Y_p[k]*w22;
                                EM_F[INDEX2(k,3,numEq)]+=36.*Y_p[k]*w22;
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

/****************************************************************************/
// PDE SYSTEM BOUNDARY
/****************************************************************************/

template<class Scalar>
void DefaultAssembler2D<Scalar>::assemblePDEBoundarySystem(
                                        AbstractSystemMatrix* mat, Data& rhs,
                                        const Data& d, const Data& y) const
{
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
    const int NE0 = m_NE[0];
    const int NE1 = m_NE[1];
    const bool addEM_S = !d.isEmpty();
    const bool addEM_F = !y.isEmpty();
    const Scalar zero = static_cast<Scalar>(0);
    rhs.requireWrite();

#pragma omp parallel
    {
        vector<Scalar> EM_S(4*4*numEq*numComp, zero);
        vector<Scalar> EM_F(4*numEq, zero);

        if (domain->m_faceOffset[0] > -1) {
            if (addEM_S)
                fill(EM_S.begin(), EM_S.end(), zero);
            if (addEM_F)
                fill(EM_F.begin(), EM_F.end(), zero);

            for (index_t k1_0=0; k1_0<2; k1_0++) { // colouring
#pragma omp for
                for (index_t k1=k1_0; k1<NE1; k1+=2) {
                    const index_t e = k1;
                    ///////////////
                    // process d //
                    ///////////////
                    if (addEM_S) {
                        const Scalar* d_p = d.getSampleDataRO(e, zero);
                        if (d.actsExpanded()) {
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const Scalar d_0 = d_p[INDEX3(k,m,0,numEq,numComp)];
                                    const Scalar d_1 = d_p[INDEX3(k,m,1,numEq,numComp)];
                                    const Scalar tmp0 = w2*(d_0 + d_1);
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,4)] = d_0*w0 + d_1*w1;
                                    EM_S[INDEX4(k,m,2,0,numEq,numComp,4)] = tmp0;
                                    EM_S[INDEX4(k,m,0,2,numEq,numComp,4)] = tmp0;
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,4)] = d_0*w1 + d_1*w0;
                                }
                             }
                        } else { // constant data
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const Scalar d_0 = d_p[INDEX2(k, m, numEq)];
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,4)] = 4.*d_0*w2;
                                    EM_S[INDEX4(k,m,2,0,numEq,numComp,4)] = 2.*d_0*w2;
                                    EM_S[INDEX4(k,m,0,2,numEq,numComp,4)] = 2.*d_0*w2;
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,4)] = 4.*d_0*w2;
                                }
                            }
                        }
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (addEM_F) {
                        const Scalar* y_p = y.getSampleDataRO(e, zero);
                        if (y.actsExpanded()) {
                            for (index_t k=0; k<numEq; k++) {
                                const Scalar y_0 = y_p[INDEX2(k, 0, numEq)];
                                const Scalar y_1 = y_p[INDEX2(k, 1, numEq)];
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
                fill(EM_S.begin(), EM_S.end(), zero);
            if (addEM_F)
                fill(EM_F.begin(), EM_F.end(), zero);

            for (index_t k1_0=0; k1_0<2; k1_0++) { // colouring            
#pragma omp for
                for (index_t k1=k1_0; k1<NE1; k1+=2) {
                    const index_t e = domain->m_faceOffset[1]+k1;
                    ///////////////
                    // process d //
                    ///////////////
                    if (addEM_S) {
                        const Scalar* d_p = d.getSampleDataRO(e, zero);
                        if (d.actsExpanded()) {
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const Scalar d_0 = d_p[INDEX3(k,m,0,numEq,numComp)];
                                    const Scalar d_1 = d_p[INDEX3(k,m,1,numEq,numComp)];
                                    const Scalar tmp0 = w2*(d_0 + d_1);
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,4)] = d_0*w0 + d_1*w1;
                                    EM_S[INDEX4(k,m,3,1,numEq,numComp,4)] = tmp0;
                                    EM_S[INDEX4(k,m,1,3,numEq,numComp,4)] = tmp0;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,4)] = d_0*w1 + d_1*w0;
                                }
                             }
                        } else { // constant data
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const Scalar d_0 = d_p[INDEX2(k, m, numEq)];
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,4)] = 4.*d_0*w2;
                                    EM_S[INDEX4(k,m,3,1,numEq,numComp,4)] = 2.*d_0*w2;
                                    EM_S[INDEX4(k,m,1,3,numEq,numComp,4)] = 2.*d_0*w2;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,4)] = 4.*d_0*w2;
                                }
                            }
                        }
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (addEM_F) {
                        const Scalar* y_p = y.getSampleDataRO(e, zero);
                        if (y.actsExpanded()) {
                            for (index_t k=0; k<numEq; k++) {
                                const Scalar y_0 = y_p[INDEX2(k, 0, numEq)];
                                const Scalar y_1 = y_p[INDEX2(k, 1, numEq)];
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
                fill(EM_S.begin(), EM_S.end(), zero);
            if (addEM_F)
                fill(EM_F.begin(), EM_F.end(), zero);

            for (index_t k0_0=0; k0_0<2; k0_0++) { // colouring
#pragma omp for
                for (index_t k0 = k0_0; k0 < NE0; k0+=2) {
                    const index_t e = domain->m_faceOffset[2]+k0;
                    ///////////////
                    // process d //
                    ///////////////
                    if (addEM_S) {
                        const Scalar* d_p = d.getSampleDataRO(e, zero);
                        if (d.actsExpanded()) {
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const Scalar d_0 = d_p[INDEX3(k,m,0,numEq,numComp)];
                                    const Scalar d_1 = d_p[INDEX3(k,m,1,numEq,numComp)];
                                    const Scalar tmp0 = w5*(d_0 + d_1);
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,4)] = d_0*w6 + d_1*w7;
                                    EM_S[INDEX4(k,m,1,0,numEq,numComp,4)] = tmp0;
                                    EM_S[INDEX4(k,m,0,1,numEq,numComp,4)] = tmp0;
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,4)] = d_0*w7 + d_1*w6;
                                }
                             }
                        } else { // constant data
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const Scalar d_0 = d_p[INDEX2(k, m, numEq)];
                                    EM_S[INDEX4(k,m,0,0,numEq,numComp,4)] = 4.*d_0*w5;
                                    EM_S[INDEX4(k,m,1,0,numEq,numComp,4)] = 2.*d_0*w5;
                                    EM_S[INDEX4(k,m,0,1,numEq,numComp,4)] = 2.*d_0*w5;
                                    EM_S[INDEX4(k,m,1,1,numEq,numComp,4)] = 4.*d_0*w5;
                                }
                            }
                        }
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (addEM_F) {
                        const Scalar* y_p = y.getSampleDataRO(e, zero);
                        if (y.actsExpanded()) {
                            for (index_t k=0; k<numEq; k++) {
                                const Scalar y_0 = y_p[INDEX2(k, 0, numEq)];
                                const Scalar y_1 = y_p[INDEX2(k, 1, numEq)];
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
                fill(EM_S.begin(), EM_S.end(), zero);
            if (addEM_F)
                fill(EM_F.begin(), EM_F.end(), zero);

            for (index_t k0_0=0; k0_0<2; k0_0++) { // colouring
#pragma omp for
                for (index_t k0 = k0_0; k0 < NE0; k0+=2) {
                    const index_t e = domain->m_faceOffset[3]+k0;
                    ///////////////
                    // process d //
                    ///////////////
                    if (addEM_S) {
                        const Scalar* d_p = d.getSampleDataRO(e, zero);
                        if (d.actsExpanded()) {
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const Scalar d_0 = d_p[INDEX3(k,m,0,numEq,numComp)];
                                    const Scalar d_1 = d_p[INDEX3(k,m,1,numEq,numComp)];
                                    const Scalar tmp0 = w5*(d_0 + d_1);
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,4)] = d_0*w6 + d_1*w7;
                                    EM_S[INDEX4(k,m,3,2,numEq,numComp,4)] = tmp0;
                                    EM_S[INDEX4(k,m,2,3,numEq,numComp,4)] = tmp0;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,4)] = d_0*w7 + d_1*w6;
                                }
                             }
                        } else { // constant data
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const Scalar d_0 = d_p[INDEX2(k, m, numEq)];
                                    EM_S[INDEX4(k,m,2,2,numEq,numComp,4)] = 4.*d_0*w5;
                                    EM_S[INDEX4(k,m,3,2,numEq,numComp,4)] = 2.*d_0*w5;
                                    EM_S[INDEX4(k,m,2,3,numEq,numComp,4)] = 2.*d_0*w5;
                                    EM_S[INDEX4(k,m,3,3,numEq,numComp,4)] = 4.*d_0*w5;
                                }
                            }
                        }
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (addEM_F) {
                        const Scalar* y_p = y.getSampleDataRO(e, zero);
                        if (y.actsExpanded()) {
                            for (index_t k=0; k<numEq; k++) {
                                const Scalar y_0 = y_p[INDEX2(k, 0, numEq)];
                                const Scalar y_1 = y_p[INDEX2(k, 1, numEq)];
                                EM_F[INDEX2(k,2,numEq)] = w8*y_0 + w9*y_1;
                                EM_F[INDEX2(k,3,numEq)] = w8*y_1 + w9*y_0;
                            }
                        } else { // constant data
                            for (index_t k=0; k<numEq; k++) {
                                EM_F[INDEX2(k,2,numEq)] = 6.*w5*y_p[k];
                                EM_F[INDEX2(k,3,numEq)] = 6.*w5*y_p[k];
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

/****************************************************************************/
// PDE SYSTEM REDUCED
/****************************************************************************/

template<class Scalar>
void DefaultAssembler2D<Scalar>::assemblePDESystemReduced(
                                    AbstractSystemMatrix* mat,
                                    Data& rhs, const Data& A, const Data& B,
                                    const Data& C, const Data& D,
                                    const Data& X, const Data& Y) const
{
    dim_t numEq, numComp;
    if (!mat)
        numEq=numComp=(rhs.isEmpty() ? 1 : rhs.getDataPointSize());
    else {
        numEq=mat->getRowBlockSize();
        numComp=mat->getColumnBlockSize();
    }

    const double w0 = 1./4;
    const double w1 = m_dx[0]/8;
    const double w2 = m_dx[1]/8;
    const double w3 = m_dx[0]*m_dx[1]/16;
    const double w4 = m_dx[0]/(4*m_dx[1]);
    const double w5 = m_dx[1]/(4*m_dx[0]);
    const int NE0 = m_NE[0];
    const int NE1 = m_NE[1];
    const bool addEM_S = (!A.isEmpty() || !B.isEmpty() || !C.isEmpty() || !D.isEmpty());
    const bool addEM_F = (!X.isEmpty() || !Y.isEmpty());
    const Scalar zero = static_cast<Scalar>(0);
    rhs.requireWrite();

#pragma omp parallel
    {
        vector<Scalar> EM_S(4*4*numEq*numComp, zero);
        vector<Scalar> EM_F(4*numEq, zero);

        for (index_t k1_0=0; k1_0<2; k1_0++) { // colouring
#pragma omp for
            for (index_t k1=k1_0; k1<NE1; k1+=2) {
                for (index_t k0=0; k0<NE0; ++k0)  {
                    if (addEM_S)
                        fill(EM_S.begin(), EM_S.end(), zero);
                    if (addEM_F)
                        fill(EM_F.begin(), EM_F.end(), zero);
                    const index_t e = k0 + NE0*k1;
                    ///////////////
                    // process A //
                    ///////////////
                    if (!A.isEmpty()) {
                        const Scalar* A_p = A.getSampleDataRO(e, zero);
                        for (index_t k=0; k<numEq; k++) {
                            for (index_t m=0; m<numComp; m++) {
                                const Scalar Aw00 = A_p[INDEX4(k,0,m,0, numEq,2, numComp)]*w5;
                                const Scalar Aw01 = A_p[INDEX4(k,0,m,1, numEq,2, numComp)]*w0;
                                const Scalar Aw10 = A_p[INDEX4(k,1,m,0, numEq,2, numComp)]*w0;
                                const Scalar Aw11 = A_p[INDEX4(k,1,m,1, numEq,2, numComp)]*w4;
                                EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+= Aw00 + Aw01 + Aw10 + Aw11;
                                EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+=-Aw00 - Aw01 + Aw10 + Aw11;
                                EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+= Aw00 + Aw01 - Aw10 - Aw11;
                                EM_S[INDEX4(k,m,3,0,numEq,numComp,4)]+=-Aw00 - Aw01 - Aw10 - Aw11;
                                EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+=-Aw00 + Aw01 - Aw10 + Aw11;
                                EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+= Aw00 - Aw01 - Aw10 + Aw11;
                                EM_S[INDEX4(k,m,2,1,numEq,numComp,4)]+=-Aw00 + Aw01 + Aw10 - Aw11;
                                EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+= Aw00 - Aw01 + Aw10 - Aw11;
                                EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+= Aw00 - Aw01 + Aw10 - Aw11;
                                EM_S[INDEX4(k,m,1,2,numEq,numComp,4)]+=-Aw00 + Aw01 + Aw10 - Aw11;
                                EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+= Aw00 - Aw01 - Aw10 + Aw11;
                                EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+=-Aw00 + Aw01 - Aw10 + Aw11;
                                EM_S[INDEX4(k,m,0,3,numEq,numComp,4)]+=-Aw00 - Aw01 - Aw10 - Aw11;
                                EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+= Aw00 + Aw01 - Aw10 - Aw11;
                                EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+=-Aw00 - Aw01 + Aw10 + Aw11;
                                EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+= Aw00 + Aw01 + Aw10 + Aw11;
                            }
                        }
                    }
                    ///////////////
                    // process B //
                    ///////////////
                    if (!B.isEmpty()) {
                        const Scalar* B_p = B.getSampleDataRO(e, zero);
                        for (index_t k=0; k<numEq; k++) {
                            for (index_t m=0; m<numComp; m++) {
                                const Scalar wB0 = B_p[INDEX3(k,0,m, numEq, 2)]*w2;
                                const Scalar wB1 = B_p[INDEX3(k,1,m, numEq, 2)]*w1;
                                EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=-wB0 - wB1;
                                EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+=-wB0 - wB1;
                                EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=-wB0 - wB1;
                                EM_S[INDEX4(k,m,0,3,numEq,numComp,4)]+=-wB0 - wB1;
                                EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+= wB0 - wB1;
                                EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+= wB0 - wB1;
                                EM_S[INDEX4(k,m,1,2,numEq,numComp,4)]+= wB0 - wB1;
                                EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+= wB0 - wB1;
                                EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=-wB0 + wB1;
                                EM_S[INDEX4(k,m,2,1,numEq,numComp,4)]+=-wB0 + wB1;
                                EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=-wB0 + wB1;
                                EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+=-wB0 + wB1;
                                EM_S[INDEX4(k,m,3,0,numEq,numComp,4)]+= wB0 + wB1;
                                EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+= wB0 + wB1;
                                EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+= wB0 + wB1;
                                EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+= wB0 + wB1;
                            }
                        }
                    }
                    ///////////////
                    // process C //
                    ///////////////
                    if (!C.isEmpty()) {
                        const Scalar* C_p = C.getSampleDataRO(e, zero);
                        for (index_t k=0; k<numEq; k++) {
                            for (index_t m=0; m<numComp; m++) {
                                const Scalar wC0 = C_p[INDEX3(k, m, 0, numEq, numComp)]*w2;
                                const Scalar wC1 = C_p[INDEX3(k, m, 1, numEq, numComp)]*w1;
                                EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=-wC0 - wC1;
                                EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+=-wC0 - wC1;
                                EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=-wC0 - wC1;
                                EM_S[INDEX4(k,m,3,0,numEq,numComp,4)]+=-wC0 - wC1;
                                EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+= wC0 - wC1;
                                EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+= wC0 - wC1;
                                EM_S[INDEX4(k,m,2,1,numEq,numComp,4)]+= wC0 - wC1;
                                EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+= wC0 - wC1;
                                EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=-wC0 + wC1;
                                EM_S[INDEX4(k,m,1,2,numEq,numComp,4)]+=-wC0 + wC1;
                                EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=-wC0 + wC1;
                                EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+=-wC0 + wC1;
                                EM_S[INDEX4(k,m,0,3,numEq,numComp,4)]+= wC0 + wC1;
                                EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+= wC0 + wC1;
                                EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+= wC0 + wC1;
                                EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+= wC0 + wC1;
                            }
                        }
                    }
                    ///////////////
                    // process D //
                    ///////////////
                    if (!D.isEmpty()) {
                        const Scalar* D_p = D.getSampleDataRO(e, zero);
                        for (index_t k=0; k<numEq; k++) {
                            for (index_t m=0; m<numComp; m++) {
                                const Scalar wD0 = D_p[INDEX2(k, m, numEq)]*w3;
                                EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=wD0;
                                EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+=wD0;
                                EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=wD0;
                                EM_S[INDEX4(k,m,3,0,numEq,numComp,4)]+=wD0;
                                EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+=wD0;
                                EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=wD0;
                                EM_S[INDEX4(k,m,2,1,numEq,numComp,4)]+=wD0;
                                EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+=wD0;
                                EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=wD0;
                                EM_S[INDEX4(k,m,1,2,numEq,numComp,4)]+=wD0;
                                EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=wD0;
                                EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+=wD0;
                                EM_S[INDEX4(k,m,0,3,numEq,numComp,4)]+=wD0;
                                EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+=wD0;
                                EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+=wD0;
                                EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=wD0;
                            }
                        }
                    }
                    ///////////////
                    // process X //
                    ///////////////
                    if (!X.isEmpty()) {
                        const Scalar* X_p = X.getSampleDataRO(e, zero);
                        for (index_t k=0; k<numEq; k++) {
                            const Scalar wX0 = 4.*X_p[INDEX2(k, 0, numEq)]*w2;
                            const Scalar wX1 = 4.*X_p[INDEX2(k, 1, numEq)]*w1;
                            EM_F[INDEX2(k,0,numEq)]+=-wX0 - wX1;
                            EM_F[INDEX2(k,1,numEq)]+= wX0 - wX1;
                            EM_F[INDEX2(k,2,numEq)]+=-wX0 + wX1;
                            EM_F[INDEX2(k,3,numEq)]+= wX0 + wX1;
                        }
                    }
                    ///////////////
                    // process Y //
                    ///////////////
                    if (!Y.isEmpty()) {
                        const Scalar* Y_p = Y.getSampleDataRO(e, zero);
                        for (index_t k=0; k<numEq; k++) {
                            EM_F[INDEX2(k,0,numEq)]+=4.*Y_p[k]*w3;
                            EM_F[INDEX2(k,1,numEq)]+=4.*Y_p[k]*w3;
                            EM_F[INDEX2(k,2,numEq)]+=4.*Y_p[k]*w3;
                            EM_F[INDEX2(k,3,numEq)]+=4.*Y_p[k]*w3;
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

/****************************************************************************/
// PDE SYSTEM REDUCED BOUNDARY
/****************************************************************************/

template<class Scalar>
void DefaultAssembler2D<Scalar>::assemblePDEBoundarySystemReduced(
                                         AbstractSystemMatrix* mat, Data& rhs,
                                         const Data& d, const Data& y) const
{
    dim_t numEq, numComp;
    if (!mat)
        numEq=numComp=(rhs.isEmpty() ? 1 : rhs.getDataPointSize());
    else {
        numEq=mat->getRowBlockSize();
        numComp=mat->getColumnBlockSize();
    }
    const double w0 = m_dx[0]/4;
    const double w1 = m_dx[1]/4;
    const int NE0 = m_NE[0];
    const int NE1 = m_NE[1];
    const bool addEM_S = !d.isEmpty();
    const bool addEM_F = !y.isEmpty();
    const Scalar zero = static_cast<Scalar>(0);
    rhs.requireWrite();

#pragma omp parallel
    {
        vector<Scalar> EM_S(4*4*numEq*numComp, zero);
        vector<Scalar> EM_F(4*numEq, zero);

        if (domain->m_faceOffset[0] > -1) {
            if (addEM_S)
                fill(EM_S.begin(), EM_S.end(), zero);
            if (addEM_F)
                fill(EM_F.begin(), EM_F.end(), zero);

            for (index_t k1_0=0; k1_0<2; k1_0++) { // colouring
#pragma omp for
                for (index_t k1=k1_0; k1<NE1; k1+=2) {
                    ///////////////
                    // process d //
                    ///////////////
                    if (addEM_S) {
                        const Scalar* d_p = d.getSampleDataRO(k1, zero);
                        for (index_t k=0; k<numEq; k++) {
                            for (index_t m=0; m<numComp; m++) {
                                const Scalar tmp0 = d_p[INDEX2(k, m, numEq)]*w1;
                                EM_S[INDEX4(k,m,0,0,numEq,numComp,4)] = tmp0;
                                EM_S[INDEX4(k,m,2,0,numEq,numComp,4)] = tmp0;
                                EM_S[INDEX4(k,m,0,2,numEq,numComp,4)] = tmp0;
                                EM_S[INDEX4(k,m,2,2,numEq,numComp,4)] = tmp0;
                            }
                        }
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (addEM_F) {
                        const Scalar* y_p = y.getSampleDataRO(k1, zero);
                        for (index_t k=0; k<numEq; k++) {
                            EM_F[INDEX2(k,0,numEq)] = 2*w1*y_p[k];
                            EM_F[INDEX2(k,2,numEq)] = 2*w1*y_p[k];
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
                fill(EM_S.begin(), EM_S.end(), zero);
            if (addEM_F)
                fill(EM_F.begin(), EM_F.end(), zero);

            for (index_t k1_0=0; k1_0<2; k1_0++) { // colouring            
#pragma omp for
                for (index_t k1=k1_0; k1<NE1; k1+=2) {
                    const index_t e = domain->m_faceOffset[1]+k1;
                    ///////////////
                    // process d //
                    ///////////////
                    if (addEM_S) {
                        const Scalar* d_p = d.getSampleDataRO(e, zero);
                        for (index_t k=0; k<numEq; k++) {
                            for (index_t m=0; m<numComp; m++) {
                                const Scalar tmp0 = d_p[INDEX2(k, m, numEq)]*w1;
                                EM_S[INDEX4(k,m,1,1,numEq,numComp,4)] = tmp0;
                                EM_S[INDEX4(k,m,3,1,numEq,numComp,4)] = tmp0;
                                EM_S[INDEX4(k,m,1,3,numEq,numComp,4)] = tmp0;
                                EM_S[INDEX4(k,m,3,3,numEq,numComp,4)] = tmp0;
                            }
                        }
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (addEM_F) {
                        const Scalar* y_p = y.getSampleDataRO(e, zero);
                        for (index_t k=0; k<numEq; k++) {
                            EM_F[INDEX2(k,1,numEq)] = 2.*w1*y_p[k];
                            EM_F[INDEX2(k,3,numEq)] = 2.*w1*y_p[k];
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
                fill(EM_S.begin(), EM_S.end(), zero);
            if (addEM_F)
                fill(EM_F.begin(), EM_F.end(), zero);

            for (index_t k0_0=0; k0_0<2; k0_0++) { // colouring
#pragma omp for
                for (index_t k0 = k0_0; k0 < NE0; k0+=2) {
                    const index_t e = domain->m_faceOffset[2]+k0;
                    ///////////////
                    // process d //
                    ///////////////
                    if (addEM_S) {
                        const Scalar* d_p = d.getSampleDataRO(e, zero);
                        for (index_t k=0; k<numEq; k++) {
                            for (index_t m=0; m<numComp; m++) {
                                const Scalar tmp0 = d_p[INDEX2(k, m, numEq)]*w0;
                                EM_S[INDEX4(k,m,0,0,numEq,numComp,4)] = tmp0;
                                EM_S[INDEX4(k,m,1,0,numEq,numComp,4)] = tmp0;
                                EM_S[INDEX4(k,m,0,1,numEq,numComp,4)] = tmp0;
                                EM_S[INDEX4(k,m,1,1,numEq,numComp,4)] = tmp0;
                            }
                        }
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (addEM_F) {
                        const Scalar* y_p = y.getSampleDataRO(e, zero);
                        for (index_t k=0; k<numEq; k++) {
                            EM_F[INDEX2(k,0,numEq)] = 2.*w0*y_p[k];
                            EM_F[INDEX2(k,1,numEq)] = 2.*w0*y_p[k];
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
                fill(EM_S.begin(), EM_S.end(), zero);
            if (addEM_F)
                fill(EM_F.begin(), EM_F.end(), zero);

            for (index_t k0_0=0; k0_0<2; k0_0++) { // colouring
#pragma omp for
                for (index_t k0 = k0_0; k0 < NE0; k0+=2) {
                    const index_t e = domain->m_faceOffset[3]+k0;
                    ///////////////
                    // process d //
                    ///////////////
                    if (addEM_S) {
                        const Scalar* d_p = d.getSampleDataRO(e, zero);
                        for (index_t k=0; k<numEq; k++) {
                            for (index_t m=0; m<numComp; m++) {
                                const Scalar tmp0 = d_p[INDEX2(k, m, numEq)]*w0;
                                EM_S[INDEX4(k,m,2,2,numEq,numComp,4)] = tmp0;
                                EM_S[INDEX4(k,m,3,2,numEq,numComp,4)] = tmp0;
                                EM_S[INDEX4(k,m,2,3,numEq,numComp,4)] = tmp0;
                                EM_S[INDEX4(k,m,3,3,numEq,numComp,4)] = tmp0;
                            }
                        }
                    }
                    ///////////////
                    // process y //
                    ///////////////
                    if (addEM_F) {
                        const Scalar* y_p = y.getSampleDataRO(e, zero);
                        for (index_t k=0; k<numEq; k++) {
                            EM_F[INDEX2(k,2,numEq)] = 2*w0*y_p[k];
                            EM_F[INDEX2(k,3,numEq)] = 2*w0*y_p[k];
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

// instantiate our two supported versions
template class DefaultAssembler2D<escript::DataTypes::real_t>;
template class DefaultAssembler2D<escript::DataTypes::cplx_t>;

} // namespace ripley

