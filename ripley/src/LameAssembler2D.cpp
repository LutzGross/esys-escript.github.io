
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
#include <ripley/LameAssembler2D.h>

using namespace std;

namespace ripley {


void LameAssembler2D::assemblePDESingle(Paso_SystemMatrix* mat,
        escript::Data& rhs, map<string, escript::Data> coefs) const
{
    throw RipleyException("assemblePDESingle not implemented in LameAssembler2D");
}

void LameAssembler2D::assemblePDEBoundarySingle(Paso_SystemMatrix* mat,
        escript::Data& rhs, map<string, escript::Data> coefs) const 
{
    throw RipleyException("assemblePDEBoundarySingle not implemented in LameAssembler2D");
}

void LameAssembler2D::assemblePDESingleReduced(Paso_SystemMatrix* mat,
            escript::Data& rhs, map<string, escript::Data> coefs) const
{
    throw RipleyException("assemblePDESingleReduced not implemented in LameAssembler2D");
}
void LameAssembler2D::assemblePDEBoundarySingleReduced(Paso_SystemMatrix* mat,
            escript::Data& rhs, map<string, escript::Data> coefs) const
{
    throw RipleyException("assemblePDEBoundarySingleReduced not implemented in LameAssembler2D");
}

void LameAssembler2D::assemblePDEBoundarySystem(Paso_SystemMatrix* mat,
            escript::Data& rhs, map<string, escript::Data> coefs) const
{
    throw RipleyException("assemblePDEBoundarySystem not implemented in LameAssembler2D");
}

void LameAssembler2D::assemblePDESystemReduced(Paso_SystemMatrix* mat,
            escript::Data& rhs, map<string, escript::Data> coefs) const
{
    throw RipleyException("assemblePDEBoundarySystem not implemented in LameAssembler2D");
}

void LameAssembler2D::assemblePDEBoundarySystemReduced(Paso_SystemMatrix* mat,
            escript::Data& rhs, map<string, escript::Data> coefs) const
{
    throw RipleyException("assemblePDEBoundarySystemReduced not implemented in LameAssembler2D");
}

void LameAssembler2D::assemblePDESystem(Paso_SystemMatrix* mat,
            escript::Data& rhs, map<string, escript::Data> coefs) const
{
    escript::Data lambda = unpackData("lame_lambda", coefs),
		  mu = unpackData("lame_mu", coefs), B = unpackData("B", coefs),
                  C = unpackData("C", coefs), D = unpackData("D", coefs),
                  X = unpackData("X", coefs), Y = unpackData("Y", coefs);
    dim_t numEq, numComp;
    if (!mat)
        numEq=numComp=(rhs.isEmpty() ? 1 : rhs.getDataPointSize());
    else {
        numEq=mat->logical_row_block_size;
        numComp=mat->logical_col_block_size;
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

    rhs.requireWrite();
#pragma omp parallel
    {
        for (index_t k1_0=0; k1_0<2; k1_0++) { // colouring
#pragma omp for
            for (index_t k1=k1_0; k1 < m_NE[1]; k1+=2) {
                for (index_t k0=0; k0 < m_NE[0]; ++k0)  {
                    bool addEM_S=false;
                    bool addEM_F=false;
                    vector<double> EM_S(4*4*numEq*numComp, 0);
                    vector<double> EM_F(4*numEq, 0);
                    const index_t e = k0 + m_NE[0]*k1;
                    ///////////////
                    // process A //
                    ///////////////
                    if (!mu.isEmpty() && !lambda.isEmpty()) {
                        addEM_S = true;
                        if (true) {
                            for (index_t k=0; k<numEq; k++) {
                                for (index_t m=0; m<numComp; m++) {
                                    const double A_00_0 = 1;//A_p[INDEX5(k,0,m,0,0,numEq,2,numComp,2)];
                                    const double A_01_0 = 1;//A_p[INDEX5(k,0,m,1,0,numEq,2,numComp,2)];
                                    const double A_10_0 = 1;//A_p[INDEX5(k,1,m,0,0,numEq,2,numComp,2)];
                                    const double A_11_0 = 1;//A_p[INDEX5(k,1,m,1,0,numEq,2,numComp,2)];
                                    const double A_00_1 = 1;//A_p[INDEX5(k,0,m,0,1,numEq,2,numComp,2)];
                                    const double A_01_1 = 1;//A_p[INDEX5(k,0,m,1,1,numEq,2,numComp,2)];
                                    const double A_10_1 = 1;//A_p[INDEX5(k,1,m,0,1,numEq,2,numComp,2)];
                                    const double A_11_1 = 1;//A_p[INDEX5(k,1,m,1,1,numEq,2,numComp,2)];
                                    const double A_00_2 = 1;//A_p[INDEX5(k,0,m,0,2,numEq,2,numComp,2)];
                                    const double A_01_2 = 1;//A_p[INDEX5(k,0,m,1,2,numEq,2,numComp,2)];
                                    const double A_10_2 = 1;//A_p[INDEX5(k,1,m,0,2,numEq,2,numComp,2)];
                                    const double A_11_2 = 1;//A_p[INDEX5(k,1,m,1,2,numEq,2,numComp,2)];
                                    const double A_00_3 = 1;//A_p[INDEX5(k,0,m,0,3,numEq,2,numComp,2)];
                                    const double A_01_3 = 1;//A_p[INDEX5(k,0,m,1,3,numEq,2,numComp,2)];
                                    const double A_10_3 = 1;//A_p[INDEX5(k,1,m,0,3,numEq,2,numComp,2)];
                                    const double A_11_3 = 1;//A_p[INDEX5(k,1,m,1,3,numEq,2,numComp,2)];
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
                                    const double A_00 = 1;//A_p[INDEX4(k,0,m,0, numEq,2, numComp)];
                                    const double A_01 = 1;//A_p[INDEX4(k,0,m,1, numEq,2, numComp)];
                                    const double A_10 = 1;//A_p[INDEX4(k,1,m,0, numEq,2, numComp)];
                                    const double A_11 = 1;//A_p[INDEX4(k,1,m,1, numEq,2, numComp)];
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
                        addEM_S=true;
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
                        addEM_S=true;
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
                        addEM_S=true;
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
                        addEM_F=true;
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
                        addEM_F=true;
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

}
