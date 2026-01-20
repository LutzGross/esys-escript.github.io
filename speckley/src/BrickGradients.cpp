
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

#include <speckley/Brick.h>

#include <escript/index.h>

namespace speckley {

template<typename Scalar>
void Brick::gradient_order2(escript::Data& out, const escript::Data& in) const
{
    const real_t lagrange_deriv_0[3] = {-1.50000000000000, -0.500000000000000, 0.500000000000000};
    const real_t lagrange_deriv_1[3] = {2.00000000000000, 0, -2.00000000000000};
    const real_t lagrange_deriv_2[3] = {-0.500000000000000, 0.500000000000000, 1.50000000000000};
    const real_t inv_jac[3] = {2/m_dx[0], 2/m_dx[1], 2/m_dx[2]}; //inverse jacobi
    const int numComp = in.getDataPointSize();
    const Scalar zero = static_cast<Scalar>(0);
    out.requireWrite();
    if (!in.actsExpanded()) {
#pragma omp parallel for
        for (int ei = 0; ei < m_NE[2]; ++ei) {
            for (int ej = 0; ej < m_NE[1]; ++ej) {
                for (int ek = 0; ek < m_NE[0]; ++ek) {
                    const Scalar* e = in.getSampleDataRO(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]), zero);
                    Scalar* grad = out.getSampleDataRW(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]), zero);
                    for (int comp = 0; comp < numComp; ++comp) {
                        const Scalar a = (lagrange_deriv_0[0] * e[comp] + lagrange_deriv_1[0] * e[comp] + lagrange_deriv_2[0] * e[comp]) * inv_jac[0];
                        const Scalar b = (lagrange_deriv_0[0] * e[comp] + lagrange_deriv_1[0] * e[comp] + lagrange_deriv_2[0] * e[comp]) * inv_jac[1];
                        const Scalar c = (lagrange_deriv_0[0] * e[comp] + lagrange_deriv_1[0] * e[comp] + lagrange_deriv_2[0] * e[comp]) * inv_jac[2];
                        for (int k = 0; k < 3; ++k) {
                            for (int j = 0; j < 3; ++j) {
                                for (int i = 0; i < 3; ++i) {
                                    const index_t ind = INDEX5(comp,0,i,j,k,numComp,3,3,3);
                                    grad[ind + 0] = a;
                                    grad[ind + 1] = b;
                                    grad[ind + 2] = c;
                                }
                            }
                        }
                    }
                }
            }
        }
    } else {
#pragma omp parallel for
        for (int ei = 0; ei < m_NE[2]; ++ei) {
            for (int ej = 0; ej < m_NE[1]; ++ej) {
                for (int ek = 0; ek < m_NE[0]; ++ek) {
                    const Scalar* e = in.getSampleDataRO(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]), zero);
                    Scalar* grad = out.getSampleDataRW(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]), zero);
                    for (int k = 0; k < 3; ++k) {
                        for (int j = 0; j < 3; ++j) {
                            for (int i = 0; i < 3; ++i) {
                                for (int comp = 0; comp < numComp; ++comp) {
                                    grad[INDEX5(comp,0,i,j,k,numComp,3,3,3)] = (lagrange_deriv_0[i] * e[INDEX4(comp,0,j,k,numComp,3,3)] + lagrange_deriv_1[i] * e[INDEX4(comp,1,j,k,numComp,3,3)] + lagrange_deriv_2[i] * e[INDEX4(comp,2,j,k,numComp,3,3)]) * inv_jac[0];
                                    grad[INDEX5(comp,1,i,j,k,numComp,3,3,3)] = (lagrange_deriv_0[j] * e[INDEX4(comp,i,0,k,numComp,3,3)] + lagrange_deriv_1[j] * e[INDEX4(comp,i,1,k,numComp,3,3)] + lagrange_deriv_2[j] * e[INDEX4(comp,i,2,k,numComp,3,3)]) * inv_jac[1];
                                    grad[INDEX5(comp,2,i,j,k,numComp,3,3,3)] = (lagrange_deriv_0[k] * e[INDEX4(comp,i,j,0,numComp,3,3)] + lagrange_deriv_1[k] * e[INDEX4(comp,i,j,1,numComp,3,3)] + lagrange_deriv_2[k] * e[INDEX4(comp,i,j,2,numComp,3,3)]) * inv_jac[2];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

template<typename Scalar>
void Brick::gradient_order3(escript::Data& out, const escript::Data& in) const {
    const real_t lagrange_deriv_0[4] = {-3.00000000000000, -0.809016994374948, 0.309016994374948, -0.500000000000000};
    const real_t lagrange_deriv_1[4] = {4.04508497187474, 4.44089209850063e-16, -1.11803398874990, 1.54508497187474};
    const real_t lagrange_deriv_2[4] = {-1.54508497187474, 1.11803398874989, 2.22044604925031e-16, -4.04508497187474};
    const real_t lagrange_deriv_3[4] = {0.500000000000000, -0.309016994374947, 0.809016994374948, 3.00000000000000};
    const real_t inv_jac[3] = {2/m_dx[0], 2/m_dx[1], 2/m_dx[2]}; //inverse jacobi
    const int numComp = in.getDataPointSize();
    const Scalar zero = static_cast<Scalar>(0);
    out.requireWrite();
    if (!in.actsExpanded()) {
#pragma omp parallel for
        for (int ei = 0; ei < m_NE[2]; ++ei) {
            for (int ej = 0; ej < m_NE[1]; ++ej) {
                for (int ek = 0; ek < m_NE[0]; ++ek) {
                    const Scalar* e = in.getSampleDataRO(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]), zero);
                    Scalar* grad = out.getSampleDataRW(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]), zero);
                    for (int comp = 0; comp < numComp; ++comp) {
                        const Scalar a = (lagrange_deriv_0[0] * e[comp] + lagrange_deriv_1[0] * e[comp] + lagrange_deriv_2[0] * e[comp] + lagrange_deriv_3[0] * e[comp]) * inv_jac[0];
                        const Scalar b = (lagrange_deriv_0[0] * e[comp] + lagrange_deriv_1[0] * e[comp] + lagrange_deriv_2[0] * e[comp] + lagrange_deriv_3[0] * e[comp]) * inv_jac[1];
                        const Scalar c = (lagrange_deriv_0[0] * e[comp] + lagrange_deriv_1[0] * e[comp] + lagrange_deriv_2[0] * e[comp] + lagrange_deriv_3[0] * e[comp]) * inv_jac[2];
                        for (int k = 0; k < 4; ++k) {
                            for (int j = 0; j < 4; ++j) {
                                for (int i = 0; i < 4; ++i) {
                                    const index_t ind = INDEX5(comp,0,i,j,k,numComp,3,4,4);
                                    grad[ind + 0] = a;
                                    grad[ind + 1] = b;
                                    grad[ind + 2] = c;
                                }
                            }
                        }
                    }
                }
            }
        }
    } else {
#pragma omp parallel for
        for (int ei = 0; ei < m_NE[2]; ++ei) {
            for (int ej = 0; ej < m_NE[1]; ++ej) {
                for (int ek = 0; ek < m_NE[0]; ++ek) {
                    const Scalar* e = in.getSampleDataRO(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]), zero);
                    Scalar* grad = out.getSampleDataRW(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]), zero);
                    for (int k = 0; k < 4; ++k) {
                        for (int j = 0; j < 4; ++j) {
                            for (int i = 0; i < 4; ++i) {
                                for (int comp = 0; comp < numComp; ++comp) {
                                    grad[INDEX5(comp,0,i,j,k,numComp,3,4,4)] = (lagrange_deriv_0[i] * e[INDEX4(comp,0,j,k,numComp,4,4)] + lagrange_deriv_1[i] * e[INDEX4(comp,1,j,k,numComp,4,4)] + lagrange_deriv_2[i] * e[INDEX4(comp,2,j,k,numComp,4,4)] + lagrange_deriv_3[i] * e[INDEX4(comp,3,j,k,numComp,4,4)]) * inv_jac[0];
                                    grad[INDEX5(comp,1,i,j,k,numComp,3,4,4)] = (lagrange_deriv_0[j] * e[INDEX4(comp,i,0,k,numComp,4,4)] + lagrange_deriv_1[j] * e[INDEX4(comp,i,1,k,numComp,4,4)] + lagrange_deriv_2[j] * e[INDEX4(comp,i,2,k,numComp,4,4)] + lagrange_deriv_3[j] * e[INDEX4(comp,i,3,k,numComp,4,4)]) * inv_jac[1];
                                    grad[INDEX5(comp,2,i,j,k,numComp,3,4,4)] = (lagrange_deriv_0[k] * e[INDEX4(comp,i,j,0,numComp,4,4)] + lagrange_deriv_1[k] * e[INDEX4(comp,i,j,1,numComp,4,4)] + lagrange_deriv_2[k] * e[INDEX4(comp,i,j,2,numComp,4,4)] + lagrange_deriv_3[k] * e[INDEX4(comp,i,j,3,numComp,4,4)]) * inv_jac[2];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

template<typename Scalar>
void Brick::gradient_order4(escript::Data& out, const escript::Data& in) const {

    const real_t lagrange_deriv_0[5] = {-4.99999999999999, -1.24099025303098, 0.374999999999999, -0.259009746969017, 0.499999999999999};
    const real_t lagrange_deriv_1[5] = {6.75650248872424, -6.66133814775094e-15, -1.33658457769545, 0.763762615825974, -1.41016417794243};
    const real_t lagrange_deriv_2[5] = {-2.66666666666667, 1.74574312188794, 1.44328993201270e-15, -1.74574312188794, 2.66666666666667};
    const real_t lagrange_deriv_3[5] = {1.41016417794243, -0.763762615825974, 1.33658457769545, 1.66533453693773e-15, -6.75650248872424};
    const real_t lagrange_deriv_4[5] = {-0.500000000000001, 0.259009746969017, -0.375000000000000, 1.24099025303098, 5.00000000000000};

    const real_t inv_jac[3] = {2/m_dx[0], 2/m_dx[1], 2/m_dx[2]}; //inverse jacobi
    const int numComp = in.getDataPointSize();
    const Scalar zero = static_cast<Scalar>(0);
    out.requireWrite();
    if (!in.actsExpanded()) {
#pragma omp parallel for
        for (int ei = 0; ei < m_NE[2]; ++ei) {
            for (int ej = 0; ej < m_NE[1]; ++ej) {
                for (int ek = 0; ek < m_NE[0]; ++ek) {
                    const Scalar* e = in.getSampleDataRO(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]), zero);
                    Scalar* grad = out.getSampleDataRW(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]), zero);
                    for (int comp = 0; comp < numComp; ++comp) {
                        const Scalar a = (lagrange_deriv_0[0] * e[comp] + lagrange_deriv_1[0] * e[comp] + lagrange_deriv_2[0] * e[comp] + lagrange_deriv_3[0] * e[comp] + lagrange_deriv_4[0] * e[comp]) * inv_jac[0];
                        const Scalar b = (lagrange_deriv_0[0] * e[comp] + lagrange_deriv_1[0] * e[comp] + lagrange_deriv_2[0] * e[comp] + lagrange_deriv_3[0] * e[comp] + lagrange_deriv_4[0] * e[comp]) * inv_jac[1];
                        const Scalar c = (lagrange_deriv_0[0] * e[comp] + lagrange_deriv_1[0] * e[comp] + lagrange_deriv_2[0] * e[comp] + lagrange_deriv_3[0] * e[comp] + lagrange_deriv_4[0] * e[comp]) * inv_jac[2];
                        for (int k = 0; k < 5; ++k) {
                            for (int j = 0; j < 5; ++j) {
                                for (int i = 0; i < 5; ++i) {
                                    const index_t ind = INDEX5(comp,0,i,j,k,numComp,3,5,5);
                                    grad[ind + 0] = a;
                                    grad[ind + 1] = b;
                                    grad[ind + 2] = c;
                                }
                            }
                        }
                    }
                }
            }
        }
    } else {
#pragma omp parallel for
        for (int ei = 0; ei < m_NE[2]; ++ei) {
            for (int ej = 0; ej < m_NE[1]; ++ej) {
                for (int ek = 0; ek < m_NE[0]; ++ek) {
                    const Scalar* e = in.getSampleDataRO(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]), zero);
                    Scalar* grad = out.getSampleDataRW(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]), zero);
                    for (int k = 0; k < 5; ++k) {
                        for (int j = 0; j < 5; ++j) {
                            for (int i = 0; i < 5; ++i) {
                                for (int comp = 0; comp < numComp; ++comp) {
                                    grad[INDEX5(comp,0,i,j,k,numComp,3,5,5)] = (lagrange_deriv_0[i] * e[INDEX4(comp,0,j,k,numComp,5,5)] + lagrange_deriv_1[i] * e[INDEX4(comp,1,j,k,numComp,5,5)] + lagrange_deriv_2[i] * e[INDEX4(comp,2,j,k,numComp,5,5)] + lagrange_deriv_3[i] * e[INDEX4(comp,3,j,k,numComp,5,5)] + lagrange_deriv_4[i] * e[INDEX4(comp,4,j,k,numComp,5,5)]) * inv_jac[0];
                                    grad[INDEX5(comp,1,i,j,k,numComp,3,5,5)] = (lagrange_deriv_0[j] * e[INDEX4(comp,i,0,k,numComp,5,5)] + lagrange_deriv_1[j] * e[INDEX4(comp,i,1,k,numComp,5,5)] + lagrange_deriv_2[j] * e[INDEX4(comp,i,2,k,numComp,5,5)] + lagrange_deriv_3[j] * e[INDEX4(comp,i,3,k,numComp,5,5)] + lagrange_deriv_4[j] * e[INDEX4(comp,i,4,k,numComp,5,5)]) * inv_jac[1];
                                    grad[INDEX5(comp,2,i,j,k,numComp,3,5,5)] = (lagrange_deriv_0[k] * e[INDEX4(comp,i,j,0,numComp,5,5)] + lagrange_deriv_1[k] * e[INDEX4(comp,i,j,1,numComp,5,5)] + lagrange_deriv_2[k] * e[INDEX4(comp,i,j,2,numComp,5,5)] + lagrange_deriv_3[k] * e[INDEX4(comp,i,j,3,numComp,5,5)] + lagrange_deriv_4[k] * e[INDEX4(comp,i,j,4,numComp,5,5)]) * inv_jac[2];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

template<typename Scalar>
void Brick::gradient_order5(escript::Data& out, const escript::Data& in) const {

    const real_t lagrange_deriv_0[6] = {-7.50000000000002, -1.78636494833911, 0.484951047853572, -0.269700610832040, 0.237781177984232, -0.500000000000002};
    const real_t lagrange_deriv_1[6] = {10.1414159363197, 2.13162820728030e-14, -1.72125695283023, 0.786356672223240, -0.653547507429800, 1.34991331419049};
    const real_t lagrange_deriv_2[6] = {-4.03618727030532, 2.52342677742945, -4.66293670342566e-15, -1.75296196636786, 1.15282815853593, -2.24468464817616};
    const real_t lagrange_deriv_3[6] = {2.24468464817616, -1.15282815853593, 1.75296196636787, -1.77635683940025e-15, -2.52342677742946, 4.03618727030535};
    const real_t lagrange_deriv_4[6] = {-1.34991331419048, 0.653547507429800, -0.786356672223242, 1.72125695283023, 2.22044604925031e-15, -10.1414159363197};
    const real_t lagrange_deriv_5[6] = {0.499999999999998, -0.237781177984231, 0.269700610832039, -0.484951047853569, 1.78636494833909, 7.50000000000000};

    const real_t inv_jac[3] = {2/m_dx[0], 2/m_dx[1], 2/m_dx[2]}; //inverse jacobi
    const int numComp = in.getDataPointSize();
    const Scalar zero = static_cast<Scalar>(0);
    out.requireWrite();
    if (!in.actsExpanded()) {
#pragma omp parallel for
        for (int ei = 0; ei < m_NE[2]; ++ei) {
            for (int ej = 0; ej < m_NE[1]; ++ej) {
                for (int ek = 0; ek < m_NE[0]; ++ek) {
                    const Scalar* e = in.getSampleDataRO(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]), zero);
                    Scalar* grad = out.getSampleDataRW(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]), zero);
                    for (int comp = 0; comp < numComp; ++comp) {
                        const Scalar a = (lagrange_deriv_0[0] * e[comp] + lagrange_deriv_1[0] * e[comp] + lagrange_deriv_2[0] * e[comp] + lagrange_deriv_3[0] * e[comp] + lagrange_deriv_4[0] * e[comp] + lagrange_deriv_5[0] * e[comp]) * inv_jac[0];
                        const Scalar b = (lagrange_deriv_0[0] * e[comp] + lagrange_deriv_1[0] * e[comp] + lagrange_deriv_2[0] * e[comp] + lagrange_deriv_3[0] * e[comp] + lagrange_deriv_4[0] * e[comp] + lagrange_deriv_5[0] * e[comp]) * inv_jac[1];
                        const Scalar c = (lagrange_deriv_0[0] * e[comp] + lagrange_deriv_1[0] * e[comp] + lagrange_deriv_2[0] * e[comp] + lagrange_deriv_3[0] * e[comp] + lagrange_deriv_4[0] * e[comp] + lagrange_deriv_5[0] * e[comp]) * inv_jac[2];
                        for (int k = 0; k < 6; ++k) {
                            for (int j = 0; j < 6; ++j) {
                                for (int i = 0; i < 6; ++i) {
                                    const index_t ind = INDEX5(comp,0,i,j,k,numComp,3,6,6);
                                    grad[ind + 0] = a;
                                    grad[ind + 1] = b;
                                    grad[ind + 2] = c;
                                }
                            }
                        }
                    }
                }
            }
        }
    } else {
#pragma omp parallel for
        for (int ei = 0; ei < m_NE[2]; ++ei) {
            for (int ej = 0; ej < m_NE[1]; ++ej) {
                for (int ek = 0; ek < m_NE[0]; ++ek) {
                    const Scalar* e = in.getSampleDataRO(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]), zero);
                    Scalar* grad = out.getSampleDataRW(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]), zero);
                    for (int k = 0; k < 6; ++k) {
                        for (int j = 0; j < 6; ++j) {
                            for (int i = 0; i < 6; ++i) {
                                for (int comp = 0; comp < numComp; ++comp) {
                                    grad[INDEX5(comp,0,i,j,k,numComp,3,6,6)] = (lagrange_deriv_0[i] * e[INDEX4(comp,0,j,k,numComp,6,6)] + lagrange_deriv_1[i] * e[INDEX4(comp,1,j,k,numComp,6,6)] + lagrange_deriv_2[i] * e[INDEX4(comp,2,j,k,numComp,6,6)] + lagrange_deriv_3[i] * e[INDEX4(comp,3,j,k,numComp,6,6)] + lagrange_deriv_4[i] * e[INDEX4(comp,4,j,k,numComp,6,6)] + lagrange_deriv_5[i] * e[INDEX4(comp,5,j,k,numComp,6,6)]) * inv_jac[0];
                                    grad[INDEX5(comp,1,i,j,k,numComp,3,6,6)] = (lagrange_deriv_0[j] * e[INDEX4(comp,i,0,k,numComp,6,6)] + lagrange_deriv_1[j] * e[INDEX4(comp,i,1,k,numComp,6,6)] + lagrange_deriv_2[j] * e[INDEX4(comp,i,2,k,numComp,6,6)] + lagrange_deriv_3[j] * e[INDEX4(comp,i,3,k,numComp,6,6)] + lagrange_deriv_4[j] * e[INDEX4(comp,i,4,k,numComp,6,6)] + lagrange_deriv_5[j] * e[INDEX4(comp,i,5,k,numComp,6,6)]) * inv_jac[1];
                                    grad[INDEX5(comp,2,i,j,k,numComp,3,6,6)] = (lagrange_deriv_0[k] * e[INDEX4(comp,i,j,0,numComp,6,6)] + lagrange_deriv_1[k] * e[INDEX4(comp,i,j,1,numComp,6,6)] + lagrange_deriv_2[k] * e[INDEX4(comp,i,j,2,numComp,6,6)] + lagrange_deriv_3[k] * e[INDEX4(comp,i,j,3,numComp,6,6)] + lagrange_deriv_4[k] * e[INDEX4(comp,i,j,4,numComp,6,6)] + lagrange_deriv_5[k] * e[INDEX4(comp,i,j,5,numComp,6,6)]) * inv_jac[2];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

template<typename Scalar>
void Brick::gradient_order6(escript::Data& out, const escript::Data& in) const {

    const real_t lagrange_deriv_0[7] = {-10.5000000000000, -2.44292601424426, 0.625256665515336, -0.312499999999997, 0.226099400942572, -0.226611870395444, 0.500000000000001};
    const real_t lagrange_deriv_1[7] = {14.2015766029198, -4.17443857259059e-14, -2.21580428316997, 0.907544471268819, -0.616390835517577, 0.602247179635785, -1.31737343570244};
    const real_t lagrange_deriv_2[7] = {-5.66898522554555, 3.45582821429430, 3.10862446895044e-15, -2.00696924058875, 1.06644190400637, -0.961339797288714, 2.04996481307676};
    const real_t lagrange_deriv_3[7] = {3.20000000000003, -1.59860668809837, 2.26669808708600, 1.33226762955019e-15, -2.26669808708599, 1.59860668809837, -3.20000000000003};
    const real_t lagrange_deriv_4[7] = {-2.04996481307676, 0.961339797288717, -1.06644190400638, 2.00696924058876, -1.77635683940025e-14, -3.45582821429431, 5.66898522554558};
    const real_t lagrange_deriv_5[7] = {1.31737343570245, -0.602247179635788, 0.616390835517580, -0.907544471268822, 2.21580428316998, 6.76125821996720e-14, -14.2015766029198};
    const real_t lagrange_deriv_6[7] = {-0.500000000000000, 0.226611870395444, -0.226099400942572, 0.312499999999997, -0.625256665515335, 2.44292601424425, 10.5000000000000};

    const real_t inv_jac[3] = {2/m_dx[0], 2/m_dx[1], 2/m_dx[2]}; //inverse jacobi
    const int numComp = in.getDataPointSize();
    const Scalar zero = static_cast<Scalar>(0);
    out.requireWrite();
    if (!in.actsExpanded()) {
#pragma omp parallel for
        for (int ei = 0; ei < m_NE[2]; ++ei) {
            for (int ej = 0; ej < m_NE[1]; ++ej) {
                for (int ek = 0; ek < m_NE[0]; ++ek) {
                    const Scalar* e = in.getSampleDataRO(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]), zero);
                    Scalar* grad = out.getSampleDataRW(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]), zero);
                    for (int comp = 0; comp < numComp; ++comp) {
                        const Scalar a = (lagrange_deriv_0[0] * e[comp] + lagrange_deriv_1[0] * e[comp] + lagrange_deriv_2[0] * e[comp] + lagrange_deriv_3[0] * e[comp] + lagrange_deriv_4[0] * e[comp] + lagrange_deriv_5[0] * e[comp] + lagrange_deriv_6[0] * e[comp]) * inv_jac[0];
                        const Scalar b = (lagrange_deriv_0[0] * e[comp] + lagrange_deriv_1[0] * e[comp] + lagrange_deriv_2[0] * e[comp] + lagrange_deriv_3[0] * e[comp] + lagrange_deriv_4[0] * e[comp] + lagrange_deriv_5[0] * e[comp] + lagrange_deriv_6[0] * e[comp]) * inv_jac[1];
                        const Scalar c = (lagrange_deriv_0[0] * e[comp] + lagrange_deriv_1[0] * e[comp] + lagrange_deriv_2[0] * e[comp] + lagrange_deriv_3[0] * e[comp] + lagrange_deriv_4[0] * e[comp] + lagrange_deriv_5[0] * e[comp] + lagrange_deriv_6[0] * e[comp]) * inv_jac[2];
                        for (int k = 0; k < 7; ++k) {
                            for (int j = 0; j < 7; ++j) {
                                for (int i = 0; i < 7; ++i) {
                                    const index_t ind = INDEX5(comp,0,i,j,k,numComp,3,7,7);
                                    grad[ind + 0] = a;
                                    grad[ind + 1] = b;
                                    grad[ind + 2] = c;
                                }
                            }
                        }
                    }
                }
            }
        }
    } else {
#pragma omp parallel for
        for (int ei = 0; ei < m_NE[2]; ++ei) {
            for (int ej = 0; ej < m_NE[1]; ++ej) {
                for (int ek = 0; ek < m_NE[0]; ++ek) {
                    const Scalar* e = in.getSampleDataRO(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]), zero);
                    Scalar* grad = out.getSampleDataRW(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]), zero);
                    for (int k = 0; k < 7; ++k) {
                        for (int j = 0; j < 7; ++j) {
                            for (int i = 0; i < 7; ++i) {
                                for (int comp = 0; comp < numComp; ++comp) {
                                    grad[INDEX5(comp,0,i,j,k,numComp,3,7,7)] = (lagrange_deriv_0[i] * e[INDEX4(comp,0,j,k,numComp,7,7)] + lagrange_deriv_1[i] * e[INDEX4(comp,1,j,k,numComp,7,7)] + lagrange_deriv_2[i] * e[INDEX4(comp,2,j,k,numComp,7,7)] + lagrange_deriv_3[i] * e[INDEX4(comp,3,j,k,numComp,7,7)] + lagrange_deriv_4[i] * e[INDEX4(comp,4,j,k,numComp,7,7)] + lagrange_deriv_5[i] * e[INDEX4(comp,5,j,k,numComp,7,7)] + lagrange_deriv_6[i] * e[INDEX4(comp,6,j,k,numComp,7,7)]) * inv_jac[0];
                                    grad[INDEX5(comp,1,i,j,k,numComp,3,7,7)] = (lagrange_deriv_0[j] * e[INDEX4(comp,i,0,k,numComp,7,7)] + lagrange_deriv_1[j] * e[INDEX4(comp,i,1,k,numComp,7,7)] + lagrange_deriv_2[j] * e[INDEX4(comp,i,2,k,numComp,7,7)] + lagrange_deriv_3[j] * e[INDEX4(comp,i,3,k,numComp,7,7)] + lagrange_deriv_4[j] * e[INDEX4(comp,i,4,k,numComp,7,7)] + lagrange_deriv_5[j] * e[INDEX4(comp,i,5,k,numComp,7,7)] + lagrange_deriv_6[j] * e[INDEX4(comp,i,6,k,numComp,7,7)]) * inv_jac[1];
                                    grad[INDEX5(comp,2,i,j,k,numComp,3,7,7)] = (lagrange_deriv_0[k] * e[INDEX4(comp,i,j,0,numComp,7,7)] + lagrange_deriv_1[k] * e[INDEX4(comp,i,j,1,numComp,7,7)] + lagrange_deriv_2[k] * e[INDEX4(comp,i,j,2,numComp,7,7)] + lagrange_deriv_3[k] * e[INDEX4(comp,i,j,3,numComp,7,7)] + lagrange_deriv_4[k] * e[INDEX4(comp,i,j,4,numComp,7,7)] + lagrange_deriv_5[k] * e[INDEX4(comp,i,j,5,numComp,7,7)] + lagrange_deriv_6[k] * e[INDEX4(comp,i,j,6,numComp,7,7)]) * inv_jac[2];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

template<typename Scalar>
void Brick::gradient_order7(escript::Data& out, const escript::Data& in) const {

    const real_t lagrange_deriv_0[8] = {-13.9999999999999, -3.20991570300295, 0.792476681320508, -0.372150435728592, 0.243330712723790, -0.203284568900591, 0.219957514771299, -0.499999999999980};
    const real_t lagrange_deriv_1[8] = {18.9375986071174, -6.23945339839338e-14, -2.80647579473643, 1.07894468879045, -0.661157350900312, 0.537039586157660, -0.573565414940254, 1.29768738832019};
    const real_t lagrange_deriv_2[8] = {-7.56928981934855, 4.54358506456659, 1.17683640610267e-14, -2.37818723351551, 1.13535801688112, -0.845022556506511, 0.869448098331479, -1.94165942554406};
    const real_t lagrange_deriv_3[8] = {4.29790816426521, -2.11206121431455, 2.87551740597250, 3.77475828372553e-15, -2.38892435915824, 1.37278583180603, -1.29423205091348, 2.81018898925786};
    const real_t lagrange_deriv_4[8] = {-2.81018898925796, 1.29423205091350, -1.37278583180602, 2.38892435915823, 7.77156117237610e-15, -2.87551740597249, 2.11206121431450, -4.29790816426503};
    const real_t lagrange_deriv_5[8] = {1.94165942554413, -0.869448098331493, 0.845022556506508, -1.13535801688111, 2.37818723351550, -2.44249065417534e-14, -4.54358506456649, 7.56928981934828};
    const real_t lagrange_deriv_6[8] = {-1.29768738832026, 0.573565414940274, -0.537039586157669, 0.661157350900321, -1.07894468879047, 2.80647579473648, -1.71085368094737e-13, -18.9375986071174};
    const real_t lagrange_deriv_7[8] = {0.500000000000021, -0.219957514771312, 0.203284568900599, -0.243330712723800, 0.372150435728609, -0.792476681320546, 3.20991570300311, 14.0000000000002};

    const real_t inv_jac[3] = {2/m_dx[0], 2/m_dx[1], 2/m_dx[2]}; //inverse jacobi
    const int numComp = in.getDataPointSize();
    const Scalar zero = static_cast<Scalar>(0);
    out.requireWrite();
    if (!in.actsExpanded()) {
#pragma omp parallel for
        for (int ei = 0; ei < m_NE[2]; ++ei) {
            for (int ej = 0; ej < m_NE[1]; ++ej) {
                for (int ek = 0; ek < m_NE[0]; ++ek) {
                    const Scalar* e = in.getSampleDataRO(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]), zero);
                    Scalar* grad = out.getSampleDataRW(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]), zero);
                    for (int comp = 0; comp < numComp; ++comp) {
                        const Scalar a = (lagrange_deriv_0[0] * e[comp] + lagrange_deriv_1[0] * e[comp] + lagrange_deriv_2[0] * e[comp] + lagrange_deriv_3[0] * e[comp] + lagrange_deriv_4[0] * e[comp] + lagrange_deriv_5[0] * e[comp] + lagrange_deriv_6[0] * e[comp] + lagrange_deriv_7[0] * e[comp]) * inv_jac[0];
                        const Scalar b = (lagrange_deriv_0[0] * e[comp] + lagrange_deriv_1[0] * e[comp] + lagrange_deriv_2[0] * e[comp] + lagrange_deriv_3[0] * e[comp] + lagrange_deriv_4[0] * e[comp] + lagrange_deriv_5[0] * e[comp] + lagrange_deriv_6[0] * e[comp] + lagrange_deriv_7[0] * e[comp]) * inv_jac[1];
                        const Scalar c = (lagrange_deriv_0[0] * e[comp] + lagrange_deriv_1[0] * e[comp] + lagrange_deriv_2[0] * e[comp] + lagrange_deriv_3[0] * e[comp] + lagrange_deriv_4[0] * e[comp] + lagrange_deriv_5[0] * e[comp] + lagrange_deriv_6[0] * e[comp] + lagrange_deriv_7[0] * e[comp]) * inv_jac[2];
                        for (int k = 0; k < 8; ++k) {
                            for (int j = 0; j < 8; ++j) {
                                for (int i = 0; i < 8; ++i) {
                                    const index_t ind = INDEX5(comp,0,i,j,k,numComp,3,8,8);
                                    grad[ind + 0] = a;
                                    grad[ind + 1] = b;
                                    grad[ind + 2] = c;
                                }
                            }
                        }
                    }
                }
            }
        }
    } else {
#pragma omp parallel for
        for (int ei = 0; ei < m_NE[2]; ++ei) {
            for (int ej = 0; ej < m_NE[1]; ++ej) {
                for (int ek = 0; ek < m_NE[0]; ++ek) {
                    const Scalar* e = in.getSampleDataRO(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]), zero);
                    Scalar* grad = out.getSampleDataRW(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]), zero);
                    for (int k = 0; k < 8; ++k) {
                        for (int j = 0; j < 8; ++j) {
                            for (int i = 0; i < 8; ++i) {
                                for (int comp = 0; comp < numComp; ++comp) {
                                    grad[INDEX5(comp,0,i,j,k,numComp,3,8,8)] = (lagrange_deriv_0[i] * e[INDEX4(comp,0,j,k,numComp,8,8)] + lagrange_deriv_1[i] * e[INDEX4(comp,1,j,k,numComp,8,8)] + lagrange_deriv_2[i] * e[INDEX4(comp,2,j,k,numComp,8,8)] + lagrange_deriv_3[i] * e[INDEX4(comp,3,j,k,numComp,8,8)] + lagrange_deriv_4[i] * e[INDEX4(comp,4,j,k,numComp,8,8)] + lagrange_deriv_5[i] * e[INDEX4(comp,5,j,k,numComp,8,8)] + lagrange_deriv_6[i] * e[INDEX4(comp,6,j,k,numComp,8,8)] + lagrange_deriv_7[i] * e[INDEX4(comp,7,j,k,numComp,8,8)]) * inv_jac[0];
                                    grad[INDEX5(comp,1,i,j,k,numComp,3,8,8)] = (lagrange_deriv_0[j] * e[INDEX4(comp,i,0,k,numComp,8,8)] + lagrange_deriv_1[j] * e[INDEX4(comp,i,1,k,numComp,8,8)] + lagrange_deriv_2[j] * e[INDEX4(comp,i,2,k,numComp,8,8)] + lagrange_deriv_3[j] * e[INDEX4(comp,i,3,k,numComp,8,8)] + lagrange_deriv_4[j] * e[INDEX4(comp,i,4,k,numComp,8,8)] + lagrange_deriv_5[j] * e[INDEX4(comp,i,5,k,numComp,8,8)] + lagrange_deriv_6[j] * e[INDEX4(comp,i,6,k,numComp,8,8)] + lagrange_deriv_7[j] * e[INDEX4(comp,i,7,k,numComp,8,8)]) * inv_jac[1];
                                    grad[INDEX5(comp,2,i,j,k,numComp,3,8,8)] = (lagrange_deriv_0[k] * e[INDEX4(comp,i,j,0,numComp,8,8)] + lagrange_deriv_1[k] * e[INDEX4(comp,i,j,1,numComp,8,8)] + lagrange_deriv_2[k] * e[INDEX4(comp,i,j,2,numComp,8,8)] + lagrange_deriv_3[k] * e[INDEX4(comp,i,j,3,numComp,8,8)] + lagrange_deriv_4[k] * e[INDEX4(comp,i,j,4,numComp,8,8)] + lagrange_deriv_5[k] * e[INDEX4(comp,i,j,5,numComp,8,8)] + lagrange_deriv_6[k] * e[INDEX4(comp,i,j,6,numComp,8,8)] + lagrange_deriv_7[k] * e[INDEX4(comp,i,j,7,numComp,8,8)]) * inv_jac[2];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

template<typename Scalar>
void Brick::gradient_order8(escript::Data& out, const escript::Data& in) const {

    const real_t lagrange_deriv_0[9] = {-18.0000000000010, -4.08701370203454, 0.985360090074639, -0.444613449281139, 0.273437500000029, -0.207734512035617, 0.189655591978376, -0.215654018702531, 0.500000000000095};
    const real_t lagrange_deriv_1[9] = {24.3497451715930, 1.34892097491957e-12, -3.48835875343438, 1.28796075006388, -0.741782397916244, 0.547300160534042, -0.492350938315503, 0.555704981283736, -1.28483063269969};
    const real_t lagrange_deriv_2[9] = {-9.73870165721010, 5.78680581663678, -3.21298543326520e-13, -2.83445891207935, 1.26941308635811, -0.855726185092640, 0.738349277190360, -0.816756381741392, 1.87444087344708};
    const real_t lagrange_deriv_3[9] = {5.54496390694879, -2.69606544031400, 3.57668094012577, -3.10862446895044e-15, -2.65931021757391, 1.37696489376050, -1.07980381128263, 1.14565373845518, -2.59074567655957};
    const real_t lagrange_deriv_4[9] = {-3.65714285714248, 1.66522164500537, -1.71783215719513, 2.85191596846290, -1.06581410364015e-14, -2.85191596846287, 1.71783215719506, -1.66522164500546, 3.65714285714316};
    const real_t lagrange_deriv_5[9] = {2.59074567655910, -1.14565373845513, 1.07980381128268, -1.37696489376052, 2.65931021757393, -3.71924713249427e-14, -3.57668094012566, 2.69606544031420, -5.54496390694988};
    const real_t lagrange_deriv_6[9] = {-1.87444087344680, 0.816756381741386, -0.738349277190418, 0.855726185092682, -1.26941308635816, 2.83445891207943, 9.27036225562006e-14, -5.78680581663756, 9.73870165721225};
    const real_t lagrange_deriv_7[9] = {1.28483063269940, -0.555704981283693, 0.492350938315507, -0.547300160534031, 0.741782397916225, -1.28796075006385, 3.48835875343430, 5.71542813077031e-13, -24.3497451715928};
    const real_t lagrange_deriv_8[9] = {-0.499999999999903, 0.215654018702479, -0.189655591978347, 0.207734512035579, -0.273437499999975, 0.444613449281046, -0.985360090074401, 4.08701370203326, 17.9999999999994};

    const real_t inv_jac[3] = {2/m_dx[0], 2/m_dx[1], 2/m_dx[2]}; //inverse jacobi
    const int numComp = in.getDataPointSize();
    const Scalar zero = static_cast<Scalar>(0);
    out.requireWrite();
    if (!in.actsExpanded()) {
#pragma omp parallel for
        for (int ei = 0; ei < m_NE[2]; ++ei) {
            for (int ej = 0; ej < m_NE[1]; ++ej) {
                for (int ek = 0; ek < m_NE[0]; ++ek) {
                    const Scalar* e = in.getSampleDataRO(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]), zero);
                    Scalar* grad = out.getSampleDataRW(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]), zero);
                    for (int comp = 0; comp < numComp; ++comp) {
                        const Scalar a = (lagrange_deriv_0[0] * e[comp] + lagrange_deriv_1[0] * e[comp] + lagrange_deriv_2[0] * e[comp] + lagrange_deriv_3[0] * e[comp] + lagrange_deriv_4[0] * e[comp] + lagrange_deriv_5[0] * e[comp] + lagrange_deriv_6[0] * e[comp] + lagrange_deriv_7[0] * e[comp] + lagrange_deriv_8[0] * e[comp]) * inv_jac[0];
                        const Scalar b = (lagrange_deriv_0[0] * e[comp] + lagrange_deriv_1[0] * e[comp] + lagrange_deriv_2[0] * e[comp] + lagrange_deriv_3[0] * e[comp] + lagrange_deriv_4[0] * e[comp] + lagrange_deriv_5[0] * e[comp] + lagrange_deriv_6[0] * e[comp] + lagrange_deriv_7[0] * e[comp] + lagrange_deriv_8[0] * e[comp]) * inv_jac[1];
                        const Scalar c = (lagrange_deriv_0[0] * e[comp] + lagrange_deriv_1[0] * e[comp] + lagrange_deriv_2[0] * e[comp] + lagrange_deriv_3[0] * e[comp] + lagrange_deriv_4[0] * e[comp] + lagrange_deriv_5[0] * e[comp] + lagrange_deriv_6[0] * e[comp] + lagrange_deriv_7[0] * e[comp] + lagrange_deriv_8[0] * e[comp]) * inv_jac[2];
                        for (int k = 0; k < 9; ++k) {
                            for (int j = 0; j < 9; ++j) {
                                for (int i = 0; i < 9; ++i) {
                                    const index_t ind = INDEX5(comp,0,i,j,k,numComp,3,10,10);
                                    grad[ind + 0] = a;
                                    grad[ind + 1] = b;
                                    grad[ind + 2] = c;
                                }
                            }
                        }
                    }
                }
            }
        }
    } else {
#pragma omp parallel for
        for (int ei = 0; ei < m_NE[2]; ++ei) {
            for (int ej = 0; ej < m_NE[1]; ++ej) {
                for (int ek = 0; ek < m_NE[0]; ++ek) {
                    const Scalar* e = in.getSampleDataRO(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]), zero);
                    Scalar* grad = out.getSampleDataRW(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]), zero);
                    for (int k = 0; k < 9; ++k) {
                        for (int j = 0; j < 9; ++j) {
                            for (int i = 0; i < 9; ++i) {
                                for (int comp = 0; comp < numComp; ++comp) {
                                    grad[INDEX5(comp,0,i,j,k,numComp,3,9,9)] = (lagrange_deriv_0[i] * e[INDEX4(comp,0,j,k,numComp,9,9)] + lagrange_deriv_1[i] * e[INDEX4(comp,1,j,k,numComp,9,9)] + lagrange_deriv_2[i] * e[INDEX4(comp,2,j,k,numComp,9,9)] + lagrange_deriv_3[i] * e[INDEX4(comp,3,j,k,numComp,9,9)] + lagrange_deriv_4[i] * e[INDEX4(comp,4,j,k,numComp,9,9)] + lagrange_deriv_5[i] * e[INDEX4(comp,5,j,k,numComp,9,9)] + lagrange_deriv_6[i] * e[INDEX4(comp,6,j,k,numComp,9,9)] + lagrange_deriv_7[i] * e[INDEX4(comp,7,j,k,numComp,9,9)] + lagrange_deriv_8[i] * e[INDEX4(comp,8,j,k,numComp,9,9)]) * inv_jac[0];
                                    grad[INDEX5(comp,1,i,j,k,numComp,3,9,9)] = (lagrange_deriv_0[j] * e[INDEX4(comp,i,0,k,numComp,9,9)] + lagrange_deriv_1[j] * e[INDEX4(comp,i,1,k,numComp,9,9)] + lagrange_deriv_2[j] * e[INDEX4(comp,i,2,k,numComp,9,9)] + lagrange_deriv_3[j] * e[INDEX4(comp,i,3,k,numComp,9,9)] + lagrange_deriv_4[j] * e[INDEX4(comp,i,4,k,numComp,9,9)] + lagrange_deriv_5[j] * e[INDEX4(comp,i,5,k,numComp,9,9)] + lagrange_deriv_6[j] * e[INDEX4(comp,i,6,k,numComp,9,9)] + lagrange_deriv_7[j] * e[INDEX4(comp,i,7,k,numComp,9,9)] + lagrange_deriv_8[j] * e[INDEX4(comp,i,8,k,numComp,9,9)]) * inv_jac[1];
                                    grad[INDEX5(comp,2,i,j,k,numComp,3,9,9)] = (lagrange_deriv_0[k] * e[INDEX4(comp,i,j,0,numComp,9,9)] + lagrange_deriv_1[k] * e[INDEX4(comp,i,j,1,numComp,9,9)] + lagrange_deriv_2[k] * e[INDEX4(comp,i,j,2,numComp,9,9)] + lagrange_deriv_3[k] * e[INDEX4(comp,i,j,3,numComp,9,9)] + lagrange_deriv_4[k] * e[INDEX4(comp,i,j,4,numComp,9,9)] + lagrange_deriv_5[k] * e[INDEX4(comp,i,j,5,numComp,9,9)] + lagrange_deriv_6[k] * e[INDEX4(comp,i,j,6,numComp,9,9)] + lagrange_deriv_7[k] * e[INDEX4(comp,i,j,7,numComp,9,9)] + lagrange_deriv_8[k] * e[INDEX4(comp,i,j,8,numComp,9,9)]) * inv_jac[2];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

template<typename Scalar>
void Brick::gradient_order9(escript::Data& out, const escript::Data& in) const {

    const real_t lagrange_deriv_0[10] = {-22.4999999999988, -5.07406470297709, 1.20335199285206, -0.528369376820220, 0.312047255608382, -0.223527944742433, 0.186645789393719, -0.180786585489230, 0.212702758009187, -0.500000000000077};
    const real_t lagrange_deriv_1[10] = {30.4381450292820, -1.52677870346452e-12, -4.25929735496529, 1.52990263818163, -0.845813573406436, 0.588082143045176, -0.483462326333953, 0.464274958908154, -0.543753738235757, 1.27595483609299};
    const real_t lagrange_deriv_2[10] = {-12.1779467074315, 7.18550286970643, 3.27293747659496e-13, -3.36412586829791, 1.44485031560171, -0.916555180336469, 0.721237312721631, -0.676797087196100, 0.783239293138005, -1.82956393190377};
    const real_t lagrange_deriv_3[10] = {6.94378848513465, -3.35166386274684, 4.36867455701003, 3.17523785042795e-14, -3.02021795819936, 1.46805550939000, -1.04618936550250, 0.936603213139437, -1.05915446364554, 2.45288417544331};
    const real_t lagrange_deriv_4[10] = {-4.59935476110357, 2.07820799403643, -2.10435017941307, 3.38731810120242, 2.22044604925031e-16, -3.02518848775198, 1.64649408398706, -1.33491548387823, 1.44494844875159, -3.29464303375000};
    const real_t lagrange_deriv_5[10] = {3.29464303374949, -1.44494844875146, 1.33491548387820, -1.64649408398705, 3.02518848775197, 6.66133814775094e-15, -3.38731810120245, 2.10435017941312, -2.07820799403661, 4.59935476110425};
    const real_t lagrange_deriv_6[10] = {-2.45288417544291, 1.05915446364544, -0.936603213139410, 1.04618936550249, -1.46805550938999, 3.02021795819934, -1.55431223447522e-15, -4.36867455701010, 3.35166386274711, -6.94378848513561};
    const real_t lagrange_deriv_7[10] = {1.82956393190345, -0.783239293137921, 0.676797087196070, -0.721237312721610, 0.916555180336448, -1.44485031560168, 3.36412586829786, -2.10609307771392e-13, -7.18550286970689, 12.1779467074328};
    const real_t lagrange_deriv_8[10] = {-1.27595483609268, 0.543753738235664, -0.464274958908104, 0.483462326333909, -0.588082143045126, 0.845813573406365, -1.52990263818151, 4.25929735496501, 2.64499533386697e-12, -30.4381450292815};
    const real_t lagrange_deriv_9[10] = {0.499999999999919, -0.212702758009135, 0.180786585489197, -0.186645789393687, 0.223527944742396, -0.312047255608330, 0.528369376820133, -1.20335199285186, 5.07406470297626, 22.4999999999976};

    const real_t inv_jac[3] = {2/m_dx[0], 2/m_dx[1], 2/m_dx[2]}; //inverse jacobi
    const int numComp = in.getDataPointSize();
    const Scalar zero = static_cast<Scalar>(0);
    out.requireWrite();
    if (!in.actsExpanded()) {
#pragma omp parallel for
        for (int ei = 0; ei < m_NE[2]; ++ei) {
            for (int ej = 0; ej < m_NE[1]; ++ej) {
                for (int ek = 0; ek < m_NE[0]; ++ek) {
                    const Scalar* e = in.getSampleDataRO(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]), zero);
                    Scalar* grad = out.getSampleDataRW(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]), zero);
                    for (int comp = 0; comp < numComp; ++comp) {
                        const Scalar a = (lagrange_deriv_0[0] * e[comp] + lagrange_deriv_1[0] * e[comp] + lagrange_deriv_2[0] * e[comp] + lagrange_deriv_3[0] * e[comp] + lagrange_deriv_4[0] * e[comp] + lagrange_deriv_5[0] * e[comp] + lagrange_deriv_6[0] * e[comp] + lagrange_deriv_7[0] * e[comp] + lagrange_deriv_8[0] * e[comp] + lagrange_deriv_9[0] * e[comp]) * inv_jac[0];
                        const Scalar b = (lagrange_deriv_0[0] * e[comp] + lagrange_deriv_1[0] * e[comp] + lagrange_deriv_2[0] * e[comp] + lagrange_deriv_3[0] * e[comp] + lagrange_deriv_4[0] * e[comp] + lagrange_deriv_5[0] * e[comp] + lagrange_deriv_6[0] * e[comp] + lagrange_deriv_7[0] * e[comp] + lagrange_deriv_8[0] * e[comp] + lagrange_deriv_9[0] * e[comp]) * inv_jac[1];
                        const Scalar c = (lagrange_deriv_0[0] * e[comp] + lagrange_deriv_1[0] * e[comp] + lagrange_deriv_2[0] * e[comp] + lagrange_deriv_3[0] * e[comp] + lagrange_deriv_4[0] * e[comp] + lagrange_deriv_5[0] * e[comp] + lagrange_deriv_6[0] * e[comp] + lagrange_deriv_7[0] * e[comp] + lagrange_deriv_8[0] * e[comp] + lagrange_deriv_9[0] * e[comp]) * inv_jac[2];
                        for (int k = 0; k < 10; ++k) {
                            for (int j = 0; j < 10; ++j) {
                                for (int i = 0; i < 10; ++i) {
                                    const index_t ind = INDEX5(comp,0,i,j,k,numComp,3,10,10);
                                    grad[ind + 0] = a;
                                    grad[ind + 1] = b;
                                    grad[ind + 2] = c;
                                }
                            }
                        }
                    }
                }
            }
        }
    } else {
#pragma omp parallel for
        for (int ei = 0; ei < m_NE[2]; ++ei) {
            for (int ej = 0; ej < m_NE[1]; ++ej) {
                for (int ek = 0; ek < m_NE[0]; ++ek) {
                    const Scalar* e = in.getSampleDataRO(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]), zero);
                    Scalar* grad = out.getSampleDataRW(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]), zero);
                    for (int k = 0; k < 10; ++k) {
                        for (int j = 0; j < 10; ++j) {
                            for (int i = 0; i < 10; ++i) {
                                for (int comp = 0; comp < numComp; ++comp) {
                                    grad[INDEX5(comp,0,i,j,k,numComp,3,10,10)] = (lagrange_deriv_0[i] * e[INDEX4(comp,0,j,k,numComp,10,10)] + lagrange_deriv_1[i] * e[INDEX4(comp,1,j,k,numComp,10,10)] + lagrange_deriv_2[i] * e[INDEX4(comp,2,j,k,numComp,10,10)] + lagrange_deriv_3[i] * e[INDEX4(comp,3,j,k,numComp,10,10)] + lagrange_deriv_4[i] * e[INDEX4(comp,4,j,k,numComp,10,10)] + lagrange_deriv_5[i] * e[INDEX4(comp,5,j,k,numComp,10,10)] + lagrange_deriv_6[i] * e[INDEX4(comp,6,j,k,numComp,10,10)] + lagrange_deriv_7[i] * e[INDEX4(comp,7,j,k,numComp,10,10)] + lagrange_deriv_8[i] * e[INDEX4(comp,8,j,k,numComp,10,10)] + lagrange_deriv_9[i] * e[INDEX4(comp,9,j,k,numComp,10,10)]) * inv_jac[0];
                                    grad[INDEX5(comp,1,i,j,k,numComp,3,10,10)] = (lagrange_deriv_0[j] * e[INDEX4(comp,i,0,k,numComp,10,10)] + lagrange_deriv_1[j] * e[INDEX4(comp,i,1,k,numComp,10,10)] + lagrange_deriv_2[j] * e[INDEX4(comp,i,2,k,numComp,10,10)] + lagrange_deriv_3[j] * e[INDEX4(comp,i,3,k,numComp,10,10)] + lagrange_deriv_4[j] * e[INDEX4(comp,i,4,k,numComp,10,10)] + lagrange_deriv_5[j] * e[INDEX4(comp,i,5,k,numComp,10,10)] + lagrange_deriv_6[j] * e[INDEX4(comp,i,6,k,numComp,10,10)] + lagrange_deriv_7[j] * e[INDEX4(comp,i,7,k,numComp,10,10)] + lagrange_deriv_8[j] * e[INDEX4(comp,i,8,k,numComp,10,10)] + lagrange_deriv_9[j] * e[INDEX4(comp,i,9,k,numComp,10,10)]) * inv_jac[1];
                                    grad[INDEX5(comp,2,i,j,k,numComp,3,10,10)] = (lagrange_deriv_0[k] * e[INDEX4(comp,i,j,0,numComp,10,10)] + lagrange_deriv_1[k] * e[INDEX4(comp,i,j,1,numComp,10,10)] + lagrange_deriv_2[k] * e[INDEX4(comp,i,j,2,numComp,10,10)] + lagrange_deriv_3[k] * e[INDEX4(comp,i,j,3,numComp,10,10)] + lagrange_deriv_4[k] * e[INDEX4(comp,i,j,4,numComp,10,10)] + lagrange_deriv_5[k] * e[INDEX4(comp,i,j,5,numComp,10,10)] + lagrange_deriv_6[k] * e[INDEX4(comp,i,j,6,numComp,10,10)] + lagrange_deriv_7[k] * e[INDEX4(comp,i,j,7,numComp,10,10)] + lagrange_deriv_8[k] * e[INDEX4(comp,i,j,8,numComp,10,10)] + lagrange_deriv_9[k] * e[INDEX4(comp,i,j,9,numComp,10,10)]) * inv_jac[2];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

template<typename Scalar>
void Brick::gradient_order10(escript::Data& out, const escript::Data& in) const {

    const real_t lagrange_deriv_0[11] = {-27.4999999999896, -6.17098569730879, 1.44617248279108, -0.622725214738251, 0.357476373116252, -0.246093749999837, 0.194287668796289, -0.172970108511720, 0.174657862947473, -0.210587346312973, 0.500000000000200};
    const real_t lagrange_deriv_1[11] = {37.2028673819635, -1.45093936865237e-11, -5.11821182477732, 1.80264679987985, -0.968487138802308, 0.646939963832195, -0.502643313012917, 0.443395636437341, -0.445313527290911, 0.535331085929629, -1.26956267628665};
    const real_t lagrange_deriv_2[11] = {-14.8873962951139, 8.73967005352577, 4.30699920173083e-12, -3.96199657704633, 1.65273663422738, -1.00650540857701, 0.747734824688410, -0.643586209360019, 0.637362056418227, -0.760400982244270, 1.79798803579712};
    const real_t lagrange_deriv_3[11] = {8.49561949463899, -4.07931619371279, 5.25066175545241, -4.06341627012807e-13, -3.45061471617355, 1.60812790372552, -1.07998725049569, 0.884587314556481, -0.852916813558380, 1.00338624297419, -2.35976991309107};
    const real_t lagrange_deriv_4[11] = {-5.64038799768972, 2.53474117868649, -2.53318340860873, 3.99079578802906, 9.19264664389630e-14, -3.30517685337838, 1.69057057046883, -1.24905529156928, 1.14606853428002, -1.31552671443958, 3.06553920088515};
    const real_t lagrange_deriv_5[11] = {4.06349206349476, -1.77190705526895, 1.61441910793459, -1.94634945457108, 3.45885134807703, 3.86357612569554e-14, -3.45885134807708, 1.94634945457117, -1.61441910793481, 1.77190705526946, -4.06349206349634};
    const real_t lagrange_deriv_6[11] = {-3.06553920088389, 1.31552671443917, -1.14606853427984, 1.24905529156920, -1.69057057046877, 3.30517685337831, -4.28546087505310e-14, -3.99079578802918, 2.53318340860904, -2.53474117868717, 5.64038799769177};
    const real_t lagrange_deriv_7[11] = {2.35976991309003, -1.00338624297385, 0.852916813558223, -0.884587314556402, 1.07998725049562, -1.60812790372545, 3.45061471617347, 5.73430192218893e-13, -5.25066175545293, 4.07931619371373, -8.49561949464169};
    const real_t lagrange_deriv_8[11] = {-1.79798803579616, 0.760400982243941, -0.637362056418051, 0.643586209359903, -0.747734824688294, 1.00650540857687, -1.65273663422718, 3.96199657704598, -3.31623617455534e-12, -8.73967005352657, 14.8873962951164};
    const real_t lagrange_deriv_9[11] = {1.26956267628575, -0.535331085929307, 0.445313527290713, -0.443395636437185, 0.502643313012754, -0.646939963831994, 0.968487138802021, -1.80264679987934, 5.11821182477613, 1.64738223062955e-11, -37.2028673819613};
    const real_t lagrange_deriv_10[11] = {-0.499999999999790, 0.210587346312823, -0.174657862947375, 0.172970108511639, -0.194287668796203, 0.246093749999731, -0.357476373116102, 0.622725214737994, -1.44617248279054, 6.17098569730708, 27.4999999999864};

    const real_t inv_jac[3] = {2/m_dx[0], 2/m_dx[1], 2/m_dx[2]}; //inverse jacobi
    const int numComp = in.getDataPointSize();
    const Scalar zero = static_cast<Scalar>(0);
    out.requireWrite();
    if (!in.actsExpanded()) {
#pragma omp parallel for
        for (int ei = 0; ei < m_NE[2]; ++ei) {
            for (int ej = 0; ej < m_NE[1]; ++ej) {
                for (int ek = 0; ek < m_NE[0]; ++ek) {
                    const Scalar* e = in.getSampleDataRO(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]), zero);
                    Scalar* grad = out.getSampleDataRW(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]), zero);
                    for (int comp = 0; comp < numComp; ++comp) {
                        const Scalar a = (lagrange_deriv_0[0] * e[comp] + lagrange_deriv_1[0] * e[comp] + lagrange_deriv_2[0] * e[comp] + lagrange_deriv_3[0] * e[comp] + lagrange_deriv_4[0] * e[comp] + lagrange_deriv_5[0] * e[comp] + lagrange_deriv_6[0] * e[comp] + lagrange_deriv_7[0] * e[comp] + lagrange_deriv_8[0] * e[comp] + lagrange_deriv_9[0] * e[comp] + lagrange_deriv_10[0] * e[comp]) * inv_jac[0];
                        const Scalar b = (lagrange_deriv_0[0] * e[comp] + lagrange_deriv_1[0] * e[comp] + lagrange_deriv_2[0] * e[comp] + lagrange_deriv_3[0] * e[comp] + lagrange_deriv_4[0] * e[comp] + lagrange_deriv_5[0] * e[comp] + lagrange_deriv_6[0] * e[comp] + lagrange_deriv_7[0] * e[comp] + lagrange_deriv_8[0] * e[comp] + lagrange_deriv_9[0] * e[comp] + lagrange_deriv_10[0] * e[comp]) * inv_jac[1];
                        const Scalar c = (lagrange_deriv_0[0] * e[comp] + lagrange_deriv_1[0] * e[comp] + lagrange_deriv_2[0] * e[comp] + lagrange_deriv_3[0] * e[comp] + lagrange_deriv_4[0] * e[comp] + lagrange_deriv_5[0] * e[comp] + lagrange_deriv_6[0] * e[comp] + lagrange_deriv_7[0] * e[comp] + lagrange_deriv_8[0] * e[comp] + lagrange_deriv_9[0] * e[comp] + lagrange_deriv_10[0] * e[comp]) * inv_jac[2];
                        for (int k = 0; k < 11; ++k) {
                            for (int j = 0; j < 11; ++j) {
                                for (int i = 0; i < 11; ++i) {
                                    const index_t ind = INDEX5(0,comp,i,j,k,3,numComp,11,11);
                                    grad[ind + 0] = a;
                                    grad[ind + 1] = b;
                                    grad[ind + 2] = c;
                                }
                            }
                        }
                    }
                }
            }
        }
    } else {
#pragma omp parallel for
        for (int ei = 0; ei < m_NE[2]; ++ei) {
            for (int ej = 0; ej < m_NE[1]; ++ej) {
                for (int ek = 0; ek < m_NE[0]; ++ek) {
                    const Scalar* e = in.getSampleDataRO(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]), zero);
                    Scalar* grad = out.getSampleDataRW(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]), zero);
                    for (int k = 0; k < 11; ++k) {
                        for (int j = 0; j < 11; ++j) {
                            for (int i = 0; i < 11; ++i) {
                                for (int comp = 0; comp < numComp; ++comp) {
                                    grad[INDEX5(comp,0,i,j,k,numComp,3,11,11)] = (lagrange_deriv_0[i] * e[INDEX4(comp,0,j,k,numComp,11,11)] + lagrange_deriv_1[i] * e[INDEX4(comp,1,j,k,numComp,11,11)] + lagrange_deriv_2[i] * e[INDEX4(comp,2,j,k,numComp,11,11)] + lagrange_deriv_3[i] * e[INDEX4(comp,3,j,k,numComp,11,11)] + lagrange_deriv_4[i] * e[INDEX4(comp,4,j,k,numComp,11,11)] + lagrange_deriv_5[i] * e[INDEX4(comp,5,j,k,numComp,11,11)] + lagrange_deriv_6[i] * e[INDEX4(comp,6,j,k,numComp,11,11)] + lagrange_deriv_7[i] * e[INDEX4(comp,7,j,k,numComp,11,11)] + lagrange_deriv_8[i] * e[INDEX4(comp,8,j,k,numComp,11,11)] + lagrange_deriv_9[i] * e[INDEX4(comp,9,j,k,numComp,11,11)] + lagrange_deriv_10[i] * e[INDEX4(comp,10,j,k,numComp,11,11)]) * inv_jac[0];
                                    grad[INDEX5(comp,1,i,j,k,numComp,3,11,11)] = (lagrange_deriv_0[j] * e[INDEX4(comp,i,0,k,numComp,11,11)] + lagrange_deriv_1[j] * e[INDEX4(comp,i,1,k,numComp,11,11)] + lagrange_deriv_2[j] * e[INDEX4(comp,i,2,k,numComp,11,11)] + lagrange_deriv_3[j] * e[INDEX4(comp,i,3,k,numComp,11,11)] + lagrange_deriv_4[j] * e[INDEX4(comp,i,4,k,numComp,11,11)] + lagrange_deriv_5[j] * e[INDEX4(comp,i,5,k,numComp,11,11)] + lagrange_deriv_6[j] * e[INDEX4(comp,i,6,k,numComp,11,11)] + lagrange_deriv_7[j] * e[INDEX4(comp,i,7,k,numComp,11,11)] + lagrange_deriv_8[j] * e[INDEX4(comp,i,8,k,numComp,11,11)] + lagrange_deriv_9[j] * e[INDEX4(comp,i,9,k,numComp,11,11)] + lagrange_deriv_10[j] * e[INDEX4(comp,i,10,k,numComp,11,11)]) * inv_jac[1];
                                    grad[INDEX5(comp,2,i,j,k,numComp,3,11,11)] = (lagrange_deriv_0[k] * e[INDEX4(comp,i,j,0,numComp,11,11)] + lagrange_deriv_1[k] * e[INDEX4(comp,i,j,1,numComp,11,11)] + lagrange_deriv_2[k] * e[INDEX4(comp,i,j,2,numComp,11,11)] + lagrange_deriv_3[k] * e[INDEX4(comp,i,j,3,numComp,11,11)] + lagrange_deriv_4[k] * e[INDEX4(comp,i,j,4,numComp,11,11)] + lagrange_deriv_5[k] * e[INDEX4(comp,i,j,5,numComp,11,11)] + lagrange_deriv_6[k] * e[INDEX4(comp,i,j,6,numComp,11,11)] + lagrange_deriv_7[k] * e[INDEX4(comp,i,j,7,numComp,11,11)] + lagrange_deriv_8[k] * e[INDEX4(comp,i,j,8,numComp,11,11)] + lagrange_deriv_9[k] * e[INDEX4(comp,i,j,9,numComp,11,11)] + lagrange_deriv_10[k] * e[INDEX4(comp,i,j,10,numComp,11,11)]) * inv_jac[2];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

// instantiate
template
void Brick::gradient_order2<real_t>(escript::Data& out,
                                    const escript::Data& in) const;
template
void Brick::gradient_order2<cplx_t>(escript::Data& out,
                                    const escript::Data& in) const;
template
void Brick::gradient_order3<real_t>(escript::Data& out,
                                    const escript::Data& in) const;
template
void Brick::gradient_order3<cplx_t>(escript::Data& out,
                                    const escript::Data& in) const;
template
void Brick::gradient_order4<real_t>(escript::Data& out,
                                    const escript::Data& in) const;
template
void Brick::gradient_order4<cplx_t>(escript::Data& out,
                                    const escript::Data& in) const;
template
void Brick::gradient_order5<real_t>(escript::Data& out,
                                    const escript::Data& in) const;
template
void Brick::gradient_order5<cplx_t>(escript::Data& out,
                                    const escript::Data& in) const;
template
void Brick::gradient_order6<real_t>(escript::Data& out,
                                    const escript::Data& in) const;
template
void Brick::gradient_order6<cplx_t>(escript::Data& out,
                                    const escript::Data& in) const;
template
void Brick::gradient_order7<real_t>(escript::Data& out,
                                    const escript::Data& in) const;
template
void Brick::gradient_order7<cplx_t>(escript::Data& out,
                                    const escript::Data& in) const;
template
void Brick::gradient_order8<real_t>(escript::Data& out,
                                    const escript::Data& in) const;
template
void Brick::gradient_order8<cplx_t>(escript::Data& out,
                                    const escript::Data& in) const;
template
void Brick::gradient_order9<real_t>(escript::Data& out,
                                    const escript::Data& in) const;
template
void Brick::gradient_order9<cplx_t>(escript::Data& out,
                                    const escript::Data& in) const;
template
void Brick::gradient_order10<real_t>(escript::Data& out,
                                     const escript::Data& in) const;
template
void Brick::gradient_order10<cplx_t>(escript::Data& out,
                                     const escript::Data& in) const;

} // namespace speckley

