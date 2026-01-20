
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
void Brick::integral_order2(std::vector<Scalar>& integrals, const escript::Data& arg) const
{
    const double weights[] = {0.333333333333, 1.33333333333, 0.333333333333};
    const int numComp = arg.getDataPointSize();
    const double volume_product = 0.125*m_dx[0]*m_dx[1]*m_dx[2];
    const Scalar zero = static_cast<Scalar>(0);
    for (int ei = 0; ei < m_NE[2]; ++ei) {
        for (int ej = 0; ej < m_NE[1]; ++ej) {
            for (int ek = 0; ek < m_NE[0]; ++ek) {
                const Scalar* e = arg.getSampleDataRO(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]), zero);
                Scalar result = zero;
                for (int comp = 0; comp < numComp; ++comp) {
                    for (int i = 0; i < 3; ++i) {
                        for (int j = 0; j < 3; ++j) {
                            for (int k = 0; k < 3; ++k) {
                                result += weights[i] * weights[j] * weights[k] * e[INDEX4(comp,i,j,k,numComp,3,3)];
                            }
                        }
                    }
                    integrals[comp] += result;
                }
            }
        }
    }
    for (int comp = 0; comp < numComp; ++comp) {
        integrals[comp] *= volume_product;
    }
}

template<typename Scalar>
void Brick::integral_order3(std::vector<Scalar>& integrals, const escript::Data& arg) const {
    const double weights[] = {0.166666666667, 0.833333333333, 0.833333333333, 0.166666666667};
    const int numComp = arg.getDataPointSize();
    const double volume_product = 0.125*m_dx[0]*m_dx[1]*m_dx[2];
    const Scalar zero = static_cast<Scalar>(0);
    for (int ei = 0; ei < m_NE[2]; ++ei) {
        for (int ej = 0; ej < m_NE[1]; ++ej) {
            for (int ek = 0; ek < m_NE[0]; ++ek) {
                const Scalar* e = arg.getSampleDataRO(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]), zero);
                Scalar result = zero;
                for (int comp = 0; comp < numComp; ++comp) {
                    for (int i = 0; i < 4; ++i) {
                        for (int j = 0; j < 4; ++j) {
                            for (int k = 0; k < 4; ++k) {
                                result += weights[i] * weights[j] * weights[k] * e[INDEX4(comp,i,j,k,numComp,4,4)];
                            }
                        }
                    }
                    integrals[comp] += result;
                }
            }
        }
    }
    for (int comp = 0; comp < numComp; ++comp) {
        integrals[comp] *= volume_product;
    }
}

template<typename Scalar>
void Brick::integral_order4(std::vector<Scalar>& integrals, const escript::Data& arg) const {
    const double weights[] = {0.1, 0.544444444444, 0.711111111111, 0.544444444444, 0.1};
    const int numComp = arg.getDataPointSize();
    const double volume_product = 0.125*m_dx[0]*m_dx[1]*m_dx[2];
    const Scalar zero = static_cast<Scalar>(0);
    for (int ei = 0; ei < m_NE[2]; ++ei) {
        for (int ej = 0; ej < m_NE[1]; ++ej) {
            for (int ek = 0; ek < m_NE[0]; ++ek) {
                const Scalar* e = arg.getSampleDataRO(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]), zero);
                Scalar result = zero;
                for (int comp = 0; comp < numComp; ++comp) {
                    for (int i = 0; i < 5; ++i) {
                        for (int j = 0; j < 5; ++j) {
                            for (int k = 0; k < 5; ++k) {
                                result += weights[i] * weights[j] * weights[k] * e[INDEX4(comp,i,j,k,numComp,5,5)];
                            }
                        }
                    }
                    integrals[comp] += result;
                }
            }
        }
    }
    for (int comp = 0; comp < numComp; ++comp) {
        integrals[comp] *= volume_product;
    }
}

template<typename Scalar>
void Brick::integral_order5(std::vector<Scalar>& integrals, const escript::Data& arg) const {
    const double weights[] = {0.0666666666667, 0.378474956298, 0.554858377035, 0.554858377035, 0.378474956298, 0.0666666666667};
    const int numComp = arg.getDataPointSize();
    const double volume_product = 0.125*m_dx[0]*m_dx[1]*m_dx[2];
    const Scalar zero = static_cast<Scalar>(0);
    for (int ei = 0; ei < m_NE[2]; ++ei) {
        for (int ej = 0; ej < m_NE[1]; ++ej) {
            for (int ek = 0; ek < m_NE[0]; ++ek) {
                const Scalar* e = arg.getSampleDataRO(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]), zero);
                Scalar result = zero;
                for (int comp = 0; comp < numComp; ++comp) {
                    for (int i = 0; i < 6; ++i) {
                        for (int j = 0; j < 6; ++j) {
                            for (int k = 0; k < 6; ++k) {
                                result += weights[i] * weights[j] * weights[k] * e[INDEX4(comp,i,j,k,numComp,6,6)];
                            }
                        }
                    }
                    integrals[comp] += result;
                }
            }
        }
    }
    for (int comp = 0; comp < numComp; ++comp) {
        integrals[comp] *= volume_product;
    }
}

template<typename Scalar>
void Brick::integral_order6(std::vector<Scalar>& integrals, const escript::Data& arg) const {
    const double weights[] = {0.047619047619, 0.276826047362, 0.43174538121, 0.487619047619, 0.43174538121, 0.276826047362, 0.047619047619};
    const int numComp = arg.getDataPointSize();
    const double volume_product = 0.125*m_dx[0]*m_dx[1]*m_dx[2];
    const Scalar zero = static_cast<Scalar>(0);
    for (int ei = 0; ei < m_NE[2]; ++ei) {
        for (int ej = 0; ej < m_NE[1]; ++ej) {
            for (int ek = 0; ek < m_NE[0]; ++ek) {
                const Scalar* e = arg.getSampleDataRO(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]), zero);
                Scalar result = zero;
                for (int comp = 0; comp < numComp; ++comp) {
                    for (int i = 0; i < 7; ++i) {
                        for (int j = 0; j < 7; ++j) {
                            for (int k = 0; k < 7; ++k) {
                                result += weights[i] * weights[j] * weights[k] * e[INDEX4(comp,i,j,k,numComp,7,7)];
                            }
                        }
                    }
                    integrals[comp] += result;
                }
            }
        }
    }
    for (int comp = 0; comp < numComp; ++comp) {
        integrals[comp] *= volume_product;
    }
}

template<typename Scalar>
void Brick::integral_order7(std::vector<Scalar>& integrals, const escript::Data& arg) const {
    const double weights[] = {0.0357142857143, 0.210704227144, 0.341122692484, 0.412458794659, 0.412458794659, 0.341122692484, 0.210704227144, 0.0357142857143};
    const int numComp = arg.getDataPointSize();
    const double volume_product = 0.125*m_dx[0]*m_dx[1]*m_dx[2];
    const Scalar zero = static_cast<Scalar>(0);
    for (int ei = 0; ei < m_NE[2]; ++ei) {
        for (int ej = 0; ej < m_NE[1]; ++ej) {
            for (int ek = 0; ek < m_NE[0]; ++ek) {
                const Scalar* e = arg.getSampleDataRO(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]), zero);
                Scalar result = zero;
                for (int comp = 0; comp < numComp; ++comp) {
                    for (int i = 0; i < 8; ++i) {
                        for (int j = 0; j < 8; ++j) {
                            for (int k = 0; k < 8; ++k) {
                                result += weights[i] * weights[j] * weights[k] * e[INDEX4(comp,i,j,k,numComp,8,8)];
                            }
                        }
                    }
                    integrals[comp] += result;
                }
            }
        }
    }
    for (int comp = 0; comp < numComp; ++comp) {
        integrals[comp] *= volume_product;
    }
}

template<typename Scalar>
void Brick::integral_order8(std::vector<Scalar>& integrals, const escript::Data& arg) const {
    const double weights[] = {0.0277777777778, 0.165495361561, 0.2745387125, 0.346428510973, 0.371519274376, 0.346428510973, 0.2745387125, 0.165495361561, 0.0277777777778};
    const int numComp = arg.getDataPointSize();
    const double volume_product = 0.125*m_dx[0]*m_dx[1]*m_dx[2];
    const Scalar zero = static_cast<Scalar>(0);
    for (int ei = 0; ei < m_NE[2]; ++ei) {
        for (int ej = 0; ej < m_NE[1]; ++ej) {
            for (int ek = 0; ek < m_NE[0]; ++ek) {
                const Scalar* e = arg.getSampleDataRO(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]), zero);
                Scalar result = zero;
                for (int comp = 0; comp < numComp; ++comp) {
                    for (int i = 0; i < 9; ++i) {
                        for (int j = 0; j < 9; ++j) {
                            for (int k = 0; k < 9; ++k) {
                                result += weights[i] * weights[j] * weights[k] * e[INDEX4(comp,i,j,k,numComp,9,9)];
                            }
                        }
                    }
                    integrals[comp] += result;
                }
            }
        }
    }
    for (int comp = 0; comp < numComp; ++comp) {
        integrals[comp] *= volume_product;
    }
}

template<typename Scalar>
void Brick::integral_order9(std::vector<Scalar>& integrals, const escript::Data& arg) const {
    const double weights[] = {0.0222222222222, 0.133305990851, 0.224889342063, 0.29204268368, 0.327539761184, 0.327539761184, 0.29204268368, 0.224889342063, 0.133305990851, 0.0222222222222};
    const int numComp = arg.getDataPointSize();
    const double volume_product = 0.125*m_dx[0]*m_dx[1]*m_dx[2];
    const Scalar zero = static_cast<Scalar>(0);
    for (int ei = 0; ei < m_NE[2]; ++ei) {
        for (int ej = 0; ej < m_NE[1]; ++ej) {
            for (int ek = 0; ek < m_NE[0]; ++ek) {
                const Scalar* e = arg.getSampleDataRO(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]), zero);
                Scalar result = zero;
                for (int comp = 0; comp < numComp; ++comp) {
                    for (int i = 0; i < 10; ++i) {
                        for (int j = 0; j < 10; ++j) {
                            for (int k = 0; k < 10; ++k) {
                                result += weights[i] * weights[j] * weights[k] * e[INDEX4(comp,i,j,k,numComp,10,10)];
                            }
                        }
                    }
                    integrals[comp] += result;
                }
            }
        }
    }
    for (int comp = 0; comp < numComp; ++comp) {
        integrals[comp] *= volume_product;
    }
}

template<typename Scalar>
void Brick::integral_order10(std::vector<Scalar>& integrals, const escript::Data& arg) const {
    const double weights[] = {0.0181818181818, 0.109612273267, 0.18716988178, 0.248048104264, 0.286879124779, 0.300217595456, 0.286879124779, 0.248048104264, 0.18716988178, 0.109612273267, 0.0181818181818};
    const int numComp = arg.getDataPointSize();
    const double volume_product = 0.125*m_dx[0]*m_dx[1]*m_dx[2];
    const Scalar zero = static_cast<Scalar>(0);
    for (int ei = 0; ei < m_NE[2]; ++ei) {
        for (int ej = 0; ej < m_NE[1]; ++ej) {
            for (int ek = 0; ek < m_NE[0]; ++ek) {
                const Scalar* e = arg.getSampleDataRO(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]), zero);
                Scalar result = zero;
                for (int comp = 0; comp < numComp; ++comp) {
                    for (int i = 0; i < 11; ++i) {
                        for (int j = 0; j < 11; ++j) {
                            for (int k = 0; k < 11; ++k) {
                                result += weights[i] * weights[j] * weights[k] * e[INDEX4(comp,i,j,k,numComp,11,11)];
                            }
                        }
                    }
                    integrals[comp] += result;
                }
            }
        }
    }
    for (int comp = 0; comp < numComp; ++comp) {
        integrals[comp] *= volume_product;
    }
}

// instantiate
template
void Brick::integral_order2<real_t>(std::vector<real_t>& integrals,
                                    const escript::Data& arg) const;
template
void Brick::integral_order2<cplx_t>(std::vector<cplx_t>& integrals,
                                    const escript::Data& arg) const;
template
void Brick::integral_order3<real_t>(std::vector<real_t>& integrals,
                                    const escript::Data& arg) const;
template
void Brick::integral_order3<cplx_t>(std::vector<cplx_t>& integrals,
                                    const escript::Data& arg) const;
template
void Brick::integral_order4<real_t>(std::vector<real_t>& integrals,
                                    const escript::Data& arg) const;
template
void Brick::integral_order4<cplx_t>(std::vector<cplx_t>& integrals,
                                    const escript::Data& arg) const;
template
void Brick::integral_order5<real_t>(std::vector<real_t>& integrals,
                                    const escript::Data& arg) const;
template
void Brick::integral_order5<cplx_t>(std::vector<cplx_t>& integrals,
                                    const escript::Data& arg) const;
template
void Brick::integral_order6<real_t>(std::vector<real_t>& integrals,
                                    const escript::Data& arg) const;
template
void Brick::integral_order6<cplx_t>(std::vector<cplx_t>& integrals,
                                    const escript::Data& arg) const;
template
void Brick::integral_order7<real_t>(std::vector<real_t>& integrals,
                                    const escript::Data& arg) const;
template
void Brick::integral_order7<cplx_t>(std::vector<cplx_t>& integrals,
                                    const escript::Data& arg) const;
template
void Brick::integral_order8<real_t>(std::vector<real_t>& integrals,
                                    const escript::Data& arg) const;
template
void Brick::integral_order8<cplx_t>(std::vector<cplx_t>& integrals,
                                    const escript::Data& arg) const;
template
void Brick::integral_order9<real_t>(std::vector<real_t>& integrals,
                                    const escript::Data& arg) const;
template
void Brick::integral_order9<cplx_t>(std::vector<cplx_t>& integrals,
                                    const escript::Data& arg) const;
template
void Brick::integral_order10<real_t>(std::vector<real_t>& integrals,
                                     const escript::Data& arg) const;
template
void Brick::integral_order10<cplx_t>(std::vector<cplx_t>& integrals,
                                     const escript::Data& arg) const;

} // namespace speckley

