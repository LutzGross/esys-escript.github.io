
/*****************************************************************************
*
* Copyright (c) 2003-2017 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

#include <speckley/Brick.h>

#include <escript/index.h>

namespace speckley {

void Brick::reduction_order2(const escript::Data& in, escript::Data& out) const
{
    const double weights[] = {0.333333333333, 1.33333333333, 0.333333333333};
    const int numComp = in.getDataPointSize();
    for (int ei = 0; ei < m_NE[2]; ++ei) {
        for (int ej = 0; ej < m_NE[1]; ++ej) {
            for (int ek = 0; ek < m_NE[0]; ++ek) {
                const double *e_in = in.getSampleDataRO(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]));
                double *e_out = out.getSampleDataRW(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]));
                for (int comp = 0; comp < numComp; ++comp) {
                    double result = 0;
                    for (int i = 0; i < 3; ++i) {
                        for (int j = 0; j < 3; ++j) {
                            for (int k = 0; k < 3; ++k) {
                                result += weights[i] * weights[j] * weights[k] * e_in[INDEX4(comp,k,j,i,numComp,3,3)];
                            }
                        }
                    }
                    e_out[comp] += result / 8.;
                }
            }
        }
    }
}

void Brick::reduction_order3(const escript::Data& in, escript::Data& out) const {
    const double weights[] = {0.166666666667, 0.833333333333, 0.833333333333, 0.166666666667};
    const int numComp = in.getDataPointSize();
    for (int ei = 0; ei < m_NE[2]; ++ei) {
        for (int ej = 0; ej < m_NE[1]; ++ej) {
            for (int ek = 0; ek < m_NE[0]; ++ek) {
                const double *e_in = in.getSampleDataRO(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]));
                double *e_out = out.getSampleDataRW(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]));
                for (int comp = 0; comp < numComp; ++comp) {
                    double result = 0;
                    for (int i = 0; i < 4; ++i) {
                        for (int j = 0; j < 4; ++j) {
                            for (int k = 0; k < 4; ++k) {
                                result += weights[i] * weights[j] * weights[k] * e_in[INDEX4(comp,k,j,i,numComp,4,4)];
                            }
                        }
                    }
                    e_out[comp] += result / 8.;
                }
            }
        }
    }
}

void Brick::reduction_order4(const escript::Data& in, escript::Data& out) const {
    const double weights[] = {0.1, 0.544444444444, 0.711111111111, 0.544444444444, 0.1};
    const int numComp = in.getDataPointSize();
    for (int ei = 0; ei < m_NE[2]; ++ei) {
        for (int ej = 0; ej < m_NE[1]; ++ej) {
            for (int ek = 0; ek < m_NE[0]; ++ek) {
                const double *e_in = in.getSampleDataRO(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]));
                double *e_out = out.getSampleDataRW(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]));
                for (int comp = 0; comp < numComp; ++comp) {
                    double result = 0;
                    for (int i = 0; i < 5; ++i) {
                        for (int j = 0; j < 5; ++j) {
                            for (int k = 0; k < 5; ++k) {
                                result += weights[i] * weights[j] * weights[k] * e_in[INDEX4(comp,k,j,i,numComp,5,5)];
                            }
                        }
                    }
                    e_out[comp] += result / 8.;
                }
            }
        }
    }
}

void Brick::reduction_order5(const escript::Data& in, escript::Data& out) const {
    const double weights[] = {0.0666666666667, 0.378474956298, 0.554858377035, 0.554858377035, 0.378474956298, 0.0666666666667};
    const int numComp = in.getDataPointSize();
    for (int ei = 0; ei < m_NE[2]; ++ei) {
        for (int ej = 0; ej < m_NE[1]; ++ej) {
            for (int ek = 0; ek < m_NE[0]; ++ek) {
                const double *e_in = in.getSampleDataRO(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]));
                double *e_out = out.getSampleDataRW(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]));
                for (int comp = 0; comp < numComp; ++comp) {
                    double result = 0;
                    for (int i = 0; i < 6; ++i) {
                        for (int j = 0; j < 6; ++j) {
                            for (int k = 0; k < 6; ++k) {
                                result += weights[i] * weights[j] * weights[k] * e_in[INDEX4(comp,k,j,i,numComp,6,6)];
                            }
                        }
                    }
                    e_out[comp] += result / 8.;
                }
            }
        }
    }
}

void Brick::reduction_order6(const escript::Data& in, escript::Data& out) const {
    const double weights[] = {0.047619047619, 0.276826047362, 0.43174538121, 0.487619047619, 0.43174538121, 0.276826047362, 0.047619047619};
    const int numComp = in.getDataPointSize();
    for (int ei = 0; ei < m_NE[2]; ++ei) {
        for (int ej = 0; ej < m_NE[1]; ++ej) {
            for (int ek = 0; ek < m_NE[0]; ++ek) {
                const double *e_in = in.getSampleDataRO(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]));
                double *e_out = out.getSampleDataRW(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]));
                for (int comp = 0; comp < numComp; ++comp) {
                    double result = 0;
                    for (int i = 0; i < 7; ++i) {
                        for (int j = 0; j < 7; ++j) {
                            for (int k = 0; k < 7; ++k) {
                                result += weights[i] * weights[j] * weights[k] * e_in[INDEX4(comp,k,j,i,numComp,7,7)];
                            }
                        }
                    }
                    e_out[comp] += result / 8.;
                }
            }
        }
    }
}

void Brick::reduction_order7(const escript::Data& in, escript::Data& out) const {
    const double weights[] = {0.0357142857143, 0.210704227144, 0.341122692484, 0.412458794659, 0.412458794659, 0.341122692484, 0.210704227144, 0.0357142857143};
    const int numComp = in.getDataPointSize();
    for (int ei = 0; ei < m_NE[2]; ++ei) {
        for (int ej = 0; ej < m_NE[1]; ++ej) {
            for (int ek = 0; ek < m_NE[0]; ++ek) {
                const double *e_in = in.getSampleDataRO(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]));
                double *e_out = out.getSampleDataRW(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]));
                for (int comp = 0; comp < numComp; ++comp) {
                    double result = 0;
                    for (int i = 0; i < 8; ++i) {
                        for (int j = 0; j < 8; ++j) {
                            for (int k = 0; k < 8; ++k) {
                                result += weights[i] * weights[j] * weights[k] * e_in[INDEX4(comp,k,j,i,numComp,8,8)];
                            }
                        }
                    }
                    e_out[comp] += result / 8.;
                }
            }
        }
    }
}

void Brick::reduction_order8(const escript::Data& in, escript::Data& out) const {
    const double weights[] = {0.0277777777778, 0.165495361561, 0.2745387125, 0.346428510973, 0.371519274376, 0.346428510973, 0.2745387125, 0.165495361561, 0.0277777777778};
    const int numComp = in.getDataPointSize();
    for (int ei = 0; ei < m_NE[2]; ++ei) {
        for (int ej = 0; ej < m_NE[1]; ++ej) {
            for (int ek = 0; ek < m_NE[0]; ++ek) {
                const double *e_in = in.getSampleDataRO(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]));
                double *e_out = out.getSampleDataRW(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]));
                for (int comp = 0; comp < numComp; ++comp) {
                    double result = 0;
                    for (int i = 0; i < 9; ++i) {
                        for (int j = 0; j < 9; ++j) {
                            for (int k = 0; k < 9; ++k) {
                                result += weights[i] * weights[j] * weights[k] * e_in[INDEX4(comp,k,j,i,numComp,9,9)];
                            }
                        }
                    }
                    e_out[comp] += result / 8.;
                }
            }
        }
    }
}

void Brick::reduction_order9(const escript::Data& in, escript::Data& out) const {
    const double weights[] = {0.0222222222222, 0.133305990851, 0.224889342063, 0.29204268368, 0.327539761184, 0.327539761184, 0.29204268368, 0.224889342063, 0.133305990851, 0.0222222222222};
    const int numComp = in.getDataPointSize();
    for (int ei = 0; ei < m_NE[2]; ++ei) {
        for (int ej = 0; ej < m_NE[1]; ++ej) {
            for (int ek = 0; ek < m_NE[0]; ++ek) {
                const double *e_in = in.getSampleDataRO(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]));
                double *e_out = out.getSampleDataRW(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]));
                for (int comp = 0; comp < numComp; ++comp) {
                    double result = 0;
                    for (int i = 0; i < 10; ++i) {
                        for (int j = 0; j < 10; ++j) {
                            for (int k = 0; k < 10; ++k) {
                                result += weights[i] * weights[j] * weights[k] * e_in[INDEX4(comp,k,j,i,numComp,10,10)];
                            }
                        }
                    }
                    e_out[comp] += result / 8.;
                }
            }
        }
    }
}

void Brick::reduction_order10(const escript::Data& in, escript::Data& out) const {
    const double weights[] = {0.0181818181818, 0.109612273267, 0.18716988178, 0.248048104264, 0.286879124779, 0.300217595456, 0.286879124779, 0.248048104264, 0.18716988178, 0.109612273267, 0.0181818181818};
    const int numComp = in.getDataPointSize();
    for (int ei = 0; ei < m_NE[2]; ++ei) {
        for (int ej = 0; ej < m_NE[1]; ++ej) {
            for (int ek = 0; ek < m_NE[0]; ++ek) {
                const double *e_in = in.getSampleDataRO(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]));
                double *e_out = out.getSampleDataRW(INDEX3(ek,ej,ei,m_NE[0],m_NE[1]));
                for (int comp = 0; comp < numComp; ++comp) {
                    double result = 0;
                    for (int i = 0; i < 11; ++i) {
                        for (int j = 0; j < 11; ++j) {
                            for (int k = 0; k < 11; ++k) {
                                result += weights[i] * weights[j] * weights[k] * e_in[INDEX4(comp,k,j,i,numComp,11,11)];
                            }
                        }
                    }
                    e_out[comp] += result / 8.;
                }
            }
        }
    }
}

}
