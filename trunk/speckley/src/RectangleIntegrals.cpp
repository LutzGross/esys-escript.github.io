#define ESNEEDPYTHON
#include "esysUtils/first.h"


#include <esysUtils/index.h>
#include <speckley/Rectangle.h>

namespace speckley {
void Rectangle::integral_order2(std::vector<double>& integrals, const escript::Data& arg) const {
    const double weights[] = {0.333333333333, 1.33333333333, 0.333333333333};
    const int numComp = arg.getDataPointSize();
    const double volume_product = 0.25*m_dx[0]*m_dx[1];
    for (int ei = 0; ei < m_NE[1]; ++ei) {
        for (int ej = 0; ej < m_NE[0]; ++ej) {
            const double *e = arg.getSampleDataRO(INDEX2(ej,ei,m_NE[0]));
            double result = 0;
            for (int comp = 0; comp < numComp; ++comp) {
                for (int i = 0; i < 3; ++i) {
                    for (int j = 0; j < 3; ++j) {
                        result += weights[i] * weights[j] * e[INDEX3(comp,i,j,numComp,3)];
                    }
                }
                integrals[comp] += result;
            }
        }
    }
    for (int comp = 0; comp < numComp; ++comp) {
        integrals[comp] *= volume_product;
    }
}

void Rectangle::integral_order3(std::vector<double>& integrals, const escript::Data& arg) const {
    const double weights[] = {0.166666666667, 0.833333333333, 0.833333333333, 0.166666666667};
    const int numComp = arg.getDataPointSize();
    const double volume_product = 0.25*m_dx[0]*m_dx[1];
    for (int ei = 0; ei < m_NE[1]; ++ei) {
        for (int ej = 0; ej < m_NE[0]; ++ej) {
            const double *e = arg.getSampleDataRO(INDEX2(ej,ei,m_NE[0]));
            double result = 0;
            for (int comp = 0; comp < numComp; ++comp) {
                for (int i = 0; i < 4; ++i) {
                    for (int j = 0; j < 4; ++j) {
                        result += weights[i] * weights[j] * e[INDEX3(comp,i,j,numComp,4)];
                    }
                }
                integrals[comp] += result;
            }
        }
    }
    for (int comp = 0; comp < numComp; ++comp) {
        integrals[comp] *= volume_product;
    }
}

void Rectangle::integral_order4(std::vector<double>& integrals, const escript::Data& arg) const {
    const double weights[] = {0.1, 0.544444444444, 0.711111111111, 0.544444444444, 0.1};
    const int numComp = arg.getDataPointSize();
    const double volume_product = 0.25*m_dx[0]*m_dx[1];
    for (int ei = 0; ei < m_NE[1]; ++ei) {
        for (int ej = 0; ej < m_NE[0]; ++ej) {
            const double *e = arg.getSampleDataRO(INDEX2(ej,ei,m_NE[0]));
            double result = 0;
            for (int comp = 0; comp < numComp; ++comp) {
                for (int i = 0; i < 5; ++i) {
                    for (int j = 0; j < 5; ++j) {
                        result += weights[i] * weights[j] * e[INDEX3(comp,i,j,numComp,5)];
                    }
                }
                integrals[comp] += result;
            }
        }
    }
    for (int comp = 0; comp < numComp; ++comp) {
        integrals[comp] *= volume_product;
    }
}

void Rectangle::integral_order5(std::vector<double>& integrals, const escript::Data& arg) const {
    const double weights[] = {0.0666666666667, 0.378474956298, 0.554858377035, 0.554858377035, 0.378474956298, 0.0666666666667};
    const int numComp = arg.getDataPointSize();
    const double volume_product = 0.25*m_dx[0]*m_dx[1];
    for (int ei = 0; ei < m_NE[1]; ++ei) {
        for (int ej = 0; ej < m_NE[0]; ++ej) {
            const double *e = arg.getSampleDataRO(INDEX2(ej,ei,m_NE[0]));
            double result = 0;
            for (int comp = 0; comp < numComp; ++comp) {
                for (int i = 0; i < 6; ++i) {
                    for (int j = 0; j < 6; ++j) {
                        result += weights[i] * weights[j] * e[INDEX3(comp,i,j,numComp,6)];
                    }
                }
                integrals[comp] += result;
            }
        }
    }
    for (int comp = 0; comp < numComp; ++comp) {
        integrals[comp] *= volume_product;
    }
}

void Rectangle::integral_order6(std::vector<double>& integrals, const escript::Data& arg) const {
    const double weights[] = {0.047619047619, 0.276826047362, 0.43174538121, 0.487619047619, 0.43174538121, 0.276826047362, 0.047619047619};
    const int numComp = arg.getDataPointSize();
    const double volume_product = 0.25*m_dx[0]*m_dx[1];
    for (int ei = 0; ei < m_NE[1]; ++ei) {
        for (int ej = 0; ej < m_NE[0]; ++ej) {
            const double *e = arg.getSampleDataRO(INDEX2(ej,ei,m_NE[0]));
            double result = 0;
            for (int comp = 0; comp < numComp; ++comp) {
                for (int i = 0; i < 7; ++i) {
                    for (int j = 0; j < 7; ++j) {
                        result += weights[i] * weights[j] * e[INDEX3(comp,i,j,numComp,7)];
                    }
                }
                integrals[comp] += result;
            }
        }
    }
    for (int comp = 0; comp < numComp; ++comp) {
        integrals[comp] *= volume_product;
    }
}

void Rectangle::integral_order7(std::vector<double>& integrals, const escript::Data& arg) const {
    const double weights[] = {0.0357142857143, 0.210704227144, 0.341122692484, 0.412458794659, 0.412458794659, 0.341122692484, 0.210704227144, 0.0357142857143};
    const int numComp = arg.getDataPointSize();
    const double volume_product = 0.25*m_dx[0]*m_dx[1];
    for (int ei = 0; ei < m_NE[1]; ++ei) {
        for (int ej = 0; ej < m_NE[0]; ++ej) {
            const double *e = arg.getSampleDataRO(INDEX2(ej,ei,m_NE[0]));
            double result = 0;
            for (int comp = 0; comp < numComp; ++comp) {
                for (int i = 0; i < 8; ++i) {
                    for (int j = 0; j < 8; ++j) {
                        result += weights[i] * weights[j] * e[INDEX3(comp,i,j,numComp,8)];
                    }
                }
                integrals[comp] += result;
            }
        }
    }
    for (int comp = 0; comp < numComp; ++comp) {
        integrals[comp] *= volume_product;
    }
}

void Rectangle::integral_order8(std::vector<double>& integrals, const escript::Data& arg) const {
    const double weights[] = {0.0277777777778, 0.165495361561, 0.2745387125, 0.346428510973, 0.371519274376, 0.346428510973, 0.2745387125, 0.165495361561, 0.0277777777778};
    const int numComp = arg.getDataPointSize();
    const double volume_product = 0.25*m_dx[0]*m_dx[1];
    for (int ei = 0; ei < m_NE[1]; ++ei) {
        for (int ej = 0; ej < m_NE[0]; ++ej) {
            const double *e = arg.getSampleDataRO(INDEX2(ej,ei,m_NE[0]));
            double result = 0;
            for (int comp = 0; comp < numComp; ++comp) {
                for (int i = 0; i < 9; ++i) {
                    for (int j = 0; j < 9; ++j) {
                        result += weights[i] * weights[j] * e[INDEX3(comp,i,j,numComp,9)];
                    }
                }
                integrals[comp] += result;
            }
        }
    }
    for (int comp = 0; comp < numComp; ++comp) {
        integrals[comp] *= volume_product;
    }
}

void Rectangle::integral_order9(std::vector<double>& integrals, const escript::Data& arg) const {
    const double weights[] = {0.0222222222222, 0.133305990851, 0.224889342063, 0.29204268368, 0.327539761184, 0.327539761184, 0.29204268368, 0.224889342063, 0.133305990851, 0.0222222222222};
    const int numComp = arg.getDataPointSize();
    const double volume_product = 0.25*m_dx[0]*m_dx[1];
    for (int ei = 0; ei < m_NE[1]; ++ei) {
        for (int ej = 0; ej < m_NE[0]; ++ej) {
            const double *e = arg.getSampleDataRO(INDEX2(ej,ei,m_NE[0]));
            double result = 0;
            for (int comp = 0; comp < numComp; ++comp) {
                for (int i = 0; i < 10; ++i) {
                    for (int j = 0; j < 10; ++j) {
                        result += weights[i] * weights[j] * e[INDEX3(comp,i,j,numComp,10)];
                    }
                }
                integrals[comp] += result;
            }
        }
    }
    for (int comp = 0; comp < numComp; ++comp) {
        integrals[comp] *= volume_product;
    }
}

void Rectangle::integral_order10(std::vector<double>& integrals, const escript::Data& arg) const {
    const double weights[] = {0.0181818181818, 0.109612273267, 0.18716988178, 0.248048104264, 0.286879124779, 0.300217595456, 0.286879124779, 0.248048104264, 0.18716988178, 0.109612273267, 0.0181818181818};
    const int numComp = arg.getDataPointSize();
    const double volume_product = 0.25*m_dx[0]*m_dx[1];
    for (int ei = 0; ei < m_NE[1]; ++ei) {
        for (int ej = 0; ej < m_NE[0]; ++ej) {
            const double *e = arg.getSampleDataRO(INDEX2(ej,ei,m_NE[0]));
            double result = 0;
            for (int comp = 0; comp < numComp; ++comp) {
                for (int i = 0; i < 11; ++i) {
                    for (int j = 0; j < 11; ++j) {
                        result += weights[i] * weights[j] * e[INDEX3(comp,i,j,numComp,11)];
                    }
                }
                integrals[comp] += result;
            }
        }
    }
    for (int comp = 0; comp < numComp; ++comp) {
        integrals[comp] *= volume_product;
    }
}

}
