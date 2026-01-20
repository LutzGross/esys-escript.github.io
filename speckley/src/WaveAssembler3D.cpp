
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

#include <speckley/WaveAssembler3D.h>
#include <speckley/domainhelpers.h>

#include <escript/index.h>

const double all_weights[][11] = {
    {0.333333333333, 1.33333333333, 0.333333333333},
    {0.166666666667, 0.833333333333, 0.833333333333, 0.166666666667},
    {0.1, 0.544444444444, 0.711111111111, 0.544444444444, 0.1},
    {0.0666666666667, 0.378474956298, 0.554858377035, 0.554858377035, 0.378474956298, 0.0666666666667},
    {0.047619047619, 0.276826047362, 0.43174538121, 0.487619047619, 0.43174538121, 0.276826047362, 0.047619047619},
    {0.0357142857143, 0.210704227144, 0.341122692484, 0.412458794659, 0.412458794659, 0.341122692484,0.210704227144, 0.0357142857143},
    {0.0277777777778, 0.165495361561, 0.2745387125, 0.346428510973, 0.371519274376, 0.346428510973, 0.2745387125, 0.165495361561, 0.0277777777778},
    {0.0222222222222, 0.133305990851, 0.224889342063, 0.29204268368, 0.327539761184, 0.327539761184, 0.29204268368, 0.224889342063, 0.133305990851, 0.0222222222222},
    {0.0181818181818, 0.109612273267, 0.18716988178, 0.248048104264, 0.286879124779, 0.300217595456, 0.286879124779, 0.248048104264, 0.18716988178, 0.109612273267, 0.0181818181818}
};

const double all_lagrange_derivs[][11][11] = {
    { // order 2
        {-1.50000000000000, -0.500000000000000, 0.500000000000000},
        {2.00000000000000, 0, -2.00000000000000},
        {-0.500000000000000, 0.500000000000000, 1.50000000000000}
    }, { //order 3
        {-3.00000000000000, -0.809016994374948, 0.309016994374948, -0.500000000000000},
        {4.04508497187474, 4.44089209850063e-16, -1.11803398874990, 1.54508497187474},
        {-1.54508497187474, 1.11803398874989, 2.22044604925031e-16, -4.04508497187474},
        {0.500000000000000, -0.309016994374947, 0.809016994374948, 3.00000000000000}
    }, { //order 4
        {-4.99999999999999, -1.24099025303098, 0.374999999999999, -0.259009746969017, 0.499999999999999},
        {6.75650248872424, -6.66133814775094e-15, -1.33658457769545, 0.763762615825974, -1.41016417794243},
        {-2.66666666666667, 1.74574312188794, 1.44328993201270e-15, -1.74574312188794, 2.66666666666667},
        {1.41016417794243, -0.763762615825974, 1.33658457769545, 1.66533453693773e-15, -6.75650248872424},
        {-0.500000000000001, 0.259009746969017, -0.375000000000000, 1.24099025303098, 5.00000000000000}
    }, { //order 5
        {-7.50000000000002, -1.78636494833911, 0.484951047853572, -0.269700610832040, 0.237781177984232, -0.500000000000002},
        {10.1414159363197, 2.13162820728030e-14, -1.72125695283023, 0.786356672223240, -0.653547507429800, 1.34991331419049},
        {-4.03618727030532, 2.52342677742945, -4.66293670342566e-15, -1.75296196636786, 1.15282815853593, -2.24468464817616},
        {2.24468464817616, -1.15282815853593, 1.75296196636787, -1.77635683940025e-15, -2.52342677742946, 4.03618727030535},
        {-1.34991331419048, 0.653547507429800, -0.786356672223242, 1.72125695283023, 2.22044604925031e-15, -10.1414159363197},
        {0.499999999999998, -0.237781177984231, 0.269700610832039, -0.484951047853569, 1.78636494833909, 7.50000000000000}
    }, { //order 6
        {-10.5000000000000, -2.44292601424426, 0.625256665515336, -0.312499999999997, 0.226099400942572, -0.226611870395444, 0.500000000000001},
        {14.2015766029198, -4.17443857259059e-14, -2.21580428316997, 0.907544471268819, -0.616390835517577, 0.602247179635785, -1.31737343570244},
        {-5.66898522554555, 3.45582821429430, 3.10862446895044e-15, -2.00696924058875, 1.06644190400637, -0.961339797288714, 2.04996481307676},
        {3.20000000000003, -1.59860668809837, 2.26669808708600, 1.33226762955019e-15, -2.26669808708599, 1.59860668809837, -3.20000000000003},
        {-2.04996481307676, 0.961339797288717, -1.06644190400638, 2.00696924058876, -1.77635683940025e-14, -3.45582821429431, 5.66898522554558},
        {1.31737343570245, -0.602247179635788, 0.616390835517580, -0.907544471268822, 2.21580428316998, 6.76125821996720e-14, -14.2015766029198},
        {-0.500000000000000, 0.226611870395444, -0.226099400942572, 0.312499999999997, -0.625256665515335, 2.44292601424425, 10.5000000000000}
    }, { //order 7
        {-13.9999999999999, -3.20991570300295, 0.792476681320508, -0.372150435728592, 0.243330712723790, -0.203284568900591, 0.219957514771299, -0.499999999999980},
        {18.9375986071174, -6.23945339839338e-14, -2.80647579473643, 1.07894468879045, -0.661157350900312, 0.537039586157660, -0.573565414940254, 1.29768738832019},
        {-7.56928981934855, 4.54358506456659, 1.17683640610267e-14, -2.37818723351551, 1.13535801688112, -0.845022556506511, 0.869448098331479, -1.94165942554406},
        {4.29790816426521, -2.11206121431455, 2.87551740597250, 3.77475828372553e-15, -2.38892435915824, 1.37278583180603, -1.29423205091348, 2.81018898925786},
        {-2.81018898925796, 1.29423205091350, -1.37278583180602, 2.38892435915823, 7.77156117237610e-15, -2.87551740597249, 2.11206121431450, -4.29790816426503},
        {1.94165942554413, -0.869448098331493, 0.845022556506508, -1.13535801688111, 2.37818723351550, -2.44249065417534e-14, -4.54358506456649, 7.56928981934828},
        {-1.29768738832026, 0.573565414940274, -0.537039586157669, 0.661157350900321, -1.07894468879047, 2.80647579473648, -1.71085368094737e-13, -18.9375986071174},
        {0.500000000000021, -0.219957514771312, 0.203284568900599, -0.243330712723800, 0.372150435728609, -0.792476681320546, 3.20991570300311, 14.0000000000002}
    }, { //order 8
        {-18.0000000000010, -4.08701370203454, 0.985360090074639, -0.444613449281139, 0.273437500000029, -0.207734512035617, 0.189655591978376, -0.215654018702531, 0.500000000000095},
        {24.3497451715930, 1.34892097491957e-12, -3.48835875343438, 1.28796075006388, -0.741782397916244, 0.547300160534042, -0.492350938315503, 0.555704981283736, -1.28483063269969},
        {-9.73870165721010, 5.78680581663678, -3.21298543326520e-13, -2.83445891207935, 1.26941308635811, -0.855726185092640, 0.738349277190360, -0.816756381741392, 1.87444087344708},
        {5.54496390694879, -2.69606544031400, 3.57668094012577, -3.10862446895044e-15, -2.65931021757391, 1.37696489376050, -1.07980381128263, 1.14565373845518, -2.59074567655957},
        {-3.65714285714248, 1.66522164500537, -1.71783215719513, 2.85191596846290, -1.06581410364015e-14, -2.85191596846287, 1.71783215719506, -1.66522164500546, 3.65714285714316},
        {2.59074567655910, -1.14565373845513, 1.07980381128268, -1.37696489376052, 2.65931021757393, -3.71924713249427e-14, -3.57668094012566, 2.69606544031420, -5.54496390694988},
        {-1.87444087344680, 0.816756381741386, -0.738349277190418, 0.855726185092682, -1.26941308635816, 2.83445891207943, 9.27036225562006e-14, -5.78680581663756, 9.73870165721225},
        {1.28483063269940, -0.555704981283693, 0.492350938315507, -0.547300160534031, 0.741782397916225, -1.28796075006385, 3.48835875343430, 5.71542813077031e-13, -24.3497451715928},
        {-0.499999999999903, 0.215654018702479, -0.189655591978347, 0.207734512035579, -0.273437499999975, 0.444613449281046, -0.985360090074401, 4.08701370203326, 17.9999999999994}
    }, { //order 9
        {-22.4999999999988, -5.07406470297709, 1.20335199285206, -0.528369376820220, 0.312047255608382, -0.223527944742433, 0.186645789393719, -0.180786585489230, 0.212702758009187, -0.500000000000077},
        {30.4381450292820, -1.52677870346452e-12, -4.25929735496529, 1.52990263818163, -0.845813573406436, 0.588082143045176, -0.483462326333953, 0.464274958908154, -0.543753738235757, 1.27595483609299},
        {-12.1779467074315, 7.18550286970643, 3.27293747659496e-13, -3.36412586829791, 1.44485031560171, -0.916555180336469, 0.721237312721631, -0.676797087196100, 0.783239293138005, -1.82956393190377},
        {6.94378848513465, -3.35166386274684, 4.36867455701003, 3.17523785042795e-14, -3.02021795819936, 1.46805550939000, -1.04618936550250, 0.936603213139437, -1.05915446364554, 2.45288417544331},
        {-4.59935476110357, 2.07820799403643, -2.10435017941307, 3.38731810120242, 2.22044604925031e-16, -3.02518848775198, 1.64649408398706, -1.33491548387823, 1.44494844875159, -3.29464303375000},
        {3.29464303374949, -1.44494844875146, 1.33491548387820, -1.64649408398705, 3.02518848775197, 6.66133814775094e-15, -3.38731810120245, 2.10435017941312, -2.07820799403661, 4.59935476110425},
        {-2.45288417544291, 1.05915446364544, -0.936603213139410, 1.04618936550249, -1.46805550938999, 3.02021795819934, -1.55431223447522e-15, -4.36867455701010, 3.35166386274711, -6.94378848513561},
        {1.82956393190345, -0.783239293137921, 0.676797087196070, -0.721237312721610, 0.916555180336448, -1.44485031560168, 3.36412586829786, -2.10609307771392e-13, -7.18550286970689, 12.1779467074328},
        {-1.27595483609268, 0.543753738235664, -0.464274958908104, 0.483462326333909, -0.588082143045126, 0.845813573406365, -1.52990263818151, 4.25929735496501, 2.64499533386697e-12, -30.4381450292815},
        {0.499999999999919, -0.212702758009135, 0.180786585489197, -0.186645789393687, 0.223527944742396, -0.312047255608330, 0.528369376820133, -1.20335199285186, 5.07406470297626, 22.4999999999976}
    }, { //order 10
        {-27.4999999999896, -6.17098569730879, 1.44617248279108, -0.622725214738251, 0.357476373116252, -0.246093749999837, 0.194287668796289, -0.172970108511720, 0.174657862947473, -0.210587346312973, 0.500000000000200},
        {37.2028673819635, -1.45093936865237e-11, -5.11821182477732, 1.80264679987985, -0.968487138802308, 0.646939963832195, -0.502643313012917, 0.443395636437341, -0.445313527290911, 0.535331085929629, -1.26956267628665},
        {-14.8873962951139, 8.73967005352577, 4.30699920173083e-12, -3.96199657704633, 1.65273663422738, -1.00650540857701, 0.747734824688410, -0.643586209360019, 0.637362056418227, -0.760400982244270, 1.79798803579712},
        {8.49561949463899, -4.07931619371279, 5.25066175545241, -4.06341627012807e-13, -3.45061471617355, 1.60812790372552, -1.07998725049569, 0.884587314556481, -0.852916813558380, 1.00338624297419, -2.35976991309107},
        {-5.64038799768972, 2.53474117868649, -2.53318340860873, 3.99079578802906, 9.19264664389630e-14, -3.30517685337838, 1.69057057046883, -1.24905529156928, 1.14606853428002, -1.31552671443958, 3.06553920088515},
        {4.06349206349476, -1.77190705526895, 1.61441910793459, -1.94634945457108, 3.45885134807703, 3.86357612569554e-14, -3.45885134807708, 1.94634945457117, -1.61441910793481, 1.77190705526946, -4.06349206349634},
        {-3.06553920088389, 1.31552671443917, -1.14606853427984, 1.24905529156920, -1.69057057046877, 3.30517685337831, -4.28546087505310e-14, -3.99079578802918, 2.53318340860904, -2.53474117868717, 5.64038799769177},
        {2.35976991309003, -1.00338624297385, 0.852916813558223, -0.884587314556402, 1.07998725049562, -1.60812790372545, 3.45061471617347, 5.73430192218893e-13, -5.25066175545293, 4.07931619371373, -8.49561949464169},
        {-1.79798803579616, 0.760400982243941, -0.637362056418051, 0.643586209359903, -0.747734824688294, 1.00650540857687, -1.65273663422718, 3.96199657704598, -3.31623617455534e-12, -8.73967005352657, 14.8873962951164},
        {1.26956267628575, -0.535331085929307, 0.445313527290713, -0.443395636437185, 0.502643313012754, -0.646939963831994, 0.968487138802021, -1.80264679987934, 5.11821182477613, 1.64738223062955e-11, -37.2028673819613},
        {-0.499999999999790, 0.210587346312823, -0.174657862947375, 0.172970108511639, -0.194287668796203, 0.246093749999731, -0.357476373116102, 0.622725214737994, -1.44617248279054, 6.17098569730708, 27.4999999999864}
    }
};

using escript::AbstractSystemMatrix;
using escript::Data;

namespace speckley {

void WaveAssembler3D::collateFunctionSpaceTypes(std::vector<int>& fsTypes,
                                   const DataMap& coefs) const
{
    if (isNotEmpty("D", coefs))
        fsTypes.push_back(coefs.find("D")->second.getFunctionSpace().getTypeCode());
    if (isNotEmpty("du", coefs))
        fsTypes.push_back(coefs.find("du")->second.getFunctionSpace().getTypeCode());
    if (isNotEmpty("Y", coefs))
        fsTypes.push_back(coefs.find("Y")->second.getFunctionSpace().getTypeCode());
}

void WaveAssembler3D::assemblePDESingle(AbstractSystemMatrix* mat,
                                        Data& rhs, const DataMap& coefs) const
{
    const Data& A = unpackData("A", coefs);
    const Data& B = unpackData("B", coefs);
    const Data& C = unpackData("C", coefs);
    const Data& D = unpackData("D", coefs);
    const Data& du = unpackData("du", coefs);
    const Data& Y = unpackData("Y", coefs);
    assemblePDESingle(mat, rhs, A, B, C, D, du, Y);

}

void WaveAssembler3D::assemblePDEBoundarySingle(AbstractSystemMatrix* mat,
                                        Data& rhs, const DataMap& coefs) const 
{
    const Data& d = unpackData("d", coefs);
    const Data& y = unpackData("y", coefs);
    assemblePDEBoundarySingle(mat, rhs, d, y);
}

void WaveAssembler3D::assemblePDESingleReduced(AbstractSystemMatrix* mat,
                                        Data& rhs, const DataMap& coefs) const
{
    const Data& A = unpackData("A", coefs);
    const Data& B = unpackData("B", coefs);
    const Data& C = unpackData("C", coefs);
    const Data& D = unpackData("D", coefs);
    const Data& du = unpackData("du", coefs);
    const Data& Y = unpackData("Y", coefs);
    assemblePDESingleReduced(mat, rhs, A, B, C, D, du, Y);
}

void WaveAssembler3D::assemblePDEBoundarySingleReduced(
                                        AbstractSystemMatrix* mat,
                                        Data& rhs, const DataMap& coefs) const
{
    const Data& d = unpackData("d", coefs);
    const Data& y = unpackData("y", coefs);
    assemblePDEBoundarySingleReduced(mat, rhs, d, y);
}

void WaveAssembler3D::assemblePDESystem(AbstractSystemMatrix* mat,
                                        Data& rhs, const DataMap& coefs) const
{
    const Data& A = unpackData("A", coefs);
    const Data& B = unpackData("B", coefs);
    const Data& C = unpackData("C", coefs);
    const Data& D = unpackData("D", coefs);
    const Data& du = unpackData("du", coefs);
    const Data& Y = unpackData("Y", coefs);
    assemblePDESystem(mat, rhs, A, B, C, D, du, Y);
}

void WaveAssembler3D::assemblePDEBoundarySystem(AbstractSystemMatrix* mat,
                                        Data& rhs, const DataMap& coefs) const
{
    const Data& d = unpackData("d", coefs);
    const Data& y = unpackData("y", coefs);
    assemblePDEBoundarySystem(mat, rhs, d, y);
}

void WaveAssembler3D::assemblePDESystemReduced(AbstractSystemMatrix* mat,
                                        Data& rhs, const DataMap& coefs) const
{
    const Data& A = unpackData("A", coefs);
    const Data& B = unpackData("B", coefs);
    const Data& C = unpackData("C", coefs);
    const Data& D = unpackData("D", coefs);
    const Data& du = unpackData("du", coefs);
    const Data& Y = unpackData("Y", coefs);
    assemblePDESystemReduced(mat, rhs, A, B, C, D, du, Y);
}

void WaveAssembler3D::assemblePDEBoundarySystemReduced(
                                        AbstractSystemMatrix* mat,
                                        Data& rhs, const DataMap& coefs) const
{
    const Data& d = unpackData("d", coefs);
    const Data& y = unpackData("y", coefs);
    assemblePDEBoundarySystemReduced(mat, rhs, d, y);
}


void WaveAssembler3D::assemblePDESystem(AbstractSystemMatrix* mat, Data& rhs,
                                           const Data& A, const Data& B,
                                           const Data& C, const Data& D,
                                           const Data& du, const Data& Y) const
{
    int order = domain->m_order;
    const double *weights = all_weights[order-2];
    const double volume_product = m_dx[0]*m_dx[1]*m_dx[2]/8.;
    const int NE0 = m_NE[0];
    const int NE1 = m_NE[1];
    const int NE2 = m_NE[2];
    const int quads = order + 1;
    const int max_x = m_NN[0];
    const int max_y = m_NN[1];
    dim_t numComp;
    if (!mat)
        numComp=(rhs.isEmpty() ? 1 : rhs.getDataPointSize());
    else {
        numComp=mat->getColumnBlockSize();
    }
    
    rhs.requireWrite();
    int size = 1, mid = 0;
    if (!Y.isEmpty()) {
        size = Y.getDataPointSize();
        mid = size/2-1;
        if (mid < 0)
            mid = 0;
    }
    const index_t y_indices[3] = {0,mid,size-1};
    size = 1;
    if (!D.isEmpty()) {
        size = D.getDataPointSize();
        mid = size/2-1;
        if (mid < 0)
            mid = 0;
    }
    const index_t d_indices[3] = {0,mid,size-1};

    if (!D.isEmpty() && (!du.isEmpty() || !Y.isEmpty()))
        throw SpeckleyException("assemblers can't deal with adding both lhs and rhs right now");
    
    for (dim_t colouring = 0; colouring < 2; colouring++) {
#pragma omp parallel for
        for (dim_t ez = colouring; ez < NE2; ez += 2) {
            for (dim_t ey = 0; ey < NE1; ey++) {
                for (dim_t ex = 0; ex < NE0; ex++) {
                    const index_t e_index = INDEX3(ex,ey,ez,NE0,NE1); //element number for Elements
                    const index_t start = order * (INDEX3(ex, ey, ez, max_x, max_y)); //node start for Nodes
                    if (!D.isEmpty()) {
                        const double *e = D.getSampleDataRO(e_index);
                        if (D.actsExpanded()) {
                            for (index_t comp = 0; comp < numComp; comp++) {
                                for (short qz = 0; qz < quads; qz++) {
                                    for (short qy = 0; qy < quads; qy++) {
                                        for (short qx = 0; qx < quads; qx++) {
                                            double *out = rhs.getSampleDataRW(start + INDEX3(qx,qy,qz,max_x,max_y));
/* D */ out[comp] += volume_product * weights[qx] * weights[qy] * weights[qz]
            * e[INDEX4(comp, qx, qy, qz, numComp, quads, quads)]; //D(i),qx,qy,qz
                                        }
                                    }
                                }
                            }
                        } else { //constant
                            for (index_t comp = 0; comp < numComp; comp++) {
                                for (short qz = 0; qz < quads; qz++) {
                                    for (short qy = 0; qy < quads; qy++) {
                                        for (short qx = 0; qx < quads; qx++) {
                                            double *out = rhs.getSampleDataRW(start + INDEX3(qx,qy,qz,max_x,max_y));
/* D */ out[comp] += volume_product * weights[qx] * weights[qy] * weights[qz]
                        * e[d_indices[comp]]; //D(i)
                                        }
                                    }
                                }
                            }
                        }
                    }
                    
                    if (!du.isEmpty()) {
                        if (isHTI) {
                            const double *du_p = du.getSampleDataRO(e_index);
                            const double c11_v = -c11.getSampleDataRO(e_index)[0];
                            const double c13_v = -c13.getSampleDataRO(e_index)[0];
                            const double c23_v = -c23.getSampleDataRO(e_index)[0];
                            const double c33_v = -c33.getSampleDataRO(e_index)[0];
                            const double c44_v = -c44.getSampleDataRO(e_index)[0];
                            const double c66_v = -c66.getSampleDataRO(e_index)[0];
                            if (du.actsExpanded()) {
                                for (short qz = 0; qz < quads; qz++) {
                                    for (short qy = 0; qy < quads; qy++) {
                                        for (short qx = 0; qx < quads; qx++) {
                                            double *out = rhs.getSampleDataRW(start + INDEX3(qx,qy,qz,max_x,max_y));
                                            double res_a[3] = {0}, res_b[3] = {0}, res_c[3] = {0};
#define SI(_x_,_y_,_qx_,_qy_,_qz_) INDEX5((_x_),(_y_),(_qx_),(_qy_),(_qz_),3,3,quads,quads)
                                            for (short k = 0; k < quads; k++) {
                                                res_a[0] += weights[k] * all_lagrange_derivs[order-2][qx][k]
                                                        * (c11_v*du_p[SI(0,0,k,qy,qz)]+c13_v*(du_p[SI(1,1,k,qy,qz)] + du_p[SI(2,2,k,qy,qz)]));
                                                res_b[0] += weights[k] * all_lagrange_derivs[order-2][qy][k] 
                                                        * c66_v*(du_p[SI(0,1,qx,k,qz)] + du_p[SI(1,0,qx,k,qz)]);
                                                res_c[0] += weights[k] * all_lagrange_derivs[order-2][qz][k]
                                                        * c66_v*(du_p[SI(0,2,qx,qy,k)] + du_p[SI(2,0,qx,qy,k)]);
                                                        
                                                res_a[1] += weights[k] * all_lagrange_derivs[order-2][qx][k]
                                                        * c66_v*(du_p[SI(1,0,k,qy,qz)] + du_p[SI(0,1,k,qy,qz)]);
                                                res_b[1] += weights[k] * all_lagrange_derivs[order-2][qy][k] 
                                                        * (c13_v*du_p[SI(0,0,qx,k,qz)] + c33_v*du_p[SI(1,1,qx,k,qz)] + c23_v*du_p[SI(2,2,qx,k,qz)]);
                                                res_c[1] += weights[k] * all_lagrange_derivs[order-2][qz][k]
                                                        * c44_v*(du_p[SI(1,2,qx,qy,k)] + du_p[SI(2,1,qx,qy,k)]);
                                                        
                                                res_a[2] += weights[k] * all_lagrange_derivs[order-2][qx][k]
                                                        * c66_v*(du_p[SI(2,0,k,qy,qz)] + du_p[SI(0,2,k,qy,qz)]);
                                                res_b[2] += weights[k] * all_lagrange_derivs[order-2][qy][k] 
                                                        * c44_v*(du_p[SI(2,1,qx,k,qz)] + du_p[SI(1,2,qx,k,qz)]);
                                                res_c[2] += weights[k] * all_lagrange_derivs[order-2][qz][k]
                                                        * (c13_v*du_p[SI(0,0,qx,qy,k)] +c23_v*du_p[SI(1,1,qx,qy,k)] + c33_v*du_p[SI(2,2,qx,qy,k)]);
                                            }
#undef SI
                                            for (short i = 0; i < 3; i++) {
                                                out[i] += 2 * volume_product * (
                                                        res_a[i] * weights[qy] * weights[qz] / m_dx[0] 
                                                      + res_b[i] * weights[qx] * weights[qz] / m_dx[1] 
                                                      + res_c[i] * weights[qx] * weights[qy] / m_dx[2]
                                                    );
                                            }
                                        }
                                    }
                                }
                            } else { //constant
                                for (short qz = 0; qz < quads; qz++) {
                                    for (short qy = 0; qy < quads; qy++) {
                                        for (short qx = 0; qx < quads; qx++) {
                                            double *out = rhs.getSampleDataRW(start + INDEX3(qx,qy,qz,max_x,max_y));
                                            double res_a[3] = {0}, res_b[3] = {0}, res_c[3] = {0};
#define SI(_x_,_y_) INDEX2((_x_),(_y_),3)
                                            for (short k = 0; k < quads; k++) {
                                                res_a[0] += weights[k] * all_lagrange_derivs[order-2][qx][k]
                                                        * (c11_v*du_p[SI(0,0)]+c13_v*(du_p[SI(1,1)] + du_p[SI(2,2)]));
                                                res_b[0] += weights[k] * all_lagrange_derivs[order-2][qy][k] 
                                                        * c66_v*(du_p[SI(0,1)] + du_p[SI(1,0)]);
                                                res_c[0] += weights[k] * all_lagrange_derivs[order-2][qz][k]
                                                        * c66_v*(du_p[SI(0,2)] + du_p[SI(2,0)]);
                                                        
                                                res_a[1] += weights[k] * all_lagrange_derivs[order-2][qx][k]
                                                        * c66_v*(du_p[SI(1,0)] + du_p[SI(0,1)]);
                                                res_b[1] += weights[k] * all_lagrange_derivs[order-2][qy][k] 
                                                        * (c13_v*du_p[SI(0,0)] + c33_v*du_p[SI(1,1)] + c23_v*du_p[SI(2,2)]);
                                                res_c[1] += weights[k] * all_lagrange_derivs[order-2][qz][k]
                                                        * c44_v*(du_p[SI(1,2)] + du_p[SI(2,1)]);
                                                        
                                                res_a[2] += weights[k] * all_lagrange_derivs[order-2][qx][k]
                                                        * c66_v*(du_p[SI(2,0)] + du_p[SI(0,2)]);
                                                res_b[2] += weights[k] * all_lagrange_derivs[order-2][qy][k] 
                                                        * c44_v*(du_p[SI(2,1)] + du_p[SI(1,2)]);
                                                res_c[2] += weights[k] * all_lagrange_derivs[order-2][qz][k]
                                                        * (c13_v*du_p[SI(0,0)] +c23_v*du_p[SI(1,1)] + c33_v*du_p[SI(2,2)]);
                                            }
#undef SI
                                            for (short comp = 0; comp < 3; comp++) {
                                                out[comp] += 2 * volume_product * (
                                                    res_a[comp] * weights[qy] * weights[qz] / m_dx[0] 
                                                  + res_b[comp] * weights[qx] * weights[qz] / m_dx[1] 
                                                  + res_c[comp] * weights[qx] * weights[qy] / m_dx[2]
                                                );
                                            }
                                        }
                                    }
                                }
                            }
                        } else { //VTI
                            const double *du_p = du.getSampleDataRO(e_index);
                            const double c11_v = -c11.getSampleDataRO(e_index)[0];
                            const double c12_v = -c12.getSampleDataRO(e_index)[0];
                            const double c13_v = -c13.getSampleDataRO(e_index)[0];
                            const double c33_v = -c33.getSampleDataRO(e_index)[0];
                            const double c44_v = -c44.getSampleDataRO(e_index)[0];
                            const double c66_v = -c66.getSampleDataRO(e_index)[0];
                            if (du.actsExpanded()) {
                                for (short qz = 0; qz < quads; qz++) {
                                    for (short qy = 0; qy < quads; qy++) {
                                        for (short qx = 0; qx < quads; qx++) {
                                            double *out = rhs.getSampleDataRW(start + INDEX3(qx,qy,qz,max_x,max_y));
                                            double res_a[3] = {0}, res_b[3] = {0}, res_c[3] = {0};
#define SI(_x_,_y_,_qx_,_qy_,_qz_) INDEX5((_x_),(_y_),(_qx_),(_qy_),(_qz_),3,3,quads,quads)
                                            for (short k = 0; k < quads; k++) {
                                                res_a[0] += weights[k] * all_lagrange_derivs[order-2][qx][k]
                                                        * (c11_v*du_p[SI(0,0,k,qy,qz)]+c12_v*du_p[SI(1,1,k,qy,qz)] + c13_v*du_p[SI(2,2,k,qy,qz)]);
                                                res_b[0] += weights[k] * all_lagrange_derivs[order-2][qy][k] 
                                                        * c66_v*(du_p[SI(0,1,qx,k,qz)] + du_p[SI(1,0,qx,k,qz)]);
                                                res_c[0] += weights[k] * all_lagrange_derivs[order-2][qz][k]
                                                        * c44_v*(du_p[SI(0,2,qx,qy,k)] + du_p[SI(2,0,qx,qy,k)]);
                                                        
                                                res_a[1] += weights[k] * all_lagrange_derivs[order-2][qx][k]
                                                        * c66_v*(du_p[SI(1,0,k,qy,qz)] + du_p[SI(0,1,k,qy,qz)]);
                                                res_b[1] += weights[k] * all_lagrange_derivs[order-2][qy][k] 
                                                        * (c12_v*du_p[SI(0,0,qx,k,qz)] + c11_v*du_p[SI(1,1,qx,k,qz)] + c13_v*du_p[SI(2,2,qx,k,qz)]);
                                                res_c[1] += weights[k] * all_lagrange_derivs[order-2][qz][k]
                                                        * c44_v*(du_p[SI(1,2,qx,qy,k)] + du_p[SI(2,1,qx,qy,k)]);
                                                        
                                                res_a[2] += weights[k] * all_lagrange_derivs[order-2][qx][k]
                                                        * c44_v*(du_p[SI(2,0,k,qy,qz)] + du_p[SI(0,2,k,qy,qz)]);
                                                res_b[2] += weights[k] * all_lagrange_derivs[order-2][qy][k] 
                                                        * c44_v*(du_p[SI(2,1,qx,k,qz)] + du_p[SI(1,2,qx,k,qz)]);
                                                res_c[2] += weights[k] * all_lagrange_derivs[order-2][qz][k]
                                                        * (c13_v*(du_p[SI(0,0,qx,qy,k)] + du_p[SI(1,1,qx,qy,k)]) + c33_v*du_p[SI(2,2,qx,qy,k)]);
                                            }
#undef SI
                                            for (short i = 0; i < 3; i++) {
                                                out[i] += 2 * volume_product * (
                                                        res_a[i] * weights[qy] * weights[qz] / m_dx[0] 
                                                      + res_b[i] * weights[qx] * weights[qz] / m_dx[1] 
                                                      + res_c[i] * weights[qx] * weights[qy] / m_dx[2]
                                                    );
                                            }
                                        }
                                    }
                                }
                            } else { //constant
                                for (short qz = 0; qz < quads; qz++) {
                                    for (short qy = 0; qy < quads; qy++) {
                                        for (short qx = 0; qx < quads; qx++) {
                                            double *out = rhs.getSampleDataRW(start + INDEX3(qx,qy,qz,max_x,max_y));
                                            double res_a[3] = {0}, res_b[3] = {0}, res_c[3] = {0};
#define SI(_x_,_y_) INDEX2((_x_),(_y_),3)
                                            for (short k = 0; k < quads; k++) {
                                                res_a[0] += weights[k] * all_lagrange_derivs[order-2][qx][k]
                                                        * (c11_v*du_p[SI(0,0)]+c12_v*du_p[SI(1,1)] + c13_v*du_p[SI(2,2)]);
                                                res_b[0] += weights[k] * all_lagrange_derivs[order-2][qy][k] 
                                                        * c66_v*(du_p[SI(0,1)] + du_p[SI(1,0)]);
                                                res_c[0] += weights[k] * all_lagrange_derivs[order-2][qz][k]
                                                        * c44_v*(du_p[SI(0,2)] + du_p[SI(2,0)]);
                                                        
                                                res_a[1] += weights[k] * all_lagrange_derivs[order-2][qx][k]
                                                        * c66_v*(du_p[SI(1,0)] + du_p[SI(0,1)]);
                                                res_b[1] += weights[k] * all_lagrange_derivs[order-2][qy][k] 
                                                        * (c12_v*du_p[SI(0,0)] + c11_v*du_p[SI(1,1)] + c13_v*du_p[SI(2,2)]);
                                                res_c[1] += weights[k] * all_lagrange_derivs[order-2][qz][k]
                                                        * c44_v*(du_p[SI(1,2)] + du_p[SI(2,1)]);
                                                        
                                                res_a[2] += weights[k] * all_lagrange_derivs[order-2][qx][k]
                                                        * c44_v*(du_p[SI(2,0)] + du_p[SI(0,2)]);
                                                res_b[2] += weights[k] * all_lagrange_derivs[order-2][qy][k] 
                                                        * c44_v*(du_p[SI(2,1)] + du_p[SI(1,2)]);
                                                res_c[2] += weights[k] * all_lagrange_derivs[order-2][qz][k]
                                                        * (c13_v*(du_p[SI(0,0)] + du_p[SI(1,1)]) + c33_v*du_p[SI(2,2)]);
                                            }
#undef SI
                                            for (short comp = 0; comp < 3; comp++) {
                                                out[comp] += 2 * volume_product * (
                                                    res_a[comp] * weights[qy] * weights[qz] / m_dx[0] 
                                                  + res_b[comp] * weights[qx] * weights[qz] / m_dx[1] 
                                                  + res_c[comp] * weights[qx] * weights[qy] / m_dx[2]
                                                );
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    
                    if (!Y.isEmpty()) {
                        const double *e = Y.getSampleDataRO(e_index);
                        if (Y.actsExpanded()) {
                            for (index_t comp = 0; comp < numComp; comp++) {
                                for (short qz = 0; qz < quads; qz++) {
                                    for (short qy = 0; qy < quads; qy++) {
                                        for (short qx = 0; qx < quads; qx++) {
                                            double *out = rhs.getSampleDataRW(start + INDEX3(qx,qy,qz,max_x,max_y));
/* Y */ out[comp] += volume_product * weights[qx] * weights[qy] * weights[qz]
                    * e[INDEX4(comp, qx, qy, qz, numComp, quads, quads)]; //Y(i),qx,qy,qz
                                        }
                                    }
                                }
                            }
                        } else { //constant
                            for (index_t comp = 0; comp < numComp; comp++) {
                                for (short qz = 0; qz < quads; qz++) {
                                    for (short qy = 0; qy < quads; qy++) {
                                        for (short qx = 0; qx < quads; qx++) {
                                            double *out = rhs.getSampleDataRW(start + INDEX3(qx,qy,qz,max_x,max_y));
/* Y */ out[comp] += volume_product * weights[qx] * weights[qy] * weights[qz]
                    * e[y_indices[comp]]; //Y(i)
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

//protected
void WaveAssembler3D::assemblePDESystemReduced(AbstractSystemMatrix* mat,
                                     Data& rhs, const Data& A, const Data& B,
                                     const Data& C, const Data& D,
                                     const Data& du, const Data& Y) const
{
    throw SpeckleyException("system reduced assemblers not implemented yet");
    
}

//protected
void WaveAssembler3D::assemblePDEBoundarySingle(AbstractSystemMatrix* mat,
                                Data& rhs, const Data& d, const Data& y) const
{
    throw SpeckleyException("boundary single assemblers not implemented yet");
}

//protected
void WaveAssembler3D::assemblePDEBoundarySingleReduced(
                                        AbstractSystemMatrix* mat, Data& rhs,
                                        const Data& d, const Data& y) const
{
    throw SpeckleyException("boundary single reduced assemblers not implemented yet");
}

void WaveAssembler3D::assemblePDEBoundarySystem(AbstractSystemMatrix* mat,
                               Data& rhs, const Data& d, const Data& y) const
{
    throw SpeckleyException("boundary system assemblers not implemented yet");
}

//protected
void WaveAssembler3D::assemblePDEBoundarySystemReduced(
                                        AbstractSystemMatrix* mat, Data& rhs,
                                        const Data& d, const Data& y) const
{
    throw SpeckleyException("boundary system reduced assemblers not implemented yet");
}

//protected
void WaveAssembler3D::assemblePDESingle(escript::AbstractSystemMatrix *mat,
        escript::Data& rhs, const escript::Data& A, const escript::Data& B,
        const escript::Data& C, const escript::Data& D,
        const escript::Data& du, const escript::Data& Y) const
{
    throw SpeckleyException("Wave assemblers do not support assemblePDESingle()");
}

//protected
void WaveAssembler3D::assemblePDESingleReduced(AbstractSystemMatrix* mat,
                                    Data& rhs, const Data& A, const Data& B,
                                    const Data& C, const Data& D,
                                    const Data& du, const Data& Y) const
{
    throw SpeckleyException("single reduced assemblers not implemented yet");
}


} // namespace speckley

