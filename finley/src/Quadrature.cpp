
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


/****************************************************************************

  Finley: quadrature schemes

*****************************************************************************/

#include "Quadrature.h"

#include <escript/index.h>

#define QUADNODES(_K_,_I_) quadNodes[INDEX2(_K_,_I_,DIM)]
#define QUADWEIGHTS(_I_) quadWeights[_I_]

namespace finley {

const QuadInfo QuadInfoList[] = {
    { PointQuad,"Point", 0, 1,  Quad_getNodesPoint, Quad_getNumNodesPoint, Quad_MacroPoint },
    { LineQuad, "Line",  1,  2, Quad_getNodesLine,  Quad_getNumNodesLine,  Quad_MacroLine },
    { TriQuad,   "Tri",  2,  3, Quad_getNodesTri,   Quad_getNumNodesTri,   Quad_MacroTri },
    { RecQuad,   "Rec",  2,  4, Quad_getNodesRec,   Quad_getNumNodesRec,   Quad_MacroRec },
    { TetQuad,   "Tet",  3,  4, Quad_getNodesTet,   Quad_getNumNodesTet,   Quad_MacroTet },
    { HexQuad,   "Hex",  3,  8, Quad_getNodesHex,   Quad_getNumNodesHex,   Quad_MacroHex },
    { NoQuad, "NoType",  0,  1, Quad_getNodesPoint, Quad_getNumNodesPoint, Quad_MacroPoint }
};

const QuadInfo* QuadInfo_getInfo(QuadTypeId id)
{
    int idx=0;
    const QuadInfo* out=NULL;
    while (QuadInfoList[idx].TypeId!=NoQuad && out==NULL) {
       if (QuadInfoList[idx].TypeId==id)
           out=&(QuadInfoList[idx]);
       idx++;
    }
    if (out==NULL) {
        throw escript::ValueError("QuadInfo_getInfo: cannot find requested quadrature scheme.");
    }
    return out;
}

/****************************************************************************/

/// get a quadrature scheme with numQuadNodes quadrature nodes for the tri
/// as a squeezed scheme on a quad [0,1]^2
void Quad_getNodesTri(int numQuadNodes, std::vector<double>& quadNodes, std::vector<double>& quadWeights)
{
#define DIM 2
    // the easy cases:

    if (numQuadNodes==1) {
        QUADNODES(0,0)=1./3.;
        QUADNODES(1,0)=1./3.;
        QUADWEIGHTS(0)=1./2.;
    } else if (numQuadNodes==3) {
        QUADNODES(0,0)=1./2.;
        QUADNODES(1,0)=0.;
        QUADWEIGHTS(0)=1./6.;
        QUADNODES(0,1)=0.;
        QUADNODES(1,1)=1./2.;
        QUADWEIGHTS(1)=1./6.;
        QUADNODES(0,2)=1./2.;
        QUADNODES(1,2)=1./2.;
        QUADWEIGHTS(2)=1./6.;
    } else if (numQuadNodes==4) {
        QUADNODES(0,0)=1./3.;
        QUADNODES(1,0)=1./3.;
        QUADWEIGHTS(0)=-27./96.;
        QUADNODES(0,1)=0.2;
        QUADNODES(1,1)=0.2;
        QUADWEIGHTS(1)=25./96.;
        QUADNODES(0,2)=0.6;
        QUADNODES(1,2)=0.2;
        QUADWEIGHTS(2)=25./96.;
        QUADNODES(0,3)=0.2;
        QUADNODES(1,3)=0.6;
        QUADWEIGHTS(3)=25./96.;
    } else if (numQuadNodes==6) {
        QUADWEIGHTS(0) = 0.109951743655322/2.;
        QUADWEIGHTS(1) = 0.109951743655322/2.;
        QUADWEIGHTS(2) = 0.109951743655322/2.;
        QUADWEIGHTS(3) = 0.223381589678011/2.;
        QUADWEIGHTS(4) = 0.223381589678011/2.;
        QUADWEIGHTS(5) = 0.223381589678011/2.;

        QUADNODES(0,0) = 0.816847572980459;
        QUADNODES(0,1) = 0.091576213509771;
        QUADNODES(0,2) = 0.091576213509771;
        QUADNODES(0,3) = 0.108103018168070;
        QUADNODES(0,4) = 0.445948490915965;
        QUADNODES(0,5) = 0.445948490915965;

        QUADNODES(1,0) = 0.091576213509771;
        QUADNODES(1,1) = 0.816847572980459;
        QUADNODES(1,2) = 0.091576213509771;
        QUADNODES(1,3) = 0.445948490915965;
        QUADNODES(1,4) = 0.108103018168070;
        QUADNODES(1,5) = 0.445948490915965;

    } else if (numQuadNodes==7) {
        QUADNODES(0,0) = 0.33333333333333333;
        QUADNODES(0,1) = 0.7974269853530872;
        QUADNODES(0,2) = 0.10128650732345633;
        QUADNODES(0,3) = 0.10128650732345633;
        QUADNODES(0,4) = 0.059715871789769809;
        QUADNODES(0,5) = 0.47014206410511505;
        QUADNODES(0,6) = 0.47014206410511505;

        QUADNODES(1,0) = 0.33333333333333333;
        QUADNODES(1,1) = 0.10128650732345633;
        QUADNODES(1,2) = 0.7974269853530872;
        QUADNODES(1,3) = 0.10128650732345633;
        QUADNODES(1,4) = 0.47014206410511505;
        QUADNODES(1,5) = 0.059715871789769809;
        QUADNODES(1,6) = 0.47014206410511505;

        QUADWEIGHTS(0) = 0.225/2.;
        QUADWEIGHTS(1) = 0.12593918054482717/2.;
        QUADWEIGHTS(2) = 0.12593918054482717/2.;
        QUADWEIGHTS(3) = 0.12593918054482717/2.;
        QUADWEIGHTS(4) = 0.13239415278850616/2.;
        QUADWEIGHTS(5) = 0.13239415278850616/2.;
        QUADWEIGHTS(6) = 0.13239415278850616/2.;

    } else if (numQuadNodes==12) {
        const double a = 0.873821971016996;
        const double b = 0.063089014491502;
        const double c = 0.501426509658179;
        const double d = 0.249286745170910;
        const double e = 0.636502499121399;
        const double f = 0.310352451033785;
        const double g = 0.053145049844816;

        const double u = 0.050844906370207/2.;
        const double v = 0.116786275726379/2.;
        const double w = 0.082851075618374/2.;

        QUADNODES(0,0) = a;
        QUADNODES(0,1) = b;
        QUADNODES(0,2) = b;
        QUADNODES(0,3) = c;
        QUADNODES(0,4) = d;
        QUADNODES(0,5) = d;
        QUADNODES(0,6) = e;
        QUADNODES(0,7) = e;
        QUADNODES(0,8) = f;
        QUADNODES(0,9) = f;
        QUADNODES(0,10) = g;
        QUADNODES(0,11) = g;

        QUADNODES(1,0) = b;
        QUADNODES(1,1) = a;
        QUADNODES(1,2) = b;
        QUADNODES(1,3) = d;
        QUADNODES(1,4) = c;
        QUADNODES(1,5) = d;
        QUADNODES(1,6) = f;
        QUADNODES(1,7) = g;
        QUADNODES(1,8) = e;
        QUADNODES(1,9) = g;
        QUADNODES(1,10) = e;
        QUADNODES(1,11) = f;

        QUADWEIGHTS(0) = u;
        QUADWEIGHTS(1) = u;
        QUADWEIGHTS(2) = u;
        QUADWEIGHTS(3) = v;
        QUADWEIGHTS(4) = v;
        QUADWEIGHTS(5) = v;
        QUADWEIGHTS(6) = w;
        QUADWEIGHTS(7) = w;
        QUADWEIGHTS(8) = w;
        QUADWEIGHTS(9) = w;
        QUADWEIGHTS(10) = w;
        QUADWEIGHTS(11) = w;

    } else if (numQuadNodes==13) {
        QUADWEIGHTS(0) =-0.149570044467670/2.;
        QUADWEIGHTS(1) = 0.175615257433204/2.;
        QUADWEIGHTS(2) = 0.175615257433204/2.;
        QUADWEIGHTS(3) = 0.175615257433204/2.;
        QUADWEIGHTS(4) = 0.053347235608839/2.;
        QUADWEIGHTS(5) = 0.053347235608839/2.;
        QUADWEIGHTS(6) = 0.053347235608839/2.;
        QUADWEIGHTS(7) = 0.077113760890257/2.;
        QUADWEIGHTS(8) = 0.077113760890257/2.;
        QUADWEIGHTS(9) = 0.077113760890257/2.;
        QUADWEIGHTS(10) = 0.077113760890257/2.;
        QUADWEIGHTS(11) = 0.077113760890257/2.;
        QUADWEIGHTS(12) = 0.077113760890257/2.;

        QUADNODES(0,0) = 0.3333333333333333;
        QUADNODES(0,1) = 0.479308067841923;
        QUADNODES(0,2) = 0.260345966079038;
        QUADNODES(0,3) = 0.260345966079038;
        QUADNODES(0,4) = 0.869739794195568;
        QUADNODES(0,5) = 0.065130102902216;
        QUADNODES(0,6) = 0.065130102902216;
        QUADNODES(0,7) = 0.638444188569809;
        QUADNODES(0,8) = 0.638444188569809;
        QUADNODES(0,9) = 0.048690315425316;
        QUADNODES(0,10) = 0.048690315425316;
        QUADNODES(0,11) = 0.312865496004875;
        QUADNODES(0,12) = 0.312865496004875;

        QUADNODES(1,0) = 0.3333333333333333;
        QUADNODES(1,1) = 0.260345966079038;
        QUADNODES(1,2) = 0.479308067841923;
        QUADNODES(1,3) = 0.260345966079038;
        QUADNODES(1,4) = 0.065130102902216;
        QUADNODES(1,5) = 0.869739794195568;
        QUADNODES(1,6) = 0.065130102902216;
        QUADNODES(1,7) = 0.048690315425316;
        QUADNODES(1,8) = 0.312865496004875;
        QUADNODES(1,9) = 0.638444188569809;
        QUADNODES(1,10) = 0.312865496004875;
        QUADNODES(1,11) = 0.638444188569809;
        QUADNODES(1,12) = 0.048690315425316;

    } else if (numQuadNodes==16) {
        QUADWEIGHTS(0) = 0.07215780;
        QUADWEIGHTS(1) = 0.04754582;
        QUADWEIGHTS(2) = 0.04754582;
        QUADWEIGHTS(3) = 0.04754582;
        QUADWEIGHTS(4) = 0.01622925;
        QUADWEIGHTS(5) = 0.01622925;
        QUADWEIGHTS(6) = 0.01622925;
        QUADWEIGHTS(7) = 0.05160869;
        QUADWEIGHTS(8) = 0.05160869;
        QUADWEIGHTS(9) = 0.05160869;
        QUADWEIGHTS(10) = 0.01361516;
        QUADWEIGHTS(11) = 0.01361516;
        QUADWEIGHTS(12) = 0.01361516;
        QUADWEIGHTS(13) = 0.01361516;
        QUADWEIGHTS(14) = 0.01361516;
        QUADWEIGHTS(15) = 0.01361516;

        QUADNODES(0,0) = 0.3333333;
        QUADNODES(0,1) = 0.08141482;
        QUADNODES(0,2) = 0.4592926;
        QUADNODES(0,3) = 0.4592926;
        QUADNODES(0,4) = 0.8989055;
        QUADNODES(0,5) = 0.05054723;
        QUADNODES(0,6) = 0.05054723;
        QUADNODES(0,7) = 0.6588614;
        QUADNODES(0,8) = 0.1705693;
        QUADNODES(0,9) = 0.1705693;
        QUADNODES(0,10) = 0.008394777;
        QUADNODES(0,11) = 0.008394777;
        QUADNODES(0,12) = 0.7284924;
        QUADNODES(0,13) = 0.7284924;
        QUADNODES(0,14) = 0.2631128;
        QUADNODES(0,15) = 0.2631128;

        QUADNODES(1,0) = 0.3333333;
        QUADNODES(1,1) = 0.4592926;
        QUADNODES(1,2) = 0.08141482;
        QUADNODES(1,3) = 0.4592926;
        QUADNODES(1,4) = 0.05054723;
        QUADNODES(1,5) = 0.8989055;
        QUADNODES(1,6) = 0.05054723;
        QUADNODES(1,7) = 0.1705693;
        QUADNODES(1,8) = 0.6588614;
        QUADNODES(1,9) = 0.1705693;
        QUADNODES(1,10) = 0.7284924;
        QUADNODES(1,11) = 0.2631128;
        QUADNODES(1,12) = 0.008394777;
        QUADNODES(1,13) = 0.2631128;
        QUADNODES(1,14) = 0.008394777;
        QUADNODES(1,15) = 0.7284924;

    } else if (numQuadNodes==19) {
        QUADWEIGHTS(0) = 0.04856790;
        QUADWEIGHTS(1) = 0.01566735;
        QUADWEIGHTS(2) = 0.01566735;
        QUADWEIGHTS(3) = 0.01566735;
        QUADWEIGHTS(4) = 0.03891377;
        QUADWEIGHTS(5) = 0.03891377;
        QUADWEIGHTS(6) = 0.03891377;
        QUADWEIGHTS(7) = 0.03982387;
        QUADWEIGHTS(8) = 0.03982387;
        QUADWEIGHTS(9) = 0.03982387;
        QUADWEIGHTS(10) = 0.01278884;
        QUADWEIGHTS(11) = 0.01278884;
        QUADWEIGHTS(12) = 0.01278884;
        QUADWEIGHTS(13) = 0.02164177;
        QUADWEIGHTS(14) = 0.02164177;
        QUADWEIGHTS(15) = 0.02164177;
        QUADWEIGHTS(16) = 0.02164177;
        QUADWEIGHTS(17) = 0.02164177;
        QUADWEIGHTS(18) = 0.02164177;

        QUADNODES(0,0) = 0.3333333;
        QUADNODES(0,1) = 0.02063496;
        QUADNODES(0,2) = 0.4896825;
        QUADNODES(0,3) = 0.4896825;
        QUADNODES(0,4) = 0.1258208;
        QUADNODES(0,5) = 0.4370896;
        QUADNODES(0,6) = 0.4370896;
        QUADNODES(0,7) = 0.6235929;
        QUADNODES(0,8) = 0.1882035;
        QUADNODES(0,9) = 0.1882035;
        QUADNODES(0,10) = 0.9105410;
        QUADNODES(0,11) = 0.04472951;
        QUADNODES(0,12) = 0.04472951;
        QUADNODES(0,13) = 0.03683841;
        QUADNODES(0,14) = 0.03683841;
        QUADNODES(0,15) = 0.7411986;
        QUADNODES(0,16) = 0.7411986;
        QUADNODES(0,17) = 0.2219630;
        QUADNODES(0,18) = 0.2219630;

        QUADNODES(1,0) = 0.3333333;
        QUADNODES(1,1) = 0.4896825;
        QUADNODES(1,2) = 0.02063496;
        QUADNODES(1,3) = 0.4896825;
        QUADNODES(1,4) = 0.4370896;
        QUADNODES(1,5) = 0.1258208;
        QUADNODES(1,6) = 0.4370896;
        QUADNODES(1,7) = 0.1882035;
        QUADNODES(1,8) = 0.6235929;
        QUADNODES(1,9) = 0.1882035;
        QUADNODES(1,10) = 0.04472951;
        QUADNODES(1,11) = 0.9105410;
        QUADNODES(1,12) = 0.04472951;
        QUADNODES(1,13) = 0.7411986;
        QUADNODES(1,14) = 0.2219630;
        QUADNODES(1,15) = 0.03683841;
        QUADNODES(1,16) = 0.2219630;
        QUADNODES(1,17) = 0.03683841;
        QUADNODES(1,18) = 0.7411986;

    } else {
        // get scheme on [0.1]^2
        Quad_getNodesRec(numQuadNodes, quadNodes, quadWeights);

        // squeeze it:
        for (int i=0; i<numQuadNodes; i++) {
            const double Q1=QUADNODES(0,i);
            const double Q2=QUADNODES(1,i);
            QUADWEIGHTS(i)=QUADWEIGHTS(i)*(1.-(1./2.)*(Q1+Q2));
            QUADNODES(0,i)=Q1*(1.-(1./2.)*Q2);
            QUADNODES(1,i)=Q2*(1.-(1./2.)*Q1);
        }
    } // numQuadNodes
#undef DIM
}

/// get a quadrature scheme with numQuadNodes quadrature nodes for the tet
/// as a squeezed scheme on a hex [0,1]^3
void Quad_getNodesTet(int numQuadNodes, std::vector<double>& quadNodes, std::vector<double>& quadWeights)
{
    const double alpha=0.58541019662496852;
    const double beta =0.1381966011250105;
#define DIM 3

    // the easy cases:
    if (numQuadNodes==1) {
        QUADNODES(0,0)=0.25;
        QUADNODES(1,0)=0.25;
        QUADNODES(2,0)=0.25;
        QUADWEIGHTS(0)=1./6.;
    } else if (numQuadNodes==4) {
        QUADNODES(0,0)=beta;
        QUADNODES(1,0)=beta;
        QUADNODES(2,0)=beta;
        QUADWEIGHTS(0)=1./24.;
        QUADNODES(0,1)=alpha;
        QUADNODES(1,1)=beta;
        QUADNODES(2,1)=beta;
        QUADWEIGHTS(1)=1./24.;
        QUADNODES(0,2)=beta;
        QUADNODES(1,2)=alpha;
        QUADNODES(2,2)=beta;
        QUADWEIGHTS(2)=1./24.;
        QUADNODES(0,3)=beta;
        QUADNODES(1,3)=beta;
        QUADNODES(2,3)=alpha;
        QUADWEIGHTS(3)=1./24.;
    } else if (numQuadNodes==5) {
        QUADNODES(0,0)=1./4.;
        QUADNODES(1,0)=1./4.;
        QUADNODES(2,0)=1./4.;
        QUADWEIGHTS(0)=-2./15.;
        QUADNODES(0,1)=1./6.;
        QUADNODES(1,1)=1./6.;
        QUADNODES(2,1)=1./6.;
        QUADWEIGHTS(1)=3./40.;
        QUADNODES(0,2)=1./2.;
        QUADNODES(1,2)=1./6.;
        QUADNODES(2,2)=1./6.;
        QUADWEIGHTS(2)=3./40.;
        QUADNODES(0,3)=1./6.;
        QUADNODES(1,3)=1./2.;
        QUADNODES(2,3)=1./6.;
        QUADWEIGHTS(3)=3./40.;
        QUADNODES(0,4)=1./6.;
        QUADNODES(1,4)=1./6.;
        QUADNODES(2,4)=1./2.;
        QUADWEIGHTS(4)=3./40.;

    } else if (numQuadNodes==11) {
        const double a =  0.25;
        const double b =  11.0/14.0;
        const double c =  1.0/14.0;
        const double d =  0.25 * (1.0 + sqrt ( 5.0 / 14.0 ) );
        const double e =  0.25 * (1.0 - sqrt ( 5.0 / 14.0 ) );
        const double f = -74.0 /  5625.0;
        const double g = 343.0 / 45000.0;
        const double h =  56.0 /  2250.0;

        QUADWEIGHTS(401-401) = f;
        QUADWEIGHTS(402-401) = g;
        QUADWEIGHTS(403-401) = g;
        QUADWEIGHTS(404-401) = g;
        QUADWEIGHTS(405-401) = g;
        QUADWEIGHTS(406-401) = h;
        QUADWEIGHTS(407-401) = h;
        QUADWEIGHTS(408-401) = h;
        QUADWEIGHTS(409-401) = h;
        QUADWEIGHTS(410-401) = h;
        QUADWEIGHTS(411-401) = h;

        QUADNODES(0,401-401) = a;
        QUADNODES(0,402-401) = b;
        QUADNODES(0,403-401) = c;
        QUADNODES(0,404-401) = c;
        QUADNODES(0,405-401) = c;
        QUADNODES(0,406-401) = d;
        QUADNODES(0,407-401) = d;
        QUADNODES(0,408-401) = d;
        QUADNODES(0,409-401) = e;
        QUADNODES(0,410-401) = e;
        QUADNODES(0,411-401) = e;

        QUADNODES(1,401-401) = a;
        QUADNODES(1,402-401) = c;
        QUADNODES(1,403-401) = b;
        QUADNODES(1,404-401) = c;
        QUADNODES(1,405-401) = c;
        QUADNODES(1,406-401) = d;
        QUADNODES(1,407-401) = e;
        QUADNODES(1,408-401) = e;
        QUADNODES(1,409-401) = d;
        QUADNODES(1,410-401) = d;
        QUADNODES(1,411-401) = e;

        QUADNODES(2,401-401) = a;
        QUADNODES(2,402-401) = c;
        QUADNODES(2,403-401) = c;
        QUADNODES(2,404-401) = b;
        QUADNODES(2,405-401) = c;
        QUADNODES(2,406-401) = e;
        QUADNODES(2,407-401) = d;
        QUADNODES(2,408-401) = e;
        QUADNODES(2,409-401) = d;
        QUADNODES(2,410-401) = e;
        QUADNODES(2,411-401) = d;

    } else if (numQuadNodes==15) {
        QUADWEIGHTS(412-412) = 0.019753086419753086;
        QUADWEIGHTS(413-412) = 0.01198951396316977;
        QUADWEIGHTS(414-412) = 0.01198951396316977;
        QUADWEIGHTS(415-412) = 0.01198951396316977;
        QUADWEIGHTS(416-412) = 0.01198951396316977;
        QUADWEIGHTS(417-412) = 0.011511367871045397;
        QUADWEIGHTS(418-412) = 0.011511367871045397;
        QUADWEIGHTS(419-412) = 0.011511367871045397;
        QUADWEIGHTS(420-412) = 0.011511367871045397;
        QUADWEIGHTS(421-412) = 0.0088183421516754845;
        QUADWEIGHTS(422-412) = 0.0088183421516754845;
        QUADWEIGHTS(423-412) = 0.0088183421516754845;
        QUADWEIGHTS(424-412) = 0.0088183421516754845;
        QUADWEIGHTS(425-412) = 0.0088183421516754845;
        QUADWEIGHTS(426-412) = 0.0088183421516754845;

        QUADNODES(0,412-412) = 0.2500000;
        QUADNODES(0,413-412) = 0.091971078052723032;
        QUADNODES(0,414-412) = 0.72408676584183096;
        QUADNODES(0,415-412) = 0.091971078052723032;
        QUADNODES(0,416-412) = 0.091971078052723032;
        QUADNODES(0,417-412) = 0.31979362782962989;
        QUADNODES(0,418-412) = 0.040619116511110234;
        QUADNODES(0,419-412) = 0.31979362782962989;
        QUADNODES(0,420-412) = 0.31979362782962989;
        QUADNODES(0,421-412) = 0.056350832689629149;
        QUADNODES(0,422-412) = 0.056350832689629149;
        QUADNODES(0,423-412) = 0.056350832689629149;
        QUADNODES(0,424-412) = 0.4436491673103708;
        QUADNODES(0,425-412) = 0.4436491673103708;
        QUADNODES(0,426-412) = 0.4436491673103708;

        QUADNODES(1,412-412) = 0.2500000;
        QUADNODES(1,413-412) = 0.091971078052723032;
        QUADNODES(1,414-412) = 0.091971078052723032;
        QUADNODES(1,415-412) = 0.72408676584183096;
        QUADNODES(1,416-412) = 0.091971078052723032;
        QUADNODES(1,417-412) = 0.31979362782962989;
        QUADNODES(1,418-412) = 0.31979362782962989;
        QUADNODES(1,419-412) = 0.040619116511110234;
        QUADNODES(1,420-412) = 0.31979362782962989;
        QUADNODES(1,421-412) = 0.056350832689629149;
        QUADNODES(1,422-412) = 0.4436491673103708;
        QUADNODES(1,423-412) = 0.4436491673103708;
        QUADNODES(1,424-412) = 0.056350832689629149;
        QUADNODES(1,425-412) = 0.056350832689629149;
        QUADNODES(1,426-412) = 0.4436491673103708;

        QUADNODES(2,412-412) = 0.2500000;
        QUADNODES(2,413-412) = 0.091971078052723032;
        QUADNODES(2,414-412) = 0.091971078052723032;
        QUADNODES(2,415-412) = 0.091971078052723032;
        QUADNODES(2,416-412) = 0.72408676584183096;
        QUADNODES(2,417-412) = 0.31979362782962989;
        QUADNODES(2,418-412) = 0.31979362782962989;
        QUADNODES(2,419-412) = 0.31979362782962989;
        QUADNODES(2,420-412) = 0.040619116511110234;
        QUADNODES(2,421-412) = 0.4436491673103708;
        QUADNODES(2,422-412) = 0.056350832689629149;
        QUADNODES(2,423-412) = 0.4436491673103708;
        QUADNODES(2,424-412) = 0.056350832689629149;
        QUADNODES(2,425-412) = 0.4436491673103708;
        QUADNODES(2,426-412) = 0.056350832689629149;

    } else if (numQuadNodes==24) {
        QUADWEIGHTS(427-427) = 0.006653792;
        QUADWEIGHTS(428-427) = 0.006653792;
        QUADWEIGHTS(429-427) = 0.006653792;
        QUADWEIGHTS(430-427) = 0.006653792;
        QUADWEIGHTS(431-427) = 0.001679535;
        QUADWEIGHTS(432-427) = 0.001679535;
        QUADWEIGHTS(433-427) = 0.001679535;
        QUADWEIGHTS(434-427) = 0.001679535;
        QUADWEIGHTS(435-427) = 0.009226197;
        QUADWEIGHTS(436-427) = 0.009226197;
        QUADWEIGHTS(437-427) = 0.009226197;
        QUADWEIGHTS(438-427) = 0.009226197;
        QUADWEIGHTS(439-427) = 0.008035714;
        QUADWEIGHTS(440-427) = 0.008035714;
        QUADWEIGHTS(441-427) = 0.008035714;
        QUADWEIGHTS(442-427) = 0.008035714;
        QUADWEIGHTS(443-427) = 0.008035714;
        QUADWEIGHTS(444-427) = 0.008035714;
        QUADWEIGHTS(445-427) = 0.008035714;
        QUADWEIGHTS(446-427) = 0.008035714;
        QUADWEIGHTS(447-427) = 0.008035714;
        QUADWEIGHTS(448-427) = 0.008035714;
        QUADWEIGHTS(449-427) = 0.008035714;
        QUADWEIGHTS(450-427) = 0.008035714;

        QUADNODES(0,427-427) = 0.3561914;
        QUADNODES(0,428-427) = 0.2146029;
        QUADNODES(0,429-427) = 0.2146029;
        QUADNODES(0,430-427) = 0.2146029;
        QUADNODES(0,431-427) = 0.8779781;
        QUADNODES(0,432-427) = 0.04067396;
        QUADNODES(0,433-427) = 0.04067396;
        QUADNODES(0,434-427) = 0.04067396;
        QUADNODES(0,435-427) = 0.03298633;
        QUADNODES(0,436-427) = 0.3223379;
        QUADNODES(0,437-427) = 0.3223379;
        QUADNODES(0,438-427) = 0.3223379;
        QUADNODES(0,439-427) = 0.6030057;
        QUADNODES(0,440-427) = 0.6030057;
        QUADNODES(0,441-427) = 0.6030057;
        QUADNODES(0,442-427) = 0.2696723;
        QUADNODES(0,443-427) = 0.2696723;
        QUADNODES(0,444-427) = 0.2696723;
        QUADNODES(0,445-427) = 0.06366100;
        QUADNODES(0,446-427) = 0.06366100;
        QUADNODES(0,447-427) = 0.06366100;
        QUADNODES(0,448-427) = 0.06366100;
        QUADNODES(0,449-427) = 0.06366100;
        QUADNODES(0,450-427) = 0.06366100;

        QUADNODES(1,427-427) = 0.2146029;
        QUADNODES(1,428-427) = 0.3561914;
        QUADNODES(1,429-427) = 0.2146029;
        QUADNODES(1,430-427) = 0.2146029;
        QUADNODES(1,431-427) = 0.04067396;
        QUADNODES(1,432-427) = 0.8779781;
        QUADNODES(1,433-427) = 0.04067396;
        QUADNODES(1,434-427) = 0.04067396;
        QUADNODES(1,435-427) = 0.3223379;
        QUADNODES(1,436-427) = 0.03298633;
        QUADNODES(1,437-427) = 0.3223379;
        QUADNODES(1,438-427) = 0.3223379;
        QUADNODES(1,439-427) = 0.2696723;
        QUADNODES(1,440-427) = 0.06366100;
        QUADNODES(1,441-427) = 0.06366100;
        QUADNODES(1,442-427) = 0.6030057;
        QUADNODES(1,443-427) = 0.06366100;
        QUADNODES(1,444-427) = 0.06366100;
        QUADNODES(1,445-427) = 0.6030057;
        QUADNODES(1,446-427) = 0.6030057;
        QUADNODES(1,447-427) = 0.2696723;
        QUADNODES(1,448-427) = 0.2696723;
        QUADNODES(1,449-427) = 0.06366100;
        QUADNODES(1,450-427) = 0.06366100;

        QUADNODES(2,427-427) = 0.2146029;
        QUADNODES(2,428-427) = 0.2146029;
        QUADNODES(2,429-427) = 0.3561914;
        QUADNODES(2,430-427) = 0.2146029;
        QUADNODES(2,431-427) = 0.04067396;
        QUADNODES(2,432-427) = 0.04067396;
        QUADNODES(2,433-427) = 0.8779781;
        QUADNODES(2,434-427) = 0.04067396;
        QUADNODES(2,435-427) = 0.3223379;
        QUADNODES(2,436-427) = 0.3223379;
        QUADNODES(2,437-427) = 0.03298633;
        QUADNODES(2,438-427) = 0.3223379;
        QUADNODES(2,439-427) = 0.06366100;
        QUADNODES(2,440-427) = 0.2696723;
        QUADNODES(2,441-427) = 0.06366100;
        QUADNODES(2,442-427) = 0.06366100;
        QUADNODES(2,443-427) = 0.6030057;
        QUADNODES(2,444-427) = 0.06366100;
        QUADNODES(2,445-427) = 0.2696723;
        QUADNODES(2,446-427) = 0.06366100;
        QUADNODES(2,447-427) = 0.6030057;
        QUADNODES(2,448-427) = 0.06366100;
        QUADNODES(2,449-427) = 0.6030057;
        QUADNODES(2,450-427) = 0.2696723;

    } else if (numQuadNodes==31) {
        QUADWEIGHTS(451-451) = 0.01826422;
        QUADWEIGHTS(452-451) = 0.01059994;
        QUADWEIGHTS(453-451) = 0.01059994;
        QUADWEIGHTS(454-451) = 0.01059994;
        QUADWEIGHTS(455-451) = 0.01059994;
        QUADWEIGHTS(456-451) =-0.06251774;
        QUADWEIGHTS(457-451) =-0.06251774;
        QUADWEIGHTS(458-451) =-0.06251774;
        QUADWEIGHTS(459-451) =-0.06251774;
        QUADWEIGHTS(460-451) = 0.004891425;
        QUADWEIGHTS(461-451) = 0.004891425;
        QUADWEIGHTS(462-451) = 0.004891425;
        QUADWEIGHTS(463-451) = 0.004891425;
        QUADWEIGHTS(464-451) = 0.0009700176;
        QUADWEIGHTS(465-451) = 0.0009700176;
        QUADWEIGHTS(466-451) = 0.0009700176;
        QUADWEIGHTS(467-451) = 0.0009700176;
        QUADWEIGHTS(468-451) = 0.0009700176;
        QUADWEIGHTS(469-451) = 0.0009700176;
        QUADWEIGHTS(470-451) = 0.02755732;
        QUADWEIGHTS(471-451) = 0.02755732;
        QUADWEIGHTS(472-451) = 0.02755732;
        QUADWEIGHTS(473-451) = 0.02755732;
        QUADWEIGHTS(474-451) = 0.02755732;
        QUADWEIGHTS(475-451) = 0.02755732;
        QUADWEIGHTS(476-451) = 0.02755732;
        QUADWEIGHTS(477-451) = 0.02755732;
        QUADWEIGHTS(478-451) = 0.02755732;
        QUADWEIGHTS(479-451) = 0.02755732;
        QUADWEIGHTS(480-451) = 0.02755732;
        QUADWEIGHTS(481-451) = 0.02755732;

        QUADNODES(0,451-451) = 0.2500000;
        QUADNODES(0,452-451) = 0.7653604;
        QUADNODES(0,453-451) = 0.07821319;
        QUADNODES(0,454-451) = 0.07821319;
        QUADNODES(0,455-451) = 0.07821319;
        QUADNODES(0,456-451) = 0.6344704;
        QUADNODES(0,457-451) = 0.1218432;
        QUADNODES(0,458-451) = 0.1218432;
        QUADNODES(0,459-451) = 0.1218432;
        QUADNODES(0,460-451) = 0.002382507;
        QUADNODES(0,461-451) = 0.3325392;
        QUADNODES(0,462-451) = 0.3325392;
        QUADNODES(0,463-451) = 0.3325392;
        QUADNODES(0,464-451) = 0.0000000;
        QUADNODES(0,465-451) = 0.0000000;
        QUADNODES(0,466-451) = 0.0000000;
        QUADNODES(0,467-451) = 0.5000000;
        QUADNODES(0,468-451) = 0.5000000;
        QUADNODES(0,469-451) = 0.5000000;
        QUADNODES(0,470-451) = 0.6000000;
        QUADNODES(0,471-451) = 0.6000000;
        QUADNODES(0,472-451) = 0.6000000;
        QUADNODES(0,473-451) = 0.2000000;
        QUADNODES(0,474-451) = 0.2000000;
        QUADNODES(0,475-451) = 0.2000000;
        QUADNODES(0,476-451) = 0.1000000;
        QUADNODES(0,477-451) = 0.1000000;
        QUADNODES(0,478-451) = 0.1000000;
        QUADNODES(0,479-451) = 0.1000000;
        QUADNODES(0,480-451) = 0.1000000;
        QUADNODES(0,481-451) = 0.1000000;

        QUADNODES(1,451-451) = 0.2500000;
        QUADNODES(1,452-451) = 0.07821319;
        QUADNODES(1,453-451) = 0.7653604;
        QUADNODES(1,454-451) = 0.07821319;
        QUADNODES(1,455-451) = 0.07821319;
        QUADNODES(1,456-451) = 0.1218432;
        QUADNODES(1,457-451) = 0.6344704;
        QUADNODES(1,458-451) = 0.1218432;
        QUADNODES(1,459-451) = 0.1218432;
        QUADNODES(1,460-451) = 0.3325392;
        QUADNODES(1,461-451) = 0.002382507;
        QUADNODES(1,462-451) = 0.3325392;
        QUADNODES(1,463-451) = 0.3325392;
        QUADNODES(1,464-451) = 0.0000000;
        QUADNODES(1,465-451) = 0.5000000;
        QUADNODES(1,466-451) = 0.5000000;
        QUADNODES(1,467-451) = 0.0000000;
        QUADNODES(1,468-451) = 0.0000000;
        QUADNODES(1,469-451) = 0.5000000;
        QUADNODES(1,470-451) = 0.2000000;
        QUADNODES(1,471-451) = 0.1000000;
        QUADNODES(1,472-451) = 0.1000000;
        QUADNODES(1,473-451) = 0.6000000;
        QUADNODES(1,474-451) = 0.1000000;
        QUADNODES(1,475-451) = 0.1000000;
        QUADNODES(1,476-451) = 0.6000000;
        QUADNODES(1,477-451) = 0.6000000;
        QUADNODES(1,478-451) = 0.2000000;
        QUADNODES(1,479-451) = 0.2000000;
        QUADNODES(1,480-451) = 0.1000000;
        QUADNODES(1,481-451) = 0.1000000;

        QUADNODES(2,451-451) = 0.2500000;
        QUADNODES(2,452-451) = 0.07821319;
        QUADNODES(2,453-451) = 0.07821319;
        QUADNODES(2,454-451) = 0.7653604;
        QUADNODES(2,455-451) = 0.07821319;
        QUADNODES(2,456-451) = 0.1218432;
        QUADNODES(2,457-451) = 0.1218432;
        QUADNODES(2,458-451) = 0.6344704;
        QUADNODES(2,459-451) = 0.1218432;
        QUADNODES(2,460-451) = 0.3325392;
        QUADNODES(2,461-451) = 0.3325392;
        QUADNODES(2,462-451) = 0.002382507;
        QUADNODES(2,463-451) = 0.3325392;
        QUADNODES(2,464-451) = 0.5000000;
        QUADNODES(2,465-451) = 0.0000000;
        QUADNODES(2,466-451) = 0.5000000;
        QUADNODES(2,467-451) = 0.0000000;
        QUADNODES(2,468-451) = 0.5000000;
        QUADNODES(2,469-451) = 0.0000000;
        QUADNODES(2,470-451) = 0.1000000;
        QUADNODES(2,471-451) = 0.2000000;
        QUADNODES(2,472-451) = 0.1000000;
        QUADNODES(2,473-451) = 0.1000000;
        QUADNODES(2,474-451) = 0.6000000;
        QUADNODES(2,475-451) = 0.1000000;
        QUADNODES(2,476-451) = 0.2000000;
        QUADNODES(2,477-451) = 0.1000000;
        QUADNODES(2,478-451) = 0.6000000;
        QUADNODES(2,479-451) = 0.1000000;
        QUADNODES(2,480-451) = 0.6000000;
        QUADNODES(2,481-451) = 0.2000000;

    } else if (numQuadNodes==45) {
        QUADWEIGHTS(482-482) =-0.03932701;
        QUADWEIGHTS(483-482) = 0.004081316;
        QUADWEIGHTS(484-482) = 0.004081316;
        QUADWEIGHTS(485-482) = 0.004081316;
        QUADWEIGHTS(486-482) = 0.004081316;
        QUADWEIGHTS(487-482) = 0.0006580868;
        QUADWEIGHTS(488-482) = 0.0006580868;
        QUADWEIGHTS(489-482) = 0.0006580868;
        QUADWEIGHTS(490-482) = 0.0006580868;
        QUADWEIGHTS(491-482) = 0.004384259;
        QUADWEIGHTS(492-482) = 0.004384259;
        QUADWEIGHTS(493-482) = 0.004384259;
        QUADWEIGHTS(494-482) = 0.004384259;
        QUADWEIGHTS(495-482) = 0.004384259;
        QUADWEIGHTS(496-482) = 0.004384259;
        QUADWEIGHTS(497-482) = 0.01383006;
        QUADWEIGHTS(498-482) = 0.01383006;
        QUADWEIGHTS(499-482) = 0.01383006;
        QUADWEIGHTS(500-482) = 0.01383006;
        QUADWEIGHTS(501-482) = 0.01383006;
        QUADWEIGHTS(502-482) = 0.01383006;
        QUADWEIGHTS(503-482) = 0.004240437;
        QUADWEIGHTS(504-482) = 0.004240437;
        QUADWEIGHTS(505-482) = 0.004240437;
        QUADWEIGHTS(506-482) = 0.004240437;
        QUADWEIGHTS(507-482) = 0.004240437;
        QUADWEIGHTS(508-482) = 0.004240437;
        QUADWEIGHTS(509-482) = 0.004240437;
        QUADWEIGHTS(510-482) = 0.004240437;
        QUADWEIGHTS(511-482) = 0.004240437;
        QUADWEIGHTS(512-482) = 0.004240437;
        QUADWEIGHTS(513-482) = 0.004240437;
        QUADWEIGHTS(514-482) = 0.004240437;
        QUADWEIGHTS(515-482) = 0.002238740;
        QUADWEIGHTS(516-482) = 0.002238740;
        QUADWEIGHTS(517-482) = 0.002238740;
        QUADWEIGHTS(518-482) = 0.002238740;
        QUADWEIGHTS(519-482) = 0.002238740;
        QUADWEIGHTS(520-482) = 0.002238740;
        QUADWEIGHTS(521-482) = 0.002238740;
        QUADWEIGHTS(522-482) = 0.002238740;
        QUADWEIGHTS(523-482) = 0.002238740;
        QUADWEIGHTS(524-482) = 0.002238740;
        QUADWEIGHTS(525-482) = 0.002238740;
        QUADWEIGHTS(526-482) = 0.002238740;

        QUADNODES(0,482-482) = 0.2500000;
        QUADNODES(0,483-482) = 0.6175872;
        QUADNODES(0,484-482) = 0.1274709;
        QUADNODES(0,485-482) = 0.1274709;
        QUADNODES(0,486-482) = 0.1274709;
        QUADNODES(0,487-482) = 0.9037635;
        QUADNODES(0,488-482) = 0.03207883;
        QUADNODES(0,489-482) = 0.03207883;
        QUADNODES(0,490-482) = 0.03207883;
        QUADNODES(0,491-482) = 0.4502229;
        QUADNODES(0,492-482) = 0.4502229;
        QUADNODES(0,493-482) = 0.4502229;
        QUADNODES(0,494-482) = 0.04977710;
        QUADNODES(0,495-482) = 0.04977710;
        QUADNODES(0,496-482) = 0.04977710;
        QUADNODES(0,497-482) = 0.3162696;
        QUADNODES(0,498-482) = 0.3162696;
        QUADNODES(0,499-482) = 0.3162696;
        QUADNODES(0,500-482) = 0.1837304;
        QUADNODES(0,501-482) = 0.1837304;
        QUADNODES(0,502-482) = 0.1837304;
        QUADNODES(0,503-482) = 0.5132800;
        QUADNODES(0,504-482) = 0.5132800;
        QUADNODES(0,505-482) = 0.5132800;
        QUADNODES(0,506-482) = 0.02291779;
        QUADNODES(0,507-482) = 0.02291779;
        QUADNODES(0,508-482) = 0.02291779;
        QUADNODES(0,509-482) = 0.2319011;
        QUADNODES(0,510-482) = 0.2319011;
        QUADNODES(0,511-482) = 0.2319011;
        QUADNODES(0,512-482) = 0.2319011;
        QUADNODES(0,513-482) = 0.2319011;
        QUADNODES(0,514-482) = 0.2319011;
        QUADNODES(0,515-482) = 0.1937465;
        QUADNODES(0,516-482) = 0.1937465;
        QUADNODES(0,517-482) = 0.1937465;
        QUADNODES(0,518-482) = 0.7303134;
        QUADNODES(0,519-482) = 0.7303134;
        QUADNODES(0,520-482) = 0.7303134;
        QUADNODES(0,521-482) = 0.03797005;
        QUADNODES(0,522-482) = 0.03797005;
        QUADNODES(0,523-482) = 0.03797005;
        QUADNODES(0,524-482) = 0.03797005;
        QUADNODES(0,525-482) = 0.03797005;
        QUADNODES(0,526-482) = 0.03797005;

        QUADNODES(1,482-482) = 0.2500000;
        QUADNODES(1,483-482) = 0.1274709;
        QUADNODES(1,484-482) = 0.6175872;
        QUADNODES(1,485-482) = 0.1274709;
        QUADNODES(1,486-482) = 0.1274709;
        QUADNODES(1,487-482) = 0.03207883;
        QUADNODES(1,488-482) = 0.9037635;
        QUADNODES(1,489-482) = 0.03207883;
        QUADNODES(1,490-482) = 0.03207883;
        QUADNODES(1,491-482) = 0.4502229;
        QUADNODES(1,492-482) = 0.04977710;
        QUADNODES(1,493-482) = 0.04977710;
        QUADNODES(1,494-482) = 0.4502229;
        QUADNODES(1,495-482) = 0.4502229;
        QUADNODES(1,496-482) = 0.04977710;
        QUADNODES(1,497-482) = 0.3162696;
        QUADNODES(1,498-482) = 0.1837304;
        QUADNODES(1,499-482) = 0.1837304;
        QUADNODES(1,500-482) = 0.3162696;
        QUADNODES(1,501-482) = 0.3162696;
        QUADNODES(1,502-482) = 0.1837304;
        QUADNODES(1,503-482) = 0.02291779;
        QUADNODES(1,504-482) = 0.2319011;
        QUADNODES(1,505-482) = 0.2319011;
        QUADNODES(1,506-482) = 0.5132800;
        QUADNODES(1,507-482) = 0.2319011;
        QUADNODES(1,508-482) = 0.2319011;
        QUADNODES(1,509-482) = 0.5132800;
        QUADNODES(1,510-482) = 0.5132800;
        QUADNODES(1,511-482) = 0.02291779;
        QUADNODES(1,512-482) = 0.02291779;
        QUADNODES(1,513-482) = 0.2319011;
        QUADNODES(1,514-482) = 0.2319011;
        QUADNODES(1,515-482) = 0.7303134;
        QUADNODES(1,516-482) = 0.03797005;
        QUADNODES(1,517-482) = 0.03797005;
        QUADNODES(1,518-482) = 0.1937465;
        QUADNODES(1,519-482) = 0.03797005;
        QUADNODES(1,520-482) = 0.03797005;
        QUADNODES(1,521-482) = 0.1937465;
        QUADNODES(1,522-482) = 0.1937465;
        QUADNODES(1,523-482) = 0.7303134;
        QUADNODES(1,524-482) = 0.7303134;
        QUADNODES(1,525-482) = 0.03797005;
        QUADNODES(1,526-482) = 0.03797005;

        QUADNODES(2,482-482) = 0.2500000;
        QUADNODES(2,483-482) = 0.1274709;
        QUADNODES(2,484-482) = 0.1274709;
        QUADNODES(2,485-482) = 0.6175872;
        QUADNODES(2,486-482) = 0.1274709;
        QUADNODES(2,487-482) = 0.03207883;
        QUADNODES(2,488-482) = 0.03207883;
        QUADNODES(2,489-482) = 0.9037635;
        QUADNODES(2,490-482) = 0.03207883;
        QUADNODES(2,491-482) = 0.04977710;
        QUADNODES(2,492-482) = 0.4502229;
        QUADNODES(2,493-482) = 0.04977710;
        QUADNODES(2,494-482) = 0.4502229;
        QUADNODES(2,495-482) = 0.04977710;
        QUADNODES(2,496-482) = 0.4502229;
        QUADNODES(2,497-482) = 0.1837304;
        QUADNODES(2,498-482) = 0.3162696;
        QUADNODES(2,499-482) = 0.1837304;
        QUADNODES(2,500-482) = 0.3162696;
        QUADNODES(2,501-482) = 0.1837304;
        QUADNODES(2,502-482) = 0.3162696;
        QUADNODES(2,503-482) = 0.2319011;
        QUADNODES(2,504-482) = 0.02291779;
        QUADNODES(2,505-482) = 0.2319011;
        QUADNODES(2,506-482) = 0.2319011;
        QUADNODES(2,507-482) = 0.5132800;
        QUADNODES(2,508-482) = 0.2319011;
        QUADNODES(2,509-482) = 0.02291779;
        QUADNODES(2,510-482) = 0.2319011;
        QUADNODES(2,511-482) = 0.5132800;
        QUADNODES(2,512-482) = 0.2319011;
        QUADNODES(2,513-482) = 0.5132800;
        QUADNODES(2,514-482) = 0.02291779;
        QUADNODES(2,515-482) = 0.03797005;
        QUADNODES(2,516-482) = 0.7303134;
        QUADNODES(2,517-482) = 0.03797005;
        QUADNODES(2,518-482) = 0.03797005;
        QUADNODES(2,519-482) = 0.1937465;
        QUADNODES(2,520-482) = 0.03797005;
        QUADNODES(2,521-482) = 0.7303134;
        QUADNODES(2,522-482) = 0.03797005;
        QUADNODES(2,523-482) = 0.1937465;
        QUADNODES(2,524-482) = 0.03797005;
        QUADNODES(2,525-482) = 0.1937465;
        QUADNODES(2,526-482) = 0.7303134;

    } else {
        // get scheme on [0.1]^3
        Quad_getNodesHex(numQuadNodes, quadNodes, quadWeights);

        // squeeze it:
        for (int i=0; i<numQuadNodes; i++) {
            const double Q1=QUADNODES(0,i);
            const double Q2=QUADNODES(1,i);
            const double Q3=QUADNODES(2,i);
            const double JA11 = (1./3.)*Q2*Q3-(1./2.)*(Q2+Q3)+1.;
            const double JA12 = (1./3.)*Q1*Q3-(1./2.)*Q1;
            const double JA13 = (1./3.)*Q1*Q2-(1./2.)*Q1;
            const double JA21 = (1./3.)*Q2*Q3-(1./2.)*Q2;
            const double JA22 = (1./3.)*Q1*Q3-(1./2.)*(Q1+Q3)+1.;
            const double JA23 = (1./3.)*Q1*Q2-(1./2.)*Q2;
            const double JA31 = (1./3.)*Q2*Q3-(1./2.)*Q3;
            const double JA32 = (1./3.)*Q1*Q3-(1./2.)*Q3;
            const double JA33 = (1./3.)*Q1*Q2-(1./2.)*(Q1+Q2)+1.;
            const double DET = JA11*JA22*JA33 + JA12*JA23*JA31 + JA13*JA21*JA32
                              -JA13*JA22*JA31 - JA11*JA23*JA32 - JA12*JA21*JA33;
            quadWeights[i]=quadWeights[i]*std::abs(DET);
            QUADNODES(0,i)=Q1*((1./3.)*Q2*Q3-(1./2.)*(Q2+Q3)+1.);
            QUADNODES(1,i)=Q2*((1./3.)*Q1*Q3-(1./2.)*(Q1+Q3)+1.);
            QUADNODES(2,i)=Q3*((1./3.)*Q1*Q2-(1./2.)*(Q1+Q2)+1.);
        }
    } // numQuadNodes
#undef DIM
}


/// get a quadrature scheme with numQuadNodes quadrature nodes for the
/// quad [0.1]^2 as a X-product of a 1D scheme
void Quad_getNodesRec(int numQuadNodes, std::vector<double>& quadNodes, std::vector<double>& quadWeights)
{
#define DIM 2
    bool set=false;
    std::vector<double> quadNodes1d(numQuadNodes);
    std::vector<double> quadWeights1d(numQuadNodes);

    // find numQuadNodes1d with numQuadNodes1d**2==numQuadNodes:
    for (int numQuadNodes1d=1; numQuadNodes1d<=MAX_numQuadNodesLine; numQuadNodes1d++) {
        if (numQuadNodes1d*numQuadNodes1d==numQuadNodes) {
            // get 1D scheme:
            Quad_getNodesLine(numQuadNodes1d, quadNodes1d, quadWeights1d);

            // make 2D scheme:
            int l=0;
            for (int i=0; i<numQuadNodes1d; i++) {
                for (int j=0; j<numQuadNodes1d; j++) {
                    QUADNODES(0,l)=quadNodes1d[i];
                    QUADNODES(1,l)=quadNodes1d[j];
                    QUADWEIGHTS(l)=quadWeights1d[i]*quadWeights1d[j];
                    l++;
                }
            }
            set=true;
            break;
        }
    }
    if (!set) {
        std::stringstream ss;
        ss << "Quad_getNodesRec: Illegal number of quadrature nodes "
           << numQuadNodes << " on hexahedron.";
        const std::string msg(ss.str());
        throw escript::ValueError(msg);
    }
#undef DIM
}

/// get a quadrature scheme with numQuadNodes quadrature nodes for the
/// hex [0.1]^3 as a X-product of a 1D scheme
void Quad_getNodesHex(int numQuadNodes, std::vector<double>& quadNodes, std::vector<double>& quadWeights)
{
#define DIM 3
    bool set=false;
    std::vector<double> quadNodes1d(numQuadNodes);
    std::vector<double> quadWeights1d(numQuadNodes);

    // find numQuadNodes1d with numQuadNodes1d**3==numQuadNodes:
    for (int numQuadNodes1d=1; numQuadNodes1d<=MAX_numQuadNodesLine; numQuadNodes1d++) {
        if (numQuadNodes1d*numQuadNodes1d*numQuadNodes1d==numQuadNodes) {
            // get 1D scheme:
            Quad_getNodesLine(numQuadNodes1d, quadNodes1d, quadWeights1d);

            // make 3D scheme:
            int l=0;
            for (int i=0; i<numQuadNodes1d; i++) {
                for (int j=0; j<numQuadNodes1d; j++) {
                    for (int k=0; k<numQuadNodes1d; k++) {
                        QUADNODES(0,l)=quadNodes1d[i];
                        QUADNODES(1,l)=quadNodes1d[j];
                        QUADNODES(2,l)=quadNodes1d[k];
                        QUADWEIGHTS(l)=quadWeights1d[i]*quadWeights1d[j]*quadWeights1d[k];
                        l++;
                    }
                }
            }
            set=true;
            break;
        }
    }
    if (!set) {
        std::stringstream ss;
        ss << "Quad_getNodesHex: Illegal number of quadrature nodes "
           << numQuadNodes << " on hexahedron.";
        const std::string msg(ss.str());
        throw escript::ValueError(msg);
    }
#undef DIM
}

/// get a quadrature scheme with numQuadNodes quadrature nodes for a point.
/// As there is no quadrature scheme for a point this function returns without
/// changing the arrays but throws an error if numQuadNodes is negative
void Quad_getNodesPoint(int numQuadNodes, std::vector<double>& quadNodes, std::vector<double>& quadWeights)
{
    if (numQuadNodes<0)
        throw escript::ValueError(
                "Quad_getNodesPoint: Illegal number of quadrature nodes.");
}

/// get a quadrature scheme with numQuadNodes quadrature nodes on the
/// line [0,1]. The nodes and weights are set from a table
void Quad_getNodesLine(int numQuadNodes, std::vector<double>& quadNodes, std::vector<double>& quadWeights)
{
    switch (numQuadNodes) {
        case 1:
            quadNodes[0]=0.5;
            quadWeights[0]=1.;
            break;

        case 2:
            quadNodes[0]=(1.-.577350269189626)/2.;
            quadNodes[1]=(1.+.577350269189626)/2.;
            quadWeights[0]=.5;
            quadWeights[1]=.5;
            break;

        case 3:
            quadNodes[0]=(1.-.774596669241483)/2.;
            quadNodes[1]=.5;
            quadNodes[2]=(1.+.774596669241483)/2.;
            quadWeights[0]=5./18.;
            quadWeights[1]=4./ 9.;
            quadWeights[2]=5./18.;
            break;

        case 4:
            quadNodes[0]=(1.-.861136311594053)/2.;
            quadNodes[1]=(1.-.339981043584856)/2.;
            quadNodes[2]=(1.+.339981043584856)/2.;
            quadNodes[3]=(1.+.861136311594053)/2.;
            quadWeights[0]=.347854845137454/2.;
            quadWeights[1]=.652145154862546/2.;
            quadWeights[2]=.652145154862546/2.;
            quadWeights[3]=.347854845137454/2.;
            break;

        case 5:
            quadNodes[0]=(1.-.906179845938664)/2.;
            quadNodes[1]=(1.-.538469310105683)/2.;
            quadNodes[2]= .5;
            quadNodes[3]=(1.+.538469310105683)/2.;
            quadNodes[4]=(1.+.906179845938664)/2.;
            quadWeights[0]=.236926885056189/2.;
            quadWeights[1]=.478628670499366/2.;
            quadWeights[2]=.568888888888889/2.;
            quadWeights[3]=.478628670499366/2.;
            quadWeights[4]=.236926885056189/2.;
            break;

        case 6:
            quadNodes[0]=(1.-.932469514203152)/2.;
            quadNodes[1]=(1.-.661209386466265)/2.;
            quadNodes[2]=(1.-.238619186083197)/2.;
            quadNodes[3]=(1.+.238619186083197)/2.;
            quadNodes[4]=(1.+.661209386466265)/2.;
            quadNodes[5]=(1.+.932469514203152)/2.;
            quadWeights[0]=.171324492379170/2.;
            quadWeights[1]=.360761573048139/2.;
            quadWeights[2]=.467913934572691/2.;
            quadWeights[3]=.467913934572691/2.;
            quadWeights[4]=.360761573048139/2.;
            quadWeights[5]=.171324492379170/2.;
            break;

        case 7:
            quadNodes[0]=(1.-.949107912342759)/2.;
            quadNodes[1]=(1.-.741531185599394)/2.;
            quadNodes[2]=(1.-.405845151377397)/2.;
            quadNodes[3]=0.5;
            quadNodes[4]=(1.+.405845151377397)/2.;
            quadNodes[5]=(1.+.741531185599394)/2.;
            quadNodes[6]=(1.+.949107912342759)/2.;
            quadWeights[0]= .129484966168870/2.;
            quadWeights[1]= .279705391489277/2.;
            quadWeights[2]= .381830050505119/2.;
            quadWeights[3]= .417959183673469/2.;
            quadWeights[4]= .381830050505119/2.;
            quadWeights[5]= .279705391489277/2.;
            quadWeights[6]= .129484966168870/2.;
            break;

        case 8:
            quadNodes[0]=(1.-.960289856497536)/2.;
            quadNodes[1]=(1.-.796666477413627)/2.;
            quadNodes[2]=(1.-.525532409916329)/2.;
            quadNodes[3]=(1.-.183434642495650)/2.;
            quadNodes[4]=(1.+.183434642495650)/2.;
            quadNodes[5]=(1.+.525532409916329)/2.;
            quadNodes[6]=(1.+.796666477413627)/2.;
            quadNodes[7]=(1.+.960289856497536)/2.;
            quadWeights[0]= .101228536290376/2.;
            quadWeights[1]= .222381034453374/2.;
            quadWeights[2]= .313706645877887/2.;
            quadWeights[3]= .362683783378362/2.;
            quadWeights[4]= .362683783378362/2.;
            quadWeights[5]= .313706645877887/2.;
            quadWeights[6]= .222381034453374/2.;
            quadWeights[7]= .101228536290376/2.;
            break;

        case 9:
            quadNodes[0]=(1.-.968160239507626)/2.;
            quadNodes[1]=(1.-.836031107326636)/2.;
            quadNodes[2]=(1.-.613371432700590)/2.;
            quadNodes[3]=(1.-.324253423403809)/2.;
            quadNodes[4]= .5;
            quadNodes[5]=(1.+.324253423403809)/2.;
            quadNodes[6]=(1.+.613371432700590)/2.;
            quadNodes[7]=(1.+.836031107326636)/2.;
            quadNodes[8]=(1.+.968160239507626)/2.;
            quadWeights[0]= .081274388361574/2.;
            quadWeights[1]= .180648160694857/2.;
            quadWeights[2]= .260610696402935/2.;
            quadWeights[3]= .312347077040003/2.;
            quadWeights[4]= .330239355001260/2.;
            quadWeights[5]= .312347077040003/2.;
            quadWeights[6]= .260610696402935/2.;
            quadWeights[7]= .180648160694857/2.;
            quadWeights[8]= .081274388361574/2.;
            break;

        case 10:
            quadNodes[0]=(1.-.973906528517172)/2.;
            quadNodes[1]=(1.-.865063366688985)/2.;
            quadNodes[2]=(1.-.679409568299024)/2.;
            quadNodes[3]=(1.-.433395394129247)/2.;
            quadNodes[4]=(1.-.148874338981631)/2.;
            quadNodes[5]=(1.+.148874338981631)/2.;
            quadNodes[6]=(1.+.433395394129247)/2.;
            quadNodes[7]=(1.+.679409568299024)/2.;
            quadNodes[8]=(1.+.865063366688985)/2.;
            quadNodes[9]=(1.+.973906528517172)/2.;
            quadWeights[0]= .066671344308688/2.;
            quadWeights[1]= .149451349150581/2.;
            quadWeights[2]= .219086362515982/2.;
            quadWeights[3]= .269266719309996/2.;
            quadWeights[4]= .295524224714753/2.;
            quadWeights[5]= .295524224714753/2.;
            quadWeights[6]= .269266719309996/2.;
            quadWeights[7]= .219086362515982/2.;
            quadWeights[8]= .149451349150581/2.;
            quadWeights[9]= .066671344308688/2.;
            break;

        default:
            throw escript::ValueError("Quad_getNodesLine: Invalid integration order.");
    }
}


// The following functions Quad_getNumNodes* return the number of quadrature
// points needed to achieve a certain accuracy. Notice that for Tet and Tri
// the order is increased to consider the accuracy reduction through the
// construction process.

int Quad_getNumNodesPoint(int order)
{
    return 1;
}

int Quad_getNumNodesLine(int order)
{
    if (order < 0) {
        throw escript::ValueError("Quad_getNumNodesLine: Negative integration order.");
    } else if (order > 2*MAX_numQuadNodesLine-1) {
        std::stringstream ss;
        ss << "Quad_getNumNodesLine: requested integration order "
           << order << " on line is too large (>"
           << 2*MAX_numQuadNodesLine-1 << ").";
        const std::string msg(ss.str());
        throw escript::ValueError(msg);
    } else {
        return order/2+1;
    }
}

int Quad_getNumNodesTri(int order)
{
    if (order<=1) {
        return 1;
    } else if (order<=2) {
        return 3;
    } else if (order<=3) {
        return 4;
    } else if (order<=4) {
        return 6;
    } else if (order<=5) {
        return 7;
    } else if (order<=6) {
        return 12;
    } else if (order<=7) {
        return 13;
    } else if (order<=8) {
        return 16;
    } else if (order<=9) {
        return 19;
    } else {
        const int numQuadNodesLine=Quad_getNumNodesLine(order+1);
        return numQuadNodesLine*numQuadNodesLine;
    }
}

int Quad_getNumNodesRec(int order)
{
    const int numQuadNodesLine=Quad_getNumNodesLine(order);
    return numQuadNodesLine*numQuadNodesLine;
}

int Quad_getNumNodesTet(int order)
{
    if (order<=1) {
        return 1;
    } else if (order<=2) {
        return 4;
    } else if (order<=3) {
        return 5;
    } else if (order<=4) {
        return 11;
    } else if (order<=5) {
        return 15;
    } else if (order<=6) {
        return 24;
    } else if (order<=7) {
        return 31;
    } else if (order<=8) {
        return 45;
    } else {
        const int numQuadNodesLine=Quad_getNumNodesLine(order+2);
        return numQuadNodesLine*numQuadNodesLine*numQuadNodesLine;
    }
}

int Quad_getNumNodesHex(int order)
{
    const int numQuadNodesLine=Quad_getNumNodesLine(order);
    return numQuadNodesLine*numQuadNodesLine*numQuadNodesLine;
}

int Quad_MacroPoint(int numSubElements, int numQuadNodes,
                    const double* quadNodes, const double* quadWeights,
                    int numF, const double* dFdv, int new_len,
                    double* new_quadNodes, double* new_quadWeights,
                    double* new_dFdv)
{
    return 0;
}

int Quad_MacroLine(int numSubElements, int numQuadNodes,
                   const double* quadNodes, const double* quadWeights,
                   int numF, const double* dFdv, int new_len,
                   double* new_quadNodes, double* new_quadWeights,
                   double* new_dFdv)
{
#define DIM 1
    if (new_len < numSubElements*numQuadNodes) {
        throw FinleyException("Quad_MacroLine: array for new quadrature scheme is too small");
    }
    const double f=1./((double)numSubElements);

    for (int q=0; q<numQuadNodes; ++q) {
        const double x0=quadNodes[INDEX2(0,q,DIM)];
        const double w=f*quadWeights[q];

        for (int s=0; s<numSubElements; ++s) {
            new_quadWeights[INDEX2(q,s, numQuadNodes)]=w;
            new_quadNodes[INDEX3(0,q,s, DIM,numQuadNodes)]=(x0+s)*f;
            for (int i=0;i<numF;i++)
                new_dFdv[INDEX4(i,0,q,s, numF, DIM,numQuadNodes)] = dFdv[INDEX3(i,0,numF, q,DIM)]*f;
        }
    }
    return numSubElements*numQuadNodes;
#undef DIM
}

#define HALF 0.5
#define TWO 2.
int Quad_MacroTri(int numSubElements, int numQuadNodes,
                  const double* quadNodes, const double* quadWeights,
                  int numF, const double* dFdv, int new_len,
                  double* new_quadNodes, double* new_quadWeights,
                  double* new_dFdv)
{
#define DIM 2
    if (new_len < numSubElements*numQuadNodes) {
        throw FinleyException("Quad_MacroTri: array for new quadrature scheme is too small");
    }

    if (numSubElements==1) {
        for (int q=0; q<numQuadNodes; ++q) {
            const double x0=quadNodes[INDEX2(0,q,DIM)];
            const double x1=quadNodes[INDEX2(1,q,DIM)];
            const double w=quadWeights[q];

            new_quadWeights[INDEX2(q,0,numQuadNodes)]=w;
            new_quadNodes[INDEX3(0,q,0,DIM,numQuadNodes)]=x0;
            new_quadNodes[INDEX3(1,q,0,DIM,numQuadNodes)]=x1;
            for (int i=0; i<numF; i++) {
                new_dFdv[INDEX4(i,0,q,0, numF,DIM,numQuadNodes)] = dFdv[INDEX3(i, 0, q, numF, DIM)];
                new_dFdv[INDEX4(i,1,q,0, numF,DIM,numQuadNodes)] = dFdv[INDEX3(i, 1, q, numF, DIM)];
            }
        }
    } else if (numSubElements==4) {
        const double f = 0.25;
        for (int q=0; q<numQuadNodes; ++q) {
            const double x0=quadNodes[INDEX2(0,q,DIM)];
            const double x1=quadNodes[INDEX2(1,q,DIM)];
            const double w=f*quadWeights[q];

            new_quadWeights[INDEX2(q,0,numQuadNodes)]=w;
            new_quadNodes[INDEX3(0,q,0,DIM,numQuadNodes)]=HALF*x0;
            new_quadNodes[INDEX3(1,q,0,DIM,numQuadNodes)]=HALF*x1;

            new_quadWeights[INDEX2(q,1,numQuadNodes)]=w;
            new_quadNodes[INDEX3(0,q,1,DIM,numQuadNodes)]=HALF*x0;
            new_quadNodes[INDEX3(1,q,1,DIM,numQuadNodes)]=HALF*(x1+1);

            new_quadWeights[INDEX2(q,2,numQuadNodes)]=w;
            new_quadNodes[INDEX3(0,q,2,DIM,numQuadNodes)]=HALF*(x0+1);
            new_quadNodes[INDEX3(1,q,2,DIM,numQuadNodes)]=HALF*x1;

            new_quadWeights[INDEX2(q,3,numQuadNodes)]=w;
            new_quadNodes[INDEX3(0,q,3,DIM,numQuadNodes)]=HALF*(1-x0);
            new_quadNodes[INDEX3(1,q,3,DIM,numQuadNodes)]=HALF*(1-x1);

            for (int i=0; i<numF; i++) {
                const double df0=dFdv[INDEX3(i,0,q, numF, DIM)]*TWO;
                const double df1=dFdv[INDEX3(i,1,q, numF, DIM)]*TWO;

                new_dFdv[INDEX4(i,0,q,0, numF,DIM,numQuadNodes)] = df0;
                new_dFdv[INDEX4(i,1,q,0, numF,DIM,numQuadNodes)] = df1;

                new_dFdv[INDEX4(i,0,q,1, numF,DIM,numQuadNodes)] = df0;
                new_dFdv[INDEX4(i,1,q,1, numF,DIM,numQuadNodes)] = df1;

                new_dFdv[INDEX4(i,0,q,2, numF,DIM,numQuadNodes)] = df0;
                new_dFdv[INDEX4(i,1,q,2, numF,DIM,numQuadNodes)] = df1;

                new_dFdv[INDEX4(i,0,q,3, numF,DIM,numQuadNodes)] = -df0;
                new_dFdv[INDEX4(i,1,q,3, numF,DIM,numQuadNodes)] = -df1;
            }
        }
    } else {
        throw escript::ValueError("Quad_MacroTri: unable to create quadrature scheme for macro element.");
    }
    return numSubElements*numQuadNodes;
#undef DIM
}

int Quad_MacroRec(int numSubElements, int numQuadNodes,
                  const double* quadNodes, const double* quadWeights,
                  int numF, const double* dFdv, int new_len,
                  double* new_quadNodes, double* new_quadWeights,
                  double* new_dFdv)
{
#define DIM 2
    if (new_len < numSubElements*numQuadNodes) {
        throw FinleyException("Quad_MacroRec: array for new quadrature scheme is too small");
    }

    if (numSubElements==1) {
        for (int q=0; q<numQuadNodes; ++q) {
            const double x0=quadNodes[INDEX2(0,q,DIM)];
            const double x1=quadNodes[INDEX2(1,q,DIM)];
            const double w=quadWeights[q];

            new_quadWeights[INDEX2(q,0,numQuadNodes)]=w;
            new_quadNodes[INDEX3(0,q,0,DIM,numQuadNodes)]=x0;
            new_quadNodes[INDEX3(1,q,0,DIM,numQuadNodes)]=x1;
            for (int i=0; i<numF; i++) {
                new_dFdv[INDEX4(i,0,q,0, numF, DIM,numQuadNodes)] = dFdv[INDEX3(i,0, q,numF, DIM)];
                new_dFdv[INDEX4(i,1,q,0, numF, DIM,numQuadNodes)] = dFdv[INDEX3(i,1, q,numF, DIM)];
            }
        }
    } else if (numSubElements==4) {
        const double f = 0.25;
        for (int q=0; q<numQuadNodes; ++q) {
            const double x0=quadNodes[INDEX2(0,q,DIM)];
            const double x1=quadNodes[INDEX2(1,q,DIM)];
            const double w=f*quadWeights[q];

            new_quadWeights[INDEX2(q,0,numQuadNodes)]=w;
            new_quadNodes[INDEX3(0,q,0,DIM,numQuadNodes)]=HALF*x0;
            new_quadNodes[INDEX3(1,q,0,DIM,numQuadNodes)]=HALF*x1;

            new_quadWeights[INDEX2(q,1,numQuadNodes)]=w;
            new_quadNodes[INDEX3(0,q,1,DIM,numQuadNodes)]=HALF*x0;
            new_quadNodes[INDEX3(1,q,1,DIM,numQuadNodes)]=HALF*(x1+1);

            new_quadWeights[INDEX2(q,2,numQuadNodes)]=w;
            new_quadNodes[INDEX3(0,q,2,DIM,numQuadNodes)]=HALF*(x0+1);
            new_quadNodes[INDEX3(1,q,2,DIM,numQuadNodes)]=HALF*x1;

            new_quadWeights[INDEX2(q,3,numQuadNodes)]=w;
            new_quadNodes[INDEX3(0,q,3,DIM,numQuadNodes)]=HALF*(x0+1);
            new_quadNodes[INDEX3(1,q,3,DIM,numQuadNodes)]=HALF*(x1+1);

            for (int i=0; i<numF; i++) {
                const double df0=dFdv[INDEX3(i,0,q, numF, DIM)]*TWO;
                const double df1=dFdv[INDEX3(i,1,q, numF, DIM)]*TWO;

                new_dFdv[INDEX4(i,0,q,0, numF,DIM,numQuadNodes)] = df0;
                new_dFdv[INDEX4(i,1,q,0, numF,DIM,numQuadNodes)] = df1;

                new_dFdv[INDEX4(i,0,q,1, numF,DIM,numQuadNodes)] = df0;
                new_dFdv[INDEX4(i,1,q,1, numF,DIM,numQuadNodes)] = df1;

                new_dFdv[INDEX4(i,0,q,2, numF,DIM,numQuadNodes)] = df0;
                new_dFdv[INDEX4(i,1,q,2, numF,DIM,numQuadNodes)] = df1;

                new_dFdv[INDEX4(i,0,q,3, numF,DIM,numQuadNodes)] = df0;
                new_dFdv[INDEX4(i,1,q,3, numF,DIM,numQuadNodes)] = df1;
            }
        }
    } else {
        throw escript::ValueError("Quad_MacroRec: unable to create quadrature scheme for macro element.");
    }
    return numSubElements*numQuadNodes;
#undef DIM
}

int Quad_MacroTet(int numSubElements, int numQuadNodes, const double* quadNodes,
                  const double* quadWeights, int numF, const double* dFdv,
                  int new_len, double* new_quadNodes, double* new_quadWeights,
                  double* new_dFdv)
{
#define DIM 3
    if (new_len < numSubElements*numQuadNodes) {
        throw FinleyException("Quad_MacroTet: array for new quadrature scheme is too small");
    }

    if (numSubElements==1) {
        for (int q=0; q<numQuadNodes; ++q) {
            const double x0=quadNodes[INDEX2(0,q,DIM)];
            const double x1=quadNodes[INDEX2(1,q,DIM)];
            const double x2=quadNodes[INDEX2(2,q,DIM)];
            const double w=quadWeights[q];

            new_quadWeights[INDEX2(q,0,numQuadNodes)]=w;
            new_quadNodes[INDEX3(0,q,0,DIM,numQuadNodes)]=x0;
            new_quadNodes[INDEX3(1,q,0,DIM,numQuadNodes)]=x1;
            new_quadNodes[INDEX3(2,q,0,DIM,numQuadNodes)]=x2;

            for (int i=0; i<numF; i++) {
                new_dFdv[INDEX4(i,0,q,0, numF, DIM,numQuadNodes)] = dFdv[INDEX3(i,0,q, numF, DIM)];
                new_dFdv[INDEX4(i,1,q,0, numF, DIM,numQuadNodes)] = dFdv[INDEX3(i,1,q, numF, DIM)];
                new_dFdv[INDEX4(i,2,q,0, numF, DIM,numQuadNodes)] = dFdv[INDEX3(i,2,q, numF, DIM)];
            }
        }
    } else if (numSubElements==8) {
        const double f = 0.125;
        for (int q=0; q<numQuadNodes; ++q) {
            const double x0=quadNodes[INDEX2(0,q,DIM)];
            const double x1=quadNodes[INDEX2(1,q,DIM)];
            const double x2=quadNodes[INDEX2(2,q,DIM)];
            const double w=f*quadWeights[q];

            new_quadWeights[INDEX2(q,0,numQuadNodes)]=w;
            new_quadNodes[INDEX3(0,q,0,DIM,numQuadNodes)]=HALF*x0;
            new_quadNodes[INDEX3(1,q,0,DIM,numQuadNodes)]=HALF*x1;
            new_quadNodes[INDEX3(2,q,0,DIM,numQuadNodes)]=HALF*x2;

            new_quadWeights[INDEX2(q,1,numQuadNodes)]=w;
            new_quadNodes[INDEX3(0,q,1,DIM,numQuadNodes)]=HALF*(x0+1);
            new_quadNodes[INDEX3(1,q,1,DIM,numQuadNodes)]=HALF*x1;
            new_quadNodes[INDEX3(2,q,1,DIM,numQuadNodes)]=HALF*x2;

            new_quadWeights[INDEX2(q,2,numQuadNodes)]=w;
            new_quadNodes[INDEX3(0,q,2,DIM,numQuadNodes)]=HALF*x0;
            new_quadNodes[INDEX3(1,q,2,DIM,numQuadNodes)]=HALF*(x1+1);
            new_quadNodes[INDEX3(2,q,2,DIM,numQuadNodes)]=HALF*x2;

            new_quadWeights[INDEX2(q,3,numQuadNodes)]=w;
            new_quadNodes[INDEX3(0,q,3,DIM,numQuadNodes)]=HALF*x0;
            new_quadNodes[INDEX3(1,q,3,DIM,numQuadNodes)]=HALF*x1;
            new_quadNodes[INDEX3(2,q,3,DIM,numQuadNodes)]=HALF*(x2+1);

            new_quadWeights[INDEX2(q,4,numQuadNodes)]=w;
            new_quadNodes[INDEX3(0,q,4,DIM,numQuadNodes)]=HALF*(1-x1);
            new_quadNodes[INDEX3(1,q,4,DIM,numQuadNodes)]=HALF*(x0+x1);
            new_quadNodes[INDEX3(2,q,4,DIM,numQuadNodes)]=HALF*x2;

            new_quadWeights[INDEX2(q,5,numQuadNodes)]=w;
            new_quadNodes[INDEX3(0,q,5,DIM,numQuadNodes)]=HALF*(1-x0-x2);
            new_quadNodes[INDEX3(1,q,5,DIM,numQuadNodes)]=HALF*(1-x1);
            new_quadNodes[INDEX3(2,q,5,DIM,numQuadNodes)]=HALF*(x0+x1);

            new_quadWeights[INDEX2(q,6,numQuadNodes)]=w;
            new_quadNodes[INDEX3(0,q,6,DIM,numQuadNodes)]=HALF*x2;
            new_quadNodes[INDEX3(1,q,6,DIM,numQuadNodes)]=HALF*(1-x0-x2);
            new_quadNodes[INDEX3(2,q,6,DIM,numQuadNodes)]=HALF*(1-x1);

            new_quadWeights[INDEX2(q,7,numQuadNodes)]=w;
            new_quadNodes[INDEX3(0,q,7,DIM,numQuadNodes)]=HALF*(x0+x2);
            new_quadNodes[INDEX3(1,q,7,DIM,numQuadNodes)]=HALF*x1;
            new_quadNodes[INDEX3(2,q,7,DIM,numQuadNodes)]=HALF*(1-x0-x1);

            for (int i=0; i<numF; i++) {
                const double df0=dFdv[INDEX3(i,0,q, numF, DIM)]*TWO;
                const double df1=dFdv[INDEX3(i,1,q, numF, DIM)]*TWO;
                const double df2=dFdv[INDEX3(i,2,q, numF, DIM)]*TWO;

                new_dFdv[INDEX4(i,0,q,0, numF,DIM,numQuadNodes)] = df0;
                new_dFdv[INDEX4(i,1,q,0, numF,DIM,numQuadNodes)] = df1;
                new_dFdv[INDEX4(i,2,q,0, numF,DIM,numQuadNodes)] = df2;

                new_dFdv[INDEX4(i,0,q,1, numF,DIM,numQuadNodes)] = df0;
                new_dFdv[INDEX4(i,1,q,1, numF,DIM,numQuadNodes)] = df1;
                new_dFdv[INDEX4(i,2,q,1, numF,DIM,numQuadNodes)] = df2;

                new_dFdv[INDEX4(i,0,q,2, numF,DIM,numQuadNodes)] = df0;
                new_dFdv[INDEX4(i,1,q,2, numF,DIM,numQuadNodes)] = df1;
                new_dFdv[INDEX4(i,2,q,2, numF,DIM,numQuadNodes)] = df2;

                new_dFdv[INDEX4(i,0,q,3, numF,DIM,numQuadNodes)] = df0;
                new_dFdv[INDEX4(i,1,q,3, numF,DIM,numQuadNodes)] = df1;
                new_dFdv[INDEX4(i,2,q,3, numF,DIM,numQuadNodes)] = df2;

                new_dFdv[INDEX4(i,0,q,4, numF,DIM,numQuadNodes)] = df0-df1;
                new_dFdv[INDEX4(i,1,q,4, numF,DIM,numQuadNodes)] = df0;
                new_dFdv[INDEX4(i,2,q,4, numF,DIM,numQuadNodes)] = df2;

                new_dFdv[INDEX4(i,0,q,5, numF,DIM,numQuadNodes)] = -df2;
                new_dFdv[INDEX4(i,1,q,5, numF,DIM,numQuadNodes)] = df0-df2-df1;
                new_dFdv[INDEX4(i,2,q,5, numF,DIM,numQuadNodes)] = df0-df2;

                new_dFdv[INDEX4(i,0,q,6, numF,DIM,numQuadNodes)] = -df0+df2;
                new_dFdv[INDEX4(i,1,q,6, numF,DIM,numQuadNodes)] = -df0;
                new_dFdv[INDEX4(i,2,q,6, numF,DIM,numQuadNodes)] = -df1;

                new_dFdv[INDEX4(i,0,q,7, numF,DIM,numQuadNodes)] = df2;
                new_dFdv[INDEX4(i,1,q,7, numF,DIM,numQuadNodes)] = -df0+df1+df2;
                new_dFdv[INDEX4(i,2,q,7, numF,DIM,numQuadNodes)] = -df0+df2;
            }
        }
    } else {
        throw escript::ValueError("Quad_MacroTet: unable to create quadrature scheme for macro element.");
    }
    return numSubElements*numQuadNodes;
#undef DIM
}

int Quad_MacroHex(int numSubElements, int numQuadNodes, const double* quadNodes,
                  const double* quadWeights, int numF, const double* dFdv,
                  int new_len, double* new_quadNodes, double* new_quadWeights,
                  double* new_dFdv)
{
#define DIM 3
    if (new_len < numSubElements*numQuadNodes) {
        throw FinleyException("Quad_MacroHex: array for new quadrature scheme is too small");
    }

    if (numSubElements==1) {
        for (int q=0; q<numQuadNodes; ++q) {
            const double x0=quadNodes[INDEX2(0,q,DIM)];
            const double x1=quadNodes[INDEX2(1,q,DIM)];
            const double x2=quadNodes[INDEX2(2,q,DIM)];
            const double w=quadWeights[q];

            new_quadWeights[INDEX2(q,0,numQuadNodes)]=w;
            new_quadNodes[INDEX3(0,q,0,DIM,numQuadNodes)]=x0;
            new_quadNodes[INDEX3(1,q,0,DIM,numQuadNodes)]=x1;
            new_quadNodes[INDEX3(2,q,0,DIM,numQuadNodes)]=x2;

            for (int i=0; i<numF; i++) {
                new_dFdv[INDEX4(i,0,q,0, numF, DIM,numQuadNodes)] = dFdv[INDEX3(i,0,q,numF, DIM)];
                new_dFdv[INDEX4(i,1,q,0, numF, DIM,numQuadNodes)] = dFdv[INDEX3(i,1,q,numF, DIM)];
                new_dFdv[INDEX4(i,2,q,0, numF, DIM,numQuadNodes)] = dFdv[INDEX3(i,2,q,numF, DIM)];
            }
        }
    } else if (numSubElements==8) {
        const double f = 0.125;
        for (int q=0; q<numQuadNodes; ++q) {
            const double x0=quadNodes[INDEX2(0,q,DIM)];
            const double x1=quadNodes[INDEX2(1,q,DIM)];
            const double x2=quadNodes[INDEX2(2,q,DIM)];
            const double w=f*quadWeights[q];

            new_quadWeights[INDEX2(q,0,numQuadNodes)]=w;
            new_quadNodes[INDEX3(0,q,0,DIM,numQuadNodes)]=HALF*x0;
            new_quadNodes[INDEX3(1,q,0,DIM,numQuadNodes)]=HALF*x1;
            new_quadNodes[INDEX3(2,q,0,DIM,numQuadNodes)]=HALF*x2;

            new_quadWeights[INDEX2(q,1,numQuadNodes)]=w;
            new_quadNodes[INDEX3(0,q,1,DIM,numQuadNodes)]=HALF*(x0+1);
            new_quadNodes[INDEX3(1,q,1,DIM,numQuadNodes)]=HALF*x1;
            new_quadNodes[INDEX3(2,q,1,DIM,numQuadNodes)]=HALF*x2;

            new_quadWeights[INDEX2(q,2,numQuadNodes)]=w;
            new_quadNodes[INDEX3(0,q,2,DIM,numQuadNodes)]=HALF*x0;
            new_quadNodes[INDEX3(1,q,2,DIM,numQuadNodes)]=HALF*(x1+1);
            new_quadNodes[INDEX3(2,q,2,DIM,numQuadNodes)]=HALF*x2;

            new_quadWeights[INDEX2(q,3,numQuadNodes)]=w;
            new_quadNodes[INDEX3(0,q,3,DIM,numQuadNodes)]=HALF*(x0+1);
            new_quadNodes[INDEX3(1,q,3,DIM,numQuadNodes)]=HALF*(x1+1);
            new_quadNodes[INDEX3(2,q,3,DIM,numQuadNodes)]=HALF*x2;

            new_quadWeights[INDEX2(q,4,numQuadNodes)]=w;
            new_quadNodes[INDEX3(0,q,4,DIM,numQuadNodes)]=HALF*x0;
            new_quadNodes[INDEX3(1,q,4,DIM,numQuadNodes)]=HALF*x1;
            new_quadNodes[INDEX3(2,q,4,DIM,numQuadNodes)]=HALF*(x2+1);

            new_quadWeights[INDEX2(q,5,numQuadNodes)]=w;
            new_quadNodes[INDEX3(0,q,5,DIM,numQuadNodes)]=HALF*(x0+1);
            new_quadNodes[INDEX3(1,q,5,DIM,numQuadNodes)]=HALF*x1;
            new_quadNodes[INDEX3(2,q,5,DIM,numQuadNodes)]=HALF*(x2+1);

            new_quadWeights[INDEX2(q,6,numQuadNodes)]=w;
            new_quadNodes[INDEX3(0,q,6,DIM,numQuadNodes)]=HALF*x0;
            new_quadNodes[INDEX3(1,q,6,DIM,numQuadNodes)]=HALF*(x1+1);
            new_quadNodes[INDEX3(2,q,6,DIM,numQuadNodes)]=HALF*(x2+1);

            new_quadWeights[INDEX2(q,7,numQuadNodes)]=w;
            new_quadNodes[INDEX3(0,q,7,DIM,numQuadNodes)]=HALF*(x0+1);
            new_quadNodes[INDEX3(1,q,7,DIM,numQuadNodes)]=HALF*(x1+1);
            new_quadNodes[INDEX3(2,q,7,DIM,numQuadNodes)]=HALF*(x2+1);

            for (int i=0; i<numF; i++) {
                const double df0=dFdv[INDEX3(i,0,q, numF, DIM)]*TWO;
                const double df1=dFdv[INDEX3(i,1,q, numF, DIM)]*TWO;
                const double df2=dFdv[INDEX3(i,2,q, numF, DIM)]*TWO;

                new_dFdv[INDEX4(i,0,q,0, numF,DIM,numQuadNodes)] = df0;
                new_dFdv[INDEX4(i,1,q,0, numF,DIM,numQuadNodes)] = df1;
                new_dFdv[INDEX4(i,2,q,0, numF,DIM,numQuadNodes)] = df2;

                new_dFdv[INDEX4(i,0,q,1, numF,DIM,numQuadNodes)] = df0;
                new_dFdv[INDEX4(i,1,q,1, numF,DIM,numQuadNodes)] = df1;
                new_dFdv[INDEX4(i,2,q,1, numF,DIM,numQuadNodes)] = df2;

                new_dFdv[INDEX4(i,0,q,2, numF,DIM,numQuadNodes)] = df0;
                new_dFdv[INDEX4(i,1,q,2, numF,DIM,numQuadNodes)] = df1;
                new_dFdv[INDEX4(i,2,q,2, numF,DIM,numQuadNodes)] = df2;

                new_dFdv[INDEX4(i,0,q,3, numF,DIM,numQuadNodes)] = df0;
                new_dFdv[INDEX4(i,1,q,3, numF,DIM,numQuadNodes)] = df1;
                new_dFdv[INDEX4(i,2,q,3, numF,DIM,numQuadNodes)] = df2;

                new_dFdv[INDEX4(i,0,q,4, numF,DIM,numQuadNodes)] = df0;
                new_dFdv[INDEX4(i,1,q,4, numF,DIM,numQuadNodes)] = df1;
                new_dFdv[INDEX4(i,2,q,4, numF,DIM,numQuadNodes)] = df2;

                new_dFdv[INDEX4(i,0,q,5, numF,DIM,numQuadNodes)] = df0;
                new_dFdv[INDEX4(i,1,q,5, numF,DIM,numQuadNodes)] = df1;
                new_dFdv[INDEX4(i,2,q,5, numF,DIM,numQuadNodes)] = df2;

                new_dFdv[INDEX4(i,0,q,6, numF,DIM,numQuadNodes)] = df0;
                new_dFdv[INDEX4(i,1,q,6, numF,DIM,numQuadNodes)] = df1;
                new_dFdv[INDEX4(i,2,q,6, numF,DIM,numQuadNodes)] = df2;

                new_dFdv[INDEX4(i,0,q,7, numF,DIM,numQuadNodes)] = df0;
                new_dFdv[INDEX4(i,1,q,7, numF,DIM,numQuadNodes)] = df1;
                new_dFdv[INDEX4(i,2,q,7, numF,DIM,numQuadNodes)] = df2;
            }
        }
    } else {
        throw escript::ValueError("Quad_MacroHex: unable to create quadrature scheme for macro element.");
    }
    return numSubElements*numQuadNodes;
#undef DIM
}

} // namespace finley

