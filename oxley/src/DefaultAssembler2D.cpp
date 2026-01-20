
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
* Development from 2019 by School of Earth and Environmental Sciences
*
*****************************************************************************/

#include <oxley/DefaultAssembler2D.h>
#include <oxley/domainhelpers.h>

#include <escript/DataTypes.h>
#include <escript/index.h>

#ifdef ESYS_HAVE_TRILINOS
#include <trilinoswrap/TrilinosMatrixAdapter.h>
#endif

using namespace std;
using escript::AbstractSystemMatrix;
using escript::Data;

namespace oxley {

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
void DefaultAssembler2D<Scalar>::assemblePDESingle(AbstractSystemMatrix* mat, Data& rhs, 
                                      const Data& A, const Data& B,
                                      const Data& C, const Data& D,
                                      const Data& X, const Data& Y) const
{
    // Find the maximum level of refinement in the mesh
    int max_level = 0;
    for(p4est_topidx_t tree = domain->p4est->first_local_tree; tree <= domain->p4est->last_local_tree; tree++) {
        p4est_tree_t * tree_t = p4est_tree_array_index(domain->p4est->trees, tree);
        max_level = tree_t->maxlevel > max_level ? tree_t->maxlevel : max_level;
    }

    // This is similar to Ripley, except that in Oxley the quads vary in size so that the 
    // weightings need to be calcuated for each level of refinement
    const double SQRT3 = 1.73205080756887719318;
    double w[32][P4EST_MAXLEVEL] = {{0}};
// #pragma omp parallel for
    for(int i = 0; i <= max_level; i++)
    {
        double m_dx[2] = {domain->forestData.m_dx[0][P4EST_MAXLEVEL-i], 
                          domain->forestData.m_dx[1][P4EST_MAXLEVEL-i]};
        w[1][i]  = 1.0/24.0;
        w[5][i]  = -SQRT3/24 + 1.0/12;
        w[2][i]  = -SQRT3/24 - 1.0/12;
        w[19][i] = -m_dx[0]/12;
        w[11][i] = w[19][i]*(SQRT3 + 3)/12;
        w[14][i] = w[19][i]*(-SQRT3 + 3)/12;
        w[16][i] = w[19][i]*(5*SQRT3 + 9)/12;
        w[17][i] = w[19][i]*(-5*SQRT3 + 9)/12;
        w[27][i] = w[19][i]*(-SQRT3 - 3)/2;
        w[28][i] = w[19][i]*(SQRT3 - 3)/2;
        w[18][i] = -m_dx[1]/12;
        w[12][i] = w[18][i]*(5*SQRT3 + 9)/12;
        w[13][i] = w[18][i]*(-5*SQRT3 + 9)/12;
        w[10][i] = w[18][i]*(SQRT3 + 3)/12;
        w[15][i] = w[18][i]*(-SQRT3 + 3)/12;
        w[25][i] = w[18][i]*(-SQRT3 - 3)/2;
        w[26][i] = w[18][i]*(SQRT3 - 3)/2;
        w[22][i] = m_dx[0]*m_dx[1]/144;
        w[20][i] = w[22][i]*(SQRT3 + 2);
        w[21][i] = w[22][i]*(-SQRT3 + 2);
        w[23][i] = w[22][i]*(4*SQRT3 + 7);
        w[24][i] = w[22][i]*(-4*SQRT3 + 7);
        w[3][i]  = m_dx[0]/(24*m_dx[1]);
        w[7][i]  = w[3][i]*(SQRT3 + 2);
        w[8][i]  = w[3][i]*(-SQRT3 + 2);
        w[6][i]  = -m_dx[1]/(24*m_dx[0]);
        w[0][i]  = w[6][i]*(SQRT3 + 2);
        w[4][i]  = w[6][i]*(-SQRT3 + 2);
    }

    // Debugging info
    // for(int i = 0; i < 28; i++)
    //     std::cout << i << ": " << w[i][0] << std::endl;

    const bool addEM_S = (!A.isEmpty() || !B.isEmpty() || !C.isEmpty() || !D.isEmpty());
    const bool addEM_F = (!X.isEmpty() || !Y.isEmpty());

    // Matrices
    const Scalar zero = static_cast<Scalar>(0);
    vector<Scalar> EM_S(4*4, zero);
    vector<Scalar> EM_F(4, zero);

    rhs.requireWrite();

// #pragma omp parallel for
    for (p4est_topidx_t t = domain->p4est->first_local_tree; t <= domain->p4est->last_local_tree; t++) // Loop over every tree
    {
        p4est_tree_t * currenttree = p4est_tree_array_index(domain->p4est->trees, t);
        sc_array_t * tquadrants = &currenttree->quadrants;
        p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
        for (int q = 0; q < Q; ++q)  
        {
            if (addEM_S)
                fill(EM_S.begin(), EM_S.end(), zero);
            if (addEM_F)
                fill(EM_F.begin(), EM_F.end(), zero);
            
            p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, q);
            int l = quad->level;
            double xy[3];
            p4est_qcoord_to_vertex(domain->p4est->connectivity, t, quad->x, quad->y, xy);
            long id = domain->getQuadID(domain->NodeIDs.find(std::make_pair(xy[0],xy[1]))->second);

            ///////////////
            // process A //
            ///////////////
            if (!A.isEmpty()) {
                const Scalar* A_p = A.getSampleDataRO(id, zero);
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
                    const Scalar tmp0  = w[3][l]*(A_11_0 + A_11_1 + A_11_2 + A_11_3);
                    const Scalar tmp1  = w[1][l]*(A_01_0 + A_01_3 - A_10_1 - A_10_2);
                    const Scalar tmp2  = w[4][l]*(A_00_2 + A_00_3);
                    const Scalar tmp3  = w[0][l]*(A_00_0 + A_00_1);
                    const Scalar tmp4  = w[5][l]*(A_01_2 - A_10_3);
                    const Scalar tmp5  = w[2][l]*(-A_01_1 + A_10_0);
                    const Scalar tmp6  = w[5][l]*(A_01_3 + A_10_0);
                    const Scalar tmp7  = w[3][l]*(-A_11_0 - A_11_1 - A_11_2 - A_11_3);
                    const Scalar tmp8  = w[6][l]*(A_00_0 + A_00_1 + A_00_2 + A_00_3);
                    const Scalar tmp9  = w[1][l]*(A_01_1 + A_01_2 + A_10_1 + A_10_2);
                    const Scalar tmp10 = w[2][l]*(-A_01_0 - A_10_3);
                    const Scalar tmp11 = w[4][l]*(A_00_0 + A_00_1);
                    const Scalar tmp12 = w[0][l]*(A_00_2 + A_00_3);
                    const Scalar tmp13 = w[5][l]*(A_01_1 - A_10_0);
                    const Scalar tmp14 = w[2][l]*(-A_01_2 + A_10_3);
                    const Scalar tmp15 = w[7][l]*(A_11_0 + A_11_2);
                    const Scalar tmp16 = w[4][l]*(-A_00_2 - A_00_3);
                    const Scalar tmp17 = w[0][l]*(-A_00_0 - A_00_1);
                    const Scalar tmp18 = w[5][l]*(A_01_3 + A_10_3);
                    const Scalar tmp19 = w[8][l]*(A_11_1 + A_11_3);
                    const Scalar tmp20 = w[2][l]*(-A_01_0 - A_10_0);
                    const Scalar tmp21 = w[7][l]*(A_11_1 + A_11_3);
                    const Scalar tmp22 = w[4][l]*(-A_00_0 - A_00_1);
                    const Scalar tmp23 = w[0][l]*(-A_00_2 - A_00_3);
                    const Scalar tmp24 = w[5][l]*(A_01_0 + A_10_0);
                    const Scalar tmp25 = w[8][l]*(A_11_0 + A_11_2);
                    const Scalar tmp26 = w[2][l]*(-A_01_3 - A_10_3);
                    const Scalar tmp27 = w[5][l]*(-A_01_1 - A_10_2);
                    const Scalar tmp28 = w[1][l]*(-A_01_0 - A_01_3 - A_10_0 - A_10_3);
                    const Scalar tmp29 = w[2][l]*(A_01_2 + A_10_1);
                    const Scalar tmp30 = w[7][l]*(-A_11_1 - A_11_3);
                    const Scalar tmp31 = w[1][l]*(-A_01_1 - A_01_2 + A_10_0 + A_10_3);
                    const Scalar tmp32 = w[5][l]*(-A_01_0 + A_10_2);
                    const Scalar tmp33 = w[8][l]*(-A_11_0 - A_11_2);
                    const Scalar tmp34 = w[6][l]*(-A_00_0 - A_00_1 - A_00_2 - A_00_3);
                    const Scalar tmp35 = w[2][l]*(A_01_3 - A_10_1);
                    const Scalar tmp36 = w[5][l]*(A_01_0 + A_10_3);
                    const Scalar tmp37 = w[2][l]*(-A_01_3 - A_10_0);
                    const Scalar tmp38 = w[7][l]*(-A_11_0 - A_11_2);
                    const Scalar tmp39 = w[5][l]*(-A_01_3 + A_10_1);
                    const Scalar tmp40 = w[8][l]*(-A_11_1 - A_11_3);
                    const Scalar tmp41 = w[2][l]*(A_01_0 - A_10_2);
                    const Scalar tmp42 = w[5][l]*(A_01_1 - A_10_3);
                    const Scalar tmp43 = w[2][l]*(-A_01_2 + A_10_0);
                    const Scalar tmp44 = w[5][l]*(A_01_2 - A_10_0);
                    const Scalar tmp45 = w[2][l]*(-A_01_1 + A_10_3);
                    const Scalar tmp46 = w[5][l]*(-A_01_0 + A_10_1);
                    const Scalar tmp47 = w[2][l]*(A_01_3 - A_10_2);
                    const Scalar tmp48 = w[5][l]*(-A_01_1 - A_10_1);
                    const Scalar tmp49 = w[2][l]*(A_01_2 + A_10_2);
                    const Scalar tmp50 = w[5][l]*(-A_01_3 + A_10_2);
                    const Scalar tmp51 = w[2][l]*(A_01_0 - A_10_1);
                    const Scalar tmp52 = w[5][l]*(-A_01_2 - A_10_1);
                    const Scalar tmp53 = w[2][l]*(A_01_1 + A_10_2);
                    const Scalar tmp54 = w[5][l]*(-A_01_2 - A_10_2);
                    const Scalar tmp55 = w[2][l]*(A_01_1 + A_10_1);
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
                    const Scalar tmp0 = 6.*w[1][l]*(A_01 - A_10);
                    const Scalar tmp1 = 6.*w[1][l]*(A_01 + A_10);
                    const Scalar tmp2 = 6.*w[1][l]*(-A_01 - A_10);
                    const Scalar tmp3 = 6.*w[1][l]*(-A_01 + A_10);
                    EM_S[INDEX2(0,0,4)]+=-8.*A_00*w[6][l] + 8.*A_11*w[3][l] + tmp1;
                    EM_S[INDEX2(0,1,4)]+= 8.*A_00*w[6][l] + 4.*A_11*w[3][l] + tmp0;
                    EM_S[INDEX2(0,2,4)]+=-4.*A_00*w[6][l] - 8.*A_11*w[3][l] + tmp3;
                    EM_S[INDEX2(0,3,4)]+= 4.*A_00*w[6][l] - 4.*A_11*w[3][l] + tmp2;
                    EM_S[INDEX2(1,0,4)]+= 8.*A_00*w[6][l] + 4.*A_11*w[3][l] + tmp3;
                    EM_S[INDEX2(1,1,4)]+=-8.*A_00*w[6][l] + 8.*A_11*w[3][l] + tmp2;
                    EM_S[INDEX2(1,2,4)]+= 4.*A_00*w[6][l] - 4.*A_11*w[3][l] + tmp1;
                    EM_S[INDEX2(1,3,4)]+=-4.*A_00*w[6][l] - 8.*A_11*w[3][l] + tmp0;
                    EM_S[INDEX2(2,0,4)]+=-4.*A_00*w[6][l] - 8.*A_11*w[3][l] + tmp0;
                    EM_S[INDEX2(2,1,4)]+= 4.*A_00*w[6][l] - 4.*A_11*w[3][l] + tmp1;
                    EM_S[INDEX2(2,2,4)]+=-8.*A_00*w[6][l] + 8.*A_11*w[3][l] + tmp2;
                    EM_S[INDEX2(2,3,4)]+= 8.*A_00*w[6][l] + 4.*A_11*w[3][l] + tmp3;
                    EM_S[INDEX2(3,0,4)]+= 4.*A_00*w[6][l] - 4.*A_11*w[3][l] + tmp2;
                    EM_S[INDEX2(3,1,4)]+=-4.*A_00*w[6][l] - 8.*A_11*w[3][l] + tmp3;
                    EM_S[INDEX2(3,2,4)]+= 8.*A_00*w[6][l] + 4.*A_11*w[3][l] + tmp0;
                    EM_S[INDEX2(3,3,4)]+=-8.*A_00*w[6][l] + 8.*A_11*w[3][l] + tmp1;
                }
            }

            ///////////////
            // process B //
            ///////////////
            if (!B.isEmpty()) {
                const Scalar* B_p = B.getSampleDataRO(id, zero);
                if (B.actsExpanded()) {
                    const Scalar B_0_0 = B_p[INDEX2(0,0,2)];
                    const Scalar B_1_0 = B_p[INDEX2(1,0,2)];
                    const Scalar B_0_1 = B_p[INDEX2(0,1,2)];
                    const Scalar B_1_1 = B_p[INDEX2(1,1,2)];
                    const Scalar B_0_2 = B_p[INDEX2(0,2,2)];
                    const Scalar B_1_2 = B_p[INDEX2(1,2,2)];
                    const Scalar B_0_3 = B_p[INDEX2(0,3,2)];
                    const Scalar B_1_3 = B_p[INDEX2(1,3,2)];
                    const Scalar tmp0  = w[11][l]*(B_1_0 + B_1_1);
                    const Scalar tmp1  = w[14][l]*(B_1_2 + B_1_3);
                    const Scalar tmp2  = w[15][l]*(-B_0_1 - B_0_3);
                    const Scalar tmp3  = w[10][l]*(-B_0_0 - B_0_2);
                    const Scalar tmp4  = w[11][l]*(B_1_2 + B_1_3);
                    const Scalar tmp5  = w[14][l]*(B_1_0 + B_1_1);
                    const Scalar tmp6  = w[11][l]*(-B_1_2 - B_1_3);
                    const Scalar tmp7  = w[14][l]*(-B_1_0 - B_1_1);
                    const Scalar tmp8  = w[11][l]*(-B_1_0 - B_1_1);
                    const Scalar tmp9  = w[14][l]*(-B_1_2 - B_1_3);
                    const Scalar tmp10 = w[10][l]*(-B_0_1 - B_0_3);
                    const Scalar tmp11 = w[15][l]*(-B_0_0 - B_0_2);
                    const Scalar tmp12 = w[15][l]*(B_0_0 + B_0_2);
                    const Scalar tmp13 = w[10][l]*(B_0_1 + B_0_3);
                    const Scalar tmp14 = w[10][l]*(B_0_0 + B_0_2);
                    const Scalar tmp15 = w[15][l]*(B_0_1 + B_0_3);
                    EM_S[INDEX2(0,0,4)]+=B_0_0*w[12][l] + B_0_1*w[10][l] + B_0_2*w[15][l] + B_0_3*w[13][l] + B_1_0*w[16][l] + B_1_1*w[14][l] + B_1_2*w[11][l] + B_1_3*w[17][l];
                    EM_S[INDEX2(0,1,4)]+=B_0_0*w[10][l] + B_0_1*w[12][l] + B_0_2*w[13][l] + B_0_3*w[15][l] + tmp0 + tmp1;
                    EM_S[INDEX2(0,2,4)]+=B_1_0*w[11][l] + B_1_1*w[17][l] + B_1_2*w[16][l] + B_1_3*w[14][l] + tmp14 + tmp15;
                    EM_S[INDEX2(0,3,4)]+=tmp12 + tmp13 + tmp4 + tmp5;
                    EM_S[INDEX2(1,0,4)]+=-B_0_0*w[12][l] - B_0_1*w[10][l] - B_0_2*w[15][l] - B_0_3*w[13][l] + tmp0 + tmp1;
                    EM_S[INDEX2(1,1,4)]+=-B_0_0*w[10][l] - B_0_1*w[12][l] - B_0_2*w[13][l] - B_0_3*w[15][l] + B_1_0*w[14][l] + B_1_1*w[16][l] + B_1_2*w[17][l] + B_1_3*w[11][l];
                    EM_S[INDEX2(1,2,4)]+=tmp2 + tmp3 + tmp4 + tmp5;
                    EM_S[INDEX2(1,3,4)]+= B_1_0*w[17][l] + B_1_1*w[11][l] + B_1_2*w[14][l] + B_1_3*w[16][l] + tmp10 + tmp11;
                    EM_S[INDEX2(2,0,4)]+=-B_1_0*w[16][l] - B_1_1*w[14][l] - B_1_2*w[11][l] - B_1_3*w[17][l] + tmp14 + tmp15;
                    EM_S[INDEX2(2,1,4)]+=tmp12 + tmp13 + tmp8 + tmp9;
                    EM_S[INDEX2(2,2,4)]+=B_0_0*w[15][l] + B_0_1*w[13][l] + B_0_2*w[12][l] + B_0_3*w[10][l] - B_1_0*w[11][l] - B_1_1*w[17][l] - B_1_2*w[16][l] - B_1_3*w[14][l];
                    EM_S[INDEX2(2,3,4)]+=B_0_0*w[13][l] + B_0_1*w[15][l] + B_0_2*w[10][l] + B_0_3*w[12][l] + tmp6 + tmp7;
                    EM_S[INDEX2(3,0,4)]+=tmp2 + tmp3 + tmp8 + tmp9;
                    EM_S[INDEX2(3,1,4)]+=-B_1_0*w[14][l] - B_1_1*w[16][l] - B_1_2*w[17][l] - B_1_3*w[11][l] + tmp10 + tmp11;
                    EM_S[INDEX2(3,2,4)]+=-B_0_0*w[15][l] - B_0_1*w[13][l] - B_0_2*w[12][l] - B_0_3*w[10][l] + tmp6 + tmp7;
                    EM_S[INDEX2(3,3,4)]+=-B_0_0*w[13][l] - B_0_1*w[15][l] - B_0_2*w[10][l] - B_0_3*w[12][l] - B_1_0*w[17][l] - B_1_1*w[11][l] - B_1_2*w[14][l] - B_1_3*w[16][l];
                } else { // constant data
                    const Scalar B_0 = B_p[0];
                    const Scalar B_1 = B_p[1];
                    EM_S[INDEX2(0,0,4)]+= 2.*B_0*w[18][l] + 2.*B_1*w[19][l];
                    EM_S[INDEX2(0,1,4)]+= 2.*B_0*w[18][l] +    B_1*w[19][l];
                    EM_S[INDEX2(0,2,4)]+=    B_0*w[18][l] + 2.*B_1*w[19][l];
                    EM_S[INDEX2(0,3,4)]+=    B_0*w[18][l] +    B_1*w[19][l];
                    EM_S[INDEX2(1,0,4)]+=-2.*B_0*w[18][l] +    B_1*w[19][l];
                    EM_S[INDEX2(1,1,4)]+=-2.*B_0*w[18][l] + 2.*B_1*w[19][l];
                    EM_S[INDEX2(1,2,4)]+=   -B_0*w[18][l] +    B_1*w[19][l];
                    EM_S[INDEX2(1,3,4)]+=   -B_0*w[18][l] + 2.*B_1*w[19][l];
                    EM_S[INDEX2(2,0,4)]+=    B_0*w[18][l] - 2.*B_1*w[19][l];
                    EM_S[INDEX2(2,1,4)]+=    B_0*w[18][l] -    B_1*w[19][l];
                    EM_S[INDEX2(2,2,4)]+= 2.*B_0*w[18][l] - 2.*B_1*w[19][l];
                    EM_S[INDEX2(2,3,4)]+= 2.*B_0*w[18][l] -    B_1*w[19][l];
                    EM_S[INDEX2(3,0,4)]+=   -B_0*w[18][l] -    B_1*w[19][l];
                    EM_S[INDEX2(3,1,4)]+=   -B_0*w[18][l] - 2.*B_1*w[19][l];
                    EM_S[INDEX2(3,2,4)]+=-2.*B_0*w[18][l] -    B_1*w[19][l];
                    EM_S[INDEX2(3,3,4)]+=-2.*B_0*w[18][l] - 2.*B_1*w[19][l];
                }
            }

            ///////////////
            // process C //
            ///////////////
            if (!C.isEmpty()) {
                const Scalar* C_p = C.getSampleDataRO(id, zero);
                if (C.actsExpanded()) {
                    const Scalar C_0_0 = C_p[INDEX2(0,0,2)];
                    const Scalar C_1_0 = C_p[INDEX2(1,0,2)];
                    const Scalar C_0_1 = C_p[INDEX2(0,1,2)];
                    const Scalar C_1_1 = C_p[INDEX2(1,1,2)];
                    const Scalar C_0_2 = C_p[INDEX2(0,2,2)];
                    const Scalar C_1_2 = C_p[INDEX2(1,2,2)];
                    const Scalar C_0_3 = C_p[INDEX2(0,3,2)];
                    const Scalar C_1_3 = C_p[INDEX2(1,3,2)];
                    const Scalar tmp0 = w[11][l]*(C_1_0 + C_1_1);
                    const Scalar tmp1 = w[14][l]*(C_1_2 + C_1_3);
                    const Scalar tmp2 = w[15][l]*(C_0_0 + C_0_2);
                    const Scalar tmp3 = w[10][l]*(C_0_1 + C_0_3);
                    const Scalar tmp4 = w[11][l]*(-C_1_0 - C_1_1);
                    const Scalar tmp5 = w[14][l]*(-C_1_2 - C_1_3);
                    const Scalar tmp6 = w[11][l]*(-C_1_2 - C_1_3);
                    const Scalar tmp7 = w[14][l]*(-C_1_0 - C_1_1);
                    const Scalar tmp8 = w[11][l]*(C_1_2 + C_1_3);
                    const Scalar tmp9 = w[14][l]*(C_1_0 + C_1_1);
                    const Scalar tmp10 = w[10][l]*(-C_0_1 - C_0_3);
                    const Scalar tmp11 = w[15][l]*(-C_0_0 - C_0_2);
                    const Scalar tmp12 = w[15][l]*(-C_0_1 - C_0_3);
                    const Scalar tmp13 = w[10][l]*(-C_0_0 - C_0_2);
                    const Scalar tmp14 = w[10][l]*(C_0_0 + C_0_2);
                    const Scalar tmp15 = w[15][l]*(C_0_1 + C_0_3);
                    EM_S[INDEX2(0,0,4)]+= C_0_0*w[12][l] + C_0_1*w[10][l] + C_0_2*w[15][l] + C_0_3*w[13][l] + C_1_0*w[16][l] + C_1_1*w[14][l] + C_1_2*w[11][l] + C_1_3*w[17][l];
                    EM_S[INDEX2(0,1,4)]+=-C_0_0*w[12][l] - C_0_1*w[10][l] - C_0_2*w[15][l] - C_0_3*w[13][l] + tmp0 + tmp1;
                    EM_S[INDEX2(0,2,4)]+=-C_1_0*w[16][l] - C_1_1*w[14][l] - C_1_2*w[11][l] - C_1_3*w[17][l] + tmp14 + tmp15;
                    EM_S[INDEX2(0,3,4)]+=tmp12 + tmp13 + tmp4 + tmp5;
                    EM_S[INDEX2(1,0,4)]+= C_0_0*w[10][l] + C_0_1*w[12][l] + C_0_2*w[13][l] + C_0_3*w[15][l] + tmp0 + tmp1;
                    EM_S[INDEX2(1,1,4)]+=-C_0_0*w[10][l] - C_0_1*w[12][l] - C_0_2*w[13][l] - C_0_3*w[15][l] + C_1_0*w[14][l] + C_1_1*w[16][l] + C_1_2*w[17][l] + C_1_3*w[11][l];
                    EM_S[INDEX2(1,2,4)]+=tmp2 + tmp3 + tmp4 + tmp5;
                    EM_S[INDEX2(1,3,4)]+=-C_1_0*w[14][l] - C_1_1*w[16][l] - C_1_2*w[17][l] - C_1_3*w[11][l] + tmp10 + tmp11;
                    EM_S[INDEX2(2,0,4)]+= C_1_0*w[11][l] + C_1_1*w[17][l] + C_1_2*w[16][l] + C_1_3*w[14][l] + tmp14 + tmp15;
                    EM_S[INDEX2(2,1,4)]+=tmp12 + tmp13 + tmp8 + tmp9;
                    EM_S[INDEX2(2,2,4)]+= C_0_0*w[15][l] + C_0_1*w[13][l] + C_0_2*w[12][l] + C_0_3*w[10][l] - C_1_0*w[11][l] - C_1_1*w[17][l] - C_1_2*w[16][l] - C_1_3*w[14][l];
                    EM_S[INDEX2(2,3,4)]+=-C_0_0*w[15][l] - C_0_1*w[13][l] - C_0_2*w[12][l] - C_0_3*w[10][l] + tmp6 + tmp7;
                    EM_S[INDEX2(3,0,4)]+=tmp2 + tmp3 + tmp8 + tmp9;
                    EM_S[INDEX2(3,1,4)]+= C_1_0*w[17][l] + C_1_1*w[11][l] + C_1_2*w[14][l] + C_1_3*w[16][l] + tmp10 + tmp11;
                    EM_S[INDEX2(3,2,4)]+= C_0_0*w[13][l] + C_0_1*w[15][l] + C_0_2*w[10][l] + C_0_3*w[12][l] + tmp6 + tmp7;
                    EM_S[INDEX2(3,3,4)]+=-C_0_0*w[13][l] - C_0_1*w[15][l] - C_0_2*w[10][l] - C_0_3*w[12][l] - C_1_0*w[17][l] - C_1_1*w[11][l] - C_1_2*w[14][l] - C_1_3*w[16][l];
                } else { // constant data
                    const Scalar C_0 = C_p[0];
                    const Scalar C_1 = C_p[1];
                    EM_S[INDEX2(0,0,4)]+= 2.*C_0*w[18][l] + 2.*C_1*w[19][l];
                    EM_S[INDEX2(0,1,4)]+=-2.*C_0*w[18][l] +    C_1*w[19][l];
                    EM_S[INDEX2(0,2,4)]+=    C_0*w[18][l] - 2.*C_1*w[19][l];
                    EM_S[INDEX2(0,3,4)]+=   -C_0*w[18][l] -    C_1*w[19][l];
                    EM_S[INDEX2(1,0,4)]+= 2.*C_0*w[18][l] +    C_1*w[19][l];
                    EM_S[INDEX2(1,1,4)]+=-2.*C_0*w[18][l] + 2.*C_1*w[19][l];
                    EM_S[INDEX2(1,2,4)]+=    C_0*w[18][l] -    C_1*w[19][l];
                    EM_S[INDEX2(1,3,4)]+=   -C_0*w[18][l] - 2.*C_1*w[19][l];
                    EM_S[INDEX2(2,0,4)]+=    C_0*w[18][l] + 2.*C_1*w[19][l];
                    EM_S[INDEX2(2,1,4)]+=   -C_0*w[18][l] +    C_1*w[19][l];
                    EM_S[INDEX2(2,2,4)]+= 2.*C_0*w[18][l] - 2.*C_1*w[19][l];
                    EM_S[INDEX2(2,3,4)]+=-2.*C_0*w[18][l] -    C_1*w[19][l];
                    EM_S[INDEX2(3,0,4)]+=    C_0*w[18][l] +    C_1*w[19][l];
                    EM_S[INDEX2(3,1,4)]+=   -C_0*w[18][l] + 2.*C_1*w[19][l];
                    EM_S[INDEX2(3,2,4)]+= 2.*C_0*w[18][l] -    C_1*w[19][l];
                    EM_S[INDEX2(3,3,4)]+=-2.*C_0*w[18][l] - 2.*C_1*w[19][l];
                }
            }

            ///////////////
            // process D //
            ///////////////
            if (!D.isEmpty()) {
                const Scalar* D_p = D.getSampleDataRO(id, zero);
                if (D.actsExpanded()) {
                    const Scalar D_0 = D_p[0];
                    const Scalar D_1 = D_p[1];
                    const Scalar D_2 = D_p[2];
                    const Scalar D_3 = D_p[3];
                    const Scalar tmp0 =  w[21][l]*(D_2 + D_3);
                    const Scalar tmp1 =  w[20][l]*(D_0 + D_1);
                    const Scalar tmp2 =  w[22][l]*(D_0 + D_1 + D_2 + D_3);
                    const Scalar tmp3 =  w[21][l]*(D_0 + D_1);
                    const Scalar tmp4 =  w[20][l]*(D_2 + D_3);
                    const Scalar tmp5 =  w[22][l]*(D_1 + D_2);
                    const Scalar tmp6 =  w[21][l]*(D_0 + D_2);
                    const Scalar tmp7 =  w[20][l]*(D_1 + D_3);
                    const Scalar tmp8 =  w[21][l]*(D_1 + D_3);
                    const Scalar tmp9 =  w[20][l]*(D_0 + D_2);
                    const Scalar tmp10 = w[22][l]*(D_0 + D_3);
                    EM_S[INDEX2(0,0,4)]+=D_0*w[23][l] + D_3*w[24][l] + tmp5;
                    EM_S[INDEX2(0,1,4)]+=tmp0 + tmp1;
                    EM_S[INDEX2(0,2,4)]+=tmp8 + tmp9;
                    EM_S[INDEX2(0,3,4)]+=tmp2;
                    EM_S[INDEX2(1,0,4)]+=tmp0 + tmp1;
                    EM_S[INDEX2(1,1,4)]+=D_1*w[23][l] + D_2*w[24][l] + tmp10;
                    EM_S[INDEX2(1,2,4)]+=tmp2;
                    EM_S[INDEX2(1,3,4)]+=tmp6 + tmp7;
                    EM_S[INDEX2(2,0,4)]+=tmp8 + tmp9;
                    EM_S[INDEX2(2,1,4)]+=tmp2;
                    EM_S[INDEX2(2,2,4)]+=D_1*w[24][l] + D_2*w[23][l] + tmp10;
                    EM_S[INDEX2(2,3,4)]+=tmp3 + tmp4;
                    EM_S[INDEX2(3,0,4)]+=tmp2;
                    EM_S[INDEX2(3,1,4)]+=tmp6 + tmp7;
                    EM_S[INDEX2(3,2,4)]+=tmp3 + tmp4;
                    EM_S[INDEX2(3,3,4)]+=D_0*w[24][l] + D_3*w[23][l] + tmp5;
                } else { // constant data
                    const Scalar D_0 = D_p[0];
                    EM_S[INDEX2(0,0,4)]+=16.*D_0*w[22][l];
                    EM_S[INDEX2(0,1,4)]+= 8.*D_0*w[22][l];
                    EM_S[INDEX2(0,2,4)]+= 8.*D_0*w[22][l];
                    EM_S[INDEX2(0,3,4)]+= 4.*D_0*w[22][l];
                    EM_S[INDEX2(1,0,4)]+= 8.*D_0*w[22][l];
                    EM_S[INDEX2(1,1,4)]+=16.*D_0*w[22][l];
                    EM_S[INDEX2(1,2,4)]+= 4.*D_0*w[22][l];
                    EM_S[INDEX2(1,3,4)]+= 8.*D_0*w[22][l];
                    EM_S[INDEX2(2,0,4)]+= 8.*D_0*w[22][l];
                    EM_S[INDEX2(2,1,4)]+= 4.*D_0*w[22][l];
                    EM_S[INDEX2(2,2,4)]+=16.*D_0*w[22][l];
                    EM_S[INDEX2(2,3,4)]+= 8.*D_0*w[22][l];
                    EM_S[INDEX2(3,0,4)]+= 4.*D_0*w[22][l];
                    EM_S[INDEX2(3,1,4)]+= 8.*D_0*w[22][l];
                    EM_S[INDEX2(3,2,4)]+= 8.*D_0*w[22][l];
                    EM_S[INDEX2(3,3,4)]+=16.*D_0*w[22][l];
                }
            }

            ///////////////
            // process X //
            ///////////////
            if (!X.isEmpty()) {
                const Scalar* X_p = X.getSampleDataRO(id, zero);
                if (X.actsExpanded()) {
                    const Scalar X_0_0 = X_p[INDEX2(0,0,2)];
                    const Scalar X_1_0 = X_p[INDEX2(1,0,2)];
                    const Scalar X_0_1 = X_p[INDEX2(0,1,2)];
                    const Scalar X_1_1 = X_p[INDEX2(1,1,2)];
                    const Scalar X_0_2 = X_p[INDEX2(0,2,2)];
                    const Scalar X_1_2 = X_p[INDEX2(1,2,2)];
                    const Scalar X_0_3 = X_p[INDEX2(0,3,2)];
                    const Scalar X_1_3 = X_p[INDEX2(1,3,2)];
                    const Scalar tmp0 = 6.*w[15][l]*(X_0_2 + X_0_3);
                    const Scalar tmp1 = 6.*w[10][l]*(X_0_0 + X_0_1);
                    const Scalar tmp2 = 6.*w[11][l]*(X_1_0 + X_1_2);
                    const Scalar tmp3 = 6.*w[14][l]*(X_1_1 + X_1_3);
                    const Scalar tmp4 = 6.*w[11][l]*(X_1_1 + X_1_3);
                    const Scalar tmp5 = w[25][l]*(X_0_0 + X_0_1);
                    const Scalar tmp6 = w[26][l]*(X_0_2 + X_0_3);
                    const Scalar tmp7 = 6.*w[14][l]*(X_1_0 + X_1_2);
                    const Scalar tmp8 =  w[27][l]*(X_1_0 + X_1_2);
                    const Scalar tmp9 =  w[28][l]*(X_1_1 + X_1_3);
                    const Scalar tmp10 = w[25][l]*(-X_0_2 - X_0_3);
                    const Scalar tmp11 = w[26][l]*(-X_0_0 - X_0_1);
                    const Scalar tmp12 = w[27][l]*(X_1_1 + X_1_3);
                    const Scalar tmp13 = w[28][l]*(X_1_0 + X_1_2);
                    const Scalar tmp14 = w[25][l]*(X_0_2 + X_0_3);
                    const Scalar tmp15 = w[26][l]*(X_0_0 + X_0_1);
                    EM_F[0]+=tmp0 + tmp1 + tmp2 + tmp3;
                    EM_F[1]+=tmp4 + tmp5 + tmp6 + tmp7;
                    EM_F[2]+=tmp10 + tmp11 + tmp8 + tmp9;
                    EM_F[3]+=tmp12 + tmp13 + tmp14 + tmp15;
                } else { // constant data
                    const Scalar X_0 = X_p[0];
                    const Scalar X_1 = X_p[1];
                    EM_F[0]+= 6.*X_0*w[18][l] + 6.*X_1*w[19][l];
                    EM_F[1]+=-6.*X_0*w[18][l] + 6.*X_1*w[19][l];
                    EM_F[2]+= 6.*X_0*w[18][l] - 6.*X_1*w[19][l];
                    EM_F[3]+=-6.*X_0*w[18][l] - 6.*X_1*w[19][l];
                }
            }

            ///////////////
            // process Y //
            ///////////////
            if (!Y.isEmpty()) {
                const Scalar* Y_p = Y.getSampleDataRO(id, zero);
                if (Y.actsExpanded()) {
                    const Scalar Y_0 = Y_p[0];
                    const Scalar Y_1 = Y_p[1];
                    const Scalar Y_2 = Y_p[2];
                    const Scalar Y_3 = Y_p[3];
                    const Scalar tmp0 = 6.*w[22][l]*(Y_1 + Y_2);
                    const Scalar tmp1 = 6.*w[22][l]*(Y_0 + Y_3);
                    EM_F[0]+=6.*Y_0*w[20][l] + 6.*Y_3*w[21][l] + tmp0;
                    EM_F[1]+=6.*Y_1*w[20][l] + 6.*Y_2*w[21][l] + tmp1;
                    EM_F[2]+=6.*Y_1*w[21][l] + 6.*Y_2*w[20][l] + tmp1;
                    EM_F[3]+=6.*Y_0*w[21][l] + 6.*Y_3*w[20][l] + tmp0;
                } else { // constant data
                    EM_F[0]+=36.*Y_p[0]*w[22][l];
                    EM_F[1]+=36.*Y_p[0]*w[22][l];
                    EM_F[2]+=36.*Y_p[0]*w[22][l];
                    EM_F[3]+=36.*Y_p[0]*w[22][l];
                }
            }
            // add to matrix (if addEM_S) and RHS (if addEM_F)
            domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, addEM_S, addEM_F, q, t);
        }
    }
} // end of parallel region

/****************************************************************************/
// PDE SINGLE BOUNDARY
/****************************************************************************/

template<class Scalar>
void DefaultAssembler2D<Scalar>::assemblePDEBoundarySingle(
                                        AbstractSystemMatrix* mat, Data& rhs,
                                        const Data& d, const Data& y) const
{
    // Find the maximum level of refinement in the mesh
    int max_level = 0;
    for(p4est_topidx_t tree = domain->p4est->first_local_tree; tree <= domain->p4est->last_local_tree; tree++) {
        p4est_tree_t * tree_t = p4est_tree_array_index(domain->p4est->trees, tree);
        max_level = tree_t->maxlevel > max_level ? tree_t->maxlevel : max_level;
    }

    // This is similar to Ripley, except that in Oxley the quads vary in size so that the 
    // weightings need to be calculated for each level of refinement
    const double SQRT3 = 1.73205080756887719318;
    double w[16][P4EST_MAXLEVEL] = {{0}};
#pragma omp parallel for
    for(int i = 0; i <= max_level; i++)
    {
        double m_dx[2] = {domain->forestData.m_dx[0][P4EST_MAXLEVEL-i], 
                          domain->forestData.m_dx[1][P4EST_MAXLEVEL-i]};

        w[5][i] = m_dx[0]/12;
        w[6][i] = w[5][i]*(SQRT3 + 2);
        w[7][i] = w[5][i]*(-SQRT3 + 2);
        w[8][i] = w[5][i]*(SQRT3 + 3);
        w[9][i] = w[5][i]*(-SQRT3 + 3);
        w[2][i] = m_dx[1]/12;
        w[0][i] = w[2][i]*(SQRT3 + 2);
        w[1][i] = w[2][i]*(-SQRT3 + 2);
        w[3][i] = w[2][i]*(SQRT3 + 3);
        w[4][i] = w[2][i]*(-SQRT3 + 3);
    }
    
    const bool addEM_S = !d.isEmpty();
    const bool addEM_F = !y.isEmpty();
    const Scalar zero = static_cast<Scalar>(0);
    rhs.requireWrite();

    vector<Scalar> EM_S(4*4);
    vector<Scalar> EM_F(4);

    if(domain->m_faceOffset[0] > -1)
    {
        if (addEM_S)
            fill(EM_S.begin(), EM_S.end(), zero);
        if (addEM_F) {
            EM_F[1] = zero;
            EM_F[3] = zero;
        }

        for (index_t k=0; k<domain->NodeIDsLeft.size()-1; k++) {
            int id = domain->m_faceOffset[0]+k;
            int l = domain->NodeIDsLeft[k].level;

            #ifdef OXLEY_ENABLE_DEBUG_ASSEMBLE_PDE_BOUNDARY
            std::cout << "Left\tid: " << id << " , level: " << l << std::endl;
            #endif

            ///////////////
            // process d //
            ///////////////
            if (addEM_S) {
                const Scalar* d_p = d.getSampleDataRO(id, zero);
                if (d.actsExpanded()) {
                    const Scalar d_0 = d_p[0];
                    const Scalar d_1 = d_p[1];
                    const Scalar tmp0 = w[2][l]*(d_0 + d_1);
                    EM_S[INDEX2(0,0,4)] = d_0*w[0][l] + d_1*w[1][l];
                    EM_S[INDEX2(2,0,4)] = tmp0;
                    EM_S[INDEX2(0,2,4)] = tmp0;
                    EM_S[INDEX2(2,2,4)] = d_0*w[1][l] + d_1*w[0][l];
                } else { // constant data
                    EM_S[INDEX2(0,0,4)] = 4.*d_p[0]*w[2][l];
                    EM_S[INDEX2(2,0,4)] = 2.*d_p[0]*w[2][l];
                    EM_S[INDEX2(0,2,4)] = 2.*d_p[0]*w[2][l];
                    EM_S[INDEX2(2,2,4)] = 4.*d_p[0]*w[2][l];
                }
            }

            ///////////////
            // process y //
            ///////////////
            if (addEM_F) {
                const Scalar* y_p = y.getSampleDataRO(id, zero);
                if (y.actsExpanded()) {
                    EM_F[0] = w[3][l]*y_p[0] + w[4][l]*y_p[1];
                    EM_F[2] = w[3][l]*y_p[1] + w[4][l]*y_p[0];
                } else { // constant data
                    EM_F[0] = 6.*w[2][l]*y_p[0];
                    EM_F[2] = 6.*w[2][l]*y_p[0];

                    #ifdef OXLEY_ENABLE_DEBUG_ASSEMBLE_PDE_BOUNDARY
                        std::cout << "w= " << w[2][l] << ", y= " << y_p[0] << std::endl;
                    #endif
                }
            }
            domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, addEM_S, addEM_F, domain->NodeIDsLeft[k]);
        }
    }

    if(domain->m_faceOffset[1] > -1)
    {
        if (addEM_S)
            fill(EM_S.begin(), EM_S.end(), zero);
        if (addEM_F) {
            EM_F[0] = zero;
            EM_F[2] = zero;
        }

        for (index_t k=0; k<domain->NodeIDsRight.size()-1; k++) {
            int id = domain->m_faceOffset[1]+k;
            int l = domain->NodeIDsRight[k].level;

            #ifdef OXLEY_ENABLE_DEBUG_ASSEMBLE_PDE_BOUNDARY
            std::cout << "Right\tid: " << id << " , level: " << l << std::endl;
            #endif

            ///////////////
            // process d //
            ///////////////
            if (addEM_S) {
                const Scalar* d_p = d.getSampleDataRO(id, zero);
                if (d.actsExpanded()) {
                    const Scalar d_0 = d_p[0];
                    const Scalar d_1 = d_p[1];
                    const Scalar tmp0 = w[2][l]*(d_0 + d_1);
                    EM_S[INDEX2(1,1,4)] = d_0*w[0][l] + d_1*w[1][l];
                    EM_S[INDEX2(3,1,4)] = tmp0;
                    EM_S[INDEX2(1,3,4)] = tmp0;
                    EM_S[INDEX2(3,3,4)] = d_0*w[1][l] + d_1*w[0][l];
                } else { // constant data
                    EM_S[INDEX2(1,1,4)] = 4.*d_p[0]*w[2][l];
                    EM_S[INDEX2(3,1,4)] = 2.*d_p[0]*w[2][l];
                    EM_S[INDEX2(1,3,4)] = 2.*d_p[0]*w[2][l];
                    EM_S[INDEX2(3,3,4)] = 4.*d_p[0]*w[2][l];
                }
            }

            ///////////////
            // process y //
            ///////////////
            if (addEM_F) {
                const Scalar* y_p = y.getSampleDataRO(id, zero);
                if (y.actsExpanded()) {
                    EM_F[1] = w[3][l]*y_p[0] + w[4][l]*y_p[1];
                    EM_F[3] = w[3][l]*y_p[1] + w[4][l]*y_p[0];
                } else { // constant data
                    EM_F[1] = 6.*w[2][l]*y_p[0];
                    EM_F[3] = 6.*w[2][l]*y_p[0];

                    #ifdef OXLEY_ENABLE_DEBUG_ASSEMBLE_PDE_BOUNDARY
                        std::cout << "w= " << w[2][l] << ", y= " << y_p[0] << std::endl;
                    #endif
                }
            }
            domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, addEM_S, addEM_F, domain->NodeIDsRight[k]);
        }
    }

    if(domain->m_faceOffset[2] > -1)
    {
        if (addEM_S)
            fill(EM_S.begin(), EM_S.end(), zero);
        if (addEM_F) {
            EM_F[2] = zero;
            EM_F[3] = zero;
        }

        for (index_t k=0; k<domain->NodeIDsBottom.size()-1; k++) {
            int id = domain->m_faceOffset[2]+k;
            int l = domain->NodeIDsBottom[k].level;

            #ifdef OXLEY_ENABLE_DEBUG_ASSEMBLE_PDE_BOUNDARY
            std::cout << "Bottom\tid: " << id << " , level: " << l << std::endl;
            #endif

            ///////////////
            // process d //
            ///////////////
            if (addEM_S) {
                const Scalar* d_p = d.getSampleDataRO(id, zero);
                if (d.actsExpanded()) {
                    const Scalar d_0 = d_p[0];
                    const Scalar d_1 = d_p[1];
                    const Scalar tmp0 = w[5][l]*(d_0 + d_1);
                    EM_S[INDEX2(0,0,4)] = d_0*w[6][l] + d_1*w[7][l];
                    EM_S[INDEX2(1,0,4)] = tmp0;
                    EM_S[INDEX2(0,1,4)] = tmp0;
                    EM_S[INDEX2(1,1,4)] = d_0*w[7][l] + d_1*w[6][l];
                } else { // constant data
                    EM_S[INDEX2(0,0,4)] = 4.*d_p[0]*w[5][l];
                    EM_S[INDEX2(1,0,4)] = 2.*d_p[0]*w[5][l];
                    EM_S[INDEX2(0,1,4)] = 2.*d_p[0]*w[5][l];
                    EM_S[INDEX2(1,1,4)] = 4.*d_p[0]*w[5][l];
                }
            }

            ///////////////
            // process y //
            ///////////////
            if (addEM_F) {
                const Scalar* y_p = y.getSampleDataRO(id, zero);
                if (y.actsExpanded()) {
                    EM_F[0] = w[8][l]*y_p[0] + w[9][l]*y_p[1];
                    EM_F[1] = w[8][l]*y_p[1] + w[9][l]*y_p[0];
                } else { // constant data
                    EM_F[0] = 6.*w[5][l]*y_p[0];
                    EM_F[1] = 6.*w[5][l]*y_p[0];

                    #ifdef OXLEY_ENABLE_DEBUG_ASSEMBLE_PDE_BOUNDARY
                        std::cout << "w= " << w[5][l] << ", y= " << y_p[0] << std::endl;
                    #endif
                }
            }
            domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, addEM_S, addEM_F, domain->NodeIDsBottom[k]);
        }
    }

    if(domain->m_faceOffset[3] > -1)
    {

        if (addEM_S)
            fill(EM_S.begin(), EM_S.end(), zero);
        if (addEM_F) {
            EM_F[0] = zero;
            EM_F[1] = zero;
        }

        for (index_t k=0; k<domain->NodeIDsTop.size()-1; k++) {
            int id = domain->m_faceOffset[3]+k;
            int l = domain->NodeIDsTop[k].level;

            #ifdef OXLEY_ENABLE_DEBUG_ASSEMBLE_PDE_BOUNDARY
            std::cout << "Top\tid: " << id << " , level: " << l << std::endl;
            #endif

            ///////////////
            // process d //
            ///////////////
            if (addEM_S) {
                const Scalar* d_p = d.getSampleDataRO(id, zero);
                if (d.actsExpanded()) {
                    const Scalar d_0 = d_p[0];
                    const Scalar d_1 = d_p[1];
                    const Scalar tmp0 = w[5][l]*(d_0 + d_1);
                    EM_S[INDEX2(2,2,4)] = d_0*w[6][l] + d_1*w[7][l];
                    EM_S[INDEX2(3,2,4)] = tmp0;
                    EM_S[INDEX2(2,3,4)] = tmp0;
                    EM_S[INDEX2(3,3,4)] = d_0*w[7][l] + d_1*w[6][l];
                } else { // constant data
                    EM_S[INDEX2(2,2,4)] = 4.*d_p[0]*w[5][l];
                    EM_S[INDEX2(3,2,4)] = 2.*d_p[0]*w[5][l];
                    EM_S[INDEX2(2,3,4)] = 2.*d_p[0]*w[5][l];
                    EM_S[INDEX2(3,3,4)] = 4.*d_p[0]*w[5][l];
                }
            }

            ///////////////
            // process y //
            ///////////////
            if (addEM_F) {
                const Scalar* y_p = y.getSampleDataRO(id, zero);
                if (y.actsExpanded()) {
                    EM_F[2] = w[8][l]*y_p[0] + w[9][l]*y_p[1];
                    EM_F[3] = w[8][l]*y_p[1] + w[9][l]*y_p[0];
                } else { // constant data
                    EM_F[2] = 6.*w[5][l]*y_p[0];
                    EM_F[3] = 6.*w[5][l]*y_p[0];

                    #ifdef OXLEY_ENABLE_DEBUG_ASSEMBLE_PDE_BOUNDARY
                        std::cout << "w= " << w[5][l] << ", y= " << y_p[0] << std::endl;
                    #endif
                }
            }
            domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, addEM_S, addEM_F, domain->NodeIDsTop[k]);
        }
    }    
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
    // Find the maximum level of refinement in the mesh
    int max_level = 0;
    for(p4est_topidx_t tree = domain->p4est->first_local_tree; tree <= domain->p4est->last_local_tree; tree++) {
        p4est_tree_t * tree_t = p4est_tree_array_index(domain->p4est->trees, tree);
        max_level = tree_t->maxlevel > max_level ? tree_t->maxlevel : max_level;
    }

    // This is similar to Ripley, except that in Oxley the quads vary in size so that the 
    // weightings need to be calculated for each level of refinement
    double w[8][P4EST_MAXLEVEL] = {{0}};
// #pragma omp parallel for
    for(int i = 0; i <= max_level; i++)
    {
        double m_dx[2] = {domain->forestData.m_dx[0][P4EST_MAXLEVEL-i], 
                          domain->forestData.m_dx[1][P4EST_MAXLEVEL-i]};

        w[0][i] = 1./4;
        w[1][i] = m_dx[0]/8;
        w[2][i] = m_dx[1]/8;
        w[3][i] = m_dx[0]*m_dx[1]/16;
        w[4][i] = m_dx[0]/(4*m_dx[1]);
        w[5][i] = m_dx[1]/(4*m_dx[0]);
    }

    const bool addEM_S = (!A.isEmpty() || !B.isEmpty() || !C.isEmpty() || !D.isEmpty());
    const bool addEM_F = (!X.isEmpty() || !Y.isEmpty());
    const Scalar zero = static_cast<Scalar>(0);
    rhs.requireWrite();

    vector<Scalar> EM_S(4*4, zero);
    vector<Scalar> EM_F(4, zero);

// #pragma omp parallel for
    for (p4est_topidx_t t = domain->p4est->first_local_tree; t <= domain->p4est->last_local_tree; t++) // Loop over every tree
    {
        p4est_tree_t * currenttree = p4est_tree_array_index(domain->p4est->trees, t);
        sc_array_t * tquadrants = &currenttree->quadrants;
        p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
        for (int q = 0; q < Q; ++q)  // Loop over the elements attached to the tree
        {                
            p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, q);
            int l = quad->level;
            double xy[3];
            p4est_qcoord_to_vertex(domain->p4est->connectivity, t, quad->x, quad->y, xy);
            long id = domain->getQuadID(domain->NodeIDs.find(std::make_pair(xy[0],xy[1]))->second);

            if (addEM_S)
                fill(EM_S.begin(), EM_S.end(), zero);
            if (addEM_F)
                fill(EM_F.begin(), EM_F.end(), zero);

            ///////////////
            // process A //
            ///////////////
            if (!A.isEmpty()) {
                const Scalar* A_p = A.getSampleDataRO(id, zero);
                const Scalar A_00 = A_p[INDEX2(0,0,2)];
                const Scalar A_10 = A_p[INDEX2(1,0,2)];
                const Scalar A_01 = A_p[INDEX2(0,1,2)];
                const Scalar A_11 = A_p[INDEX2(1,1,2)];
                const Scalar tmp0 = (A_01 + A_10)*w[0][l];
                const Scalar tmp1 = A_00*w[5][l];
                const Scalar tmp2 = A_01*w[0][l];
                const Scalar tmp3 = A_10*w[0][l];
                const Scalar tmp4 = A_11*w[4][l];
                EM_S[INDEX2(0,0,4)]+= tmp4 + tmp0 + tmp1;
                EM_S[INDEX2(1,0,4)]+= tmp4 - tmp1 + tmp3 - tmp2;
                EM_S[INDEX2(2,0,4)]+= tmp2 - tmp3 - tmp4 + tmp1;
                EM_S[INDEX2(3,0,4)]+=-tmp1 - tmp4 - tmp0;
                EM_S[INDEX2(0,1,4)]+= tmp4 - tmp1 + tmp2 - tmp3;
                EM_S[INDEX2(1,1,4)]+= tmp4 + tmp1 - tmp0;
                EM_S[INDEX2(2,1,4)]+=-tmp1 + tmp0 - tmp4;
                EM_S[INDEX2(3,1,4)]+=-tmp4 + tmp1 + tmp3 - tmp2;
                EM_S[INDEX2(0,2,4)]+=-tmp4 + tmp1 + tmp3 - tmp2;
                EM_S[INDEX2(1,2,4)]+=-tmp1 + tmp0 - tmp4;
                EM_S[INDEX2(2,2,4)]+= tmp4 + tmp1 - tmp0;
                EM_S[INDEX2(3,2,4)]+= tmp4 - tmp1 + tmp2 - tmp3;
                EM_S[INDEX2(0,3,4)]+=-tmp1 - tmp4 - tmp0;
                EM_S[INDEX2(1,3,4)]+= tmp2 - tmp3 - tmp4 + tmp1;
                EM_S[INDEX2(2,3,4)]+= tmp4 - tmp1 + tmp3 - tmp2;
                EM_S[INDEX2(3,3,4)]+= tmp4 + tmp0 + tmp1;
            }

            ///////////////
            // process B //
            ///////////////
            if (!B.isEmpty()) {
                const Scalar* B_p = B.getSampleDataRO(id, zero);
                const Scalar tmp0 = B_p[0]*w[2][l];
                const Scalar tmp1 = B_p[1]*w[1][l];
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
                const Scalar* C_p = C.getSampleDataRO(id, zero);
                const Scalar tmp0 = C_p[0]*w[2][l];
                const Scalar tmp1 = C_p[1]*w[1][l];
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
                const Scalar* D_p = D.getSampleDataRO(id, zero);
                EM_S[INDEX2(0,0,4)]+=D_p[0]*w[3][l];
                EM_S[INDEX2(1,0,4)]+=D_p[0]*w[3][l];
                EM_S[INDEX2(2,0,4)]+=D_p[0]*w[3][l];
                EM_S[INDEX2(3,0,4)]+=D_p[0]*w[3][l];
                EM_S[INDEX2(0,1,4)]+=D_p[0]*w[3][l];
                EM_S[INDEX2(1,1,4)]+=D_p[0]*w[3][l];
                EM_S[INDEX2(2,1,4)]+=D_p[0]*w[3][l];
                EM_S[INDEX2(3,1,4)]+=D_p[0]*w[3][l];
                EM_S[INDEX2(0,2,4)]+=D_p[0]*w[3][l];
                EM_S[INDEX2(1,2,4)]+=D_p[0]*w[3][l];
                EM_S[INDEX2(2,2,4)]+=D_p[0]*w[3][l];
                EM_S[INDEX2(3,2,4)]+=D_p[0]*w[3][l];
                EM_S[INDEX2(0,3,4)]+=D_p[0]*w[3][l];
                EM_S[INDEX2(1,3,4)]+=D_p[0]*w[3][l];
                EM_S[INDEX2(2,3,4)]+=D_p[0]*w[3][l];
                EM_S[INDEX2(3,3,4)]+=D_p[0]*w[3][l];
            }

            ///////////////
            // process X //
            ///////////////
            if (!X.isEmpty()) {
                const Scalar* X_p = X.getSampleDataRO(id, zero);
                const Scalar wX0 = 4.*X_p[0]*w[2][l];
                const Scalar wX1 = 4.*X_p[1]*w[1][l];
                EM_F[0]+=-wX0 - wX1;
                EM_F[1]+=-wX1 + wX0;
                EM_F[2]+=-wX0 + wX1;
                EM_F[3]+= wX0 + wX1;
            }

            ///////////////
            // process Y //
            ///////////////
            if (!Y.isEmpty()) {
                const Scalar* Y_p = Y.getSampleDataRO(id, zero);
                EM_F[0]+=4.*Y_p[0]*w[3][l];
                EM_F[1]+=4.*Y_p[0]*w[3][l];
                EM_F[2]+=4.*Y_p[0]*w[3][l];
                EM_F[3]+=4.*Y_p[0]*w[3][l];
            }

            // add to matrix (if addEM_S) and RHS (if addEM_F)
            domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, addEM_S, addEM_F, q, t);
        } 
    }
}

/****************************************************************************/
// PDE SINGLE REDUCED BOUNDARY
/****************************************************************************/

template<class Scalar>
void DefaultAssembler2D<Scalar>::assemblePDEBoundarySingleReduced(
                                        AbstractSystemMatrix* mat, Data& rhs,
                                        const Data& d, const Data& y) const
{
    // Find the maximum level of refinement in the mesh
    int max_level = 0;
    for(p4est_topidx_t tree = domain->p4est->first_local_tree; tree <= domain->p4est->last_local_tree; tree++) {
        p4est_tree_t * tree_t = p4est_tree_array_index(domain->p4est->trees, tree);
        max_level = tree_t->maxlevel > max_level ? tree_t->maxlevel : max_level;
    }

    // This is similar to Ripley, except that in Oxley the quads vary in size so that the 
    // weightings need to be calcuated for each level of refinement
    double w[2][P4EST_MAXLEVEL] = {{0}};
// #pragma omp parallel for
    for(int i = 0; i <= max_level; i++)
    {
        double m_dx[2] = {domain->forestData.m_dx[0][P4EST_MAXLEVEL-i], 
                          domain->forestData.m_dx[1][P4EST_MAXLEVEL-i]};

        w[0][i] = m_dx[0]/4;
        w[1][i] = m_dx[1]/4;
    }

    const bool addEM_S = !d.isEmpty();
    const bool addEM_F = !y.isEmpty();
    const Scalar zero = static_cast<Scalar>(0);
    rhs.requireWrite();

    vector<Scalar> EM_S(4*4, zero);
    vector<Scalar> EM_F(4, zero);

    if(domain->m_faceOffset[0]>-1)
    {
        if (addEM_S)
            fill(EM_S.begin(), EM_S.end(), zero);
        if (addEM_F) {
            EM_F[1] = zero;
            EM_F[3] = zero;
        }

        for (index_t k=0; k<domain->NodeIDsLeft.size()-1; k++) {
            int id = k;
            int l = domain->NodeIDsLeft[k].level;

            ///////////////
            // process d //
            ///////////////
            if (addEM_S) {
                const Scalar* d_p = d.getSampleDataRO(id, zero);
                EM_S[INDEX2(0,0,4)] = d_p[0]*w[1][l];
                EM_S[INDEX2(2,0,4)] = d_p[0]*w[1][l];
                EM_S[INDEX2(0,2,4)] = d_p[0]*w[1][l];
                EM_S[INDEX2(2,2,4)] = d_p[0]*w[1][l];
            }

            ///////////////
            // process y //
            ///////////////
            if (addEM_F) {
                const Scalar* y_p = y.getSampleDataRO(id, zero);
                EM_F[0] = 2.*w[1][l]*y_p[0];
                EM_F[2] = 2.*w[1][l]*y_p[0];
            }
            domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, addEM_S, addEM_F, domain->NodeIDsLeft[k]);
        } 
    }

    if(domain->m_faceOffset[1]>-1)
    {
        if (addEM_S)
            fill(EM_S.begin(), EM_S.end(), zero);
        if (addEM_F) {
            EM_F[0] = zero;
            EM_F[2] = zero;
        }

        for (index_t k=0; k<domain->NodeIDsRight.size()-1; k++) {
            int id = domain->m_faceOffset[1]+k;
            int l = domain->NodeIDsRight[k].level;

            ///////////////
            // process d //
            ///////////////
            if (addEM_S) {
                const Scalar* d_p = d.getSampleDataRO(id, zero);
                EM_S[INDEX2(1,1,4)] = d_p[0]*w[1][l];
                EM_S[INDEX2(3,1,4)] = d_p[0]*w[1][l];
                EM_S[INDEX2(1,3,4)] = d_p[0]*w[1][l];
                EM_S[INDEX2(3,3,4)] = d_p[0]*w[1][l];
            }

            ///////////////
            // process y //
            ///////////////
            if (addEM_F) {
                const Scalar* y_p = y.getSampleDataRO(id, zero);
                EM_F[1] = 2.*w[1][l]*y_p[0];
                EM_F[3] = 2.*w[1][l]*y_p[0];
            }
            domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, addEM_S, addEM_F, domain->NodeIDsRight[k]);
        } 
    }

    if(domain->m_faceOffset[2]>-1)
    {
        if (addEM_S)
            fill(EM_S.begin(), EM_S.end(), zero);
        if (addEM_F) {
            EM_F[2] = zero;
            EM_F[3] = zero;
        }

        for (index_t k=0; k<domain->NodeIDsBottom.size()-1; k++) {
            int id = domain->m_faceOffset[2]+k;
            int l = domain->NodeIDsBottom[k].level;

            ///////////////
            // process d //
            ///////////////
            if (addEM_S) {
                const Scalar* d_p = d.getSampleDataRO(id, zero);
                EM_S[INDEX2(0,0,4)] = d_p[0]*w[0][l];
                EM_S[INDEX2(1,0,4)] = d_p[0]*w[0][l];
                EM_S[INDEX2(0,1,4)] = d_p[0]*w[0][l];
                EM_S[INDEX2(1,1,4)] = d_p[0]*w[0][l];
            }

            ///////////////
            // process y //
            ///////////////
            if (addEM_F) {
                const Scalar* y_p = y.getSampleDataRO(id, zero);
                EM_F[0] = 2.*w[0][l]*y_p[0];
                EM_F[1] = 2.*w[0][l]*y_p[0];
            }
            domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, addEM_S, addEM_F, domain->NodeIDsBottom[k]);
        }
    }

    if(domain->m_faceOffset[3]>-1)
    {
        if (addEM_S)
            fill(EM_S.begin(), EM_S.end(), zero);
        if (addEM_F) {
            EM_F[0] = zero;
            EM_F[1] = zero;
        }
        
        for (index_t k=0; k<domain->NodeIDsTop.size()-1; k++) {
            int id = domain->m_faceOffset[3]+k;
            int l = domain->NodeIDsTop[k].level;

            ///////////////
            // process d //
            ///////////////
            if (addEM_S) {
                const Scalar* d_p = d.getSampleDataRO(id, zero);
                EM_S[INDEX2(2,2,4)] = d_p[0]*w[0][l];
                EM_S[INDEX2(3,2,4)] = d_p[0]*w[0][l];
                EM_S[INDEX2(2,3,4)] = d_p[0]*w[0][l];
                EM_S[INDEX2(3,3,4)] = d_p[0]*w[0][l];
            }

            ///////////////
            // process y //
            ///////////////
            if (addEM_F) {
                const Scalar* y_p = y.getSampleDataRO(id, zero);
                EM_F[2] = 2.*w[0][l]*y_p[0];
                EM_F[3] = 2.*w[0][l]*y_p[0];
            }
            domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, addEM_S, addEM_F, domain->NodeIDsTop[k]);
        } 
    }
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

    // Find the maximum level of refinement in the mesh
    int max_level = 0;
    for(p4est_topidx_t tree = domain->p4est->first_local_tree; tree <= domain->p4est->last_local_tree; tree++) {
        p4est_tree_t * tree_t = p4est_tree_array_index(domain->p4est->trees, tree);
        max_level = tree_t->maxlevel > max_level ? tree_t->maxlevel : max_level;
    }

    // This is similar to Ripley, except that in Oxley the quads vary in size so that the 
    // weightings need to be calcuated for each level of refinement
    const double SQRT3 = 1.73205080756887719318;
    double w[32][P4EST_MAXLEVEL] = {{0}};
// #pragma omp parallel for
    for(int i = 0; i <= max_level; i++)
    {
        double m_dx[2] = {domain->forestData.m_dx[0][P4EST_MAXLEVEL-i], 
                          domain->forestData.m_dx[1][P4EST_MAXLEVEL-i]};

        w[1][i]  = 1.0/24;
        w[5][i]  = -SQRT3/24 + 1.0/12;
        w[2][i]  = -SQRT3/24 - 1.0/12;
        w[19][i] = -m_dx[0]/12;
        w[11][i] = w[19][i]*(SQRT3 + 3)/12;
        w[14][i] = w[19][i]*(-SQRT3 + 3)/12;
        w[16][i] = w[19][i]*(5*SQRT3 + 9)/12;
        w[17][i] = w[19][i]*(-5*SQRT3 + 9)/12;
        w[27][i] = w[19][i]*(-SQRT3 - 3)/2;
        w[28][i] = w[19][i]*(SQRT3 - 3)/2;
        w[18][i] = -m_dx[1]/12;
        w[10][i] = w[18][i]*(SQRT3 + 3)/12;
        w[15][i] = w[18][i]*(-SQRT3 + 3)/12;
        w[12][i] = w[18][i]*(5*SQRT3 + 9)/12;
        w[13][i] = w[18][i]*(-5*SQRT3 + 9)/12;
        w[25][i] = w[18][i]*(-SQRT3 - 3)/2;
        w[26][i] = w[18][i]*(SQRT3 - 3)/2;
        w[22][i] = m_dx[0]*m_dx[1]/144;
        w[20][i] = w[22][i]*(SQRT3 + 2);
        w[21][i] = w[22][i]*(-SQRT3 + 2);
        w[23][i] = w[22][i]*(4*SQRT3 + 7);
        w[24][i] = w[22][i]*(-4*SQRT3 + 7);
        w[3][i]  = m_dx[0]/(24*m_dx[1]);
        w[7][i]  = w[3][i]*(SQRT3 + 2);
        w[8][i]  = w[3][i]*(-SQRT3 + 2);
        w[6][i]  = -m_dx[1]/(24*m_dx[0]);
        w[0][i]  = w[6][i]*(SQRT3 + 2);
        w[4][i]  = w[6][i]*(-SQRT3 + 2);
    }

#ifdef OXLEY_ENABLE_DEBUG
    std::cout << "w" << std::endl;
    for(int i = 0; i < 29; i++)
        std::cout << i << ": " << w[i][0] << std::endl;
    std::cout << std::endl;
#endif


    const bool addEM_S = (!A.isEmpty() || !B.isEmpty() || !C.isEmpty() || !D.isEmpty());
    const bool addEM_F = (!X.isEmpty() || !Y.isEmpty());
    const Scalar zero = static_cast<Scalar>(0);
    rhs.requireWrite();

    vector<Scalar> EM_S(4*4*numEq*numComp, zero);
    vector<Scalar> EM_F(4*numEq, zero);

// // #pragma omp parallel for
    for (p4est_topidx_t t = domain->p4est->first_local_tree; t <= domain->p4est->last_local_tree; t++) // Loop over every tree
    {
        p4est_tree_t * currenttree = p4est_tree_array_index(domain->p4est->trees, t);
        sc_array_t * tquadrants = &currenttree->quadrants;
        p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
        for (int q = 0; q < Q; ++q)  // Loop over the elements attached to the tree
        {                
            p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, q);
            int l = quad->level;
            double xy[3];
            p4est_qcoord_to_vertex(domain->p4est->connectivity, t, quad->x, quad->y, xy);
            long id = domain->getQuadID(domain->NodeIDs.find(std::make_pair(xy[0],xy[1]))->second);

            if (addEM_S)
                fill(EM_S.begin(), EM_S.end(), zero);
            if (addEM_F)
                fill(EM_F.begin(), EM_F.end(), zero);

            ///////////////
            // process A //
            ///////////////
            if (!A.isEmpty()) {
                const Scalar* A_p = A.getSampleDataRO(id, zero);
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
                            const Scalar tmp0  = w[3][l]*(A_11_0 + A_11_1 + A_11_2 + A_11_3);
                            const Scalar tmp1  = w[1][l]*(A_01_0 + A_01_3 - A_10_1 - A_10_2);
                            const Scalar tmp2  = w[4][l]*(A_00_2 + A_00_3);
                            const Scalar tmp3  = w[0][l]*(A_00_0 + A_00_1);
                            const Scalar tmp4  = w[5][l]*(A_01_2 - A_10_3);
                            const Scalar tmp5  = w[2][l]*(-A_01_1 + A_10_0);
                            const Scalar tmp6  = w[5][l]*(A_01_3 + A_10_0);
                            const Scalar tmp7  = w[3][l]*(-A_11_0 - A_11_1 - A_11_2 - A_11_3);
                            const Scalar tmp8  = w[6][l]*(A_00_0 + A_00_1 + A_00_2 + A_00_3);
                            const Scalar tmp9  = w[1][l]*(A_01_1 + A_01_2 + A_10_1 + A_10_2);
                            const Scalar tmp10 = w[2][l]*(-A_01_0 - A_10_3);
                            const Scalar tmp11 = w[4][l]*(A_00_0 + A_00_1);
                            const Scalar tmp12 = w[0][l]*(A_00_2 + A_00_3);
                            const Scalar tmp13 = w[5][l]*(A_01_1 - A_10_0);
                            const Scalar tmp14 = w[2][l]*(-A_01_2 + A_10_3);
                            const Scalar tmp15 = w[7][l]*(A_11_0 + A_11_2);
                            const Scalar tmp16 = w[4][l]*(-A_00_2 - A_00_3);
                            const Scalar tmp17 = w[0][l]*(-A_00_0 - A_00_1);
                            const Scalar tmp18 = w[5][l]*(A_01_3 + A_10_3);
                            const Scalar tmp19 = w[8][l]*(A_11_1 + A_11_3);
                            const Scalar tmp20 = w[2][l]*(-A_01_0 - A_10_0);
                            const Scalar tmp21 = w[7][l]*(A_11_1 + A_11_3);
                            const Scalar tmp22 = w[4][l]*(-A_00_0 - A_00_1);
                            const Scalar tmp23 = w[0][l]*(-A_00_2 - A_00_3);
                            const Scalar tmp24 = w[5][l]*(A_01_0 + A_10_0);
                            const Scalar tmp25 = w[8][l]*(A_11_0 + A_11_2);
                            const Scalar tmp26 = w[2][l]*(-A_01_3 - A_10_3);
                            const Scalar tmp27 = w[5][l]*(-A_01_1 - A_10_2);
                            const Scalar tmp28 = w[1][l]*(-A_01_0 - A_01_3 - A_10_0 - A_10_3);
                            const Scalar tmp29 = w[2][l]*(A_01_2 + A_10_1);
                            const Scalar tmp30 = w[7][l]*(-A_11_1 - A_11_3);
                            const Scalar tmp31 = w[1][l]*(-A_01_1 - A_01_2 + A_10_0 + A_10_3);
                            const Scalar tmp32 = w[5][l]*(-A_01_0 + A_10_2);
                            const Scalar tmp33 = w[8][l]*(-A_11_0 - A_11_2);
                            const Scalar tmp34 = w[6][l]*(-A_00_0 - A_00_1 - A_00_2 - A_00_3);
                            const Scalar tmp35 = w[2][l]*(A_01_3 - A_10_1);
                            const Scalar tmp36 = w[5][l]*(A_01_0 + A_10_3);
                            const Scalar tmp37 = w[2][l]*(-A_01_3 - A_10_0);
                            const Scalar tmp38 = w[7][l]*(-A_11_0 - A_11_2);
                            const Scalar tmp39 = w[5][l]*(-A_01_3 + A_10_1);
                            const Scalar tmp40 = w[8][l]*(-A_11_1 - A_11_3);
                            const Scalar tmp41 = w[2][l]*(A_01_0 - A_10_2);
                            const Scalar tmp42 = w[5][l]*(A_01_1 - A_10_3);
                            const Scalar tmp43 = w[2][l]*(-A_01_2 + A_10_0);
                            const Scalar tmp44 = w[5][l]*(A_01_2 - A_10_0);
                            const Scalar tmp45 = w[2][l]*(-A_01_1 + A_10_3);
                            const Scalar tmp46 = w[5][l]*(-A_01_0 + A_10_1);
                            const Scalar tmp47 = w[2][l]*(A_01_3 - A_10_2);
                            const Scalar tmp48 = w[5][l]*(-A_01_1 - A_10_1);
                            const Scalar tmp49 = w[2][l]*(A_01_2 + A_10_2);
                            const Scalar tmp50 = w[5][l]*(-A_01_3 + A_10_2);
                            const Scalar tmp51 = w[2][l]*(A_01_0 - A_10_1);
                            const Scalar tmp52 = w[5][l]*(-A_01_2 - A_10_1);
                            const Scalar tmp53 = w[2][l]*(A_01_1 + A_10_2);
                            const Scalar tmp54 = w[5][l]*(-A_01_2 - A_10_2);
                            const Scalar tmp55 = w[2][l]*(A_01_1 + A_10_1);
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
                            const Scalar tmp0 = 6.*w[1][l]*(A_01 - A_10);
                            const Scalar tmp1 = 6.*w[1][l]*(A_01 + A_10);
                            const Scalar tmp2 = 6.*w[1][l]*(-A_01 - A_10);
                            const Scalar tmp3 = 6.*w[1][l]*(-A_01 + A_10);
                            EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=-8.*A_00*w[6][l] + 8.*A_11*w[3][l] + tmp1;
                            EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+= 8.*A_00*w[6][l] + 4.*A_11*w[3][l] + tmp0;
                            EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=-4.*A_00*w[6][l] - 8.*A_11*w[3][l] + tmp3;
                            EM_S[INDEX4(k,m,0,3,numEq,numComp,4)]+= 4.*A_00*w[6][l] - 4.*A_11*w[3][l] + tmp2;
                            EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+= 8.*A_00*w[6][l] + 4.*A_11*w[3][l] + tmp3;
                            EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=-8.*A_00*w[6][l] + 8.*A_11*w[3][l] + tmp2;
                            EM_S[INDEX4(k,m,1,2,numEq,numComp,4)]+= 4.*A_00*w[6][l] - 4.*A_11*w[3][l] + tmp1;
                            EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+=-4.*A_00*w[6][l] - 8.*A_11*w[3][l] + tmp0;
                            EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=-4.*A_00*w[6][l] - 8.*A_11*w[3][l] + tmp0;
                            EM_S[INDEX4(k,m,2,1,numEq,numComp,4)]+= 4.*A_00*w[6][l] - 4.*A_11*w[3][l] + tmp1;
                            EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=-8.*A_00*w[6][l] + 8.*A_11*w[3][l] + tmp2;
                            EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+= 8.*A_00*w[6][l] + 4.*A_11*w[3][l] + tmp3;
                            EM_S[INDEX4(k,m,3,0,numEq,numComp,4)]+= 4.*A_00*w[6][l] - 4.*A_11*w[3][l] + tmp2;
                            EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+=-4.*A_00*w[6][l] - 8.*A_11*w[3][l] + tmp3;
                            EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+= 8.*A_00*w[6][l] + 4.*A_11*w[3][l] + tmp0;
                            EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=-8.*A_00*w[6][l] + 8.*A_11*w[3][l] + tmp1;
                        }
                    }
                }
            }

            ///////////////
            // process B //
            ///////////////
            if (!B.isEmpty()) {
                const Scalar* B_p = B.getSampleDataRO(id, zero);
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
                            const Scalar tmp0  = w[11][l]*(B_1_0 + B_1_1);
                            const Scalar tmp1  = w[14][l]*(B_1_2 + B_1_3);
                            const Scalar tmp2  = w[15][l]*(-B_0_1 - B_0_3);
                            const Scalar tmp3  = w[10][l]*(-B_0_0 - B_0_2);
                            const Scalar tmp4  = w[11][l]*(B_1_2 + B_1_3);
                            const Scalar tmp5  = w[14][l]*(B_1_0 + B_1_1);
                            const Scalar tmp6  = w[11][l]*(-B_1_2 - B_1_3);
                            const Scalar tmp7  = w[14][l]*(-B_1_0 - B_1_1);
                            const Scalar tmp8  = w[11][l]*(-B_1_0 - B_1_1);
                            const Scalar tmp9  = w[14][l]*(-B_1_2 - B_1_3);
                            const Scalar tmp10 = w[10][l]*(-B_0_1 - B_0_3);
                            const Scalar tmp11 = w[15][l]*(-B_0_0 - B_0_2);
                            const Scalar tmp12 = w[15][l]*(B_0_0 + B_0_2);
                            const Scalar tmp13 = w[10][l]*(B_0_1 + B_0_3);
                            const Scalar tmp14 = w[10][l]*(B_0_0 + B_0_2);
                            const Scalar tmp15 = w[15][l]*(B_0_1 + B_0_3);
                            EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+= B_0_0*w[12][l] + B_0_1*w[10][l] + B_0_2*w[15][l] + B_0_3*w[13][l] + B_1_0*w[16][l] + B_1_1*w[14][l] + B_1_2*w[11][l] + B_1_3*w[17][l];
                            EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+= B_0_0*w[10][l] + B_0_1*w[12][l] + B_0_2*w[13][l] + B_0_3*w[15][l] + tmp0 + tmp1;
                            EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+= B_1_0*w[11][l] + B_1_1*w[17][l] + B_1_2*w[16][l] + B_1_3*w[14][l] + tmp14 + tmp15;
                            EM_S[INDEX4(k,m,0,3,numEq,numComp,4)]+=tmp12 + tmp13 + tmp4 + tmp5;
                            EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+=-B_0_0*w[12][l] - B_0_1*w[10][l] - B_0_2*w[15][l] - B_0_3*w[13][l] + tmp0 + tmp1;
                            EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=-B_0_0*w[10][l] - B_0_1*w[12][l] - B_0_2*w[13][l] - B_0_3*w[15][l] + B_1_0*w[14][l] + B_1_1*w[16][l] + B_1_2*w[17][l] + B_1_3*w[11][l];
                            EM_S[INDEX4(k,m,1,2,numEq,numComp,4)]+=tmp2 + tmp3 + tmp4 + tmp5;
                            EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+= B_1_0*w[17][l] + B_1_1*w[11][l] + B_1_2*w[14][l] + B_1_3*w[16][l] + tmp10 + tmp11;
                            EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=-B_1_0*w[16][l] - B_1_1*w[14][l] - B_1_2*w[11][l] - B_1_3*w[17][l] + tmp14 + tmp15;
                            EM_S[INDEX4(k,m,2,1,numEq,numComp,4)]+=tmp12 + tmp13 + tmp8 + tmp9;
                            EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+= B_0_0*w[15][l] + B_0_1*w[13][l] + B_0_2*w[12][l] + B_0_3*w[10][l] - B_1_0*w[11][l] - B_1_1*w[17][l] - B_1_2*w[16][l] - B_1_3*w[14][l];
                            EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+= B_0_0*w[13][l] + B_0_1*w[15][l] + B_0_2*w[10][l] + B_0_3*w[12][l] + tmp6 + tmp7;
                            EM_S[INDEX4(k,m,3,0,numEq,numComp,4)]+=tmp2 + tmp3 + tmp8 + tmp9;
                            EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+=-B_1_0*w[14][l] - B_1_1*w[16][l] - B_1_2*w[17][l] - B_1_3*w[11][l] + tmp10 + tmp11;
                            EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+=-B_0_0*w[15][l] - B_0_1*w[13][l] - B_0_2*w[12][l] - B_0_3*w[10][l] + tmp6 + tmp7;
                            EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=-B_0_0*w[13][l] - B_0_1*w[15][l] - B_0_2*w[10][l] - B_0_3*w[12][l] - B_1_0*w[17][l] - B_1_1*w[11][l] - B_1_2*w[14][l] - B_1_3*w[16][l];
                        }
                    }
                } else { // constant data
                    for (index_t k=0; k<numEq; k++) {
                        for (index_t m=0; m<numComp; m++) {
                            const Scalar wB0 = B_p[INDEX3(k,0,m,numEq,2)]*w[18][l];
                            const Scalar wB1 = B_p[INDEX3(k,1,m,numEq,2)]*w[19][l];
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
                const Scalar* C_p = C.getSampleDataRO(id, zero);
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
                            const Scalar tmp0  = w[11][l]*(C_1_0 + C_1_1);
                            const Scalar tmp1  = w[14][l]*(C_1_2 + C_1_3);
                            const Scalar tmp2  = w[15][l]*(C_0_0 + C_0_2);
                            const Scalar tmp3  = w[10][l]*(C_0_1 + C_0_3);
                            const Scalar tmp4  = w[11][l]*(-C_1_0 - C_1_1);
                            const Scalar tmp5  = w[14][l]*(-C_1_2 - C_1_3);
                            const Scalar tmp6  = w[11][l]*(-C_1_2 - C_1_3);
                            const Scalar tmp7  = w[14][l]*(-C_1_0 - C_1_1);
                            const Scalar tmp8  = w[11][l]*(C_1_2 + C_1_3);
                            const Scalar tmp9  = w[14][l]*(C_1_0 + C_1_1);
                            const Scalar tmp10 = w[10][l]*(-C_0_1 - C_0_3);
                            const Scalar tmp11 = w[15][l]*(-C_0_0 - C_0_2);
                            const Scalar tmp12 = w[15][l]*(-C_0_1 - C_0_3);
                            const Scalar tmp13 = w[10][l]*(-C_0_0 - C_0_2);
                            const Scalar tmp14 = w[10][l]*(C_0_0 + C_0_2);
                            const Scalar tmp15 = w[15][l]*(C_0_1 + C_0_3);
                            EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+= C_0_0*w[12][l] + C_0_1*w[10][l] + C_0_2*w[15][l] + C_0_3*w[13][l] + C_1_0*w[16][l] + C_1_1*w[14][l] + C_1_2*w[11][l] + C_1_3*w[17][l];
                            EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+=-C_0_0*w[12][l] - C_0_1*w[10][l] - C_0_2*w[15][l] - C_0_3*w[13][l] + tmp0 + tmp1;
                            EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=-C_1_0*w[16][l] - C_1_1*w[14][l] - C_1_2*w[11][l] - C_1_3*w[17][l] + tmp14 + tmp15;
                            EM_S[INDEX4(k,m,0,3,numEq,numComp,4)]+=tmp12 + tmp13 + tmp4 + tmp5;
                            EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+= C_0_0*w[10][l] + C_0_1*w[12][l] + C_0_2*w[13][l] + C_0_3*w[15][l] + tmp0 + tmp1;
                            EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=-C_0_0*w[10][l] - C_0_1*w[12][l] - C_0_2*w[13][l] - C_0_3*w[15][l] + C_1_0*w[14][l] + C_1_1*w[16][l] + C_1_2*w[17][l] + C_1_3*w[11][l];
                            EM_S[INDEX4(k,m,1,2,numEq,numComp,4)]+=tmp2 + tmp3 + tmp4 + tmp5;
                            EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+=-C_1_0*w[14][l] - C_1_1*w[16][l] - C_1_2*w[17][l] - C_1_3*w[11][l] + tmp10 + tmp11;
                            EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+= C_1_0*w[11][l] + C_1_1*w[17][l] + C_1_2*w[16][l] + C_1_3*w[14][l] + tmp14 + tmp15;
                            EM_S[INDEX4(k,m,2,1,numEq,numComp,4)]+=tmp12 + tmp13 + tmp8 + tmp9;
                            EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+= C_0_0*w[15][l] + C_0_1*w[13][l] + C_0_2*w[12][l] + C_0_3*w[10][l] - C_1_0*w[11][l] - C_1_1*w[17][l] - C_1_2*w[16][l] - C_1_3*w[14][l];
                            EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+=-C_0_0*w[15][l] - C_0_1*w[13][l] - C_0_2*w[12][l] - C_0_3*w[10][l] + tmp6 + tmp7;
                            EM_S[INDEX4(k,m,3,0,numEq,numComp,4)]+=tmp2 + tmp3 + tmp8 + tmp9;
                            EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+= C_1_0*w[17][l] + C_1_1*w[11][l] + C_1_2*w[14][l] + C_1_3*w[16][l] + tmp10 + tmp11;
                            EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+= C_0_0*w[13][l] + C_0_1*w[15][l] + C_0_2*w[10][l] + C_0_3*w[12][l] + tmp6 + tmp7;
                            EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=-C_0_0*w[13][l] - C_0_1*w[15][l] - C_0_2*w[10][l] - C_0_3*w[12][l] - C_1_0*w[17][l] - C_1_1*w[11][l] - C_1_2*w[14][l] - C_1_3*w[16][l];
                        }
                    }
                } else { // constant data
                    for (index_t k=0; k<numEq; k++) {
                        for (index_t m=0; m<numComp; m++) {
                            const Scalar wC0 = C_p[INDEX3(k,m,0,numEq,numComp)]*w[18][l];
                            const Scalar wC1 = C_p[INDEX3(k,m,1,numEq,numComp)]*w[19][l];
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
                const Scalar* D_p = D.getSampleDataRO(id, zero);
                if (D.actsExpanded()) {
                    for (index_t k=0; k<numEq; k++) {
                        for (index_t m=0; m<numComp; m++) {
                            const Scalar D_0 = D_p[INDEX3(k,m,0,numEq,numComp)];
                            const Scalar D_1 = D_p[INDEX3(k,m,1,numEq,numComp)];
                            const Scalar D_2 = D_p[INDEX3(k,m,2,numEq,numComp)];
                            const Scalar D_3 = D_p[INDEX3(k,m,3,numEq,numComp)];
                            const Scalar tmp0  = w[21][l]*(D_2 + D_3);
                            const Scalar tmp1  = w[20][l]*(D_0 + D_1);
                            const Scalar tmp2  = w[22][l]*(D_0 + D_1 + D_2 + D_3);
                            const Scalar tmp3  = w[21][l]*(D_0 + D_1);
                            const Scalar tmp4  = w[20][l]*(D_2 + D_3);
                            const Scalar tmp5  = w[22][l]*(D_1 + D_2);
                            const Scalar tmp6  = w[21][l]*(D_0 + D_2);
                            const Scalar tmp7  = w[20][l]*(D_1 + D_3);
                            const Scalar tmp8  = w[21][l]*(D_1 + D_3);
                            const Scalar tmp9  = w[20][l]*(D_0 + D_2);
                            const Scalar tmp10 = w[22][l]*(D_0 + D_3);
                            EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=D_0*w[23][l] + D_3*w[24][l] + tmp5;
                            EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+=tmp0 + tmp1;
                            EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+=tmp8 + tmp9;
                            EM_S[INDEX4(k,m,0,3,numEq,numComp,4)]+=tmp2;
                            EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+=tmp0 + tmp1;
                            EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=D_1*w[23][l] + D_2*w[24][l] + tmp10;
                            EM_S[INDEX4(k,m,1,2,numEq,numComp,4)]+=tmp2;
                            EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+=tmp6 + tmp7;
                            EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+=tmp8 + tmp9;
                            EM_S[INDEX4(k,m,2,1,numEq,numComp,4)]+=tmp2;
                            EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=D_1*w[24][l] + D_2*w[23][l] + tmp10;
                            EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+=tmp3 + tmp4;
                            EM_S[INDEX4(k,m,3,0,numEq,numComp,4)]+=tmp2;
                            EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+=tmp6 + tmp7;
                            EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+=tmp3 + tmp4;
                            EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=D_0*w[24][l] + D_3*w[23][l] + tmp5;
                        }
                     }
                } else { // constant data
                    for (index_t k=0; k<numEq; k++) {
                        for (index_t m=0; m<numComp; m++) {
                            const Scalar D_0 = D_p[INDEX2(k, m, numEq)];
                            EM_S[INDEX4(k,m,0,0,numEq,numComp,4)]+=16.*D_0*w[22][l];
                            EM_S[INDEX4(k,m,0,1,numEq,numComp,4)]+= 8.*D_0*w[22][l];
                            EM_S[INDEX4(k,m,0,2,numEq,numComp,4)]+= 8.*D_0*w[22][l];
                            EM_S[INDEX4(k,m,0,3,numEq,numComp,4)]+= 4.*D_0*w[22][l];
                            EM_S[INDEX4(k,m,1,0,numEq,numComp,4)]+= 8.*D_0*w[22][l];
                            EM_S[INDEX4(k,m,1,1,numEq,numComp,4)]+=16.*D_0*w[22][l];
                            EM_S[INDEX4(k,m,1,2,numEq,numComp,4)]+= 4.*D_0*w[22][l];
                            EM_S[INDEX4(k,m,1,3,numEq,numComp,4)]+= 8.*D_0*w[22][l];
                            EM_S[INDEX4(k,m,2,0,numEq,numComp,4)]+= 8.*D_0*w[22][l];
                            EM_S[INDEX4(k,m,2,1,numEq,numComp,4)]+= 4.*D_0*w[22][l];
                            EM_S[INDEX4(k,m,2,2,numEq,numComp,4)]+=16.*D_0*w[22][l];
                            EM_S[INDEX4(k,m,2,3,numEq,numComp,4)]+= 8.*D_0*w[22][l];
                            EM_S[INDEX4(k,m,3,0,numEq,numComp,4)]+= 4.*D_0*w[22][l];
                            EM_S[INDEX4(k,m,3,1,numEq,numComp,4)]+= 8.*D_0*w[22][l];
                            EM_S[INDEX4(k,m,3,2,numEq,numComp,4)]+= 8.*D_0*w[22][l];
                            EM_S[INDEX4(k,m,3,3,numEq,numComp,4)]+=16.*D_0*w[22][l];
                        }
                    }
                }
            }

            ///////////////
            // process X //
            ///////////////
            if (!X.isEmpty()) {
                const Scalar* X_p = X.getSampleDataRO(id, zero);
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
                        const Scalar tmp0 = 6.*w[15][l]*(X_0_2 + X_0_3);
                        const Scalar tmp1 = 6.*w[10][l]*(X_0_0 + X_0_1);
                        const Scalar tmp2 = 6.*w[11][l]*(X_1_0 + X_1_2);
                        const Scalar tmp3 = 6.*w[14][l]*(X_1_1 + X_1_3);
                        const Scalar tmp4 = 6.*w[11][l]*(X_1_1 + X_1_3);
                        const Scalar tmp5 =    w[25][l]*(X_0_0 + X_0_1);
                        const Scalar tmp6 =    w[26][l]*(X_0_2 + X_0_3);
                        const Scalar tmp7 = 6.*w[14][l]*(X_1_0 + X_1_2);
                        const Scalar tmp8 =    w[27][l]*(X_1_0 + X_1_2);
                        const Scalar tmp9 =    w[28][l]*(X_1_1 + X_1_3);
                        const Scalar tmp10 =   w[25][l]*(-X_0_2- X_0_3);
                        const Scalar tmp11 =   w[26][l]*(-X_0_0- X_0_1);
                        const Scalar tmp12 =   w[27][l]*(X_1_1 + X_1_3);
                        const Scalar tmp13 =   w[28][l]*(X_1_0 + X_1_2);
                        const Scalar tmp14 =   w[25][l]*(X_0_2 + X_0_3);
                        const Scalar tmp15 =   w[26][l]*(X_0_0 + X_0_1);
                        EM_F[INDEX2(k,0,numEq)]+=tmp0 + tmp1 + tmp2 + tmp3;
                        EM_F[INDEX2(k,1,numEq)]+=tmp4 + tmp5 + tmp6 + tmp7;
                        EM_F[INDEX2(k,2,numEq)]+=tmp10 + tmp11 + tmp8 + tmp9;
                        EM_F[INDEX2(k,3,numEq)]+=tmp12 + tmp13 + tmp14 + tmp15;
                    }
                } else { // constant data
                    for (index_t k=0; k<numEq; k++) {
                        const Scalar wX0 = X_p[INDEX2(k, 0, numEq)]*w[18][l];
                        const Scalar wX1 = X_p[INDEX2(k, 1, numEq)]*w[19][l];
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
                const Scalar* Y_p = Y.getSampleDataRO(id, zero);
                if (Y.actsExpanded()) {
                    for (index_t k=0; k<numEq; k++) {
                        const Scalar Y_0 = Y_p[INDEX2(k, 0, numEq)];
                        const Scalar Y_1 = Y_p[INDEX2(k, 1, numEq)];
                        const Scalar Y_2 = Y_p[INDEX2(k, 2, numEq)];
                        const Scalar Y_3 = Y_p[INDEX2(k, 3, numEq)];
                        const Scalar tmp0 = 6.*w[22][l]*(Y_1 + Y_2);
                        const Scalar tmp1 = 6.*w[22][l]*(Y_0 + Y_3);
                        EM_F[INDEX2(k,0,numEq)]+=6.*Y_0*w[20][l] + 6.*Y_3*w[21][l] + tmp0;
                        EM_F[INDEX2(k,1,numEq)]+=6.*Y_1*w[20][l] + 6.*Y_2*w[21][l] + tmp1;
                        EM_F[INDEX2(k,2,numEq)]+=6.*Y_1*w[21][l] + 6.*Y_2*w[20][l] + tmp1;
                        EM_F[INDEX2(k,3,numEq)]+=6.*Y_0*w[21][l] + 6.*Y_3*w[20][l] + tmp0;
                    }
                } else { // constant data
                    for (index_t k=0; k<numEq; k++) {
                        EM_F[INDEX2(k,0,numEq)]+=36.*Y_p[k]*w[22][l];
                        EM_F[INDEX2(k,1,numEq)]+=36.*Y_p[k]*w[22][l];
                        EM_F[INDEX2(k,2,numEq)]+=36.*Y_p[k]*w[22][l];
                        EM_F[INDEX2(k,3,numEq)]+=36.*Y_p[k]*w[22][l];
                    }
                }
            }

            // add to matrix (if addEM_S) and RHS (if addEM_F)
            // const index_t firstNode=m_NN[0]*k1+k0;
            // domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, addEM_S, addEM_F, q, t, numEq, numComp);
            domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, addEM_S, addEM_F, q, t, numEq, numComp);
        }
    }
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

    // Find the maximum level of refinement in the mesh
    int max_level = 0;
    for(p4est_topidx_t tree = domain->p4est->first_local_tree; tree <= domain->p4est->last_local_tree; tree++) {
        p4est_tree_t * tree_t = p4est_tree_array_index(domain->p4est->trees, tree);
        max_level = tree_t->maxlevel > max_level ? tree_t->maxlevel : max_level;
    }

    // This is similar to Ripley, except that in Oxley the quads vary in size so that the 
    // weightings need to be calcuated for each level of refinement
    const double SQRT3 = 1.73205080756887719318;
    double w[16][P4EST_MAXLEVEL] = {{0}};
// #pragma omp parallel for
    for(int i = 0; i <= max_level; i++)
    {
        double m_dx[2] = {domain->forestData.m_dx[0][P4EST_MAXLEVEL-i], 
                          domain->forestData.m_dx[1][P4EST_MAXLEVEL-i]};

        w[5][i] = m_dx[0]/12;
        w[6][i] = w[5][i]*( SQRT3 + 2);
        w[7][i] = w[5][i]*(-SQRT3 + 2);
        w[8][i] = w[5][i]*( SQRT3 + 3);
        w[9][i] = w[5][i]*(-SQRT3 + 3);
        w[2][i] = m_dx[1]/12;
        w[0][i] = w[2][i]*( SQRT3 + 2);
        w[1][i] = w[2][i]*(-SQRT3 + 2);
        w[3][i] = w[2][i]*( SQRT3 + 3);
        w[4][i] = w[2][i]*(-SQRT3 + 3);
    }

    const bool addEM_S = !d.isEmpty();
    const bool addEM_F = !y.isEmpty();
    const Scalar zero = static_cast<Scalar>(0);
    rhs.requireWrite();

    vector<Scalar> EM_S(4*4*numEq*numComp, zero);
    vector<Scalar> EM_F(4*numEq, zero);

    if(domain->m_faceOffset[0]>-1)
    {

        if (addEM_S)
            fill(EM_S.begin(), EM_S.end(), zero);
        if (addEM_F)
            fill(EM_F.begin(), EM_F.end(), zero);

        for (index_t k=0; k<domain->NodeIDsLeft.size()-1; k++) {
            int id = k;
            int l = domain->NodeIDsLeft[k].level;

            ///////////////
            // process d //
            ///////////////
            if (addEM_S) {
                const Scalar* d_p = d.getSampleDataRO(id, zero);
                if (d.actsExpanded()) {
                    for (index_t k=0; k<numEq; k++) {
                        for (index_t m=0; m<numComp; m++) {
                            const Scalar d_0 = d_p[INDEX3(k,m,0,numEq,numComp)];
                            const Scalar d_1 = d_p[INDEX3(k,m,1,numEq,numComp)];
                            const Scalar tmp0 = w[2][l]*(d_0 + d_1);
                            EM_S[INDEX4(k,m,0,0,numEq,numComp,4)] = d_0*w[0][l] + d_1*w[1][l];
                            EM_S[INDEX4(k,m,2,0,numEq,numComp,4)] = tmp0;
                            EM_S[INDEX4(k,m,0,2,numEq,numComp,4)] = tmp0;
                            EM_S[INDEX4(k,m,2,2,numEq,numComp,4)] = d_0*w[1][l] + d_1*w[0][l];
                        }
                     }
                } else { // constant data
                    for (index_t k=0; k<numEq; k++) {
                        for (index_t m=0; m<numComp; m++) {
                            const Scalar d_0 = d_p[INDEX2(k, m, numEq)];
                            EM_S[INDEX4(k,m,0,0,numEq,numComp,4)] = 4.*d_0*w[2][l];
                            EM_S[INDEX4(k,m,2,0,numEq,numComp,4)] = 2.*d_0*w[2][l];
                            EM_S[INDEX4(k,m,0,2,numEq,numComp,4)] = 2.*d_0*w[2][l];
                            EM_S[INDEX4(k,m,2,2,numEq,numComp,4)] = 4.*d_0*w[2][l];
                        }
                    }
                }
            }

            ///////////////
            // process y //
            ///////////////
            if (addEM_F) {
                const Scalar* y_p = y.getSampleDataRO(id, zero);
                if (y.actsExpanded()) {
                    for (index_t k=0; k<numEq; k++) {
                        const Scalar y_0 = y_p[INDEX2(k, 0, numEq)];
                        const Scalar y_1 = y_p[INDEX2(k, 1, numEq)];
                        EM_F[INDEX2(k,0,numEq)] = w[3][l]*y_0 + w[4][l]*y_1;
                        EM_F[INDEX2(k,2,numEq)] = w[3][l]*y_1 + w[4][l]*y_0;
                    }
                } else { // constant data
                    for (index_t k=0; k<numEq; k++) {
                        EM_F[INDEX2(k,0,numEq)] = 6*w[2][l]*y_p[k];
                        EM_F[INDEX2(k,2,numEq)] = 6*w[2][l]*y_p[k];
                    }
                }
            }
            domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, addEM_S, addEM_F, domain->NodeIDsLeft[k], numEq, numComp);
        } 
    }

    if(domain->m_faceOffset[1]>-1)
    {
        if (addEM_S)
            fill(EM_S.begin(), EM_S.end(), zero);
        if (addEM_F)
            fill(EM_F.begin(), EM_F.end(), zero);

        for (index_t k=0; k<domain->NodeIDsRight.size()-1; k++) {
            int id = domain->m_faceOffset[1]+k;
            int l = domain->NodeIDsRight[k].level;

            ///////////////
            // process d //
            ///////////////
            if (addEM_S) {
                const Scalar* d_p = d.getSampleDataRO(id, zero);
                if (d.actsExpanded()) {
                    for (index_t k=0; k<numEq; k++) {
                        for (index_t m=0; m<numComp; m++) {
                            const Scalar d_0 = d_p[INDEX3(k,m,0,numEq,numComp)];
                            const Scalar d_1 = d_p[INDEX3(k,m,1,numEq,numComp)];
                            const Scalar tmp0 = w[2][l]*(d_0 + d_1);
                            EM_S[INDEX4(k,m,1,1,numEq,numComp,4)] = d_0*w[0][l] + d_1*w[1][l];
                            EM_S[INDEX4(k,m,3,1,numEq,numComp,4)] = tmp0;
                            EM_S[INDEX4(k,m,1,3,numEq,numComp,4)] = tmp0;
                            EM_S[INDEX4(k,m,3,3,numEq,numComp,4)] = d_0*w[1][l] + d_1*w[0][l];
                        }
                     }
                } else { // constant data
                    for (index_t k=0; k<numEq; k++) {
                        for (index_t m=0; m<numComp; m++) {
                            const Scalar d_0 = d_p[INDEX2(k, m, numEq)];
                            EM_S[INDEX4(k,m,1,1,numEq,numComp,4)] = 4.*d_0*w[2][l];
                            EM_S[INDEX4(k,m,3,1,numEq,numComp,4)] = 2.*d_0*w[2][l];
                            EM_S[INDEX4(k,m,1,3,numEq,numComp,4)] = 2.*d_0*w[2][l];
                            EM_S[INDEX4(k,m,3,3,numEq,numComp,4)] = 4.*d_0*w[2][l];
                        }
                    }
                }
            }

            ///////////////
            // process y //
            ///////////////
            if (addEM_F) {
                const Scalar* y_p = y.getSampleDataRO(id, zero);
                if (y.actsExpanded()) {
                    for (index_t k=0; k<numEq; k++) {
                        const Scalar y_0 = y_p[INDEX2(k, 0, numEq)];
                        const Scalar y_1 = y_p[INDEX2(k, 1, numEq)];
                        EM_F[INDEX2(k,1,numEq)] = w[3][l]*y_0 + w[4][l]*y_1;
                        EM_F[INDEX2(k,3,numEq)] = w[3][l]*y_1 + w[4][l]*y_0;
                    }
                } else { // constant data
                    for (index_t k=0; k<numEq; k++) {
                        EM_F[INDEX2(k,1,numEq)] = 6*w[2][l]*y_p[k];
                        EM_F[INDEX2(k,3,numEq)] = 6*w[2][l]*y_p[k];
                    }
                }
            }
            domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, addEM_S, addEM_F, domain->NodeIDsRight[k], numEq, numComp);
        }
    }

    if(domain->m_faceOffset[2]>-1)
    {

        if (addEM_S)
            fill(EM_S.begin(), EM_S.end(), zero);
        if (addEM_F)
            fill(EM_F.begin(), EM_F.end(), zero);

        for (index_t k=0; k<domain->NodeIDsBottom.size()-1; k++) {
            int id = domain->m_faceOffset[2]+k;
            int l = domain->NodeIDsBottom[k].level;

            ///////////////
            // process d //
            ///////////////
            if (addEM_S) {
                const Scalar* d_p = d.getSampleDataRO(id, zero);
                if (d.actsExpanded()) {
                    for (index_t k=0; k<numEq; k++) {
                        for (index_t m=0; m<numComp; m++) {
                            const Scalar d_0 = d_p[INDEX3(k,m,0,numEq,numComp)];
                            const Scalar d_1 = d_p[INDEX3(k,m,1,numEq,numComp)];
                            const Scalar tmp0 = w[5][l]*(d_0 + d_1);
                            EM_S[INDEX4(k,m,0,0,numEq,numComp,4)] = d_0*w[6][l] + d_1*w[7][l];
                            EM_S[INDEX4(k,m,1,0,numEq,numComp,4)] = tmp0;
                            EM_S[INDEX4(k,m,0,1,numEq,numComp,4)] = tmp0;
                            EM_S[INDEX4(k,m,1,1,numEq,numComp,4)] = d_0*w[7][l] + d_1*w[6][l];
                        }
                     }
                } else { // constant data
                    for (index_t k=0; k<numEq; k++) {
                        for (index_t m=0; m<numComp; m++) {
                            const Scalar d_0 = d_p[INDEX2(k, m, numEq)];
                            EM_S[INDEX4(k,m,0,0,numEq,numComp,4)] = 4.*d_0*w[5][l];
                            EM_S[INDEX4(k,m,1,0,numEq,numComp,4)] = 2.*d_0*w[5][l];
                            EM_S[INDEX4(k,m,0,1,numEq,numComp,4)] = 2.*d_0*w[5][l];
                            EM_S[INDEX4(k,m,1,1,numEq,numComp,4)] = 4.*d_0*w[5][l];
                        }
                    }
                }
            }

            ///////////////
            // process y //
            ///////////////
            if (addEM_F) {
                const Scalar* y_p = y.getSampleDataRO(id, zero);
                if (y.actsExpanded()) {
                    for (index_t k=0; k<numEq; k++) {
                        const Scalar y_0 = y_p[INDEX2(k, 0, numEq)];
                        const Scalar y_1 = y_p[INDEX2(k, 1, numEq)];
                        EM_F[INDEX2(k,0,numEq)] = w[8][l]*y_0 + w[9][l]*y_1;
                        EM_F[INDEX2(k,1,numEq)] = w[8][l]*y_1 + w[9][l]*y_0;
                    }
                } else { // constant data
                    for (index_t k=0; k<numEq; k++) {
                        EM_F[INDEX2(k,0,numEq)] = 6*w[5][l]*y_p[k];
                        EM_F[INDEX2(k,1,numEq)] = 6*w[5][l]*y_p[k];
                    }
                }
            }
            domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, addEM_S, addEM_F, domain->NodeIDsBottom[k], numEq, numComp);
        } 
    }

    if(domain->m_faceOffset[3]>-1)
    {

        if (addEM_S)
            fill(EM_S.begin(), EM_S.end(), zero);
        if (addEM_F)
            fill(EM_F.begin(), EM_F.end(), zero);

        for (index_t k=0; k<domain->NodeIDsTop.size()-1; k++) {
            int id = domain->m_faceOffset[3]+k;
            int l = domain->NodeIDsTop[k].level;

            ///////////////
            // process d //
            ///////////////
            if (addEM_S) {
                const Scalar* d_p = d.getSampleDataRO(id, zero);
                if (d.actsExpanded()) {
                    for (index_t k=0; k<numEq; k++) {
                        for (index_t m=0; m<numComp; m++) {
                            const Scalar d_0 = d_p[INDEX3(k,m,0,numEq,numComp)];
                            const Scalar d_1 = d_p[INDEX3(k,m,1,numEq,numComp)];
                            const Scalar tmp0 = w[5][l]*(d_0 + d_1);
                            EM_S[INDEX4(k,m,2,2,numEq,numComp,4)] = d_0*w[6][l] + d_1*w[7][l];
                            EM_S[INDEX4(k,m,3,2,numEq,numComp,4)] = tmp0;
                            EM_S[INDEX4(k,m,2,3,numEq,numComp,4)] = tmp0;
                            EM_S[INDEX4(k,m,3,3,numEq,numComp,4)] = d_0*w[7][l] + d_1*w[6][l];
                        }
                     }
                } else { // constant data
                    for (index_t k=0; k<numEq; k++) {
                        for (index_t m=0; m<numComp; m++) {
                            const Scalar d_0 = d_p[INDEX2(k, m, numEq)];
                            EM_S[INDEX4(k,m,2,2,numEq,numComp,4)] = 4.*d_0*w[5][l];
                            EM_S[INDEX4(k,m,3,2,numEq,numComp,4)] = 2.*d_0*w[5][l];
                            EM_S[INDEX4(k,m,2,3,numEq,numComp,4)] = 2.*d_0*w[5][l];
                            EM_S[INDEX4(k,m,3,3,numEq,numComp,4)] = 4.*d_0*w[5][l];
                        }
                    }
                }
            }

            ///////////////
            // process y //
            ///////////////
            if (addEM_F) {
                const Scalar* y_p = y.getSampleDataRO(id, zero);
                if (y.actsExpanded()) {
                    for (index_t k=0; k<numEq; k++) {
                        const Scalar y_0 = y_p[INDEX2(k, 0, numEq)];
                        const Scalar y_1 = y_p[INDEX2(k, 1, numEq)];
                        EM_F[INDEX2(k,2,numEq)] = w[8][l]*y_0 + w[9][l]*y_1;
                        EM_F[INDEX2(k,3,numEq)] = w[8][l]*y_1 + w[9][l]*y_0;
                    }
                } else { // constant data
                    for (index_t k=0; k<numEq; k++) {
                        EM_F[INDEX2(k,2,numEq)] = 6.*w[5][l]*y_p[k];
                        EM_F[INDEX2(k,3,numEq)] = 6.*w[5][l]*y_p[k];
                    }
                }
            }
            domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, addEM_S, addEM_F, domain->NodeIDsTop[k], numEq, numComp);
        }
    }
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

    // Find the maximum level of refinement in the mesh
    int max_level = 0;
    for(p4est_topidx_t tree = domain->p4est->first_local_tree; tree <= domain->p4est->last_local_tree; tree++) {
        p4est_tree_t * tree_t = p4est_tree_array_index(domain->p4est->trees, tree);
        max_level = tree_t->maxlevel > max_level ? tree_t->maxlevel : max_level;
    }

    // This is similar to Ripley, except that in Oxley the quads vary in size so that the 
    // weightings need to be calcuated for each level of refinement
    double w[12][P4EST_MAXLEVEL] = {{0}};
// #pragma omp parallel for
    for(int i = 0; i <= max_level; i++)
    {
        double m_dx[2] = {domain->forestData.m_dx[0][P4EST_MAXLEVEL-i], 
                          domain->forestData.m_dx[1][P4EST_MAXLEVEL-i]};

        w[0][i] = 1./4;
        w[1][i] = m_dx[0]/8;
        w[2][i] = m_dx[1]/8;
        w[3][i] = m_dx[0] *  m_dx[1]/16;
        w[4][i] = m_dx[0]/(4*m_dx[1]);
        w[5][i] = m_dx[1]/(4*m_dx[0]);
    }

    const bool addEM_S = (!A.isEmpty() || !B.isEmpty() || !C.isEmpty() || !D.isEmpty());
    const bool addEM_F = (!X.isEmpty() || !Y.isEmpty());
    const Scalar zero = static_cast<Scalar>(0);
    rhs.requireWrite();

    vector<Scalar> EM_S(4*4*numEq*numComp, zero);
    vector<Scalar> EM_F(4*numEq, zero);

// #pragma omp parallel for
    for (p4est_topidx_t t = domain->p4est->first_local_tree; t <= domain->p4est->last_local_tree; t++) // Loop over every tree
    {
        p4est_tree_t * currenttree = p4est_tree_array_index(domain->p4est->trees, t);
        sc_array_t * tquadrants = &currenttree->quadrants;
        p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
        for (int q = 0; q < Q; ++q)  // Loop over the elements attached to the tree
        {                
            p4est_quadrant_t * quad = p4est_quadrant_array_index(tquadrants, q);
            int l = quad->level;
            double xy[3];
            p4est_qcoord_to_vertex(domain->p4est->connectivity, t, quad->x, quad->y, xy);
            long id = domain->getQuadID(domain->NodeIDs.find(std::make_pair(xy[0],xy[1]))->second);

                if (addEM_S)
                    fill(EM_S.begin(), EM_S.end(), zero);
                if (addEM_F)
                    fill(EM_F.begin(), EM_F.end(), zero);

                ///////////////
                // process A //
                ///////////////
                if (!A.isEmpty()) {
                    const Scalar* A_p = A.getSampleDataRO(id, zero);
                    for (index_t k=0; k<numEq; k++) {
                        for (index_t m=0; m<numComp; m++) {
                            const Scalar Aw00 = A_p[INDEX4(k,0,m,0, numEq,2, numComp)]*w[5][l];
                            const Scalar Aw01 = A_p[INDEX4(k,0,m,1, numEq,2, numComp)]*w[0][l];
                            const Scalar Aw10 = A_p[INDEX4(k,1,m,0, numEq,2, numComp)]*w[0][l];
                            const Scalar Aw11 = A_p[INDEX4(k,1,m,1, numEq,2, numComp)]*w[4][l];
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
                    const Scalar* B_p = B.getSampleDataRO(id, zero);
                    for (index_t k=0; k<numEq; k++) {
                        for (index_t m=0; m<numComp; m++) {
                            const Scalar wB0 = B_p[INDEX3(k,0,m, numEq, 2)]*w[2][l];
                            const Scalar wB1 = B_p[INDEX3(k,1,m, numEq, 2)]*w[1][l];
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
                    const Scalar* C_p = C.getSampleDataRO(id, zero);
                    for (index_t k=0; k<numEq; k++) {
                        for (index_t m=0; m<numComp; m++) {
                            const Scalar wC0 = C_p[INDEX3(k, m, 0, numEq, numComp)]*w[2][l];
                            const Scalar wC1 = C_p[INDEX3(k, m, 1, numEq, numComp)]*w[1][l];
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
                    const Scalar* D_p = D.getSampleDataRO(id, zero);
                    for (index_t k=0; k<numEq; k++) {
                        for (index_t m=0; m<numComp; m++) {
                            const Scalar wD0 = D_p[INDEX2(k, m, numEq)]*w[3][l];
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
                    const Scalar* X_p = X.getSampleDataRO(id, zero);
                    for (index_t k=0; k<numEq; k++) {
                        const Scalar wX0 = 4.*X_p[INDEX2(k, 0, numEq)]*w[2][l];
                        const Scalar wX1 = 4.*X_p[INDEX2(k, 1, numEq)]*w[1][l];
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
                    const Scalar* Y_p = Y.getSampleDataRO(id, zero);
                    for (index_t k=0; k<numEq; k++) {
                        EM_F[INDEX2(k,0,numEq)]+=4.*Y_p[k]*w[3][l];
                        EM_F[INDEX2(k,1,numEq)]+=4.*Y_p[k]*w[3][l];
                        EM_F[INDEX2(k,2,numEq)]+=4.*Y_p[k]*w[3][l];
                        EM_F[INDEX2(k,3,numEq)]+=4.*Y_p[k]*w[3][l];
                    }
                }

                // add to matrix (if addEM_S) and RHS (if addEM_F)
                // domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, addEM_S, addEM_F, id, numEq, numComp);
                domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, addEM_S, addEM_F, q, t, numEq, numComp);
        } 
    } 
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

    // Find the maximum level of refinement in the mesh
    int max_level = 0;
    for(p4est_topidx_t tree = domain->p4est->first_local_tree; tree <= domain->p4est->last_local_tree; tree++) {
        p4est_tree_t * tree_t = p4est_tree_array_index(domain->p4est->trees, tree);
        max_level = tree_t->maxlevel > max_level ? tree_t->maxlevel : max_level;
    }

    // This is similar to Ripley, except that in Oxley the quads vary in size so that the 
    // weightings need to be calcuated for each level of refinement
    double w[2][P4EST_MAXLEVEL] = {{0}};
// #pragma omp parallel for
    for(int i = 0; i <= max_level; i++)
    {
        double m_dx[2] = {domain->forestData.m_dx[0][P4EST_MAXLEVEL-i], 
                          domain->forestData.m_dx[1][P4EST_MAXLEVEL-i]};

        w[0][i] = m_dx[0]/4;
        w[1][i] = m_dx[1]/4;
    }

    const bool addEM_S = !d.isEmpty();
    const bool addEM_F = !y.isEmpty();
    const Scalar zero = static_cast<Scalar>(0);
    rhs.requireWrite();

    vector<Scalar> EM_S(4*4*numEq*numComp, zero);
    vector<Scalar> EM_F(4*numEq, zero);

    if(domain->m_faceOffset[0]>-1)
    {
        if (addEM_S)
            fill(EM_S.begin(), EM_S.end(), zero);
        if (addEM_F)
            fill(EM_F.begin(), EM_F.end(), zero);

        for (index_t k=0; k<domain->NodeIDsLeft.size()-1; k++) {
            int id = k;
            int l = domain->NodeIDsLeft[k].level;

            ///////////////
            // process d //
            ///////////////
            if (addEM_S) {
                const Scalar* d_p = d.getSampleDataRO(id, zero);
                for (index_t k=0; k<numEq; k++) {
                    for (index_t m=0; m<numComp; m++) {
                        const Scalar tmp0 = d_p[INDEX2(k, m, numEq)]*w[1][l];
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
                const Scalar* y_p = y.getSampleDataRO(id, zero);
                for (index_t k=0; k<numEq; k++) {
                    EM_F[INDEX2(k,0,numEq)] = 2*w[1][l]*y_p[k];
                    EM_F[INDEX2(k,2,numEq)] = 2*w[1][l]*y_p[k];
                }
            }
            domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, addEM_S, addEM_F, domain->NodeIDsLeft[k], numEq, numComp);
        }
    }

    if(domain->m_faceOffset[1]>-1)
    {
        if (addEM_S)
            fill(EM_S.begin(), EM_S.end(), zero);
        if (addEM_F)
            fill(EM_F.begin(), EM_F.end(), zero);

        for (index_t k=0; k<domain->NodeIDsRight.size()-1; k++) {
            int id = domain->m_faceOffset[1]+k;
            int l = domain->NodeIDsRight[k].level;

            ///////////////
            // process d //
            ///////////////
            if (addEM_S) {
                const Scalar* d_p = d.getSampleDataRO(id, zero);
                for (index_t k=0; k<numEq; k++) {
                    for (index_t m=0; m<numComp; m++) {
                        const Scalar tmp0 = d_p[INDEX2(k, m, numEq)]*w[1][l];
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
                const Scalar* y_p = y.getSampleDataRO(id, zero);
                for (index_t k=0; k<numEq; k++) {
                    EM_F[INDEX2(k,1,numEq)] = 2.*w[1][l]*y_p[k];
                    EM_F[INDEX2(k,3,numEq)] = 2.*w[1][l]*y_p[k];
                }
            }
            domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, addEM_S, addEM_F, domain->NodeIDsRight[k], numEq, numComp);
        }
    }

    if(domain->m_faceOffset[2]>-1)
    {
        if (addEM_S)
            fill(EM_S.begin(), EM_S.end(), zero);
        if (addEM_F)
            fill(EM_F.begin(), EM_F.end(), zero);

        for (index_t k=0; k<domain->NodeIDsBottom.size()-1; k++) {
            int id = domain->m_faceOffset[2]+k;
            int l = domain->NodeIDsBottom[k].level;

            ///////////////
            // process d //
            ///////////////
            if (addEM_S) {
                const Scalar* d_p = d.getSampleDataRO(id, zero);
                for (index_t k=0; k<numEq; k++) {
                    for (index_t m=0; m<numComp; m++) {
                        const Scalar tmp0 = d_p[INDEX2(k, m, numEq)]*w[0][l];
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
                const Scalar* y_p = y.getSampleDataRO(id, zero);
                for (index_t k=0; k<numEq; k++) {
                    EM_F[INDEX2(k,0,numEq)] = 2.*w[0][l]*y_p[k];
                    EM_F[INDEX2(k,1,numEq)] = 2.*w[0][l]*y_p[k];
                }
            }
            domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, addEM_S, addEM_F, domain->NodeIDsBottom[k], numEq, numComp);
        }
    }

    if(domain->m_faceOffset[3]>-1)
    {
        if (addEM_S)
            fill(EM_S.begin(), EM_S.end(), zero);
        if (addEM_F)
            fill(EM_F.begin(), EM_F.end(), zero);

        for (index_t k=0; k<domain->NodeIDsTop.size()-1; k++) {
            int id = domain->m_faceOffset[3]+k;
            int l = domain->NodeIDsTop[k].level;

            ///////////////
            // process d //
            ///////////////
            if (addEM_S) {
                const Scalar* d_p = d.getSampleDataRO(id, zero);
                for (index_t k=0; k<numEq; k++) {
                    for (index_t m=0; m<numComp; m++) {
                        const Scalar tmp0 = d_p[INDEX2(k, m, numEq)]*w[0][l];
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
                const Scalar* y_p = y.getSampleDataRO(id, zero);
                for (index_t k=0; k<numEq; k++) {
                    EM_F[INDEX2(k,2,numEq)] = 2*w[0][l]*y_p[k];
                    EM_F[INDEX2(k,3,numEq)] = 2*w[0][l]*y_p[k];
                }
            }
            domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, addEM_S, addEM_F, domain->NodeIDsTop[k], numEq, numComp);
        }
    }
}

/****************************************************************************/
// PDE HANGING NODES
/****************************************************************************/
#ifdef ESYS_HAVE_TRILINOS
template<class Scalar>
void DefaultAssembler2D<Scalar>::assemblePDEHanging(
                    Teuchos::RCP<Tpetra::CrsMatrix<double,esys_trilinos::LO,esys_trilinos::GO,esys_trilinos::NT>>* IZ) const
{
    // const Tpetra::Vector<>::scalar_type one  = static_cast<Tpetra::Vector<>::scalar_type> (1.0);
    // const Tpetra::Vector<>::scalar_type half = static_cast<Tpetra::Vector<>::scalar_type> (0.5);

     // A(0, 0:1) = [2, -1]
        //     if (gblRow == 0) {
        //         A->insertGlobalValues (gblRow,
        //             tuple<global_ordinal_type> (gblRow, gblRow + 1),
        //             tuple<scalar_type> (two, negOne));

    // Create I
    // for(esys_trilinos::GO i = 0; i < domain->getNumNodes()-0.5*domain->num_hanging; i++)
    //     IZ->insertGlobalValues(i,i,one);

    // Loop over hanging nodes
    // std::vector<LongPair> hanging_faces=domain->hanging_faces;

    // for(int i = 0; i < domain->num_hanging; i++)
    // {
    //     int a, b;
    //     for(int j = 0; j < hanging_faces.size(); j++)
    //     {
    //         if(hanging_faces[j].first > hanging_faces[i].first)
    //         {
    //             a=j;
    //             break;
    //         }
    //     }
    //     for(int j = 0; j < hanging_faces.size(); j++)
    //     {
    //         if(hanging_faces[j].first > hanging_faces[i].first)
    //         {
    //             b=j;
    //             break;
    //         }
    //     }

        // IZ->insertGlobalValues(a,b,half);
        // IZ->insertGlobalValues(b,a,half);
    // }

    // IZ->fillComplete();

}

// instantiate the supported templates
template class DefaultAssembler2D<escript::DataTypes::real_t>;
template class DefaultAssembler2D<escript::DataTypes::cplx_t>;

#endif //ESYS_HAVE_TRILINOS

} // namespace oxley

