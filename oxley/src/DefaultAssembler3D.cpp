
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

#include <oxley/DefaultAssembler3D.h>
#include <oxley/domainhelpers.h>

#include <escript/DataTypes.h>
#include <escript/index.h>

using namespace std;
using escript::AbstractSystemMatrix;
using escript::Data;

namespace oxley {

template<class Scalar>
void DefaultAssembler3D<Scalar>::collateFunctionSpaceTypes(
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
void DefaultAssembler3D<Scalar>::assemblePDESingle(AbstractSystemMatrix* mat,
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
void DefaultAssembler3D<Scalar>::assemblePDEBoundarySingle(
                                        AbstractSystemMatrix* mat, Data& rhs,
                                        const DataMap& coefs) const 
{
    const Data& d = unpackData("d", coefs);
    const Data& y = unpackData("y", coefs);
    assemblePDEBoundarySingle(mat, rhs, d, y);
}

template<class Scalar>
void DefaultAssembler3D<Scalar>::assemblePDESingleReduced(
                                        AbstractSystemMatrix* mat, Data& rhs,
                                        const DataMap& coefs) const
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
void DefaultAssembler3D<Scalar>::assemblePDEBoundarySingleReduced(
                                        AbstractSystemMatrix* mat, Data& rhs,
                                        const DataMap& coefs) const
{
    const Data& d = unpackData("d", coefs);
    const Data& y = unpackData("y", coefs);
    assemblePDEBoundarySingleReduced(mat, rhs, d, y);
}

template<class Scalar>
void DefaultAssembler3D<Scalar>::assemblePDESystem(AbstractSystemMatrix* mat,
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
void DefaultAssembler3D<Scalar>::assemblePDEBoundarySystem(
                                        AbstractSystemMatrix* mat, Data& rhs,
                                        const DataMap& coefs) const
{
    const Data& d = unpackData("d", coefs);
    const Data& y = unpackData("y", coefs);
    assemblePDEBoundarySystem(mat, rhs, d, y);
}

template<class Scalar>
void DefaultAssembler3D<Scalar>::assemblePDESystemReduced(
                                        AbstractSystemMatrix* mat, Data& rhs,
                                        const DataMap& coefs) const
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
void DefaultAssembler3D<Scalar>::assemblePDEBoundarySystemReduced(
                                        AbstractSystemMatrix* mat,
                                        Data& rhs, const DataMap& coefs) const
{
    const Data& d = unpackData("d", coefs);
    const Data& y = unpackData("y", coefs);
    assemblePDEBoundarySystemReduced(mat, rhs, d, y);
}

#ifdef ESYS_HAVE_TRILINOS
template<class Scalar>
void DefaultAssembler3D<Scalar>::assemblePDEHanging(Teuchos::RCP<Tpetra::CrsMatrix<double,esys_trilinos::LO,esys_trilinos::GO,esys_trilinos::NT>>* mat) const
{

}
#endif

/****************************************************************************/
// PDE SINGLE
/****************************************************************************/

template<class Scalar>
void DefaultAssembler3D<Scalar>::assemblePDESingle(AbstractSystemMatrix* mat,
                                       Data& rhs, const Data& A, const Data& B,
                                       const Data& C, const Data& D,
                                       const Data& X, const Data& Y) const
{
    // Find the maximum level of refinement in the mesh
    int max_level = 0;
    for(p8est_topidx_t tree = domain->p8est->first_local_tree; tree <= domain->p8est->last_local_tree; tree++) 
    {
        p8est_tree_t * tree_t = p8est_tree_array_index(domain->p8est->trees, tree);
        max_level = tree_t->maxlevel > max_level ? tree_t->maxlevel : max_level;
    }

    const double SQRT3 = 1.73205080756887719318;
    double w[74][P4EST_MAXLEVEL] = {{0}};
#pragma omp parallel for
    for(int i = 0; i < max_level; i++)
    {
        double m_dx[3] = {domain->m_NX[0]*domain->forestData->m_dx[0][P4EST_MAXLEVEL-i], 
                          domain->m_NX[1]*domain->forestData->m_dx[1][P4EST_MAXLEVEL-i],
                          domain->m_NX[2]*domain->forestData->m_dx[2][P4EST_MAXLEVEL-i]};

        w[10][i] = -m_dx[0]/288;
        w[6][i]  = w[10][i]*(SQRT3 - 2);
        w[12][i] = w[10][i]*(-SQRT3 - 2);
        w[4][i]  = w[10][i]*(-4*SQRT3 + 7);
        w[18][i] = w[10][i]*(-4*SQRT3 - 7);
        w[11][i] = m_dx[1]/288;
        w[5][i]  = w[11][i]*(-SQRT3 + 2);
        w[15][i] = w[11][i]*(SQRT3 + 2);
        w[2][i]  = w[11][i]*(4*SQRT3 - 7);
        w[17][i] = w[11][i]*(4*SQRT3 + 7);
        w[8][i]  = m_dx[2]/288;
        w[1][i]  = w[8][i]*(-SQRT3 + 2);
        w[16][i] = w[8][i]*(SQRT3 + 2);
        w[20][i] = w[8][i]*(4*SQRT3 - 7);
        w[21][i] = w[8][i]*(-4*SQRT3 - 7);
        w[50][i] =  m_dx[0]*m_dx[1]/72;
        w[65][i] = -m_dx[0]*m_dx[1]/48;
        w[35][i] = w[65][i]*(-SQRT3 - 3)/36;
        w[42][i] = w[65][i]*(SQRT3 - 3)/36;
        w[32][i] = w[65][i]*(5*SQRT3 - 9)/36;
        w[43][i] = w[65][i]*(-5*SQRT3 - 9)/36;
        w[40][i] = w[65][i]*(-19*SQRT3 - 33)/36;
        w[41][i] = w[65][i]*(19*SQRT3 - 33)/36;
        w[63][i] = w[65][i]*(SQRT3 + 2);
        w[67][i] = w[65][i]*(-SQRT3 + 2);
        w[51][i] = -m_dx[0]*m_dx[2]/72;
        w[64][i] = -m_dx[0]*m_dx[2]/48;
        w[34][i] = w[64][i]*(-SQRT3 - 3)/36;
        w[37][i] = w[64][i]*(SQRT3 - 3)/36;
        w[31][i] = w[64][i]*(5*SQRT3 - 9)/36;
        w[39][i] = w[64][i]*(-5*SQRT3 - 9)/36;
        w[44][i] = w[64][i]*(19*SQRT3 + 33)/36;
        w[45][i] = w[64][i]*(-19*SQRT3 + 33)/36;
        w[62][i] = w[64][i]*(SQRT3 + 2);
        w[68][i] = w[64][i]*(-SQRT3 + 2);
        w[53][i] = -m_dx[1]*m_dx[2]/72;
        w[66][i] = -m_dx[1]*m_dx[2]/48;
        w[33][i] = w[66][i]*(SQRT3 - 3)/36;
        w[36][i] = w[66][i]*(-SQRT3 - 3)/36;
        w[30][i] = w[66][i]*(5*SQRT3 - 9)/36;
        w[38][i] = w[66][i]*(-5*SQRT3 - 9)/36;
        w[46][i] = w[66][i]*(19*SQRT3 - 33)/36;
        w[47][i] = w[66][i]*(-19*SQRT3 - 33)/36;
        w[61][i] = w[66][i]*(SQRT3 + 2);
        w[69][i] = w[66][i]*(-SQRT3 + 2);
        w[55][i] = m_dx[0]*m_dx[1]*m_dx[2]/1728;
        w[57][i] = w[55][i]*(-SQRT3 + 2);
        w[58][i] = w[55][i]*(SQRT3 + 2);
        w[54][i] = w[55][i]*(-4*SQRT3 + 7);
        w[56][i] = w[55][i]*(4*SQRT3 + 7);
        w[59][i] = w[55][i]*(15*SQRT3 + 26);
        w[60][i] = w[55][i]*(-15*SQRT3 + 26);
        w[71][i] = w[55][i]*6*(SQRT3 + 3);
        w[72][i] = w[55][i]*6*(-SQRT3 + 3);
        w[70][i] = w[55][i]*6*(5*SQRT3 + 9);
        w[73][i] = w[55][i]*6*(-5*SQRT3 + 9);
        w[13][i] = -m_dx[0]*m_dx[1]/(288*m_dx[2]);
        w[23][i] = w[13][i]*(SQRT3 - 2);
        w[25][i] = w[13][i]*(-SQRT3 - 2);
        w[7][i]  = w[13][i]*(-4*SQRT3 + 7);
        w[19][i] = w[13][i]*(4*SQRT3 + 7);
        w[22][i] = -m_dx[0]*m_dx[2]/(288*m_dx[1]);
        w[3][i]  = w[22][i]*(SQRT3 - 2);
        w[9][i]  = w[22][i]*(-SQRT3 - 2);
        w[24][i] = w[22][i]*(4*SQRT3 + 7);
        w[26][i] = w[22][i]*(-4*SQRT3 + 7);
        w[27][i] = -m_dx[1]*m_dx[2]/(288*m_dx[0]);
        w[0][i]  = w[27][i]*(SQRT3 - 2);
        w[14][i] = w[27][i]*(-SQRT3 - 2);
        w[28][i] = w[27][i]*(-4*SQRT3 + 7);
        w[29][i] = w[27][i]*(4*SQRT3 + 7);
    }

    // const dim_t NE0 = m_NE[0];
    // const dim_t NE1 = m_NE[1];
    // const dim_t NE2 = m_NE[2];

    const bool add_EM_S = (!A.isEmpty() || !B.isEmpty() || !C.isEmpty() || !D.isEmpty());
    const bool add_EM_F = (!X.isEmpty() || !Y.isEmpty());
    const Scalar zero = static_cast<Scalar>(0);

    std::vector<Scalar> EM_S(8*8);
    std::vector<Scalar> EM_F(8);
    if (add_EM_S)
        fill(EM_S.begin(), EM_S.end(), zero);
    if (add_EM_F)
        fill(EM_F.begin(), EM_F.end(), zero);
    rhs.requireWrite();

    for (p8est_topidx_t t = domain->p8est->first_local_tree; t <= domain->p8est->last_local_tree; t++) // Loop over every tree
    {
        p8est_tree_t * currenttree = p8est_tree_array_index(domain->p8est->trees, t);
        sc_array_t * tquadrants = &currenttree->quadrants;
        p4est_locidx_t Q = (p4est_locidx_t) tquadrants->elem_count;
#pragma omp parallel for
        for (int q = 0; q < Q; ++q)  // Loop over the elements attached to the tree
        {
            if (add_EM_S)
                fill(EM_S.begin(), EM_S.end(), zero);
            if (add_EM_F)
                fill(EM_F.begin(), EM_F.end(), zero);
            
            p8est_quadrant_t * quad = p8est_quadrant_array_index(tquadrants, q);
            int l = quad->level;
            double xyz[3];
            p8est_qcoord_to_vertex(domain->p8est->connectivity, t, quad->x, quad->y, quad->z, xyz);
            long id = domain->NodeIDs.find(std::make_tuple(xyz[0],xyz[1],xyz[2]))->second;
                        
            ///////////////
            // process A //
            ///////////////
            if (!A.isEmpty()) {
                const Scalar* A_p = A.getSampleDataRO(id, zero);
                if (A.actsExpanded()) {
                    const Scalar A_00_0 = A_p[INDEX3(0,0,0,3,3)];
                    const Scalar A_01_0 = A_p[INDEX3(0,1,0,3,3)];
                    const Scalar A_02_0 = A_p[INDEX3(0,2,0,3,3)];
                    const Scalar A_10_0 = A_p[INDEX3(1,0,0,3,3)];
                    const Scalar A_11_0 = A_p[INDEX3(1,1,0,3,3)];
                    const Scalar A_12_0 = A_p[INDEX3(1,2,0,3,3)];
                    const Scalar A_20_0 = A_p[INDEX3(2,0,0,3,3)];
                    const Scalar A_21_0 = A_p[INDEX3(2,1,0,3,3)];
                    const Scalar A_22_0 = A_p[INDEX3(2,2,0,3,3)];
                    const Scalar A_00_1 = A_p[INDEX3(0,0,1,3,3)];
                    const Scalar A_01_1 = A_p[INDEX3(0,1,1,3,3)];
                    const Scalar A_02_1 = A_p[INDEX3(0,2,1,3,3)];
                    const Scalar A_10_1 = A_p[INDEX3(1,0,1,3,3)];
                    const Scalar A_11_1 = A_p[INDEX3(1,1,1,3,3)];
                    const Scalar A_12_1 = A_p[INDEX3(1,2,1,3,3)];
                    const Scalar A_20_1 = A_p[INDEX3(2,0,1,3,3)];
                    const Scalar A_21_1 = A_p[INDEX3(2,1,1,3,3)];
                    const Scalar A_22_1 = A_p[INDEX3(2,2,1,3,3)];
                    const Scalar A_00_2 = A_p[INDEX3(0,0,2,3,3)];
                    const Scalar A_01_2 = A_p[INDEX3(0,1,2,3,3)];
                    const Scalar A_02_2 = A_p[INDEX3(0,2,2,3,3)];
                    const Scalar A_10_2 = A_p[INDEX3(1,0,2,3,3)];
                    const Scalar A_11_2 = A_p[INDEX3(1,1,2,3,3)];
                    const Scalar A_12_2 = A_p[INDEX3(1,2,2,3,3)];
                    const Scalar A_20_2 = A_p[INDEX3(2,0,2,3,3)];
                    const Scalar A_21_2 = A_p[INDEX3(2,1,2,3,3)];
                    const Scalar A_22_2 = A_p[INDEX3(2,2,2,3,3)];
                    const Scalar A_00_3 = A_p[INDEX3(0,0,3,3,3)];
                    const Scalar A_01_3 = A_p[INDEX3(0,1,3,3,3)];
                    const Scalar A_02_3 = A_p[INDEX3(0,2,3,3,3)];
                    const Scalar A_10_3 = A_p[INDEX3(1,0,3,3,3)];
                    const Scalar A_11_3 = A_p[INDEX3(1,1,3,3,3)];
                    const Scalar A_12_3 = A_p[INDEX3(1,2,3,3,3)];
                    const Scalar A_20_3 = A_p[INDEX3(2,0,3,3,3)];
                    const Scalar A_21_3 = A_p[INDEX3(2,1,3,3,3)];
                    const Scalar A_22_3 = A_p[INDEX3(2,2,3,3,3)];
                    const Scalar A_00_4 = A_p[INDEX3(0,0,4,3,3)];
                    const Scalar A_01_4 = A_p[INDEX3(0,1,4,3,3)];
                    const Scalar A_02_4 = A_p[INDEX3(0,2,4,3,3)];
                    const Scalar A_10_4 = A_p[INDEX3(1,0,4,3,3)];
                    const Scalar A_11_4 = A_p[INDEX3(1,1,4,3,3)];
                    const Scalar A_12_4 = A_p[INDEX3(1,2,4,3,3)];
                    const Scalar A_20_4 = A_p[INDEX3(2,0,4,3,3)];
                    const Scalar A_21_4 = A_p[INDEX3(2,1,4,3,3)];
                    const Scalar A_22_4 = A_p[INDEX3(2,2,4,3,3)];
                    const Scalar A_00_5 = A_p[INDEX3(0,0,5,3,3)];
                    const Scalar A_01_5 = A_p[INDEX3(0,1,5,3,3)];
                    const Scalar A_02_5 = A_p[INDEX3(0,2,5,3,3)];
                    const Scalar A_10_5 = A_p[INDEX3(1,0,5,3,3)];
                    const Scalar A_11_5 = A_p[INDEX3(1,1,5,3,3)];
                    const Scalar A_12_5 = A_p[INDEX3(1,2,5,3,3)];
                    const Scalar A_20_5 = A_p[INDEX3(2,0,5,3,3)];
                    const Scalar A_21_5 = A_p[INDEX3(2,1,5,3,3)];
                    const Scalar A_22_5 = A_p[INDEX3(2,2,5,3,3)];
                    const Scalar A_00_6 = A_p[INDEX3(0,0,6,3,3)];
                    const Scalar A_01_6 = A_p[INDEX3(0,1,6,3,3)];
                    const Scalar A_02_6 = A_p[INDEX3(0,2,6,3,3)];
                    const Scalar A_10_6 = A_p[INDEX3(1,0,6,3,3)];
                    const Scalar A_11_6 = A_p[INDEX3(1,1,6,3,3)];
                    const Scalar A_12_6 = A_p[INDEX3(1,2,6,3,3)];
                    const Scalar A_20_6 = A_p[INDEX3(2,0,6,3,3)];
                    const Scalar A_21_6 = A_p[INDEX3(2,1,6,3,3)];
                    const Scalar A_22_6 = A_p[INDEX3(2,2,6,3,3)];
                    const Scalar A_00_7 = A_p[INDEX3(0,0,7,3,3)];
                    const Scalar A_01_7 = A_p[INDEX3(0,1,7,3,3)];
                    const Scalar A_02_7 = A_p[INDEX3(0,2,7,3,3)];
                    const Scalar A_10_7 = A_p[INDEX3(1,0,7,3,3)];
                    const Scalar A_11_7 = A_p[INDEX3(1,1,7,3,3)];
                    const Scalar A_12_7 = A_p[INDEX3(1,2,7,3,3)];
                    const Scalar A_20_7 = A_p[INDEX3(2,0,7,3,3)];
                    const Scalar A_21_7 = A_p[INDEX3(2,1,7,3,3)];
                    const Scalar A_22_7 = A_p[INDEX3(2,2,7,3,3)];
                    const Scalar tmp0 = w[18][l]*(-A_12_7 + A_21_3);
                    const Scalar tmp1 = w[13][l]*(A_22_1 + A_22_2 + A_22_5 + A_22_6);
                    const Scalar tmp2 = w[11][l]*(-A_02_2 - A_02_5 + A_20_1 + A_20_6);
                    const Scalar tmp3 = w[14][l]*(A_00_2 + A_00_3 + A_00_6 + A_00_7);
                    const Scalar tmp4 = w[7][l]*(A_22_0 + A_22_4);
                    const Scalar tmp5 = w[10][l]*(A_12_1 + A_12_6 - A_21_2 - A_21_5);
                    const Scalar tmp6 = w[3][l]*(A_11_0 + A_11_2 + A_11_4 + A_11_6);
                    const Scalar tmp7 = w[1][l]*(A_01_0 + A_01_4 + A_10_0 + A_10_4);
                    const Scalar tmp8 = w[4][l]*(A_12_0 - A_21_4);
                    const Scalar tmp9 = w[15][l]*(-A_02_3 - A_02_6 + A_20_2 + A_20_7);
                    const Scalar tmp10 = w[0][l]*(A_00_0 + A_00_1 + A_00_4 + A_00_5);
                    const Scalar tmp11 = w[16][l]*(A_01_3 + A_01_7 + A_10_3 + A_10_7);
                    const Scalar tmp12 = w[9][l]*(A_11_1 + A_11_3 + A_11_5 + A_11_7);
                    const Scalar tmp13 = w[12][l]*(-A_12_3 - A_12_5 + A_21_1 + A_21_7);
                    const Scalar tmp14 = w[5][l]*(-A_02_1 - A_02_4 + A_20_0 + A_20_5);
                    const Scalar tmp15 = w[8][l]*(A_01_1 + A_01_2 + A_01_5 + A_01_6 + A_10_1 + A_10_2 + A_10_5 + A_10_6);
                    const Scalar tmp16 = w[6][l]*(-A_12_2 - A_12_4 + A_21_0 + A_21_6);
                    const Scalar tmp17 = w[19][l]*(A_22_3 + A_22_7);
                    const Scalar tmp18 = w[17][l]*(-A_02_7 + A_20_3);
                    const Scalar tmp19 = w[2][l]*(A_02_0 - A_20_4);
                    const Scalar tmp20 = w[13][l]*(-A_22_0 - A_22_1 - A_22_2 - A_22_3 - A_22_4 - A_22_5 - A_22_6 - A_22_7);
                    const Scalar tmp21 = w[11][l]*(-A_02_1 - A_02_3 - A_02_4 - A_02_6 + A_20_0 + A_20_2 + A_20_5 + A_20_7);
                    const Scalar tmp22 = w[14][l]*(-A_00_4 - A_00_5 - A_00_6 - A_00_7);
                    const Scalar tmp23 = w[20][l]*(A_01_2 + A_10_1);
                    const Scalar tmp24 = w[10][l]*(A_12_2 + A_12_3 + A_12_4 + A_12_5 - A_21_0 - A_21_1 - A_21_6 - A_21_7);
                    const Scalar tmp25 = w[3][l]*(-A_11_0 - A_11_1 - A_11_2 - A_11_3);
                    const Scalar tmp26 = w[1][l]*(-A_01_0 - A_01_3 - A_10_0 - A_10_3);
                    const Scalar tmp27 = w[15][l]*(-A_02_5 - A_02_7 + A_20_4 + A_20_6);
                    const Scalar tmp28 = w[0][l]*(-A_00_0 - A_00_1 - A_00_2 - A_00_3);
                    const Scalar tmp29 = w[16][l]*(-A_01_4 - A_01_7 - A_10_4 - A_10_7);
                    const Scalar tmp30 = w[9][l]*(-A_11_4 - A_11_5 - A_11_6 - A_11_7);
                    const Scalar tmp31 = w[21][l]*(A_01_5 + A_10_6);
                    const Scalar tmp32 = w[12][l]*(-A_12_6 - A_12_7 + A_21_4 + A_21_5);
                    const Scalar tmp33 = w[5][l]*(-A_02_0 - A_02_2 + A_20_1 + A_20_3);
                    const Scalar tmp34 = w[8][l]*(-A_01_1 - A_01_6 - A_10_2 - A_10_5);
                    const Scalar tmp35 = w[6][l]*(-A_12_0 - A_12_1 + A_21_2 + A_21_3);
                    const Scalar tmp36 = w[20][l]*(-A_01_6 + A_10_4);
                    const Scalar tmp37 = w[18][l]*(A_12_3 - A_21_1);
                    const Scalar tmp38 = w[11][l]*(-A_02_0 - A_02_2 - A_02_5 - A_02_7 - A_20_0 - A_20_2 - A_20_5 - A_20_7);
                    const Scalar tmp39 = w[14][l]*(A_00_0 + A_00_1 + A_00_2 + A_00_3);
                    const Scalar tmp40 = w[26][l]*(A_11_4 + A_11_6);
                    const Scalar tmp41 = w[0][l]*(A_00_4 + A_00_5 + A_00_6 + A_00_7);
                    const Scalar tmp42 = w[10][l]*(-A_12_2 - A_12_5 + A_21_0 + A_21_7);
                    const Scalar tmp43 = w[22][l]*(A_11_0 + A_11_2 + A_11_5 + A_11_7);
                    const Scalar tmp44 = w[1][l]*(A_01_4 + A_01_7 - A_10_5 - A_10_6);
                    const Scalar tmp45 = w[25][l]*(A_22_1 + A_22_3 + A_22_5 + A_22_7);
                    const Scalar tmp46 = w[4][l]*(-A_12_4 + A_21_6);
                    const Scalar tmp47 = w[15][l]*(-A_02_1 - A_02_3 - A_20_1 - A_20_3);
                    const Scalar tmp48 = w[21][l]*(-A_01_1 + A_10_3);
                    const Scalar tmp49 = w[16][l]*(A_01_0 + A_01_3 - A_10_1 - A_10_2);
                    const Scalar tmp50 = w[5][l]*(-A_02_4 - A_02_6 - A_20_4 - A_20_6);
                    const Scalar tmp51 = w[12][l]*(A_12_1 + A_12_7 - A_21_3 - A_21_5);
                    const Scalar tmp52 = w[24][l]*(A_11_1 + A_11_3);
                    const Scalar tmp53 = w[8][l]*(A_01_2 + A_01_5 - A_10_0 - A_10_7);
                    const Scalar tmp54 = w[6][l]*(A_12_0 + A_12_6 - A_21_2 - A_21_4);
                    const Scalar tmp55 = w[23][l]*(A_22_0 + A_22_2 + A_22_4 + A_22_6);
                    const Scalar tmp56 = w[18][l]*(A_12_4 - A_21_6);
                    const Scalar tmp57 = w[14][l]*(A_00_4 + A_00_5 + A_00_6 + A_00_7);
                    const Scalar tmp58 = w[26][l]*(A_11_1 + A_11_3);
                    const Scalar tmp59 = w[20][l]*(-A_01_1 + A_10_3);
                    const Scalar tmp60 = w[1][l]*(A_01_0 + A_01_3 - A_10_1 - A_10_2);
                    const Scalar tmp61 = w[25][l]*(A_22_0 + A_22_2 + A_22_4 + A_22_6);
                    const Scalar tmp62 = w[4][l]*(-A_12_3 + A_21_1);
                    const Scalar tmp63 = w[15][l]*(-A_02_4 - A_02_6 - A_20_4 - A_20_6);
                    const Scalar tmp64 = w[0][l]*(A_00_0 + A_00_1 + A_00_2 + A_00_3);
                    const Scalar tmp65 = w[16][l]*(A_01_4 + A_01_7 - A_10_5 - A_10_6);
                    const Scalar tmp66 = w[24][l]*(A_11_4 + A_11_6);
                    const Scalar tmp67 = w[21][l]*(-A_01_6 + A_10_4);
                    const Scalar tmp68 = w[12][l]*(A_12_0 + A_12_6 - A_21_2 - A_21_4);
                    const Scalar tmp69 = w[5][l]*(-A_02_1 - A_02_3 - A_20_1 - A_20_3);
                    const Scalar tmp70 = w[6][l]*(A_12_1 + A_12_7 - A_21_3 - A_21_5);
                    const Scalar tmp71 = w[23][l]*(A_22_1 + A_22_3 + A_22_5 + A_22_7);
                    const Scalar tmp72 = w[20][l]*(A_01_5 + A_10_6);
                    const Scalar tmp73 = w[14][l]*(-A_00_0 - A_00_1 - A_00_2 - A_00_3);
                    const Scalar tmp74 = w[0][l]*(-A_00_4 - A_00_5 - A_00_6 - A_00_7);
                    const Scalar tmp75 = w[3][l]*(-A_11_4 - A_11_5 - A_11_6 - A_11_7);
                    const Scalar tmp76 = w[1][l]*(-A_01_4 - A_01_7 - A_10_4 - A_10_7);
                    const Scalar tmp77 = w[15][l]*(-A_02_0 - A_02_2 + A_20_1 + A_20_3);
                    const Scalar tmp78 = w[21][l]*(A_01_2 + A_10_1);
                    const Scalar tmp79 = w[16][l]*(-A_01_0 - A_01_3 - A_10_0 - A_10_3);
                    const Scalar tmp80 = w[9][l]*(-A_11_0 - A_11_1 - A_11_2 - A_11_3);
                    const Scalar tmp81 = w[12][l]*(-A_12_0 - A_12_1 + A_21_2 + A_21_3);
                    const Scalar tmp82 = w[5][l]*(-A_02_5 - A_02_7 + A_20_4 + A_20_6);
                    const Scalar tmp83 = w[6][l]*(-A_12_6 - A_12_7 + A_21_4 + A_21_5);
                    const Scalar tmp84 = w[6][l]*(-A_12_2 - A_12_3 - A_21_2 - A_21_3);
                    const Scalar tmp85 = w[11][l]*(A_02_1 + A_02_6 - A_20_0 - A_20_7);
                    const Scalar tmp86 = w[20][l]*(A_01_3 - A_10_2);
                    const Scalar tmp87 = w[10][l]*(A_12_0 + A_12_1 + A_12_6 + A_12_7 + A_21_0 + A_21_1 + A_21_6 + A_21_7);
                    const Scalar tmp88 = w[3][l]*(A_11_0 + A_11_1 + A_11_2 + A_11_3);
                    const Scalar tmp89 = w[23][l]*(A_22_2 + A_22_3 + A_22_6 + A_22_7);
                    const Scalar tmp90 = w[1][l]*(-A_01_1 - A_01_2 + A_10_0 + A_10_3);
                    const Scalar tmp91 = w[25][l]*(A_22_0 + A_22_1 + A_22_4 + A_22_5);
                    const Scalar tmp92 = w[15][l]*(A_02_0 + A_02_5 - A_20_1 - A_20_4);
                    const Scalar tmp93 = w[21][l]*(A_01_4 - A_10_5);
                    const Scalar tmp94 = w[16][l]*(-A_01_5 - A_01_6 + A_10_4 + A_10_7);
                    const Scalar tmp95 = w[28][l]*(A_00_2 + A_00_3);
                    const Scalar tmp96 = w[12][l]*(-A_12_4 - A_12_5 - A_21_4 - A_21_5);
                    const Scalar tmp97 = w[29][l]*(A_00_4 + A_00_5);
                    const Scalar tmp98 = w[5][l]*(A_02_2 + A_02_7 - A_20_3 - A_20_6);
                    const Scalar tmp99 = w[8][l]*(-A_01_0 - A_01_7 + A_10_1 + A_10_6);
                    const Scalar tmp100 = w[9][l]*(A_11_4 + A_11_5 + A_11_6 + A_11_7);
                    const Scalar tmp101 = w[27][l]*(A_00_0 + A_00_1 + A_00_6 + A_00_7);
                    const Scalar tmp102 = w[17][l]*(A_02_4 - A_20_5);
                    const Scalar tmp103 = w[2][l]*(-A_02_3 + A_20_2);
                    const Scalar tmp104 = w[13][l]*(A_22_0 + A_22_1 + A_22_2 + A_22_3 + A_22_4 + A_22_5 + A_22_6 + A_22_7);
                    const Scalar tmp105 = w[6][l]*(-A_12_4 - A_12_5 - A_21_2 - A_21_3);
                    const Scalar tmp106 = w[22][l]*(A_11_0 + A_11_1 + A_11_2 + A_11_3 + A_11_4 + A_11_5 + A_11_6 + A_11_7);
                    const Scalar tmp107 = w[1][l]*(-A_01_2 - A_01_6 - A_10_1 - A_10_5);
                    const Scalar tmp108 = w[15][l]*(-A_02_1 - A_02_3 - A_20_4 - A_20_6);
                    const Scalar tmp109 = w[16][l]*(-A_01_1 - A_01_5 - A_10_2 - A_10_6);
                    const Scalar tmp110 = w[12][l]*(-A_12_2 - A_12_3 - A_21_4 - A_21_5);
                    const Scalar tmp111 = w[5][l]*(-A_02_4 - A_02_6 - A_20_1 - A_20_3);
                    const Scalar tmp112 = w[8][l]*(-A_01_0 - A_01_3 - A_01_4 - A_01_7 - A_10_0 - A_10_3 - A_10_4 - A_10_7);
                    const Scalar tmp113 = w[27][l]*(A_00_0 + A_00_1 + A_00_2 + A_00_3 + A_00_4 + A_00_5 + A_00_6 + A_00_7);
                    const Scalar tmp114 = w[11][l]*(A_02_0 + A_02_2 + A_02_5 + A_02_7 - A_20_1 - A_20_3 - A_20_4 - A_20_6);
                    const Scalar tmp115 = w[21][l]*(-A_01_4 - A_10_7);
                    const Scalar tmp116 = w[20][l]*(-A_01_3 - A_10_0);
                    const Scalar tmp117 = w[15][l]*(A_02_4 + A_02_6 - A_20_5 - A_20_7);
                    const Scalar tmp118 = w[16][l]*(A_01_5 + A_01_6 + A_10_5 + A_10_6);
                    const Scalar tmp119 = w[5][l]*(A_02_1 + A_02_3 - A_20_0 - A_20_2);
                    const Scalar tmp120 = w[8][l]*(A_01_0 + A_01_7 + A_10_3 + A_10_4);
                    const Scalar tmp121 = w[1][l]*(A_01_1 + A_01_2 + A_10_1 + A_10_2);
                    const Scalar tmp122 = w[18][l]*(A_12_2 - A_21_6);
                    const Scalar tmp123 = w[13][l]*(A_22_0 + A_22_3 + A_22_4 + A_22_7);
                    const Scalar tmp124 = w[11][l]*(-A_02_0 - A_02_7 + A_20_3 + A_20_4);
                    const Scalar tmp125 = w[7][l]*(A_22_1 + A_22_5);
                    const Scalar tmp126 = w[10][l]*(-A_12_3 - A_12_4 + A_21_0 + A_21_7);
                    const Scalar tmp127 = w[3][l]*(A_11_1 + A_11_3 + A_11_5 + A_11_7);
                    const Scalar tmp128 = w[1][l]*(-A_01_1 - A_01_5 - A_10_1 - A_10_5);
                    const Scalar tmp129 = w[4][l]*(-A_12_5 + A_21_1);
                    const Scalar tmp130 = w[16][l]*(-A_01_2 - A_01_6 - A_10_2 - A_10_6);
                    const Scalar tmp131 = w[9][l]*(A_11_0 + A_11_2 + A_11_4 + A_11_6);
                    const Scalar tmp132 = w[19][l]*(A_22_2 + A_22_6);
                    const Scalar tmp133 = w[17][l]*(-A_02_2 + A_20_6);
                    const Scalar tmp134 = w[2][l]*(A_02_5 - A_20_1);
                    const Scalar tmp135 = w[11][l]*(A_02_1 + A_02_3 + A_02_4 + A_02_6 + A_20_1 + A_20_3 + A_20_4 + A_20_6);
                    const Scalar tmp136 = w[1][l]*(A_01_3 + A_01_7 + A_10_0 + A_10_4);
                    const Scalar tmp137 = w[15][l]*(A_02_0 + A_02_2 + A_20_5 + A_20_7);
                    const Scalar tmp138 = w[16][l]*(A_01_0 + A_01_4 + A_10_3 + A_10_7);
                    const Scalar tmp139 = w[5][l]*(A_02_5 + A_02_7 + A_20_0 + A_20_2);
                    const Scalar tmp140 = w[18][l]*(A_12_5 - A_21_1);
                    const Scalar tmp141 = w[14][l]*(A_00_0 + A_00_1 + A_00_4 + A_00_5);
                    const Scalar tmp142 = w[7][l]*(A_22_2 + A_22_6);
                    const Scalar tmp143 = w[1][l]*(-A_01_2 - A_01_6 - A_10_2 - A_10_6);
                    const Scalar tmp144 = w[4][l]*(-A_12_2 + A_21_6);
                    const Scalar tmp145 = w[15][l]*(-A_02_1 - A_02_4 + A_20_0 + A_20_5);
                    const Scalar tmp146 = w[0][l]*(A_00_2 + A_00_3 + A_00_6 + A_00_7);
                    const Scalar tmp147 = w[16][l]*(-A_01_1 - A_01_5 - A_10_1 - A_10_5);
                    const Scalar tmp148 = w[5][l]*(-A_02_3 - A_02_6 + A_20_2 + A_20_7);
                    const Scalar tmp149 = w[19][l]*(A_22_1 + A_22_5);
                    const Scalar tmp150 = w[17][l]*(-A_02_5 + A_20_1);
                    const Scalar tmp151 = w[2][l]*(A_02_2 - A_20_6);
                    const Scalar tmp152 = w[18][l]*(A_12_3 - A_21_7);
                    const Scalar tmp153 = w[11][l]*(A_02_1 + A_02_6 - A_20_2 - A_20_5);
                    const Scalar tmp154 = w[10][l]*(-A_12_2 - A_12_5 + A_21_1 + A_21_6);
                    const Scalar tmp155 = w[4][l]*(-A_12_4 + A_21_0);
                    const Scalar tmp156 = w[15][l]*(A_02_2 + A_02_7 - A_20_3 - A_20_6);
                    const Scalar tmp157 = w[5][l]*(A_02_0 + A_02_5 - A_20_1 - A_20_4);
                    const Scalar tmp158 = w[17][l]*(A_02_3 - A_20_7);
                    const Scalar tmp159 = w[2][l]*(-A_02_4 + A_20_0);
                    const Scalar tmp160 = w[6][l]*(A_12_6 + A_12_7 + A_21_0 + A_21_1);
                    const Scalar tmp161 = w[10][l]*(-A_12_2 - A_12_3 - A_12_4 - A_12_5 - A_21_2 - A_21_3 - A_21_4 - A_21_5);
                    const Scalar tmp162 = w[1][l]*(A_01_0 + A_01_4 + A_10_3 + A_10_7);
                    const Scalar tmp163 = w[16][l]*(A_01_3 + A_01_7 + A_10_0 + A_10_4);
                    const Scalar tmp164 = w[12][l]*(A_12_0 + A_12_1 + A_21_6 + A_21_7);
                    const Scalar tmp165 = w[20][l]*(A_01_6 + A_10_5);
                    const Scalar tmp166 = w[10][l]*(-A_12_0 - A_12_1 - A_12_6 - A_12_7 + A_21_2 + A_21_3 + A_21_4 + A_21_5);
                    const Scalar tmp167 = w[15][l]*(A_02_1 + A_02_3 - A_20_0 - A_20_2);
                    const Scalar tmp168 = w[21][l]*(A_01_1 + A_10_2);
                    const Scalar tmp169 = w[12][l]*(A_12_2 + A_12_3 - A_21_0 - A_21_1);
                    const Scalar tmp170 = w[5][l]*(A_02_4 + A_02_6 - A_20_5 - A_20_7);
                    const Scalar tmp171 = w[8][l]*(-A_01_2 - A_01_5 - A_10_1 - A_10_6);
                    const Scalar tmp172 = w[6][l]*(A_12_4 + A_12_5 - A_21_6 - A_21_7);
                    const Scalar tmp173 = w[2][l]*(A_02_1 + A_20_4);
                    const Scalar tmp174 = w[11][l]*(-A_02_3 - A_02_4 - A_20_1 - A_20_6);
                    const Scalar tmp175 = w[14][l]*(-A_00_2 - A_00_3 - A_00_6 - A_00_7);
                    const Scalar tmp176 = w[22][l]*(-A_11_0 - A_11_1 - A_11_2 - A_11_3 - A_11_4 - A_11_5 - A_11_6 - A_11_7);
                    const Scalar tmp177 = w[1][l]*(A_01_1 + A_01_5 - A_10_0 - A_10_4);
                    const Scalar tmp178 = w[25][l]*(-A_22_2 - A_22_3 - A_22_6 - A_22_7);
                    const Scalar tmp179 = w[15][l]*(-A_02_2 - A_02_7 - A_20_2 - A_20_7);
                    const Scalar tmp180 = w[0][l]*(-A_00_0 - A_00_1 - A_00_4 - A_00_5);
                    const Scalar tmp181 = w[16][l]*(A_01_2 + A_01_6 - A_10_3 - A_10_7);
                    const Scalar tmp182 = w[12][l]*(-A_12_6 - A_12_7 + A_21_2 + A_21_3);
                    const Scalar tmp183 = w[5][l]*(-A_02_0 - A_02_5 - A_20_0 - A_20_5);
                    const Scalar tmp184 = w[8][l]*(A_01_0 + A_01_3 + A_01_4 + A_01_7 - A_10_1 - A_10_2 - A_10_5 - A_10_6);
                    const Scalar tmp185 = w[6][l]*(-A_12_0 - A_12_1 + A_21_4 + A_21_5);
                    const Scalar tmp186 = w[17][l]*(-A_02_6 - A_20_3);
                    const Scalar tmp187 = w[23][l]*(-A_22_0 - A_22_1 - A_22_4 - A_22_5);
                    const Scalar tmp188 = w[18][l]*(A_12_4 - A_21_0);
                    const Scalar tmp189 = w[7][l]*(A_22_3 + A_22_7);
                    const Scalar tmp190 = w[1][l]*(A_01_3 + A_01_7 + A_10_3 + A_10_7);
                    const Scalar tmp191 = w[4][l]*(-A_12_3 + A_21_7);
                    const Scalar tmp192 = w[16][l]*(A_01_0 + A_01_4 + A_10_0 + A_10_4);
                    const Scalar tmp193 = w[19][l]*(A_22_0 + A_22_4);
                    const Scalar tmp194 = w[17][l]*(A_02_4 - A_20_0);
                    const Scalar tmp195 = w[2][l]*(-A_02_3 + A_20_7);
                    const Scalar tmp196 = w[20][l]*(-A_01_7 - A_10_4);
                    const Scalar tmp197 = w[21][l]*(-A_01_0 - A_10_3);
                    const Scalar tmp198 = w[16][l]*(A_01_1 + A_01_2 + A_10_1 + A_10_2);
                    const Scalar tmp199 = w[8][l]*(A_01_3 + A_01_4 + A_10_0 + A_10_7);
                    const Scalar tmp200 = w[1][l]*(A_01_5 + A_01_6 + A_10_5 + A_10_6);
                    const Scalar tmp201 = w[27][l]*(A_00_2 + A_00_3 + A_00_4 + A_00_5);
                    const Scalar tmp202 = w[11][l]*(-A_02_2 - A_02_5 + A_20_3 + A_20_4);
                    const Scalar tmp203 = w[20][l]*(A_01_0 - A_10_1);
                    const Scalar tmp204 = w[23][l]*(A_22_0 + A_22_1 + A_22_4 + A_22_5);
                    const Scalar tmp205 = w[25][l]*(A_22_2 + A_22_3 + A_22_6 + A_22_7);
                    const Scalar tmp206 = w[21][l]*(A_01_7 - A_10_6);
                    const Scalar tmp207 = w[12][l]*(A_12_6 + A_12_7 + A_21_6 + A_21_7);
                    const Scalar tmp208 = w[28][l]*(A_00_0 + A_00_1);
                    const Scalar tmp209 = w[29][l]*(A_00_6 + A_00_7);
                    const Scalar tmp210 = w[8][l]*(-A_01_3 - A_01_4 + A_10_2 + A_10_5);
                    const Scalar tmp211 = w[6][l]*(A_12_0 + A_12_1 + A_21_0 + A_21_1);
                    const Scalar tmp212 = w[17][l]*(-A_02_7 + A_20_6);
                    const Scalar tmp213 = w[2][l]*(A_02_0 - A_20_1);
                    const Scalar tmp214 = w[13][l]*(-A_22_1 - A_22_2 - A_22_5 - A_22_6);
                    const Scalar tmp215 = w[22][l]*(-A_11_0 - A_11_2 - A_11_5 - A_11_7);
                    const Scalar tmp216 = w[8][l]*(A_01_0 + A_01_7 + A_10_0 + A_10_7);
                    const Scalar tmp217 = w[27][l]*(-A_00_0 - A_00_1 - A_00_6 - A_00_7);
                    const Scalar tmp218 = w[17][l]*(-A_02_3 - A_20_3);
                    const Scalar tmp219 = w[2][l]*(A_02_4 + A_20_4);
                    const Scalar tmp220 = w[11][l]*(-A_02_1 - A_02_6 - A_20_1 - A_20_6);
                    const Scalar tmp221 = w[26][l]*(-A_11_4 - A_11_6);
                    const Scalar tmp222 = w[10][l]*(A_12_2 + A_12_5 + A_21_2 + A_21_5);
                    const Scalar tmp223 = w[20][l]*(-A_01_4 - A_10_4);
                    const Scalar tmp224 = w[21][l]*(-A_01_3 - A_10_3);
                    const Scalar tmp225 = w[6][l]*(-A_12_0 - A_12_6 - A_21_0 - A_21_6);
                    const Scalar tmp226 = w[7][l]*(-A_22_0 - A_22_4);
                    const Scalar tmp227 = w[24][l]*(-A_11_1 - A_11_3);
                    const Scalar tmp228 = w[19][l]*(-A_22_3 - A_22_7);
                    const Scalar tmp229 = w[18][l]*(-A_12_3 - A_21_3);
                    const Scalar tmp230 = w[4][l]*(A_12_4 + A_21_4);
                    const Scalar tmp231 = w[28][l]*(-A_00_4 - A_00_5);
                    const Scalar tmp232 = w[12][l]*(-A_12_1 - A_12_7 - A_21_1 - A_21_7);
                    const Scalar tmp233 = w[29][l]*(-A_00_2 - A_00_3);
                    const Scalar tmp234 = w[20][l]*(-A_01_5 + A_10_7);
                    const Scalar tmp235 = w[18][l]*(-A_12_0 + A_21_2);
                    const Scalar tmp236 = w[26][l]*(A_11_5 + A_11_7);
                    const Scalar tmp237 = w[10][l]*(A_12_1 + A_12_6 - A_21_3 - A_21_4);
                    const Scalar tmp238 = w[22][l]*(A_11_1 + A_11_3 + A_11_4 + A_11_6);
                    const Scalar tmp239 = w[4][l]*(A_12_7 - A_21_5);
                    const Scalar tmp240 = w[15][l]*(A_02_0 + A_02_2 + A_20_0 + A_20_2);
                    const Scalar tmp241 = w[21][l]*(-A_01_2 + A_10_0);
                    const Scalar tmp242 = w[5][l]*(A_02_5 + A_02_7 + A_20_5 + A_20_7);
                    const Scalar tmp243 = w[12][l]*(-A_12_2 - A_12_4 + A_21_0 + A_21_6);
                    const Scalar tmp244 = w[24][l]*(A_11_0 + A_11_2);
                    const Scalar tmp245 = w[8][l]*(A_01_1 + A_01_6 - A_10_3 - A_10_4);
                    const Scalar tmp246 = w[6][l]*(-A_12_3 - A_12_5 + A_21_1 + A_21_7);
                    const Scalar tmp247 = w[11][l]*(A_02_3 + A_02_4 - A_20_2 - A_20_5);
                    const Scalar tmp248 = w[20][l]*(-A_01_1 + A_10_0);
                    const Scalar tmp249 = w[21][l]*(-A_01_6 + A_10_7);
                    const Scalar tmp250 = w[8][l]*(A_01_2 + A_01_5 - A_10_3 - A_10_4);
                    const Scalar tmp251 = w[17][l]*(A_02_6 - A_20_7);
                    const Scalar tmp252 = w[2][l]*(-A_02_1 + A_20_0);
                    const Scalar tmp253 = w[17][l]*(-A_02_4 - A_20_4);
                    const Scalar tmp254 = w[2][l]*(A_02_3 + A_20_3);
                    const Scalar tmp255 = w[26][l]*(-A_11_1 - A_11_3);
                    const Scalar tmp256 = w[20][l]*(-A_01_3 - A_10_3);
                    const Scalar tmp257 = w[21][l]*(-A_01_4 - A_10_4);
                    const Scalar tmp258 = w[6][l]*(-A_12_1 - A_12_7 - A_21_1 - A_21_7);
                    const Scalar tmp259 = w[7][l]*(-A_22_3 - A_22_7);
                    const Scalar tmp260 = w[15][l]*(-A_02_0 - A_02_5 - A_20_0 - A_20_5);
                    const Scalar tmp261 = w[24][l]*(-A_11_4 - A_11_6);
                    const Scalar tmp262 = w[19][l]*(-A_22_0 - A_22_4);
                    const Scalar tmp263 = w[18][l]*(-A_12_4 - A_21_4);
                    const Scalar tmp264 = w[4][l]*(A_12_3 + A_21_3);
                    const Scalar tmp265 = w[28][l]*(-A_00_2 - A_00_3);
                    const Scalar tmp266 = w[12][l]*(-A_12_0 - A_12_6 - A_21_0 - A_21_6);
                    const Scalar tmp267 = w[5][l]*(-A_02_2 - A_02_7 - A_20_2 - A_20_7);
                    const Scalar tmp268 = w[29][l]*(-A_00_4 - A_00_5);
                    const Scalar tmp269 = w[11][l]*(A_02_2 + A_02_5 + A_20_0 + A_20_7);
                    const Scalar tmp270 = w[1][l]*(-A_01_0 - A_01_4 + A_10_1 + A_10_5);
                    const Scalar tmp271 = w[15][l]*(A_02_3 + A_02_6 + A_20_3 + A_20_6);
                    const Scalar tmp272 = w[16][l]*(-A_01_3 - A_01_7 + A_10_2 + A_10_6);
                    const Scalar tmp273 = w[5][l]*(A_02_1 + A_02_4 + A_20_1 + A_20_4);
                    const Scalar tmp274 = w[8][l]*(-A_01_1 - A_01_2 - A_01_5 - A_01_6 + A_10_0 + A_10_3 + A_10_4 + A_10_7);
                    const Scalar tmp275 = w[17][l]*(A_02_7 + A_20_2);
                    const Scalar tmp276 = w[2][l]*(-A_02_0 - A_20_5);
                    const Scalar tmp277 = w[18][l]*(-A_12_1 + A_21_5);
                    const Scalar tmp278 = w[11][l]*(A_02_3 + A_02_4 - A_20_0 - A_20_7);
                    const Scalar tmp279 = w[10][l]*(A_12_0 + A_12_7 - A_21_3 - A_21_4);
                    const Scalar tmp280 = w[4][l]*(A_12_6 - A_21_2);
                    const Scalar tmp281 = w[17][l]*(A_02_1 - A_20_5);
                    const Scalar tmp282 = w[2][l]*(-A_02_6 + A_20_2);
                    const Scalar tmp283 = w[11][l]*(A_02_0 + A_02_7 + A_20_2 + A_20_5);
                    const Scalar tmp284 = w[12][l]*(A_12_2 + A_12_3 - A_21_6 - A_21_7);
                    const Scalar tmp285 = w[6][l]*(A_12_4 + A_12_5 - A_21_0 - A_21_1);
                    const Scalar tmp286 = w[17][l]*(A_02_2 + A_20_7);
                    const Scalar tmp287 = w[2][l]*(-A_02_5 - A_20_0);
                    const Scalar tmp288 = w[13][l]*(-A_22_0 - A_22_3 - A_22_4 - A_22_7);
                    const Scalar tmp289 = w[22][l]*(-A_11_1 - A_11_3 - A_11_4 - A_11_6);
                    const Scalar tmp290 = w[8][l]*(-A_01_1 - A_01_6 - A_10_1 - A_10_6);
                    const Scalar tmp291 = w[17][l]*(A_02_2 + A_20_2);
                    const Scalar tmp292 = w[2][l]*(-A_02_5 - A_20_5);
                    const Scalar tmp293 = w[11][l]*(A_02_0 + A_02_7 + A_20_0 + A_20_7);
                    const Scalar tmp294 = w[26][l]*(-A_11_5 - A_11_7);
                    const Scalar tmp295 = w[10][l]*(A_12_3 + A_12_4 + A_21_3 + A_21_4);
                    const Scalar tmp296 = w[20][l]*(A_01_5 + A_10_5);
                    const Scalar tmp297 = w[21][l]*(A_01_2 + A_10_2);
                    const Scalar tmp298 = w[7][l]*(-A_22_1 - A_22_5);
                    const Scalar tmp299 = w[24][l]*(-A_11_0 - A_11_2);
                    const Scalar tmp300 = w[19][l]*(-A_22_2 - A_22_6);
                    const Scalar tmp301 = w[18][l]*(-A_12_2 - A_21_2);
                    const Scalar tmp302 = w[4][l]*(A_12_5 + A_21_5);
                    const Scalar tmp303 = w[8][l]*(A_01_3 + A_01_4 + A_10_3 + A_10_4);
                    const Scalar tmp304 = w[27][l]*(-A_00_2 - A_00_3 - A_00_4 - A_00_5);
                    const Scalar tmp305 = w[17][l]*(A_02_7 + A_20_7);
                    const Scalar tmp306 = w[2][l]*(-A_02_0 - A_20_0);
                    const Scalar tmp307 = w[11][l]*(A_02_2 + A_02_5 + A_20_2 + A_20_5);
                    const Scalar tmp308 = w[26][l]*(-A_11_0 - A_11_2);
                    const Scalar tmp309 = w[10][l]*(-A_12_1 - A_12_6 - A_21_1 - A_21_6);
                    const Scalar tmp310 = w[20][l]*(-A_01_0 - A_10_0);
                    const Scalar tmp311 = w[21][l]*(-A_01_7 - A_10_7);
                    const Scalar tmp312 = w[6][l]*(A_12_2 + A_12_4 + A_21_2 + A_21_4);
                    const Scalar tmp313 = w[24][l]*(-A_11_5 - A_11_7);
                    const Scalar tmp314 = w[18][l]*(A_12_7 + A_21_7);
                    const Scalar tmp315 = w[4][l]*(-A_12_0 - A_21_0);
                    const Scalar tmp316 = w[28][l]*(-A_00_0 - A_00_1);
                    const Scalar tmp317 = w[12][l]*(A_12_3 + A_12_5 + A_21_3 + A_21_5);
                    const Scalar tmp318 = w[29][l]*(-A_00_6 - A_00_7);
                    const Scalar tmp319 = w[18][l]*(-A_12_7 + A_21_5);
                    const Scalar tmp320 = w[26][l]*(A_11_0 + A_11_2);
                    const Scalar tmp321 = w[21][l]*(-A_01_5 + A_10_7);
                    const Scalar tmp322 = w[20][l]*(-A_01_2 + A_10_0);
                    const Scalar tmp323 = w[4][l]*(A_12_0 - A_21_2);
                    const Scalar tmp324 = w[15][l]*(A_02_5 + A_02_7 + A_20_5 + A_20_7);
                    const Scalar tmp325 = w[24][l]*(A_11_5 + A_11_7);
                    const Scalar tmp326 = w[5][l]*(A_02_0 + A_02_2 + A_20_0 + A_20_2);
                    const Scalar tmp327 = w[18][l]*(A_12_7 + A_21_1);
                    const Scalar tmp328 = w[10][l]*(-A_12_1 - A_12_6 - A_21_0 - A_21_7);
                    const Scalar tmp329 = w[3][l]*(-A_11_0 - A_11_2 - A_11_4 - A_11_6);
                    const Scalar tmp330 = w[1][l]*(A_01_2 + A_01_6 - A_10_0 - A_10_4);
                    const Scalar tmp331 = w[4][l]*(-A_12_0 - A_21_6);
                    const Scalar tmp332 = w[25][l]*(-A_22_1 - A_22_3 - A_22_5 - A_22_7);
                    const Scalar tmp333 = w[15][l]*(-A_02_5 - A_02_7 + A_20_1 + A_20_3);
                    const Scalar tmp334 = w[16][l]*(A_01_1 + A_01_5 - A_10_3 - A_10_7);
                    const Scalar tmp335 = w[9][l]*(-A_11_1 - A_11_3 - A_11_5 - A_11_7);
                    const Scalar tmp336 = w[5][l]*(-A_02_0 - A_02_2 + A_20_4 + A_20_6);
                    const Scalar tmp337 = w[27][l]*(-A_00_0 - A_00_1 - A_00_2 - A_00_3 - A_00_4 - A_00_5 - A_00_6 - A_00_7);
                    const Scalar tmp338 = w[23][l]*(-A_22_0 - A_22_2 - A_22_4 - A_22_6);
                    const Scalar tmp339 = w[14][l]*(-A_00_0 - A_00_1 - A_00_4 - A_00_5);
                    const Scalar tmp340 = w[23][l]*(-A_22_2 - A_22_3 - A_22_6 - A_22_7);
                    const Scalar tmp341 = w[1][l]*(A_01_2 + A_01_6 - A_10_3 - A_10_7);
                    const Scalar tmp342 = w[25][l]*(-A_22_0 - A_22_1 - A_22_4 - A_22_5);
                    const Scalar tmp343 = w[15][l]*(A_02_1 + A_02_4 + A_20_1 + A_20_4);
                    const Scalar tmp344 = w[0][l]*(-A_00_2 - A_00_3 - A_00_6 - A_00_7);
                    const Scalar tmp345 = w[16][l]*(A_01_1 + A_01_5 - A_10_0 - A_10_4);
                    const Scalar tmp346 = w[12][l]*(A_12_4 + A_12_5 - A_21_0 - A_21_1);
                    const Scalar tmp347 = w[5][l]*(A_02_3 + A_02_6 + A_20_3 + A_20_6);
                    const Scalar tmp348 = w[6][l]*(A_12_2 + A_12_3 - A_21_6 - A_21_7);
                    const Scalar tmp349 = w[17][l]*(A_02_5 + A_20_0);
                    const Scalar tmp350 = w[2][l]*(-A_02_2 - A_20_7);
                    const Scalar tmp351 = w[8][l]*(-A_01_2 - A_01_5 - A_10_2 - A_10_5);
                    const Scalar tmp352 = w[17][l]*(-A_02_1 - A_20_1);
                    const Scalar tmp353 = w[2][l]*(A_02_6 + A_20_6);
                    const Scalar tmp354 = w[11][l]*(-A_02_3 - A_02_4 - A_20_3 - A_20_4);
                    const Scalar tmp355 = w[10][l]*(-A_12_0 - A_12_7 - A_21_0 - A_21_7);
                    const Scalar tmp356 = w[20][l]*(A_01_6 + A_10_6);
                    const Scalar tmp357 = w[21][l]*(A_01_1 + A_10_1);
                    const Scalar tmp358 = w[7][l]*(-A_22_2 - A_22_6);
                    const Scalar tmp359 = w[19][l]*(-A_22_1 - A_22_5);
                    const Scalar tmp360 = w[18][l]*(A_12_1 + A_21_1);
                    const Scalar tmp361 = w[4][l]*(-A_12_6 - A_21_6);
                    const Scalar tmp362 = w[28][l]*(-A_00_6 - A_00_7);
                    const Scalar tmp363 = w[29][l]*(-A_00_0 - A_00_1);
                    const Scalar tmp364 = w[2][l]*(A_02_4 + A_20_1);
                    const Scalar tmp365 = w[11][l]*(-A_02_1 - A_02_6 - A_20_3 - A_20_4);
                    const Scalar tmp366 = w[17][l]*(-A_02_3 - A_20_6);
                    const Scalar tmp367 = w[2][l]*(A_02_5 - A_20_4);
                    const Scalar tmp368 = w[6][l]*(-A_12_4 - A_12_5 - A_21_4 - A_21_5);
                    const Scalar tmp369 = w[11][l]*(-A_02_0 - A_02_7 + A_20_1 + A_20_6);
                    const Scalar tmp370 = w[20][l]*(-A_01_5 + A_10_4);
                    const Scalar tmp371 = w[3][l]*(A_11_4 + A_11_5 + A_11_6 + A_11_7);
                    const Scalar tmp372 = w[12][l]*(-A_12_2 - A_12_3 - A_21_2 - A_21_3);
                    const Scalar tmp373 = w[21][l]*(-A_01_2 + A_10_3);
                    const Scalar tmp374 = w[9][l]*(A_11_0 + A_11_1 + A_11_2 + A_11_3);
                    const Scalar tmp375 = w[29][l]*(A_00_2 + A_00_3);
                    const Scalar tmp376 = w[8][l]*(A_01_1 + A_01_6 - A_10_0 - A_10_7);
                    const Scalar tmp377 = w[28][l]*(A_00_4 + A_00_5);
                    const Scalar tmp378 = w[17][l]*(-A_02_2 + A_20_3);
                    const Scalar tmp379 = w[17][l]*(A_02_0 + A_20_0);
                    const Scalar tmp380 = w[2][l]*(-A_02_7 - A_20_7);
                    const Scalar tmp381 = w[20][l]*(-A_01_7 - A_10_7);
                    const Scalar tmp382 = w[21][l]*(-A_01_0 - A_10_0);
                    const Scalar tmp383 = w[6][l]*(A_12_3 + A_12_5 + A_21_3 + A_21_5);
                    const Scalar tmp384 = w[18][l]*(A_12_0 + A_21_0);
                    const Scalar tmp385 = w[4][l]*(-A_12_7 - A_21_7);
                    const Scalar tmp386 = w[12][l]*(A_12_2 + A_12_4 + A_21_2 + A_21_4);
                    const Scalar tmp387 = w[17][l]*(-A_02_6 - A_20_6);
                    const Scalar tmp388 = w[2][l]*(A_02_1 + A_20_1);
                    const Scalar tmp389 = w[20][l]*(A_01_1 + A_10_1);
                    const Scalar tmp390 = w[21][l]*(A_01_6 + A_10_6);
                    const Scalar tmp391 = w[18][l]*(A_12_6 + A_21_6);
                    const Scalar tmp392 = w[4][l]*(-A_12_1 - A_21_1);
                    const Scalar tmp393 = w[2][l]*(A_02_3 + A_20_6);
                    const Scalar tmp394 = w[1][l]*(-A_01_3 - A_01_7 + A_10_2 + A_10_6);
                    const Scalar tmp395 = w[16][l]*(-A_01_0 - A_01_4 + A_10_1 + A_10_5);
                    const Scalar tmp396 = w[17][l]*(-A_02_4 - A_20_1);
                    const Scalar tmp397 = w[18][l]*(-A_12_5 - A_21_3);
                    const Scalar tmp398 = w[10][l]*(A_12_3 + A_12_4 + A_21_2 + A_21_5);
                    const Scalar tmp399 = w[1][l]*(-A_01_0 - A_01_4 + A_10_2 + A_10_6);
                    const Scalar tmp400 = w[4][l]*(A_12_2 + A_21_4);
                    const Scalar tmp401 = w[16][l]*(-A_01_3 - A_01_7 + A_10_1 + A_10_5);
                    const Scalar tmp402 = w[20][l]*(-A_01_2 + A_10_3);
                    const Scalar tmp403 = w[21][l]*(-A_01_5 + A_10_4);
                    const Scalar tmp404 = w[17][l]*(-A_02_5 + A_20_4);
                    const Scalar tmp405 = w[2][l]*(A_02_2 - A_20_3);
                    const Scalar tmp406 = w[18][l]*(-A_12_0 + A_21_4);
                    const Scalar tmp407 = w[4][l]*(A_12_7 - A_21_3);
                    const Scalar tmp408 = w[17][l]*(-A_02_0 + A_20_4);
                    const Scalar tmp409 = w[2][l]*(A_02_7 - A_20_3);
                    const Scalar tmp410 = w[17][l]*(A_02_5 + A_20_5);
                    const Scalar tmp411 = w[2][l]*(-A_02_2 - A_20_2);
                    const Scalar tmp412 = w[20][l]*(A_01_2 + A_10_2);
                    const Scalar tmp413 = w[21][l]*(A_01_5 + A_10_5);
                    const Scalar tmp414 = w[18][l]*(-A_12_5 - A_21_5);
                    const Scalar tmp415 = w[4][l]*(A_12_2 + A_21_2);
                    const Scalar tmp416 = w[12][l]*(-A_12_0 - A_12_1 + A_21_4 + A_21_5);
                    const Scalar tmp417 = w[6][l]*(-A_12_6 - A_12_7 + A_21_2 + A_21_3);
                    const Scalar tmp418 = w[17][l]*(A_02_0 + A_20_5);
                    const Scalar tmp419 = w[2][l]*(-A_02_7 - A_20_2);
                    const Scalar tmp420 = w[18][l]*(-A_12_4 - A_21_2);
                    const Scalar tmp421 = w[10][l]*(A_12_2 + A_12_5 + A_21_3 + A_21_4);
                    const Scalar tmp422 = w[3][l]*(-A_11_1 - A_11_3 - A_11_5 - A_11_7);
                    const Scalar tmp423 = w[1][l]*(A_01_1 + A_01_5 - A_10_3 - A_10_7);
                    const Scalar tmp424 = w[25][l]*(-A_22_0 - A_22_2 - A_22_4 - A_22_6);
                    const Scalar tmp425 = w[4][l]*(A_12_3 + A_21_5);
                    const Scalar tmp426 = w[15][l]*(A_02_4 + A_02_6 - A_20_0 - A_20_2);
                    const Scalar tmp427 = w[16][l]*(A_01_2 + A_01_6 - A_10_0 - A_10_4);
                    const Scalar tmp428 = w[9][l]*(-A_11_0 - A_11_2 - A_11_4 - A_11_6);
                    const Scalar tmp429 = w[5][l]*(A_02_1 + A_02_3 - A_20_5 - A_20_7);
                    const Scalar tmp430 = w[23][l]*(-A_22_1 - A_22_3 - A_22_5 - A_22_7);
                    const Scalar tmp431 = w[18][l]*(A_12_5 - A_21_7);
                    const Scalar tmp432 = w[10][l]*(-A_12_3 - A_12_4 + A_21_1 + A_21_6);
                    const Scalar tmp433 = w[21][l]*(A_01_7 - A_10_5);
                    const Scalar tmp434 = w[20][l]*(A_01_0 - A_10_2);
                    const Scalar tmp435 = w[4][l]*(-A_12_2 + A_21_0);
                    const Scalar tmp436 = w[8][l]*(-A_01_3 - A_01_4 + A_10_1 + A_10_6);
                    const Scalar tmp437 = w[2][l]*(-A_02_4 + A_20_5);
                    const Scalar tmp438 = w[20][l]*(A_01_4 - A_10_5);
                    const Scalar tmp439 = w[21][l]*(A_01_3 - A_10_2);
                    const Scalar tmp440 = w[16][l]*(-A_01_1 - A_01_2 + A_10_0 + A_10_3);
                    const Scalar tmp441 = w[1][l]*(-A_01_5 - A_01_6 + A_10_4 + A_10_7);
                    const Scalar tmp442 = w[17][l]*(A_02_3 - A_20_2);
                    const Scalar tmp443 = w[20][l]*(-A_01_4 - A_10_7);
                    const Scalar tmp444 = w[21][l]*(-A_01_3 - A_10_0);
                    const Scalar tmp445 = w[18][l]*(A_12_6 + A_21_0);
                    const Scalar tmp446 = w[10][l]*(-A_12_0 - A_12_7 - A_21_1 - A_21_6);
                    const Scalar tmp447 = w[1][l]*(-A_01_3 - A_01_7 + A_10_1 + A_10_5);
                    const Scalar tmp448 = w[4][l]*(-A_12_1 - A_21_7);
                    const Scalar tmp449 = w[16][l]*(-A_01_0 - A_01_4 + A_10_2 + A_10_6);
                    const Scalar tmp450 = w[2][l]*(A_02_7 - A_20_6);
                    const Scalar tmp451 = w[6][l]*(A_12_6 + A_12_7 + A_21_6 + A_21_7);
                    const Scalar tmp452 = w[20][l]*(A_01_7 - A_10_6);
                    const Scalar tmp453 = w[21][l]*(A_01_0 - A_10_1);
                    const Scalar tmp454 = w[12][l]*(A_12_0 + A_12_1 + A_21_0 + A_21_1);
                    const Scalar tmp455 = w[29][l]*(A_00_0 + A_00_1);
                    const Scalar tmp456 = w[28][l]*(A_00_6 + A_00_7);
                    const Scalar tmp457 = w[17][l]*(-A_02_0 + A_20_1);
                    const Scalar tmp458 = w[21][l]*(-A_01_7 - A_10_4);
                    const Scalar tmp459 = w[20][l]*(-A_01_0 - A_10_3);
                    const Scalar tmp460 = w[12][l]*(A_12_4 + A_12_5 - A_21_6 - A_21_7);
                    const Scalar tmp461 = w[6][l]*(A_12_2 + A_12_3 - A_21_0 - A_21_1);
                    const Scalar tmp462 = w[18][l]*(A_12_1 + A_21_7);
                    const Scalar tmp463 = w[4][l]*(-A_12_6 - A_21_0);
                    const Scalar tmp464 = w[15][l]*(A_02_1 + A_02_3 - A_20_5 - A_20_7);
                    const Scalar tmp465 = w[5][l]*(A_02_4 + A_02_6 - A_20_0 - A_20_2);
                    const Scalar tmp466 = w[2][l]*(-A_02_6 + A_20_7);
                    const Scalar tmp467 = w[20][l]*(-A_01_6 + A_10_7);
                    const Scalar tmp468 = w[21][l]*(-A_01_1 + A_10_0);
                    const Scalar tmp469 = w[17][l]*(A_02_1 - A_20_0);
                    const Scalar tmp470 = w[6][l]*(-A_12_2 - A_12_3 - A_21_4 - A_21_5);
                    const Scalar tmp471 = w[1][l]*(-A_01_1 - A_01_5 - A_10_2 - A_10_6);
                    const Scalar tmp472 = w[15][l]*(-A_02_4 - A_02_6 - A_20_1 - A_20_3);
                    const Scalar tmp473 = w[16][l]*(-A_01_2 - A_01_6 - A_10_1 - A_10_5);
                    const Scalar tmp474 = w[12][l]*(-A_12_4 - A_12_5 - A_21_2 - A_21_3);
                    const Scalar tmp475 = w[5][l]*(-A_02_1 - A_02_3 - A_20_4 - A_20_6);
                    const Scalar tmp476 = w[18][l]*(-A_12_6 + A_21_4);
                    const Scalar tmp477 = w[20][l]*(A_01_3 - A_10_1);
                    const Scalar tmp478 = w[10][l]*(A_12_0 + A_12_7 - A_21_2 - A_21_5);
                    const Scalar tmp479 = w[4][l]*(A_12_1 - A_21_3);
                    const Scalar tmp480 = w[21][l]*(A_01_4 - A_10_6);
                    const Scalar tmp481 = w[8][l]*(-A_01_0 - A_01_7 + A_10_2 + A_10_5);
                    const Scalar tmp482 = w[6][l]*(A_12_0 + A_12_1 + A_21_6 + A_21_7);
                    const Scalar tmp483 = w[12][l]*(A_12_6 + A_12_7 + A_21_0 + A_21_1);
                    const Scalar tmp484 = w[15][l]*(A_02_5 + A_02_7 + A_20_0 + A_20_2);
                    const Scalar tmp485 = w[5][l]*(A_02_0 + A_02_2 + A_20_5 + A_20_7);
                    const Scalar tmp486 = w[18][l]*(-A_12_1 + A_21_3);
                    const Scalar tmp487 = w[20][l]*(A_01_4 - A_10_6);
                    const Scalar tmp488 = w[4][l]*(A_12_6 - A_21_4);
                    const Scalar tmp489 = w[21][l]*(A_01_3 - A_10_1);
                    const Scalar tmp490 = w[20][l]*(A_01_7 - A_10_5);
                    const Scalar tmp491 = w[18][l]*(A_12_2 - A_21_0);
                    const Scalar tmp492 = w[4][l]*(-A_12_5 + A_21_7);
                    const Scalar tmp493 = w[21][l]*(A_01_0 - A_10_2);
                    const Scalar tmp494 = w[20][l]*(A_01_1 + A_10_2);
                    const Scalar tmp495 = w[21][l]*(A_01_6 + A_10_5);
                    const Scalar tmp496 = w[18][l]*(-A_12_2 - A_21_4);
                    const Scalar tmp497 = w[4][l]*(A_12_5 + A_21_3);
                    const Scalar tmp498 = w[15][l]*(-A_02_0 - A_02_2 + A_20_4 + A_20_6);
                    const Scalar tmp499 = w[5][l]*(-A_02_5 - A_02_7 + A_20_1 + A_20_3);
                    const Scalar tmp500 = w[18][l]*(-A_12_6 + A_21_2);
                    const Scalar tmp501 = w[4][l]*(A_12_1 - A_21_5);
                    const Scalar tmp502 = w[17][l]*(A_02_6 - A_20_2);
                    const Scalar tmp503 = w[2][l]*(-A_02_1 + A_20_5);
                    const Scalar tmp504 = w[18][l]*(-A_12_3 - A_21_5);
                    const Scalar tmp505 = w[4][l]*(A_12_4 + A_21_2);
                    const Scalar tmp506 = w[2][l]*(A_02_6 + A_20_3);
                    const Scalar tmp507 = w[17][l]*(-A_02_1 - A_20_4);
                    const Scalar tmp508 = w[18][l]*(A_12_0 + A_21_6);
                    const Scalar tmp509 = w[4][l]*(-A_12_7 - A_21_1);
                    EM_S[INDEX2(0,0,8)]+=tmp198 + tmp200 + tmp214 + tmp259 + tmp262 + tmp289 + tmp294 + tmp299 + tmp303 + tmp304 + tmp307 + tmp309 + tmp343 + tmp347 + tmp362 + tmp363 + tmp379 + tmp380 + tmp381 + tmp382 + tmp383 + tmp384 + tmp385 + tmp386;
                    EM_S[INDEX2(1,0,8)]+=tmp145 + tmp148 + tmp161 + tmp201 + tmp202 + tmp210 + tmp371 + tmp374 + tmp440 + tmp441 + tmp450 + tmp451 + tmp452 + tmp453 + tmp454 + tmp455 + tmp456 + tmp457 + tmp89 + tmp91;
                    EM_S[INDEX2(2,0,8)]+=tmp135 + tmp234 + tmp235 + tmp236 + tmp237 + tmp238 + tmp239 + tmp240 + tmp241 + tmp242 + tmp243 + tmp244 + tmp245 + tmp246 + tmp39 + tmp41 + tmp44 + tmp49 + tmp61 + tmp71;
                    EM_S[INDEX2(3,0,8)]+=tmp20 + tmp21 + tmp24 + tmp34 + tmp72 + tmp73 + tmp74 + tmp75 + tmp76 + tmp77 + tmp78 + tmp79 + tmp80 + tmp81 + tmp82 + tmp83;
                    EM_S[INDEX2(4,0,8)]+=tmp1 + tmp127 + tmp131 + tmp141 + tmp146 + tmp15 + tmp153 + tmp154 + tmp188 + tmp189 + tmp190 + tmp191 + tmp192 + tmp193 + tmp194 + tmp195 + tmp68 + tmp70 + tmp92 + tmp98;
                    EM_S[INDEX2(5,0,8)]+=tmp166 + tmp176 + tmp260 + tmp267 + tmp274 + tmp339 + tmp340 + tmp342 + tmp344 + tmp346 + tmp348 + tmp365 + tmp393 + tmp394 + tmp395 + tmp396;
                    EM_S[INDEX2(6,0,8)]+=tmp114 + tmp184 + tmp258 + tmp266 + tmp337 + tmp420 + tmp421 + tmp422 + tmp423 + tmp424 + tmp425 + tmp426 + tmp427 + tmp428 + tmp429 + tmp430;
                    EM_S[INDEX2(7,0,8)]+=tmp104 + tmp106 + tmp112 + tmp113 + tmp38 + tmp470 + tmp471 + tmp472 + tmp473 + tmp474 + tmp475 + tmp87;
                    EM_S[INDEX2(0,1,8)]+=tmp161 + tmp201 + tmp247 + tmp250 + tmp371 + tmp374 + tmp44 + tmp451 + tmp454 + tmp455 + tmp456 + tmp466 + tmp467 + tmp468 + tmp469 + tmp49 + tmp89 + tmp91 + tmp92 + tmp98;
                    EM_S[INDEX2(1,1,8)]+=tmp215 + tmp221 + tmp227 + tmp260 + tmp267 + tmp288 + tmp304 + tmp312 + tmp317 + tmp351 + tmp352 + tmp353 + tmp354 + tmp355 + tmp356 + tmp357 + tmp358 + tmp359 + tmp360 + tmp361 + tmp362 + tmp363 + tmp76 + tmp79;
                    EM_S[INDEX2(2,1,8)]+=tmp114 + tmp120 + tmp167 + tmp170 + tmp198 + tmp20 + tmp200 + tmp24 + tmp443 + tmp444 + tmp73 + tmp74 + tmp75 + tmp80 + tmp81 + tmp83;
                    EM_S[INDEX2(3,1,8)]+=tmp13 + tmp16 + tmp38 + tmp39 + tmp40 + tmp41 + tmp43 + tmp440 + tmp441 + tmp45 + tmp47 + tmp478 + tmp481 + tmp486 + tmp487 + tmp488 + tmp489 + tmp50 + tmp52 + tmp55;
                    EM_S[INDEX2(4,1,8)]+=tmp166 + tmp176 + tmp184 + tmp283 + tmp339 + tmp340 + tmp341 + tmp342 + tmp343 + tmp344 + tmp345 + tmp346 + tmp347 + tmp348 + tmp349 + tmp350;
                    EM_S[INDEX2(5,1,8)]+=tmp112 + tmp12 + tmp123 + tmp124 + tmp126 + tmp140 + tmp141 + tmp142 + tmp143 + tmp144 + tmp145 + tmp146 + tmp147 + tmp148 + tmp149 + tmp150 + tmp151 + tmp51 + tmp54 + tmp6;
                    EM_S[INDEX2(6,1,8)]+=tmp104 + tmp106 + tmp113 + tmp135 + tmp15 + tmp162 + tmp163 + tmp470 + tmp474 + tmp484 + tmp485 + tmp87;
                    EM_S[INDEX2(7,1,8)]+=tmp21 + tmp225 + tmp232 + tmp274 + tmp329 + tmp332 + tmp333 + tmp335 + tmp336 + tmp337 + tmp338 + tmp397 + tmp398 + tmp399 + tmp400 + tmp401;
                    EM_S[INDEX2(0,2,8)]+=tmp135 + tmp236 + tmp238 + tmp240 + tmp242 + tmp244 + tmp39 + tmp41 + tmp432 + tmp436 + tmp440 + tmp441 + tmp490 + tmp491 + tmp492 + tmp493 + tmp61 + tmp68 + tmp70 + tmp71;
                    EM_S[INDEX2(1,2,8)]+=tmp166 + tmp169 + tmp172 + tmp196 + tmp197 + tmp198 + tmp199 + tmp20 + tmp200 + tmp21 + tmp73 + tmp74 + tmp75 + tmp77 + tmp80 + tmp82;
                    EM_S[INDEX2(2,2,8)]+=tmp217 + tmp231 + tmp233 + tmp258 + tmp266 + tmp271 + tmp273 + tmp288 + tmp289 + tmp290 + tmp291 + tmp292 + tmp293 + tmp294 + tmp295 + tmp296 + tmp297 + tmp298 + tmp299 + tmp300 + tmp301 + tmp302 + tmp76 + tmp79;
                    EM_S[INDEX2(3,2,8)]+=tmp101 + tmp14 + tmp204 + tmp205 + tmp367 + tmp368 + tmp369 + tmp370 + tmp371 + tmp372 + tmp373 + tmp374 + tmp375 + tmp376 + tmp377 + tmp378 + tmp44 + tmp49 + tmp87 + tmp9;
                    EM_S[INDEX2(4,2,8)]+=tmp114 + tmp274 + tmp337 + tmp383 + tmp386 + tmp422 + tmp424 + tmp426 + tmp428 + tmp429 + tmp430 + tmp445 + tmp446 + tmp447 + tmp448 + tmp449;
                    EM_S[INDEX2(5,2,8)]+=tmp104 + tmp106 + tmp113 + tmp136 + tmp138 + tmp15 + tmp161 + tmp38 + tmp472 + tmp475 + tmp482 + tmp483;
                    EM_S[INDEX2(6,2,8)]+=tmp10 + tmp112 + tmp123 + tmp125 + tmp127 + tmp128 + tmp130 + tmp131 + tmp132 + tmp156 + tmp157 + tmp243 + tmp246 + tmp278 + tmp279 + tmp3 + tmp500 + tmp501 + tmp502 + tmp503;
                    EM_S[INDEX2(7,2,8)]+=tmp173 + tmp174 + tmp175 + tmp176 + tmp177 + tmp178 + tmp179 + tmp180 + tmp181 + tmp182 + tmp183 + tmp184 + tmp185 + tmp186 + tmp187 + tmp24;
                    EM_S[INDEX2(0,3,8)]+=tmp114 + tmp165 + tmp166 + tmp167 + tmp168 + tmp169 + tmp170 + tmp171 + tmp172 + tmp20 + tmp73 + tmp74 + tmp75 + tmp76 + tmp79 + tmp80;
                    EM_S[INDEX2(1,3,8)]+=tmp36 + tmp37 + tmp38 + tmp39 + tmp40 + tmp41 + tmp42 + tmp43 + tmp44 + tmp45 + tmp46 + tmp47 + tmp48 + tmp49 + tmp50 + tmp51 + tmp52 + tmp53 + tmp54 + tmp55;
                    EM_S[INDEX2(2,3,8)]+=tmp101 + tmp156 + tmp157 + tmp204 + tmp205 + tmp368 + tmp371 + tmp372 + tmp374 + tmp375 + tmp377 + tmp437 + tmp438 + tmp439 + tmp440 + tmp441 + tmp442 + tmp85 + tmp87 + tmp99;
                    EM_S[INDEX2(3,3,8)]+=tmp179 + tmp183 + tmp198 + tmp200 + tmp214 + tmp215 + tmp216 + tmp217 + tmp218 + tmp219 + tmp220 + tmp221 + tmp222 + tmp223 + tmp224 + tmp225 + tmp226 + tmp227 + tmp228 + tmp229 + tmp230 + tmp231 + tmp232 + tmp233;
                    EM_S[INDEX2(4,3,8)]+=tmp104 + tmp106 + tmp107 + tmp109 + tmp112 + tmp113 + tmp135 + tmp161 + tmp482 + tmp483 + tmp484 + tmp485;
                    EM_S[INDEX2(5,3,8)]+=tmp184 + tmp21 + tmp312 + tmp317 + tmp327 + tmp328 + tmp329 + tmp330 + tmp331 + tmp332 + tmp333 + tmp334 + tmp335 + tmp336 + tmp337 + tmp338;
                    EM_S[INDEX2(6,3,8)]+=tmp175 + tmp176 + tmp178 + tmp180 + tmp182 + tmp185 + tmp187 + tmp24 + tmp269 + tmp270 + tmp271 + tmp272 + tmp273 + tmp274 + tmp275 + tmp276;
                    EM_S[INDEX2(7,3,8)]+=tmp0 + tmp1 + tmp10 + tmp11 + tmp12 + tmp13 + tmp14 + tmp15 + tmp16 + tmp17 + tmp18 + tmp19 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6 + tmp7 + tmp8 + tmp9;
                    EM_S[INDEX2(0,4,8)]+=tmp1 + tmp127 + tmp131 + tmp141 + tmp145 + tmp146 + tmp148 + tmp15 + tmp189 + tmp190 + tmp192 + tmp193 + tmp2 + tmp243 + tmp246 + tmp406 + tmp407 + tmp408 + tmp409 + tmp5;
                    EM_S[INDEX2(1,4,8)]+=tmp176 + tmp24 + tmp269 + tmp274 + tmp339 + tmp340 + tmp342 + tmp343 + tmp344 + tmp347 + tmp394 + tmp395 + tmp416 + tmp417 + tmp418 + tmp419;
                    EM_S[INDEX2(2,4,8)]+=tmp184 + tmp21 + tmp328 + tmp337 + tmp383 + tmp386 + tmp422 + tmp423 + tmp424 + tmp427 + tmp428 + tmp430 + tmp498 + tmp499 + tmp508 + tmp509;
                    EM_S[INDEX2(3,4,8)]+=tmp104 + tmp106 + tmp112 + tmp113 + tmp135 + tmp137 + tmp139 + tmp160 + tmp161 + tmp164 + tmp471 + tmp473;
                    EM_S[INDEX2(4,4,8)]+=tmp118 + tmp121 + tmp214 + tmp215 + tmp216 + tmp217 + tmp220 + tmp222 + tmp253 + tmp254 + tmp255 + tmp256 + tmp257 + tmp258 + tmp259 + tmp260 + tmp261 + tmp262 + tmp263 + tmp264 + tmp265 + tmp266 + tmp267 + tmp268;
                    EM_S[INDEX2(5,4,8)]+=tmp100 + tmp101 + tmp102 + tmp103 + tmp84 + tmp85 + tmp86 + tmp87 + tmp88 + tmp89 + tmp90 + tmp91 + tmp92 + tmp93 + tmp94 + tmp95 + tmp96 + tmp97 + tmp98 + tmp99;
                    EM_S[INDEX2(6,4,8)]+=tmp38 + tmp42 + tmp43 + tmp53 + tmp56 + tmp57 + tmp58 + tmp59 + tmp60 + tmp61 + tmp62 + tmp63 + tmp64 + tmp65 + tmp66 + tmp67 + tmp68 + tmp69 + tmp70 + tmp71;
                    EM_S[INDEX2(7,4,8)]+=tmp114 + tmp117 + tmp119 + tmp166 + tmp171 + tmp20 + tmp22 + tmp25 + tmp26 + tmp28 + tmp29 + tmp30 + tmp460 + tmp461 + tmp494 + tmp495;
                    EM_S[INDEX2(0,5,8)]+=tmp174 + tmp176 + tmp184 + tmp24 + tmp260 + tmp267 + tmp339 + tmp340 + tmp341 + tmp342 + tmp344 + tmp345 + tmp416 + tmp417 + tmp506 + tmp507;
                    EM_S[INDEX2(1,5,8)]+=tmp112 + tmp12 + tmp123 + tmp13 + tmp141 + tmp142 + tmp143 + tmp146 + tmp147 + tmp149 + tmp16 + tmp277 + tmp278 + tmp279 + tmp280 + tmp281 + tmp282 + tmp6 + tmp92 + tmp98;
                    EM_S[INDEX2(2,5,8)]+=tmp104 + tmp106 + tmp108 + tmp111 + tmp113 + tmp15 + tmp160 + tmp161 + tmp162 + tmp163 + tmp164 + tmp38;
                    EM_S[INDEX2(3,5,8)]+=tmp114 + tmp274 + tmp312 + tmp317 + tmp329 + tmp332 + tmp335 + tmp337 + tmp338 + tmp399 + tmp401 + tmp446 + tmp462 + tmp463 + tmp464 + tmp465;
                    EM_S[INDEX2(4,5,8)]+=tmp100 + tmp101 + tmp145 + tmp148 + tmp369 + tmp376 + tmp402 + tmp403 + tmp404 + tmp405 + tmp60 + tmp65 + tmp84 + tmp87 + tmp88 + tmp89 + tmp91 + tmp95 + tmp96 + tmp97;
                    EM_S[INDEX2(5,5,8)]+=tmp217 + tmp225 + tmp232 + tmp26 + tmp265 + tmp268 + tmp288 + tmp289 + tmp29 + tmp290 + tmp293 + tmp295 + tmp308 + tmp313 + tmp343 + tmp347 + tmp358 + tmp359 + tmp410 + tmp411 + tmp412 + tmp413 + tmp414 + tmp415;
                    EM_S[INDEX2(6,5,8)]+=tmp118 + tmp121 + tmp166 + tmp199 + tmp20 + tmp21 + tmp22 + tmp25 + tmp27 + tmp28 + tmp30 + tmp33 + tmp458 + tmp459 + tmp460 + tmp461;
                    EM_S[INDEX2(7,5,8)]+=tmp135 + tmp238 + tmp320 + tmp324 + tmp325 + tmp326 + tmp431 + tmp432 + tmp433 + tmp434 + tmp435 + tmp436 + tmp45 + tmp51 + tmp54 + tmp55 + tmp57 + tmp64 + tmp90 + tmp94;
                    EM_S[INDEX2(0,6,8)]+=tmp21 + tmp258 + tmp266 + tmp274 + tmp337 + tmp398 + tmp422 + tmp424 + tmp428 + tmp430 + tmp447 + tmp449 + tmp496 + tmp497 + tmp498 + tmp499;
                    EM_S[INDEX2(1,6,8)]+=tmp104 + tmp105 + tmp106 + tmp110 + tmp113 + tmp135 + tmp136 + tmp137 + tmp138 + tmp139 + tmp15 + tmp87;
                    EM_S[INDEX2(2,6,8)]+=tmp10 + tmp112 + tmp122 + tmp123 + tmp124 + tmp125 + tmp126 + tmp127 + tmp128 + tmp129 + tmp130 + tmp131 + tmp132 + tmp133 + tmp134 + tmp14 + tmp3 + tmp68 + tmp70 + tmp9;
                    EM_S[INDEX2(3,6,8)]+=tmp166 + tmp175 + tmp176 + tmp177 + tmp178 + tmp180 + tmp181 + tmp184 + tmp187 + tmp271 + tmp273 + tmp283 + tmp284 + tmp285 + tmp286 + tmp287;
                    EM_S[INDEX2(4,6,8)]+=tmp243 + tmp246 + tmp38 + tmp43 + tmp476 + tmp477 + tmp478 + tmp479 + tmp480 + tmp481 + tmp57 + tmp58 + tmp61 + tmp63 + tmp64 + tmp66 + tmp69 + tmp71 + tmp90 + tmp94;
                    EM_S[INDEX2(5,6,8)]+=tmp114 + tmp115 + tmp116 + tmp117 + tmp118 + tmp119 + tmp120 + tmp121 + tmp20 + tmp22 + tmp24 + tmp25 + tmp28 + tmp30 + tmp32 + tmp35;
                    EM_S[INDEX2(6,6,8)]+=tmp179 + tmp183 + tmp215 + tmp255 + tmp26 + tmp261 + tmp288 + tmp29 + tmp298 + tmp300 + tmp304 + tmp316 + tmp318 + tmp351 + tmp354 + tmp355 + tmp383 + tmp386 + tmp387 + tmp388 + tmp389 + tmp390 + tmp391 + tmp392;
                    EM_S[INDEX2(7,6,8)]+=tmp100 + tmp156 + tmp157 + tmp161 + tmp201 + tmp204 + tmp205 + tmp207 + tmp208 + tmp209 + tmp211 + tmp247 + tmp248 + tmp249 + tmp250 + tmp251 + tmp252 + tmp60 + tmp65 + tmp88;
                    EM_S[INDEX2(0,7,8)]+=tmp104 + tmp105 + tmp106 + tmp107 + tmp108 + tmp109 + tmp110 + tmp111 + tmp112 + tmp113 + tmp38 + tmp87;
                    EM_S[INDEX2(1,7,8)]+=tmp114 + tmp184 + tmp225 + tmp232 + tmp329 + tmp330 + tmp332 + tmp334 + tmp335 + tmp337 + tmp338 + tmp421 + tmp464 + tmp465 + tmp504 + tmp505;
                    EM_S[INDEX2(2,7,8)]+=tmp166 + tmp175 + tmp176 + tmp178 + tmp179 + tmp180 + tmp183 + tmp187 + tmp270 + tmp272 + tmp274 + tmp284 + tmp285 + tmp364 + tmp365 + tmp366;
                    EM_S[INDEX2(3,7,8)]+=tmp1 + tmp10 + tmp11 + tmp12 + tmp15 + tmp152 + tmp153 + tmp154 + tmp155 + tmp156 + tmp157 + tmp158 + tmp159 + tmp17 + tmp3 + tmp4 + tmp51 + tmp54 + tmp6 + tmp7;
                    EM_S[INDEX2(4,7,8)]+=tmp20 + tmp21 + tmp22 + tmp23 + tmp24 + tmp25 + tmp26 + tmp27 + tmp28 + tmp29 + tmp30 + tmp31 + tmp32 + tmp33 + tmp34 + tmp35;
                    EM_S[INDEX2(5,7,8)]+=tmp13 + tmp135 + tmp16 + tmp237 + tmp238 + tmp245 + tmp319 + tmp320 + tmp321 + tmp322 + tmp323 + tmp324 + tmp325 + tmp326 + tmp45 + tmp55 + tmp57 + tmp60 + tmp64 + tmp65;
                    EM_S[INDEX2(6,7,8)]+=tmp100 + tmp14 + tmp161 + tmp201 + tmp202 + tmp203 + tmp204 + tmp205 + tmp206 + tmp207 + tmp208 + tmp209 + tmp210 + tmp211 + tmp212 + tmp213 + tmp88 + tmp9 + tmp90 + tmp94;
                    EM_S[INDEX2(7,7,8)]+=tmp118 + tmp121 + tmp214 + tmp226 + tmp228 + tmp271 + tmp273 + tmp289 + tmp303 + tmp304 + tmp305 + tmp306 + tmp307 + tmp308 + tmp309 + tmp310 + tmp311 + tmp312 + tmp313 + tmp314 + tmp315 + tmp316 + tmp317 + tmp318;
                } else { // constant data
                    const Scalar Aw00 =  8.*A_p[INDEX2(0,0,3)]*w[27][l];
                    const Scalar Aw01 = 12.*A_p[INDEX2(0,1,3)]*w[8][l];
                    const Scalar Aw02 = 12.*A_p[INDEX2(0,2,3)]*w[11][l];
                    const Scalar Aw10 = 12.*A_p[INDEX2(1,0,3)]*w[8][l];
                    const Scalar Aw11 =  8.*A_p[INDEX2(1,1,3)]*w[22][l];
                    const Scalar Aw12 = 12.*A_p[INDEX2(1,2,3)]*w[10][l];
                    const Scalar Aw20 = 12.*A_p[INDEX2(2,0,3)]*w[11][l];
                    const Scalar Aw21 = 12.*A_p[INDEX2(2,1,3)]*w[10][l];
                    const Scalar Aw22 =  8.*A_p[INDEX2(2,2,3)]*w[13][l];
                    const Scalar tmp0 = Aw01 + Aw10;
                    const Scalar tmp1 = Aw01 - Aw10;
                    const Scalar tmp2 = -Aw01 - Aw10;
                    const Scalar tmp3 = -Aw01 + Aw10;
                    const Scalar tmp4 = Aw02 + Aw20;
                    const Scalar tmp5 = Aw02 - Aw20;
                    const Scalar tmp6 = -Aw02 - Aw20;
                    const Scalar tmp7 = -Aw02 + Aw20;
                    const Scalar tmp8 = Aw12 + Aw21;
                    const Scalar tmp9 = Aw12 - Aw21;
                    const Scalar tmp10 = -Aw12 - Aw21;
                    const Scalar tmp11 = -Aw12 + Aw21;
                    EM_S[INDEX2(0,0,8)]+=-4.*Aw00 - 4.*Aw11 - 4.*Aw22 + 2.*tmp0 + 2.*tmp4 + 2.*tmp10;
                    EM_S[INDEX2(1,0,8)]+= 4.*Aw00 - 2.*Aw11 - 2.*Aw22 + 2.*tmp7 + 2.*tmp3 + tmp10;
                    EM_S[INDEX2(2,0,8)]+=-2.*Aw00 + 4.*Aw11 - 2.*Aw22 + 2.*tmp9 + tmp4 + 2.*tmp1;
                    EM_S[INDEX2(3,0,8)]+= 2.*Aw00 + 2.*Aw11 -    Aw22 + tmp7 + tmp9 + 2.*tmp2;
                    EM_S[INDEX2(4,0,8)]+=-2.*Aw00 - 2.*Aw11 + 4.*Aw22 + tmp0 + 2.*tmp11 + 2.*tmp5;
                    EM_S[INDEX2(5,0,8)]+= 2.*Aw00 -    Aw11 + 2.*Aw22 + tmp11 + 2.*tmp6 + tmp3;
                    EM_S[INDEX2(6,0,8)]+=   -Aw00 + 2.*Aw11 + 2.*Aw22 + tmp5 + tmp1 + 2.*tmp8;
                    EM_S[INDEX2(7,0,8)]+=    Aw00 +    Aw11 +    Aw22 + tmp8 + tmp2 + tmp6;
                    EM_S[INDEX2(0,1,8)]+= 4.*Aw00 - 2.*Aw11 - 2.*Aw22 + tmp10 + 2.*tmp1 + 2.*tmp5;
                    EM_S[INDEX2(1,1,8)]+=-4.*Aw00 - 4.*Aw11 - 4.*Aw22 + 2.*tmp6 + 2.*tmp10 + 2.*tmp2;
                    EM_S[INDEX2(2,1,8)]+= 2.*Aw00 + 2.*Aw11 -    Aw22 + tmp5 + 2.*tmp0 + tmp9;
                    EM_S[INDEX2(3,1,8)]+=-2.*Aw00 + 4.*Aw11 - 2.*Aw22 + 2.*tmp9 + 2.*tmp3 + tmp6;
                    EM_S[INDEX2(4,1,8)]+= 2.*Aw00 -    Aw11 + 2.*Aw22 + tmp11 + tmp1 + 2.*tmp4;
                    EM_S[INDEX2(5,1,8)]+=-2.*Aw00 - 2.*Aw11 + 4.*Aw22 + 2.*tmp7 + tmp2 + 2.*tmp11;
                    EM_S[INDEX2(6,1,8)]+=    Aw00 +    Aw11 +    Aw22 + tmp8 + tmp4 + tmp0;
                    EM_S[INDEX2(7,1,8)]+=   -Aw00 + 2.*Aw11 + 2.*Aw22 + 2.*tmp8 + tmp3 + tmp7;
                    EM_S[INDEX2(0,2,8)]+=-2.*Aw00 + 4.*Aw11 - 2.*Aw22 + 2.*tmp3 + tmp4 + 2.*tmp11;
                    EM_S[INDEX2(1,2,8)]+= 2.*Aw00 + 2.*Aw11 -    Aw22 + 2.*tmp0 + tmp11 + tmp7;
                    EM_S[INDEX2(2,2,8)]+=-4.*Aw00 - 4.*Aw11 - 4.*Aw22 + 2.*tmp8 + 2.*tmp4 + 2.*tmp2;
                    EM_S[INDEX2(3,2,8)]+= 4.*Aw00 - 2.*Aw11 - 2.*Aw22 + 2.*tmp7 + tmp8 + 2.*tmp1;
                    EM_S[INDEX2(4,2,8)]+=   -Aw00 + 2.*Aw11 + 2.*Aw22 + tmp5 + tmp3 + 2.*tmp10;
                    EM_S[INDEX2(5,2,8)]+=    Aw00 +    Aw11 +    Aw22 + tmp10 + tmp0 + tmp6;
                    EM_S[INDEX2(6,2,8)]+=-2.*Aw00 - 2.*Aw11 + 4.*Aw22 + 2.*tmp9 + tmp2 + 2.*tmp5;
                    EM_S[INDEX2(7,2,8)]+= 2.*Aw00 -    Aw11 + 2.*Aw22 + 2.*tmp6 + tmp1 + tmp9;
                    EM_S[INDEX2(0,3,8)]+= 2.*Aw00 + 2.*Aw11 -    Aw22 + tmp5 + tmp11 + 2.*tmp2;
                    EM_S[INDEX2(1,3,8)]+=-2.*Aw00 + 4.*Aw11 - 2.*Aw22 + tmp6 + 2.*tmp11 + 2.*tmp1;
                    EM_S[INDEX2(2,3,8)]+= 4.*Aw00 - 2.*Aw11 - 2.*Aw22 + tmp8 + 2.*tmp3 + 2.*tmp5;
                    EM_S[INDEX2(3,3,8)]+=-4.*Aw00 - 4.*Aw11 - 4.*Aw22 + 2.*tmp0 + 2.*tmp6 + 2.*tmp8;
                    EM_S[INDEX2(4,3,8)]+=    Aw00 +    Aw11 +    Aw22 + tmp2 + tmp4 + tmp10;
                    EM_S[INDEX2(5,3,8)]+=   -Aw00 + 2.*Aw11 + 2.*Aw22 + tmp1 + 2.*tmp10 + tmp7;
                    EM_S[INDEX2(6,3,8)]+= 2.*Aw00 -    Aw11 + 2.*Aw22 + 2.*tmp4 + tmp3 + tmp9;
                    EM_S[INDEX2(7,3,8)]+=-2.*Aw00 - 2.*Aw11 + 4.*Aw22 + 2.*tmp7 + 2.*tmp9 + tmp0;
                    EM_S[INDEX2(0,4,8)]+=-2.*Aw00 - 2.*Aw11 + 4.*Aw22 + 2.*tmp7 + 2.*tmp9 + tmp0;
                    EM_S[INDEX2(1,4,8)]+= 2.*Aw00 -    Aw11 + 2.*Aw22 + 2.*tmp4 + tmp3 + tmp9;
                    EM_S[INDEX2(2,4,8)]+=   -Aw00 + 2.*Aw11 + 2.*Aw22 + tmp1 + 2.*tmp10 + tmp7;
                    EM_S[INDEX2(3,4,8)]+=    Aw00 +    Aw11 +    Aw22 + tmp2 + tmp4 + tmp10;
                    EM_S[INDEX2(4,4,8)]+=-4.*Aw00 - 4.*Aw11 - 4.*Aw22 + 2.*tmp0 + 2.*tmp6 + 2.*tmp8;
                    EM_S[INDEX2(5,4,8)]+= 4.*Aw00 - 2.*Aw11 - 2.*Aw22 + tmp8 + 2.*tmp3 + 2.*tmp5;
                    EM_S[INDEX2(6,4,8)]+=-2.*Aw00 + 4.*Aw11 - 2.*Aw22 + tmp6 + 2.*tmp11 + 2.*tmp1;
                    EM_S[INDEX2(7,4,8)]+= 2.*Aw00 + 2.*Aw11 -    Aw22 + tmp5 + tmp11 + 2.*tmp2;
                    EM_S[INDEX2(0,5,8)]+= 2.*Aw00 -    Aw11 + 2.*Aw22 + 2.*tmp6 + tmp1 + tmp9;
                    EM_S[INDEX2(1,5,8)]+=-2.*Aw00 - 2.*Aw11 + 4.*Aw22 + 2.*tmp9 + tmp2 + 2.*tmp5;
                    EM_S[INDEX2(2,5,8)]+=    Aw00 +    Aw11 +    Aw22 + tmp10 + tmp0 + tmp6;
                    EM_S[INDEX2(3,5,8)]+=   -Aw00 + 2.*Aw11 + 2.*Aw22 + tmp5 + tmp3 + 2.*tmp10;
                    EM_S[INDEX2(4,5,8)]+= 4.*Aw00 - 2.*Aw11 - 2.*Aw22 + 2.*tmp7 + tmp8 + 2.*tmp1;
                    EM_S[INDEX2(5,5,8)]+=-4.*Aw00 - 4.*Aw11 - 4.*Aw22 + 2.*tmp8 + 2.*tmp4 + 2.*tmp2;
                    EM_S[INDEX2(6,5,8)]+= 2.*Aw00 + 2.*Aw11 -    Aw22 + 2.*tmp0 + tmp11 + tmp7;
                    EM_S[INDEX2(7,5,8)]+=-2.*Aw00 + 4.*Aw11 - 2.*Aw22 + 2.*tmp3 + tmp4 + 2.*tmp11;
                    EM_S[INDEX2(0,6,8)]+=   -Aw00 + 2.*Aw11 + 2.*Aw22 + 2.*tmp8 + tmp3 + tmp7;
                    EM_S[INDEX2(1,6,8)]+=    Aw00 +    Aw11 +    Aw22 + tmp8 + tmp4 + tmp0;
                    EM_S[INDEX2(2,6,8)]+=-2.*Aw00 - 2.*Aw11 + 4.*Aw22 + 2.*tmp7 + tmp2 + 2.*tmp11;
                    EM_S[INDEX2(3,6,8)]+= 2.*Aw00 -    Aw11 + 2.*Aw22 + tmp11 + tmp1 + 2.*tmp4;
                    EM_S[INDEX2(4,6,8)]+=-2.*Aw00 + 4.*Aw11 - 2.*Aw22 + 2.*tmp9 + 2.*tmp3 + tmp6;
                    EM_S[INDEX2(5,6,8)]+= 2.*Aw00 + 2.*Aw11 -    Aw22 + tmp5 + 2.*tmp0 + tmp9;
                    EM_S[INDEX2(6,6,8)]+=-4.*Aw00 - 4.*Aw11 - 4.*Aw22 + 2.*tmp6 + 2.*tmp10 + 2.*tmp2;
                    EM_S[INDEX2(7,6,8)]+= 4.*Aw00 - 2.*Aw11 - 2.*Aw22 + tmp10 + 2.*tmp1 + 2.*tmp5;
                    EM_S[INDEX2(0,7,8)]+=    Aw00 +    Aw11 +    Aw22 + tmp8 + tmp2 + tmp6;
                    EM_S[INDEX2(1,7,8)]+=   -Aw00 + 2.*Aw11 + 2.*Aw22 + tmp5 + tmp1 + 2.*tmp8;
                    EM_S[INDEX2(2,7,8)]+= 2.*Aw00 -    Aw11 + 2.*Aw22 + tmp11 + 2.*tmp6 + tmp3;
                    EM_S[INDEX2(3,7,8)]+=-2.*Aw00 - 2.*Aw11 + 4.*Aw22 + tmp0 + 2.*tmp11 + 2.*tmp5;
                    EM_S[INDEX2(4,7,8)]+= 2.*Aw00 + 2.*Aw11 -    Aw22 + tmp7 + tmp9 + 2.*tmp2;
                    EM_S[INDEX2(5,7,8)]+=-2.*Aw00 + 4.*Aw11 - 2.*Aw22 + 2.*tmp9 + tmp4 + 2.*tmp1;
                    EM_S[INDEX2(6,7,8)]+= 4.*Aw00 - 2.*Aw11 - 2.*Aw22 + 2.*tmp7 + 2.*tmp3 + tmp10;
                    EM_S[INDEX2(7,7,8)]+=-4.*Aw00 - 4.*Aw11 - 4.*Aw22 + 2.*tmp0 + 2.*tmp4 + 2.*tmp10;
                }
            }

            ///////////////
            // process B //
            ///////////////
            if (!B.isEmpty()) {
                const Scalar* B_p = B.getSampleDataRO(id, zero);
                if (B.actsExpanded()) {
                    const Scalar B_0_0 = B_p[INDEX2(0,0,3)];
                    const Scalar B_1_0 = B_p[INDEX2(1,0,3)];
                    const Scalar B_2_0 = B_p[INDEX2(2,0,3)];
                    const Scalar B_0_1 = B_p[INDEX2(0,1,3)];
                    const Scalar B_1_1 = B_p[INDEX2(1,1,3)];
                    const Scalar B_2_1 = B_p[INDEX2(2,1,3)];
                    const Scalar B_0_2 = B_p[INDEX2(0,2,3)];
                    const Scalar B_1_2 = B_p[INDEX2(1,2,3)];
                    const Scalar B_2_2 = B_p[INDEX2(2,2,3)];
                    const Scalar B_0_3 = B_p[INDEX2(0,3,3)];
                    const Scalar B_1_3 = B_p[INDEX2(1,3,3)];
                    const Scalar B_2_3 = B_p[INDEX2(2,3,3)];
                    const Scalar B_0_4 = B_p[INDEX2(0,4,3)];
                    const Scalar B_1_4 = B_p[INDEX2(1,4,3)];
                    const Scalar B_2_4 = B_p[INDEX2(2,4,3)];
                    const Scalar B_0_5 = B_p[INDEX2(0,5,3)];
                    const Scalar B_1_5 = B_p[INDEX2(1,5,3)];
                    const Scalar B_2_5 = B_p[INDEX2(2,5,3)];
                    const Scalar B_0_6 = B_p[INDEX2(0,6,3)];
                    const Scalar B_1_6 = B_p[INDEX2(1,6,3)];
                    const Scalar B_2_6 = B_p[INDEX2(2,6,3)];
                    const Scalar B_0_7 = B_p[INDEX2(0,7,3)];
                    const Scalar B_1_7 = B_p[INDEX2(1,7,3)];
                    const Scalar B_2_7 = B_p[INDEX2(2,7,3)];
                    const Scalar tmp0   = w[38][l]*(B_0_3 + B_0_7);
                    const Scalar tmp1   = w[31][l]*(B_1_0 + B_1_4);
                    const Scalar tmp2   = w[42][l]*(B_2_5 + B_2_6);
                    const Scalar tmp3   = w[35][l]*(B_2_1 + B_2_2);
                    const Scalar tmp4   = w[37][l]*(B_1_2 + B_1_6);
                    const Scalar tmp5   = w[39][l]*(B_1_3 + B_1_7);
                    const Scalar tmp6   = w[36][l]*(B_0_2 + B_0_6);
                    const Scalar tmp7   = w[33][l]*(B_0_1 + B_0_5);
                    const Scalar tmp8   = w[30][l]*(B_0_0 + B_0_4);
                    const Scalar tmp9   = w[34][l]*(B_1_1 + B_1_5);
                    const Scalar tmp10  = w[38][l]*(-B_0_5 - B_0_7);
                    const Scalar tmp11  = w[31][l]*(-B_1_0 - B_1_1);
                    const Scalar tmp12  = w[42][l]*(B_2_0 + B_2_1 + B_2_2 + B_2_3);
                    const Scalar tmp13  = w[35][l]*(B_2_4 + B_2_5 + B_2_6 + B_2_7);
                    const Scalar tmp14  = w[37][l]*(-B_1_2 - B_1_3);
                    const Scalar tmp15  = w[39][l]*(-B_1_6 - B_1_7);
                    const Scalar tmp16  = w[36][l]*(-B_0_4 - B_0_6);
                    const Scalar tmp17  = w[33][l]*(-B_0_1 - B_0_3);
                    const Scalar tmp18  = w[30][l]*(-B_0_0 - B_0_2);
                    const Scalar tmp19  = w[34][l]*(-B_1_4 - B_1_5);
                    const Scalar tmp20  = w[38][l]*(B_0_1 + B_0_3);
                    const Scalar tmp21  = w[42][l]*(-B_2_0 - B_2_2);
                    const Scalar tmp22  = w[35][l]*(-B_2_5 - B_2_7);
                    const Scalar tmp23  = w[37][l]*(-B_1_0 - B_1_5);
                    const Scalar tmp24  = w[32][l]*(-B_2_4 - B_2_6);
                    const Scalar tmp25  = w[36][l]*(B_0_0 + B_0_2);
                    const Scalar tmp26  = w[33][l]*(B_0_5 + B_0_7);
                    const Scalar tmp27  = w[30][l]*(B_0_4 + B_0_6);
                    const Scalar tmp28  = w[43][l]*(-B_2_1 - B_2_3);
                    const Scalar tmp29  = w[34][l]*(-B_1_2 - B_1_7);
                    const Scalar tmp30  = w[38][l]*(-B_0_4 - B_0_6);
                    const Scalar tmp31  = w[42][l]*(B_2_5 + B_2_7);
                    const Scalar tmp32  = w[35][l]*(B_2_0 + B_2_2);
                    const Scalar tmp33  = w[37][l]*(B_1_2 + B_1_7);
                    const Scalar tmp34  = w[32][l]*(B_2_1 + B_2_3);
                    const Scalar tmp35  = w[36][l]*(-B_0_5 - B_0_7);
                    const Scalar tmp36  = w[33][l]*(-B_0_0 - B_0_2);
                    const Scalar tmp37  = w[30][l]*(-B_0_1 - B_0_3);
                    const Scalar tmp38  = w[43][l]*(B_2_4 + B_2_6);
                    const Scalar tmp39  = w[34][l]*(B_1_0 + B_1_5);
                    const Scalar tmp40  = w[38][l]*(B_0_0 + B_0_2);
                    const Scalar tmp41  = w[31][l]*(B_1_6 + B_1_7);
                    const Scalar tmp42  = w[42][l]*(-B_2_4 - B_2_5 - B_2_6 - B_2_7);
                    const Scalar tmp43  = w[35][l]*(-B_2_0 - B_2_1 - B_2_2 - B_2_3);
                    const Scalar tmp44  = w[37][l]*(B_1_4 + B_1_5);
                    const Scalar tmp45  = w[39][l]*(B_1_0 + B_1_1);
                    const Scalar tmp46  = w[36][l]*(B_0_1 + B_0_3);
                    const Scalar tmp47  = w[33][l]*(B_0_4 + B_0_6);
                    const Scalar tmp48  = w[30][l]*(B_0_5 + B_0_7);
                    const Scalar tmp49  = w[34][l]*(B_1_2 + B_1_3);
                    const Scalar tmp50  = w[31][l]*(-B_1_2 - B_1_3);
                    const Scalar tmp51  = w[42][l]*(B_2_6 + B_2_7);
                    const Scalar tmp52  = w[35][l]*(B_2_0 + B_2_1);
                    const Scalar tmp53  = w[37][l]*(-B_1_0 - B_1_1);
                    const Scalar tmp54  = w[32][l]*(B_2_2 + B_2_3);
                    const Scalar tmp55  = w[39][l]*(-B_1_4 - B_1_5);
                    const Scalar tmp56  = w[36][l]*(B_0_0 + B_0_6);
                    const Scalar tmp57  = w[33][l]*(B_0_1 + B_0_7);
                    const Scalar tmp58  = w[43][l]*(B_2_4 + B_2_5);
                    const Scalar tmp59  = w[34][l]*(-B_1_6 - B_1_7);
                    const Scalar tmp60  = w[42][l]*(-B_2_0 - B_2_1 - B_2_2 - B_2_3);
                    const Scalar tmp61  = w[35][l]*(-B_2_4 - B_2_5 - B_2_6 - B_2_7);
                    const Scalar tmp62  = w[37][l]*(-B_1_0 - B_1_1 - B_1_4 - B_1_5);
                    const Scalar tmp63  = w[36][l]*(-B_0_1 - B_0_3 - B_0_5 - B_0_7);
                    const Scalar tmp64  = w[33][l]*(-B_0_0 - B_0_2 - B_0_4 - B_0_6);
                    const Scalar tmp65  = w[34][l]*(-B_1_2 - B_1_3 - B_1_6 - B_1_7);
                    const Scalar tmp66  = w[38][l]*(B_0_4 + B_0_6);
                    const Scalar tmp67  = w[36][l]*(B_0_5 + B_0_7);
                    const Scalar tmp68  = w[33][l]*(B_0_0 + B_0_2);
                    const Scalar tmp69  = w[30][l]*(B_0_1 + B_0_3);
                    const Scalar tmp70  = w[38][l]*(-B_0_2 - B_0_6);
                    const Scalar tmp71  = w[31][l]*(B_1_1 + B_1_5);
                    const Scalar tmp72  = w[42][l]*(-B_2_0 - B_2_3);
                    const Scalar tmp73  = w[35][l]*(-B_2_4 - B_2_7);
                    const Scalar tmp74  = w[37][l]*(B_1_3 + B_1_7);
                    const Scalar tmp75  = w[39][l]*(B_1_2 + B_1_6);
                    const Scalar tmp76  = w[36][l]*(-B_0_3 - B_0_7);
                    const Scalar tmp77  = w[33][l]*(-B_0_0 - B_0_4);
                    const Scalar tmp78  = w[30][l]*(-B_0_1 - B_0_5);
                    const Scalar tmp79  = w[34][l]*(B_1_0 + B_1_4);
                    const Scalar tmp80  = w[36][l]*(B_0_0 + B_0_2 + B_0_4 + B_0_6);
                    const Scalar tmp81  = w[33][l]*(B_0_1 + B_0_3 + B_0_5 + B_0_7);
                    const Scalar tmp82  = w[38][l]*(B_0_1 + B_0_5);
                    const Scalar tmp83  = w[31][l]*(-B_1_2 - B_1_6);
                    const Scalar tmp84  = w[42][l]*(B_2_4 + B_2_7);
                    const Scalar tmp85  = w[35][l]*(B_2_0 + B_2_3);
                    const Scalar tmp86  = w[37][l]*(-B_1_0 - B_1_4);
                    const Scalar tmp87  = w[39][l]*(-B_1_1 - B_1_5);
                    const Scalar tmp88  = w[36][l]*(B_0_0 + B_0_4);
                    const Scalar tmp89  = w[33][l]*(B_0_3 + B_0_7);
                    const Scalar tmp90  = w[30][l]*(B_0_2 + B_0_6);
                    const Scalar tmp91  = w[34][l]*(-B_1_3 - B_1_7);
                    const Scalar tmp92  = w[42][l]*(-B_2_1 - B_2_2);
                    const Scalar tmp93  = w[35][l]*(-B_2_5 - B_2_6);
                    const Scalar tmp94  = w[37][l]*(B_1_2 + B_1_3 + B_1_6 + B_1_7);
                    const Scalar tmp95  = w[34][l]*(B_1_0 + B_1_1 + B_1_4 + B_1_5);
                    const Scalar tmp96  = w[38][l]*(-B_0_1 - B_0_3);
                    const Scalar tmp97  = w[31][l]*(-B_1_4 - B_1_5);
                    const Scalar tmp98  = w[37][l]*(-B_1_6 - B_1_7);
                    const Scalar tmp99  = w[39][l]*(-B_1_2 - B_1_3);
                    const Scalar tmp100 = w[36][l]*(-B_0_0 - B_0_2);
                    const Scalar tmp101 = w[33][l]*(-B_0_5 - B_0_7);
                    const Scalar tmp102 = w[30][l]*(-B_0_4 - B_0_6);
                    const Scalar tmp103 = w[34][l]*(-B_1_0 - B_1_1);
                    const Scalar tmp104 = w[38][l]*(B_0_2 + B_0_6);
                    const Scalar tmp105 = w[42][l]*(B_2_0 + B_2_1);
                    const Scalar tmp106 = w[35][l]*(B_2_6 + B_2_7);
                    const Scalar tmp107 = w[37][l]*(B_1_0 + B_1_1 + B_1_4 + B_1_5);
                    const Scalar tmp108 = w[32][l]*(B_2_4 + B_2_5);
                    const Scalar tmp109 = w[36][l]*(B_0_3 + B_0_7);
                    const Scalar tmp110 = w[33][l]*(B_0_0 + B_0_4);
                    const Scalar tmp111 = w[30][l]*(B_0_1 + B_0_5);
                    const Scalar tmp112 = w[43][l]*(B_2_2 + B_2_3);
                    const Scalar tmp113 = w[34][l]*(B_1_2 + B_1_3 + B_1_6 + B_1_7);
                    const Scalar tmp114 = w[38][l]*(-B_0_0 - B_0_4);
                    const Scalar tmp115 = w[31][l]*(-B_1_3 - B_1_7);
                    const Scalar tmp116 = w[37][l]*(-B_1_1 - B_1_5);
                    const Scalar tmp117 = w[39][l]*(-B_1_0 - B_1_4);
                    const Scalar tmp118 = w[36][l]*(-B_0_1 - B_0_5);
                    const Scalar tmp119 = w[33][l]*(-B_0_2 - B_0_6);
                    const Scalar tmp120 = w[30][l]*(-B_0_3 - B_0_7);
                    const Scalar tmp121 = w[34][l]*(-B_1_2 - B_1_6);
                    const Scalar tmp122 = w[31][l]*(B_1_0 + B_1_1);
                    const Scalar tmp123 = w[42][l]*(B_2_4 + B_2_5);
                    const Scalar tmp124 = w[35][l]*(B_2_2 + B_2_3);
                    const Scalar tmp125 = w[37][l]*(B_1_2 + B_1_3);
                    const Scalar tmp126 = w[32][l]*(B_2_0 + B_2_1);
                    const Scalar tmp127 = w[39][l]*(B_1_6 + B_1_7);
                    const Scalar tmp128 = w[36][l]*(-B_0_3 - B_0_5);
                    const Scalar tmp129 = w[33][l]*(-B_0_2 - B_0_4);
                    const Scalar tmp130 = w[43][l]*(B_2_6 + B_2_7);
                    const Scalar tmp131 = w[34][l]*(B_1_4 + B_1_5);
                    const Scalar tmp132 = w[42][l]*(-B_2_5 - B_2_6);
                    const Scalar tmp133 = w[35][l]*(-B_2_1 - B_2_2);
                    const Scalar tmp134 = w[37][l]*(B_1_0 + B_1_5);
                    const Scalar tmp135 = w[36][l]*(B_0_1 + B_0_7);
                    const Scalar tmp136 = w[33][l]*(B_0_0 + B_0_6);
                    const Scalar tmp137 = w[34][l]*(B_1_2 + B_1_7);
                    const Scalar tmp138 = w[38][l]*(-B_0_0 - B_0_2);
                    const Scalar tmp139 = w[42][l]*(-B_2_1 - B_2_3);
                    const Scalar tmp140 = w[35][l]*(-B_2_4 - B_2_6);
                    const Scalar tmp141 = w[37][l]*(B_1_3 + B_1_6);
                    const Scalar tmp142 = w[32][l]*(-B_2_5 - B_2_7);
                    const Scalar tmp143 = w[36][l]*(-B_0_1 - B_0_3);
                    const Scalar tmp144 = w[33][l]*(-B_0_4 - B_0_6);
                    const Scalar tmp145 = w[30][l]*(-B_0_5 - B_0_7);
                    const Scalar tmp146 = w[43][l]*(-B_2_0 - B_2_2);
                    const Scalar tmp147 = w[34][l]*(B_1_1 + B_1_4);
                    const Scalar tmp148 = w[36][l]*(B_0_2 + B_0_4);
                    const Scalar tmp149 = w[33][l]*(B_0_3 + B_0_5);
                    const Scalar tmp150 = w[42][l]*(B_2_1 + B_2_2);
                    const Scalar tmp151 = w[35][l]*(B_2_5 + B_2_6);
                    const Scalar tmp152 = w[37][l]*(-B_1_2 - B_1_7);
                    const Scalar tmp153 = w[36][l]*(-B_0_0 - B_0_6);
                    const Scalar tmp154 = w[33][l]*(-B_0_1 - B_0_7);
                    const Scalar tmp155 = w[34][l]*(-B_1_0 - B_1_5);
                    const Scalar tmp156 = w[38][l]*(-B_0_3 - B_0_7);
                    const Scalar tmp157 = w[36][l]*(-B_0_2 - B_0_6);
                    const Scalar tmp158 = w[33][l]*(-B_0_1 - B_0_5);
                    const Scalar tmp159 = w[30][l]*(-B_0_0 - B_0_4);
                    const Scalar tmp160 = w[42][l]*(-B_2_4 - B_2_5);
                    const Scalar tmp161 = w[35][l]*(-B_2_2 - B_2_3);
                    const Scalar tmp162 = w[32][l]*(-B_2_0 - B_2_1);
                    const Scalar tmp163 = w[43][l]*(-B_2_6 - B_2_7);
                    const Scalar tmp164 = w[42][l]*(-B_2_4 - B_2_7);
                    const Scalar tmp165 = w[35][l]*(-B_2_0 - B_2_3);
                    const Scalar tmp166 = w[37][l]*(B_1_1 + B_1_4);
                    const Scalar tmp167 = w[34][l]*(B_1_3 + B_1_6);
                    const Scalar tmp168 = w[36][l]*(B_0_3 + B_0_5);
                    const Scalar tmp169 = w[33][l]*(B_0_2 + B_0_4);
                    const Scalar tmp170 = w[38][l]*(B_0_5 + B_0_7);
                    const Scalar tmp171 = w[42][l]*(B_2_4 + B_2_6);
                    const Scalar tmp172 = w[35][l]*(B_2_1 + B_2_3);
                    const Scalar tmp173 = w[37][l]*(-B_1_1 - B_1_4);
                    const Scalar tmp174 = w[32][l]*(B_2_0 + B_2_2);
                    const Scalar tmp175 = w[36][l]*(B_0_4 + B_0_6);
                    const Scalar tmp176 = w[33][l]*(B_0_1 + B_0_3);
                    const Scalar tmp177 = w[30][l]*(B_0_0 + B_0_2);
                    const Scalar tmp178 = w[43][l]*(B_2_5 + B_2_7);
                    const Scalar tmp179 = w[34][l]*(-B_1_3 - B_1_6);
                    const Scalar tmp180 = w[31][l]*(-B_1_0 - B_1_4);
                    const Scalar tmp181 = w[42][l]*(B_2_0 + B_2_2);
                    const Scalar tmp182 = w[35][l]*(B_2_5 + B_2_7);
                    const Scalar tmp183 = w[37][l]*(-B_1_2 - B_1_6);
                    const Scalar tmp184 = w[32][l]*(B_2_4 + B_2_6);
                    const Scalar tmp185 = w[39][l]*(-B_1_3 - B_1_7);
                    const Scalar tmp186 = w[36][l]*(B_0_1 + B_0_3 + B_0_5 + B_0_7);
                    const Scalar tmp187 = w[33][l]*(B_0_0 + B_0_2 + B_0_4 + B_0_6);
                    const Scalar tmp188 = w[43][l]*(B_2_1 + B_2_3);
                    const Scalar tmp189 = w[34][l]*(-B_1_1 - B_1_5);
                    const Scalar tmp190 = w[38][l]*(-B_0_1 - B_0_5);
                    const Scalar tmp191 = w[42][l]*(B_2_2 + B_2_3);
                    const Scalar tmp192 = w[35][l]*(B_2_4 + B_2_5);
                    const Scalar tmp193 = w[37][l]*(-B_1_2 - B_1_3 - B_1_6 - B_1_7);
                    const Scalar tmp194 = w[32][l]*(B_2_6 + B_2_7);
                    const Scalar tmp195 = w[36][l]*(-B_0_0 - B_0_4);
                    const Scalar tmp196 = w[33][l]*(-B_0_3 - B_0_7);
                    const Scalar tmp197 = w[30][l]*(-B_0_2 - B_0_6);
                    const Scalar tmp198 = w[43][l]*(B_2_0 + B_2_1);
                    const Scalar tmp199 = w[34][l]*(-B_1_0 - B_1_1 - B_1_4 - B_1_5);
                    const Scalar tmp200 = w[31][l]*(B_1_4 + B_1_5);
                    const Scalar tmp201 = w[42][l]*(-B_2_0 - B_2_1);
                    const Scalar tmp202 = w[35][l]*(-B_2_6 - B_2_7);
                    const Scalar tmp203 = w[37][l]*(B_1_6 + B_1_7);
                    const Scalar tmp204 = w[32][l]*(-B_2_4 - B_2_5);
                    const Scalar tmp205 = w[39][l]*(B_1_2 + B_1_3);
                    const Scalar tmp206 = w[43][l]*(-B_2_2 - B_2_3);
                    const Scalar tmp207 = w[34][l]*(B_1_0 + B_1_1);
                    const Scalar tmp208 = w[37][l]*(-B_1_3 - B_1_6);
                    const Scalar tmp209 = w[36][l]*(-B_0_2 - B_0_4);
                    const Scalar tmp210 = w[33][l]*(-B_0_3 - B_0_5);
                    const Scalar tmp211 = w[34][l]*(-B_1_1 - B_1_4);
                    const Scalar tmp212 = w[42][l]*(B_2_0 + B_2_3);
                    const Scalar tmp213 = w[35][l]*(B_2_4 + B_2_7);
                    const Scalar tmp214 = w[38][l]*(B_0_0 + B_0_4);
                    const Scalar tmp215 = w[36][l]*(B_0_1 + B_0_5);
                    const Scalar tmp216 = w[33][l]*(B_0_2 + B_0_6);
                    const Scalar tmp217 = w[30][l]*(B_0_3 + B_0_7);
                    const Scalar tmp218 = w[31][l]*(B_1_2 + B_1_6);
                    const Scalar tmp219 = w[37][l]*(B_1_0 + B_1_4);
                    const Scalar tmp220 = w[39][l]*(B_1_1 + B_1_5);
                    const Scalar tmp221 = w[34][l]*(B_1_3 + B_1_7);
                    const Scalar tmp222 = w[36][l]*(-B_0_1 - B_0_7);
                    const Scalar tmp223 = w[33][l]*(-B_0_0 - B_0_6);
                    const Scalar tmp224 = w[42][l]*(-B_2_6 - B_2_7);
                    const Scalar tmp225 = w[35][l]*(-B_2_0 - B_2_1);
                    const Scalar tmp226 = w[32][l]*(-B_2_2 - B_2_3);
                    const Scalar tmp227 = w[43][l]*(-B_2_4 - B_2_5);
                    const Scalar tmp228 = w[31][l]*(B_1_3 + B_1_7);
                    const Scalar tmp229 = w[42][l]*(B_2_1 + B_2_3);
                    const Scalar tmp230 = w[35][l]*(B_2_4 + B_2_6);
                    const Scalar tmp231 = w[37][l]*(B_1_1 + B_1_5);
                    const Scalar tmp232 = w[32][l]*(B_2_5 + B_2_7);
                    const Scalar tmp233 = w[39][l]*(B_1_0 + B_1_4);
                    const Scalar tmp234 = w[36][l]*(-B_0_0 - B_0_2 - B_0_4 - B_0_6);
                    const Scalar tmp235 = w[33][l]*(-B_0_1 - B_0_3 - B_0_5 - B_0_7);
                    const Scalar tmp236 = w[43][l]*(B_2_0 + B_2_2);
                    const Scalar tmp237 = w[34][l]*(B_1_2 + B_1_6);
                    const Scalar tmp238 = w[31][l]*(-B_1_1 - B_1_5);
                    const Scalar tmp239 = w[37][l]*(-B_1_3 - B_1_7);
                    const Scalar tmp240 = w[39][l]*(-B_1_2 - B_1_6);
                    const Scalar tmp241 = w[34][l]*(-B_1_0 - B_1_4);
                    const Scalar tmp242 = w[31][l]*(-B_1_6 - B_1_7);
                    const Scalar tmp243 = w[42][l]*(-B_2_2 - B_2_3);
                    const Scalar tmp244 = w[35][l]*(-B_2_4 - B_2_5);
                    const Scalar tmp245 = w[37][l]*(-B_1_4 - B_1_5);
                    const Scalar tmp246 = w[32][l]*(-B_2_6 - B_2_7);
                    const Scalar tmp247 = w[39][l]*(-B_1_0 - B_1_1);
                    const Scalar tmp248 = w[43][l]*(-B_2_0 - B_2_1);
                    const Scalar tmp249 = w[34][l]*(-B_1_2 - B_1_3);
                    const Scalar tmp250 = w[31][l]*(B_1_2 + B_1_3);
                    const Scalar tmp251 = w[37][l]*(B_1_0 + B_1_1);
                    const Scalar tmp252 = w[39][l]*(B_1_4 + B_1_5);
                    const Scalar tmp253 = w[34][l]*(B_1_6 + B_1_7);
                    const Scalar tmp254 = w[42][l]*(-B_2_4 - B_2_6);
                    const Scalar tmp255 = w[35][l]*(-B_2_1 - B_2_3);
                    const Scalar tmp256 = w[32][l]*(-B_2_0 - B_2_2);
                    const Scalar tmp257 = w[43][l]*(-B_2_5 - B_2_7);
                    const Scalar tmp258 = w[42][l]*(B_2_4 + B_2_5 + B_2_6 + B_2_7);
                    const Scalar tmp259 = w[35][l]*(B_2_0 + B_2_1 + B_2_2 + B_2_3);
                    const Scalar tmp260 = w[42][l]*(-B_2_5 - B_2_7);
                    const Scalar tmp261 = w[35][l]*(-B_2_0 - B_2_2);
                    const Scalar tmp262 = w[32][l]*(-B_2_1 - B_2_3);
                    const Scalar tmp263 = w[43][l]*(-B_2_4 - B_2_6);
                    EM_S[INDEX2(0,0,8)]+=-B_0_0*w[47][l] - B_0_1*w[38][l] - B_0_6*w[30][l] - B_0_7*w[46][l] + B_1_0*w[44][l] - B_1_2*w[39][l] - B_1_5*w[31][l] + B_1_7*w[45][l] - B_2_0*w[40][l] - B_2_3*w[32][l] - B_2_4*w[43][l] - B_2_7*w[41][l] + tmp132 + tmp133 + tmp208 + tmp209 + tmp210 + tmp211;
                    EM_S[INDEX2(1,0,8)]+= B_0_0*w[47][l] + B_0_1*w[38][l] + B_0_6*w[30][l] + B_0_7*w[46][l] + tmp148 + tmp149 + tmp242 + tmp243 + tmp244 + tmp245 + tmp246 + tmp247 + tmp248 + tmp249;
                    EM_S[INDEX2(2,0,8)]+=-B_1_0*w[44][l] + B_1_2*w[39][l] + B_1_5*w[31][l] - B_1_7*w[45][l] + tmp138 + tmp139 + tmp140 + tmp141 + tmp142 + tmp143 + tmp144 + tmp145 + tmp146 + tmp147;
                    EM_S[INDEX2(3,0,8)]+=tmp40 + tmp41 + tmp42 + tmp43 + tmp44 + tmp45 + tmp46 + tmp47 + tmp48 + tmp49;
                    EM_S[INDEX2(4,0,8)]+=B_2_0*w[40][l] + B_2_3*w[32][l] + B_2_4*w[43][l] + B_2_7*w[41][l] + tmp114 + tmp115 + tmp116 + tmp117 + tmp118 + tmp119 + tmp120 + tmp121 + tmp2 + tmp3;
                    EM_S[INDEX2(5,0,8)]+=tmp191 + tmp192 + tmp193 + tmp194 + tmp198 + tmp199 + tmp214 + tmp215 + tmp216 + tmp217;
                    EM_S[INDEX2(6,0,8)]+=tmp228 + tmp229 + tmp230 + tmp231 + tmp232 + tmp233 + tmp234 + tmp235 + tmp236 + tmp237;
                    EM_S[INDEX2(7,0,8)]+=tmp258 + tmp259 + tmp80 + tmp81 + tmp94 + tmp95;
                    EM_S[INDEX2(0,1,8)]+=-B_0_0*w[38][l] - B_0_1*w[47][l] - B_0_6*w[46][l] - B_0_7*w[30][l] + tmp128 + tmp129 + tmp242 + tmp243 + tmp244 + tmp245 + tmp246 + tmp247 + tmp248 + tmp249;
                    EM_S[INDEX2(1,1,8)]+= B_0_0*w[38][l] + B_0_1*w[47][l] + B_0_6*w[46][l] + B_0_7*w[30][l] + B_1_1*w[44][l] - B_1_3*w[39][l] - B_1_4*w[31][l] + B_1_6*w[45][l] - B_2_1*w[40][l] - B_2_2*w[32][l] - B_2_5*w[43][l] - B_2_6*w[41][l] + tmp152 + tmp155 + tmp164 + tmp165 + tmp168 + tmp169;
                    EM_S[INDEX2(2,1,8)]+=tmp100 + tmp101 + tmp102 + tmp41 + tmp42 + tmp43 + tmp44 + tmp45 + tmp49 + tmp96;
                    EM_S[INDEX2(3,1,8)]+=-B_1_1*w[44][l] + B_1_3*w[39][l] + B_1_4*w[31][l] - B_1_6*w[45][l] + tmp20 + tmp21 + tmp22 + tmp24 + tmp25 + tmp26 + tmp27 + tmp28 + tmp33 + tmp39;
                    EM_S[INDEX2(4,1,8)]+=tmp190 + tmp191 + tmp192 + tmp193 + tmp194 + tmp195 + tmp196 + tmp197 + tmp198 + tmp199;
                    EM_S[INDEX2(5,1,8)]+= B_2_1*w[40][l] + B_2_2*w[32][l] + B_2_5*w[43][l] + B_2_6*w[41][l] + tmp82 + tmp83 + tmp84 + tmp85 + tmp86 + tmp87 + tmp88 + tmp89 + tmp90 + tmp91;
                    EM_S[INDEX2(6,1,8)]+=tmp258 + tmp259 + tmp63 + tmp64 + tmp94 + tmp95;
                    EM_S[INDEX2(7,1,8)]+=tmp181 + tmp182 + tmp184 + tmp186 + tmp187 + tmp188 + tmp218 + tmp219 + tmp220 + tmp221;
                    EM_S[INDEX2(0,2,8)]+=-B_1_0*w[39][l] + B_1_2*w[44][l] + B_1_5*w[45][l] - B_1_7*w[31][l] + tmp138 + tmp139 + tmp140 + tmp142 + tmp143 + tmp144 + tmp145 + tmp146 + tmp173 + tmp179;
                    EM_S[INDEX2(1,2,8)]+=tmp103 + tmp40 + tmp42 + tmp43 + tmp46 + tmp47 + tmp48 + tmp97 + tmp98 + tmp99;
                    EM_S[INDEX2(2,2,8)]+=-B_0_2*w[47][l] - B_0_3*w[38][l] - B_0_4*w[30][l] - B_0_5*w[46][l] + B_1_0*w[39][l] - B_1_2*w[44][l] - B_1_5*w[45][l] + B_1_7*w[31][l] - B_2_1*w[32][l] - B_2_2*w[40][l] - B_2_5*w[41][l] - B_2_6*w[43][l] + tmp153 + tmp154 + tmp164 + tmp165 + tmp166 + tmp167;
                    EM_S[INDEX2(3,2,8)]+= B_0_2*w[47][l] + B_0_3*w[38][l] + B_0_4*w[30][l] + B_0_5*w[46][l] + tmp200 + tmp201 + tmp202 + tmp203 + tmp204 + tmp205 + tmp206 + tmp207 + tmp56 + tmp57;
                    EM_S[INDEX2(4,2,8)]+=tmp229 + tmp230 + tmp232 + tmp234 + tmp235 + tmp236 + tmp238 + tmp239 + tmp240 + tmp241;
                    EM_S[INDEX2(5,2,8)]+=tmp258 + tmp259 + tmp62 + tmp65 + tmp80 + tmp81;
                    EM_S[INDEX2(6,2,8)]+= B_2_1*w[32][l] + B_2_2*w[40][l] + B_2_5*w[41][l] + B_2_6*w[43][l] + tmp70 + tmp71 + tmp74 + tmp75 + tmp76 + tmp77 + tmp78 + tmp79 + tmp84 + tmp85;
                    EM_S[INDEX2(7,2,8)]+=tmp104 + tmp105 + tmp106 + tmp107 + tmp108 + tmp109 + tmp110 + tmp111 + tmp112 + tmp113;
                    EM_S[INDEX2(0,3,8)]+=tmp100 + tmp101 + tmp102 + tmp103 + tmp42 + tmp43 + tmp96 + tmp97 + tmp98 + tmp99;
                    EM_S[INDEX2(1,3,8)]+=-B_1_1*w[39][l] + B_1_3*w[44][l] + B_1_4*w[45][l] - B_1_6*w[31][l] + tmp20 + tmp21 + tmp22 + tmp23 + tmp24 + tmp25 + tmp26 + tmp27 + tmp28 + tmp29;
                    EM_S[INDEX2(2,3,8)]+=-B_0_2*w[38][l] - B_0_3*w[47][l] - B_0_4*w[46][l] - B_0_5*w[30][l] + tmp200 + tmp201 + tmp202 + tmp203 + tmp204 + tmp205 + tmp206 + tmp207 + tmp222 + tmp223;
                    EM_S[INDEX2(3,3,8)]+= B_0_2*w[38][l] + B_0_3*w[47][l] + B_0_4*w[46][l] + B_0_5*w[30][l] + B_1_1*w[39][l] - B_1_3*w[44][l] - B_1_4*w[45][l] + B_1_6*w[31][l] - B_2_0*w[32][l] - B_2_3*w[40][l] - B_2_4*w[41][l] - B_2_7*w[43][l] + tmp132 + tmp133 + tmp134 + tmp135 + tmp136 + tmp137;
                    EM_S[INDEX2(4,3,8)]+=tmp258 + tmp259 + tmp62 + tmp63 + tmp64 + tmp65;
                    EM_S[INDEX2(5,3,8)]+=tmp180 + tmp181 + tmp182 + tmp183 + tmp184 + tmp185 + tmp186 + tmp187 + tmp188 + tmp189;
                    EM_S[INDEX2(6,3,8)]+=tmp105 + tmp106 + tmp107 + tmp108 + tmp112 + tmp113 + tmp156 + tmp157 + tmp158 + tmp159;
                    EM_S[INDEX2(7,3,8)]+= B_2_0*w[32][l] + B_2_3*w[40][l] + B_2_4*w[41][l] + B_2_7*w[43][l] + tmp0 + tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6 + tmp7 + tmp8 + tmp9;
                    EM_S[INDEX2(0,4,8)]+=-B_2_0*w[43][l] - B_2_3*w[41][l] - B_2_4*w[40][l] - B_2_7*w[32][l] + tmp114 + tmp115 + tmp116 + tmp117 + tmp118 + tmp119 + tmp120 + tmp121 + tmp92 + tmp93;
                    EM_S[INDEX2(1,4,8)]+=tmp193 + tmp199 + tmp214 + tmp215 + tmp216 + tmp217 + tmp224 + tmp225 + tmp226 + tmp227;
                    EM_S[INDEX2(2,4,8)]+=tmp228 + tmp231 + tmp233 + tmp234 + tmp235 + tmp237 + tmp260 + tmp261 + tmp262 + tmp263;
                    EM_S[INDEX2(3,4,8)]+=tmp60 + tmp61 + tmp80 + tmp81 + tmp94 + tmp95;
                    EM_S[INDEX2(4,4,8)]+=-B_0_2*w[30][l] - B_0_3*w[46][l] - B_0_4*w[47][l] - B_0_5*w[38][l] - B_1_1*w[31][l] + B_1_3*w[45][l] + B_1_4*w[44][l] - B_1_6*w[39][l] + B_2_0*w[43][l] + B_2_3*w[41][l] + B_2_4*w[40][l] + B_2_7*w[32][l] + tmp150 + tmp151 + tmp152 + tmp153 + tmp154 + tmp155;
                    EM_S[INDEX2(5,4,8)]+= B_0_2*w[30][l] + B_0_3*w[46][l] + B_0_4*w[47][l] + B_0_5*w[38][l] + tmp50 + tmp51 + tmp52 + tmp53 + tmp54 + tmp55 + tmp56 + tmp57 + tmp58 + tmp59;
                    EM_S[INDEX2(6,4,8)]+= B_1_1*w[31][l] - B_1_3*w[45][l] - B_1_4*w[44][l] + B_1_6*w[39][l] + tmp30 + tmp31 + tmp32 + tmp33 + tmp34 + tmp35 + tmp36 + tmp37 + tmp38 + tmp39;
                    EM_S[INDEX2(7,4,8)]+=tmp12 + tmp13 + tmp250 + tmp251 + tmp252 + tmp253 + tmp66 + tmp67 + tmp68 + tmp69;
                    EM_S[INDEX2(0,5,8)]+=tmp190 + tmp193 + tmp195 + tmp196 + tmp197 + tmp199 + tmp224 + tmp225 + tmp226 + tmp227;
                    EM_S[INDEX2(1,5,8)]+=-B_2_1*w[43][l] - B_2_2*w[41][l] - B_2_5*w[40][l] - B_2_6*w[32][l] + tmp72 + tmp73 + tmp82 + tmp83 + tmp86 + tmp87 + tmp88 + tmp89 + tmp90 + tmp91;
                    EM_S[INDEX2(2,5,8)]+=tmp60 + tmp61 + tmp63 + tmp64 + tmp94 + tmp95;
                    EM_S[INDEX2(3,5,8)]+=tmp186 + tmp187 + tmp218 + tmp219 + tmp220 + tmp221 + tmp254 + tmp255 + tmp256 + tmp257;
                    EM_S[INDEX2(4,5,8)]+=-B_0_2*w[46][l] - B_0_3*w[30][l] - B_0_4*w[38][l] - B_0_5*w[47][l] + tmp222 + tmp223 + tmp50 + tmp51 + tmp52 + tmp53 + tmp54 + tmp55 + tmp58 + tmp59;
                    EM_S[INDEX2(5,5,8)]+= B_0_2*w[46][l] + B_0_3*w[30][l] + B_0_4*w[38][l] + B_0_5*w[47][l] - B_1_0*w[31][l] + B_1_2*w[45][l] + B_1_5*w[44][l] - B_1_7*w[39][l] + B_2_1*w[43][l] + B_2_2*w[41][l] + B_2_5*w[40][l] + B_2_6*w[32][l] + tmp135 + tmp136 + tmp208 + tmp211 + tmp212 + tmp213;
                    EM_S[INDEX2(6,5,8)]+=tmp10 + tmp12 + tmp13 + tmp16 + tmp17 + tmp18 + tmp250 + tmp251 + tmp252 + tmp253;
                    EM_S[INDEX2(7,5,8)]+= B_1_0*w[31][l] - B_1_2*w[45][l] - B_1_5*w[44][l] + B_1_7*w[39][l] + tmp141 + tmp147 + tmp170 + tmp171 + tmp172 + tmp174 + tmp175 + tmp176 + tmp177 + tmp178;
                    EM_S[INDEX2(0,6,8)]+=tmp234 + tmp235 + tmp238 + tmp239 + tmp240 + tmp241 + tmp260 + tmp261 + tmp262 + tmp263;
                    EM_S[INDEX2(1,6,8)]+=tmp60 + tmp61 + tmp62 + tmp65 + tmp80 + tmp81;
                    EM_S[INDEX2(2,6,8)]+=-B_2_1*w[41][l] - B_2_2*w[43][l] - B_2_5*w[32][l] - B_2_6*w[40][l] + tmp70 + tmp71 + tmp72 + tmp73 + tmp74 + tmp75 + tmp76 + tmp77 + tmp78 + tmp79;
                    EM_S[INDEX2(3,6,8)]+=tmp104 + tmp107 + tmp109 + tmp110 + tmp111 + tmp113 + tmp160 + tmp161 + tmp162 + tmp163;
                    EM_S[INDEX2(4,6,8)]+= B_1_1*w[45][l] - B_1_3*w[31][l] - B_1_4*w[39][l] + B_1_6*w[44][l] + tmp23 + tmp29 + tmp30 + tmp31 + tmp32 + tmp34 + tmp35 + tmp36 + tmp37 + tmp38;
                    EM_S[INDEX2(5,6,8)]+=tmp11 + tmp12 + tmp13 + tmp14 + tmp15 + tmp19 + tmp66 + tmp67 + tmp68 + tmp69;
                    EM_S[INDEX2(6,6,8)]+=-B_0_0*w[30][l] - B_0_1*w[46][l] - B_0_6*w[47][l] - B_0_7*w[38][l] - B_1_1*w[45][l] + B_1_3*w[31][l] + B_1_4*w[39][l] - B_1_6*w[44][l] + B_2_1*w[41][l] + B_2_2*w[43][l] + B_2_5*w[32][l] + B_2_6*w[40][l] + tmp134 + tmp137 + tmp209 + tmp210 + tmp212 + tmp213;
                    EM_S[INDEX2(7,6,8)]+= B_0_0*w[30][l] + B_0_1*w[46][l] + B_0_6*w[47][l] + B_0_7*w[38][l] + tmp122 + tmp123 + tmp124 + tmp125 + tmp126 + tmp127 + tmp130 + tmp131 + tmp148 + tmp149;
                    EM_S[INDEX2(0,7,8)]+=tmp60 + tmp61 + tmp62 + tmp63 + tmp64 + tmp65;
                    EM_S[INDEX2(1,7,8)]+=tmp180 + tmp183 + tmp185 + tmp186 + tmp187 + tmp189 + tmp254 + tmp255 + tmp256 + tmp257;
                    EM_S[INDEX2(2,7,8)]+=tmp107 + tmp113 + tmp156 + tmp157 + tmp158 + tmp159 + tmp160 + tmp161 + tmp162 + tmp163;
                    EM_S[INDEX2(3,7,8)]+=-B_2_0*w[41][l] - B_2_3*w[43][l] - B_2_4*w[32][l] - B_2_7*w[40][l] + tmp0 + tmp1 + tmp4 + tmp5 + tmp6 + tmp7 + tmp8 + tmp9 + tmp92 + tmp93;
                    EM_S[INDEX2(4,7,8)]+=tmp10 + tmp11 + tmp12 + tmp13 + tmp14 + tmp15 + tmp16 + tmp17 + tmp18 + tmp19;
                    EM_S[INDEX2(5,7,8)]+= B_1_0*w[45][l] - B_1_2*w[31][l] - B_1_5*w[39][l] + B_1_7*w[44][l] + tmp170 + tmp171 + tmp172 + tmp173 + tmp174 + tmp175 + tmp176 + tmp177 + tmp178 + tmp179;
                    EM_S[INDEX2(6,7,8)]+=-B_0_0*w[46][l] - B_0_1*w[30][l] - B_0_6*w[38][l] - B_0_7*w[47][l] + tmp122 + tmp123 + tmp124 + tmp125 + tmp126 + tmp127 + tmp128 + tmp129 + tmp130 + tmp131;
                    EM_S[INDEX2(7,7,8)]+= B_0_0*w[46][l] + B_0_1*w[30][l] + B_0_6*w[38][l] + B_0_7*w[47][l] - B_1_0*w[45][l]+ B_1_2*w[31][l] + B_1_5*w[39][l] - B_1_7*w[44][l] + B_2_0*w[41][l] + B_2_3*w[43][l] + B_2_4*w[32][l] + B_2_7*w[40][l] + tmp150 + tmp151 + tmp166 + tmp167 + tmp168 + tmp169;
                } else { // constant data
                    const Scalar wB0 = B_p[0]*w[53][l];
                    const Scalar wB1 = B_p[1]*w[51][l];
                    const Scalar wB2 = B_p[2]*w[50][l];
                    EM_S[INDEX2(0,0,8)]+= 4.*wB0 + 4.*wB1 - 4.*wB2;
                    EM_S[INDEX2(1,0,8)]+=-4.*wB0 + 2.*wB1 - 2.*wB2;
                    EM_S[INDEX2(2,0,8)]+= 2.*wB0 - 4.*wB1 - 2.*wB2;
                    EM_S[INDEX2(3,0,8)]+=-2.*wB0 - 2.*wB1 - wB2;
                    EM_S[INDEX2(4,0,8)]+= 2.*wB0 + 2.*wB1 + 4.*wB2;
                    EM_S[INDEX2(5,0,8)]+=-2.*wB0 +    wB1 + 2.*wB2;
                    EM_S[INDEX2(6,0,8)]+=    wB0 - 2.*wB1 + 2.*wB2;
                    EM_S[INDEX2(7,0,8)]+=   -wB0 -    wB1 + wB2;
                    EM_S[INDEX2(0,1,8)]+= 4.*wB0 + 2.*wB1 - 2.*wB2;
                    EM_S[INDEX2(1,1,8)]+=-4.*wB0 + 4.*wB1 - 4.*wB2;
                    EM_S[INDEX2(2,1,8)]+= 2.*wB0 - 2.*wB1 - wB2;
                    EM_S[INDEX2(3,1,8)]+=-2.*wB0 - 4.*wB1 - 2.*wB2;
                    EM_S[INDEX2(4,1,8)]+= 2.*wB0 +    wB1 + 2.*wB2;
                    EM_S[INDEX2(5,1,8)]+=-2.*wB0 + 2.*wB1 + 4.*wB2;
                    EM_S[INDEX2(6,1,8)]+=    wB0 -    wB1 + wB2;
                    EM_S[INDEX2(7,1,8)]+=   -wB0 - 2.*wB1 + 2.*wB2;
                    EM_S[INDEX2(0,2,8)]+= 2.*wB0 + 4.*wB1 - 2.*wB2;
                    EM_S[INDEX2(1,2,8)]+=-2.*wB0 + 2.*wB1 - wB2;
                    EM_S[INDEX2(2,2,8)]+= 4.*wB0 - 4.*wB1 - 4.*wB2;
                    EM_S[INDEX2(3,2,8)]+=-4.*wB0 - 2.*wB1 - 2.*wB2;
                    EM_S[INDEX2(4,2,8)]+=    wB0 + 2.*wB1 + 2.*wB2;
                    EM_S[INDEX2(5,2,8)]+=   -wB0 +    wB1 + wB2;
                    EM_S[INDEX2(6,2,8)]+= 2.*wB0 - 2.*wB1 + 4.*wB2;
                    EM_S[INDEX2(7,2,8)]+=-2.*wB0 -    wB1 + 2.*wB2;
                    EM_S[INDEX2(0,3,8)]+= 2.*wB0 + 2.*wB1 - wB2;
                    EM_S[INDEX2(1,3,8)]+=-2.*wB0 + 4.*wB1 - 2.*wB2;
                    EM_S[INDEX2(2,3,8)]+= 4.*wB0 - 2.*wB1 - 2.*wB2;
                    EM_S[INDEX2(3,3,8)]+=-4.*wB0 - 4.*wB1 - 4.*wB2;
                    EM_S[INDEX2(4,3,8)]+=    wB0 +    wB1 + wB2;
                    EM_S[INDEX2(5,3,8)]+=   -wB0 + 2.*wB1 + 2.*wB2;
                    EM_S[INDEX2(6,3,8)]+= 2.*wB0 -    wB1 + 2.*wB2;
                    EM_S[INDEX2(7,3,8)]+=-2.*wB0 - 2.*wB1 + 4.*wB2;
                    EM_S[INDEX2(0,4,8)]+= 2.*wB0 + 2.*wB1 - 4.*wB2;
                    EM_S[INDEX2(1,4,8)]+=-2.*wB0 +    wB1 - 2.*wB2;
                    EM_S[INDEX2(2,4,8)]+=    wB0 - 2.*wB1 - 2.*wB2;
                    EM_S[INDEX2(3,4,8)]+=   -wB0 -    wB1 - wB2;
                    EM_S[INDEX2(4,4,8)]+= 4.*wB0 + 4.*wB1 + 4.*wB2;
                    EM_S[INDEX2(5,4,8)]+=-4.*wB0 + 2.*wB1 + 2.*wB2;
                    EM_S[INDEX2(6,4,8)]+= 2.*wB0 - 4.*wB1 + 2.*wB2;
                    EM_S[INDEX2(7,4,8)]+=-2.*wB0 - 2.*wB1 + wB2;
                    EM_S[INDEX2(0,5,8)]+= 2.*wB0 +    wB1 - 2.*wB2;
                    EM_S[INDEX2(1,5,8)]+=-2.*wB0 + 2.*wB1 - 4.*wB2;
                    EM_S[INDEX2(2,5,8)]+=    wB0 -    wB1 - wB2;
                    EM_S[INDEX2(3,5,8)]+=   -wB0 - 2.*wB1 - 2.*wB2;
                    EM_S[INDEX2(4,5,8)]+= 4.*wB0 + 2.*wB1 + 2.*wB2;
                    EM_S[INDEX2(5,5,8)]+=-4.*wB0 + 4.*wB1 + 4.*wB2;
                    EM_S[INDEX2(6,5,8)]+= 2.*wB0 - 2.*wB1 + wB2;
                    EM_S[INDEX2(7,5,8)]+=-2.*wB0 - 4.*wB1 + 2.*wB2;
                    EM_S[INDEX2(0,6,8)]+=    wB0 + 2.*wB1 - 2.*wB2;
                    EM_S[INDEX2(1,6,8)]+=   -wB0 +    wB1 - wB2;
                    EM_S[INDEX2(2,6,8)]+= 2.*wB0 - 2.*wB1 - 4.*wB2;
                    EM_S[INDEX2(3,6,8)]+=-2.*wB0 -    wB1 - 2.*wB2;
                    EM_S[INDEX2(4,6,8)]+= 2.*wB0 + 4.*wB1 + 2.*wB2;
                    EM_S[INDEX2(5,6,8)]+=-2.*wB0 + 2.*wB1 + wB2;
                    EM_S[INDEX2(6,6,8)]+= 4.*wB0 - 4.*wB1 + 4.*wB2;
                    EM_S[INDEX2(7,6,8)]+=-4.*wB0 - 2.*wB1 + 2.*wB2;
                    EM_S[INDEX2(0,7,8)]+=    wB0 +    wB1 - wB2;
                    EM_S[INDEX2(1,7,8)]+=   -wB0 + 2.*wB1 - 2.*wB2;
                    EM_S[INDEX2(2,7,8)]+= 2.*wB0 -    wB1 - 2.*wB2;
                    EM_S[INDEX2(3,7,8)]+=-2.*wB0 - 2.*wB1 - 4.*wB2;
                    EM_S[INDEX2(4,7,8)]+= 2.*wB0 + 2.*wB1 + wB2;
                    EM_S[INDEX2(5,7,8)]+=-2.*wB0 + 4.*wB1 + 2.*wB2;
                    EM_S[INDEX2(6,7,8)]+= 4.*wB0 - 2.*wB1 + 2.*wB2;
                    EM_S[INDEX2(7,7,8)]+=-4.*wB0 - 4.*wB1 + 4.*wB2;
                }
            }
            ///////////////
            // process C //
            ///////////////
            if (!C.isEmpty()) {
                const Scalar* C_p = C.getSampleDataRO(id, zero);
                if (C.actsExpanded()) {
                    const Scalar C_0_0 = C_p[INDEX2(0,0,3)];
                    const Scalar C_1_0 = C_p[INDEX2(1,0,3)];
                    const Scalar C_2_0 = C_p[INDEX2(2,0,3)];
                    const Scalar C_0_1 = C_p[INDEX2(0,1,3)];
                    const Scalar C_1_1 = C_p[INDEX2(1,1,3)];
                    const Scalar C_2_1 = C_p[INDEX2(2,1,3)];
                    const Scalar C_0_2 = C_p[INDEX2(0,2,3)];
                    const Scalar C_1_2 = C_p[INDEX2(1,2,3)];
                    const Scalar C_2_2 = C_p[INDEX2(2,2,3)];
                    const Scalar C_0_3 = C_p[INDEX2(0,3,3)];
                    const Scalar C_1_3 = C_p[INDEX2(1,3,3)];
                    const Scalar C_2_3 = C_p[INDEX2(2,3,3)];
                    const Scalar C_0_4 = C_p[INDEX2(0,4,3)];
                    const Scalar C_1_4 = C_p[INDEX2(1,4,3)];
                    const Scalar C_2_4 = C_p[INDEX2(2,4,3)];
                    const Scalar C_0_5 = C_p[INDEX2(0,5,3)];
                    const Scalar C_1_5 = C_p[INDEX2(1,5,3)];
                    const Scalar C_2_5 = C_p[INDEX2(2,5,3)];
                    const Scalar C_0_6 = C_p[INDEX2(0,6,3)];
                    const Scalar C_1_6 = C_p[INDEX2(1,6,3)];
                    const Scalar C_2_6 = C_p[INDEX2(2,6,3)];
                    const Scalar C_0_7 = C_p[INDEX2(0,7,3)];
                    const Scalar C_1_7 = C_p[INDEX2(1,7,3)];
                    const Scalar C_2_7 = C_p[INDEX2(2,7,3)];
                    const Scalar tmp0 =  w[38][l]*(C_0_3 + C_0_7);
                    const Scalar tmp1 =  w[31][l]*(C_1_0 + C_1_4);
                    const Scalar tmp2 =  w[42][l]*(-C_2_1 - C_2_2);
                    const Scalar tmp3 =  w[35][l]*(-C_2_5 - C_2_6);
                    const Scalar tmp4 =  w[37][l]*(C_1_2 + C_1_6);
                    const Scalar tmp5 =  w[39][l]*(C_1_3 + C_1_7);
                    const Scalar tmp6 =  w[36][l]*(C_0_2 + C_0_6);
                    const Scalar tmp7 =  w[33][l]*(C_0_1 + C_0_5);
                    const Scalar tmp8 =  w[30][l]*(C_0_0 + C_0_4);
                    const Scalar tmp9 =  w[34][l]*(C_1_1 + C_1_5);
                    const Scalar tmp10 = w[38][l]*(C_0_4 + C_0_6);
                    const Scalar tmp11 = w[31][l]*(C_1_2 + C_1_3);
                    const Scalar tmp12 = w[42][l]*(C_2_0 + C_2_1 + C_2_2 + C_2_3);
                    const Scalar tmp13 = w[35][l]*(C_2_4 + C_2_5 + C_2_6 + C_2_7);
                    const Scalar tmp14 = w[37][l]*(C_1_0 + C_1_1);
                    const Scalar tmp15 = w[39][l]*(C_1_4 + C_1_5);
                    const Scalar tmp16 = w[36][l]*(C_0_5 + C_0_7);
                    const Scalar tmp17 = w[33][l]*(C_0_0 + C_0_2);
                    const Scalar tmp18 = w[30][l]*(C_0_1 + C_0_3);
                    const Scalar tmp19 = w[34][l]*(C_1_6 + C_1_7);
                    const Scalar tmp20 = w[38][l]*(C_0_1 + C_0_3);
                    const Scalar tmp21 = w[42][l]*(-C_2_0 - C_2_2);
                    const Scalar tmp22 = w[35][l]*(-C_2_5 - C_2_7);
                    const Scalar tmp23 = w[37][l]*(C_1_2 + C_1_7);
                    const Scalar tmp24 = w[32][l]*(-C_2_4 - C_2_6);
                    const Scalar tmp25 = w[36][l]*(C_0_0 + C_0_2);
                    const Scalar tmp26 = w[33][l]*(C_0_5 + C_0_7);
                    const Scalar tmp27 = w[30][l]*(C_0_4 + C_0_6);
                    const Scalar tmp28 = w[43][l]*(-C_2_1 - C_2_3);
                    const Scalar tmp29 = w[34][l]*(C_1_0 + C_1_5);
                    const Scalar tmp30 = w[38][l]*(-C_0_4 - C_0_6);
                    const Scalar tmp31 = w[42][l]*(C_2_5 + C_2_7);
                    const Scalar tmp32 = w[35][l]*(C_2_0 + C_2_2);
                    const Scalar tmp33 = w[37][l]*(-C_1_0 - C_1_5);
                    const Scalar tmp34 = w[32][l]*(C_2_1 + C_2_3);
                    const Scalar tmp35 = w[36][l]*(-C_0_5 - C_0_7);
                    const Scalar tmp36 = w[33][l]*(-C_0_0 - C_0_2);
                    const Scalar tmp37 = w[30][l]*(-C_0_1 - C_0_3);
                    const Scalar tmp38 = w[43][l]*(C_2_4 + C_2_6);
                    const Scalar tmp39 = w[34][l]*(-C_1_2 - C_1_7);
                    const Scalar tmp40 = w[38][l]*(-C_0_1 - C_0_3);
                    const Scalar tmp41 = w[31][l]*(-C_1_4 - C_1_5);
                    const Scalar tmp42 = w[42][l]*(-C_2_4 - C_2_5 - C_2_6 - C_2_7);
                    const Scalar tmp43 = w[35][l]*(-C_2_0 - C_2_1 - C_2_2 - C_2_3);
                    const Scalar tmp44 = w[37][l]*(-C_1_6 - C_1_7);
                    const Scalar tmp45 = w[39][l]*(-C_1_2 - C_1_3);
                    const Scalar tmp46 = w[36][l]*(-C_0_0 - C_0_2);
                    const Scalar tmp47 = w[33][l]*(-C_0_5 - C_0_7);
                    const Scalar tmp48 = w[30][l]*(-C_0_4 - C_0_6);
                    const Scalar tmp49 = w[34][l]*(-C_1_0 - C_1_1);
                    const Scalar tmp50 = w[31][l]*(-C_1_2 - C_1_3);
                    const Scalar tmp51 = w[42][l]*(C_2_6 + C_2_7);
                    const Scalar tmp52 = w[35][l]*(C_2_0 + C_2_1);
                    const Scalar tmp53 = w[37][l]*(-C_1_0 - C_1_1);
                    const Scalar tmp54 = w[32][l]*(C_2_2 + C_2_3);
                    const Scalar tmp55 = w[39][l]*(-C_1_4 - C_1_5);
                    const Scalar tmp56 = w[36][l]*(-C_0_1 - C_0_7);
                    const Scalar tmp57 = w[33][l]*(-C_0_0 - C_0_6);
                    const Scalar tmp58 = w[43][l]*(C_2_4 + C_2_5);
                    const Scalar tmp59 = w[34][l]*(-C_1_6 - C_1_7);
                    const Scalar tmp60 = w[42][l]*(C_2_4 + C_2_5 + C_2_6 + C_2_7);
                    const Scalar tmp61 = w[35][l]*(C_2_0 + C_2_1 + C_2_2 + C_2_3);
                    const Scalar tmp62 = w[37][l]*(C_1_2 + C_1_3 + C_1_6 + C_1_7);
                    const Scalar tmp63 = w[36][l]*(C_0_0 + C_0_2 + C_0_4 + C_0_6);
                    const Scalar tmp64 = w[33][l]*(C_0_1 + C_0_3 + C_0_5 + C_0_7);
                    const Scalar tmp65 = w[34][l]*(C_1_0 + C_1_1 + C_1_4 + C_1_5);
                    const Scalar tmp66 = w[38][l]*(-C_0_5 - C_0_7);
                    const Scalar tmp67 = w[36][l]*(-C_0_4 - C_0_6);
                    const Scalar tmp68 = w[33][l]*(-C_0_1 - C_0_3);
                    const Scalar tmp69 = w[30][l]*(-C_0_0 - C_0_2);
                    const Scalar tmp70 = w[38][l]*(-C_0_2 - C_0_6);
                    const Scalar tmp71 = w[31][l]*(C_1_1 + C_1_5);
                    const Scalar tmp72 = w[42][l]*(C_2_4 + C_2_7);
                    const Scalar tmp73 = w[35][l]*(C_2_0 + C_2_3);
                    const Scalar tmp74 = w[37][l]*(C_1_3 + C_1_7);
                    const Scalar tmp75 = w[39][l]*(C_1_2 + C_1_6);
                    const Scalar tmp76 = w[36][l]*(-C_0_3 - C_0_7);
                    const Scalar tmp77 = w[33][l]*(-C_0_0 - C_0_4);
                    const Scalar tmp78 = w[30][l]*(-C_0_1 - C_0_5);
                    const Scalar tmp79 = w[34][l]*(C_1_0 + C_1_4);
                    const Scalar tmp80 = w[36][l]*(-C_0_1 - C_0_3 - C_0_5 - C_0_7);
                    const Scalar tmp81 = w[33][l]*(-C_0_0 - C_0_2 - C_0_4 - C_0_6);
                    const Scalar tmp82 = w[38][l]*(C_0_1 + C_0_5);
                    const Scalar tmp83 = w[31][l]*(-C_1_2 - C_1_6);
                    const Scalar tmp84 = w[42][l]*(-C_2_0 - C_2_3);
                    const Scalar tmp85 = w[35][l]*(-C_2_4 - C_2_7);
                    const Scalar tmp86 = w[37][l]*(-C_1_0 - C_1_4);
                    const Scalar tmp87 = w[39][l]*(-C_1_1 - C_1_5);
                    const Scalar tmp88 = w[36][l]*(C_0_0 + C_0_4);
                    const Scalar tmp89 = w[33][l]*(C_0_3 + C_0_7);
                    const Scalar tmp90 = w[30][l]*(C_0_2 + C_0_6);
                    const Scalar tmp91 = w[34][l]*(-C_1_3 - C_1_7);
                    const Scalar tmp92 = w[42][l]*(C_2_5 + C_2_6);
                    const Scalar tmp93 = w[35][l]*(C_2_1 + C_2_2);
                    const Scalar tmp94 = w[37][l]*(-C_1_0 - C_1_1 - C_1_4 - C_1_5);
                    const Scalar tmp95 = w[34][l]*(-C_1_2 - C_1_3 - C_1_6 - C_1_7);
                    const Scalar tmp96 = w[38][l]*(C_0_0 + C_0_2);
                    const Scalar tmp97 = w[31][l]*(C_1_6 + C_1_7);
                    const Scalar tmp98 = w[37][l]*(C_1_4 + C_1_5);
                    const Scalar tmp99 = w[39][l]*(C_1_0 + C_1_1);
                    const Scalar tmp100 = w[36][l]*(C_0_1 + C_0_3);
                    const Scalar tmp101 = w[33][l]*(C_0_4 + C_0_6);
                    const Scalar tmp102 = w[30][l]*(C_0_5 + C_0_7);
                    const Scalar tmp103 = w[34][l]*(C_1_2 + C_1_3);
                    const Scalar tmp104 = w[38][l]*(-C_0_3 - C_0_7);
                    const Scalar tmp105 = w[42][l]*(-C_2_4 - C_2_5);
                    const Scalar tmp106 = w[35][l]*(-C_2_2 - C_2_3);
                    const Scalar tmp107 = w[37][l]*(C_1_0 + C_1_1 + C_1_4 + C_1_5);
                    const Scalar tmp108 = w[32][l]*(-C_2_0 - C_2_1);
                    const Scalar tmp109 = w[36][l]*(-C_0_2 - C_0_6);
                    const Scalar tmp110 = w[33][l]*(-C_0_1 - C_0_5);
                    const Scalar tmp111 = w[30][l]*(-C_0_0 - C_0_4);
                    const Scalar tmp112 = w[43][l]*(-C_2_6 - C_2_7);
                    const Scalar tmp113 = w[34][l]*(C_1_2 + C_1_3 + C_1_6 + C_1_7);
                    const Scalar tmp114 = w[38][l]*(-C_0_0 - C_0_4);
                    const Scalar tmp115 = w[31][l]*(-C_1_3 - C_1_7);
                    const Scalar tmp116 = w[37][l]*(-C_1_1 - C_1_5);
                    const Scalar tmp117 = w[39][l]*(-C_1_0 - C_1_4);
                    const Scalar tmp118 = w[36][l]*(-C_0_1 - C_0_5);
                    const Scalar tmp119 = w[33][l]*(-C_0_2 - C_0_6);
                    const Scalar tmp120 = w[30][l]*(-C_0_3 - C_0_7);
                    const Scalar tmp121 = w[34][l]*(-C_1_2 - C_1_6);
                    const Scalar tmp122 = w[31][l]*(C_1_0 + C_1_1);
                    const Scalar tmp123 = w[42][l]*(C_2_4 + C_2_5);
                    const Scalar tmp124 = w[35][l]*(C_2_2 + C_2_3);
                    const Scalar tmp125 = w[37][l]*(C_1_2 + C_1_3);
                    const Scalar tmp126 = w[32][l]*(C_2_0 + C_2_1);
                    const Scalar tmp127 = w[39][l]*(C_1_6 + C_1_7);
                    const Scalar tmp128 = w[36][l]*(C_0_2 + C_0_4);
                    const Scalar tmp129 = w[33][l]*(C_0_3 + C_0_5);
                    const Scalar tmp130 = w[43][l]*(C_2_6 + C_2_7);
                    const Scalar tmp131 = w[34][l]*(C_1_4 + C_1_5);
                    const Scalar tmp132 = w[42][l]*(-C_2_5 - C_2_6);
                    const Scalar tmp133 = w[35][l]*(-C_2_1 - C_2_2);
                    const Scalar tmp134 = w[37][l]*(C_1_0 + C_1_5);
                    const Scalar tmp135 = w[36][l]*(C_0_1 + C_0_7);
                    const Scalar tmp136 = w[33][l]*(C_0_0 + C_0_6);
                    const Scalar tmp137 = w[34][l]*(C_1_2 + C_1_7);
                    const Scalar tmp138 = w[38][l]*(-C_0_0 - C_0_2);
                    const Scalar tmp139 = w[42][l]*(-C_2_1 - C_2_3);
                    const Scalar tmp140 = w[35][l]*(-C_2_4 - C_2_6);
                    const Scalar tmp141 = w[37][l]*(-C_1_1 - C_1_4);
                    const Scalar tmp142 = w[32][l]*(-C_2_5 - C_2_7);
                    const Scalar tmp143 = w[36][l]*(-C_0_1 - C_0_3);
                    const Scalar tmp144 = w[33][l]*(-C_0_4 - C_0_6);
                    const Scalar tmp145 = w[30][l]*(-C_0_5 - C_0_7);
                    const Scalar tmp146 = w[43][l]*(-C_2_0 - C_2_2);
                    const Scalar tmp147 = w[34][l]*(-C_1_3 - C_1_6);
                    const Scalar tmp148 = w[36][l]*(-C_0_3 - C_0_5);
                    const Scalar tmp149 = w[33][l]*(-C_0_2 - C_0_4);
                    const Scalar tmp150 = w[42][l]*(C_2_1 + C_2_2);
                    const Scalar tmp151 = w[35][l]*(C_2_5 + C_2_6);
                    const Scalar tmp152 = w[37][l]*(-C_1_2 - C_1_7);
                    const Scalar tmp153 = w[36][l]*(-C_0_0 - C_0_6);
                    const Scalar tmp154 = w[33][l]*(-C_0_1 - C_0_7);
                    const Scalar tmp155 = w[34][l]*(-C_1_0 - C_1_5);
                    const Scalar tmp156 = w[38][l]*(C_0_2 + C_0_6);
                    const Scalar tmp157 = w[36][l]*(C_0_3 + C_0_7);
                    const Scalar tmp158 = w[33][l]*(C_0_0 + C_0_4);
                    const Scalar tmp159 = w[30][l]*(C_0_1 + C_0_5);
                    const Scalar tmp160 = w[42][l]*(C_2_0 + C_2_1);
                    const Scalar tmp161 = w[35][l]*(C_2_6 + C_2_7);
                    const Scalar tmp162 = w[32][l]*(C_2_4 + C_2_5);
                    const Scalar tmp163 = w[43][l]*(C_2_2 + C_2_3);
                    const Scalar tmp164 = w[42][l]*(-C_2_4 - C_2_7);
                    const Scalar tmp165 = w[35][l]*(-C_2_0 - C_2_3);
                    const Scalar tmp166 = w[37][l]*(C_1_1 + C_1_4);
                    const Scalar tmp167 = w[34][l]*(C_1_3 + C_1_6);
                    const Scalar tmp168 = w[36][l]*(C_0_3 + C_0_5);
                    const Scalar tmp169 = w[33][l]*(C_0_2 + C_0_4);
                    const Scalar tmp170 = w[38][l]*(C_0_5 + C_0_7);
                    const Scalar tmp171 = w[42][l]*(C_2_4 + C_2_6);
                    const Scalar tmp172 = w[35][l]*(C_2_1 + C_2_3);
                    const Scalar tmp173 = w[37][l]*(C_1_3 + C_1_6);
                    const Scalar tmp174 = w[32][l]*(C_2_0 + C_2_2);
                    const Scalar tmp175 = w[36][l]*(C_0_4 + C_0_6);
                    const Scalar tmp176 = w[33][l]*(C_0_1 + C_0_3);
                    const Scalar tmp177 = w[30][l]*(C_0_0 + C_0_2);
                    const Scalar tmp178 = w[43][l]*(C_2_5 + C_2_7);
                    const Scalar tmp179 = w[34][l]*(C_1_1 + C_1_4);
                    const Scalar tmp180 = w[31][l]*(C_1_2 + C_1_6);
                    const Scalar tmp181 = w[42][l]*(-C_2_4 - C_2_6);
                    const Scalar tmp182 = w[35][l]*(-C_2_1 - C_2_3);
                    const Scalar tmp183 = w[37][l]*(C_1_0 + C_1_4);
                    const Scalar tmp184 = w[32][l]*(-C_2_0 - C_2_2);
                    const Scalar tmp185 = w[39][l]*(C_1_1 + C_1_5);
                    const Scalar tmp186 = w[36][l]*(C_0_1 + C_0_3 + C_0_5 + C_0_7);
                    const Scalar tmp187 = w[33][l]*(C_0_0 + C_0_2 + C_0_4 + C_0_6);
                    const Scalar tmp188 = w[43][l]*(-C_2_5 - C_2_7);
                    const Scalar tmp189 = w[34][l]*(C_1_3 + C_1_7);
                    const Scalar tmp190 = w[38][l]*(C_0_0 + C_0_4);
                    const Scalar tmp191 = w[42][l]*(-C_2_6 - C_2_7);
                    const Scalar tmp192 = w[35][l]*(-C_2_0 - C_2_1);
                    const Scalar tmp193 = w[37][l]*(-C_1_2 - C_1_3 - C_1_6 - C_1_7);
                    const Scalar tmp194 = w[32][l]*(-C_2_2 - C_2_3);
                    const Scalar tmp195 = w[36][l]*(C_0_1 + C_0_5);
                    const Scalar tmp196 = w[33][l]*(C_0_2 + C_0_6);
                    const Scalar tmp197 = w[30][l]*(C_0_3 + C_0_7);
                    const Scalar tmp198 = w[43][l]*(-C_2_4 - C_2_5);
                    const Scalar tmp199 = w[34][l]*(-C_1_0 - C_1_1 - C_1_4 - C_1_5);
                    const Scalar tmp200 = w[31][l]*(C_1_4 + C_1_5);
                    const Scalar tmp201 = w[42][l]*(-C_2_0 - C_2_1);
                    const Scalar tmp202 = w[35][l]*(-C_2_6 - C_2_7);
                    const Scalar tmp203 = w[37][l]*(C_1_6 + C_1_7);
                    const Scalar tmp204 = w[32][l]*(-C_2_4 - C_2_5);
                    const Scalar tmp205 = w[39][l]*(C_1_2 + C_1_3);
                    const Scalar tmp206 = w[43][l]*(-C_2_2 - C_2_3);
                    const Scalar tmp207 = w[34][l]*(C_1_0 + C_1_1);
                    const Scalar tmp208 = w[37][l]*(-C_1_3 - C_1_6);
                    const Scalar tmp209 = w[36][l]*(-C_0_2 - C_0_4);
                    const Scalar tmp210 = w[33][l]*(-C_0_3 - C_0_5);
                    const Scalar tmp211 = w[34][l]*(-C_1_1 - C_1_4);
                    const Scalar tmp212 = w[42][l]*(C_2_0 + C_2_3);
                    const Scalar tmp213 = w[35][l]*(C_2_4 + C_2_7);
                    const Scalar tmp214 = w[38][l]*(-C_0_1 - C_0_5);
                    const Scalar tmp215 = w[36][l]*(-C_0_0 - C_0_4);
                    const Scalar tmp216 = w[33][l]*(-C_0_3 - C_0_7);
                    const Scalar tmp217 = w[30][l]*(-C_0_2 - C_0_6);
                    const Scalar tmp218 = w[31][l]*(-C_1_0 - C_1_4);
                    const Scalar tmp219 = w[37][l]*(-C_1_2 - C_1_6);
                    const Scalar tmp220 = w[39][l]*(-C_1_3 - C_1_7);
                    const Scalar tmp221 = w[34][l]*(-C_1_1 - C_1_5);
                    const Scalar tmp222 = w[36][l]*(C_0_0 + C_0_6);
                    const Scalar tmp223 = w[33][l]*(C_0_1 + C_0_7);
                    const Scalar tmp224 = w[42][l]*(C_2_2 + C_2_3);
                    const Scalar tmp225 = w[35][l]*(C_2_4 + C_2_5);
                    const Scalar tmp226 = w[32][l]*(C_2_6 + C_2_7);
                    const Scalar tmp227 = w[43][l]*(C_2_0 + C_2_1);
                    const Scalar tmp228 = w[31][l]*(-C_1_1 - C_1_5);
                    const Scalar tmp229 = w[42][l]*(-C_2_5 - C_2_7);
                    const Scalar tmp230 = w[35][l]*(-C_2_0 - C_2_2);
                    const Scalar tmp231 = w[37][l]*(-C_1_3 - C_1_7);
                    const Scalar tmp232 = w[32][l]*(-C_2_1 - C_2_3);
                    const Scalar tmp233 = w[39][l]*(-C_1_2 - C_1_6);
                    const Scalar tmp234 = w[36][l]*(-C_0_0 - C_0_2 - C_0_4 - C_0_6);
                    const Scalar tmp235 = w[33][l]*(-C_0_1 - C_0_3 - C_0_5 - C_0_7);
                    const Scalar tmp236 = w[43][l]*(-C_2_4 - C_2_6);
                    const Scalar tmp237 = w[34][l]*(-C_1_0 - C_1_4);
                    const Scalar tmp238 = w[31][l]*(C_1_3 + C_1_7);
                    const Scalar tmp239 = w[37][l]*(C_1_1 + C_1_5);
                    const Scalar tmp240 = w[39][l]*(C_1_0 + C_1_4);
                    const Scalar tmp241 = w[34][l]*(C_1_2 + C_1_6);
                    const Scalar tmp242 = w[31][l]*(-C_1_6 - C_1_7);
                    const Scalar tmp243 = w[42][l]*(-C_2_2 - C_2_3);
                    const Scalar tmp244 = w[35][l]*(-C_2_4 - C_2_5);
                    const Scalar tmp245 = w[37][l]*(-C_1_4 - C_1_5);
                    const Scalar tmp246 = w[32][l]*(-C_2_6 - C_2_7);
                    const Scalar tmp247 = w[39][l]*(-C_1_0 - C_1_1);
                    const Scalar tmp248 = w[43][l]*(-C_2_0 - C_2_1);
                    const Scalar tmp249 = w[34][l]*(-C_1_2 - C_1_3);
                    const Scalar tmp250 = w[31][l]*(-C_1_0 - C_1_1);
                    const Scalar tmp251 = w[37][l]*(-C_1_2 - C_1_3);
                    const Scalar tmp252 = w[39][l]*(-C_1_6 - C_1_7);
                    const Scalar tmp253 = w[34][l]*(-C_1_4 - C_1_5);
                    const Scalar tmp254 = w[42][l]*(C_2_0 + C_2_2);
                    const Scalar tmp255 = w[35][l]*(C_2_5 + C_2_7);
                    const Scalar tmp256 = w[32][l]*(C_2_4 + C_2_6);
                    const Scalar tmp257 = w[43][l]*(C_2_1 + C_2_3);
                    const Scalar tmp258 = w[42][l]*(-C_2_0 - C_2_1 - C_2_2 - C_2_3);
                    const Scalar tmp259 = w[35][l]*(-C_2_4 - C_2_5 - C_2_6 - C_2_7);
                    const Scalar tmp260 = w[42][l]*(C_2_1 + C_2_3);
                    const Scalar tmp261 = w[35][l]*(C_2_4 + C_2_6);
                    const Scalar tmp262 = w[32][l]*(C_2_5 + C_2_7);
                    const Scalar tmp263 = w[43][l]*(C_2_0 + C_2_2);
                    EM_S[INDEX2(0,0,8)]+=-C_0_0*w[47][l] - C_0_1*w[38][l] - C_0_6*w[30][l] - C_0_7*w[46][l] + C_1_0*w[44][l] - C_1_2*w[39][l] - C_1_5*w[31][l] + C_1_7*w[45][l] - C_2_0*w[40][l] - C_2_3*w[32][l] - C_2_4*w[43][l] - C_2_7*w[41][l] + tmp132 + tmp133 + tmp208 + tmp209 + tmp210 + tmp211;
                    EM_S[INDEX2(1,0,8)]+=-C_0_0*w[38][l] - C_0_1*w[47][l] - C_0_6*w[46][l] - C_0_7*w[30][l] + tmp148 + tmp149 + tmp242 + tmp243 + tmp244 + tmp245 + tmp246 + tmp247 + tmp248 + tmp249;
                    EM_S[INDEX2(2,0,8)]+=-C_1_0*w[39][l] + C_1_2*w[44][l] + C_1_5*w[45][l] - C_1_7*w[31][l] + tmp138 + tmp139 + tmp140 + tmp141 + tmp142 + tmp143 + tmp144 + tmp145 + tmp146 + tmp147;
                    EM_S[INDEX2(3,0,8)]+=tmp40 + tmp41 + tmp42 + tmp43 + tmp44 + tmp45 + tmp46 + tmp47 + tmp48 + tmp49;
                    EM_S[INDEX2(4,0,8)]+=-C_2_0*w[43][l] - C_2_3*w[41][l] - C_2_4*w[40][l] - C_2_7*w[32][l] + tmp114 + tmp115 + tmp116 + tmp117 + tmp118 + tmp119 + tmp120 + tmp121 + tmp2 + tmp3;
                    EM_S[INDEX2(5,0,8)]+=tmp191 + tmp192 + tmp193 + tmp194 + tmp198 + tmp199 + tmp214 + tmp215 + tmp216 + tmp217;
                    EM_S[INDEX2(6,0,8)]+=tmp228 + tmp229 + tmp230 + tmp231 + tmp232 + tmp233 + tmp234 + tmp235 + tmp236 + tmp237;
                    EM_S[INDEX2(7,0,8)]+=tmp258 + tmp259 + tmp80 + tmp81 + tmp94 + tmp95;
                    EM_S[INDEX2(0,1,8)]+= C_0_0*w[47][l] + C_0_1*w[38][l] + C_0_6*w[30][l] + C_0_7*w[46][l] + tmp128 + tmp129 + tmp242 + tmp243 + tmp244 + tmp245 + tmp246 + tmp247 + tmp248 + tmp249;
                    EM_S[INDEX2(1,1,8)]+= C_0_0*w[38][l] + C_0_1*w[47][l] + C_0_6*w[46][l] + C_0_7*w[30][l] + C_1_1*w[44][l] - C_1_3*w[39][l] - C_1_4*w[31][l] + C_1_6*w[45][l] - C_2_1*w[40][l] - C_2_2*w[32][l] - C_2_5*w[43][l] - C_2_6*w[41][l] + tmp152 + tmp155 + tmp164 + tmp165 + tmp168 + tmp169;
                    EM_S[INDEX2(2,1,8)]+=tmp100 + tmp101 + tmp102 + tmp41 + tmp42 + tmp43 + tmp44 + tmp45 + tmp49 + tmp96;
                    EM_S[INDEX2(3,1,8)]+=-C_1_1*w[39][l] + C_1_3*w[44][l] + C_1_4*w[45][l] - C_1_6*w[31][l] + tmp20 + tmp21 + tmp22 + tmp24 + tmp25 + tmp26 + tmp27 + tmp28 + tmp33 + tmp39;
                    EM_S[INDEX2(4,1,8)]+=tmp190 + tmp191 + tmp192 + tmp193 + tmp194 + tmp195 + tmp196 + tmp197 + tmp198 + tmp199;
                    EM_S[INDEX2(5,1,8)]+=-C_2_1*w[43][l] - C_2_2*w[41][l] - C_2_5*w[40][l] - C_2_6*w[32][l] + tmp82 + tmp83 + tmp84 + tmp85 + tmp86 + tmp87 + tmp88 + tmp89 + tmp90 + tmp91;
                    EM_S[INDEX2(6,1,8)]+=tmp258 + tmp259 + tmp63 + tmp64 + tmp94 + tmp95;
                    EM_S[INDEX2(7,1,8)]+=tmp181 + tmp182 + tmp184 + tmp186 + tmp187 + tmp188 + tmp218 + tmp219 + tmp220 + tmp221;
                    EM_S[INDEX2(0,2,8)]+=-C_1_0*w[44][l] + C_1_2*w[39][l] + C_1_5*w[31][l] - C_1_7*w[45][l] + tmp138 + tmp139 + tmp140 + tmp142 + tmp143 + tmp144 + tmp145 + tmp146 + tmp173 + tmp179;
                    EM_S[INDEX2(1,2,8)]+=tmp103 + tmp40 + tmp42 + tmp43 + tmp46 + tmp47 + tmp48 + tmp97 + tmp98 + tmp99;
                    EM_S[INDEX2(2,2,8)]+=-C_0_2*w[47][l] - C_0_3*w[38][l] - C_0_4*w[30][l] - C_0_5*w[46][l] + C_1_0*w[39][l] - C_1_2*w[44][l] - C_1_5*w[45][l] + C_1_7*w[31][l] - C_2_1*w[32][l] - C_2_2*w[40][l] - C_2_5*w[41][l] - C_2_6*w[43][l] + tmp153 + tmp154 + tmp164 + tmp165 + tmp166 + tmp167;
                    EM_S[INDEX2(3,2,8)]+=-C_0_2*w[38][l] - C_0_3*w[47][l] - C_0_4*w[46][l] - C_0_5*w[30][l] + tmp200 + tmp201 + tmp202 + tmp203 + tmp204 + tmp205 + tmp206 + tmp207 + tmp56 + tmp57;
                    EM_S[INDEX2(4,2,8)]+=tmp229 + tmp230 + tmp232 + tmp234 + tmp235 + tmp236 + tmp238 + tmp239 + tmp240 + tmp241;
                    EM_S[INDEX2(5,2,8)]+=tmp258 + tmp259 + tmp62 + tmp65 + tmp80 + tmp81;
                    EM_S[INDEX2(6,2,8)]+=-C_2_1*w[41][l] - C_2_2*w[43][l] - C_2_5*w[32][l] - C_2_6*w[40][l] + tmp70 + tmp71 + tmp74 + tmp75 + tmp76 + tmp77 + tmp78 + tmp79 + tmp84 + tmp85;
                    EM_S[INDEX2(7,2,8)]+=tmp104 + tmp105 + tmp106 + tmp107 + tmp108 + tmp109 + tmp110 + tmp111 + tmp112 + tmp113;
                    EM_S[INDEX2(0,3,8)]+=tmp100 + tmp101 + tmp102 + tmp103 + tmp42 + tmp43 + tmp96 + tmp97 + tmp98 + tmp99;
                    EM_S[INDEX2(1,3,8)]+=-C_1_1*w[44][l] + C_1_3*w[39][l] + C_1_4*w[31][l] - C_1_6*w[45][l] + tmp20 + tmp21 + tmp22 + tmp23 + tmp24 + tmp25 + tmp26 + tmp27 + tmp28 + tmp29;
                    EM_S[INDEX2(2,3,8)]+= C_0_2*w[47][l] + C_0_3*w[38][l] + C_0_4*w[30][l] + C_0_5*w[46][l] + tmp200 + tmp201 + tmp202 + tmp203 + tmp204 + tmp205 + tmp206 + tmp207 + tmp222 + tmp223;
                    EM_S[INDEX2(3,3,8)]+= C_0_2*w[38][l] + C_0_3*w[47][l] + C_0_4*w[46][l] + C_0_5*w[30][l] + C_1_1*w[39][l] - C_1_3*w[44][l] - C_1_4*w[45][l] + C_1_6*w[31][l] - C_2_0*w[32][l] - C_2_3*w[40][l] - C_2_4*w[41][l] - C_2_7*w[43][l] + tmp132 + tmp133 + tmp134 + tmp135 + tmp136 + tmp137;
                    EM_S[INDEX2(4,3,8)]+=tmp258 + tmp259 + tmp62 + tmp63 + tmp64 + tmp65;
                    EM_S[INDEX2(5,3,8)]+=tmp180 + tmp181 + tmp182 + tmp183 + tmp184 + tmp185 + tmp186 + tmp187 + tmp188 + tmp189;
                    EM_S[INDEX2(6,3,8)]+=tmp105 + tmp106 + tmp107 + tmp108 + tmp112 + tmp113 + tmp156 + tmp157 + tmp158 + tmp159;
                    EM_S[INDEX2(7,3,8)]+=-C_2_0*w[41][l] - C_2_3*w[43][l] - C_2_4*w[32][l] - C_2_7*w[40][l] + tmp0 + tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6 + tmp7 + tmp8 + tmp9;
                    EM_S[INDEX2(0,4,8)]+= C_2_0*w[40][l] + C_2_3*w[32][l] + C_2_4*w[43][l] + C_2_7*w[41][l] + tmp114 + tmp115 + tmp116 + tmp117 + tmp118 + tmp119 + tmp120 + tmp121 + tmp92 + tmp93;
                    EM_S[INDEX2(1,4,8)]+=tmp193 + tmp199 + tmp214 + tmp215 + tmp216 + tmp217 + tmp224 + tmp225 + tmp226 + tmp227;
                    EM_S[INDEX2(2,4,8)]+=tmp228 + tmp231 + tmp233 + tmp234 + tmp235 + tmp237 + tmp260 + tmp261 + tmp262 + tmp263;
                    EM_S[INDEX2(3,4,8)]+=tmp60 + tmp61 + tmp80 + tmp81 + tmp94 + tmp95;
                    EM_S[INDEX2(4,4,8)]+=-C_0_2*w[30][l] - C_0_3*w[46][l] - C_0_4*w[47][l] - C_0_5*w[38][l] - C_1_1*w[31][l] + C_1_3*w[45][l] + C_1_4*w[44][l] - C_1_6*w[39][l] + C_2_0*w[43][l] + C_2_3*w[41][l] + C_2_4*w[40][l] + C_2_7*w[32][l] + tmp150 + tmp151 + tmp152 + tmp153 + tmp154 + tmp155;
                    EM_S[INDEX2(5,4,8)]+=-C_0_2*w[46][l] - C_0_3*w[30][l] - C_0_4*w[38][l] - C_0_5*w[47][l] + tmp50 + tmp51 + tmp52 + tmp53 + tmp54 + tmp55 + tmp56 + tmp57 + tmp58 + tmp59;
                    EM_S[INDEX2(6,4,8)]+= C_1_1*w[45][l] - C_1_3*w[31][l] - C_1_4*w[39][l] + C_1_6*w[44][l] + tmp30 + tmp31 + tmp32 + tmp33 + tmp34 + tmp35 + tmp36 + tmp37 + tmp38 + tmp39;
                    EM_S[INDEX2(7,4,8)]+=tmp12 + tmp13 + tmp250 + tmp251 + tmp252 + tmp253 + tmp66 + tmp67 + tmp68 + tmp69;
                    EM_S[INDEX2(0,5,8)]+=tmp190 + tmp193 + tmp195 + tmp196 + tmp197 + tmp199 + tmp224 + tmp225 + tmp226 + tmp227;
                    EM_S[INDEX2(1,5,8)]+= C_2_1*w[40][l] + C_2_2*w[32][l] + C_2_5*w[43][l] + C_2_6*w[41][l] + tmp72 + tmp73 + tmp82 + tmp83 + tmp86 + tmp87 + tmp88 + tmp89 + tmp90 + tmp91;
                    EM_S[INDEX2(2,5,8)]+=tmp60 + tmp61 + tmp63 + tmp64 + tmp94 + tmp95;
                    EM_S[INDEX2(3,5,8)]+=tmp186 + tmp187 + tmp218 + tmp219 + tmp220 + tmp221 + tmp254 + tmp255 + tmp256 + tmp257;
                    EM_S[INDEX2(4,5,8)]+= C_0_2*w[30][l] + C_0_3*w[46][l] + C_0_4*w[47][l] + C_0_5*w[38][l] + tmp222 + tmp223 + tmp50 + tmp51 + tmp52 + tmp53 + tmp54 + tmp55 + tmp58 + tmp59;
                    EM_S[INDEX2(5,5,8)]+= C_0_2*w[46][l] + C_0_3*w[30][l] + C_0_4*w[38][l] + C_0_5*w[47][l] - C_1_0*w[31][l] + C_1_2*w[45][l] + C_1_5*w[44][l] - C_1_7*w[39][l] + C_2_1*w[43][l] + C_2_2*w[41][l] + C_2_5*w[40][l] + C_2_6*w[32][l] + tmp135 + tmp136 + tmp208 + tmp211 + tmp212 + tmp213;
                    EM_S[INDEX2(6,5,8)]+=tmp10 + tmp12 + tmp13 + tmp16 + tmp17 + tmp18 + tmp250 + tmp251 + tmp252 + tmp253;
                    EM_S[INDEX2(7,5,8)]+= C_1_0*w[45][l] - C_1_2*w[31][l] - C_1_5*w[39][l] + C_1_7*w[44][l] + tmp141 + tmp147 + tmp170 + tmp171 + tmp172 + tmp174 + tmp175 + tmp176 + tmp177 + tmp178;
                    EM_S[INDEX2(0,6,8)]+=tmp234 + tmp235 + tmp238 + tmp239 + tmp240 + tmp241 + tmp260 + tmp261 + tmp262 + tmp263;
                    EM_S[INDEX2(1,6,8)]+=tmp60 + tmp61 + tmp62 + tmp65 + tmp80 + tmp81;
                    EM_S[INDEX2(2,6,8)]+= C_2_1*w[32][l] + C_2_2*w[40][l] + C_2_5*w[41][l] + C_2_6*w[43][l] + tmp70 + tmp71 + tmp72 + tmp73 + tmp74 + tmp75 + tmp76 + tmp77 + tmp78 + tmp79;
                    EM_S[INDEX2(3,6,8)]+=tmp104 + tmp107 + tmp109 + tmp110 + tmp111 + tmp113 + tmp160 + tmp161 + tmp162 + tmp163;
                    EM_S[INDEX2(4,6,8)]+= C_1_1*w[31][l] - C_1_3*w[45][l] - C_1_4*w[44][l] + C_1_6*w[39][l] + tmp23 + tmp29 + tmp30 + tmp31 + tmp32 + tmp34 + tmp35 + tmp36 + tmp37 + tmp38;
                    EM_S[INDEX2(5,6,8)]+=tmp11 + tmp12 + tmp13 + tmp14 + tmp15 + tmp19 + tmp66 + tmp67 + tmp68 + tmp69;
                    EM_S[INDEX2(6,6,8)]+=-C_0_0*w[30][l] - C_0_1*w[46][l] - C_0_6*w[47][l] - C_0_7*w[38][l] - C_1_1*w[45][l] + C_1_3*w[31][l] + C_1_4*w[39][l] - C_1_6*w[44][l] + C_2_1*w[41][l] + C_2_2*w[43][l] + C_2_5*w[32][l] + C_2_6*w[40][l] + tmp134 + tmp137 + tmp209 + tmp210 + tmp212 + tmp213;
                    EM_S[INDEX2(7,6,8)]+=-C_0_0*w[46][l] - C_0_1*w[30][l] - C_0_6*w[38][l] - C_0_7*w[47][l] + tmp122 + tmp123 + tmp124 + tmp125 + tmp126 + tmp127 + tmp130 + tmp131 + tmp148 + tmp149;
                    EM_S[INDEX2(0,7,8)]+=tmp60 + tmp61 + tmp62 + tmp63 + tmp64 + tmp65;
                    EM_S[INDEX2(1,7,8)]+=tmp180 + tmp183 + tmp185 + tmp186 + tmp187 + tmp189 + tmp254 + tmp255 + tmp256 + tmp257;
                    EM_S[INDEX2(2,7,8)]+=tmp107 + tmp113 + tmp156 + tmp157 + tmp158 + tmp159 + tmp160 + tmp161 + tmp162 + tmp163;
                    EM_S[INDEX2(3,7,8)]+= C_2_0*w[32][l] + C_2_3*w[40][l] + C_2_4*w[41][l] + C_2_7*w[43][l] + tmp0 + tmp1 + tmp4 + tmp5 + tmp6 + tmp7 + tmp8 + tmp9 + tmp92 + tmp93;
                    EM_S[INDEX2(4,7,8)]+=tmp10 + tmp11 + tmp12 + tmp13 + tmp14 + tmp15 + tmp16 + tmp17 + tmp18 + tmp19;
                    EM_S[INDEX2(5,7,8)]+= C_1_0*w[31][l] - C_1_2*w[45][l] - C_1_5*w[44][l] + C_1_7*w[39][l] + tmp170 + tmp171 + tmp172 + tmp173 + tmp174 + tmp175 + tmp176 + tmp177 + tmp178 + tmp179;
                    EM_S[INDEX2(6,7,8)]+= C_0_0*w[30][l] + C_0_1*w[46][l] + C_0_6*w[47][l] + C_0_7*w[38][l] + tmp122 + tmp123 + tmp124 + tmp125 + tmp126 + tmp127 + tmp128 + tmp129 + tmp130 + tmp131;
                    EM_S[INDEX2(7,7,8)]+= C_0_0*w[46][l] + C_0_1*w[30][l] + C_0_6*w[38][l] + C_0_7*w[47][l] - C_1_0*w[45][l] + C_1_2*w[31][l] + C_1_5*w[39][l] - C_1_7*w[44][l] + C_2_0*w[41][l] + C_2_3*w[43][l] + C_2_4*w[32][l] + C_2_7*w[40][l] + tmp150 + tmp151 + tmp166 + tmp167 + tmp168 + tmp169;
                } else { // constant data
                    const Scalar wC0 = C_p[0]*w[53][l];
                    const Scalar wC1 = C_p[1]*w[51][l];
                    const Scalar wC2 = C_p[2]*w[50][l];
                    EM_S[INDEX2(0,0,8)]+= 4.*wC0 + 4.*wC1 - 4.*wC2;
                    EM_S[INDEX2(1,0,8)]+= 4.*wC0 + 2.*wC1 - 2.*wC2;
                    EM_S[INDEX2(2,0,8)]+= 2.*wC0 + 4.*wC1 - 2.*wC2;
                    EM_S[INDEX2(3,0,8)]+= 2.*wC0 + 2.*wC1 -    wC2;
                    EM_S[INDEX2(4,0,8)]+= 2.*wC0 + 2.*wC1 - 4.*wC2;
                    EM_S[INDEX2(5,0,8)]+= 2.*wC0 +    wC1 - 2.*wC2;
                    EM_S[INDEX2(6,0,8)]+=    wC0 + 2.*wC1 - 2.*wC2;
                    EM_S[INDEX2(7,0,8)]+=    wC0 +    wC1 -    wC2;
                    EM_S[INDEX2(0,1,8)]+=-4.*wC0 + 2.*wC1 - 2.*wC2;
                    EM_S[INDEX2(1,1,8)]+=-4.*wC0 + 4.*wC1 - 4.*wC2;
                    EM_S[INDEX2(2,1,8)]+=-2.*wC0 + 2.*wC1 -    wC2;
                    EM_S[INDEX2(3,1,8)]+=-2.*wC0 + 4.*wC1 - 2.*wC2;
                    EM_S[INDEX2(4,1,8)]+=-2.*wC0 +    wC1 - 2.*wC2;
                    EM_S[INDEX2(5,1,8)]+=-2.*wC0 + 2.*wC1 - 4.*wC2;
                    EM_S[INDEX2(6,1,8)]+=   -wC0 +    wC1 -    wC2;
                    EM_S[INDEX2(7,1,8)]+=   -wC0 + 2.*wC1 - 2.*wC2;
                    EM_S[INDEX2(0,2,8)]+= 2.*wC0 - 4.*wC1 - 2.*wC2;
                    EM_S[INDEX2(1,2,8)]+= 2.*wC0 - 2.*wC1 -    wC2;
                    EM_S[INDEX2(2,2,8)]+= 4.*wC0 - 4.*wC1 - 4.*wC2;
                    EM_S[INDEX2(3,2,8)]+= 4.*wC0 - 2.*wC1 - 2.*wC2;
                    EM_S[INDEX2(4,2,8)]+=    wC0 - 2.*wC1 - 2.*wC2;
                    EM_S[INDEX2(5,2,8)]+=    wC0 -    wC1 -    wC2;
                    EM_S[INDEX2(6,2,8)]+= 2.*wC0 - 2.*wC1 - 4.*wC2;
                    EM_S[INDEX2(7,2,8)]+= 2.*wC0 -    wC1 - 2.*wC2;
                    EM_S[INDEX2(0,3,8)]+=-2.*wC0 - 2.*wC1 -    wC2;
                    EM_S[INDEX2(1,3,8)]+=-2.*wC0 - 4.*wC1 - 2.*wC2;
                    EM_S[INDEX2(2,3,8)]+=-4.*wC0 - 2.*wC1 - 2.*wC2;
                    EM_S[INDEX2(3,3,8)]+=-4.*wC0 - 4.*wC1 - 4.*wC2;
                    EM_S[INDEX2(4,3,8)]+=   -wC0 -    wC1 -    wC2;
                    EM_S[INDEX2(5,3,8)]+=   -wC0 - 2.*wC1 - 2.*wC2;
                    EM_S[INDEX2(6,3,8)]+=-2.*wC0 -    wC1 - 2.*wC2;
                    EM_S[INDEX2(7,3,8)]+=-2.*wC0 - 2.*wC1 - 4.*wC2;
                    EM_S[INDEX2(0,4,8)]+= 2.*wC0 + 2.*wC1 + 4.*wC2;
                    EM_S[INDEX2(1,4,8)]+= 2.*wC0 +    wC1 + 2.*wC2;
                    EM_S[INDEX2(2,4,8)]+=    wC0 + 2.*wC1 + 2.*wC2;
                    EM_S[INDEX2(3,4,8)]+=    wC0 +    wC1 +    wC2;
                    EM_S[INDEX2(4,4,8)]+= 4.*wC0 + 4.*wC1 + 4.*wC2;
                    EM_S[INDEX2(5,4,8)]+= 4.*wC0 + 2.*wC1 + 2.*wC2;
                    EM_S[INDEX2(6,4,8)]+= 2.*wC0 + 4.*wC1 + 2.*wC2;
                    EM_S[INDEX2(7,4,8)]+= 2.*wC0 + 2.*wC1 +    wC2;
                    EM_S[INDEX2(0,5,8)]+=-2.*wC0 +    wC1 + 2.*wC2;
                    EM_S[INDEX2(1,5,8)]+=-2.*wC0 + 2.*wC1 + 4.*wC2;
                    EM_S[INDEX2(2,5,8)]+=   -wC0 +    wC1 +    wC2;
                    EM_S[INDEX2(3,5,8)]+=   -wC0 + 2.*wC1 + 2.*wC2;
                    EM_S[INDEX2(4,5,8)]+=-4.*wC0 + 2.*wC1 + 2.*wC2;
                    EM_S[INDEX2(5,5,8)]+=-4.*wC0 + 4.*wC1 + 4.*wC2;
                    EM_S[INDEX2(6,5,8)]+=-2.*wC0 + 2.*wC1 +    wC2;
                    EM_S[INDEX2(7,5,8)]+=-2.*wC0 + 4.*wC1 + 2.*wC2;
                    EM_S[INDEX2(0,6,8)]+=    wC0 - 2.*wC1 + 2.*wC2;
                    EM_S[INDEX2(1,6,8)]+=    wC0 -    wC1 +    wC2;
                    EM_S[INDEX2(2,6,8)]+= 2.*wC0 - 2.*wC1 + 4.*wC2;
                    EM_S[INDEX2(3,6,8)]+= 2.*wC0 -    wC1 + 2.*wC2;
                    EM_S[INDEX2(4,6,8)]+= 2.*wC0 - 4.*wC1 + 2.*wC2;
                    EM_S[INDEX2(5,6,8)]+= 2.*wC0 - 2.*wC1 +    wC2;
                    EM_S[INDEX2(6,6,8)]+= 4.*wC0 - 4.*wC1 + 4.*wC2;
                    EM_S[INDEX2(7,6,8)]+= 4.*wC0 - 2.*wC1 + 2.*wC2;
                    EM_S[INDEX2(0,7,8)]+=   -wC0 -    wC1 +    wC2;
                    EM_S[INDEX2(1,7,8)]+=   -wC0 - 2.*wC1 + 2.*wC2;
                    EM_S[INDEX2(2,7,8)]+=-2.*wC0 -    wC1 + 2.*wC2;
                    EM_S[INDEX2(3,7,8)]+=-2.*wC0 - 2.*wC1 + 4.*wC2;
                    EM_S[INDEX2(4,7,8)]+=-2.*wC0 - 2.*wC1 +    wC2;
                    EM_S[INDEX2(5,7,8)]+=-2.*wC0 - 4.*wC1 + 2.*wC2;
                    EM_S[INDEX2(6,7,8)]+=-4.*wC0 - 2.*wC1 + 2.*wC2;
                    EM_S[INDEX2(7,7,8)]+=-4.*wC0 - 4.*wC1 + 4.*wC2;
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
                    const Scalar D_4 = D_p[4];
                    const Scalar D_5 = D_p[5];
                    const Scalar D_6 = D_p[6];
                    const Scalar D_7 = D_p[7];
                    const Scalar tmp0  = w[54][l]*(D_0 + D_4);
                    const Scalar tmp1  = w[55][l]*(D_1 + D_2 + D_5 + D_6);
                    const Scalar tmp2  = w[56][l]*(D_3 + D_7);
                    const Scalar tmp3  = w[57][l]*(D_0 + D_1 + D_2 + D_3);
                    const Scalar tmp4  = w[58][l]*(D_4 + D_5 + D_6 + D_7);
                    const Scalar tmp5  = w[54][l]*(D_4 + D_6);
                    const Scalar tmp6  = w[55][l]*(D_0 + D_2 + D_5 + D_7);
                    const Scalar tmp7  = w[56][l]*(D_1 + D_3);
                    const Scalar tmp8  = w[54][l]*(D_1 + D_3);
                    const Scalar tmp9  = w[56][l]*(D_4 + D_6);
                    const Scalar tmp10 = w[57][l]*(D_4 + D_5 + D_6 + D_7);
                    const Scalar tmp11 = w[58][l]*(D_0 + D_1 + D_2 + D_3);
                    const Scalar tmp12 = w[54][l]*(D_2 + D_3);
                    const Scalar tmp13 = w[55][l]*(D_0 + D_1 + D_6 + D_7);
                    const Scalar tmp14 = w[56][l]*(D_4 + D_5);
                    const Scalar tmp15 = w[55][l]*(D_0 + D_1 + D_2 + D_3 + D_4 + D_5 + D_6 + D_7);
                    const Scalar tmp16 = w[54][l]*(D_1 + D_5);
                    const Scalar tmp17 = w[55][l]*(D_0 + D_3 + D_4 + D_7);
                    const Scalar tmp18 = w[56][l]*(D_2 + D_6);
                    const Scalar tmp19 = w[54][l]*(D_2 + D_6);
                    const Scalar tmp20 = w[56][l]*(D_1 + D_5);
                    const Scalar tmp21 = w[57][l]*(D_0 + D_1 + D_4 + D_5);
                    const Scalar tmp22 = w[58][l]*(D_2 + D_3 + D_6 + D_7);
                    const Scalar tmp23 = w[54][l]*(D_3 + D_7);
                    const Scalar tmp24 = w[56][l]*(D_0 + D_4);
                    const Scalar tmp25 = w[54][l]*(D_0 + D_1);
                    const Scalar tmp26 = w[55][l]*(D_2 + D_3 + D_4 + D_5);
                    const Scalar tmp27 = w[56][l]*(D_6 + D_7);
                    const Scalar tmp28 = w[57][l]*(D_0 + D_5 + D_6);
                    const Scalar tmp29 = w[58][l]*(D_1 + D_2 + D_7);
                    const Scalar tmp30 = w[54][l]*(D_5 + D_7);
                    const Scalar tmp31 = w[55][l]*(D_1 + D_3 + D_4 + D_6);
                    const Scalar tmp32 = w[56][l]*(D_0 + D_2);
                    const Scalar tmp33 = w[57][l]*(D_1 + D_2 + D_7);
                    const Scalar tmp34 = w[58][l]*(D_0 + D_5 + D_6);
                    const Scalar tmp35 = w[57][l]*(D_1 + D_4 + D_7);
                    const Scalar tmp36 = w[58][l]*(D_0 + D_3 + D_6);
                    const Scalar tmp37 = w[57][l]*(D_1 + D_2 + D_4);
                    const Scalar tmp38 = w[58][l]*(D_3 + D_5 + D_6);
                    const Scalar tmp39 = w[54][l]*(D_0 + D_2);
                    const Scalar tmp40 = w[56][l]*(D_5 + D_7);
                    const Scalar tmp41 = w[57][l]*(D_0 + D_2 + D_4 + D_6);
                    const Scalar tmp42 = w[58][l]*(D_1 + D_3 + D_5 + D_7);
                    const Scalar tmp43 = w[57][l]*(D_2 + D_3 + D_6 + D_7);
                    const Scalar tmp44 = w[58][l]*(D_0 + D_1 + D_4 + D_5);
                    const Scalar tmp45 = w[57][l]*(D_2 + D_4 + D_7);
                    const Scalar tmp46 = w[58][l]*(D_0 + D_3 + D_5);
                    const Scalar tmp47 = w[54][l]*(D_4 + D_5);
                    const Scalar tmp48 = w[56][l]*(D_2 + D_3);
                    const Scalar tmp49 = w[57][l]*(D_3 + D_5 + D_6);
                    const Scalar tmp50 = w[58][l]*(D_1 + D_2 + D_4);
                    const Scalar tmp51 = w[57][l]*(D_0 + D_3 + D_5);
                    const Scalar tmp52 = w[58][l]*(D_2 + D_4 + D_7);
                    const Scalar tmp53 = w[57][l]*(D_0 + D_3 + D_6);
                    const Scalar tmp54 = w[58][l]*(D_1 + D_4 + D_7);
                    const Scalar tmp55 = w[57][l]*(D_1 + D_3 + D_5 + D_7);
                    const Scalar tmp56 = w[58][l]*(D_0 + D_2 + D_4 + D_6);
                    const Scalar tmp57 = w[54][l]*(D_6 + D_7);
                    const Scalar tmp58 = w[56][l]*(D_0 + D_1);
                    EM_S[INDEX2(0,0,8)]+=D_0*w[59][l] + D_7*w[60][l] + tmp49 + tmp50;
                    EM_S[INDEX2(1,0,8)]+=tmp26 + tmp57 + tmp58;
                    EM_S[INDEX2(2,0,8)]+=tmp30 + tmp31 + tmp32;
                    EM_S[INDEX2(3,0,8)]+=tmp10 + tmp11;
                    EM_S[INDEX2(4,0,8)]+=tmp1 + tmp23 + tmp24;
                    EM_S[INDEX2(5,0,8)]+=tmp43 + tmp44;
                    EM_S[INDEX2(6,0,8)]+=tmp55 + tmp56;
                    EM_S[INDEX2(7,0,8)]+=tmp15;
                    EM_S[INDEX2(0,1,8)]+=tmp26 + tmp57 + tmp58;
                    EM_S[INDEX2(1,1,8)]+=D_1*w[59][l] + D_6*w[60][l] + tmp45 + tmp46;
                    EM_S[INDEX2(2,1,8)]+=tmp10 + tmp11;
                    EM_S[INDEX2(3,1,8)]+=tmp5 + tmp6 + tmp7;
                    EM_S[INDEX2(4,1,8)]+=tmp43 + tmp44;
                    EM_S[INDEX2(5,1,8)]+=tmp17 + tmp19 + tmp20;
                    EM_S[INDEX2(6,1,8)]+=tmp15;
                    EM_S[INDEX2(7,1,8)]+=tmp41 + tmp42;
                    EM_S[INDEX2(0,2,8)]+=tmp30 + tmp31 + tmp32;
                    EM_S[INDEX2(1,2,8)]+=tmp10 + tmp11;
                    EM_S[INDEX2(2,2,8)]+=D_2*w[59][l] + D_5*w[60][l] + tmp35 + tmp36;
                    EM_S[INDEX2(3,2,8)]+=tmp13 + tmp47 + tmp48;
                    EM_S[INDEX2(4,2,8)]+=tmp55 + tmp56;
                    EM_S[INDEX2(5,2,8)]+=tmp15;
                    EM_S[INDEX2(6,2,8)]+=tmp16 + tmp17 + tmp18;
                    EM_S[INDEX2(7,2,8)]+=tmp21 + tmp22;
                    EM_S[INDEX2(0,3,8)]+=tmp10 + tmp11;
                    EM_S[INDEX2(1,3,8)]+=tmp5 + tmp6 + tmp7;
                    EM_S[INDEX2(2,3,8)]+=tmp13 + tmp47 + tmp48;
                    EM_S[INDEX2(3,3,8)]+=D_3*w[59][l] + D_4*w[60][l] + tmp28 + tmp29;
                    EM_S[INDEX2(4,3,8)]+=tmp15;
                    EM_S[INDEX2(5,3,8)]+=tmp41 + tmp42;
                    EM_S[INDEX2(6,3,8)]+=tmp21 + tmp22;
                    EM_S[INDEX2(7,3,8)]+=tmp0 + tmp1 + tmp2;
                    EM_S[INDEX2(0,4,8)]+=tmp1 + tmp23 + tmp24;
                    EM_S[INDEX2(1,4,8)]+=tmp43 + tmp44;
                    EM_S[INDEX2(2,4,8)]+=tmp55 + tmp56;
                    EM_S[INDEX2(3,4,8)]+=tmp15;
                    EM_S[INDEX2(4,4,8)]+=D_3*w[60][l] + D_4*w[59][l] + tmp33 + tmp34;
                    EM_S[INDEX2(5,4,8)]+=tmp12 + tmp13 + tmp14;
                    EM_S[INDEX2(6,4,8)]+=tmp6 + tmp8 + tmp9;
                    EM_S[INDEX2(7,4,8)]+=tmp3 + tmp4;
                    EM_S[INDEX2(0,5,8)]+=tmp43 + tmp44;
                    EM_S[INDEX2(1,5,8)]+=tmp17 + tmp19 + tmp20;
                    EM_S[INDEX2(2,5,8)]+=tmp15;
                    EM_S[INDEX2(3,5,8)]+=tmp41 + tmp42;
                    EM_S[INDEX2(4,5,8)]+=tmp12 + tmp13 + tmp14;
                    EM_S[INDEX2(5,5,8)]+=D_2*w[60][l] + D_5*w[59][l] + tmp53 + tmp54;
                    EM_S[INDEX2(6,5,8)]+=tmp3 + tmp4;
                    EM_S[INDEX2(7,5,8)]+=tmp31 + tmp39 + tmp40;
                    EM_S[INDEX2(0,6,8)]+=tmp55 + tmp56;
                    EM_S[INDEX2(1,6,8)]+=tmp15;
                    EM_S[INDEX2(2,6,8)]+=tmp16 + tmp17 + tmp18;
                    EM_S[INDEX2(3,6,8)]+=tmp21 + tmp22;
                    EM_S[INDEX2(4,6,8)]+=tmp6 + tmp8 + tmp9;
                    EM_S[INDEX2(5,6,8)]+=tmp3 + tmp4;
                    EM_S[INDEX2(6,6,8)]+=D_1*w[60][l] + D_6*w[59][l] + tmp51 + tmp52;
                    EM_S[INDEX2(7,6,8)]+=tmp25 + tmp26 + tmp27;
                    EM_S[INDEX2(0,7,8)]+=tmp15;
                    EM_S[INDEX2(1,7,8)]+=tmp41 + tmp42;
                    EM_S[INDEX2(2,7,8)]+=tmp21 + tmp22;
                    EM_S[INDEX2(3,7,8)]+=tmp0 + tmp1 + tmp2;
                    EM_S[INDEX2(4,7,8)]+=tmp3 + tmp4;
                    EM_S[INDEX2(5,7,8)]+=tmp31 + tmp39 + tmp40;
                    EM_S[INDEX2(6,7,8)]+=tmp25 + tmp26 + tmp27;
                    EM_S[INDEX2(7,7,8)]+=D_0*w[60][l] + D_7*w[59][l] + tmp37 + tmp38;
                } else { // constant data
                    const Scalar wD0 = 8.*D_p[0]*w[55][l];
                    EM_S[INDEX2(0,0,8)]+=8.*wD0;
                    EM_S[INDEX2(1,0,8)]+=4.*wD0;
                    EM_S[INDEX2(2,0,8)]+=4.*wD0;
                    EM_S[INDEX2(3,0,8)]+=2.*wD0;
                    EM_S[INDEX2(4,0,8)]+=4.*wD0;
                    EM_S[INDEX2(5,0,8)]+=2.*wD0;
                    EM_S[INDEX2(6,0,8)]+=2.*wD0;
                    EM_S[INDEX2(7,0,8)]+=   wD0;
                    EM_S[INDEX2(0,1,8)]+=4.*wD0;
                    EM_S[INDEX2(1,1,8)]+=8.*wD0;
                    EM_S[INDEX2(2,1,8)]+=2.*wD0;
                    EM_S[INDEX2(3,1,8)]+=4.*wD0;
                    EM_S[INDEX2(4,1,8)]+=2.*wD0;
                    EM_S[INDEX2(5,1,8)]+=4.*wD0;
                    EM_S[INDEX2(6,1,8)]+=   wD0;
                    EM_S[INDEX2(7,1,8)]+=2.*wD0;
                    EM_S[INDEX2(0,2,8)]+=4.*wD0;
                    EM_S[INDEX2(1,2,8)]+=2.*wD0;
                    EM_S[INDEX2(2,2,8)]+=8.*wD0;
                    EM_S[INDEX2(3,2,8)]+=4.*wD0;
                    EM_S[INDEX2(4,2,8)]+=2.*wD0;
                    EM_S[INDEX2(5,2,8)]+=   wD0;
                    EM_S[INDEX2(6,2,8)]+=4.*wD0;
                    EM_S[INDEX2(7,2,8)]+=2.*wD0;
                    EM_S[INDEX2(0,3,8)]+=2.*wD0;
                    EM_S[INDEX2(1,3,8)]+=4.*wD0;
                    EM_S[INDEX2(2,3,8)]+=4.*wD0;
                    EM_S[INDEX2(3,3,8)]+=8.*wD0;
                    EM_S[INDEX2(4,3,8)]+=   wD0;
                    EM_S[INDEX2(5,3,8)]+=2.*wD0;
                    EM_S[INDEX2(6,3,8)]+=2.*wD0;
                    EM_S[INDEX2(7,3,8)]+=4.*wD0;
                    EM_S[INDEX2(0,4,8)]+=4.*wD0;
                    EM_S[INDEX2(1,4,8)]+=2.*wD0;
                    EM_S[INDEX2(2,4,8)]+=2.*wD0;
                    EM_S[INDEX2(3,4,8)]+=   wD0;
                    EM_S[INDEX2(4,4,8)]+=8.*wD0;
                    EM_S[INDEX2(5,4,8)]+=4.*wD0;
                    EM_S[INDEX2(6,4,8)]+=4.*wD0;
                    EM_S[INDEX2(7,4,8)]+=2.*wD0;
                    EM_S[INDEX2(0,5,8)]+=2.*wD0;
                    EM_S[INDEX2(1,5,8)]+=4.*wD0;
                    EM_S[INDEX2(2,5,8)]+=   wD0;
                    EM_S[INDEX2(3,5,8)]+=2.*wD0;
                    EM_S[INDEX2(4,5,8)]+=4.*wD0;
                    EM_S[INDEX2(5,5,8)]+=8.*wD0;
                    EM_S[INDEX2(6,5,8)]+=2.*wD0;
                    EM_S[INDEX2(7,5,8)]+=4.*wD0;
                    EM_S[INDEX2(0,6,8)]+=2.*wD0;
                    EM_S[INDEX2(1,6,8)]+=   wD0;
                    EM_S[INDEX2(2,6,8)]+=4.*wD0;
                    EM_S[INDEX2(3,6,8)]+=2.*wD0;
                    EM_S[INDEX2(4,6,8)]+=4.*wD0;
                    EM_S[INDEX2(5,6,8)]+=2.*wD0;
                    EM_S[INDEX2(6,6,8)]+=8.*wD0;
                    EM_S[INDEX2(7,6,8)]+=4.*wD0;
                    EM_S[INDEX2(0,7,8)]+=   wD0;
                    EM_S[INDEX2(1,7,8)]+=2.*wD0;
                    EM_S[INDEX2(2,7,8)]+=2.*wD0;
                    EM_S[INDEX2(3,7,8)]+=4.*wD0;
                    EM_S[INDEX2(4,7,8)]+=2.*wD0;
                    EM_S[INDEX2(5,7,8)]+=4.*wD0;
                    EM_S[INDEX2(6,7,8)]+=4.*wD0;
                    EM_S[INDEX2(7,7,8)]+=8.*wD0;
                }
            }
            ///////////////
            // process X //
            ///////////////
            if (!X.isEmpty()) {
                const Scalar* X_p = X.getSampleDataRO(id, zero);
                if (X.actsExpanded()) {
                    const Scalar X_0_0 = X_p[INDEX2(0,0,3)];
                    const Scalar X_1_0 = X_p[INDEX2(1,0,3)];
                    const Scalar X_2_0 = X_p[INDEX2(2,0,3)];
                    const Scalar X_0_1 = X_p[INDEX2(0,1,3)];
                    const Scalar X_1_1 = X_p[INDEX2(1,1,3)];
                    const Scalar X_2_1 = X_p[INDEX2(2,1,3)];
                    const Scalar X_0_2 = X_p[INDEX2(0,2,3)];
                    const Scalar X_1_2 = X_p[INDEX2(1,2,3)];
                    const Scalar X_2_2 = X_p[INDEX2(2,2,3)];
                    const Scalar X_0_3 = X_p[INDEX2(0,3,3)];
                    const Scalar X_1_3 = X_p[INDEX2(1,3,3)];
                    const Scalar X_2_3 = X_p[INDEX2(2,3,3)];
                    const Scalar X_0_4 = X_p[INDEX2(0,4,3)];
                    const Scalar X_1_4 = X_p[INDEX2(1,4,3)];
                    const Scalar X_2_4 = X_p[INDEX2(2,4,3)];
                    const Scalar X_0_5 = X_p[INDEX2(0,5,3)];
                    const Scalar X_1_5 = X_p[INDEX2(1,5,3)];
                    const Scalar X_2_5 = X_p[INDEX2(2,5,3)];
                    const Scalar X_0_6 = X_p[INDEX2(0,6,3)];
                    const Scalar X_1_6 = X_p[INDEX2(1,6,3)];
                    const Scalar X_2_6 = X_p[INDEX2(2,6,3)];
                    const Scalar X_0_7 = X_p[INDEX2(0,7,3)];
                    const Scalar X_1_7 = X_p[INDEX2(1,7,3)];
                    const Scalar X_2_7 = X_p[INDEX2(2,7,3)];
                    const Scalar tmp0  = w[66][l]*(X_0_2 + X_0_3 + X_0_4 + X_0_5);
                    const Scalar tmp1  = w[64][l]*(X_1_1 + X_1_3 + X_1_4 + X_1_6);
                    const Scalar tmp2  = w[61][l]*(X_0_0 + X_0_1);
                    const Scalar tmp3  = w[68][l]*(X_1_5 + X_1_7);
                    const Scalar tmp4  = w[65][l]*(X_2_1 + X_2_2 + X_2_5 + X_2_6);
                    const Scalar tmp5  = w[63][l]*(X_2_0 + X_2_4);
                    const Scalar tmp6  = w[67][l]*(X_2_3 + X_2_7);
                    const Scalar tmp7  = w[69][l]*(X_0_6 + X_0_7);
                    const Scalar tmp8  = w[62][l]*(X_1_0 + X_1_2);
                    const Scalar tmp9  = w[66][l]*(-X_0_2 - X_0_3 - X_0_4 - X_0_5);
                    const Scalar tmp10 = w[64][l]*(X_1_0 + X_1_2 + X_1_5 + X_1_7);
                    const Scalar tmp11 = w[61][l]*(-X_0_0 - X_0_1);
                    const Scalar tmp12 = w[68][l]*(X_1_4 + X_1_6);
                    const Scalar tmp13 = w[65][l]*(X_2_0 + X_2_3 + X_2_4 + X_2_7);
                    const Scalar tmp14 = w[63][l]*(X_2_1 + X_2_5);
                    const Scalar tmp15 = w[67][l]*(X_2_2 + X_2_6);
                    const Scalar tmp16 = w[69][l]*(-X_0_6 - X_0_7);
                    const Scalar tmp17 = w[62][l]*(X_1_1 + X_1_3);
                    const Scalar tmp18 = w[66][l]*(X_0_0 + X_0_1 + X_0_6 + X_0_7);
                    const Scalar tmp19 = w[64][l]*(-X_1_1 - X_1_3 - X_1_4 - X_1_6);
                    const Scalar tmp20 = w[61][l]*(X_0_2 + X_0_3);
                    const Scalar tmp21 = w[68][l]*(-X_1_5 - X_1_7);
                    const Scalar tmp22 = w[63][l]*(X_2_2 + X_2_6);
                    const Scalar tmp23 = w[67][l]*(X_2_1 + X_2_5);
                    const Scalar tmp24 = w[69][l]*(X_0_4 + X_0_5);
                    const Scalar tmp25 = w[62][l]*(-X_1_0 - X_1_2);
                    const Scalar tmp26 = w[66][l]*(-X_0_0 - X_0_1 - X_0_6 - X_0_7);
                    const Scalar tmp27 = w[64][l]*(-X_1_0 - X_1_2 - X_1_5 - X_1_7);
                    const Scalar tmp28 = w[61][l]*(-X_0_2 - X_0_3);
                    const Scalar tmp29 = w[68][l]*(-X_1_4 - X_1_6);
                    const Scalar tmp30 = w[63][l]*(X_2_3 + X_2_7);
                    const Scalar tmp31 = w[67][l]*(X_2_0 + X_2_4);
                    const Scalar tmp32 = w[69][l]*(-X_0_4 - X_0_5);
                    const Scalar tmp33 = w[62][l]*(-X_1_1 - X_1_3);
                    const Scalar tmp34 = w[61][l]*(X_0_4 + X_0_5);
                    const Scalar tmp35 = w[68][l]*(X_1_1 + X_1_3);
                    const Scalar tmp36 = w[65][l]*(-X_2_1 - X_2_2 - X_2_5 - X_2_6);
                    const Scalar tmp37 = w[63][l]*(-X_2_0 - X_2_4);
                    const Scalar tmp38 = w[67][l]*(-X_2_3 - X_2_7);
                    const Scalar tmp39 = w[69][l]*(X_0_2 + X_0_3);
                    const Scalar tmp40 = w[62][l]*(X_1_4 + X_1_6);
                    const Scalar tmp41 = w[61][l]*(-X_0_4 - X_0_5);
                    const Scalar tmp42 = w[68][l]*(X_1_0 + X_1_2);
                    const Scalar tmp43 = w[65][l]*(-X_2_0 - X_2_3 - X_2_4 - X_2_7);
                    const Scalar tmp44 = w[63][l]*(-X_2_1 - X_2_5);
                    const Scalar tmp45 = w[67][l]*(-X_2_2 - X_2_6);
                    const Scalar tmp46 = w[69][l]*(-X_0_2 - X_0_3);
                    const Scalar tmp47 = w[62][l]*(X_1_5 + X_1_7);
                    const Scalar tmp48 = w[61][l]*(X_0_6 + X_0_7);
                    const Scalar tmp49 = w[68][l]*(-X_1_1 - X_1_3);
                    const Scalar tmp50 = w[63][l]*(-X_2_2 - X_2_6);
                    const Scalar tmp51 = w[67][l]*(-X_2_1 - X_2_5);
                    const Scalar tmp52 = w[69][l]*(X_0_0 + X_0_1);
                    const Scalar tmp53 = w[62][l]*(-X_1_4 - X_1_6);
                    const Scalar tmp54 = w[61][l]*(-X_0_6 - X_0_7);
                    const Scalar tmp55 = w[68][l]*(-X_1_0 - X_1_2);
                    const Scalar tmp56 = w[63][l]*(-X_2_3 - X_2_7);
                    const Scalar tmp57 = w[67][l]*(-X_2_0 - X_2_4);
                    const Scalar tmp58 = w[69][l]*(-X_0_0 - X_0_1);
                    const Scalar tmp59 = w[62][l]*(-X_1_5 - X_1_7);
                    EM_F[0]+=tmp0 + tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6 + tmp7 + tmp8;
                    EM_F[1]+=tmp10 + tmp11 + tmp12 + tmp13 + tmp14 + tmp15 + tmp16 + tmp17 + tmp9;
                    EM_F[2]+=tmp13 + tmp18 + tmp19 + tmp20 + tmp21 + tmp22 + tmp23 + tmp24 + tmp25;
                    EM_F[3]+=tmp26 + tmp27 + tmp28 + tmp29 + tmp30 + tmp31 + tmp32 + tmp33 + tmp4;
                    EM_F[4]+=tmp10 + tmp18 + tmp34 + tmp35 + tmp36 + tmp37 + tmp38 + tmp39 + tmp40;
                    EM_F[5]+=tmp1 + tmp26 + tmp41 + tmp42 + tmp43 + tmp44 + tmp45 + tmp46 + tmp47;
                    EM_F[6]+=tmp0 + tmp27 + tmp43 + tmp48 + tmp49 + tmp50 + tmp51 + tmp52 + tmp53;
                    EM_F[7]+=tmp19 + tmp36 + tmp54 + tmp55 + tmp56 + tmp57 + tmp58 + tmp59 + tmp9;
                } else { // constant data
                    const Scalar wX0 = 12.*X_p[0]*w[66][l];
                    const Scalar wX1 = 12.*X_p[1]*w[64][l];
                    const Scalar wX2 = 18.*X_p[2]*w[50][l];
                    EM_F[0]+= wX0 + wX1 - wX2;
                    EM_F[1]+=-wX0 + wX1 - wX2;
                    EM_F[2]+= wX0 - wX1 - wX2;
                    EM_F[3]+=-wX0 - wX1 - wX2;
                    EM_F[4]+= wX0 + wX1 + wX2;
                    EM_F[5]+=-wX0 + wX1 + wX2;
                    EM_F[6]+= wX0 - wX1 + wX2;
                    EM_F[7]+=-wX0 - wX1 + wX2;
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
                    const Scalar Y_4 = Y_p[4];
                    const Scalar Y_5 = Y_p[5];
                    const Scalar Y_6 = Y_p[6];
                    const Scalar Y_7 = Y_p[7];
                    const Scalar tmp0  = w[72][l]*(Y_3 + Y_5 + Y_6);
                    const Scalar tmp1  = w[71][l]*(Y_1 + Y_2 + Y_4);
                    const Scalar tmp2  = w[72][l]*(Y_2 + Y_4 + Y_7);
                    const Scalar tmp3  = w[71][l]*(Y_0 + Y_3 + Y_5);
                    const Scalar tmp4  = w[72][l]*(Y_1 + Y_4 + Y_7);
                    const Scalar tmp5  = w[71][l]*(Y_0 + Y_3 + Y_6);
                    const Scalar tmp6  = w[72][l]*(Y_0 + Y_5 + Y_6);
                    const Scalar tmp7  = w[71][l]*(Y_1 + Y_2 + Y_7);
                    const Scalar tmp8  = w[72][l]*(Y_1 + Y_2 + Y_7);
                    const Scalar tmp9  = w[71][l]*(Y_0 + Y_5 + Y_6);
                    const Scalar tmp10 = w[72][l]*(Y_0 + Y_3 + Y_6);
                    const Scalar tmp11 = w[71][l]*(Y_1 + Y_4 + Y_7);
                    const Scalar tmp12 = w[72][l]*(Y_0 + Y_3 + Y_5);
                    const Scalar tmp13 = w[71][l]*(Y_2 + Y_4 + Y_7);
                    const Scalar tmp14 = w[72][l]*(Y_1 + Y_2 + Y_4);
                    const Scalar tmp15 = w[71][l]*(Y_3 + Y_5 + Y_6);
                    EM_F[0]+=Y_0*w[70][l] + Y_7*w[73][l] + tmp0 + tmp1;
                    EM_F[1]+=Y_1*w[70][l] + Y_6*w[73][l] + tmp2 + tmp3;
                    EM_F[2]+=Y_2*w[70][l] + Y_5*w[73][l] + tmp4 + tmp5;
                    EM_F[3]+=Y_3*w[70][l] + Y_4*w[73][l] + tmp6 + tmp7;
                    EM_F[4]+=Y_3*w[73][l] + Y_4*w[70][l] + tmp8 + tmp9;
                    EM_F[5]+=Y_2*w[73][l] + Y_5*w[70][l] + tmp10 + tmp11;
                    EM_F[6]+=Y_1*w[73][l] + Y_6*w[70][l] + tmp12 + tmp13;
                    EM_F[7]+=Y_0*w[73][l] + Y_7*w[70][l] + tmp14 + tmp15;
                } else { // constant data
                    EM_F[0]+=216.*Y_p[0]*w[55][l];
                    EM_F[1]+=216.*Y_p[0]*w[55][l];
                    EM_F[2]+=216.*Y_p[0]*w[55][l];
                    EM_F[3]+=216.*Y_p[0]*w[55][l];
                    EM_F[4]+=216.*Y_p[0]*w[55][l];
                    EM_F[5]+=216.*Y_p[0]*w[55][l];
                    EM_F[6]+=216.*Y_p[0]*w[55][l];
                    EM_F[7]+=216.*Y_p[0]*w[55][l];
                }
            }

            // add to matrix (if add_EM_S) and RHS (if add_EM_F)
            // const index_t firstNode=m_NN[0]*m_NN[1]*k2+m_NN[0]*k1+k0;
            domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, add_EM_S, add_EM_F, q, t);
        } // end of colouring
    } // end of parallel region
}

/****************************************************************************/
// PDE SINGLE BOUNDARY
/****************************************************************************/

template<class Scalar>
void DefaultAssembler3D<Scalar>::assemblePDEBoundarySingle(
                                        AbstractSystemMatrix* mat, Data& rhs,
                                        const Data& d, const Data& y) const
{

    // Find the maximum level of refinement in the mesh
    int max_level = 0;
    for(p8est_topidx_t tree = domain->p8est->first_local_tree; tree <= domain->p8est->last_local_tree; tree++) 
    {
        p8est_tree_t * tree_t = p8est_tree_array_index(domain->p8est->trees, tree);
        max_level = tree_t->maxlevel > max_level ? tree_t->maxlevel : max_level;
    }

    const double SQRT3 = 1.73205080756887719318;
    double w[16][P4EST_MAXLEVEL] = {{0}};
#pragma omp parallel for
    for(int i = 0; i < max_level; i++)
    {
        double m_dx[3] = {domain->m_NX[0]*domain->forestData->m_dx[0][P4EST_MAXLEVEL-i], 
                          domain->m_NX[1]*domain->forestData->m_dx[1][P4EST_MAXLEVEL-i],
                          domain->m_NX[2]*domain->forestData->m_dx[2][P4EST_MAXLEVEL-i]};
        w[12][i] = m_dx[0]*m_dx[1]/144;
        w[10][i] = w[12][i]*(-SQRT3 + 2);
        w[11][i] = w[12][i]*(SQRT3 + 2);
        w[13][i] = w[12][i]*(-4*SQRT3 + 7);
        w[14][i] = w[12][i]*(4*SQRT3 + 7);
        w[7][i] = m_dx[0]*m_dx[2]/144;
        w[5][i] = w[7][i]*(-SQRT3 + 2);
        w[6][i] = w[7][i]*(SQRT3 + 2);
        w[8][i] = w[7][i]*(-4*SQRT3 + 7);
        w[9][i] = w[7][i]*(4*SQRT3 + 7);
        w[2][i] = m_dx[1]*m_dx[2]/144;
        w[0][i] = w[2][i]*(-SQRT3 + 2);
        w[1][i] = w[2][i]*(SQRT3 + 2);
        w[3][i] = w[2][i]*(-4*SQRT3 + 7);
        w[4][i] = w[2][i]*(4*SQRT3 + 7);
    }



//     const dim_t NE0 = m_NE[0];
//     const dim_t NE1 = m_NE[1];
//     const dim_t NE2 = m_NE[2];

    const bool add_EM_S = !d.isEmpty();
    const bool add_EM_F = !y.isEmpty();
    const Scalar zero = static_cast<Scalar>(0);
    rhs.requireWrite();

    vector<Scalar> EM_S(8*8);
    vector<Scalar> EM_F(8);

    for(p8est_topidx_t t = domain->p8est->first_local_tree; t <= domain->p8est->last_local_tree; t++) 
    {
        p8est_tree_t * currenttree = p8est_tree_array_index(domain->p8est->trees, t);
        sc_array_t * tquadrants = &currenttree->quadrants;
        p4est_qcoord_t Q = (p4est_qcoord_t) tquadrants->elem_count;

#pragma omp parallel for
        for (int q = 0; q < Q; ++q)  // Loop over the elements attached to the tree
        {
            if (add_EM_S)
                fill(EM_S.begin(), EM_S.end(), zero);
            if (add_EM_F)
                fill(EM_F.begin(), EM_F.end(), zero);
            
            p8est_quadrant_t * quad = p8est_quadrant_array_index(tquadrants, q);
            octantData * quaddata = (octantData *) quad->p.user_data;  

            int l = quad->level;
            double xyz[3];
            p8est_qcoord_to_vertex(domain->p8est->connectivity, t, quad->x, quad->y, quad->z, xyz);
            long e = domain->NodeIDs.find(std::make_tuple(xyz[0],xyz[1],xyz[2]))->second;

            if(quaddata->m_faceOffset[0]) {
                if (add_EM_S)
                    fill(EM_S.begin(), EM_S.end(), zero);
                if (add_EM_F) {
                    EM_F[1] = zero;
                    EM_F[3] = zero;
                    EM_F[5] = zero;
                    EM_F[7] = zero;
                }
                
                ///////////////
                // process d //
                ///////////////
                if (add_EM_S) {
                    const Scalar* d_p = d.getSampleDataRO(e, zero);
                    if (d.actsExpanded()) {
                        const Scalar d_0 = d_p[0];
                        const Scalar d_1 = d_p[1];
                        const Scalar d_2 = d_p[2];
                        const Scalar d_3 = d_p[3];
                        const Scalar tmp0 = w[0][l]*(d_0 + d_1);
                        const Scalar tmp1 = w[1][l]*(d_2 + d_3);
                        const Scalar tmp2 = w[0][l]*(d_0 + d_2);
                        const Scalar tmp3 = w[1][l]*(d_1 + d_3);
                        const Scalar tmp4 = w[0][l]*(d_1 + d_3);
                        const Scalar tmp5 = w[1][l]*(d_0 + d_2);
                        const Scalar tmp6 = w[0][l]*(d_2 + d_3);
                        const Scalar tmp7 = w[1][l]*(d_0 + d_1);
                        const Scalar tmp8 = w[2][l]*(d_0 + d_3);
                        const Scalar tmp9 = w[2][l]*(d_1 + d_2);
                        const Scalar tmp10 = w[2][l]*(d_0 + d_1 + d_2 + d_3);
                        EM_S[INDEX2(0,0,8)] = d_0*w[4][l] + d_3*w[3][l] + tmp9;
                        EM_S[INDEX2(2,0,8)] = tmp6 + tmp7;
                        EM_S[INDEX2(4,0,8)] = tmp4 + tmp5;
                        EM_S[INDEX2(6,0,8)] = tmp10;
                        EM_S[INDEX2(0,2,8)] = tmp6 + tmp7;
                        EM_S[INDEX2(2,2,8)] = d_1*w[4][l] + d_2*w[3][l] + tmp8;
                        EM_S[INDEX2(4,2,8)] = tmp10;
                        EM_S[INDEX2(6,2,8)] = tmp2 + tmp3;
                        EM_S[INDEX2(0,4,8)] = tmp4 + tmp5;
                        EM_S[INDEX2(2,4,8)] = tmp10;
                        EM_S[INDEX2(4,4,8)] = d_1*w[3][l] + d_2*w[4][l] + tmp8;
                        EM_S[INDEX2(6,4,8)] = tmp0 + tmp1;
                        EM_S[INDEX2(0,6,8)] = tmp10;
                        EM_S[INDEX2(2,6,8)] = tmp2 + tmp3;
                        EM_S[INDEX2(4,6,8)] = tmp0 + tmp1;
                        EM_S[INDEX2(6,6,8)] = d_0*w[3][l] + d_3*w[4][l] + tmp9;
                    } else { // constant data
                        const Scalar wd0 = 4.*d_p[0]*w[2][l];
                        EM_S[INDEX2(0,0,8)] = 4.*wd0;
                        EM_S[INDEX2(2,0,8)] = 2.*wd0;
                        EM_S[INDEX2(4,0,8)] = 2.*wd0;
                        EM_S[INDEX2(6,0,8)] =    wd0;
                        EM_S[INDEX2(0,2,8)] = 2.*wd0;
                        EM_S[INDEX2(2,2,8)] = 4.*wd0;
                        EM_S[INDEX2(4,2,8)] =    wd0;
                        EM_S[INDEX2(6,2,8)] = 2.*wd0;
                        EM_S[INDEX2(0,4,8)] = 2.*wd0;
                        EM_S[INDEX2(2,4,8)] =    wd0;
                        EM_S[INDEX2(4,4,8)] = 4.*wd0;
                        EM_S[INDEX2(6,4,8)] = 2.*wd0;
                        EM_S[INDEX2(0,6,8)] =    wd0;
                        EM_S[INDEX2(2,6,8)] = 2.*wd0;
                        EM_S[INDEX2(4,6,8)] = 2.*wd0;
                        EM_S[INDEX2(6,6,8)] = 4.*wd0;
                    }
                }
                ///////////////
                // process y //
                ///////////////
                if (add_EM_F) {
                    const Scalar* y_p = y.getSampleDataRO(e, zero);
                    if (y.actsExpanded()) {
                        const Scalar y_0 = y_p[0];
                        const Scalar y_1 = y_p[1];
                        const Scalar y_2 = y_p[2];
                        const Scalar y_3 = y_p[3];
                        const Scalar tmp0 = 6.*w[2][l]*(y_1 + y_2);
                        const Scalar tmp1 = 6.*w[2][l]*(y_0 + y_3);
                        EM_F[0] = tmp0 + 6.*w[0][l]*y_3 + 6.*w[1][l]*y_0;
                        EM_F[2] = tmp1 + 6.*w[0][l]*y_2 + 6.*w[1][l]*y_1;
                        EM_F[4] = tmp1 + 6.*w[0][l]*y_1 + 6.*w[1][l]*y_2;
                        EM_F[6] = tmp0 + 6.*w[0][l]*y_0 + 6.*w[1][l]*y_3;
                    } else { // constant data
                        EM_F[0] = 36.*w[2][l]*y_p[0];
                        EM_F[2] = 36.*w[2][l]*y_p[0];
                        EM_F[4] = 36.*w[2][l]*y_p[0];
                        EM_F[6] = 36.*w[2][l]*y_p[0];
                    }
                }
                
                domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, add_EM_S, add_EM_F, q, t);
            } // face 0

            if(quaddata->m_faceOffset[1]) {
                if (add_EM_S)
                    fill(EM_S.begin(), EM_S.end(), zero);
                if (add_EM_F) {
                    EM_F[0] = zero;
                    EM_F[2] = zero;
                    EM_F[4] = zero;
                    EM_F[6] = zero;
                }
                ///////////////
                // process d //
                ///////////////
                if (add_EM_S) {
                    const Scalar* d_p = d.getSampleDataRO(e, zero);
                    if (d.actsExpanded()) {
                        const Scalar d_0 = d_p[0];
                        const Scalar d_1 = d_p[1];
                        const Scalar d_2 = d_p[2];
                        const Scalar d_3 = d_p[3];
                        const Scalar tmp0  = w[0][l]*(d_0 + d_2);
                        const Scalar tmp1  = w[1][l]*(d_1 + d_3);
                        const Scalar tmp2  = w[0][l]*(d_2 + d_3);
                        const Scalar tmp3  = w[1][l]*(d_0 + d_1);
                        const Scalar tmp4  = w[0][l]*(d_1 + d_3);
                        const Scalar tmp5  = w[1][l]*(d_0 + d_2);
                        const Scalar tmp6  = w[2][l]*(d_0 + d_3);
                        const Scalar tmp7  = w[2][l]*(d_1 + d_2);
                        const Scalar tmp8  = w[0][l]*(d_0 + d_1);
                        const Scalar tmp9  = w[1][l]*(d_2 + d_3);
                        const Scalar tmp10 = w[2][l]*(d_0 + d_1 + d_2 + d_3);
                        EM_S[INDEX2(1,1,8)] = d_0*w[4][l] + d_3*w[3][l] + tmp7;
                        EM_S[INDEX2(3,1,8)] = tmp2 + tmp3;
                        EM_S[INDEX2(5,1,8)] = tmp4 + tmp5;
                        EM_S[INDEX2(7,1,8)] = tmp10;
                        EM_S[INDEX2(1,3,8)] = tmp2 + tmp3;
                        EM_S[INDEX2(3,3,8)] = d_1*w[4][l] + d_2*w[3][l] + tmp6;
                        EM_S[INDEX2(5,3,8)] = tmp10;
                        EM_S[INDEX2(7,3,8)] = tmp0 + tmp1;
                        EM_S[INDEX2(1,5,8)] = tmp4 + tmp5;
                        EM_S[INDEX2(3,5,8)] = tmp10;
                        EM_S[INDEX2(5,5,8)] = d_1*w[3][l] + d_2*w[4][l] + tmp6;
                        EM_S[INDEX2(7,5,8)] = tmp8 + tmp9;
                        EM_S[INDEX2(1,7,8)] = tmp10;
                        EM_S[INDEX2(3,7,8)] = tmp0 + tmp1;
                        EM_S[INDEX2(5,7,8)] = tmp8 + tmp9;
                        EM_S[INDEX2(7,7,8)] = d_0*w[3][l] + d_3*w[4][l] + tmp7;
                    } else { // constant data
                        const Scalar wd0 = 4.*d_p[0]*w[2][l];
                        EM_S[INDEX2(1,1,8)] = 4.*wd0;
                        EM_S[INDEX2(3,1,8)] = 2.*wd0;
                        EM_S[INDEX2(5,1,8)] = 2.*wd0;
                        EM_S[INDEX2(7,1,8)] =    wd0;
                        EM_S[INDEX2(1,3,8)] = 2.*wd0;
                        EM_S[INDEX2(3,3,8)] = 4.*wd0;
                        EM_S[INDEX2(5,3,8)] =    wd0;
                        EM_S[INDEX2(7,3,8)] = 2.*wd0;
                        EM_S[INDEX2(1,5,8)] = 2.*wd0;
                        EM_S[INDEX2(3,5,8)] =    wd0;
                        EM_S[INDEX2(5,5,8)] = 4.*wd0;
                        EM_S[INDEX2(7,5,8)] = 2.*wd0;
                        EM_S[INDEX2(1,7,8)] =    wd0;
                        EM_S[INDEX2(3,7,8)] = 2.*wd0;
                        EM_S[INDEX2(5,7,8)] = 2.*wd0;
                        EM_S[INDEX2(7,7,8)] = 4.*wd0;
                    }
                }
                ///////////////
                // process y //
                ///////////////
                if (add_EM_F) {
                    const Scalar* y_p = y.getSampleDataRO(e, zero);
                    if (y.actsExpanded()) {
                        const Scalar y_0 = y_p[0];
                        const Scalar y_1 = y_p[1];
                        const Scalar y_2 = y_p[2];
                        const Scalar y_3 = y_p[3];
                        const Scalar tmp0 = 6.*w[2][l]*(y_1 + y_2);
                        const Scalar tmp1 = 6.*w[2][l]*(y_0 + y_3);
                        EM_F[1] = tmp0 + 6.*w[0][l]*y_3 + 6.*w[1][l]*y_0;
                        EM_F[3] = tmp1 + 6.*w[0][l]*y_2 + 6.*w[1][l]*y_1;
                        EM_F[5] = tmp1 + 6.*w[0][l]*y_1 + 6.*w[1][l]*y_2;
                        EM_F[7] = tmp0 + 6.*w[0][l]*y_0 + 6.*w[1][l]*y_3;
                    } else { // constant data
                        EM_F[1] = 36.*w[2][l]*y_p[0];
                        EM_F[3] = 36.*w[2][l]*y_p[0];
                        EM_F[5] = 36.*w[2][l]*y_p[0];
                        EM_F[7] = 36.*w[2][l]*y_p[0];
                    }
                }
                domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, add_EM_S, add_EM_F, q, t);
            } // face 1

            if(quaddata->m_faceOffset[2]) {
                if (add_EM_S)
                    fill(EM_S.begin(), EM_S.end(), zero);
                if (add_EM_F) {
                    EM_F[2] = zero;
                    EM_F[3] = zero;
                    EM_F[6] = zero;
                    EM_F[7] = zero;
                }
                ///////////////
                // process d //
                ///////////////
                if (add_EM_S) {
                    const Scalar* d_p = d.getSampleDataRO(e, zero);
                    if (d.actsExpanded()) {
                        const Scalar d_0 = d_p[0];
                        const Scalar d_1 = d_p[1];
                        const Scalar d_2 = d_p[2];
                        const Scalar d_3 = d_p[3];
                        const Scalar tmp0  = w[5][l]*(d_0 + d_1);
                        const Scalar tmp1  = w[6][l]*(d_2 + d_3);
                        const Scalar tmp2  = w[5][l]*(d_0 + d_2);
                        const Scalar tmp3  = w[6][l]*(d_1 + d_3);
                        const Scalar tmp4  = w[5][l]*(d_1 + d_3);
                        const Scalar tmp5  = w[6][l]*(d_0 + d_2);
                        const Scalar tmp6  = w[7][l]*(d_0 + d_3);
                        const Scalar tmp7  = w[7][l]*(d_0 + d_1 + d_2 + d_3);
                        const Scalar tmp8  = w[7][l]*(d_1 + d_2);
                        const Scalar tmp9  = w[5][l]*(d_2 + d_3);
                        const Scalar tmp10 = w[6][l]*(d_0 + d_1);
                        EM_S[INDEX2(0,0,8)] = d_0*w[9][l] + d_3*w[8][l] + tmp8;
                        EM_S[INDEX2(1,0,8)] = tmp10 + tmp9;
                        EM_S[INDEX2(4,0,8)] = tmp4 + tmp5;
                        EM_S[INDEX2(5,0,8)] = tmp7;
                        EM_S[INDEX2(0,1,8)] = tmp10 + tmp9;
                        EM_S[INDEX2(1,1,8)] = d_1*w[9][l] + d_2*w[8][l] + tmp6;
                        EM_S[INDEX2(4,1,8)] = tmp7;
                        EM_S[INDEX2(5,1,8)] = tmp2 + tmp3;
                        EM_S[INDEX2(0,4,8)] = tmp4 + tmp5;
                        EM_S[INDEX2(1,4,8)] = tmp7;
                        EM_S[INDEX2(4,4,8)] = d_1*w[8][l] + d_2*w[9][l] + tmp6;
                        EM_S[INDEX2(5,4,8)] = tmp0 + tmp1;
                        EM_S[INDEX2(0,5,8)] = tmp7;
                        EM_S[INDEX2(1,5,8)] = tmp2 + tmp3;
                        EM_S[INDEX2(4,5,8)] = tmp0 + tmp1;
                        EM_S[INDEX2(5,5,8)] = d_0*w[8][l] + d_3*w[9][l] + tmp8;
                    } else { // constant data
                        const Scalar wd0 = 4.*d_p[0]*w[7][l];
                        EM_S[INDEX2(0,0,8)] = 4.*wd0;
                        EM_S[INDEX2(1,0,8)] = 2.*wd0;
                        EM_S[INDEX2(4,0,8)] = 2.*wd0;
                        EM_S[INDEX2(5,0,8)] =    wd0;
                        EM_S[INDEX2(0,1,8)] = 2.*wd0;
                        EM_S[INDEX2(1,1,8)] = 4.*wd0;
                        EM_S[INDEX2(4,1,8)] =    wd0;
                        EM_S[INDEX2(5,1,8)] = 2.*wd0;
                        EM_S[INDEX2(0,4,8)] = 2.*wd0;
                        EM_S[INDEX2(1,4,8)] =    wd0;
                        EM_S[INDEX2(4,4,8)] = 4.*wd0;
                        EM_S[INDEX2(5,4,8)] = 2.*wd0;
                        EM_S[INDEX2(0,5,8)] =    wd0;
                        EM_S[INDEX2(1,5,8)] = 2.*wd0;
                        EM_S[INDEX2(4,5,8)] = 2.*wd0;
                        EM_S[INDEX2(5,5,8)] = 4.*wd0;
                    }
                }
                ///////////////
                // process y //
                ///////////////
                if (add_EM_F) {
                    const Scalar* y_p = y.getSampleDataRO(e, zero);
                    if (y.actsExpanded()) {
                        const Scalar y_0 = y_p[0];
                        const Scalar y_1 = y_p[1];
                        const Scalar y_2 = y_p[2];
                        const Scalar y_3 = y_p[3];
                        const Scalar tmp0 = 6.*w[7][l]*(y_1 + y_2);
                        const Scalar tmp1 = 6.*w[7][l]*(y_0 + y_3);
                        EM_F[0] = tmp0 + 6.*w[5][l]*y_3 + 6.*w[6][l]*y_0;
                        EM_F[1] = tmp1 + 6.*w[5][l]*y_2 + 6.*w[6][l]*y_1;
                        EM_F[4] = tmp1 + 6.*w[5][l]*y_1 + 6.*w[6][l]*y_2;
                        EM_F[5] = tmp0 + 6.*w[5][l]*y_0 + 6.*w[6][l]*y_3;
                    } else { // constant data
                        EM_F[0] = 36.*w[7][l]*y_p[0];
                        EM_F[1] = 36.*w[7][l]*y_p[0];
                        EM_F[4] = 36.*w[7][l]*y_p[0];
                        EM_F[5] = 36.*w[7][l]*y_p[0];
                    }
                }
                domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, add_EM_S, add_EM_F, q, t);
            } // face 2

            if(quaddata->m_faceOffset[3]) {
                if (add_EM_S)
                    fill(EM_S.begin(), EM_S.end(), zero);
                if (add_EM_F) {
                    EM_F[0] = zero;
                    EM_F[1] = zero;
                    EM_F[4] = zero;
                    EM_F[5] = zero;
                }
                ///////////////
                // process d //
                ///////////////
                if (add_EM_S) {
                    const Scalar* d_p = d.getSampleDataRO(e, zero);
                    if (d.actsExpanded()) {
                        const Scalar d_0 = d_p[0];
                        const Scalar d_1 = d_p[1];
                        const Scalar d_2 = d_p[2];
                        const Scalar d_3 = d_p[3];
                        const Scalar tmp0  = w[5][l]*(d_0 + d_2);
                        const Scalar tmp1  = w[6][l]*(d_1 + d_3);
                        const Scalar tmp2  = w[5][l]*(d_1 + d_3);
                        const Scalar tmp3  = w[6][l]*(d_0 + d_2);
                        const Scalar tmp4  = w[7][l]*(d_0 + d_1 + d_2 + d_3);
                        const Scalar tmp5  = w[5][l]*(d_0 + d_1);
                        const Scalar tmp6  = w[6][l]*(d_2 + d_3);
                        const Scalar tmp7  = w[7][l]*(d_0 + d_3);
                        const Scalar tmp8  = w[7][l]*(d_1 + d_2);
                        const Scalar tmp9  = w[5][l]*(d_2 + d_3);
                        const Scalar tmp10 = w[6][l]*(d_0 + d_1);
                        EM_S[INDEX2(2,2,8)] = d_0*w[9][l] + d_3*w[8][l] + tmp8;
                        EM_S[INDEX2(3,2,8)] = tmp10 + tmp9;
                        EM_S[INDEX2(6,2,8)] = tmp2 + tmp3;
                        EM_S[INDEX2(7,2,8)] = tmp4;
                        EM_S[INDEX2(2,3,8)] = tmp10 + tmp9;
                        EM_S[INDEX2(3,3,8)] = d_1*w[9][l] + d_2*w[8][l] + tmp7;
                        EM_S[INDEX2(6,3,8)] = tmp4;
                        EM_S[INDEX2(7,3,8)] = tmp0 + tmp1;
                        EM_S[INDEX2(2,6,8)] = tmp2 + tmp3;
                        EM_S[INDEX2(3,6,8)] = tmp4;
                        EM_S[INDEX2(6,6,8)] = d_1*w[8][l] + d_2*w[9][l] + tmp7;
                        EM_S[INDEX2(7,6,8)] = tmp5 + tmp6;
                        EM_S[INDEX2(2,7,8)] = tmp4;
                        EM_S[INDEX2(3,7,8)] = tmp0 + tmp1;
                        EM_S[INDEX2(6,7,8)] = tmp5 + tmp6;
                        EM_S[INDEX2(7,7,8)] = d_0*w[8][l] + d_3*w[9][l] + tmp8;
                    } else { // constant data
                        const Scalar wd0 = 4.*d_p[0]*w[7][l];
                        EM_S[INDEX2(2,2,8)] = 4.*wd0;
                        EM_S[INDEX2(3,2,8)] = 2.*wd0;
                        EM_S[INDEX2(6,2,8)] = 2.*wd0;
                        EM_S[INDEX2(7,2,8)] =    wd0;
                        EM_S[INDEX2(2,3,8)] = 2.*wd0;
                        EM_S[INDEX2(3,3,8)] = 4.*wd0;
                        EM_S[INDEX2(6,3,8)] =    wd0;
                        EM_S[INDEX2(7,3,8)] = 2.*wd0;
                        EM_S[INDEX2(2,6,8)] = 2.*wd0;
                        EM_S[INDEX2(3,6,8)] =    wd0;
                        EM_S[INDEX2(6,6,8)] = 4.*wd0;
                        EM_S[INDEX2(7,6,8)] = 2.*wd0;
                        EM_S[INDEX2(2,7,8)] =    wd0;
                        EM_S[INDEX2(3,7,8)] = 2.*wd0;
                        EM_S[INDEX2(6,7,8)] = 2.*wd0;
                        EM_S[INDEX2(7,7,8)] = 4.*wd0;
                    }
                }
                ///////////////
                // process y //
                ///////////////
                if (add_EM_F) {
                    const Scalar* y_p = y.getSampleDataRO(e, zero);
                    if (y.actsExpanded()) {
                        const Scalar y_0 = y_p[0];
                        const Scalar y_1 = y_p[1];
                        const Scalar y_2 = y_p[2];
                        const Scalar y_3 = y_p[3];
                        const Scalar tmp0 = 6.*w[7][l]*(y_1 + y_2);
                        const Scalar tmp1 = 6.*w[7][l]*(y_0 + y_3);
                        EM_F[2] = tmp0 + 6.*w[5][l]*y_3 + 6.*w[6][l]*y_0;
                        EM_F[3] = tmp1 + 6.*w[5][l]*y_2 + 6.*w[6][l]*y_1;
                        EM_F[6] = tmp1 + 6.*w[5][l]*y_1 + 6.*w[6][l]*y_2;
                        EM_F[7] = tmp0 + 6.*w[5][l]*y_0 + 6.*w[6][l]*y_3;
                    } else { // constant data
                        EM_F[2] = 36.*w[7][l]*y_p[0];
                        EM_F[3] = 36.*w[7][l]*y_p[0];
                        EM_F[6] = 36.*w[7][l]*y_p[0];
                        EM_F[7] = 36.*w[7][l]*y_p[0];
                    }
                }
                domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, add_EM_S, add_EM_F, q, t);
            } // face 3

            if(quaddata->m_faceOffset[4]) {
                if (add_EM_S)
                    fill(EM_S.begin(), EM_S.end(), zero);
                if (add_EM_F) {
                    EM_F[4] = zero;
                    EM_F[5] = zero;
                    EM_F[6] = zero;
                    EM_F[7] = zero;
                }
                ///////////////
                // process d //
                ///////////////
                if (add_EM_S) {
                    const Scalar* d_p = d.getSampleDataRO(e, zero);
                    if (d.actsExpanded()) {
                        const Scalar d_0 = d_p[0];
                        const Scalar d_1 = d_p[1];
                        const Scalar d_2 = d_p[2];
                        const Scalar d_3 = d_p[3];
                        const Scalar tmp0  = w[10][l]*(d_0 + d_2);
                        const Scalar tmp1  = w[11][l]*(d_1 + d_3);
                        const Scalar tmp2  = w[12][l]*(d_0 + d_1 + d_2 + d_3);
                        const Scalar tmp3  = w[12][l]*(d_1 + d_2);
                        const Scalar tmp4  = w[10][l]*(d_1 + d_3);
                        const Scalar tmp5  = w[11][l]*(d_0 + d_2);
                        const Scalar tmp6  = w[12][l]*(d_0 + d_3);
                        const Scalar tmp7  = w[10][l]*(d_0 + d_1);
                        const Scalar tmp8  = w[11][l]*(d_2 + d_3);
                        const Scalar tmp9  = w[10][l]*(d_2 + d_3);
                        const Scalar tmp10 = w[11][l]*(d_0 + d_1);
                        EM_S[INDEX2(0,0,8)] = d_0*w[14][l] + d_3*w[13][l] + tmp3;
                        EM_S[INDEX2(1,0,8)] = tmp10 + tmp9;
                        EM_S[INDEX2(2,0,8)] = tmp4 + tmp5;
                        EM_S[INDEX2(3,0,8)] = tmp2;
                        EM_S[INDEX2(0,1,8)] = tmp10 + tmp9;
                        EM_S[INDEX2(1,1,8)] = d_1*w[14][l] + d_2*w[13][l] + tmp6;
                        EM_S[INDEX2(2,1,8)] = tmp2;
                        EM_S[INDEX2(3,1,8)] = tmp0 + tmp1;
                        EM_S[INDEX2(0,2,8)] = tmp4 + tmp5;
                        EM_S[INDEX2(1,2,8)] = tmp2;
                        EM_S[INDEX2(2,2,8)] = d_1*w[13][l] + d_2*w[14][l] + tmp6;
                        EM_S[INDEX2(3,2,8)] = tmp7 + tmp8;
                        EM_S[INDEX2(0,3,8)] = tmp2;
                        EM_S[INDEX2(1,3,8)] = tmp0 + tmp1;
                        EM_S[INDEX2(2,3,8)] = tmp7 + tmp8;
                        EM_S[INDEX2(3,3,8)] = d_0*w[13][l] + d_3*w[14][l] + tmp3;
                    } else { // constant data
                        const Scalar wd0 = 4.*d_p[0]*w[12][l];
                        EM_S[INDEX2(0,0,8)] = 4.*wd0;
                        EM_S[INDEX2(1,0,8)] = 2.*wd0;
                        EM_S[INDEX2(2,0,8)] = 2.*wd0;
                        EM_S[INDEX2(3,0,8)] =    wd0;
                        EM_S[INDEX2(0,1,8)] = 2.*wd0;
                        EM_S[INDEX2(1,1,8)] = 4.*wd0;
                        EM_S[INDEX2(2,1,8)] =    wd0;
                        EM_S[INDEX2(3,1,8)] = 2.*wd0;
                        EM_S[INDEX2(0,2,8)] = 2.*wd0;
                        EM_S[INDEX2(1,2,8)] =    wd0;
                        EM_S[INDEX2(2,2,8)] = 4.*wd0;
                        EM_S[INDEX2(3,2,8)] = 2.*wd0;
                        EM_S[INDEX2(0,3,8)] =    wd0;
                        EM_S[INDEX2(1,3,8)] = 2.*wd0;
                        EM_S[INDEX2(2,3,8)] = 2.*wd0;
                        EM_S[INDEX2(3,3,8)] = 4.*wd0;
                    }
                }
                ///////////////
                // process y //
                ///////////////
                if (add_EM_F) {
                    const Scalar* y_p = y.getSampleDataRO(e, zero);
                    if (y.actsExpanded()) {
                        const Scalar y_0 = y_p[0];
                        const Scalar y_1 = y_p[1];
                        const Scalar y_2 = y_p[2];
                        const Scalar y_3 = y_p[3];
                        const Scalar tmp0 = 6.*w[12][l]*(y_1 + y_2);
                        const Scalar tmp1 = 6.*w[12][l]*(y_0 + y_3);
                        EM_F[0] = tmp0 + 6.*w[10][l]*y_3 + 6.*w[11][l]*y_0;
                        EM_F[1] = tmp1 + 6.*w[10][l]*y_2 + 6.*w[11][l]*y_1;
                        EM_F[2] = tmp1 + 6.*w[10][l]*y_1 + 6.*w[11][l]*y_2;
                        EM_F[3] = tmp0 + 6.*w[10][l]*y_0 + 6.*w[11][l]*y_3;
                    } else { // constant data
                        EM_F[0] = 36.*w[12][l]*y_p[0];
                        EM_F[1] = 36.*w[12][l]*y_p[0];
                        EM_F[2] = 36.*w[12][l]*y_p[0];
                        EM_F[3] = 36.*w[12][l]*y_p[0];
                    }
                }
                domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, add_EM_S, add_EM_F, q, t);
            } // face 4

            if(quaddata->m_faceOffset[5]) {
                if (add_EM_S)
                    fill(EM_S.begin(), EM_S.end(), zero);
                if (add_EM_F) {
                    EM_F[0] = zero;
                    EM_F[1] = zero;
                    EM_F[2] = zero;
                    EM_F[3] = zero;
                }
                ///////////////
                // process d //
                ///////////////
                if (add_EM_S) {
                    const Scalar* d_p = d.getSampleDataRO(e, zero);
                    if (d.actsExpanded()) {
                        const Scalar d_0 = d_p[0];
                        const Scalar d_1 = d_p[1];
                        const Scalar d_2 = d_p[2];
                        const Scalar d_3 = d_p[3];
                        const Scalar tmp0  = w[12][l]*(d_0 + d_1 + d_2 + d_3);
                        const Scalar tmp1  = w[10][l]*(d_1 + d_3);
                        const Scalar tmp2  = w[11][l]*(d_0 + d_2);
                        const Scalar tmp3  = w[10][l]*(d_2 + d_3);
                        const Scalar tmp4  = w[11][l]*(d_0 + d_1);
                        const Scalar tmp5  = w[10][l]*(d_0 + d_1);
                        const Scalar tmp6  = w[11][l]*(d_2 + d_3);
                        const Scalar tmp7  = w[12][l]*(d_1 + d_2);
                        const Scalar tmp8  = w[10][l]*(d_0 + d_2);
                        const Scalar tmp9  = w[11][l]*(d_1 + d_3);
                        const Scalar tmp10 = w[12][l]*(d_0 + d_3);
                        EM_S[INDEX2(4,4,8)] = d_0*w[14][l] + d_3*w[13][l] + tmp7;
                        EM_S[INDEX2(5,4,8)] = tmp3 + tmp4;
                        EM_S[INDEX2(6,4,8)] = tmp1 + tmp2;
                        EM_S[INDEX2(7,4,8)] = tmp0;
                        EM_S[INDEX2(4,5,8)] = tmp3 + tmp4;
                        EM_S[INDEX2(5,5,8)] = d_1*w[14][l] + d_2*w[13][l] + tmp10;
                        EM_S[INDEX2(6,5,8)] = tmp0;
                        EM_S[INDEX2(7,5,8)] = tmp8 + tmp9;
                        EM_S[INDEX2(4,6,8)] = tmp1 + tmp2;
                        EM_S[INDEX2(5,6,8)] = tmp0;
                        EM_S[INDEX2(6,6,8)] = d_1*w[13][l] + d_2*w[14][l] + tmp10;
                        EM_S[INDEX2(7,6,8)] = tmp5 + tmp6;
                        EM_S[INDEX2(4,7,8)] = tmp0;
                        EM_S[INDEX2(5,7,8)] = tmp8 + tmp9;
                        EM_S[INDEX2(6,7,8)] = tmp5 + tmp6;
                        EM_S[INDEX2(7,7,8)] = d_0*w[13][l] + d_3*w[14][l] + tmp7;
                    } else { // constant data
                        const Scalar wd0 = 4.*d_p[0]*w[12][l];
                        EM_S[INDEX2(4,4,8)] = 4.*wd0;
                        EM_S[INDEX2(5,4,8)] = 2.*wd0;
                        EM_S[INDEX2(6,4,8)] = 2.*wd0;
                        EM_S[INDEX2(7,4,8)] =    wd0;
                        EM_S[INDEX2(4,5,8)] = 2.*wd0;
                        EM_S[INDEX2(5,5,8)] = 4.*wd0;
                        EM_S[INDEX2(6,5,8)] =    wd0;
                        EM_S[INDEX2(7,5,8)] = 2.*wd0;
                        EM_S[INDEX2(4,6,8)] = 2.*wd0;
                        EM_S[INDEX2(5,6,8)] =    wd0;
                        EM_S[INDEX2(6,6,8)] = 4.*wd0;
                        EM_S[INDEX2(7,6,8)] = 2.*wd0;
                        EM_S[INDEX2(4,7,8)] =    wd0;
                        EM_S[INDEX2(5,7,8)] = 2.*wd0;
                        EM_S[INDEX2(6,7,8)] = 2.*wd0;
                        EM_S[INDEX2(7,7,8)] = 4.*wd0;
                    }
                }
                ///////////////
                // process y //
                ///////////////
                if (add_EM_F) {
                    const Scalar* y_p = y.getSampleDataRO(e, zero);
                    if (y.actsExpanded()) {
                        const Scalar y_0 = y_p[0];
                        const Scalar y_1 = y_p[1];
                        const Scalar y_2 = y_p[2];
                        const Scalar y_3 = y_p[3];
                        const Scalar tmp0 = 6.*w[12][l]*(y_1 + y_2);
                        const Scalar tmp1 = 6.*w[12][l]*(y_0 + y_3);
                        EM_F[4] = tmp0 + 6.*w[10][l]*y_3 + 6.*w[11][l]*y_0;
                        EM_F[5] = tmp1 + 6.*w[10][l]*y_2 + 6.*w[11][l]*y_1;
                        EM_F[6] = tmp1 + 6.*w[10][l]*y_1 + 6.*w[11][l]*y_2;
                        EM_F[7] = tmp0 + 6.*w[10][l]*y_0 + 6.*w[11][l]*y_3;
                    } else { // constant data
                        EM_F[4] = 36.*w[12][l]*y_p[0];
                        EM_F[5] = 36.*w[12][l]*y_p[0];
                        EM_F[6] = 36.*w[12][l]*y_p[0];
                        EM_F[7] = 36.*w[12][l]*y_p[0];
                    }
                }
                domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F, add_EM_S, add_EM_F, q, t);
            } // face 5
        }
    }
}

/****************************************************************************/
// PDE SINGLE REDUCED
/****************************************************************************/

template<class Scalar>
void DefaultAssembler3D<Scalar>::assemblePDESingleReduced(
                                        AbstractSystemMatrix* mat, Data& rhs,
                                        const Data& A, const Data& B,
                                        const Data& C, const Data& D,
                                        const Data& X, const Data& Y) const
{
    // TODO 
//     const double w6 = m_dx[0]/16;
//     const double w5 = m_dx[1]/16;
//     const double w1 = m_dx[2]/16;
//     const double w14 = m_dx[0]*m_dx[1]/32;
//     const double w13 = m_dx[0]*m_dx[2]/32;
//     const double w12 = m_dx[1]*m_dx[2]/32;
//     const double w18 = m_dx[0]*m_dx[1]*m_dx[2]/64;
//     const double w11 = m_dx[0]*m_dx[1]/(16*m_dx[2]);
//     const double w3 = m_dx[0]*m_dx[2]/(16*m_dx[1]);
//     const double w0 = m_dx[1]*m_dx[2]/(16*m_dx[0]);
//     const dim_t NE0 = m_NE[0];
//     const dim_t NE1 = m_NE[1];
//     const dim_t NE2 = m_NE[2];
//     const bool add_EM_S = (!A.isEmpty() || !B.isEmpty() || !C.isEmpty() || !D.isEmpty());
//     const bool add_EM_F = (!X.isEmpty() || !Y.isEmpty());
//     const Scalar zero = static_cast<Scalar>(0);
//     rhs.requireWrite();

// #pragma omp parallel
//     {
//         vector<Scalar> EM_S(8*8, zero);
//         vector<Scalar> EM_F(8, zero);

//         for (index_t k2_0=0; k2_0<2; k2_0++) { // colouring
// #pragma omp for
//             for (index_t k2=k2_0; k2<NE2; k2+=2) {
//                 for (index_t k1=0; k1<NE1; ++k1) {
//                     for (index_t k0=0; k0<NE0; ++k0)  {
//                         const index_t e = k0 + NE0*k1 + NE0*NE1*k2;
//                         if (add_EM_S)
//                             fill(EM_S.begin(), EM_S.end(), zero);
//                         if (add_EM_F)
//                             fill(EM_F.begin(), EM_F.end(), zero);

//                         ///////////////
//                         // process A //
//                         ///////////////
//                         if (!A.isEmpty()) {
//                             const Scalar* A_p = A.getSampleDataRO(id, zero);
//                             const Scalar Aw00 = A_p[INDEX2(0,0,3)]*w0;
//                             const Scalar Aw10 = A_p[INDEX2(1,0,3)]*w1;
//                             const Scalar Aw20 = A_p[INDEX2(2,0,3)]*w5;
//                             const Scalar Aw01 = A_p[INDEX2(0,1,3)]*w1;
//                             const Scalar Aw11 = A_p[INDEX2(1,1,3)]*w3;
//                             const Scalar Aw21 = A_p[INDEX2(2,1,3)]*w6;
//                             const Scalar Aw02 = A_p[INDEX2(0,2,3)]*w5;
//                             const Scalar Aw12 = A_p[INDEX2(1,2,3)]*w6;
//                             const Scalar Aw22 = A_p[INDEX2(2,2,3)]*w11;
//                             EM_S[INDEX2(0,0,8)]+= Aw00 + Aw01 + Aw02 + Aw10 + Aw11 + Aw12 + Aw20 + Aw21 + Aw22;
//                             EM_S[INDEX2(1,0,8)]+=-Aw00 - Aw01 - Aw02 + Aw10 + Aw11 + Aw12 + Aw20 + Aw21 + Aw22;
//                             EM_S[INDEX2(2,0,8)]+= Aw00 + Aw01 + Aw02 - Aw10 - Aw11 - Aw12 + Aw20 + Aw21 + Aw22;
//                             EM_S[INDEX2(3,0,8)]+=-Aw00 - Aw01 - Aw02 - Aw10 - Aw11 - Aw12 + Aw20 + Aw21 + Aw22;
//                             EM_S[INDEX2(4,0,8)]+= Aw00 + Aw01 + Aw02 + Aw10 + Aw11 + Aw12 - Aw20 - Aw21 - Aw22;
//                             EM_S[INDEX2(5,0,8)]+=-Aw00 - Aw01 - Aw02 + Aw10 + Aw11 + Aw12 - Aw20 - Aw21 - Aw22;
//                             EM_S[INDEX2(6,0,8)]+= Aw00 + Aw01 + Aw02 - Aw10 - Aw11 - Aw12 - Aw20 - Aw21 - Aw22;
//                             EM_S[INDEX2(7,0,8)]+=-Aw00 - Aw01 - Aw02 - Aw10 - Aw11 - Aw12 - Aw20 - Aw21 - Aw22;
//                             EM_S[INDEX2(0,1,8)]+=-Aw00 + Aw01 + Aw02 - Aw10 + Aw11 + Aw12 - Aw20 + Aw21 + Aw22;
//                             EM_S[INDEX2(1,1,8)]+= Aw00 - Aw01 - Aw02 - Aw10 + Aw11 + Aw12 - Aw20 + Aw21 + Aw22;
//                             EM_S[INDEX2(2,1,8)]+=-Aw00 + Aw01 + Aw02 + Aw10 - Aw11 - Aw12 - Aw20 + Aw21 + Aw22;
//                             EM_S[INDEX2(3,1,8)]+= Aw00 - Aw01 - Aw02 + Aw10 - Aw11 - Aw12 - Aw20 + Aw21 + Aw22;
//                             EM_S[INDEX2(4,1,8)]+=-Aw00 + Aw01 + Aw02 - Aw10 + Aw11 + Aw12 + Aw20 - Aw21 - Aw22;
//                             EM_S[INDEX2(5,1,8)]+= Aw00 - Aw01 - Aw02 - Aw10 + Aw11 + Aw12 + Aw20 - Aw21 - Aw22;
//                             EM_S[INDEX2(6,1,8)]+=-Aw00 + Aw01 + Aw02 + Aw10 - Aw11 - Aw12 + Aw20 - Aw21 - Aw22;
//                             EM_S[INDEX2(7,1,8)]+= Aw00 - Aw01 - Aw02 + Aw10 - Aw11 - Aw12 + Aw20 - Aw21 - Aw22;
//                             EM_S[INDEX2(0,2,8)]+= Aw00 - Aw01 + Aw02 + Aw10 - Aw11 + Aw12 + Aw20 - Aw21 + Aw22;
//                             EM_S[INDEX2(1,2,8)]+=-Aw00 + Aw01 - Aw02 + Aw10 - Aw11 + Aw12 + Aw20 - Aw21 + Aw22;
//                             EM_S[INDEX2(2,2,8)]+= Aw00 - Aw01 + Aw02 - Aw10 + Aw11 - Aw12 + Aw20 - Aw21 + Aw22;
//                             EM_S[INDEX2(3,2,8)]+=-Aw00 + Aw01 - Aw02 - Aw10 + Aw11 - Aw12 + Aw20 - Aw21 + Aw22;
//                             EM_S[INDEX2(4,2,8)]+= Aw00 - Aw01 + Aw02 + Aw10 - Aw11 + Aw12 - Aw20 + Aw21 - Aw22;
//                             EM_S[INDEX2(5,2,8)]+=-Aw00 + Aw01 - Aw02 + Aw10 - Aw11 + Aw12 - Aw20 + Aw21 - Aw22;
//                             EM_S[INDEX2(6,2,8)]+= Aw00 - Aw01 + Aw02 - Aw10 + Aw11 - Aw12 - Aw20 + Aw21 - Aw22;
//                             EM_S[INDEX2(7,2,8)]+=-Aw00 + Aw01 - Aw02 - Aw10 + Aw11 - Aw12 - Aw20 + Aw21 - Aw22;
//                             EM_S[INDEX2(0,3,8)]+=-Aw00 - Aw01 + Aw02 - Aw10 - Aw11 + Aw12 - Aw20 - Aw21 + Aw22;
//                             EM_S[INDEX2(1,3,8)]+= Aw00 + Aw01 - Aw02 - Aw10 - Aw11 + Aw12 - Aw20 - Aw21 + Aw22;
//                             EM_S[INDEX2(2,3,8)]+=-Aw00 - Aw01 + Aw02 + Aw10 + Aw11 - Aw12 - Aw20 - Aw21 + Aw22;
//                             EM_S[INDEX2(3,3,8)]+= Aw00 + Aw01 - Aw02 + Aw10 + Aw11 - Aw12 - Aw20 - Aw21 + Aw22;
//                             EM_S[INDEX2(4,3,8)]+=-Aw00 - Aw01 + Aw02 - Aw10 - Aw11 + Aw12 + Aw20 + Aw21 - Aw22;
//                             EM_S[INDEX2(5,3,8)]+= Aw00 + Aw01 - Aw02 - Aw10 - Aw11 + Aw12 + Aw20 + Aw21 - Aw22;
//                             EM_S[INDEX2(6,3,8)]+=-Aw00 - Aw01 + Aw02 + Aw10 + Aw11 - Aw12 + Aw20 + Aw21 - Aw22;
//                             EM_S[INDEX2(7,3,8)]+= Aw00 + Aw01 - Aw02 + Aw10 + Aw11 - Aw12 + Aw20 + Aw21 - Aw22;
//                             EM_S[INDEX2(0,4,8)]+= Aw00 + Aw01 - Aw02 + Aw10 + Aw11 - Aw12 + Aw20 + Aw21 - Aw22;
//                             EM_S[INDEX2(1,4,8)]+=-Aw00 - Aw01 + Aw02 + Aw10 + Aw11 - Aw12 + Aw20 + Aw21 - Aw22;
//                             EM_S[INDEX2(2,4,8)]+= Aw00 + Aw01 - Aw02 - Aw10 - Aw11 + Aw12 + Aw20 + Aw21 - Aw22;
//                             EM_S[INDEX2(3,4,8)]+=-Aw00 - Aw01 + Aw02 - Aw10 - Aw11 + Aw12 + Aw20 + Aw21 - Aw22;
//                             EM_S[INDEX2(4,4,8)]+= Aw00 + Aw01 - Aw02 + Aw10 + Aw11 - Aw12 - Aw20 - Aw21 + Aw22;
//                             EM_S[INDEX2(5,4,8)]+=-Aw00 - Aw01 + Aw02 + Aw10 + Aw11 - Aw12 - Aw20 - Aw21 + Aw22;
//                             EM_S[INDEX2(6,4,8)]+= Aw00 + Aw01 - Aw02 - Aw10 - Aw11 + Aw12 - Aw20 - Aw21 + Aw22;
//                             EM_S[INDEX2(7,4,8)]+=-Aw00 - Aw01 + Aw02 - Aw10 - Aw11 + Aw12 - Aw20 - Aw21 + Aw22;
//                             EM_S[INDEX2(0,5,8)]+=-Aw00 + Aw01 - Aw02 - Aw10 + Aw11 - Aw12 - Aw20 + Aw21 - Aw22;
//                             EM_S[INDEX2(1,5,8)]+= Aw00 - Aw01 + Aw02 - Aw10 + Aw11 - Aw12 - Aw20 + Aw21 - Aw22;
//                             EM_S[INDEX2(2,5,8)]+=-Aw00 + Aw01 - Aw02 + Aw10 - Aw11 + Aw12 - Aw20 + Aw21 - Aw22;
//                             EM_S[INDEX2(3,5,8)]+= Aw00 - Aw01 + Aw02 + Aw10 - Aw11 + Aw12 - Aw20 + Aw21 - Aw22;
//                             EM_S[INDEX2(4,5,8)]+=-Aw00 + Aw01 - Aw02 - Aw10 + Aw11 - Aw12 + Aw20 - Aw21 + Aw22;
//                             EM_S[INDEX2(5,5,8)]+= Aw00 - Aw01 + Aw02 - Aw10 + Aw11 - Aw12 + Aw20 - Aw21 + Aw22;
//                             EM_S[INDEX2(6,5,8)]+=-Aw00 + Aw01 - Aw02 + Aw10 - Aw11 + Aw12 + Aw20 - Aw21 + Aw22;
//                             EM_S[INDEX2(7,5,8)]+= Aw00 - Aw01 + Aw02 + Aw10 - Aw11 + Aw12 + Aw20 - Aw21 + Aw22;
//                             EM_S[INDEX2(0,6,8)]+= Aw00 - Aw01 - Aw02 + Aw10 - Aw11 - Aw12 + Aw20 - Aw21 - Aw22;
//                             EM_S[INDEX2(1,6,8)]+=-Aw00 + Aw01 + Aw02 + Aw10 - Aw11 - Aw12 + Aw20 - Aw21 - Aw22;
//                             EM_S[INDEX2(2,6,8)]+= Aw00 - Aw01 - Aw02 - Aw10 + Aw11 + Aw12 + Aw20 - Aw21 - Aw22;
//                             EM_S[INDEX2(3,6,8)]+=-Aw00 + Aw01 + Aw02 - Aw10 + Aw11 + Aw12 + Aw20 - Aw21 - Aw22;
//                             EM_S[INDEX2(4,6,8)]+= Aw00 - Aw01 - Aw02 + Aw10 - Aw11 - Aw12 - Aw20 + Aw21 + Aw22;
//                             EM_S[INDEX2(5,6,8)]+=-Aw00 + Aw01 + Aw02 + Aw10 - Aw11 - Aw12 - Aw20 + Aw21 + Aw22;
//                             EM_S[INDEX2(6,6,8)]+= Aw00 - Aw01 - Aw02 - Aw10 + Aw11 + Aw12 - Aw20 + Aw21 + Aw22;
//                             EM_S[INDEX2(7,6,8)]+=-Aw00 + Aw01 + Aw02 - Aw10 + Aw11 + Aw12 - Aw20 + Aw21 + Aw22;
//                             EM_S[INDEX2(0,7,8)]+=-Aw00 - Aw01 - Aw02 - Aw10 - Aw11 - Aw12 - Aw20 - Aw21 - Aw22;
//                             EM_S[INDEX2(1,7,8)]+= Aw00 + Aw01 + Aw02 - Aw10 - Aw11 - Aw12 - Aw20 - Aw21 - Aw22;
//                             EM_S[INDEX2(2,7,8)]+=-Aw00 - Aw01 - Aw02 + Aw10 + Aw11 + Aw12 - Aw20 - Aw21 - Aw22;
//                             EM_S[INDEX2(3,7,8)]+= Aw00 + Aw01 + Aw02 + Aw10 + Aw11 + Aw12 - Aw20 - Aw21 - Aw22;
//                             EM_S[INDEX2(4,7,8)]+=-Aw00 - Aw01 - Aw02 - Aw10 - Aw11 - Aw12 + Aw20 + Aw21 + Aw22;
//                             EM_S[INDEX2(5,7,8)]+= Aw00 + Aw01 + Aw02 - Aw10 - Aw11 - Aw12 + Aw20 + Aw21 + Aw22;
//                             EM_S[INDEX2(6,7,8)]+=-Aw00 - Aw01 - Aw02 + Aw10 + Aw11 + Aw12 + Aw20 + Aw21 + Aw22;
//                             EM_S[INDEX2(7,7,8)]+= Aw00 + Aw01 + Aw02 + Aw10 + Aw11 + Aw12 + Aw20 + Aw21 + Aw22;
//                         }
//                         ///////////////
//                         // process B //
//                         ///////////////
//                         if (!B.isEmpty()) {
//                             const Scalar* B_p = B.getSampleDataRO(id, zero);
//                             const Scalar wB0 = B_p[0]*w12;
//                             const Scalar wB1 = B_p[1]*w13;
//                             const Scalar wB2 = B_p[2]*w14;
//                             EM_S[INDEX2(0,0,8)]+=-wB0 - wB1 - wB2;
//                             EM_S[INDEX2(1,0,8)]+= wB0 - wB1 - wB2;
//                             EM_S[INDEX2(2,0,8)]+=-wB0 + wB1 - wB2;
//                             EM_S[INDEX2(3,0,8)]+= wB0 + wB1 - wB2;
//                             EM_S[INDEX2(4,0,8)]+=-wB0 - wB1 + wB2;
//                             EM_S[INDEX2(5,0,8)]+= wB0 - wB1 + wB2;
//                             EM_S[INDEX2(6,0,8)]+=-wB0 + wB1 + wB2;
//                             EM_S[INDEX2(7,0,8)]+= wB0 + wB1 + wB2;
//                             EM_S[INDEX2(0,1,8)]+=-wB0 - wB1 - wB2;
//                             EM_S[INDEX2(1,1,8)]+= wB0 - wB1 - wB2;
//                             EM_S[INDEX2(2,1,8)]+=-wB0 + wB1 - wB2;
//                             EM_S[INDEX2(3,1,8)]+= wB0 + wB1 - wB2;
//                             EM_S[INDEX2(4,1,8)]+=-wB0 - wB1 + wB2;
//                             EM_S[INDEX2(5,1,8)]+= wB0 - wB1 + wB2;
//                             EM_S[INDEX2(6,1,8)]+=-wB0 + wB1 + wB2;
//                             EM_S[INDEX2(7,1,8)]+= wB0 + wB1 + wB2;
//                             EM_S[INDEX2(0,2,8)]+=-wB0 - wB1 - wB2;
//                             EM_S[INDEX2(1,2,8)]+= wB0 - wB1 - wB2;
//                             EM_S[INDEX2(2,2,8)]+=-wB0 + wB1 - wB2;
//                             EM_S[INDEX2(3,2,8)]+= wB0 + wB1 - wB2;
//                             EM_S[INDEX2(4,2,8)]+=-wB0 - wB1 + wB2;
//                             EM_S[INDEX2(5,2,8)]+= wB0 - wB1 + wB2;
//                             EM_S[INDEX2(6,2,8)]+=-wB0 + wB1 + wB2;
//                             EM_S[INDEX2(7,2,8)]+= wB0 + wB1 + wB2;
//                             EM_S[INDEX2(0,3,8)]+=-wB0 - wB1 - wB2;
//                             EM_S[INDEX2(1,3,8)]+= wB0 - wB1 - wB2;
//                             EM_S[INDEX2(2,3,8)]+=-wB0 + wB1 - wB2;
//                             EM_S[INDEX2(3,3,8)]+= wB0 + wB1 - wB2;
//                             EM_S[INDEX2(4,3,8)]+=-wB0 - wB1 + wB2;
//                             EM_S[INDEX2(5,3,8)]+= wB0 - wB1 + wB2;
//                             EM_S[INDEX2(6,3,8)]+=-wB0 + wB1 + wB2;
//                             EM_S[INDEX2(7,3,8)]+= wB0 + wB1 + wB2;
//                             EM_S[INDEX2(0,4,8)]+=-wB0 - wB1 - wB2;
//                             EM_S[INDEX2(1,4,8)]+= wB0 - wB1 - wB2;
//                             EM_S[INDEX2(2,4,8)]+=-wB0 + wB1 - wB2;
//                             EM_S[INDEX2(3,4,8)]+= wB0 + wB1 - wB2;
//                             EM_S[INDEX2(4,4,8)]+=-wB0 - wB1 + wB2;
//                             EM_S[INDEX2(5,4,8)]+= wB0 - wB1 + wB2;
//                             EM_S[INDEX2(6,4,8)]+=-wB0 + wB1 + wB2;
//                             EM_S[INDEX2(7,4,8)]+= wB0 + wB1 + wB2;
//                             EM_S[INDEX2(0,5,8)]+=-wB0 - wB1 - wB2;
//                             EM_S[INDEX2(1,5,8)]+= wB0 - wB1 - wB2;
//                             EM_S[INDEX2(2,5,8)]+=-wB0 + wB1 - wB2;
//                             EM_S[INDEX2(3,5,8)]+= wB0 + wB1 - wB2;
//                             EM_S[INDEX2(4,5,8)]+=-wB0 - wB1 + wB2;
//                             EM_S[INDEX2(5,5,8)]+= wB0 - wB1 + wB2;
//                             EM_S[INDEX2(6,5,8)]+=-wB0 + wB1 + wB2;
//                             EM_S[INDEX2(7,5,8)]+= wB0 + wB1 + wB2;
//                             EM_S[INDEX2(0,6,8)]+=-wB0 - wB1 - wB2;
//                             EM_S[INDEX2(1,6,8)]+= wB0 - wB1 - wB2;
//                             EM_S[INDEX2(2,6,8)]+=-wB0 + wB1 - wB2;
//                             EM_S[INDEX2(3,6,8)]+= wB0 + wB1 - wB2;
//                             EM_S[INDEX2(4,6,8)]+=-wB0 - wB1 + wB2;
//                             EM_S[INDEX2(5,6,8)]+= wB0 - wB1 + wB2;
//                             EM_S[INDEX2(6,6,8)]+=-wB0 + wB1 + wB2;
//                             EM_S[INDEX2(7,6,8)]+= wB0 + wB1 + wB2;
//                             EM_S[INDEX2(0,7,8)]+=-wB0 - wB1 - wB2;
//                             EM_S[INDEX2(1,7,8)]+= wB0 - wB1 - wB2;
//                             EM_S[INDEX2(2,7,8)]+=-wB0 + wB1 - wB2;
//                             EM_S[INDEX2(3,7,8)]+= wB0 + wB1 - wB2;
//                             EM_S[INDEX2(4,7,8)]+=-wB0 - wB1 + wB2;
//                             EM_S[INDEX2(5,7,8)]+= wB0 - wB1 + wB2;
//                             EM_S[INDEX2(6,7,8)]+=-wB0 + wB1 + wB2;
//                             EM_S[INDEX2(7,7,8)]+= wB0 + wB1 + wB2;
//                         }
//                         ///////////////
//                         // process C //
//                         ///////////////
//                         if (!C.isEmpty()) {
//                             const Scalar* C_p = C.getSampleDataRO(id, zero);
//                             const Scalar wC0 = C_p[0]*w12;
//                             const Scalar wC1 = C_p[1]*w13;
//                             const Scalar wC2 = C_p[2]*w14;
//                             EM_S[INDEX2(0,0,8)]+=-wC0 - wC1 - wC2;
//                             EM_S[INDEX2(1,0,8)]+=-wC0 - wC1 - wC2;
//                             EM_S[INDEX2(2,0,8)]+=-wC0 - wC1 - wC2;
//                             EM_S[INDEX2(3,0,8)]+=-wC0 - wC1 - wC2;
//                             EM_S[INDEX2(4,0,8)]+=-wC0 - wC1 - wC2;
//                             EM_S[INDEX2(5,0,8)]+=-wC0 - wC1 - wC2;
//                             EM_S[INDEX2(6,0,8)]+=-wC0 - wC1 - wC2;
//                             EM_S[INDEX2(7,0,8)]+=-wC0 - wC1 - wC2;
//                             EM_S[INDEX2(0,1,8)]+= wC0 - wC1 - wC2;
//                             EM_S[INDEX2(1,1,8)]+= wC0 - wC1 - wC2;
//                             EM_S[INDEX2(2,1,8)]+= wC0 - wC1 - wC2;
//                             EM_S[INDEX2(3,1,8)]+= wC0 - wC1 - wC2;
//                             EM_S[INDEX2(4,1,8)]+= wC0 - wC1 - wC2;
//                             EM_S[INDEX2(5,1,8)]+= wC0 - wC1 - wC2;
//                             EM_S[INDEX2(6,1,8)]+= wC0 - wC1 - wC2;
//                             EM_S[INDEX2(7,1,8)]+= wC0 - wC1 - wC2;
//                             EM_S[INDEX2(0,2,8)]+=-wC0 + wC1 - wC2;
//                             EM_S[INDEX2(1,2,8)]+=-wC0 + wC1 - wC2;
//                             EM_S[INDEX2(2,2,8)]+=-wC0 + wC1 - wC2;
//                             EM_S[INDEX2(3,2,8)]+=-wC0 + wC1 - wC2;
//                             EM_S[INDEX2(4,2,8)]+=-wC0 + wC1 - wC2;
//                             EM_S[INDEX2(5,2,8)]+=-wC0 + wC1 - wC2;
//                             EM_S[INDEX2(6,2,8)]+=-wC0 + wC1 - wC2;
//                             EM_S[INDEX2(7,2,8)]+=-wC0 + wC1 - wC2;
//                             EM_S[INDEX2(0,3,8)]+= wC0 + wC1 - wC2;
//                             EM_S[INDEX2(1,3,8)]+= wC0 + wC1 - wC2;
//                             EM_S[INDEX2(2,3,8)]+= wC0 + wC1 - wC2;
//                             EM_S[INDEX2(3,3,8)]+= wC0 + wC1 - wC2;
//                             EM_S[INDEX2(4,3,8)]+= wC0 + wC1 - wC2;
//                             EM_S[INDEX2(5,3,8)]+= wC0 + wC1 - wC2;
//                             EM_S[INDEX2(6,3,8)]+= wC0 + wC1 - wC2;
//                             EM_S[INDEX2(7,3,8)]+= wC0 + wC1 - wC2;
//                             EM_S[INDEX2(0,4,8)]+=-wC0 - wC1 + wC2;
//                             EM_S[INDEX2(1,4,8)]+=-wC0 - wC1 + wC2;
//                             EM_S[INDEX2(2,4,8)]+=-wC0 - wC1 + wC2;
//                             EM_S[INDEX2(3,4,8)]+=-wC0 - wC1 + wC2;
//                             EM_S[INDEX2(4,4,8)]+=-wC0 - wC1 + wC2;
//                             EM_S[INDEX2(5,4,8)]+=-wC0 - wC1 + wC2;
//                             EM_S[INDEX2(6,4,8)]+=-wC0 - wC1 + wC2;
//                             EM_S[INDEX2(7,4,8)]+=-wC0 - wC1 + wC2;
//                             EM_S[INDEX2(0,5,8)]+= wC0 - wC1 + wC2;
//                             EM_S[INDEX2(1,5,8)]+= wC0 - wC1 + wC2;
//                             EM_S[INDEX2(2,5,8)]+= wC0 - wC1 + wC2;
//                             EM_S[INDEX2(3,5,8)]+= wC0 - wC1 + wC2;
//                             EM_S[INDEX2(4,5,8)]+= wC0 - wC1 + wC2;
//                             EM_S[INDEX2(5,5,8)]+= wC0 - wC1 + wC2;
//                             EM_S[INDEX2(6,5,8)]+= wC0 - wC1 + wC2;
//                             EM_S[INDEX2(7,5,8)]+= wC0 - wC1 + wC2;
//                             EM_S[INDEX2(0,6,8)]+=-wC0 + wC1 + wC2;
//                             EM_S[INDEX2(1,6,8)]+=-wC0 + wC1 + wC2;
//                             EM_S[INDEX2(2,6,8)]+=-wC0 + wC1 + wC2;
//                             EM_S[INDEX2(3,6,8)]+=-wC0 + wC1 + wC2;
//                             EM_S[INDEX2(4,6,8)]+=-wC0 + wC1 + wC2;
//                             EM_S[INDEX2(5,6,8)]+=-wC0 + wC1 + wC2;
//                             EM_S[INDEX2(6,6,8)]+=-wC0 + wC1 + wC2;
//                             EM_S[INDEX2(7,6,8)]+=-wC0 + wC1 + wC2;
//                             EM_S[INDEX2(0,7,8)]+= wC0 + wC1 + wC2;
//                             EM_S[INDEX2(1,7,8)]+= wC0 + wC1 + wC2;
//                             EM_S[INDEX2(2,7,8)]+= wC0 + wC1 + wC2;
//                             EM_S[INDEX2(3,7,8)]+= wC0 + wC1 + wC2;
//                             EM_S[INDEX2(4,7,8)]+= wC0 + wC1 + wC2;
//                             EM_S[INDEX2(5,7,8)]+= wC0 + wC1 + wC2;
//                             EM_S[INDEX2(6,7,8)]+= wC0 + wC1 + wC2;
//                             EM_S[INDEX2(7,7,8)]+= wC0 + wC1 + wC2;
//                         }
//                         ///////////////
//                         // process D //
//                         ///////////////
//                         if (!D.isEmpty()) {
//                             const Scalar* D_p = D.getSampleDataRO(id, zero);
//                             EM_S[INDEX2(0,0,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(1,0,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(2,0,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(3,0,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(4,0,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(5,0,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(6,0,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(7,0,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(0,1,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(1,1,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(2,1,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(3,1,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(4,1,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(5,1,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(6,1,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(7,1,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(0,2,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(1,2,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(2,2,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(3,2,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(4,2,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(5,2,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(6,2,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(7,2,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(0,3,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(1,3,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(2,3,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(3,3,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(4,3,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(5,3,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(6,3,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(7,3,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(0,4,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(1,4,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(2,4,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(3,4,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(4,4,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(5,4,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(6,4,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(7,4,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(0,5,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(1,5,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(2,5,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(3,5,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(4,5,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(5,5,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(6,5,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(7,5,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(0,6,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(1,6,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(2,6,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(3,6,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(4,6,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(5,6,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(6,6,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(7,6,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(0,7,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(1,7,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(2,7,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(3,7,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(4,7,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(5,7,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(6,7,8)]+=D_p[0]*w18;
//                             EM_S[INDEX2(7,7,8)]+=D_p[0]*w18;
//                         }
//                         ///////////////
//                         // process X //
//                         ///////////////
//                         if (!X.isEmpty()) {
//                             const Scalar* X_p = X.getSampleDataRO(id, zero);
//                             const Scalar wX0 = 8.*X_p[0]*w12;
//                             const Scalar wX1 = 8.*X_p[1]*w13;
//                             const Scalar wX2 = 8.*X_p[2]*w14;
//                             EM_F[0]+=-wX0 - wX1 - wX2;
//                             EM_F[1]+= wX0 - wX1 - wX2;
//                             EM_F[2]+=-wX0 + wX1 - wX2;
//                             EM_F[3]+= wX0 + wX1 - wX2;
//                             EM_F[4]+=-wX0 - wX1 + wX2;
//                             EM_F[5]+= wX0 - wX1 + wX2;
//                             EM_F[6]+=-wX0 + wX1 + wX2;
//                             EM_F[7]+= wX0 + wX1 + wX2;
//                         }
//                         ///////////////
//                         // process Y //
//                         ///////////////
//                         if (!Y.isEmpty()) {
//                             const Scalar* Y_p = Y.getSampleDataRO(id, zero);
//                             EM_F[0]+=8.*Y_p[0]*w18;
//                             EM_F[1]+=8.*Y_p[0]*w18;
//                             EM_F[2]+=8.*Y_p[0]*w18;
//                             EM_F[3]+=8.*Y_p[0]*w18;
//                             EM_F[4]+=8.*Y_p[0]*w18;
//                             EM_F[5]+=8.*Y_p[0]*w18;
//                             EM_F[6]+=8.*Y_p[0]*w18;
//                             EM_F[7]+=8.*Y_p[0]*w18;
//                         }

//                         // add to matrix (if add_EM_S) and RHS (if add_EM_F)
//                         const index_t firstNode=m_NN[0]*m_NN[1]*k2+m_NN[0]*k1+k0;
//                         domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F,
//                                                add_EM_S, add_EM_F, firstNode);
//                     } // end k0 loop
//                 } // end k1 loop
//             } // end k2 loop
//         } // end of colouring
//     } // end of parallel region
}

/****************************************************************************/
// PDE SINGLE REDUCED BOUNDARY
/****************************************************************************/

template<class Scalar>
void DefaultAssembler3D<Scalar>::assemblePDEBoundarySingleReduced(
                                        AbstractSystemMatrix* mat, Data& rhs,
                                        const Data& d, const Data& y) const
{
//     const double w0 = m_dx[0]*m_dx[1]/16;
//     const double w1 = m_dx[0]*m_dx[2]/16;
//     const double w2 = m_dx[1]*m_dx[2]/16;
//     const dim_t NE0 = m_NE[0];
//     const dim_t NE1 = m_NE[1];
//     const dim_t NE2 = m_NE[2];
//     const bool add_EM_S = !d.isEmpty();
//     const bool add_EM_F = !y.isEmpty();
//     const Scalar zero = static_cast<Scalar>(0);
//     rhs.requireWrite();

// #pragma omp parallel
//     {
//         vector<Scalar> EM_S(8*8);
//         vector<Scalar> EM_F(8);

//         if (domain->m_faceOffset[0] > -1) {
//             if (add_EM_S)
//                 fill(EM_S.begin(), EM_S.end(), zero);
//             if (add_EM_F) {
//                 EM_F[1] = zero;
//                 EM_F[3] = zero;
//                 EM_F[5] = zero;
//                 EM_F[7] = zero;
//             }

//             for (index_t k2_0=0; k2_0<2; k2_0++) { // colouring
// #pragma omp for
//                 for (index_t k2=k2_0; k2<NE2; k2+=2) {
//                     for (index_t k1=0; k1<NE1; ++k1) {
//                         const index_t e = INDEX2(k1,k2,NE1);
//                         ///////////////
//                         // process d //
//                         ///////////////
//                         if (add_EM_S) {
//                             const Scalar* d_p = d.getSampleDataRO(id, zero);
//                             EM_S[INDEX2(0,0,8)] = d_p[0]*w2;
//                             EM_S[INDEX2(2,0,8)] = d_p[0]*w2;
//                             EM_S[INDEX2(4,0,8)] = d_p[0]*w2;
//                             EM_S[INDEX2(6,0,8)] = d_p[0]*w2;
//                             EM_S[INDEX2(0,2,8)] = d_p[0]*w2;
//                             EM_S[INDEX2(2,2,8)] = d_p[0]*w2;
//                             EM_S[INDEX2(4,2,8)] = d_p[0]*w2;
//                             EM_S[INDEX2(6,2,8)] = d_p[0]*w2;
//                             EM_S[INDEX2(0,4,8)] = d_p[0]*w2;
//                             EM_S[INDEX2(2,4,8)] = d_p[0]*w2;
//                             EM_S[INDEX2(4,4,8)] = d_p[0]*w2;
//                             EM_S[INDEX2(6,4,8)] = d_p[0]*w2;
//                             EM_S[INDEX2(0,6,8)] = d_p[0]*w2;
//                             EM_S[INDEX2(2,6,8)] = d_p[0]*w2;
//                             EM_S[INDEX2(4,6,8)] = d_p[0]*w2;
//                             EM_S[INDEX2(6,6,8)] = d_p[0]*w2;
//                         }
//                         ///////////////
//                         // process y //
//                         ///////////////
//                         if (add_EM_F) {
//                             const Scalar* y_p = y.getSampleDataRO(id, zero);
//                             EM_F[0] = 4.*w2*y_p[0];
//                             EM_F[2] = 4.*w2*y_p[0];
//                             EM_F[4] = 4.*w2*y_p[0];
//                             EM_F[6] = 4.*w2*y_p[0];
//                         }
//                         const index_t firstNode=m_NN[0]*m_NN[1]*k2+m_NN[0]*k1;
//                         domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F,
//                                                add_EM_S, add_EM_F, firstNode);
//                     } // k1 loop
//                 } // k2 loop
//             } // colouring
//         } // face 0

//         if (domain->m_faceOffset[1] > -1) {
//             if (add_EM_S)
//                 fill(EM_S.begin(), EM_S.end(), zero);
//             if (add_EM_F) {
//                 EM_F[0] = zero;
//                 EM_F[2] = zero;
//                 EM_F[4] = zero;
//                 EM_F[6] = zero;
//             }

//             for (index_t k2_0=0; k2_0<2; k2_0++) { // colouring
// #pragma omp for
//                 for (index_t k2=k2_0; k2<NE2; k2+=2) {
//                     for (index_t k1=0; k1<NE1; ++k1) {
//                         const index_t e = domain->m_faceOffset[1]+INDEX2(k1,k2,NE1);
//                         ///////////////
//                         // process d //
//                         ///////////////
//                         if (add_EM_S) {
//                             const Scalar* d_p = d.getSampleDataRO(id, zero);
//                             EM_S[INDEX2(1,1,8)] = d_p[0]*w2;
//                             EM_S[INDEX2(3,1,8)] = d_p[0]*w2;
//                             EM_S[INDEX2(5,1,8)] = d_p[0]*w2;
//                             EM_S[INDEX2(7,1,8)] = d_p[0]*w2;
//                             EM_S[INDEX2(1,3,8)] = d_p[0]*w2;
//                             EM_S[INDEX2(3,3,8)] = d_p[0]*w2;
//                             EM_S[INDEX2(5,3,8)] = d_p[0]*w2;
//                             EM_S[INDEX2(7,3,8)] = d_p[0]*w2;
//                             EM_S[INDEX2(1,5,8)] = d_p[0]*w2;
//                             EM_S[INDEX2(3,5,8)] = d_p[0]*w2;
//                             EM_S[INDEX2(5,5,8)] = d_p[0]*w2;
//                             EM_S[INDEX2(7,5,8)] = d_p[0]*w2;
//                             EM_S[INDEX2(1,7,8)] = d_p[0]*w2;
//                             EM_S[INDEX2(3,7,8)] = d_p[0]*w2;
//                             EM_S[INDEX2(5,7,8)] = d_p[0]*w2;
//                             EM_S[INDEX2(7,7,8)] = d_p[0]*w2;
//                         }
//                         ///////////////
//                         // process y //
//                         ///////////////
//                         if (add_EM_F) {
//                             const Scalar* y_p = y.getSampleDataRO(id, zero);
//                             EM_F[1] = 4.*w2*y_p[0];
//                             EM_F[3] = 4.*w2*y_p[0];
//                             EM_F[5] = 4.*w2*y_p[0];
//                             EM_F[7] = 4.*w2*y_p[0];
//                         }
//                         const index_t firstNode=m_NN[0]*m_NN[1]*k2+m_NN[0]*(k1+1)-2;
//                         domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F,
//                                                add_EM_S, add_EM_F, firstNode);
//                     } // k1 loop
//                 } // k2 loop
//             } // colouring
//         } // face 1

//         if (domain->m_faceOffset[2] > -1) {
//             if (add_EM_S)
//                 fill(EM_S.begin(), EM_S.end(), zero);
//             if (add_EM_F) {
//                 EM_F[2] = zero;
//                 EM_F[3] = zero;
//                 EM_F[6] = zero;
//                 EM_F[7] = zero;
//             }

//             for (index_t k2_0=0; k2_0<2; k2_0++) { // colouring
// #pragma omp for
//                 for (index_t k2=k2_0; k2<NE2; k2+=2) {
//                     for (index_t k0=0; k0<NE0; ++k0) {
//                         const index_t e = domain->m_faceOffset[2]+INDEX2(k0,k2,NE0);
//                         ///////////////
//                         // process d //
//                         ///////////////
//                         if (add_EM_S) {
//                             const Scalar* d_p = d.getSampleDataRO(id, zero);
//                             EM_S[INDEX2(0,0,8)] = d_p[0]*w1;
//                             EM_S[INDEX2(1,0,8)] = d_p[0]*w1;
//                             EM_S[INDEX2(4,0,8)] = d_p[0]*w1;
//                             EM_S[INDEX2(5,0,8)] = d_p[0]*w1;
//                             EM_S[INDEX2(0,1,8)] = d_p[0]*w1;
//                             EM_S[INDEX2(1,1,8)] = d_p[0]*w1;
//                             EM_S[INDEX2(4,1,8)] = d_p[0]*w1;
//                             EM_S[INDEX2(5,1,8)] = d_p[0]*w1;
//                             EM_S[INDEX2(0,4,8)] = d_p[0]*w1;
//                             EM_S[INDEX2(1,4,8)] = d_p[0]*w1;
//                             EM_S[INDEX2(4,4,8)] = d_p[0]*w1;
//                             EM_S[INDEX2(5,4,8)] = d_p[0]*w1;
//                             EM_S[INDEX2(0,5,8)] = d_p[0]*w1;
//                             EM_S[INDEX2(1,5,8)] = d_p[0]*w1;
//                             EM_S[INDEX2(4,5,8)] = d_p[0]*w1;
//                             EM_S[INDEX2(5,5,8)] = d_p[0]*w1;
//                         }
//                         ///////////////
//                         // process y //
//                         ///////////////
//                         if (add_EM_F) {
//                             const Scalar* y_p = y.getSampleDataRO(id, zero);
//                             EM_F[0] = 4.*w1*y_p[0];
//                             EM_F[1] = 4.*w1*y_p[0];
//                             EM_F[4] = 4.*w1*y_p[0];
//                             EM_F[5] = 4.*w1*y_p[0];
//                         }
//                         const index_t firstNode=m_NN[0]*m_NN[1]*k2+k0;
//                         domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F,
//                                                add_EM_S, add_EM_F, firstNode);
//                     } // k0 loop
//                 } // k2 loop
//             } // colouring
//         } // face 2

//         if (domain->m_faceOffset[3] > -1) {
//             if (add_EM_S)
//                 fill(EM_S.begin(), EM_S.end(), zero);
//             if (add_EM_F) {
//                 EM_F[0] = zero;
//                 EM_F[1] = zero;
//                 EM_F[4] = zero;
//                 EM_F[5] = zero;
//             }

//             for (index_t k2_0=0; k2_0<2; k2_0++) { // colouring
// #pragma omp for
//                 for (index_t k2=k2_0; k2<NE2; k2+=2) {
//                     for (index_t k0=0; k0<NE0; ++k0) {
//                         const index_t e = domain->m_faceOffset[3]+INDEX2(k0,k2,NE0);
//                         ///////////////
//                         // process d //
//                         ///////////////
//                         if (add_EM_S) {
//                             const Scalar* d_p = d.getSampleDataRO(id, zero);
//                             EM_S[INDEX2(2,2,8)] = d_p[0]*w1;
//                             EM_S[INDEX2(3,2,8)] = d_p[0]*w1;
//                             EM_S[INDEX2(6,2,8)] = d_p[0]*w1;
//                             EM_S[INDEX2(7,2,8)] = d_p[0]*w1;
//                             EM_S[INDEX2(2,3,8)] = d_p[0]*w1;
//                             EM_S[INDEX2(3,3,8)] = d_p[0]*w1;
//                             EM_S[INDEX2(6,3,8)] = d_p[0]*w1;
//                             EM_S[INDEX2(7,3,8)] = d_p[0]*w1;
//                             EM_S[INDEX2(2,6,8)] = d_p[0]*w1;
//                             EM_S[INDEX2(3,6,8)] = d_p[0]*w1;
//                             EM_S[INDEX2(6,6,8)] = d_p[0]*w1;
//                             EM_S[INDEX2(7,6,8)] = d_p[0]*w1;
//                             EM_S[INDEX2(2,7,8)] = d_p[0]*w1;
//                             EM_S[INDEX2(3,7,8)] = d_p[0]*w1;
//                             EM_S[INDEX2(6,7,8)] = d_p[0]*w1;
//                             EM_S[INDEX2(7,7,8)] = d_p[0]*w1;
//                         }
//                         ///////////////
//                         // process y //
//                         ///////////////
//                         if (add_EM_F) {
//                             const Scalar* y_p = y.getSampleDataRO(id, zero);
//                             EM_F[2] = 4.*w1*y_p[0];
//                             EM_F[3] = 4.*w1*y_p[0];
//                             EM_F[6] = 4.*w1*y_p[0];
//                             EM_F[7] = 4.*w1*y_p[0];
//                         }
//                         const index_t firstNode=m_NN[0]*m_NN[1]*k2+m_NN[0]*(m_NN[1]-2)+k0;
//                         domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F,
//                                                add_EM_S, add_EM_F, firstNode);
//                     } // k0 loop
//                 } // k2 loop
//             } // colouring
//         } // face 3

//         if (domain->m_faceOffset[4] > -1) {
//             if (add_EM_S)
//                 fill(EM_S.begin(), EM_S.end(), zero);
//             if (add_EM_F) {
//                 EM_F[4] = zero;
//                 EM_F[5] = zero;
//                 EM_F[6] = zero;
//                 EM_F[7] = zero;
//             }

//             for (index_t k1_0=0; k1_0<2; k1_0++) { // colouring
// #pragma omp for
//                 for (index_t k1=k1_0; k1<NE1; k1+=2) {
//                     for (index_t k0=0; k0<NE0; ++k0) {
//                         const index_t e = domain->m_faceOffset[4]+INDEX2(k0,k1,NE0);
//                         ///////////////
//                         // process d //
//                         ///////////////
//                         if (add_EM_S) {
//                             const Scalar* d_p = d.getSampleDataRO(id, zero);
//                             EM_S[INDEX2(0,0,8)] = d_p[0]*w0;
//                             EM_S[INDEX2(1,0,8)] = d_p[0]*w0;
//                             EM_S[INDEX2(2,0,8)] = d_p[0]*w0;
//                             EM_S[INDEX2(3,0,8)] = d_p[0]*w0;
//                             EM_S[INDEX2(0,1,8)] = d_p[0]*w0;
//                             EM_S[INDEX2(1,1,8)] = d_p[0]*w0;
//                             EM_S[INDEX2(2,1,8)] = d_p[0]*w0;
//                             EM_S[INDEX2(3,1,8)] = d_p[0]*w0;
//                             EM_S[INDEX2(0,2,8)] = d_p[0]*w0;
//                             EM_S[INDEX2(1,2,8)] = d_p[0]*w0;
//                             EM_S[INDEX2(2,2,8)] = d_p[0]*w0;
//                             EM_S[INDEX2(3,2,8)] = d_p[0]*w0;
//                             EM_S[INDEX2(0,3,8)] = d_p[0]*w0;
//                             EM_S[INDEX2(1,3,8)] = d_p[0]*w0;
//                             EM_S[INDEX2(2,3,8)] = d_p[0]*w0;
//                             EM_S[INDEX2(3,3,8)] = d_p[0]*w0;
//                         }
//                         ///////////////
//                         // process y //
//                         ///////////////
//                         if (add_EM_F) {
//                             const Scalar* y_p = y.getSampleDataRO(id, zero);
//                             EM_F[0] = 4.*w0*y_p[0];
//                             EM_F[1] = 4.*w0*y_p[0];
//                             EM_F[2] = 4.*w0*y_p[0];
//                             EM_F[3] = 4.*w0*y_p[0];
//                         }
//                         const index_t firstNode=m_NN[0]*k1+k0;
//                         domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F,
//                                                add_EM_S, add_EM_F, firstNode);
//                     } // k0 loop
//                 } // k1 loop
//             } // colouring
//         } // face 4

//         if (domain->m_faceOffset[5] > -1) {
//             if (add_EM_S)
//                 fill(EM_S.begin(), EM_S.end(), zero);
//             if (add_EM_F) {
//                 EM_F[0] = zero;
//                 EM_F[1] = zero;
//                 EM_F[2] = zero;
//                 EM_F[3] = zero;
//             }

//             for (index_t k1_0=0; k1_0<2; k1_0++) { // colouring
// #pragma omp for
//                 for (index_t k1=k1_0; k1<NE1; k1+=2) {
//                     for (index_t k0=0; k0<NE0; ++k0) {
//                         const index_t e = domain->m_faceOffset[5]+INDEX2(k0,k1,NE0);
//                         ///////////////
//                         // process d //
//                         ///////////////
//                         if (add_EM_S) {
//                             const Scalar* d_p = d.getSampleDataRO(id, zero);
//                             EM_S[INDEX2(4,4,8)] = d_p[0]*w0;
//                             EM_S[INDEX2(5,4,8)] = d_p[0]*w0;
//                             EM_S[INDEX2(6,4,8)] = d_p[0]*w0;
//                             EM_S[INDEX2(7,4,8)] = d_p[0]*w0;
//                             EM_S[INDEX2(4,5,8)] = d_p[0]*w0;
//                             EM_S[INDEX2(5,5,8)] = d_p[0]*w0;
//                             EM_S[INDEX2(6,5,8)] = d_p[0]*w0;
//                             EM_S[INDEX2(7,5,8)] = d_p[0]*w0;
//                             EM_S[INDEX2(4,6,8)] = d_p[0]*w0;
//                             EM_S[INDEX2(5,6,8)] = d_p[0]*w0;
//                             EM_S[INDEX2(6,6,8)] = d_p[0]*w0;
//                             EM_S[INDEX2(7,6,8)] = d_p[0]*w0;
//                             EM_S[INDEX2(4,7,8)] = d_p[0]*w0;
//                             EM_S[INDEX2(5,7,8)] = d_p[0]*w0;
//                             EM_S[INDEX2(6,7,8)] = d_p[0]*w0;
//                             EM_S[INDEX2(7,7,8)] = d_p[0]*w0;
//                         }
//                         ///////////////
//                         // process y //
//                         ///////////////
//                         if (add_EM_F) {
//                             const Scalar* y_p = y.getSampleDataRO(id, zero);
//                             EM_F[4] = 4.*w0*y_p[0];
//                             EM_F[5] = 4.*w0*y_p[0];
//                             EM_F[6] = 4.*w0*y_p[0];
//                             EM_F[7] = 4.*w0*y_p[0];
//                         }
//                         const index_t firstNode=m_NN[0]*m_NN[1]*(m_NN[2]-2)+m_NN[0]*k1+k0;
//                         domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F,
//                                                add_EM_S, add_EM_F, firstNode);
//                     } // k0 loop
//                 } // k1 loop
//             } // colouring
//         } // face 5
//     } // end of parallel region
}

/****************************************************************************/
// PDE SYSTEM
/****************************************************************************/

template<class Scalar>
void DefaultAssembler3D<Scalar>::assemblePDESystem(AbstractSystemMatrix* mat,
                                       Data& rhs, const Data& A, const Data& B,
                                       const Data& C, const Data& D,
                                       const Data& X, const Data& Y) const
{
    // TODO
//     dim_t numEq, numComp;
//     if (!mat)
//         numEq=numComp=(rhs.isEmpty() ? 1 : rhs.getDataPointSize());
//     else {
//         numEq=mat->getRowBlockSize();
//         numComp=mat->getColumnBlockSize();
//     }

//     const double SQRT3 = 1.73205080756887719318;
//     const double w10 = -m_dx[0]/288;
//     const double w12 = w10*(-SQRT3 - 2);
//     const double w6 = w10*(SQRT3 - 2);
//     const double w18 = w10*(-4*SQRT3 - 7);
//     const double w4 = w10*(-4*SQRT3 + 7);
//     const double w11 = m_dx[1]/288;
//     const double w15 = w11*(SQRT3 + 2);
//     const double w5 = w11*(-SQRT3 + 2);
//     const double w2 = w11*(4*SQRT3 - 7);
//     const double w17 = w11*(4*SQRT3 + 7);
//     const double w8 = m_dx[2]/288;
//     const double w16 = w8*(SQRT3 + 2);
//     const double w1 = w8*(-SQRT3 + 2);
//     const double w20 = w8*(4*SQRT3 - 7);
//     const double w21 = w8*(-4*SQRT3 - 7);
//     const double w54 = -m_dx[0]*m_dx[1]/72;
//     const double w68 = -m_dx[0]*m_dx[1]/48;
//     const double w38 = w68*(-SQRT3 - 3)/36;
//     const double w45 = w68*(SQRT3 - 3)/36;
//     const double w35 = w68*(5*SQRT3 - 9)/36;
//     const double w46 = w68*(-5*SQRT3 - 9)/36;
//     const double w43 = w68*(-19*SQRT3 - 33)/36;
//     const double w44 = w68*(19*SQRT3 - 33)/36;
//     const double w66 = w68*(SQRT3 + 2);
//     const double w70 = w68*(-SQRT3 + 2);
//     const double w56 = -m_dx[0]*m_dx[2]/72;
//     const double w67 = -m_dx[0]*m_dx[2]/48;
//     const double w37 = w67*(-SQRT3 - 3)/36;
//     const double w40 = w67*(SQRT3 - 3)/36;
//     const double w34 = w67*(5*SQRT3 - 9)/36;
//     const double w42 = w67*(-5*SQRT3 - 9)/36;
//     const double w47 = w67*(19*SQRT3 + 33)/36;
//     const double w48 = w67*(-19*SQRT3 + 33)/36;
//     const double w65 = w67*(SQRT3 + 2);
//     const double w71 = w67*(-SQRT3 + 2);
//     const double w55 = -m_dx[1]*m_dx[2]/72;
//     const double w69 = -m_dx[1]*m_dx[2]/48;
//     const double w36 = w69*(SQRT3 - 3)/36;
//     const double w39 = w69*(-SQRT3 - 3)/36;
//     const double w33 = w69*(5*SQRT3 - 9)/36;
//     const double w41 = w69*(-5*SQRT3 - 9)/36;
//     const double w49 = w69*(19*SQRT3 - 33)/36;
//     const double w50 = w69*(-19*SQRT3 - 33)/36;
//     const double w64 = w69*(SQRT3 + 2);
//     const double w72 = w69*(-SQRT3 + 2);
//     const double w58 = m_dx[0]*m_dx[1]*m_dx[2]/1728;
//     const double w60 = w58*(-SQRT3 + 2);
//     const double w61 = w58*(SQRT3 + 2);
//     const double w57 = w58*(-4*SQRT3 + 7);
//     const double w59 = w58*(4*SQRT3 + 7);
//     const double w62 = w58*(15*SQRT3 + 26);
//     const double w63 = w58*(-15*SQRT3 + 26);
//     const double w75 = w58*6*(SQRT3 + 3);
//     const double w76 = w58*6*(-SQRT3 + 3);
//     const double w74 = w58*6*(5*SQRT3 + 9);
//     const double w77 = w58*6*(-5*SQRT3 + 9);
//     const double w13 = -m_dx[0]*m_dx[1]/(288*m_dx[2]);
//     const double w19 = w13*(4*SQRT3 + 7);
//     const double w7 = w13*(-4*SQRT3 + 7);
//     const double w23 = w13*(+SQRT3 - 2);
//     const double w25 = w13*(-SQRT3 - 2);
//     const double w22 = -m_dx[0]*m_dx[2]/(288*m_dx[1]);
//     const double w3 = w22*(SQRT3 - 2);
//     const double w9 = w22*(-SQRT3 - 2);
//     const double w24 = w22*(4*SQRT3 + 7);
//     const double w26 = w22*(-4*SQRT3 + 7);
//     const double w27 = -m_dx[1]*m_dx[2]/(288*m_dx[0]);
//     const double w0 = w27*(SQRT3 - 2);
//     const double w14 = w27*(-SQRT3 - 2);
//     const double w28 = w27*(-4*SQRT3 + 7);
//     const double w29 = w27*(4*SQRT3 + 7);
//     const dim_t NE0 = m_NE[0];
//     const dim_t NE1 = m_NE[1];
//     const dim_t NE2 = m_NE[2];
//     const bool add_EM_S = (!A.isEmpty() || !B.isEmpty() || !C.isEmpty() || !D.isEmpty());
//     const bool add_EM_F = (!X.isEmpty() || !Y.isEmpty());
//     const Scalar zero = static_cast<Scalar>(0);
//     rhs.requireWrite();

// #pragma omp parallel
//     {
//         vector<Scalar> EM_S(8*8*numEq*numComp, zero);
//         vector<Scalar> EM_F(8*numEq, zero);

//         for (index_t k2_0=0; k2_0<2; k2_0++) { // colouring
// #pragma omp for
//             for (index_t k2=k2_0; k2<NE2; k2+=2) {
//                 for (index_t k1=0; k1<NE1; ++k1) {
//                     for (index_t k0=0; k0<NE0; ++k0)  {
//                         const index_t e = k0 + NE0*k1 + NE0*NE1*k2;
//                         if (add_EM_S)
//                             fill(EM_S.begin(), EM_S.end(), zero);
//                         if (add_EM_F)
//                             fill(EM_F.begin(), EM_F.end(), zero);

//                         ///////////////
//                         // process A //
//                         ///////////////
//                         if (!A.isEmpty()) {
//                             const Scalar* A_p = A.getSampleDataRO(id, zero);
//                             if (A.actsExpanded()) {
//                                 for (index_t k=0; k<numEq; k++) {
//                                     for (index_t m=0; m<numComp; m++) {
//                                         const Scalar A_00_0 = A_p[INDEX5(k,0,m,0,0,numEq,3,numComp,3)];
//                                         const Scalar A_01_0 = A_p[INDEX5(k,0,m,1,0,numEq,3,numComp,3)];
//                                         const Scalar A_02_0 = A_p[INDEX5(k,0,m,2,0,numEq,3,numComp,3)];
//                                         const Scalar A_10_0 = A_p[INDEX5(k,1,m,0,0,numEq,3,numComp,3)];
//                                         const Scalar A_11_0 = A_p[INDEX5(k,1,m,1,0,numEq,3,numComp,3)];
//                                         const Scalar A_12_0 = A_p[INDEX5(k,1,m,2,0,numEq,3,numComp,3)];
//                                         const Scalar A_20_0 = A_p[INDEX5(k,2,m,0,0,numEq,3,numComp,3)];
//                                         const Scalar A_21_0 = A_p[INDEX5(k,2,m,1,0,numEq,3,numComp,3)];
//                                         const Scalar A_22_0 = A_p[INDEX5(k,2,m,2,0,numEq,3,numComp,3)];
//                                         const Scalar A_00_1 = A_p[INDEX5(k,0,m,0,1,numEq,3,numComp,3)];
//                                         const Scalar A_01_1 = A_p[INDEX5(k,0,m,1,1,numEq,3,numComp,3)];
//                                         const Scalar A_02_1 = A_p[INDEX5(k,0,m,2,1,numEq,3,numComp,3)];
//                                         const Scalar A_10_1 = A_p[INDEX5(k,1,m,0,1,numEq,3,numComp,3)];
//                                         const Scalar A_11_1 = A_p[INDEX5(k,1,m,1,1,numEq,3,numComp,3)];
//                                         const Scalar A_12_1 = A_p[INDEX5(k,1,m,2,1,numEq,3,numComp,3)];
//                                         const Scalar A_20_1 = A_p[INDEX5(k,2,m,0,1,numEq,3,numComp,3)];
//                                         const Scalar A_21_1 = A_p[INDEX5(k,2,m,1,1,numEq,3,numComp,3)];
//                                         const Scalar A_22_1 = A_p[INDEX5(k,2,m,2,1,numEq,3,numComp,3)];
//                                         const Scalar A_00_2 = A_p[INDEX5(k,0,m,0,2,numEq,3,numComp,3)];
//                                         const Scalar A_01_2 = A_p[INDEX5(k,0,m,1,2,numEq,3,numComp,3)];
//                                         const Scalar A_02_2 = A_p[INDEX5(k,0,m,2,2,numEq,3,numComp,3)];
//                                         const Scalar A_10_2 = A_p[INDEX5(k,1,m,0,2,numEq,3,numComp,3)];
//                                         const Scalar A_11_2 = A_p[INDEX5(k,1,m,1,2,numEq,3,numComp,3)];
//                                         const Scalar A_12_2 = A_p[INDEX5(k,1,m,2,2,numEq,3,numComp,3)];
//                                         const Scalar A_20_2 = A_p[INDEX5(k,2,m,0,2,numEq,3,numComp,3)];
//                                         const Scalar A_21_2 = A_p[INDEX5(k,2,m,1,2,numEq,3,numComp,3)];
//                                         const Scalar A_22_2 = A_p[INDEX5(k,2,m,2,2,numEq,3,numComp,3)];
//                                         const Scalar A_00_3 = A_p[INDEX5(k,0,m,0,3,numEq,3,numComp,3)];
//                                         const Scalar A_01_3 = A_p[INDEX5(k,0,m,1,3,numEq,3,numComp,3)];
//                                         const Scalar A_02_3 = A_p[INDEX5(k,0,m,2,3,numEq,3,numComp,3)];
//                                         const Scalar A_10_3 = A_p[INDEX5(k,1,m,0,3,numEq,3,numComp,3)];
//                                         const Scalar A_11_3 = A_p[INDEX5(k,1,m,1,3,numEq,3,numComp,3)];
//                                         const Scalar A_12_3 = A_p[INDEX5(k,1,m,2,3,numEq,3,numComp,3)];
//                                         const Scalar A_20_3 = A_p[INDEX5(k,2,m,0,3,numEq,3,numComp,3)];
//                                         const Scalar A_21_3 = A_p[INDEX5(k,2,m,1,3,numEq,3,numComp,3)];
//                                         const Scalar A_22_3 = A_p[INDEX5(k,2,m,2,3,numEq,3,numComp,3)];
//                                         const Scalar A_00_4 = A_p[INDEX5(k,0,m,0,4,numEq,3,numComp,3)];
//                                         const Scalar A_01_4 = A_p[INDEX5(k,0,m,1,4,numEq,3,numComp,3)];
//                                         const Scalar A_02_4 = A_p[INDEX5(k,0,m,2,4,numEq,3,numComp,3)];
//                                         const Scalar A_10_4 = A_p[INDEX5(k,1,m,0,4,numEq,3,numComp,3)];
//                                         const Scalar A_11_4 = A_p[INDEX5(k,1,m,1,4,numEq,3,numComp,3)];
//                                         const Scalar A_12_4 = A_p[INDEX5(k,1,m,2,4,numEq,3,numComp,3)];
//                                         const Scalar A_20_4 = A_p[INDEX5(k,2,m,0,4,numEq,3,numComp,3)];
//                                         const Scalar A_21_4 = A_p[INDEX5(k,2,m,1,4,numEq,3,numComp,3)];
//                                         const Scalar A_22_4 = A_p[INDEX5(k,2,m,2,4,numEq,3,numComp,3)];
//                                         const Scalar A_00_5 = A_p[INDEX5(k,0,m,0,5,numEq,3,numComp,3)];
//                                         const Scalar A_01_5 = A_p[INDEX5(k,0,m,1,5,numEq,3,numComp,3)];
//                                         const Scalar A_02_5 = A_p[INDEX5(k,0,m,2,5,numEq,3,numComp,3)];
//                                         const Scalar A_10_5 = A_p[INDEX5(k,1,m,0,5,numEq,3,numComp,3)];
//                                         const Scalar A_11_5 = A_p[INDEX5(k,1,m,1,5,numEq,3,numComp,3)];
//                                         const Scalar A_12_5 = A_p[INDEX5(k,1,m,2,5,numEq,3,numComp,3)];
//                                         const Scalar A_20_5 = A_p[INDEX5(k,2,m,0,5,numEq,3,numComp,3)];
//                                         const Scalar A_21_5 = A_p[INDEX5(k,2,m,1,5,numEq,3,numComp,3)];
//                                         const Scalar A_22_5 = A_p[INDEX5(k,2,m,2,5,numEq,3,numComp,3)];
//                                         const Scalar A_00_6 = A_p[INDEX5(k,0,m,0,6,numEq,3,numComp,3)];
//                                         const Scalar A_01_6 = A_p[INDEX5(k,0,m,1,6,numEq,3,numComp,3)];
//                                         const Scalar A_02_6 = A_p[INDEX5(k,0,m,2,6,numEq,3,numComp,3)];
//                                         const Scalar A_10_6 = A_p[INDEX5(k,1,m,0,6,numEq,3,numComp,3)];
//                                         const Scalar A_11_6 = A_p[INDEX5(k,1,m,1,6,numEq,3,numComp,3)];
//                                         const Scalar A_12_6 = A_p[INDEX5(k,1,m,2,6,numEq,3,numComp,3)];
//                                         const Scalar A_20_6 = A_p[INDEX5(k,2,m,0,6,numEq,3,numComp,3)];
//                                         const Scalar A_21_6 = A_p[INDEX5(k,2,m,1,6,numEq,3,numComp,3)];
//                                         const Scalar A_22_6 = A_p[INDEX5(k,2,m,2,6,numEq,3,numComp,3)];
//                                         const Scalar A_00_7 = A_p[INDEX5(k,0,m,0,7,numEq,3,numComp,3)];
//                                         const Scalar A_01_7 = A_p[INDEX5(k,0,m,1,7,numEq,3,numComp,3)];
//                                         const Scalar A_02_7 = A_p[INDEX5(k,0,m,2,7,numEq,3,numComp,3)];
//                                         const Scalar A_10_7 = A_p[INDEX5(k,1,m,0,7,numEq,3,numComp,3)];
//                                         const Scalar A_11_7 = A_p[INDEX5(k,1,m,1,7,numEq,3,numComp,3)];
//                                         const Scalar A_12_7 = A_p[INDEX5(k,1,m,2,7,numEq,3,numComp,3)];
//                                         const Scalar A_20_7 = A_p[INDEX5(k,2,m,0,7,numEq,3,numComp,3)];
//                                         const Scalar A_21_7 = A_p[INDEX5(k,2,m,1,7,numEq,3,numComp,3)];
//                                         const Scalar A_22_7 = A_p[INDEX5(k,2,m,2,7,numEq,3,numComp,3)];
//                                         const Scalar tmp0 = w18*(-A_12_7 + A_21_3);
//                                         const Scalar tmp1 = w13*(A_22_1 + A_22_2 + A_22_5 + A_22_6);
//                                         const Scalar tmp2 = w11*(-A_02_2 - A_02_5 + A_20_1 + A_20_6);
//                                         const Scalar tmp3 = w14*(A_00_2 + A_00_3 + A_00_6 + A_00_7);
//                                         const Scalar tmp4 = w7*(A_22_0 + A_22_4);
//                                         const Scalar tmp5 = w10*(A_12_1 + A_12_6 - A_21_2 - A_21_5);
//                                         const Scalar tmp6 = w3*(A_11_0 + A_11_2 + A_11_4 + A_11_6);
//                                         const Scalar tmp7 = w1*(A_01_0 + A_01_4 + A_10_0 + A_10_4);
//                                         const Scalar tmp8 = w4*(A_12_0 - A_21_4);
//                                         const Scalar tmp9 = w15*(-A_02_3 - A_02_6 + A_20_2 + A_20_7);
//                                         const Scalar tmp10 = w0*(A_00_0 + A_00_1 + A_00_4 + A_00_5);
//                                         const Scalar tmp11 = w16*(A_01_3 + A_01_7 + A_10_3 + A_10_7);
//                                         const Scalar tmp12 = w9*(A_11_1 + A_11_3 + A_11_5 + A_11_7);
//                                         const Scalar tmp13 = w12*(-A_12_3 - A_12_5 + A_21_1 + A_21_7);
//                                         const Scalar tmp14 = w5*(-A_02_1 - A_02_4 + A_20_0 + A_20_5);
//                                         const Scalar tmp15 = w8*(A_01_1 + A_01_2 + A_01_5 + A_01_6 + A_10_1 + A_10_2 + A_10_5 + A_10_6);
//                                         const Scalar tmp16 = w6*(-A_12_2 - A_12_4 + A_21_0 + A_21_6);
//                                         const Scalar tmp17 = w19*(A_22_3 + A_22_7);
//                                         const Scalar tmp18 = w17*(-A_02_7 + A_20_3);
//                                         const Scalar tmp19 = w2*(A_02_0 - A_20_4);
//                                         const Scalar tmp20 = w13*(-A_22_0 - A_22_1 - A_22_2 - A_22_3 - A_22_4 - A_22_5 - A_22_6 - A_22_7);
//                                         const Scalar tmp21 = w11*(-A_02_1 - A_02_3 - A_02_4 - A_02_6 + A_20_0 + A_20_2 + A_20_5 + A_20_7);
//                                         const Scalar tmp22 = w14*(-A_00_4 - A_00_5 - A_00_6 - A_00_7);
//                                         const Scalar tmp23 = w20*(A_01_2 + A_10_1);
//                                         const Scalar tmp24 = w10*(A_12_2 + A_12_3 + A_12_4 + A_12_5 - A_21_0 - A_21_1 - A_21_6 - A_21_7);
//                                         const Scalar tmp25 = w3*(-A_11_0 - A_11_1 - A_11_2 - A_11_3);
//                                         const Scalar tmp26 = w1*(-A_01_0 - A_01_3 - A_10_0 - A_10_3);
//                                         const Scalar tmp27 = w15*(-A_02_5 - A_02_7 + A_20_4 + A_20_6);
//                                         const Scalar tmp28 = w0*(-A_00_0 - A_00_1 - A_00_2 - A_00_3);
//                                         const Scalar tmp29 = w16*(-A_01_4 - A_01_7 - A_10_4 - A_10_7);
//                                         const Scalar tmp30 = w9*(-A_11_4 - A_11_5 - A_11_6 - A_11_7);
//                                         const Scalar tmp31 = w21*(A_01_5 + A_10_6);
//                                         const Scalar tmp32 = w12*(-A_12_6 - A_12_7 + A_21_4 + A_21_5);
//                                         const Scalar tmp33 = w5*(-A_02_0 - A_02_2 + A_20_1 + A_20_3);
//                                         const Scalar tmp34 = w8*(-A_01_1 - A_01_6 - A_10_2 - A_10_5);
//                                         const Scalar tmp35 = w6*(-A_12_0 - A_12_1 + A_21_2 + A_21_3);
//                                         const Scalar tmp36 = w20*(-A_01_6 + A_10_4);
//                                         const Scalar tmp37 = w18*(A_12_3 - A_21_1);
//                                         const Scalar tmp38 = w11*(-A_02_0 - A_02_2 - A_02_5 - A_02_7 - A_20_0 - A_20_2 - A_20_5 - A_20_7);
//                                         const Scalar tmp39 = w14*(A_00_0 + A_00_1 + A_00_2 + A_00_3);
//                                         const Scalar tmp40 = w26*(A_11_4 + A_11_6);
//                                         const Scalar tmp41 = w0*(A_00_4 + A_00_5 + A_00_6 + A_00_7);
//                                         const Scalar tmp42 = w10*(-A_12_2 - A_12_5 + A_21_0 + A_21_7);
//                                         const Scalar tmp43 = w22*(A_11_0 + A_11_2 + A_11_5 + A_11_7);
//                                         const Scalar tmp44 = w1*(A_01_4 + A_01_7 - A_10_5 - A_10_6);
//                                         const Scalar tmp45 = w25*(A_22_1 + A_22_3 + A_22_5 + A_22_7);
//                                         const Scalar tmp46 = w4*(-A_12_4 + A_21_6);
//                                         const Scalar tmp47 = w15*(-A_02_1 - A_02_3 - A_20_1 - A_20_3);
//                                         const Scalar tmp48 = w21*(-A_01_1 + A_10_3);
//                                         const Scalar tmp49 = w16*(A_01_0 + A_01_3 - A_10_1 - A_10_2);
//                                         const Scalar tmp50 = w5*(-A_02_4 - A_02_6 - A_20_4 - A_20_6);
//                                         const Scalar tmp51 = w12*(A_12_1 + A_12_7 - A_21_3 - A_21_5);
//                                         const Scalar tmp52 = w24*(A_11_1 + A_11_3);
//                                         const Scalar tmp53 = w8*(A_01_2 + A_01_5 - A_10_0 - A_10_7);
//                                         const Scalar tmp54 = w6*(A_12_0 + A_12_6 - A_21_2 - A_21_4);
//                                         const Scalar tmp55 = w23*(A_22_0 + A_22_2 + A_22_4 + A_22_6);
//                                         const Scalar tmp56 = w18*(A_12_4 - A_21_6);
//                                         const Scalar tmp57 = w14*(A_00_4 + A_00_5 + A_00_6 + A_00_7);
//                                         const Scalar tmp58 = w26*(A_11_1 + A_11_3);
//                                         const Scalar tmp59 = w20*(-A_01_1 + A_10_3);
//                                         const Scalar tmp60 = w1*(A_01_0 + A_01_3 - A_10_1 - A_10_2);
//                                         const Scalar tmp61 = w25*(A_22_0 + A_22_2 + A_22_4 + A_22_6);
//                                         const Scalar tmp62 = w4*(-A_12_3 + A_21_1);
//                                         const Scalar tmp63 = w15*(-A_02_4 - A_02_6 - A_20_4 - A_20_6);
//                                         const Scalar tmp64 = w0*(A_00_0 + A_00_1 + A_00_2 + A_00_3);
//                                         const Scalar tmp65 = w16*(A_01_4 + A_01_7 - A_10_5 - A_10_6);
//                                         const Scalar tmp66 = w24*(A_11_4 + A_11_6);
//                                         const Scalar tmp67 = w21*(-A_01_6 + A_10_4);
//                                         const Scalar tmp68 = w12*(A_12_0 + A_12_6 - A_21_2 - A_21_4);
//                                         const Scalar tmp69 = w5*(-A_02_1 - A_02_3 - A_20_1 - A_20_3);
//                                         const Scalar tmp70 = w6*(A_12_1 + A_12_7 - A_21_3 - A_21_5);
//                                         const Scalar tmp71 = w23*(A_22_1 + A_22_3 + A_22_5 + A_22_7);
//                                         const Scalar tmp72 = w20*(A_01_5 + A_10_6);
//                                         const Scalar tmp73 = w14*(-A_00_0 - A_00_1 - A_00_2 - A_00_3);
//                                         const Scalar tmp74 = w0*(-A_00_4 - A_00_5 - A_00_6 - A_00_7);
//                                         const Scalar tmp75 = w3*(-A_11_4 - A_11_5 - A_11_6 - A_11_7);
//                                         const Scalar tmp76 = w1*(-A_01_4 - A_01_7 - A_10_4 - A_10_7);
//                                         const Scalar tmp77 = w15*(-A_02_0 - A_02_2 + A_20_1 + A_20_3);
//                                         const Scalar tmp78 = w21*(A_01_2 + A_10_1);
//                                         const Scalar tmp79 = w16*(-A_01_0 - A_01_3 - A_10_0 - A_10_3);
//                                         const Scalar tmp80 = w9*(-A_11_0 - A_11_1 - A_11_2 - A_11_3);
//                                         const Scalar tmp81 = w12*(-A_12_0 - A_12_1 + A_21_2 + A_21_3);
//                                         const Scalar tmp82 = w5*(-A_02_5 - A_02_7 + A_20_4 + A_20_6);
//                                         const Scalar tmp83 = w6*(-A_12_6 - A_12_7 + A_21_4 + A_21_5);
//                                         const Scalar tmp84 = w6*(-A_12_2 - A_12_3 - A_21_2 - A_21_3);
//                                         const Scalar tmp85 = w11*(A_02_1 + A_02_6 - A_20_0 - A_20_7);
//                                         const Scalar tmp86 = w20*(A_01_3 - A_10_2);
//                                         const Scalar tmp87 = w10*(A_12_0 + A_12_1 + A_12_6 + A_12_7 + A_21_0 + A_21_1 + A_21_6 + A_21_7);
//                                         const Scalar tmp88 = w3*(A_11_0 + A_11_1 + A_11_2 + A_11_3);
//                                         const Scalar tmp89 = w23*(A_22_2 + A_22_3 + A_22_6 + A_22_7);
//                                         const Scalar tmp90 = w1*(-A_01_1 - A_01_2 + A_10_0 + A_10_3);
//                                         const Scalar tmp91 = w25*(A_22_0 + A_22_1 + A_22_4 + A_22_5);
//                                         const Scalar tmp92 = w15*(A_02_0 + A_02_5 - A_20_1 - A_20_4);
//                                         const Scalar tmp93 = w21*(A_01_4 - A_10_5);
//                                         const Scalar tmp94 = w16*(-A_01_5 - A_01_6 + A_10_4 + A_10_7);
//                                         const Scalar tmp95 = w28*(A_00_2 + A_00_3);
//                                         const Scalar tmp96 = w12*(-A_12_4 - A_12_5 - A_21_4 - A_21_5);
//                                         const Scalar tmp97 = w29*(A_00_4 + A_00_5);
//                                         const Scalar tmp98 = w5*(A_02_2 + A_02_7 - A_20_3 - A_20_6);
//                                         const Scalar tmp99 = w8*(-A_01_0 - A_01_7 + A_10_1 + A_10_6);
//                                         const Scalar tmp100 = w9*(A_11_4 + A_11_5 + A_11_6 + A_11_7);
//                                         const Scalar tmp101 = w27*(A_00_0 + A_00_1 + A_00_6 + A_00_7);
//                                         const Scalar tmp102 = w17*(A_02_4 - A_20_5);
//                                         const Scalar tmp103 = w2*(-A_02_3 + A_20_2);
//                                         const Scalar tmp104 = w13*(A_22_0 + A_22_1 + A_22_2 + A_22_3 + A_22_4 + A_22_5 + A_22_6 + A_22_7);
//                                         const Scalar tmp105 = w6*(-A_12_4 - A_12_5 - A_21_2 - A_21_3);
//                                         const Scalar tmp106 = w22*(A_11_0 + A_11_1 + A_11_2 + A_11_3 + A_11_4 + A_11_5 + A_11_6 + A_11_7);
//                                         const Scalar tmp107 = w1*(-A_01_2 - A_01_6 - A_10_1 - A_10_5);
//                                         const Scalar tmp108 = w15*(-A_02_1 - A_02_3 - A_20_4 - A_20_6);
//                                         const Scalar tmp109 = w16*(-A_01_1 - A_01_5 - A_10_2 - A_10_6);
//                                         const Scalar tmp110 = w12*(-A_12_2 - A_12_3 - A_21_4 - A_21_5);
//                                         const Scalar tmp111 = w5*(-A_02_4 - A_02_6 - A_20_1 - A_20_3);
//                                         const Scalar tmp112 = w8*(-A_01_0 - A_01_3 - A_01_4 - A_01_7 - A_10_0 - A_10_3 - A_10_4 - A_10_7);
//                                         const Scalar tmp113 = w27*(A_00_0 + A_00_1 + A_00_2 + A_00_3 + A_00_4 + A_00_5 + A_00_6 + A_00_7);
//                                         const Scalar tmp114 = w11*(A_02_0 + A_02_2 + A_02_5 + A_02_7 - A_20_1 - A_20_3 - A_20_4 - A_20_6);
//                                         const Scalar tmp115 = w21*(-A_01_4 - A_10_7);
//                                         const Scalar tmp116 = w20*(-A_01_3 - A_10_0);
//                                         const Scalar tmp117 = w15*(A_02_4 + A_02_6 - A_20_5 - A_20_7);
//                                         const Scalar tmp118 = w16*(A_01_5 + A_01_6 + A_10_5 + A_10_6);
//                                         const Scalar tmp119 = w5*(A_02_1 + A_02_3 - A_20_0 - A_20_2);
//                                         const Scalar tmp120 = w8*(A_01_0 + A_01_7 + A_10_3 + A_10_4);
//                                         const Scalar tmp121 = w1*(A_01_1 + A_01_2 + A_10_1 + A_10_2);
//                                         const Scalar tmp122 = w18*(A_12_2 - A_21_6);
//                                         const Scalar tmp123 = w13*(A_22_0 + A_22_3 + A_22_4 + A_22_7);
//                                         const Scalar tmp124 = w11*(-A_02_0 - A_02_7 + A_20_3 + A_20_4);
//                                         const Scalar tmp125 = w7*(A_22_1 + A_22_5);
//                                         const Scalar tmp126 = w10*(-A_12_3 - A_12_4 + A_21_0 + A_21_7);
//                                         const Scalar tmp127 = w3*(A_11_1 + A_11_3 + A_11_5 + A_11_7);
//                                         const Scalar tmp128 = w1*(-A_01_1 - A_01_5 - A_10_1 - A_10_5);
//                                         const Scalar tmp129 = w4*(-A_12_5 + A_21_1);
//                                         const Scalar tmp130 = w16*(-A_01_2 - A_01_6 - A_10_2 - A_10_6);
//                                         const Scalar tmp131 = w9*(A_11_0 + A_11_2 + A_11_4 + A_11_6);
//                                         const Scalar tmp132 = w19*(A_22_2 + A_22_6);
//                                         const Scalar tmp133 = w17*(-A_02_2 + A_20_6);
//                                         const Scalar tmp134 = w2*(A_02_5 - A_20_1);
//                                         const Scalar tmp135 = w11*(A_02_1 + A_02_3 + A_02_4 + A_02_6 + A_20_1 + A_20_3 + A_20_4 + A_20_6);
//                                         const Scalar tmp136 = w1*(A_01_3 + A_01_7 + A_10_0 + A_10_4);
//                                         const Scalar tmp137 = w15*(A_02_0 + A_02_2 + A_20_5 + A_20_7);
//                                         const Scalar tmp138 = w16*(A_01_0 + A_01_4 + A_10_3 + A_10_7);
//                                         const Scalar tmp139 = w5*(A_02_5 + A_02_7 + A_20_0 + A_20_2);
//                                         const Scalar tmp140 = w18*(A_12_5 - A_21_1);
//                                         const Scalar tmp141 = w14*(A_00_0 + A_00_1 + A_00_4 + A_00_5);
//                                         const Scalar tmp142 = w7*(A_22_2 + A_22_6);
//                                         const Scalar tmp143 = w1*(-A_01_2 - A_01_6 - A_10_2 - A_10_6);
//                                         const Scalar tmp144 = w4*(-A_12_2 + A_21_6);
//                                         const Scalar tmp145 = w15*(-A_02_1 - A_02_4 + A_20_0 + A_20_5);
//                                         const Scalar tmp146 = w0*(A_00_2 + A_00_3 + A_00_6 + A_00_7);
//                                         const Scalar tmp147 = w16*(-A_01_1 - A_01_5 - A_10_1 - A_10_5);
//                                         const Scalar tmp148 = w5*(-A_02_3 - A_02_6 + A_20_2 + A_20_7);
//                                         const Scalar tmp149 = w19*(A_22_1 + A_22_5);
//                                         const Scalar tmp150 = w17*(-A_02_5 + A_20_1);
//                                         const Scalar tmp151 = w2*(A_02_2 - A_20_6);
//                                         const Scalar tmp152 = w18*(A_12_3 - A_21_7);
//                                         const Scalar tmp153 = w11*(A_02_1 + A_02_6 - A_20_2 - A_20_5);
//                                         const Scalar tmp154 = w10*(-A_12_2 - A_12_5 + A_21_1 + A_21_6);
//                                         const Scalar tmp155 = w4*(-A_12_4 + A_21_0);
//                                         const Scalar tmp156 = w15*(A_02_2 + A_02_7 - A_20_3 - A_20_6);
//                                         const Scalar tmp157 = w5*(A_02_0 + A_02_5 - A_20_1 - A_20_4);
//                                         const Scalar tmp158 = w17*(A_02_3 - A_20_7);
//                                         const Scalar tmp159 = w2*(-A_02_4 + A_20_0);
//                                         const Scalar tmp160 = w6*(A_12_6 + A_12_7 + A_21_0 + A_21_1);
//                                         const Scalar tmp161 = w10*(-A_12_2 - A_12_3 - A_12_4 - A_12_5 - A_21_2 - A_21_3 - A_21_4 - A_21_5);
//                                         const Scalar tmp162 = w1*(A_01_0 + A_01_4 + A_10_3 + A_10_7);
//                                         const Scalar tmp163 = w16*(A_01_3 + A_01_7 + A_10_0 + A_10_4);
//                                         const Scalar tmp164 = w12*(A_12_0 + A_12_1 + A_21_6 + A_21_7);
//                                         const Scalar tmp165 = w20*(A_01_6 + A_10_5);
//                                         const Scalar tmp166 = w10*(-A_12_0 - A_12_1 - A_12_6 - A_12_7 + A_21_2 + A_21_3 + A_21_4 + A_21_5);
//                                         const Scalar tmp167 = w15*(A_02_1 + A_02_3 - A_20_0 - A_20_2);
//                                         const Scalar tmp168 = w21*(A_01_1 + A_10_2);
//                                         const Scalar tmp169 = w12*(A_12_2 + A_12_3 - A_21_0 - A_21_1);
//                                         const Scalar tmp170 = w5*(A_02_4 + A_02_6 - A_20_5 - A_20_7);
//                                         const Scalar tmp171 = w8*(-A_01_2 - A_01_5 - A_10_1 - A_10_6);
//                                         const Scalar tmp172 = w6*(A_12_4 + A_12_5 - A_21_6 - A_21_7);
//                                         const Scalar tmp173 = w2*(A_02_1 + A_20_4);
//                                         const Scalar tmp174 = w11*(-A_02_3 - A_02_4 - A_20_1 - A_20_6);
//                                         const Scalar tmp175 = w14*(-A_00_2 - A_00_3 - A_00_6 - A_00_7);
//                                         const Scalar tmp176 = w22*(-A_11_0 - A_11_1 - A_11_2 - A_11_3 - A_11_4 - A_11_5 - A_11_6 - A_11_7);
//                                         const Scalar tmp177 = w1*(A_01_1 + A_01_5 - A_10_0 - A_10_4);
//                                         const Scalar tmp178 = w25*(-A_22_2 - A_22_3 - A_22_6 - A_22_7);
//                                         const Scalar tmp179 = w15*(-A_02_2 - A_02_7 - A_20_2 - A_20_7);
//                                         const Scalar tmp180 = w0*(-A_00_0 - A_00_1 - A_00_4 - A_00_5);
//                                         const Scalar tmp181 = w16*(A_01_2 + A_01_6 - A_10_3 - A_10_7);
//                                         const Scalar tmp182 = w12*(-A_12_6 - A_12_7 + A_21_2 + A_21_3);
//                                         const Scalar tmp183 = w5*(-A_02_0 - A_02_5 - A_20_0 - A_20_5);
//                                         const Scalar tmp184 = w8*(A_01_0 + A_01_3 + A_01_4 + A_01_7 - A_10_1 - A_10_2 - A_10_5 - A_10_6);
//                                         const Scalar tmp185 = w6*(-A_12_0 - A_12_1 + A_21_4 + A_21_5);
//                                         const Scalar tmp186 = w17*(-A_02_6 - A_20_3);
//                                         const Scalar tmp187 = w23*(-A_22_0 - A_22_1 - A_22_4 - A_22_5);
//                                         const Scalar tmp188 = w18*(A_12_4 - A_21_0);
//                                         const Scalar tmp189 = w7*(A_22_3 + A_22_7);
//                                         const Scalar tmp190 = w1*(A_01_3 + A_01_7 + A_10_3 + A_10_7);
//                                         const Scalar tmp191 = w4*(-A_12_3 + A_21_7);
//                                         const Scalar tmp192 = w16*(A_01_0 + A_01_4 + A_10_0 + A_10_4);
//                                         const Scalar tmp193 = w19*(A_22_0 + A_22_4);
//                                         const Scalar tmp194 = w17*(A_02_4 - A_20_0);
//                                         const Scalar tmp195 = w2*(-A_02_3 + A_20_7);
//                                         const Scalar tmp196 = w20*(-A_01_7 - A_10_4);
//                                         const Scalar tmp197 = w21*(-A_01_0 - A_10_3);
//                                         const Scalar tmp198 = w16*(A_01_1 + A_01_2 + A_10_1 + A_10_2);
//                                         const Scalar tmp199 = w8*(A_01_3 + A_01_4 + A_10_0 + A_10_7);
//                                         const Scalar tmp200 = w1*(A_01_5 + A_01_6 + A_10_5 + A_10_6);
//                                         const Scalar tmp201 = w27*(A_00_2 + A_00_3 + A_00_4 + A_00_5);
//                                         const Scalar tmp202 = w11*(-A_02_2 - A_02_5 + A_20_3 + A_20_4);
//                                         const Scalar tmp203 = w20*(A_01_0 - A_10_1);
//                                         const Scalar tmp204 = w23*(A_22_0 + A_22_1 + A_22_4 + A_22_5);
//                                         const Scalar tmp205 = w25*(A_22_2 + A_22_3 + A_22_6 + A_22_7);
//                                         const Scalar tmp206 = w21*(A_01_7 - A_10_6);
//                                         const Scalar tmp207 = w12*(A_12_6 + A_12_7 + A_21_6 + A_21_7);
//                                         const Scalar tmp208 = w28*(A_00_0 + A_00_1);
//                                         const Scalar tmp209 = w29*(A_00_6 + A_00_7);
//                                         const Scalar tmp210 = w8*(-A_01_3 - A_01_4 + A_10_2 + A_10_5);
//                                         const Scalar tmp211 = w6*(A_12_0 + A_12_1 + A_21_0 + A_21_1);
//                                         const Scalar tmp212 = w17*(-A_02_7 + A_20_6);
//                                         const Scalar tmp213 = w2*(A_02_0 - A_20_1);
//                                         const Scalar tmp214 = w13*(-A_22_1 - A_22_2 - A_22_5 - A_22_6);
//                                         const Scalar tmp215 = w22*(-A_11_0 - A_11_2 - A_11_5 - A_11_7);
//                                         const Scalar tmp216 = w8*(A_01_0 + A_01_7 + A_10_0 + A_10_7);
//                                         const Scalar tmp217 = w27*(-A_00_0 - A_00_1 - A_00_6 - A_00_7);
//                                         const Scalar tmp218 = w17*(-A_02_3 - A_20_3);
//                                         const Scalar tmp219 = w2*(A_02_4 + A_20_4);
//                                         const Scalar tmp220 = w11*(-A_02_1 - A_02_6 - A_20_1 - A_20_6);
//                                         const Scalar tmp221 = w26*(-A_11_4 - A_11_6);
//                                         const Scalar tmp222 = w10*(A_12_2 + A_12_5 + A_21_2 + A_21_5);
//                                         const Scalar tmp223 = w20*(-A_01_4 - A_10_4);
//                                         const Scalar tmp224 = w21*(-A_01_3 - A_10_3);
//                                         const Scalar tmp225 = w6*(-A_12_0 - A_12_6 - A_21_0 - A_21_6);
//                                         const Scalar tmp226 = w7*(-A_22_0 - A_22_4);
//                                         const Scalar tmp227 = w24*(-A_11_1 - A_11_3);
//                                         const Scalar tmp228 = w19*(-A_22_3 - A_22_7);
//                                         const Scalar tmp229 = w18*(-A_12_3 - A_21_3);
//                                         const Scalar tmp230 = w4*(A_12_4 + A_21_4);
//                                         const Scalar tmp231 = w28*(-A_00_4 - A_00_5);
//                                         const Scalar tmp232 = w12*(-A_12_1 - A_12_7 - A_21_1 - A_21_7);
//                                         const Scalar tmp233 = w29*(-A_00_2 - A_00_3);
//                                         const Scalar tmp234 = w20*(-A_01_5 + A_10_7);
//                                         const Scalar tmp235 = w18*(-A_12_0 + A_21_2);
//                                         const Scalar tmp236 = w26*(A_11_5 + A_11_7);
//                                         const Scalar tmp237 = w10*(A_12_1 + A_12_6 - A_21_3 - A_21_4);
//                                         const Scalar tmp238 = w22*(A_11_1 + A_11_3 + A_11_4 + A_11_6);
//                                         const Scalar tmp239 = w4*(A_12_7 - A_21_5);
//                                         const Scalar tmp240 = w15*(A_02_0 + A_02_2 + A_20_0 + A_20_2);
//                                         const Scalar tmp241 = w21*(-A_01_2 + A_10_0);
//                                         const Scalar tmp242 = w5*(A_02_5 + A_02_7 + A_20_5 + A_20_7);
//                                         const Scalar tmp243 = w12*(-A_12_2 - A_12_4 + A_21_0 + A_21_6);
//                                         const Scalar tmp244 = w24*(A_11_0 + A_11_2);
//                                         const Scalar tmp245 = w8*(A_01_1 + A_01_6 - A_10_3 - A_10_4);
//                                         const Scalar tmp246 = w6*(-A_12_3 - A_12_5 + A_21_1 + A_21_7);
//                                         const Scalar tmp247 = w11*(A_02_3 + A_02_4 - A_20_2 - A_20_5);
//                                         const Scalar tmp248 = w20*(-A_01_1 + A_10_0);
//                                         const Scalar tmp249 = w21*(-A_01_6 + A_10_7);
//                                         const Scalar tmp250 = w8*(A_01_2 + A_01_5 - A_10_3 - A_10_4);
//                                         const Scalar tmp251 = w17*(A_02_6 - A_20_7);
//                                         const Scalar tmp252 = w2*(-A_02_1 + A_20_0);
//                                         const Scalar tmp253 = w17*(-A_02_4 - A_20_4);
//                                         const Scalar tmp254 = w2*(A_02_3 + A_20_3);
//                                         const Scalar tmp255 = w26*(-A_11_1 - A_11_3);
//                                         const Scalar tmp256 = w20*(-A_01_3 - A_10_3);
//                                         const Scalar tmp257 = w21*(-A_01_4 - A_10_4);
//                                         const Scalar tmp258 = w6*(-A_12_1 - A_12_7 - A_21_1 - A_21_7);
//                                         const Scalar tmp259 = w7*(-A_22_3 - A_22_7);
//                                         const Scalar tmp260 = w15*(-A_02_0 - A_02_5 - A_20_0 - A_20_5);
//                                         const Scalar tmp261 = w24*(-A_11_4 - A_11_6);
//                                         const Scalar tmp262 = w19*(-A_22_0 - A_22_4);
//                                         const Scalar tmp263 = w18*(-A_12_4 - A_21_4);
//                                         const Scalar tmp264 = w4*(A_12_3 + A_21_3);
//                                         const Scalar tmp265 = w28*(-A_00_2 - A_00_3);
//                                         const Scalar tmp266 = w12*(-A_12_0 - A_12_6 - A_21_0 - A_21_6);
//                                         const Scalar tmp267 = w5*(-A_02_2 - A_02_7 - A_20_2 - A_20_7);
//                                         const Scalar tmp268 = w29*(-A_00_4 - A_00_5);
//                                         const Scalar tmp269 = w11*(A_02_2 + A_02_5 + A_20_0 + A_20_7);
//                                         const Scalar tmp270 = w1*(-A_01_0 - A_01_4 + A_10_1 + A_10_5);
//                                         const Scalar tmp271 = w15*(A_02_3 + A_02_6 + A_20_3 + A_20_6);
//                                         const Scalar tmp272 = w16*(-A_01_3 - A_01_7 + A_10_2 + A_10_6);
//                                         const Scalar tmp273 = w5*(A_02_1 + A_02_4 + A_20_1 + A_20_4);
//                                         const Scalar tmp274 = w8*(-A_01_1 - A_01_2 - A_01_5 - A_01_6 + A_10_0 + A_10_3 + A_10_4 + A_10_7);
//                                         const Scalar tmp275 = w17*(A_02_7 + A_20_2);
//                                         const Scalar tmp276 = w2*(-A_02_0 - A_20_5);
//                                         const Scalar tmp277 = w18*(-A_12_1 + A_21_5);
//                                         const Scalar tmp278 = w11*(A_02_3 + A_02_4 - A_20_0 - A_20_7);
//                                         const Scalar tmp279 = w10*(A_12_0 + A_12_7 - A_21_3 - A_21_4);
//                                         const Scalar tmp280 = w4*(A_12_6 - A_21_2);
//                                         const Scalar tmp281 = w17*(A_02_1 - A_20_5);
//                                         const Scalar tmp282 = w2*(-A_02_6 + A_20_2);
//                                         const Scalar tmp283 = w11*(A_02_0 + A_02_7 + A_20_2 + A_20_5);
//                                         const Scalar tmp284 = w12*(A_12_2 + A_12_3 - A_21_6 - A_21_7);
//                                         const Scalar tmp285 = w6*(A_12_4 + A_12_5 - A_21_0 - A_21_1);
//                                         const Scalar tmp286 = w17*(A_02_2 + A_20_7);
//                                         const Scalar tmp287 = w2*(-A_02_5 - A_20_0);
//                                         const Scalar tmp288 = w13*(-A_22_0 - A_22_3 - A_22_4 - A_22_7);
//                                         const Scalar tmp289 = w22*(-A_11_1 - A_11_3 - A_11_4 - A_11_6);
//                                         const Scalar tmp290 = w8*(-A_01_1 - A_01_6 - A_10_1 - A_10_6);
//                                         const Scalar tmp291 = w17*(A_02_2 + A_20_2);
//                                         const Scalar tmp292 = w2*(-A_02_5 - A_20_5);
//                                         const Scalar tmp293 = w11*(A_02_0 + A_02_7 + A_20_0 + A_20_7);
//                                         const Scalar tmp294 = w26*(-A_11_5 - A_11_7);
//                                         const Scalar tmp295 = w10*(A_12_3 + A_12_4 + A_21_3 + A_21_4);
//                                         const Scalar tmp296 = w20*(A_01_5 + A_10_5);
//                                         const Scalar tmp297 = w21*(A_01_2 + A_10_2);
//                                         const Scalar tmp298 = w7*(-A_22_1 - A_22_5);
//                                         const Scalar tmp299 = w24*(-A_11_0 - A_11_2);
//                                         const Scalar tmp300 = w19*(-A_22_2 - A_22_6);
//                                         const Scalar tmp301 = w18*(-A_12_2 - A_21_2);
//                                         const Scalar tmp302 = w4*(A_12_5 + A_21_5);
//                                         const Scalar tmp303 = w8*(A_01_3 + A_01_4 + A_10_3 + A_10_4);
//                                         const Scalar tmp304 = w27*(-A_00_2 - A_00_3 - A_00_4 - A_00_5);
//                                         const Scalar tmp305 = w17*(A_02_7 + A_20_7);
//                                         const Scalar tmp306 = w2*(-A_02_0 - A_20_0);
//                                         const Scalar tmp307 = w11*(A_02_2 + A_02_5 + A_20_2 + A_20_5);
//                                         const Scalar tmp308 = w26*(-A_11_0 - A_11_2);
//                                         const Scalar tmp309 = w10*(-A_12_1 - A_12_6 - A_21_1 - A_21_6);
//                                         const Scalar tmp310 = w20*(-A_01_0 - A_10_0);
//                                         const Scalar tmp311 = w21*(-A_01_7 - A_10_7);
//                                         const Scalar tmp312 = w6*(A_12_2 + A_12_4 + A_21_2 + A_21_4);
//                                         const Scalar tmp313 = w24*(-A_11_5 - A_11_7);
//                                         const Scalar tmp314 = w18*(A_12_7 + A_21_7);
//                                         const Scalar tmp315 = w4*(-A_12_0 - A_21_0);
//                                         const Scalar tmp316 = w28*(-A_00_0 - A_00_1);
//                                         const Scalar tmp317 = w12*(A_12_3 + A_12_5 + A_21_3 + A_21_5);
//                                         const Scalar tmp318 = w29*(-A_00_6 - A_00_7);
//                                         const Scalar tmp319 = w18*(-A_12_7 + A_21_5);
//                                         const Scalar tmp320 = w26*(A_11_0 + A_11_2);
//                                         const Scalar tmp321 = w21*(-A_01_5 + A_10_7);
//                                         const Scalar tmp322 = w20*(-A_01_2 + A_10_0);
//                                         const Scalar tmp323 = w4*(A_12_0 - A_21_2);
//                                         const Scalar tmp324 = w15*(A_02_5 + A_02_7 + A_20_5 + A_20_7);
//                                         const Scalar tmp325 = w24*(A_11_5 + A_11_7);
//                                         const Scalar tmp326 = w5*(A_02_0 + A_02_2 + A_20_0 + A_20_2);
//                                         const Scalar tmp327 = w18*(A_12_7 + A_21_1);
//                                         const Scalar tmp328 = w10*(-A_12_1 - A_12_6 - A_21_0 - A_21_7);
//                                         const Scalar tmp329 = w3*(-A_11_0 - A_11_2 - A_11_4 - A_11_6);
//                                         const Scalar tmp330 = w1*(A_01_2 + A_01_6 - A_10_0 - A_10_4);
//                                         const Scalar tmp331 = w4*(-A_12_0 - A_21_6);
//                                         const Scalar tmp332 = w25*(-A_22_1 - A_22_3 - A_22_5 - A_22_7);
//                                         const Scalar tmp333 = w15*(-A_02_5 - A_02_7 + A_20_1 + A_20_3);
//                                         const Scalar tmp334 = w16*(A_01_1 + A_01_5 - A_10_3 - A_10_7);
//                                         const Scalar tmp335 = w9*(-A_11_1 - A_11_3 - A_11_5 - A_11_7);
//                                         const Scalar tmp336 = w5*(-A_02_0 - A_02_2 + A_20_4 + A_20_6);
//                                         const Scalar tmp337 = w27*(-A_00_0 - A_00_1 - A_00_2 - A_00_3 - A_00_4 - A_00_5 - A_00_6 - A_00_7);
//                                         const Scalar tmp338 = w23*(-A_22_0 - A_22_2 - A_22_4 - A_22_6);
//                                         const Scalar tmp339 = w14*(-A_00_0 - A_00_1 - A_00_4 - A_00_5);
//                                         const Scalar tmp340 = w23*(-A_22_2 - A_22_3 - A_22_6 - A_22_7);
//                                         const Scalar tmp341 = w1*(A_01_2 + A_01_6 - A_10_3 - A_10_7);
//                                         const Scalar tmp342 = w25*(-A_22_0 - A_22_1 - A_22_4 - A_22_5);
//                                         const Scalar tmp343 = w15*(A_02_1 + A_02_4 + A_20_1 + A_20_4);
//                                         const Scalar tmp344 = w0*(-A_00_2 - A_00_3 - A_00_6 - A_00_7);
//                                         const Scalar tmp345 = w16*(A_01_1 + A_01_5 - A_10_0 - A_10_4);
//                                         const Scalar tmp346 = w12*(A_12_4 + A_12_5 - A_21_0 - A_21_1);
//                                         const Scalar tmp347 = w5*(A_02_3 + A_02_6 + A_20_3 + A_20_6);
//                                         const Scalar tmp348 = w6*(A_12_2 + A_12_3 - A_21_6 - A_21_7);
//                                         const Scalar tmp349 = w17*(A_02_5 + A_20_0);
//                                         const Scalar tmp350 = w2*(-A_02_2 - A_20_7);
//                                         const Scalar tmp351 = w8*(-A_01_2 - A_01_5 - A_10_2 - A_10_5);
//                                         const Scalar tmp352 = w17*(-A_02_1 - A_20_1);
//                                         const Scalar tmp353 = w2*(A_02_6 + A_20_6);
//                                         const Scalar tmp354 = w11*(-A_02_3 - A_02_4 - A_20_3 - A_20_4);
//                                         const Scalar tmp355 = w10*(-A_12_0 - A_12_7 - A_21_0 - A_21_7);
//                                         const Scalar tmp356 = w20*(A_01_6 + A_10_6);
//                                         const Scalar tmp357 = w21*(A_01_1 + A_10_1);
//                                         const Scalar tmp358 = w7*(-A_22_2 - A_22_6);
//                                         const Scalar tmp359 = w19*(-A_22_1 - A_22_5);
//                                         const Scalar tmp360 = w18*(A_12_1 + A_21_1);
//                                         const Scalar tmp361 = w4*(-A_12_6 - A_21_6);
//                                         const Scalar tmp362 = w28*(-A_00_6 - A_00_7);
//                                         const Scalar tmp363 = w29*(-A_00_0 - A_00_1);
//                                         const Scalar tmp364 = w2*(A_02_4 + A_20_1);
//                                         const Scalar tmp365 = w11*(-A_02_1 - A_02_6 - A_20_3 - A_20_4);
//                                         const Scalar tmp366 = w17*(-A_02_3 - A_20_6);
//                                         const Scalar tmp367 = w2*(A_02_5 - A_20_4);
//                                         const Scalar tmp368 = w6*(-A_12_4 - A_12_5 - A_21_4 - A_21_5);
//                                         const Scalar tmp369 = w11*(-A_02_0 - A_02_7 + A_20_1 + A_20_6);
//                                         const Scalar tmp370 = w20*(-A_01_5 + A_10_4);
//                                         const Scalar tmp371 = w3*(A_11_4 + A_11_5 + A_11_6 + A_11_7);
//                                         const Scalar tmp372 = w12*(-A_12_2 - A_12_3 - A_21_2 - A_21_3);
//                                         const Scalar tmp373 = w21*(-A_01_2 + A_10_3);
//                                         const Scalar tmp374 = w9*(A_11_0 + A_11_1 + A_11_2 + A_11_3);
//                                         const Scalar tmp375 = w29*(A_00_2 + A_00_3);
//                                         const Scalar tmp376 = w8*(A_01_1 + A_01_6 - A_10_0 - A_10_7);
//                                         const Scalar tmp377 = w28*(A_00_4 + A_00_5);
//                                         const Scalar tmp378 = w17*(-A_02_2 + A_20_3);
//                                         const Scalar tmp379 = w17*(A_02_0 + A_20_0);
//                                         const Scalar tmp380 = w2*(-A_02_7 - A_20_7);
//                                         const Scalar tmp381 = w20*(-A_01_7 - A_10_7);
//                                         const Scalar tmp382 = w21*(-A_01_0 - A_10_0);
//                                         const Scalar tmp383 = w6*(A_12_3 + A_12_5 + A_21_3 + A_21_5);
//                                         const Scalar tmp384 = w18*(A_12_0 + A_21_0);
//                                         const Scalar tmp385 = w4*(-A_12_7 - A_21_7);
//                                         const Scalar tmp386 = w12*(A_12_2 + A_12_4 + A_21_2 + A_21_4);
//                                         const Scalar tmp387 = w17*(-A_02_6 - A_20_6);
//                                         const Scalar tmp388 = w2*(A_02_1 + A_20_1);
//                                         const Scalar tmp389 = w20*(A_01_1 + A_10_1);
//                                         const Scalar tmp390 = w21*(A_01_6 + A_10_6);
//                                         const Scalar tmp391 = w18*(A_12_6 + A_21_6);
//                                         const Scalar tmp392 = w4*(-A_12_1 - A_21_1);
//                                         const Scalar tmp393 = w2*(A_02_3 + A_20_6);
//                                         const Scalar tmp394 = w1*(-A_01_3 - A_01_7 + A_10_2 + A_10_6);
//                                         const Scalar tmp395 = w16*(-A_01_0 - A_01_4 + A_10_1 + A_10_5);
//                                         const Scalar tmp396 = w17*(-A_02_4 - A_20_1);
//                                         const Scalar tmp397 = w18*(-A_12_5 - A_21_3);
//                                         const Scalar tmp398 = w10*(A_12_3 + A_12_4 + A_21_2 + A_21_5);
//                                         const Scalar tmp399 = w1*(-A_01_0 - A_01_4 + A_10_2 + A_10_6);
//                                         const Scalar tmp400 = w4*(A_12_2 + A_21_4);
//                                         const Scalar tmp401 = w16*(-A_01_3 - A_01_7 + A_10_1 + A_10_5);
//                                         const Scalar tmp402 = w20*(-A_01_2 + A_10_3);
//                                         const Scalar tmp403 = w21*(-A_01_5 + A_10_4);
//                                         const Scalar tmp404 = w17*(-A_02_5 + A_20_4);
//                                         const Scalar tmp405 = w2*(A_02_2 - A_20_3);
//                                         const Scalar tmp406 = w18*(-A_12_0 + A_21_4);
//                                         const Scalar tmp407 = w4*(A_12_7 - A_21_3);
//                                         const Scalar tmp408 = w17*(-A_02_0 + A_20_4);
//                                         const Scalar tmp409 = w2*(A_02_7 - A_20_3);
//                                         const Scalar tmp410 = w17*(A_02_5 + A_20_5);
//                                         const Scalar tmp411 = w2*(-A_02_2 - A_20_2);
//                                         const Scalar tmp412 = w20*(A_01_2 + A_10_2);
//                                         const Scalar tmp413 = w21*(A_01_5 + A_10_5);
//                                         const Scalar tmp414 = w18*(-A_12_5 - A_21_5);
//                                         const Scalar tmp415 = w4*(A_12_2 + A_21_2);
//                                         const Scalar tmp416 = w12*(-A_12_0 - A_12_1 + A_21_4 + A_21_5);
//                                         const Scalar tmp417 = w6*(-A_12_6 - A_12_7 + A_21_2 + A_21_3);
//                                         const Scalar tmp418 = w17*(A_02_0 + A_20_5);
//                                         const Scalar tmp419 = w2*(-A_02_7 - A_20_2);
//                                         const Scalar tmp420 = w18*(-A_12_4 - A_21_2);
//                                         const Scalar tmp421 = w10*(A_12_2 + A_12_5 + A_21_3 + A_21_4);
//                                         const Scalar tmp422 = w3*(-A_11_1 - A_11_3 - A_11_5 - A_11_7);
//                                         const Scalar tmp423 = w1*(A_01_1 + A_01_5 - A_10_3 - A_10_7);
//                                         const Scalar tmp424 = w25*(-A_22_0 - A_22_2 - A_22_4 - A_22_6);
//                                         const Scalar tmp425 = w4*(A_12_3 + A_21_5);
//                                         const Scalar tmp426 = w15*(A_02_4 + A_02_6 - A_20_0 - A_20_2);
//                                         const Scalar tmp427 = w16*(A_01_2 + A_01_6 - A_10_0 - A_10_4);
//                                         const Scalar tmp428 = w9*(-A_11_0 - A_11_2 - A_11_4 - A_11_6);
//                                         const Scalar tmp429 = w5*(A_02_1 + A_02_3 - A_20_5 - A_20_7);
//                                         const Scalar tmp430 = w23*(-A_22_1 - A_22_3 - A_22_5 - A_22_7);
//                                         const Scalar tmp431 = w18*(A_12_5 - A_21_7);
//                                         const Scalar tmp432 = w10*(-A_12_3 - A_12_4 + A_21_1 + A_21_6);
//                                         const Scalar tmp433 = w21*(A_01_7 - A_10_5);
//                                         const Scalar tmp434 = w20*(A_01_0 - A_10_2);
//                                         const Scalar tmp435 = w4*(-A_12_2 + A_21_0);
//                                         const Scalar tmp436 = w8*(-A_01_3 - A_01_4 + A_10_1 + A_10_6);
//                                         const Scalar tmp437 = w2*(-A_02_4 + A_20_5);
//                                         const Scalar tmp438 = w20*(A_01_4 - A_10_5);
//                                         const Scalar tmp439 = w21*(A_01_3 - A_10_2);
//                                         const Scalar tmp440 = w16*(-A_01_1 - A_01_2 + A_10_0 + A_10_3);
//                                         const Scalar tmp441 = w1*(-A_01_5 - A_01_6 + A_10_4 + A_10_7);
//                                         const Scalar tmp442 = w17*(A_02_3 - A_20_2);
//                                         const Scalar tmp443 = w20*(-A_01_4 - A_10_7);
//                                         const Scalar tmp444 = w21*(-A_01_3 - A_10_0);
//                                         const Scalar tmp445 = w18*(A_12_6 + A_21_0);
//                                         const Scalar tmp446 = w10*(-A_12_0 - A_12_7 - A_21_1 - A_21_6);
//                                         const Scalar tmp447 = w1*(-A_01_3 - A_01_7 + A_10_1 + A_10_5);
//                                         const Scalar tmp448 = w4*(-A_12_1 - A_21_7);
//                                         const Scalar tmp449 = w16*(-A_01_0 - A_01_4 + A_10_2 + A_10_6);
//                                         const Scalar tmp450 = w2*(A_02_7 - A_20_6);
//                                         const Scalar tmp451 = w6*(A_12_6 + A_12_7 + A_21_6 + A_21_7);
//                                         const Scalar tmp452 = w20*(A_01_7 - A_10_6);
//                                         const Scalar tmp453 = w21*(A_01_0 - A_10_1);
//                                         const Scalar tmp454 = w12*(A_12_0 + A_12_1 + A_21_0 + A_21_1);
//                                         const Scalar tmp455 = w29*(A_00_0 + A_00_1);
//                                         const Scalar tmp456 = w28*(A_00_6 + A_00_7);
//                                         const Scalar tmp457 = w17*(-A_02_0 + A_20_1);
//                                         const Scalar tmp458 = w21*(-A_01_7 - A_10_4);
//                                         const Scalar tmp459 = w20*(-A_01_0 - A_10_3);
//                                         const Scalar tmp460 = w12*(A_12_4 + A_12_5 - A_21_6 - A_21_7);
//                                         const Scalar tmp461 = w6*(A_12_2 + A_12_3 - A_21_0 - A_21_1);
//                                         const Scalar tmp462 = w18*(A_12_1 + A_21_7);
//                                         const Scalar tmp463 = w4*(-A_12_6 - A_21_0);
//                                         const Scalar tmp464 = w15*(A_02_1 + A_02_3 - A_20_5 - A_20_7);
//                                         const Scalar tmp465 = w5*(A_02_4 + A_02_6 - A_20_0 - A_20_2);
//                                         const Scalar tmp466 = w2*(-A_02_6 + A_20_7);
//                                         const Scalar tmp467 = w20*(-A_01_6 + A_10_7);
//                                         const Scalar tmp468 = w21*(-A_01_1 + A_10_0);
//                                         const Scalar tmp469 = w17*(A_02_1 - A_20_0);
//                                         const Scalar tmp470 = w6*(-A_12_2 - A_12_3 - A_21_4 - A_21_5);
//                                         const Scalar tmp471 = w1*(-A_01_1 - A_01_5 - A_10_2 - A_10_6);
//                                         const Scalar tmp472 = w15*(-A_02_4 - A_02_6 - A_20_1 - A_20_3);
//                                         const Scalar tmp473 = w16*(-A_01_2 - A_01_6 - A_10_1 - A_10_5);
//                                         const Scalar tmp474 = w12*(-A_12_4 - A_12_5 - A_21_2 - A_21_3);
//                                         const Scalar tmp475 = w5*(-A_02_1 - A_02_3 - A_20_4 - A_20_6);
//                                         const Scalar tmp476 = w18*(-A_12_6 + A_21_4);
//                                         const Scalar tmp477 = w20*(A_01_3 - A_10_1);
//                                         const Scalar tmp478 = w10*(A_12_0 + A_12_7 - A_21_2 - A_21_5);
//                                         const Scalar tmp479 = w4*(A_12_1 - A_21_3);
//                                         const Scalar tmp480 = w21*(A_01_4 - A_10_6);
//                                         const Scalar tmp481 = w8*(-A_01_0 - A_01_7 + A_10_2 + A_10_5);
//                                         const Scalar tmp482 = w6*(A_12_0 + A_12_1 + A_21_6 + A_21_7);
//                                         const Scalar tmp483 = w12*(A_12_6 + A_12_7 + A_21_0 + A_21_1);
//                                         const Scalar tmp484 = w15*(A_02_5 + A_02_7 + A_20_0 + A_20_2);
//                                         const Scalar tmp485 = w5*(A_02_0 + A_02_2 + A_20_5 + A_20_7);
//                                         const Scalar tmp486 = w18*(-A_12_1 + A_21_3);
//                                         const Scalar tmp487 = w20*(A_01_4 - A_10_6);
//                                         const Scalar tmp488 = w4*(A_12_6 - A_21_4);
//                                         const Scalar tmp489 = w21*(A_01_3 - A_10_1);
//                                         const Scalar tmp490 = w20*(A_01_7 - A_10_5);
//                                         const Scalar tmp491 = w18*(A_12_2 - A_21_0);
//                                         const Scalar tmp492 = w4*(-A_12_5 + A_21_7);
//                                         const Scalar tmp493 = w21*(A_01_0 - A_10_2);
//                                         const Scalar tmp494 = w20*(A_01_1 + A_10_2);
//                                         const Scalar tmp495 = w21*(A_01_6 + A_10_5);
//                                         const Scalar tmp496 = w18*(-A_12_2 - A_21_4);
//                                         const Scalar tmp497 = w4*(A_12_5 + A_21_3);
//                                         const Scalar tmp498 = w15*(-A_02_0 - A_02_2 + A_20_4 + A_20_6);
//                                         const Scalar tmp499 = w5*(-A_02_5 - A_02_7 + A_20_1 + A_20_3);
//                                         const Scalar tmp500 = w18*(-A_12_6 + A_21_2);
//                                         const Scalar tmp501 = w4*(A_12_1 - A_21_5);
//                                         const Scalar tmp502 = w17*(A_02_6 - A_20_2);
//                                         const Scalar tmp503 = w2*(-A_02_1 + A_20_5);
//                                         const Scalar tmp504 = w18*(-A_12_3 - A_21_5);
//                                         const Scalar tmp505 = w4*(A_12_4 + A_21_2);
//                                         const Scalar tmp506 = w2*(A_02_6 + A_20_3);
//                                         const Scalar tmp507 = w17*(-A_02_1 - A_20_4);
//                                         const Scalar tmp508 = w18*(A_12_0 + A_21_6);
//                                         const Scalar tmp509 = w4*(-A_12_7 - A_21_1);
//                                         EM_S[INDEX4(k,m,0,0,numEq,numComp,8)]+=tmp198 + tmp200 + tmp214 + tmp259 + tmp262 + tmp289 + tmp294 + tmp299 + tmp303 + tmp304 + tmp307 + tmp309 + tmp343 + tmp347 + tmp362 + tmp363 + tmp379 + tmp380 + tmp381 + tmp382 + tmp383 + tmp384 + tmp385 + tmp386;
//                                         EM_S[INDEX4(k,m,0,1,numEq,numComp,8)]+=tmp161 + tmp201 + tmp247 + tmp250 + tmp371 + tmp374 + tmp44 + tmp451 + tmp454 + tmp455 + tmp456 + tmp466 + tmp467 + tmp468 + tmp469 + tmp49 + tmp89 + tmp91 + tmp92 + tmp98;
//                                         EM_S[INDEX4(k,m,0,2,numEq,numComp,8)]+=tmp135 + tmp236 + tmp238 + tmp240 + tmp242 + tmp244 + tmp39 + tmp41 + tmp432 + tmp436 + tmp440 + tmp441 + tmp490 + tmp491 + tmp492 + tmp493 + tmp61 + tmp68 + tmp70 + tmp71;
//                                         EM_S[INDEX4(k,m,0,3,numEq,numComp,8)]+=tmp114 + tmp165 + tmp166 + tmp167 + tmp168 + tmp169 + tmp170 + tmp171 + tmp172 + tmp20 + tmp73 + tmp74 + tmp75 + tmp76 + tmp79 + tmp80;
//                                         EM_S[INDEX4(k,m,0,4,numEq,numComp,8)]+=tmp1 + tmp127 + tmp131 + tmp141 + tmp145 + tmp146 + tmp148 + tmp15 + tmp189 + tmp190 + tmp192 + tmp193 + tmp2 + tmp243 + tmp246 + tmp406 + tmp407 + tmp408 + tmp409 + tmp5;
//                                         EM_S[INDEX4(k,m,0,5,numEq,numComp,8)]+=tmp174 + tmp176 + tmp184 + tmp24 + tmp260 + tmp267 + tmp339 + tmp340 + tmp341 + tmp342 + tmp344 + tmp345 + tmp416 + tmp417 + tmp506 + tmp507;
//                                         EM_S[INDEX4(k,m,0,6,numEq,numComp,8)]+=tmp21 + tmp258 + tmp266 + tmp274 + tmp337 + tmp398 + tmp422 + tmp424 + tmp428 + tmp430 + tmp447 + tmp449 + tmp496 + tmp497 + tmp498 + tmp499;
//                                         EM_S[INDEX4(k,m,0,7,numEq,numComp,8)]+=tmp104 + tmp105 + tmp106 + tmp107 + tmp108 + tmp109 + tmp110 + tmp111 + tmp112 + tmp113 + tmp38 + tmp87;
//                                         EM_S[INDEX4(k,m,1,0,numEq,numComp,8)]+=tmp145 + tmp148 + tmp161 + tmp201 + tmp202 + tmp210 + tmp371 + tmp374 + tmp440 + tmp441 + tmp450 + tmp451 + tmp452 + tmp453 + tmp454 + tmp455 + tmp456 + tmp457 + tmp89 + tmp91;
//                                         EM_S[INDEX4(k,m,1,1,numEq,numComp,8)]+=tmp215 + tmp221 + tmp227 + tmp260 + tmp267 + tmp288 + tmp304 + tmp312 + tmp317 + tmp351 + tmp352 + tmp353 + tmp354 + tmp355 + tmp356 + tmp357 + tmp358 + tmp359 + tmp360 + tmp361 + tmp362 + tmp363 + tmp76 + tmp79;
//                                         EM_S[INDEX4(k,m,1,2,numEq,numComp,8)]+=tmp166 + tmp169 + tmp172 + tmp196 + tmp197 + tmp198 + tmp199 + tmp20 + tmp200 + tmp21 + tmp73 + tmp74 + tmp75 + tmp77 + tmp80 + tmp82;
//                                         EM_S[INDEX4(k,m,1,3,numEq,numComp,8)]+=tmp36 + tmp37 + tmp38 + tmp39 + tmp40 + tmp41 + tmp42 + tmp43 + tmp44 + tmp45 + tmp46 + tmp47 + tmp48 + tmp49 + tmp50 + tmp51 + tmp52 + tmp53 + tmp54 + tmp55;
//                                         EM_S[INDEX4(k,m,1,4,numEq,numComp,8)]+=tmp176 + tmp24 + tmp269 + tmp274 + tmp339 + tmp340 + tmp342 + tmp343 + tmp344 + tmp347 + tmp394 + tmp395 + tmp416 + tmp417 + tmp418 + tmp419;
//                                         EM_S[INDEX4(k,m,1,5,numEq,numComp,8)]+=tmp112 + tmp12 + tmp123 + tmp13 + tmp141 + tmp142 + tmp143 + tmp146 + tmp147 + tmp149 + tmp16 + tmp277 + tmp278 + tmp279 + tmp280 + tmp281 + tmp282 + tmp6 + tmp92 + tmp98;
//                                         EM_S[INDEX4(k,m,1,6,numEq,numComp,8)]+=tmp104 + tmp105 + tmp106 + tmp110 + tmp113 + tmp135 + tmp136 + tmp137 + tmp138 + tmp139 + tmp15 + tmp87;
//                                         EM_S[INDEX4(k,m,1,7,numEq,numComp,8)]+=tmp114 + tmp184 + tmp225 + tmp232 + tmp329 + tmp330 + tmp332 + tmp334 + tmp335 + tmp337 + tmp338 + tmp421 + tmp464 + tmp465 + tmp504 + tmp505;
//                                         EM_S[INDEX4(k,m,2,0,numEq,numComp,8)]+=tmp135 + tmp234 + tmp235 + tmp236 + tmp237 + tmp238 + tmp239 + tmp240 + tmp241 + tmp242 + tmp243 + tmp244 + tmp245 + tmp246 + tmp39 + tmp41 + tmp44 + tmp49 + tmp61 + tmp71;
//                                         EM_S[INDEX4(k,m,2,1,numEq,numComp,8)]+=tmp114 + tmp120 + tmp167 + tmp170 + tmp198 + tmp20 + tmp200 + tmp24 + tmp443 + tmp444 + tmp73 + tmp74 + tmp75 + tmp80 + tmp81 + tmp83;
//                                         EM_S[INDEX4(k,m,2,2,numEq,numComp,8)]+=tmp217 + tmp231 + tmp233 + tmp258 + tmp266 + tmp271 + tmp273 + tmp288 + tmp289 + tmp290 + tmp291 + tmp292 + tmp293 + tmp294 + tmp295 + tmp296 + tmp297 + tmp298 + tmp299 + tmp300 + tmp301 + tmp302 + tmp76 + tmp79;
//                                         EM_S[INDEX4(k,m,2,3,numEq,numComp,8)]+=tmp101 + tmp156 + tmp157 + tmp204 + tmp205 + tmp368 + tmp371 + tmp372 + tmp374 + tmp375 + tmp377 + tmp437 + tmp438 + tmp439 + tmp440 + tmp441 + tmp442 + tmp85 + tmp87 + tmp99;
//                                         EM_S[INDEX4(k,m,2,4,numEq,numComp,8)]+=tmp184 + tmp21 + tmp328 + tmp337 + tmp383 + tmp386 + tmp422 + tmp423 + tmp424 + tmp427 + tmp428 + tmp430 + tmp498 + tmp499 + tmp508 + tmp509;
//                                         EM_S[INDEX4(k,m,2,5,numEq,numComp,8)]+=tmp104 + tmp106 + tmp108 + tmp111 + tmp113 + tmp15 + tmp160 + tmp161 + tmp162 + tmp163 + tmp164 + tmp38;
//                                         EM_S[INDEX4(k,m,2,6,numEq,numComp,8)]+=tmp10 + tmp112 + tmp122 + tmp123 + tmp124 + tmp125 + tmp126 + tmp127 + tmp128 + tmp129 + tmp130 + tmp131 + tmp132 + tmp133 + tmp134 + tmp14 + tmp3 + tmp68 + tmp70 + tmp9;
//                                         EM_S[INDEX4(k,m,2,7,numEq,numComp,8)]+=tmp166 + tmp175 + tmp176 + tmp178 + tmp179 + tmp180 + tmp183 + tmp187 + tmp270 + tmp272 + tmp274 + tmp284 + tmp285 + tmp364 + tmp365 + tmp366;
//                                         EM_S[INDEX4(k,m,3,0,numEq,numComp,8)]+=tmp20 + tmp21 + tmp24 + tmp34 + tmp72 + tmp73 + tmp74 + tmp75 + tmp76 + tmp77 + tmp78 + tmp79 + tmp80 + tmp81 + tmp82 + tmp83;
//                                         EM_S[INDEX4(k,m,3,1,numEq,numComp,8)]+=tmp13 + tmp16 + tmp38 + tmp39 + tmp40 + tmp41 + tmp43 + tmp440 + tmp441 + tmp45 + tmp47 + tmp478 + tmp481 + tmp486 + tmp487 + tmp488 + tmp489 + tmp50 + tmp52 + tmp55;
//                                         EM_S[INDEX4(k,m,3,2,numEq,numComp,8)]+=tmp101 + tmp14 + tmp204 + tmp205 + tmp367 + tmp368 + tmp369 + tmp370 + tmp371 + tmp372 + tmp373 + tmp374 + tmp375 + tmp376 + tmp377 + tmp378 + tmp44 + tmp49 + tmp87 + tmp9;
//                                         EM_S[INDEX4(k,m,3,3,numEq,numComp,8)]+=tmp179 + tmp183 + tmp198 + tmp200 + tmp214 + tmp215 + tmp216 + tmp217 + tmp218 + tmp219 + tmp220 + tmp221 + tmp222 + tmp223 + tmp224 + tmp225 + tmp226 + tmp227 + tmp228 + tmp229 + tmp230 + tmp231 + tmp232 + tmp233;
//                                         EM_S[INDEX4(k,m,3,4,numEq,numComp,8)]+=tmp104 + tmp106 + tmp112 + tmp113 + tmp135 + tmp137 + tmp139 + tmp160 + tmp161 + tmp164 + tmp471 + tmp473;
//                                         EM_S[INDEX4(k,m,3,5,numEq,numComp,8)]+=tmp114 + tmp274 + tmp312 + tmp317 + tmp329 + tmp332 + tmp335 + tmp337 + tmp338 + tmp399 + tmp401 + tmp446 + tmp462 + tmp463 + tmp464 + tmp465;
//                                         EM_S[INDEX4(k,m,3,6,numEq,numComp,8)]+=tmp166 + tmp175 + tmp176 + tmp177 + tmp178 + tmp180 + tmp181 + tmp184 + tmp187 + tmp271 + tmp273 + tmp283 + tmp284 + tmp285 + tmp286 + tmp287;
//                                         EM_S[INDEX4(k,m,3,7,numEq,numComp,8)]+=tmp1 + tmp10 + tmp11 + tmp12 + tmp15 + tmp152 + tmp153 + tmp154 + tmp155 + tmp156 + tmp157 + tmp158 + tmp159 + tmp17 + tmp3 + tmp4 + tmp51 + tmp54 + tmp6 + tmp7;
//                                         EM_S[INDEX4(k,m,4,0,numEq,numComp,8)]+=tmp1 + tmp127 + tmp131 + tmp141 + tmp146 + tmp15 + tmp153 + tmp154 + tmp188 + tmp189 + tmp190 + tmp191 + tmp192 + tmp193 + tmp194 + tmp195 + tmp68 + tmp70 + tmp92 + tmp98;
//                                         EM_S[INDEX4(k,m,4,1,numEq,numComp,8)]+=tmp166 + tmp176 + tmp184 + tmp283 + tmp339 + tmp340 + tmp341 + tmp342 + tmp343 + tmp344 + tmp345 + tmp346 + tmp347 + tmp348 + tmp349 + tmp350;
//                                         EM_S[INDEX4(k,m,4,2,numEq,numComp,8)]+=tmp114 + tmp274 + tmp337 + tmp383 + tmp386 + tmp422 + tmp424 + tmp426 + tmp428 + tmp429 + tmp430 + tmp445 + tmp446 + tmp447 + tmp448 + tmp449;
//                                         EM_S[INDEX4(k,m,4,3,numEq,numComp,8)]+=tmp104 + tmp106 + tmp107 + tmp109 + tmp112 + tmp113 + tmp135 + tmp161 + tmp482 + tmp483 + tmp484 + tmp485;
//                                         EM_S[INDEX4(k,m,4,4,numEq,numComp,8)]+=tmp118 + tmp121 + tmp214 + tmp215 + tmp216 + tmp217 + tmp220 + tmp222 + tmp253 + tmp254 + tmp255 + tmp256 + tmp257 + tmp258 + tmp259 + tmp260 + tmp261 + tmp262 + tmp263 + tmp264 + tmp265 + tmp266 + tmp267 + tmp268;
//                                         EM_S[INDEX4(k,m,4,5,numEq,numComp,8)]+=tmp100 + tmp101 + tmp145 + tmp148 + tmp369 + tmp376 + tmp402 + tmp403 + tmp404 + tmp405 + tmp60 + tmp65 + tmp84 + tmp87 + tmp88 + tmp89 + tmp91 + tmp95 + tmp96 + tmp97;
//                                         EM_S[INDEX4(k,m,4,6,numEq,numComp,8)]+=tmp243 + tmp246 + tmp38 + tmp43 + tmp476 + tmp477 + tmp478 + tmp479 + tmp480 + tmp481 + tmp57 + tmp58 + tmp61 + tmp63 + tmp64 + tmp66 + tmp69 + tmp71 + tmp90 + tmp94;
//                                         EM_S[INDEX4(k,m,4,7,numEq,numComp,8)]+=tmp20 + tmp21 + tmp22 + tmp23 + tmp24 + tmp25 + tmp26 + tmp27 + tmp28 + tmp29 + tmp30 + tmp31 + tmp32 + tmp33 + tmp34 + tmp35;
//                                         EM_S[INDEX4(k,m,5,0,numEq,numComp,8)]+=tmp166 + tmp176 + tmp260 + tmp267 + tmp274 + tmp339 + tmp340 + tmp342 + tmp344 + tmp346 + tmp348 + tmp365 + tmp393 + tmp394 + tmp395 + tmp396;
//                                         EM_S[INDEX4(k,m,5,1,numEq,numComp,8)]+=tmp112 + tmp12 + tmp123 + tmp124 + tmp126 + tmp140 + tmp141 + tmp142 + tmp143 + tmp144 + tmp145 + tmp146 + tmp147 + tmp148 + tmp149 + tmp150 + tmp151 + tmp51 + tmp54 + tmp6;
//                                         EM_S[INDEX4(k,m,5,2,numEq,numComp,8)]+=tmp104 + tmp106 + tmp113 + tmp136 + tmp138 + tmp15 + tmp161 + tmp38 + tmp472 + tmp475 + tmp482 + tmp483;
//                                         EM_S[INDEX4(k,m,5,3,numEq,numComp,8)]+=tmp184 + tmp21 + tmp312 + tmp317 + tmp327 + tmp328 + tmp329 + tmp330 + tmp331 + tmp332 + tmp333 + tmp334 + tmp335 + tmp336 + tmp337 + tmp338;
//                                         EM_S[INDEX4(k,m,5,4,numEq,numComp,8)]+=tmp100 + tmp101 + tmp102 + tmp103 + tmp84 + tmp85 + tmp86 + tmp87 + tmp88 + tmp89 + tmp90 + tmp91 + tmp92 + tmp93 + tmp94 + tmp95 + tmp96 + tmp97 + tmp98 + tmp99;
//                                         EM_S[INDEX4(k,m,5,5,numEq,numComp,8)]+=tmp217 + tmp225 + tmp232 + tmp26 + tmp265 + tmp268 + tmp288 + tmp289 + tmp29 + tmp290 + tmp293 + tmp295 + tmp308 + tmp313 + tmp343 + tmp347 + tmp358 + tmp359 + tmp410 + tmp411 + tmp412 + tmp413 + tmp414 + tmp415;
//                                         EM_S[INDEX4(k,m,5,6,numEq,numComp,8)]+=tmp114 + tmp115 + tmp116 + tmp117 + tmp118 + tmp119 + tmp120 + tmp121 + tmp20 + tmp22 + tmp24 + tmp25 + tmp28 + tmp30 + tmp32 + tmp35;
//                                         EM_S[INDEX4(k,m,5,7,numEq,numComp,8)]+=tmp13 + tmp135 + tmp16 + tmp237 + tmp238 + tmp245 + tmp319 + tmp320 + tmp321 + tmp322 + tmp323 + tmp324 + tmp325 + tmp326 + tmp45 + tmp55 + tmp57 + tmp60 + tmp64 + tmp65;
//                                         EM_S[INDEX4(k,m,6,0,numEq,numComp,8)]+=tmp114 + tmp184 + tmp258 + tmp266 + tmp337 + tmp420 + tmp421 + tmp422 + tmp423 + tmp424 + tmp425 + tmp426 + tmp427 + tmp428 + tmp429 + tmp430;
//                                         EM_S[INDEX4(k,m,6,1,numEq,numComp,8)]+=tmp104 + tmp106 + tmp113 + tmp135 + tmp15 + tmp162 + tmp163 + tmp470 + tmp474 + tmp484 + tmp485 + tmp87;
//                                         EM_S[INDEX4(k,m,6,2,numEq,numComp,8)]+=tmp10 + tmp112 + tmp123 + tmp125 + tmp127 + tmp128 + tmp130 + tmp131 + tmp132 + tmp156 + tmp157 + tmp243 + tmp246 + tmp278 + tmp279 + tmp3 + tmp500 + tmp501 + tmp502 + tmp503;
//                                         EM_S[INDEX4(k,m,6,3,numEq,numComp,8)]+=tmp175 + tmp176 + tmp178 + tmp180 + tmp182 + tmp185 + tmp187 + tmp24 + tmp269 + tmp270 + tmp271 + tmp272 + tmp273 + tmp274 + tmp275 + tmp276;
//                                         EM_S[INDEX4(k,m,6,4,numEq,numComp,8)]+=tmp38 + tmp42 + tmp43 + tmp53 + tmp56 + tmp57 + tmp58 + tmp59 + tmp60 + tmp61 + tmp62 + tmp63 + tmp64 + tmp65 + tmp66 + tmp67 + tmp68 + tmp69 + tmp70 + tmp71;
//                                         EM_S[INDEX4(k,m,6,5,numEq,numComp,8)]+=tmp118 + tmp121 + tmp166 + tmp199 + tmp20 + tmp21 + tmp22 + tmp25 + tmp27 + tmp28 + tmp30 + tmp33 + tmp458 + tmp459 + tmp460 + tmp461;
//                                         EM_S[INDEX4(k,m,6,6,numEq,numComp,8)]+=tmp179 + tmp183 + tmp215 + tmp255 + tmp26 + tmp261 + tmp288 + tmp29 + tmp298 + tmp300 + tmp304 + tmp316 + tmp318 + tmp351 + tmp354 + tmp355 + tmp383 + tmp386 + tmp387 + tmp388 + tmp389 + tmp390 + tmp391 + tmp392;
//                                         EM_S[INDEX4(k,m,6,7,numEq,numComp,8)]+=tmp100 + tmp14 + tmp161 + tmp201 + tmp202 + tmp203 + tmp204 + tmp205 + tmp206 + tmp207 + tmp208 + tmp209 + tmp210 + tmp211 + tmp212 + tmp213 + tmp88 + tmp9 + tmp90 + tmp94;
//                                         EM_S[INDEX4(k,m,7,0,numEq,numComp,8)]+=tmp104 + tmp106 + tmp112 + tmp113 + tmp38 + tmp470 + tmp471 + tmp472 + tmp473 + tmp474 + tmp475 + tmp87;
//                                         EM_S[INDEX4(k,m,7,1,numEq,numComp,8)]+=tmp21 + tmp225 + tmp232 + tmp274 + tmp329 + tmp332 + tmp333 + tmp335 + tmp336 + tmp337 + tmp338 + tmp397 + tmp398 + tmp399 + tmp400 + tmp401;
//                                         EM_S[INDEX4(k,m,7,2,numEq,numComp,8)]+=tmp173 + tmp174 + tmp175 + tmp176 + tmp177 + tmp178 + tmp179 + tmp180 + tmp181 + tmp182 + tmp183 + tmp184 + tmp185 + tmp186 + tmp187 + tmp24;
//                                         EM_S[INDEX4(k,m,7,3,numEq,numComp,8)]+=tmp0 + tmp1 + tmp10 + tmp11 + tmp12 + tmp13 + tmp14 + tmp15 + tmp16 + tmp17 + tmp18 + tmp19 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6 + tmp7 + tmp8 + tmp9;
//                                         EM_S[INDEX4(k,m,7,4,numEq,numComp,8)]+=tmp114 + tmp117 + tmp119 + tmp166 + tmp171 + tmp20 + tmp22 + tmp25 + tmp26 + tmp28 + tmp29 + tmp30 + tmp460 + tmp461 + tmp494 + tmp495;
//                                         EM_S[INDEX4(k,m,7,5,numEq,numComp,8)]+=tmp135 + tmp238 + tmp320 + tmp324 + tmp325 + tmp326 + tmp431 + tmp432 + tmp433 + tmp434 + tmp435 + tmp436 + tmp45 + tmp51 + tmp54 + tmp55 + tmp57 + tmp64 + tmp90 + tmp94;
//                                         EM_S[INDEX4(k,m,7,6,numEq,numComp,8)]+=tmp100 + tmp156 + tmp157 + tmp161 + tmp201 + tmp204 + tmp205 + tmp207 + tmp208 + tmp209 + tmp211 + tmp247 + tmp248 + tmp249 + tmp250 + tmp251 + tmp252 + tmp60 + tmp65 + tmp88;
//                                         EM_S[INDEX4(k,m,7,7,numEq,numComp,8)]+=tmp118 + tmp121 + tmp214 + tmp226 + tmp228 + tmp271 + tmp273 + tmp289 + tmp303 + tmp304 + tmp305 + tmp306 + tmp307 + tmp308 + tmp309 + tmp310 + tmp311 + tmp312 + tmp313 + tmp314 + tmp315 + tmp316 + tmp317 + tmp318;
//                                     }
//                                 }
//                             } else { // constant data
//                                 for (index_t k=0; k<numEq; k++) {
//                                     for (index_t m=0; m<numComp; m++) {
//                                         const Scalar Aw00 =  8.*A_p[INDEX4(k,0,m,0, numEq,3, numComp)]*w27;
//                                         const Scalar Aw01 = 12.*A_p[INDEX4(k,0,m,1, numEq,3, numComp)]*w8;
//                                         const Scalar Aw02 = 12.*A_p[INDEX4(k,0,m,2, numEq,3, numComp)]*w11;
//                                         const Scalar Aw10 = 12.*A_p[INDEX4(k,1,m,0, numEq,3, numComp)]*w8;
//                                         const Scalar Aw11 =  8.*A_p[INDEX4(k,1,m,1, numEq,3, numComp)]*w22;
//                                         const Scalar Aw12 = 12.*A_p[INDEX4(k,1,m,2, numEq,3, numComp)]*w10;
//                                         const Scalar Aw20 = 12.*A_p[INDEX4(k,2,m,0, numEq,3, numComp)]*w11;
//                                         const Scalar Aw21 = 12.*A_p[INDEX4(k,2,m,1, numEq,3, numComp)]*w10;
//                                         const Scalar Aw22 =  8.*A_p[INDEX4(k,2,m,2, numEq,3, numComp)]*w13;
//                                         const Scalar tmp0 = Aw01 + Aw10;
//                                         const Scalar tmp1 = Aw01 - Aw10;
//                                         const Scalar tmp2 = Aw02 + Aw20;
//                                         const Scalar tmp3 = Aw02 - Aw20;
//                                         const Scalar tmp4 = Aw12 + Aw21;
//                                         const Scalar tmp5 = Aw12 - Aw21;
//                                         EM_S[INDEX4(k,m,0,0,numEq,numComp,8)]+=-4.*Aw00 - 4.*Aw11 - 4.*Aw22 + 2.*tmp0 + 2.*tmp2 - 2.*tmp4;
//                                         EM_S[INDEX4(k,m,0,1,numEq,numComp,8)]+= 4.*Aw00 - 2.*Aw11 - 2.*Aw22 - tmp4 + 2.*tmp1 + 2.*tmp3;
//                                         EM_S[INDEX4(k,m,0,2,numEq,numComp,8)]+=-2.*Aw00 + 4.*Aw11 - 2.*Aw22 - 2.*tmp1 + tmp2 - 2.*tmp5;
//                                         EM_S[INDEX4(k,m,0,3,numEq,numComp,8)]+= 2.*Aw00 + 2.*Aw11 -    Aw22 - 2.*tmp0 + tmp3 - tmp5;
//                                         EM_S[INDEX4(k,m,0,4,numEq,numComp,8)]+=-2.*Aw00 - 2.*Aw11 + 4.*Aw22 - 2.*tmp3 + 2.*tmp5 + tmp0;
//                                         EM_S[INDEX4(k,m,0,5,numEq,numComp,8)]+= 2.*Aw00 -    Aw11 + 2.*Aw22 - 2.*tmp2 + tmp1 + tmp5;
//                                         EM_S[INDEX4(k,m,0,6,numEq,numComp,8)]+=   -Aw00 + 2.*Aw11 + 2.*Aw22 + 2.*tmp4 - tmp1 - tmp3;
//                                         EM_S[INDEX4(k,m,0,7,numEq,numComp,8)]+=    Aw00 +    Aw11 +    Aw22 + tmp4 - tmp0 - tmp2;
//                                         EM_S[INDEX4(k,m,1,0,numEq,numComp,8)]+= 4.*Aw00 - 2.*Aw11 - 2.*Aw22 - 2.*tmp3 - 2.*tmp1 - tmp4;
//                                         EM_S[INDEX4(k,m,1,1,numEq,numComp,8)]+=-4.*Aw00 - 4.*Aw11 - 4.*Aw22 - 2.*tmp2 - 2.*tmp4 - 2.*tmp0;
//                                         EM_S[INDEX4(k,m,1,2,numEq,numComp,8)]+= 2.*Aw00 + 2.*Aw11 -    Aw22 + 2.*tmp0 - tmp5 - tmp3;
//                                         EM_S[INDEX4(k,m,1,3,numEq,numComp,8)]+=-2.*Aw00 + 4.*Aw11 - 2.*Aw22 - tmp2 - 2.*tmp5 + 2.*tmp1;
//                                         EM_S[INDEX4(k,m,1,4,numEq,numComp,8)]+= 2.*Aw00 -    Aw11 + 2.*Aw22 + 2.*tmp2 - tmp1 + tmp5;
//                                         EM_S[INDEX4(k,m,1,5,numEq,numComp,8)]+=-2.*Aw00 - 2.*Aw11 + 4.*Aw22 + 2.*tmp5 - tmp0 + 2.*tmp3;
//                                         EM_S[INDEX4(k,m,1,6,numEq,numComp,8)]+=    Aw00 +    Aw11 +    Aw22 + tmp4 + tmp2 + tmp0;
//                                         EM_S[INDEX4(k,m,1,7,numEq,numComp,8)]+=   -Aw00 + 2.*Aw11 + 2.*Aw22 + tmp3 + tmp1 + 2.*tmp4;
//                                         EM_S[INDEX4(k,m,2,0,numEq,numComp,8)]+=-2.*Aw00 + 4.*Aw11 - 2.*Aw22 + 2.*tmp5 + tmp2 + 2.*tmp1;
//                                         EM_S[INDEX4(k,m,2,1,numEq,numComp,8)]+= 2.*Aw00 + 2.*Aw11 -    Aw22 + tmp3 + 2.*tmp0 + tmp5;
//                                         EM_S[INDEX4(k,m,2,2,numEq,numComp,8)]+=-4.*Aw00 - 4.*Aw11 - 4.*Aw22 + 2.*tmp4 + 2.*tmp2 - 2.*tmp0;
//                                         EM_S[INDEX4(k,m,2,3,numEq,numComp,8)]+= 4.*Aw00 - 2.*Aw11 - 2.*Aw22 + tmp4 - 2.*tmp1 + 2.*tmp3;
//                                         EM_S[INDEX4(k,m,2,4,numEq,numComp,8)]+=   -Aw00 + 2.*Aw11 + 2.*Aw22 + tmp1 - 2.*tmp4 - tmp3;
//                                         EM_S[INDEX4(k,m,2,5,numEq,numComp,8)]+=    Aw00 +    Aw11 +    Aw22 - tmp4 + tmp0 - tmp2;
//                                         EM_S[INDEX4(k,m,2,6,numEq,numComp,8)]+=-2.*Aw00 - 2.*Aw11 + 4.*Aw22 - 2.*tmp3 - tmp0 - 2.*tmp5;
//                                         EM_S[INDEX4(k,m,2,7,numEq,numComp,8)]+= 2.*Aw00 -    Aw11 + 2.*Aw22 - tmp5 - 2.*tmp2 - tmp1;
//                                         EM_S[INDEX4(k,m,3,0,numEq,numComp,8)]+= 2.*Aw00 + 2.*Aw11 -    Aw22 - tmp3 + tmp5 - 2.*tmp0;
//                                         EM_S[INDEX4(k,m,3,1,numEq,numComp,8)]+=-2.*Aw00 + 4.*Aw11 - 2.*Aw22 + 2.*tmp5 - 2.*tmp1 - tmp2;
//                                         EM_S[INDEX4(k,m,3,2,numEq,numComp,8)]+= 4.*Aw00 - 2.*Aw11 - 2.*Aw22 - 2.*tmp3 + tmp4 + 2.*tmp1;
//                                         EM_S[INDEX4(k,m,3,3,numEq,numComp,8)]+=-4.*Aw00 - 4.*Aw11 - 4.*Aw22 + 2.*tmp0 - 2.*tmp2 + 2.*tmp4;
//                                         EM_S[INDEX4(k,m,3,4,numEq,numComp,8)]+=    Aw00 +    Aw11 +    Aw22 - tmp0 + tmp2 - tmp4;
//                                         EM_S[INDEX4(k,m,3,5,numEq,numComp,8)]+=   -Aw00 + 2.*Aw11 + 2.*Aw22 + tmp3 - tmp1 - 2.*tmp4;
//                                         EM_S[INDEX4(k,m,3,6,numEq,numComp,8)]+= 2.*Aw00 -    Aw11 + 2.*Aw22 - tmp5 + tmp1 + 2.*tmp2;
//                                         EM_S[INDEX4(k,m,3,7,numEq,numComp,8)]+=-2.*Aw00 - 2.*Aw11 + 4.*Aw22 + tmp0 - 2.*tmp5 + 2.*tmp3;
//                                         EM_S[INDEX4(k,m,4,0,numEq,numComp,8)]+=-2.*Aw00 - 2.*Aw11 + 4.*Aw22 + tmp0 - 2.*tmp5 + 2.*tmp3;
//                                         EM_S[INDEX4(k,m,4,1,numEq,numComp,8)]+= 2.*Aw00 -    Aw11 + 2.*Aw22 - tmp5 + tmp1 + 2.*tmp2;
//                                         EM_S[INDEX4(k,m,4,2,numEq,numComp,8)]+=   -Aw00 + 2.*Aw11 + 2.*Aw22 + tmp3 - tmp1 - 2.*tmp4;
//                                         EM_S[INDEX4(k,m,4,3,numEq,numComp,8)]+=    Aw00 +    Aw11 +    Aw22 - tmp0 + tmp2 - tmp4;
//                                         EM_S[INDEX4(k,m,4,4,numEq,numComp,8)]+=-4.*Aw00 - 4.*Aw11 - 4.*Aw22 + 2.*tmp0 - 2.*tmp2 + 2.*tmp4;
//                                         EM_S[INDEX4(k,m,4,5,numEq,numComp,8)]+= 4.*Aw00 - 2.*Aw11 - 2.*Aw22 - 2.*tmp3 + tmp4 + 2.*tmp1;
//                                         EM_S[INDEX4(k,m,4,6,numEq,numComp,8)]+=-2.*Aw00 + 4.*Aw11 - 2.*Aw22 + 2.*tmp5 - 2.*tmp1 - tmp2;
//                                         EM_S[INDEX4(k,m,4,7,numEq,numComp,8)]+= 2.*Aw00 + 2.*Aw11 -    Aw22 - tmp3 + tmp5 - 2.*tmp0;
//                                         EM_S[INDEX4(k,m,5,0,numEq,numComp,8)]+= 2.*Aw00 -    Aw11 + 2.*Aw22 - tmp5 - 2.*tmp2 - tmp1;
//                                         EM_S[INDEX4(k,m,5,1,numEq,numComp,8)]+=-2.*Aw00 - 2.*Aw11 + 4.*Aw22 - 2.*tmp3 - tmp0 - 2.*tmp5;
//                                         EM_S[INDEX4(k,m,5,2,numEq,numComp,8)]+=    Aw00 +    Aw11 +    Aw22 - tmp4 + tmp0 - tmp2;
//                                         EM_S[INDEX4(k,m,5,3,numEq,numComp,8)]+=   -Aw00 + 2.*Aw11 + 2.*Aw22 + tmp1 - 2.*tmp4 - tmp3;
//                                         EM_S[INDEX4(k,m,5,4,numEq,numComp,8)]+= 4.*Aw00 - 2.*Aw11 - 2.*Aw22 + tmp4 - 2.*tmp1 + 2.*tmp3;
//                                         EM_S[INDEX4(k,m,5,5,numEq,numComp,8)]+=-4.*Aw00 - 4.*Aw11 - 4.*Aw22 + 2.*tmp4 + 2.*tmp2 - 2.*tmp0;
//                                         EM_S[INDEX4(k,m,5,6,numEq,numComp,8)]+= 2.*Aw00 + 2.*Aw11 -    Aw22 + tmp3 + 2.*tmp0 + tmp5;
//                                         EM_S[INDEX4(k,m,5,7,numEq,numComp,8)]+=-2.*Aw00 + 4.*Aw11 - 2.*Aw22 + 2.*tmp5 + tmp2 + 2.*tmp1;
//                                         EM_S[INDEX4(k,m,6,0,numEq,numComp,8)]+=   -Aw00 + 2.*Aw11 + 2.*Aw22 + tmp3 + tmp1 + 2.*tmp4;
//                                         EM_S[INDEX4(k,m,6,1,numEq,numComp,8)]+=    Aw00 +    Aw11 +    Aw22 + tmp4 + tmp2 + tmp0;
//                                         EM_S[INDEX4(k,m,6,2,numEq,numComp,8)]+=-2.*Aw00 - 2.*Aw11 + 4.*Aw22 + 2.*tmp5 - tmp0 + 2.*tmp3;
//                                         EM_S[INDEX4(k,m,6,3,numEq,numComp,8)]+= 2.*Aw00 -    Aw11 + 2.*Aw22 + 2.*tmp2 - tmp1 + tmp5;
//                                         EM_S[INDEX4(k,m,6,4,numEq,numComp,8)]+=-2.*Aw00 + 4.*Aw11 - 2.*Aw22 - tmp2 - 2.*tmp5 + 2.*tmp1;
//                                         EM_S[INDEX4(k,m,6,5,numEq,numComp,8)]+= 2.*Aw00 + 2.*Aw11 -    Aw22 + 2.*tmp0 - tmp5 - tmp3;
//                                         EM_S[INDEX4(k,m,6,6,numEq,numComp,8)]+=-4.*Aw00 - 4.*Aw11 - 4.*Aw22 - 2.*tmp2 - 2.*tmp4 - 2.*tmp0;
//                                         EM_S[INDEX4(k,m,6,7,numEq,numComp,8)]+= 4.*Aw00 - 2.*Aw11 - 2.*Aw22 - 2.*tmp3 - 2.*tmp1 - tmp4;
//                                         EM_S[INDEX4(k,m,7,0,numEq,numComp,8)]+=    Aw00 +    Aw11 +    Aw22 + tmp4 - tmp0 - tmp2;
//                                         EM_S[INDEX4(k,m,7,1,numEq,numComp,8)]+=   -Aw00 + 2.*Aw11 + 2.*Aw22 + 2.*tmp4 - tmp1 - tmp3;
//                                         EM_S[INDEX4(k,m,7,2,numEq,numComp,8)]+= 2.*Aw00 -    Aw11 + 2.*Aw22 - 2.*tmp2 + tmp1 + tmp5;
//                                         EM_S[INDEX4(k,m,7,3,numEq,numComp,8)]+=-2.*Aw00 - 2.*Aw11 + 4.*Aw22 - 2.*tmp3 + 2.*tmp5 + tmp0;
//                                         EM_S[INDEX4(k,m,7,4,numEq,numComp,8)]+= 2.*Aw00 + 2.*Aw11 -    Aw22 + tmp3 - tmp5 - 2.*tmp0;
//                                         EM_S[INDEX4(k,m,7,5,numEq,numComp,8)]+=-2.*Aw00 + 4.*Aw11 - 2.*Aw22 - 2.*tmp1 + tmp2 - 2.*tmp5;
//                                         EM_S[INDEX4(k,m,7,6,numEq,numComp,8)]+= 4.*Aw00 - 2.*Aw11 - 2.*Aw22 - tmp4 + 2.*tmp1 + 2.*tmp3;
//                                         EM_S[INDEX4(k,m,7,7,numEq,numComp,8)]+=-4.*Aw00 - 4.*Aw11 - 4.*Aw22 + 2.*tmp0 + 2.*tmp2 - 2.*tmp4;
//                                     }
//                                 }
//                             }
//                         }
//                         ///////////////
//                         // process B //
//                         ///////////////
//                         if (!B.isEmpty()) {
//                             const Scalar* B_p = B.getSampleDataRO(id, zero);
//                             if (B.actsExpanded()) {
//                                 for (index_t k=0; k<numEq; k++) {
//                                     for (index_t m=0; m<numComp; m++) {
//                                         const Scalar B_0_0 = B_p[INDEX4(k,0,m,0, numEq,3,numComp)];
//                                         const Scalar B_1_0 = B_p[INDEX4(k,1,m,0, numEq,3,numComp)];
//                                         const Scalar B_2_0 = B_p[INDEX4(k,2,m,0, numEq,3,numComp)];
//                                         const Scalar B_0_1 = B_p[INDEX4(k,0,m,1, numEq,3,numComp)];
//                                         const Scalar B_1_1 = B_p[INDEX4(k,1,m,1, numEq,3,numComp)];
//                                         const Scalar B_2_1 = B_p[INDEX4(k,2,m,1, numEq,3,numComp)];
//                                         const Scalar B_0_2 = B_p[INDEX4(k,0,m,2, numEq,3,numComp)];
//                                         const Scalar B_1_2 = B_p[INDEX4(k,1,m,2, numEq,3,numComp)];
//                                         const Scalar B_2_2 = B_p[INDEX4(k,2,m,2, numEq,3,numComp)];
//                                         const Scalar B_0_3 = B_p[INDEX4(k,0,m,3, numEq,3,numComp)];
//                                         const Scalar B_1_3 = B_p[INDEX4(k,1,m,3, numEq,3,numComp)];
//                                         const Scalar B_2_3 = B_p[INDEX4(k,2,m,3, numEq,3,numComp)];
//                                         const Scalar B_0_4 = B_p[INDEX4(k,0,m,4, numEq,3,numComp)];
//                                         const Scalar B_1_4 = B_p[INDEX4(k,1,m,4, numEq,3,numComp)];
//                                         const Scalar B_2_4 = B_p[INDEX4(k,2,m,4, numEq,3,numComp)];
//                                         const Scalar B_0_5 = B_p[INDEX4(k,0,m,5, numEq,3,numComp)];
//                                         const Scalar B_1_5 = B_p[INDEX4(k,1,m,5, numEq,3,numComp)];
//                                         const Scalar B_2_5 = B_p[INDEX4(k,2,m,5, numEq,3,numComp)];
//                                         const Scalar B_0_6 = B_p[INDEX4(k,0,m,6, numEq,3,numComp)];
//                                         const Scalar B_1_6 = B_p[INDEX4(k,1,m,6, numEq,3,numComp)];
//                                         const Scalar B_2_6 = B_p[INDEX4(k,2,m,6, numEq,3,numComp)];
//                                         const Scalar B_0_7 = B_p[INDEX4(k,0,m,7, numEq,3,numComp)];
//                                         const Scalar B_1_7 = B_p[INDEX4(k,1,m,7, numEq,3,numComp)];
//                                         const Scalar B_2_7 = B_p[INDEX4(k,2,m,7, numEq,3,numComp)];
//                                         const Scalar tmp0 = w38*(B_2_1 + B_2_2);
//                                         const Scalar tmp1 = w42*(B_1_3 + B_1_7);
//                                         const Scalar tmp2 = w41*(B_0_3 + B_0_7);
//                                         const Scalar tmp3 = w37*(B_1_1 + B_1_5);
//                                         const Scalar tmp4 = w39*(B_0_2 + B_0_6);
//                                         const Scalar tmp5 = w45*(B_2_5 + B_2_6);
//                                         const Scalar tmp6 = w36*(B_0_1 + B_0_5);
//                                         const Scalar tmp7 = w40*(B_1_2 + B_1_6);
//                                         const Scalar tmp8 = w33*(B_0_0 + B_0_4);
//                                         const Scalar tmp9 = w34*(B_1_0 + B_1_4);
//                                         const Scalar tmp10 = w38*(B_2_4 + B_2_5 + B_2_6 + B_2_7);
//                                         const Scalar tmp11 = w42*(-B_1_6 - B_1_7);
//                                         const Scalar tmp12 = w41*(-B_0_5 - B_0_7);
//                                         const Scalar tmp13 = w37*(-B_1_4 - B_1_5);
//                                         const Scalar tmp14 = w39*(-B_0_4 - B_0_6);
//                                         const Scalar tmp15 = w45*(B_2_0 + B_2_1 + B_2_2 + B_2_3);
//                                         const Scalar tmp16 = w36*(-B_0_1 - B_0_3);
//                                         const Scalar tmp17 = w40*(-B_1_2 - B_1_3);
//                                         const Scalar tmp18 = w33*(-B_0_0 - B_0_2);
//                                         const Scalar tmp19 = w34*(-B_1_0 - B_1_1);
//                                         const Scalar tmp20 = w38*(-B_2_5 - B_2_7);
//                                         const Scalar tmp21 = w35*(-B_2_4 - B_2_6);
//                                         const Scalar tmp22 = w41*(B_0_1 + B_0_3);
//                                         const Scalar tmp23 = w37*(-B_1_2 - B_1_7);
//                                         const Scalar tmp24 = w39*(B_0_0 + B_0_2);
//                                         const Scalar tmp25 = w45*(-B_2_0 - B_2_2);
//                                         const Scalar tmp26 = w36*(B_0_5 + B_0_7);
//                                         const Scalar tmp27 = w40*(-B_1_0 - B_1_5);
//                                         const Scalar tmp28 = w33*(B_0_4 + B_0_6);
//                                         const Scalar tmp29 = w46*(-B_2_1 - B_2_3);
//                                         const Scalar tmp30 = w38*(B_2_0 + B_2_2);
//                                         const Scalar tmp31 = w35*(B_2_1 + B_2_3);
//                                         const Scalar tmp32 = w41*(-B_0_4 - B_0_6);
//                                         const Scalar tmp33 = w37*(B_1_0 + B_1_5);
//                                         const Scalar tmp34 = w39*(-B_0_5 - B_0_7);
//                                         const Scalar tmp35 = w45*(B_2_5 + B_2_7);
//                                         const Scalar tmp36 = w36*(-B_0_0 - B_0_2);
//                                         const Scalar tmp37 = w40*(B_1_2 + B_1_7);
//                                         const Scalar tmp38 = w33*(-B_0_1 - B_0_3);
//                                         const Scalar tmp39 = w46*(B_2_4 + B_2_6);
//                                         const Scalar tmp40 = w38*(-B_2_0 - B_2_1 - B_2_2 - B_2_3);
//                                         const Scalar tmp41 = w42*(B_1_0 + B_1_1);
//                                         const Scalar tmp42 = w41*(B_0_0 + B_0_2);
//                                         const Scalar tmp43 = w37*(B_1_2 + B_1_3);
//                                         const Scalar tmp44 = w39*(B_0_1 + B_0_3);
//                                         const Scalar tmp45 = w45*(-B_2_4 - B_2_5 - B_2_6 - B_2_7);
//                                         const Scalar tmp46 = w36*(B_0_4 + B_0_6);
//                                         const Scalar tmp47 = w40*(B_1_4 + B_1_5);
//                                         const Scalar tmp48 = w33*(B_0_5 + B_0_7);
//                                         const Scalar tmp49 = w34*(B_1_6 + B_1_7);
//                                         const Scalar tmp50 = w38*(B_2_0 + B_2_1);
//                                         const Scalar tmp51 = w42*(-B_1_4 - B_1_5);
//                                         const Scalar tmp52 = w35*(B_2_2 + B_2_3);
//                                         const Scalar tmp53 = w37*(-B_1_6 - B_1_7);
//                                         const Scalar tmp54 = w39*(B_0_0 + B_0_6);
//                                         const Scalar tmp55 = w45*(B_2_6 + B_2_7);
//                                         const Scalar tmp56 = w36*(B_0_1 + B_0_7);
//                                         const Scalar tmp57 = w40*(-B_1_0 - B_1_1);
//                                         const Scalar tmp58 = w46*(B_2_4 + B_2_5);
//                                         const Scalar tmp59 = w34*(-B_1_2 - B_1_3);
//                                         const Scalar tmp60 = w38*(-B_2_4 - B_2_5 - B_2_6 - B_2_7);
//                                         const Scalar tmp61 = w37*(-B_1_2 - B_1_3 - B_1_6 - B_1_7);
//                                         const Scalar tmp62 = w39*(-B_0_1 - B_0_3 - B_0_5 - B_0_7);
//                                         const Scalar tmp63 = w45*(-B_2_0 - B_2_1 - B_2_2 - B_2_3);
//                                         const Scalar tmp64 = w36*(-B_0_0 - B_0_2 - B_0_4 - B_0_6);
//                                         const Scalar tmp65 = w40*(-B_1_0 - B_1_1 - B_1_4 - B_1_5);
//                                         const Scalar tmp66 = w41*(B_0_4 + B_0_6);
//                                         const Scalar tmp67 = w39*(B_0_5 + B_0_7);
//                                         const Scalar tmp68 = w36*(B_0_0 + B_0_2);
//                                         const Scalar tmp69 = w33*(B_0_1 + B_0_3);
//                                         const Scalar tmp70 = w38*(-B_2_4 - B_2_7);
//                                         const Scalar tmp71 = w42*(B_1_2 + B_1_6);
//                                         const Scalar tmp72 = w41*(-B_0_2 - B_0_6);
//                                         const Scalar tmp73 = w37*(B_1_0 + B_1_4);
//                                         const Scalar tmp74 = w39*(-B_0_3 - B_0_7);
//                                         const Scalar tmp75 = w45*(-B_2_0 - B_2_3);
//                                         const Scalar tmp76 = w36*(-B_0_0 - B_0_4);
//                                         const Scalar tmp77 = w40*(B_1_3 + B_1_7);
//                                         const Scalar tmp78 = w33*(-B_0_1 - B_0_5);
//                                         const Scalar tmp79 = w34*(B_1_1 + B_1_5);
//                                         const Scalar tmp80 = w39*(B_0_0 + B_0_2 + B_0_4 + B_0_6);
//                                         const Scalar tmp81 = w36*(B_0_1 + B_0_3 + B_0_5 + B_0_7);
//                                         const Scalar tmp82 = w38*(B_2_0 + B_2_3);
//                                         const Scalar tmp83 = w42*(-B_1_1 - B_1_5);
//                                         const Scalar tmp84 = w41*(B_0_1 + B_0_5);
//                                         const Scalar tmp85 = w37*(-B_1_3 - B_1_7);
//                                         const Scalar tmp86 = w39*(B_0_0 + B_0_4);
//                                         const Scalar tmp87 = w45*(B_2_4 + B_2_7);
//                                         const Scalar tmp88 = w36*(B_0_3 + B_0_7);
//                                         const Scalar tmp89 = w40*(-B_1_0 - B_1_4);
//                                         const Scalar tmp90 = w33*(B_0_2 + B_0_6);
//                                         const Scalar tmp91 = w34*(-B_1_2 - B_1_6);
//                                         const Scalar tmp92 = w38*(-B_2_5 - B_2_6);
//                                         const Scalar tmp93 = w45*(-B_2_1 - B_2_2);
//                                         const Scalar tmp94 = w37*(B_1_0 + B_1_1 + B_1_4 + B_1_5);
//                                         const Scalar tmp95 = w40*(B_1_2 + B_1_3 + B_1_6 + B_1_7);
//                                         const Scalar tmp96 = w42*(-B_1_2 - B_1_3);
//                                         const Scalar tmp97 = w41*(-B_0_1 - B_0_3);
//                                         const Scalar tmp98 = w37*(-B_1_0 - B_1_1);
//                                         const Scalar tmp99 = w39*(-B_0_0 - B_0_2);
//                                         const Scalar tmp100 = w36*(-B_0_5 - B_0_7);
//                                         const Scalar tmp101 = w40*(-B_1_6 - B_1_7);
//                                         const Scalar tmp102 = w33*(-B_0_4 - B_0_6);
//                                         const Scalar tmp103 = w34*(-B_1_4 - B_1_5);
//                                         const Scalar tmp104 = w38*(B_2_6 + B_2_7);
//                                         const Scalar tmp105 = w35*(B_2_4 + B_2_5);
//                                         const Scalar tmp106 = w41*(B_0_2 + B_0_6);
//                                         const Scalar tmp107 = w37*(B_1_2 + B_1_3 + B_1_6 + B_1_7);
//                                         const Scalar tmp108 = w39*(B_0_3 + B_0_7);
//                                         const Scalar tmp109 = w45*(B_2_0 + B_2_1);
//                                         const Scalar tmp110 = w36*(B_0_0 + B_0_4);
//                                         const Scalar tmp111 = w40*(B_1_0 + B_1_1 + B_1_4 + B_1_5);
//                                         const Scalar tmp112 = w33*(B_0_1 + B_0_5);
//                                         const Scalar tmp113 = w46*(B_2_2 + B_2_3);
//                                         const Scalar tmp114 = w42*(-B_1_0 - B_1_4);
//                                         const Scalar tmp115 = w41*(-B_0_0 - B_0_4);
//                                         const Scalar tmp116 = w37*(-B_1_2 - B_1_6);
//                                         const Scalar tmp117 = w39*(-B_0_1 - B_0_5);
//                                         const Scalar tmp118 = w36*(-B_0_2 - B_0_6);
//                                         const Scalar tmp119 = w40*(-B_1_1 - B_1_5);
//                                         const Scalar tmp120 = w33*(-B_0_3 - B_0_7);
//                                         const Scalar tmp121 = w34*(-B_1_3 - B_1_7);
//                                         const Scalar tmp122 = w38*(B_2_2 + B_2_3);
//                                         const Scalar tmp123 = w42*(B_1_6 + B_1_7);
//                                         const Scalar tmp124 = w35*(B_2_0 + B_2_1);
//                                         const Scalar tmp125 = w37*(B_1_4 + B_1_5);
//                                         const Scalar tmp126 = w39*(-B_0_3 - B_0_5);
//                                         const Scalar tmp127 = w45*(B_2_4 + B_2_5);
//                                         const Scalar tmp128 = w36*(-B_0_2 - B_0_4);
//                                         const Scalar tmp129 = w40*(B_1_2 + B_1_3);
//                                         const Scalar tmp130 = w46*(B_2_6 + B_2_7);
//                                         const Scalar tmp131 = w34*(B_1_0 + B_1_1);
//                                         const Scalar tmp132 = w38*(-B_2_1 - B_2_2);
//                                         const Scalar tmp133 = w37*(B_1_2 + B_1_7);
//                                         const Scalar tmp134 = w39*(B_0_1 + B_0_7);
//                                         const Scalar tmp135 = w36*(B_0_0 + B_0_6);
//                                         const Scalar tmp136 = w40*(B_1_0 + B_1_5);
//                                         const Scalar tmp137 = w45*(-B_2_5 - B_2_6);
//                                         const Scalar tmp138 = w38*(-B_2_4 - B_2_6);
//                                         const Scalar tmp139 = w35*(-B_2_5 - B_2_7);
//                                         const Scalar tmp140 = w41*(-B_0_0 - B_0_2);
//                                         const Scalar tmp141 = w37*(B_1_1 + B_1_4);
//                                         const Scalar tmp142 = w39*(-B_0_1 - B_0_3);
//                                         const Scalar tmp143 = w45*(-B_2_1 - B_2_3);
//                                         const Scalar tmp144 = w36*(-B_0_4 - B_0_6);
//                                         const Scalar tmp145 = w40*(B_1_3 + B_1_6);
//                                         const Scalar tmp146 = w33*(-B_0_5 - B_0_7);
//                                         const Scalar tmp147 = w46*(-B_2_0 - B_2_2);
//                                         const Scalar tmp148 = w39*(B_0_2 + B_0_4);
//                                         const Scalar tmp149 = w36*(B_0_3 + B_0_5);
//                                         const Scalar tmp150 = w38*(B_2_5 + B_2_6);
//                                         const Scalar tmp151 = w37*(-B_1_0 - B_1_5);
//                                         const Scalar tmp152 = w39*(-B_0_0 - B_0_6);
//                                         const Scalar tmp153 = w45*(B_2_1 + B_2_2);
//                                         const Scalar tmp154 = w36*(-B_0_1 - B_0_7);
//                                         const Scalar tmp155 = w40*(-B_1_2 - B_1_7);
//                                         const Scalar tmp156 = w41*(-B_0_3 - B_0_7);
//                                         const Scalar tmp157 = w39*(-B_0_2 - B_0_6);
//                                         const Scalar tmp158 = w36*(-B_0_1 - B_0_5);
//                                         const Scalar tmp159 = w33*(-B_0_0 - B_0_4);
//                                         const Scalar tmp160 = w38*(-B_2_2 - B_2_3);
//                                         const Scalar tmp161 = w35*(-B_2_0 - B_2_1);
//                                         const Scalar tmp162 = w45*(-B_2_4 - B_2_5);
//                                         const Scalar tmp163 = w46*(-B_2_6 - B_2_7);
//                                         const Scalar tmp164 = w38*(-B_2_0 - B_2_3);
//                                         const Scalar tmp165 = w37*(B_1_3 + B_1_6);
//                                         const Scalar tmp166 = w40*(B_1_1 + B_1_4);
//                                         const Scalar tmp167 = w45*(-B_2_4 - B_2_7);
//                                         const Scalar tmp168 = w39*(B_0_3 + B_0_5);
//                                         const Scalar tmp169 = w36*(B_0_2 + B_0_4);
//                                         const Scalar tmp170 = w38*(B_2_1 + B_2_3);
//                                         const Scalar tmp171 = w35*(B_2_0 + B_2_2);
//                                         const Scalar tmp172 = w41*(B_0_5 + B_0_7);
//                                         const Scalar tmp173 = w37*(-B_1_3 - B_1_6);
//                                         const Scalar tmp174 = w39*(B_0_4 + B_0_6);
//                                         const Scalar tmp175 = w45*(B_2_4 + B_2_6);
//                                         const Scalar tmp176 = w36*(B_0_1 + B_0_3);
//                                         const Scalar tmp177 = w40*(-B_1_1 - B_1_4);
//                                         const Scalar tmp178 = w33*(B_0_0 + B_0_2);
//                                         const Scalar tmp179 = w46*(B_2_5 + B_2_7);
//                                         const Scalar tmp180 = w38*(B_2_5 + B_2_7);
//                                         const Scalar tmp181 = w42*(-B_1_3 - B_1_7);
//                                         const Scalar tmp182 = w35*(B_2_4 + B_2_6);
//                                         const Scalar tmp183 = w37*(-B_1_1 - B_1_5);
//                                         const Scalar tmp184 = w39*(B_0_1 + B_0_3 + B_0_5 + B_0_7);
//                                         const Scalar tmp185 = w45*(B_2_0 + B_2_2);
//                                         const Scalar tmp186 = w36*(B_0_0 + B_0_2 + B_0_4 + B_0_6);
//                                         const Scalar tmp187 = w40*(-B_1_2 - B_1_6);
//                                         const Scalar tmp188 = w46*(B_2_1 + B_2_3);
//                                         const Scalar tmp189 = w34*(-B_1_0 - B_1_4);
//                                         const Scalar tmp190 = w38*(B_2_4 + B_2_5);
//                                         const Scalar tmp191 = w35*(B_2_6 + B_2_7);
//                                         const Scalar tmp192 = w41*(-B_0_1 - B_0_5);
//                                         const Scalar tmp193 = w37*(-B_1_0 - B_1_1 - B_1_4 - B_1_5);
//                                         const Scalar tmp194 = w39*(-B_0_0 - B_0_4);
//                                         const Scalar tmp195 = w45*(B_2_2 + B_2_3);
//                                         const Scalar tmp196 = w36*(-B_0_3 - B_0_7);
//                                         const Scalar tmp197 = w40*(-B_1_2 - B_1_3 - B_1_6 - B_1_7);
//                                         const Scalar tmp198 = w33*(-B_0_2 - B_0_6);
//                                         const Scalar tmp199 = w46*(B_2_0 + B_2_1);
//                                         const Scalar tmp200 = w38*(-B_2_6 - B_2_7);
//                                         const Scalar tmp201 = w42*(B_1_2 + B_1_3);
//                                         const Scalar tmp202 = w35*(-B_2_4 - B_2_5);
//                                         const Scalar tmp203 = w37*(B_1_0 + B_1_1);
//                                         const Scalar tmp204 = w45*(-B_2_0 - B_2_1);
//                                         const Scalar tmp205 = w40*(B_1_6 + B_1_7);
//                                         const Scalar tmp206 = w46*(-B_2_2 - B_2_3);
//                                         const Scalar tmp207 = w34*(B_1_4 + B_1_5);
//                                         const Scalar tmp208 = w37*(-B_1_1 - B_1_4);
//                                         const Scalar tmp209 = w39*(-B_0_2 - B_0_4);
//                                         const Scalar tmp210 = w36*(-B_0_3 - B_0_5);
//                                         const Scalar tmp211 = w40*(-B_1_3 - B_1_6);
//                                         const Scalar tmp212 = w38*(B_2_4 + B_2_7);
//                                         const Scalar tmp213 = w45*(B_2_0 + B_2_3);
//                                         const Scalar tmp214 = w41*(B_0_0 + B_0_4);
//                                         const Scalar tmp215 = w39*(B_0_1 + B_0_5);
//                                         const Scalar tmp216 = w36*(B_0_2 + B_0_6);
//                                         const Scalar tmp217 = w33*(B_0_3 + B_0_7);
//                                         const Scalar tmp218 = w42*(B_1_1 + B_1_5);
//                                         const Scalar tmp219 = w37*(B_1_3 + B_1_7);
//                                         const Scalar tmp220 = w40*(B_1_0 + B_1_4);
//                                         const Scalar tmp221 = w34*(B_1_2 + B_1_6);
//                                         const Scalar tmp222 = w39*(-B_0_1 - B_0_7);
//                                         const Scalar tmp223 = w36*(-B_0_0 - B_0_6);
//                                         const Scalar tmp224 = w38*(-B_2_0 - B_2_1);
//                                         const Scalar tmp225 = w35*(-B_2_2 - B_2_3);
//                                         const Scalar tmp226 = w45*(-B_2_6 - B_2_7);
//                                         const Scalar tmp227 = w46*(-B_2_4 - B_2_5);
//                                         const Scalar tmp228 = w38*(B_2_4 + B_2_6);
//                                         const Scalar tmp229 = w42*(B_1_0 + B_1_4);
//                                         const Scalar tmp230 = w35*(B_2_5 + B_2_7);
//                                         const Scalar tmp231 = w37*(B_1_2 + B_1_6);
//                                         const Scalar tmp232 = w39*(-B_0_0 - B_0_2 - B_0_4 - B_0_6);
//                                         const Scalar tmp233 = w45*(B_2_1 + B_2_3);
//                                         const Scalar tmp234 = w36*(-B_0_1 - B_0_3 - B_0_5 - B_0_7);
//                                         const Scalar tmp235 = w40*(B_1_1 + B_1_5);
//                                         const Scalar tmp236 = w46*(B_2_0 + B_2_2);
//                                         const Scalar tmp237 = w34*(B_1_3 + B_1_7);
//                                         const Scalar tmp238 = w42*(-B_1_2 - B_1_6);
//                                         const Scalar tmp239 = w37*(-B_1_0 - B_1_4);
//                                         const Scalar tmp240 = w40*(-B_1_3 - B_1_7);
//                                         const Scalar tmp241 = w34*(-B_1_1 - B_1_5);
//                                         const Scalar tmp242 = w38*(-B_2_4 - B_2_5);
//                                         const Scalar tmp243 = w42*(-B_1_0 - B_1_1);
//                                         const Scalar tmp244 = w35*(-B_2_6 - B_2_7);
//                                         const Scalar tmp245 = w37*(-B_1_2 - B_1_3);
//                                         const Scalar tmp246 = w45*(-B_2_2 - B_2_3);
//                                         const Scalar tmp247 = w40*(-B_1_4 - B_1_5);
//                                         const Scalar tmp248 = w46*(-B_2_0 - B_2_1);
//                                         const Scalar tmp249 = w34*(-B_1_6 - B_1_7);
//                                         const Scalar tmp250 = w42*(B_1_4 + B_1_5);
//                                         const Scalar tmp251 = w37*(B_1_6 + B_1_7);
//                                         const Scalar tmp252 = w40*(B_1_0 + B_1_1);
//                                         const Scalar tmp253 = w34*(B_1_2 + B_1_3);
//                                         const Scalar tmp254 = w38*(-B_2_1 - B_2_3);
//                                         const Scalar tmp255 = w35*(-B_2_0 - B_2_2);
//                                         const Scalar tmp256 = w45*(-B_2_4 - B_2_6);
//                                         const Scalar tmp257 = w46*(-B_2_5 - B_2_7);
//                                         const Scalar tmp258 = w38*(B_2_0 + B_2_1 + B_2_2 + B_2_3);
//                                         const Scalar tmp259 = w45*(B_2_4 + B_2_5 + B_2_6 + B_2_7);
//                                         const Scalar tmp260 = w38*(-B_2_0 - B_2_2);
//                                         const Scalar tmp261 = w35*(-B_2_1 - B_2_3);
//                                         const Scalar tmp262 = w45*(-B_2_5 - B_2_7);
//                                         const Scalar tmp263 = w46*(-B_2_4 - B_2_6);
//                                         EM_S[INDEX4(k,m,0,0,numEq,numComp,8)]+=-B_0_0*w50 - B_0_1*w41 - B_0_6*w33 - B_0_7*w49 + B_1_0*w47 - B_1_2*w42 - B_1_5*w34 + B_1_7*w48 - B_2_0*w43 - B_2_3*w35 - B_2_4*w46 - B_2_7*w44 + tmp132 + tmp137 + tmp208 + tmp209 + tmp210 + tmp211;
//                                         EM_S[INDEX4(k,m,0,1,numEq,numComp,8)]+=-B_0_0*w41 - B_0_1*w50 - B_0_6*w49 - B_0_7*w33 + tmp126 + tmp128 + tmp242 + tmp243 + tmp244 + tmp245 + tmp246 + tmp247 + tmp248 + tmp249;
//                                         EM_S[INDEX4(k,m,0,2,numEq,numComp,8)]+=-B_1_0*w42 + B_1_2*w47 + B_1_5*w48 - B_1_7*w34 + tmp138 + tmp139 + tmp140 + tmp142 + tmp143 + tmp144 + tmp146 + tmp147 + tmp173 + tmp177;
//                                         EM_S[INDEX4(k,m,0,3,numEq,numComp,8)]+=tmp100 + tmp101 + tmp102 + tmp103 + tmp40 + tmp45 + tmp96 + tmp97 + tmp98 + tmp99;
//                                         EM_S[INDEX4(k,m,0,4,numEq,numComp,8)]+=-B_2_0*w46 - B_2_3*w44 - B_2_4*w43 - B_2_7*w35 + tmp114 + tmp115 + tmp116 + tmp117 + tmp118 + tmp119 + tmp120 + tmp121 + tmp92 + tmp93;
//                                         EM_S[INDEX4(k,m,0,5,numEq,numComp,8)]+=tmp192 + tmp193 + tmp194 + tmp196 + tmp197 + tmp198 + tmp224 + tmp225 + tmp226 + tmp227;
//                                         EM_S[INDEX4(k,m,0,6,numEq,numComp,8)]+=tmp232 + tmp234 + tmp238 + tmp239 + tmp240 + tmp241 + tmp260 + tmp261 + tmp262 + tmp263;
//                                         EM_S[INDEX4(k,m,0,7,numEq,numComp,8)]+=tmp60 + tmp61 + tmp62 + tmp63 + tmp64 + tmp65;
//                                         EM_S[INDEX4(k,m,1,0,numEq,numComp,8)]+=B_0_0*w50 + B_0_1*w41 + B_0_6*w33 + B_0_7*w49 + tmp148 + tmp149 + tmp242 + tmp243 + tmp244 + tmp245 + tmp246 + tmp247 + tmp248 + tmp249;
//                                         EM_S[INDEX4(k,m,1,1,numEq,numComp,8)]+=B_0_0*w41 + B_0_1*w50 + B_0_6*w49 + B_0_7*w33 + B_1_1*w47 - B_1_3*w42 - B_1_4*w34 + B_1_6*w48 - B_2_1*w43 - B_2_2*w35 - B_2_5*w46 - B_2_6*w44 + tmp151 + tmp155 + tmp164 + tmp167 + tmp168 + tmp169;
//                                         EM_S[INDEX4(k,m,1,2,numEq,numComp,8)]+=tmp101 + tmp103 + tmp40 + tmp42 + tmp44 + tmp45 + tmp46 + tmp48 + tmp96 + tmp98;
//                                         EM_S[INDEX4(k,m,1,3,numEq,numComp,8)]+=-B_1_1*w42 + B_1_3*w47 + B_1_4*w48 - B_1_6*w34 + tmp20 + tmp21 + tmp22 + tmp23 + tmp24 + tmp25 + tmp26 + tmp27 + tmp28 + tmp29;
//                                         EM_S[INDEX4(k,m,1,4,numEq,numComp,8)]+=tmp193 + tmp197 + tmp214 + tmp215 + tmp216 + tmp217 + tmp224 + tmp225 + tmp226 + tmp227;
//                                         EM_S[INDEX4(k,m,1,5,numEq,numComp,8)]+=-B_2_1*w46 - B_2_2*w44 - B_2_5*w43 - B_2_6*w35 + tmp70 + tmp75 + tmp83 + tmp84 + tmp85 + tmp86 + tmp88 + tmp89 + tmp90 + tmp91;
//                                         EM_S[INDEX4(k,m,1,6,numEq,numComp,8)]+=tmp60 + tmp61 + tmp63 + tmp65 + tmp80 + tmp81;
//                                         EM_S[INDEX4(k,m,1,7,numEq,numComp,8)]+=tmp181 + tmp183 + tmp184 + tmp186 + tmp187 + tmp189 + tmp254 + tmp255 + tmp256 + tmp257;
//                                         EM_S[INDEX4(k,m,2,0,numEq,numComp,8)]+=-B_1_0*w47 + B_1_2*w42 + B_1_5*w34 - B_1_7*w48 + tmp138 + tmp139 + tmp140 + tmp141 + tmp142 + tmp143 + tmp144 + tmp145 + tmp146 + tmp147;
//                                         EM_S[INDEX4(k,m,2,1,numEq,numComp,8)]+=tmp100 + tmp102 + tmp40 + tmp41 + tmp43 + tmp45 + tmp47 + tmp49 + tmp97 + tmp99;
//                                         EM_S[INDEX4(k,m,2,2,numEq,numComp,8)]+=-B_0_2*w50 - B_0_3*w41 - B_0_4*w33 - B_0_5*w49 + B_1_0*w42 - B_1_2*w47 - B_1_5*w48 + B_1_7*w34 - B_2_1*w35 - B_2_2*w43 - B_2_5*w44 - B_2_6*w46 + tmp152 + tmp154 + tmp164 + tmp165 + tmp166 + tmp167;
//                                         EM_S[INDEX4(k,m,2,3,numEq,numComp,8)]+=-B_0_2*w41 - B_0_3*w50 - B_0_4*w49 - B_0_5*w33 + tmp200 + tmp201 + tmp202 + tmp203 + tmp204 + tmp205 + tmp206 + tmp207 + tmp222 + tmp223;
//                                         EM_S[INDEX4(k,m,2,4,numEq,numComp,8)]+=tmp229 + tmp231 + tmp232 + tmp234 + tmp235 + tmp237 + tmp260 + tmp261 + tmp262 + tmp263;
//                                         EM_S[INDEX4(k,m,2,5,numEq,numComp,8)]+=tmp60 + tmp62 + tmp63 + tmp64 + tmp94 + tmp95;
//                                         EM_S[INDEX4(k,m,2,6,numEq,numComp,8)]+=-B_2_1*w44 - B_2_2*w46 - B_2_5*w35 - B_2_6*w43 + tmp70 + tmp71 + tmp72 + tmp73 + tmp74 + tmp75 + tmp76 + tmp77 + tmp78 + tmp79;
//                                         EM_S[INDEX4(k,m,2,7,numEq,numComp,8)]+=tmp107 + tmp111 + tmp156 + tmp157 + tmp158 + tmp159 + tmp160 + tmp161 + tmp162 + tmp163;
//                                         EM_S[INDEX4(k,m,3,0,numEq,numComp,8)]+=tmp40 + tmp41 + tmp42 + tmp43 + tmp44 + tmp45 + tmp46 + tmp47 + tmp48 + tmp49;
//                                         EM_S[INDEX4(k,m,3,1,numEq,numComp,8)]+=-B_1_1*w47 + B_1_3*w42 + B_1_4*w34 - B_1_6*w48 + tmp20 + tmp21 + tmp22 + tmp24 + tmp25 + tmp26 + tmp28 + tmp29 + tmp33 + tmp37;
//                                         EM_S[INDEX4(k,m,3,2,numEq,numComp,8)]+=B_0_2*w50 + B_0_3*w41 + B_0_4*w33 + B_0_5*w49 + tmp200 + tmp201 + tmp202 + tmp203 + tmp204 + tmp205 + tmp206 + tmp207 + tmp54 + tmp56;
//                                         EM_S[INDEX4(k,m,3,3,numEq,numComp,8)]+=B_0_2*w41 + B_0_3*w50 + B_0_4*w49 + B_0_5*w33 + B_1_1*w42 - B_1_3*w47 - B_1_4*w48 + B_1_6*w34 - B_2_0*w35 - B_2_3*w43 - B_2_4*w44 - B_2_7*w46 + tmp132 + tmp133 + tmp134 + tmp135 + tmp136 + tmp137;
//                                         EM_S[INDEX4(k,m,3,4,numEq,numComp,8)]+=tmp60 + tmp63 + tmp80 + tmp81 + tmp94 + tmp95;
//                                         EM_S[INDEX4(k,m,3,5,numEq,numComp,8)]+=tmp184 + tmp186 + tmp218 + tmp219 + tmp220 + tmp221 + tmp254 + tmp255 + tmp256 + tmp257;
//                                         EM_S[INDEX4(k,m,3,6,numEq,numComp,8)]+=tmp106 + tmp107 + tmp108 + tmp110 + tmp111 + tmp112 + tmp160 + tmp161 + tmp162 + tmp163;
//                                         EM_S[INDEX4(k,m,3,7,numEq,numComp,8)]+=-B_2_0*w44 - B_2_3*w46 - B_2_4*w35 - B_2_7*w43 + tmp1 + tmp2 + tmp3 + tmp4 + tmp6 + tmp7 + tmp8 + tmp9 + tmp92 + tmp93;
//                                         EM_S[INDEX4(k,m,4,0,numEq,numComp,8)]+=B_2_0*w43 + B_2_3*w35 + B_2_4*w46 + B_2_7*w44 + tmp0 + tmp114 + tmp115 + tmp116 + tmp117 + tmp118 + tmp119 + tmp120 + tmp121 + tmp5;
//                                         EM_S[INDEX4(k,m,4,1,numEq,numComp,8)]+=tmp190 + tmp191 + tmp192 + tmp193 + tmp194 + tmp195 + tmp196 + tmp197 + tmp198 + tmp199;
//                                         EM_S[INDEX4(k,m,4,2,numEq,numComp,8)]+=tmp228 + tmp230 + tmp232 + tmp233 + tmp234 + tmp236 + tmp238 + tmp239 + tmp240 + tmp241;
//                                         EM_S[INDEX4(k,m,4,3,numEq,numComp,8)]+=tmp258 + tmp259 + tmp61 + tmp62 + tmp64 + tmp65;
//                                         EM_S[INDEX4(k,m,4,4,numEq,numComp,8)]+=-B_0_2*w33 - B_0_3*w49 - B_0_4*w50 - B_0_5*w41 - B_1_1*w34 + B_1_3*w48 + B_1_4*w47 - B_1_6*w42 + B_2_0*w46 + B_2_3*w44 + B_2_4*w43 + B_2_7*w35 + tmp150 + tmp151 + tmp152 + tmp153 + tmp154 + tmp155;
//                                         EM_S[INDEX4(k,m,4,5,numEq,numComp,8)]+=-B_0_2*w49 - B_0_3*w33 - B_0_4*w41 - B_0_5*w50 + tmp222 + tmp223 + tmp50 + tmp51 + tmp52 + tmp53 + tmp55 + tmp57 + tmp58 + tmp59;
//                                         EM_S[INDEX4(k,m,4,6,numEq,numComp,8)]+=B_1_1*w48 - B_1_3*w34 - B_1_4*w42 + B_1_6*w47 + tmp23 + tmp27 + tmp30 + tmp31 + tmp32 + tmp34 + tmp35 + tmp36 + tmp38 + tmp39;
//                                         EM_S[INDEX4(k,m,4,7,numEq,numComp,8)]+=tmp10 + tmp11 + tmp12 + tmp13 + tmp14 + tmp15 + tmp16 + tmp17 + tmp18 + tmp19;
//                                         EM_S[INDEX4(k,m,5,0,numEq,numComp,8)]+=tmp190 + tmp191 + tmp193 + tmp195 + tmp197 + tmp199 + tmp214 + tmp215 + tmp216 + tmp217;
//                                         EM_S[INDEX4(k,m,5,1,numEq,numComp,8)]+=B_2_1*w43 + B_2_2*w35 + B_2_5*w46 + B_2_6*w44 + tmp82 + tmp83 + tmp84 + tmp85 + tmp86 + tmp87 + tmp88 + tmp89 + tmp90 + tmp91;
//                                         EM_S[INDEX4(k,m,5,2,numEq,numComp,8)]+=tmp258 + tmp259 + tmp61 + tmp65 + tmp80 + tmp81;
//                                         EM_S[INDEX4(k,m,5,3,numEq,numComp,8)]+=tmp180 + tmp181 + tmp182 + tmp183 + tmp184 + tmp185 + tmp186 + tmp187 + tmp188 + tmp189;
//                                         EM_S[INDEX4(k,m,5,4,numEq,numComp,8)]+=B_0_2*w33 + B_0_3*w49 + B_0_4*w50 + B_0_5*w41 + tmp50 + tmp51 + tmp52 + tmp53 + tmp54 + tmp55 + tmp56 + tmp57 + tmp58 + tmp59;
//                                         EM_S[INDEX4(k,m,5,5,numEq,numComp,8)]+=B_0_2*w49 + B_0_3*w33 + B_0_4*w41 + B_0_5*w50 - B_1_0*w34 + B_1_2*w48 + B_1_5*w47 - B_1_7*w42 + B_2_1*w46 + B_2_2*w44 + B_2_5*w43 + B_2_6*w35 + tmp134 + tmp135 + tmp208 + tmp211 + tmp212 + tmp213;
//                                         EM_S[INDEX4(k,m,5,6,numEq,numComp,8)]+=tmp10 + tmp11 + tmp13 + tmp15 + tmp17 + tmp19 + tmp66 + tmp67 + tmp68 + tmp69;
//                                         EM_S[INDEX4(k,m,5,7,numEq,numComp,8)]+=B_1_0*w48 - B_1_2*w34 - B_1_5*w42 + B_1_7*w47 + tmp170 + tmp171 + tmp172 + tmp173 + tmp174 + tmp175 + tmp176 + tmp177 + tmp178 + tmp179;
//                                         EM_S[INDEX4(k,m,6,0,numEq,numComp,8)]+=tmp228 + tmp229 + tmp230 + tmp231 + tmp232 + tmp233 + tmp234 + tmp235 + tmp236 + tmp237;
//                                         EM_S[INDEX4(k,m,6,1,numEq,numComp,8)]+=tmp258 + tmp259 + tmp62 + tmp64 + tmp94 + tmp95;
//                                         EM_S[INDEX4(k,m,6,2,numEq,numComp,8)]+=B_2_1*w35 + B_2_2*w43 + B_2_5*w44 + B_2_6*w46 + tmp71 + tmp72 + tmp73 + tmp74 + tmp76 + tmp77 + tmp78 + tmp79 + tmp82 + tmp87;
//                                         EM_S[INDEX4(k,m,6,3,numEq,numComp,8)]+=tmp104 + tmp105 + tmp107 + tmp109 + tmp111 + tmp113 + tmp156 + tmp157 + tmp158 + tmp159;
//                                         EM_S[INDEX4(k,m,6,4,numEq,numComp,8)]+=B_1_1*w34 - B_1_3*w48 - B_1_4*w47 + B_1_6*w42 + tmp30 + tmp31 + tmp32 + tmp33 + tmp34 + tmp35 + tmp36 + tmp37 + tmp38 + tmp39;
//                                         EM_S[INDEX4(k,m,6,5,numEq,numComp,8)]+=tmp10 + tmp12 + tmp14 + tmp15 + tmp16 + tmp18 + tmp250 + tmp251 + tmp252 + tmp253;
//                                         EM_S[INDEX4(k,m,6,6,numEq,numComp,8)]+=-B_0_0*w33 - B_0_1*w49 - B_0_6*w50 - B_0_7*w41 - B_1_1*w48 + B_1_3*w34 + B_1_4*w42 - B_1_6*w47 + B_2_1*w44 + B_2_2*w46 + B_2_5*w35 + B_2_6*w43 + tmp133 + tmp136 + tmp209 + tmp210 + tmp212 + tmp213;
//                                         EM_S[INDEX4(k,m,6,7,numEq,numComp,8)]+=-B_0_0*w49 - B_0_1*w33 - B_0_6*w41 - B_0_7*w50 + tmp122 + tmp123 + tmp124 + tmp125 + tmp126 + tmp127 + tmp128 + tmp129 + tmp130 + tmp131;
//                                         EM_S[INDEX4(k,m,7,0,numEq,numComp,8)]+=tmp258 + tmp259 + tmp80 + tmp81 + tmp94 + tmp95;
//                                         EM_S[INDEX4(k,m,7,1,numEq,numComp,8)]+=tmp180 + tmp182 + tmp184 + tmp185 + tmp186 + tmp188 + tmp218 + tmp219 + tmp220 + tmp221;
//                                         EM_S[INDEX4(k,m,7,2,numEq,numComp,8)]+=tmp104 + tmp105 + tmp106 + tmp107 + tmp108 + tmp109 + tmp110 + tmp111 + tmp112 + tmp113;
//                                         EM_S[INDEX4(k,m,7,3,numEq,numComp,8)]+=B_2_0*w35 + B_2_3*w43 + B_2_4*w44 + B_2_7*w46 + tmp0 + tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6 + tmp7 + tmp8 + tmp9;
//                                         EM_S[INDEX4(k,m,7,4,numEq,numComp,8)]+=tmp10 + tmp15 + tmp250 + tmp251 + tmp252 + tmp253 + tmp66 + tmp67 + tmp68 + tmp69;
//                                         EM_S[INDEX4(k,m,7,5,numEq,numComp,8)]+=B_1_0*w34 - B_1_2*w48 - B_1_5*w47 + B_1_7*w42 + tmp141 + tmp145 + tmp170 + tmp171 + tmp172 + tmp174 + tmp175 + tmp176 + tmp178 + tmp179;
//                                         EM_S[INDEX4(k,m,7,6,numEq,numComp,8)]+=B_0_0*w33 + B_0_1*w49 + B_0_6*w50 + B_0_7*w41 + tmp122 + tmp123 + tmp124 + tmp125 + tmp127 + tmp129 + tmp130 + tmp131 + tmp148 + tmp149;
//                                         EM_S[INDEX4(k,m,7,7,numEq,numComp,8)]+=B_0_0*w49 + B_0_1*w33 + B_0_6*w41 + B_0_7*w50 - B_1_0*w48 + B_1_2*w34 + B_1_5*w42 - B_1_7*w47 + B_2_0*w44 + B_2_3*w46 + B_2_4*w35 + B_2_7*w43 + tmp150 + tmp153 + tmp165 + tmp166 + tmp168 + tmp169;
//                                     }
//                                 }
//                             } else { // constant data
//                                 for (index_t k=0; k<numEq; k++) {
//                                     for (index_t m=0; m<numComp; m++) {
//                                         const Scalar wB0 = B_p[INDEX3(k,0,m,numEq,3)]*w55;
//                                         const Scalar wB1 = B_p[INDEX3(k,1,m,numEq,3)]*w56;
//                                         const Scalar wB2 = B_p[INDEX3(k,2,m,numEq,3)]*w54;
//                                         EM_S[INDEX4(k,m,0,0,numEq,numComp,8)]+= 4.*wB0 + 4.*wB1 + 4.*wB2;
//                                         EM_S[INDEX4(k,m,0,1,numEq,numComp,8)]+= 4.*wB0 + 2.*wB1 + 2.*wB2;
//                                         EM_S[INDEX4(k,m,0,2,numEq,numComp,8)]+= 2.*wB0 + 4.*wB1 + 2.*wB2;
//                                         EM_S[INDEX4(k,m,0,3,numEq,numComp,8)]+= 2.*wB0 + 2.*wB1 +    wB2;
//                                         EM_S[INDEX4(k,m,0,4,numEq,numComp,8)]+= 2.*wB0 + 2.*wB1 + 4.*wB2;
//                                         EM_S[INDEX4(k,m,0,5,numEq,numComp,8)]+= 2.*wB0 +    wB1 + 2.*wB2;
//                                         EM_S[INDEX4(k,m,0,6,numEq,numComp,8)]+=    wB0 + 2.*wB1 + 2.*wB2;
//                                         EM_S[INDEX4(k,m,0,7,numEq,numComp,8)]+=    wB0 +    wB1 +    wB2;
//                                         EM_S[INDEX4(k,m,1,0,numEq,numComp,8)]+=-4.*wB0 + 2.*wB1 + 2.*wB2;
//                                         EM_S[INDEX4(k,m,1,1,numEq,numComp,8)]+=-4.*wB0 + 4.*wB1 + 4.*wB2;
//                                         EM_S[INDEX4(k,m,1,2,numEq,numComp,8)]+=-2.*wB0 + 2.*wB1 +    wB2;
//                                         EM_S[INDEX4(k,m,1,3,numEq,numComp,8)]+=-2.*wB0 + 4.*wB1 + 2.*wB2;
//                                         EM_S[INDEX4(k,m,1,4,numEq,numComp,8)]+=-2.*wB0 +    wB1 + 2.*wB2;
//                                         EM_S[INDEX4(k,m,1,5,numEq,numComp,8)]+=-2.*wB0 + 2.*wB1 + 4.*wB2;
//                                         EM_S[INDEX4(k,m,1,6,numEq,numComp,8)]+=   -wB0 +    wB1 +    wB2;
//                                         EM_S[INDEX4(k,m,1,7,numEq,numComp,8)]+=   -wB0 + 2.*wB1 + 2.*wB2;
//                                         EM_S[INDEX4(k,m,2,0,numEq,numComp,8)]+= 2.*wB0 - 4.*wB1 + 2.*wB2;
//                                         EM_S[INDEX4(k,m,2,1,numEq,numComp,8)]+= 2.*wB0 - 2.*wB1 +    wB2;
//                                         EM_S[INDEX4(k,m,2,2,numEq,numComp,8)]+= 4.*wB0 - 4.*wB1 + 4.*wB2;
//                                         EM_S[INDEX4(k,m,2,3,numEq,numComp,8)]+= 4.*wB0 - 2.*wB1 + 2.*wB2;
//                                         EM_S[INDEX4(k,m,2,4,numEq,numComp,8)]+=    wB0 - 2.*wB1 + 2.*wB2;
//                                         EM_S[INDEX4(k,m,2,5,numEq,numComp,8)]+=    wB0 -    wB1 +    wB2;
//                                         EM_S[INDEX4(k,m,2,6,numEq,numComp,8)]+= 2.*wB0 - 2.*wB1 + 4.*wB2;
//                                         EM_S[INDEX4(k,m,2,7,numEq,numComp,8)]+= 2.*wB0 -    wB1 + 2.*wB2;
//                                         EM_S[INDEX4(k,m,3,0,numEq,numComp,8)]+=-2.*wB0 - 2.*wB1 +    wB2;
//                                         EM_S[INDEX4(k,m,3,1,numEq,numComp,8)]+=-2.*wB0 - 4.*wB1 + 2.*wB2;
//                                         EM_S[INDEX4(k,m,3,2,numEq,numComp,8)]+=-4.*wB0 - 2.*wB1 + 2.*wB2;
//                                         EM_S[INDEX4(k,m,3,3,numEq,numComp,8)]+=-4.*wB0 - 4.*wB1 + 4.*wB2;
//                                         EM_S[INDEX4(k,m,3,4,numEq,numComp,8)]+=   -wB0 -    wB1 +    wB2;
//                                         EM_S[INDEX4(k,m,3,5,numEq,numComp,8)]+=   -wB0 - 2.*wB1 + 2.*wB2;
//                                         EM_S[INDEX4(k,m,3,6,numEq,numComp,8)]+=-2.*wB0 -    wB1 + 2.*wB2;
//                                         EM_S[INDEX4(k,m,3,7,numEq,numComp,8)]+=-2.*wB0 - 2.*wB1 + 4.*wB2;
//                                         EM_S[INDEX4(k,m,4,0,numEq,numComp,8)]+= 2.*wB0 + 2.*wB1 - 4.*wB2;
//                                         EM_S[INDEX4(k,m,4,1,numEq,numComp,8)]+= 2.*wB0 +    wB1 - 2.*wB2;
//                                         EM_S[INDEX4(k,m,4,2,numEq,numComp,8)]+=    wB0 + 2.*wB1 - 2.*wB2;
//                                         EM_S[INDEX4(k,m,4,3,numEq,numComp,8)]+=    wB0 +    wB1 -    wB2;
//                                         EM_S[INDEX4(k,m,4,4,numEq,numComp,8)]+= 4.*wB0 + 4.*wB1 - 4.*wB2;
//                                         EM_S[INDEX4(k,m,4,5,numEq,numComp,8)]+= 4.*wB0 + 2.*wB1 - 2.*wB2;
//                                         EM_S[INDEX4(k,m,4,6,numEq,numComp,8)]+= 2.*wB0 + 4.*wB1 - 2.*wB2;
//                                         EM_S[INDEX4(k,m,4,7,numEq,numComp,8)]+= 2.*wB0 + 2.*wB1 -    wB2;
//                                         EM_S[INDEX4(k,m,5,0,numEq,numComp,8)]+=-2.*wB0 +    wB1 - 2.*wB2;
//                                         EM_S[INDEX4(k,m,5,1,numEq,numComp,8)]+=-2.*wB0 + 2.*wB1 - 4.*wB2;
//                                         EM_S[INDEX4(k,m,5,2,numEq,numComp,8)]+=   -wB0 +    wB1 -    wB2;
//                                         EM_S[INDEX4(k,m,5,3,numEq,numComp,8)]+=   -wB0 + 2.*wB1 - 2.*wB2;
//                                         EM_S[INDEX4(k,m,5,4,numEq,numComp,8)]+=-4.*wB0 + 2.*wB1 - 2.*wB2;
//                                         EM_S[INDEX4(k,m,5,5,numEq,numComp,8)]+=-4.*wB0 + 4.*wB1 - 4.*wB2;
//                                         EM_S[INDEX4(k,m,5,6,numEq,numComp,8)]+=-2.*wB0 + 2.*wB1 -    wB2;
//                                         EM_S[INDEX4(k,m,5,7,numEq,numComp,8)]+=-2.*wB0 + 4.*wB1 - 2.*wB2;
//                                         EM_S[INDEX4(k,m,6,0,numEq,numComp,8)]+=    wB0 - 2.*wB1 - 2.*wB2;
//                                         EM_S[INDEX4(k,m,6,1,numEq,numComp,8)]+=    wB0 -    wB1 -    wB2;
//                                         EM_S[INDEX4(k,m,6,2,numEq,numComp,8)]+= 2.*wB0 - 2.*wB1 - 4.*wB2;
//                                         EM_S[INDEX4(k,m,6,3,numEq,numComp,8)]+= 2.*wB0 -    wB1 - 2.*wB2;
//                                         EM_S[INDEX4(k,m,6,4,numEq,numComp,8)]+= 2.*wB0 - 4.*wB1 - 2.*wB2;
//                                         EM_S[INDEX4(k,m,6,5,numEq,numComp,8)]+= 2.*wB0 - 2.*wB1 -    wB2;
//                                         EM_S[INDEX4(k,m,6,6,numEq,numComp,8)]+= 4.*wB0 - 4.*wB1 - 4.*wB2;
//                                         EM_S[INDEX4(k,m,6,7,numEq,numComp,8)]+= 4.*wB0 - 2.*wB1 - 2.*wB2;
//                                         EM_S[INDEX4(k,m,7,0,numEq,numComp,8)]+=   -wB0 -    wB1 -    wB2;
//                                         EM_S[INDEX4(k,m,7,1,numEq,numComp,8)]+=   -wB0 - 2.*wB1 - 2.*wB2;
//                                         EM_S[INDEX4(k,m,7,2,numEq,numComp,8)]+=-2.*wB0 -    wB1 - 2.*wB2;
//                                         EM_S[INDEX4(k,m,7,3,numEq,numComp,8)]+=-2.*wB0 - 2.*wB1 - 4.*wB2;
//                                         EM_S[INDEX4(k,m,7,4,numEq,numComp,8)]+=-2.*wB0 - 2.*wB1 -    wB2;
//                                         EM_S[INDEX4(k,m,7,5,numEq,numComp,8)]+=-2.*wB0 - 4.*wB1 - 2.*wB2;
//                                         EM_S[INDEX4(k,m,7,6,numEq,numComp,8)]+=-4.*wB0 - 2.*wB1 - 2.*wB2;
//                                         EM_S[INDEX4(k,m,7,7,numEq,numComp,8)]+=-4.*wB0 - 4.*wB1 - 4.*wB2;
//                                     }
//                                 }
//                             }
//                         }
//                         ///////////////
//                         // process C //
//                         ///////////////
//                         if (!C.isEmpty()) {
//                             const Scalar* C_p = C.getSampleDataRO(id, zero);
//                             if (C.actsExpanded()) {
//                                 for (index_t k=0; k<numEq; k++) {
//                                     for (index_t m=0; m<numComp; m++) {
//                                         const Scalar C_0_0 = C_p[INDEX4(k,m,0, 0, numEq,numComp,3)];
//                                         const Scalar C_1_0 = C_p[INDEX4(k,m,1, 0, numEq,numComp,3)];
//                                         const Scalar C_2_0 = C_p[INDEX4(k,m,2, 0, numEq,numComp,3)];
//                                         const Scalar C_0_1 = C_p[INDEX4(k,m,0, 1, numEq,numComp,3)];
//                                         const Scalar C_1_1 = C_p[INDEX4(k,m,1, 1, numEq,numComp,3)];
//                                         const Scalar C_2_1 = C_p[INDEX4(k,m,2, 1, numEq,numComp,3)];
//                                         const Scalar C_0_2 = C_p[INDEX4(k,m,0, 2, numEq,numComp,3)];
//                                         const Scalar C_1_2 = C_p[INDEX4(k,m,1, 2, numEq,numComp,3)];
//                                         const Scalar C_2_2 = C_p[INDEX4(k,m,2, 2, numEq,numComp,3)];
//                                         const Scalar C_0_3 = C_p[INDEX4(k,m,0, 3, numEq,numComp,3)];
//                                         const Scalar C_1_3 = C_p[INDEX4(k,m,1, 3, numEq,numComp,3)];
//                                         const Scalar C_2_3 = C_p[INDEX4(k,m,2, 3, numEq,numComp,3)];
//                                         const Scalar C_0_4 = C_p[INDEX4(k,m,0, 4, numEq,numComp,3)];
//                                         const Scalar C_1_4 = C_p[INDEX4(k,m,1, 4, numEq,numComp,3)];
//                                         const Scalar C_2_4 = C_p[INDEX4(k,m,2, 4, numEq,numComp,3)];
//                                         const Scalar C_0_5 = C_p[INDEX4(k,m,0, 5, numEq,numComp,3)];
//                                         const Scalar C_1_5 = C_p[INDEX4(k,m,1, 5, numEq,numComp,3)];
//                                         const Scalar C_2_5 = C_p[INDEX4(k,m,2, 5, numEq,numComp,3)];
//                                         const Scalar C_0_6 = C_p[INDEX4(k,m,0, 6, numEq,numComp,3)];
//                                         const Scalar C_1_6 = C_p[INDEX4(k,m,1, 6, numEq,numComp,3)];
//                                         const Scalar C_2_6 = C_p[INDEX4(k,m,2, 6, numEq,numComp,3)];
//                                         const Scalar C_0_7 = C_p[INDEX4(k,m,0, 7, numEq,numComp,3)];
//                                         const Scalar C_1_7 = C_p[INDEX4(k,m,1, 7, numEq,numComp,3)];
//                                         const Scalar C_2_7 = C_p[INDEX4(k,m,2, 7, numEq,numComp,3)];
//                                         const Scalar tmp0 = w38*(-C_2_5 - C_2_6);
//                                         const Scalar tmp1 = w42*(C_1_3 + C_1_7);
//                                         const Scalar tmp2 = w41*(C_0_3 + C_0_7);
//                                         const Scalar tmp3 = w37*(C_1_1 + C_1_5);
//                                         const Scalar tmp4 = w39*(C_0_2 + C_0_6);
//                                         const Scalar tmp5 = w45*(-C_2_1 - C_2_2);
//                                         const Scalar tmp6 = w36*(C_0_1 + C_0_5);
//                                         const Scalar tmp7 = w40*(C_1_2 + C_1_6);
//                                         const Scalar tmp8 = w33*(C_0_0 + C_0_4);
//                                         const Scalar tmp9 = w34*(C_1_0 + C_1_4);
//                                         const Scalar tmp10 = w38*(C_2_4 + C_2_5 + C_2_6 + C_2_7);
//                                         const Scalar tmp11 = w42*(C_1_4 + C_1_5);
//                                         const Scalar tmp12 = w41*(C_0_4 + C_0_6);
//                                         const Scalar tmp13 = w37*(C_1_6 + C_1_7);
//                                         const Scalar tmp14 = w39*(C_0_5 + C_0_7);
//                                         const Scalar tmp15 = w45*(C_2_0 + C_2_1 + C_2_2 + C_2_3);
//                                         const Scalar tmp16 = w36*(C_0_0 + C_0_2);
//                                         const Scalar tmp17 = w40*(C_1_0 + C_1_1);
//                                         const Scalar tmp18 = w33*(C_0_1 + C_0_3);
//                                         const Scalar tmp19 = w34*(C_1_2 + C_1_3);
//                                         const Scalar tmp20 = w38*(-C_2_5 - C_2_7);
//                                         const Scalar tmp21 = w35*(-C_2_4 - C_2_6);
//                                         const Scalar tmp22 = w41*(C_0_1 + C_0_3);
//                                         const Scalar tmp23 = w37*(C_1_0 + C_1_5);
//                                         const Scalar tmp24 = w39*(C_0_0 + C_0_2);
//                                         const Scalar tmp25 = w45*(-C_2_0 - C_2_2);
//                                         const Scalar tmp26 = w36*(C_0_5 + C_0_7);
//                                         const Scalar tmp27 = w40*(C_1_2 + C_1_7);
//                                         const Scalar tmp28 = w33*(C_0_4 + C_0_6);
//                                         const Scalar tmp29 = w46*(-C_2_1 - C_2_3);
//                                         const Scalar tmp30 = w38*(C_2_0 + C_2_2);
//                                         const Scalar tmp31 = w35*(C_2_1 + C_2_3);
//                                         const Scalar tmp32 = w41*(-C_0_4 - C_0_6);
//                                         const Scalar tmp33 = w37*(-C_1_2 - C_1_7);
//                                         const Scalar tmp34 = w39*(-C_0_5 - C_0_7);
//                                         const Scalar tmp35 = w45*(C_2_5 + C_2_7);
//                                         const Scalar tmp36 = w36*(-C_0_0 - C_0_2);
//                                         const Scalar tmp37 = w40*(-C_1_0 - C_1_5);
//                                         const Scalar tmp38 = w33*(-C_0_1 - C_0_3);
//                                         const Scalar tmp39 = w46*(C_2_4 + C_2_6);
//                                         const Scalar tmp40 = w38*(-C_2_0 - C_2_1 - C_2_2 - C_2_3);
//                                         const Scalar tmp41 = w42*(-C_1_2 - C_1_3);
//                                         const Scalar tmp42 = w41*(-C_0_1 - C_0_3);
//                                         const Scalar tmp43 = w37*(-C_1_0 - C_1_1);
//                                         const Scalar tmp44 = w39*(-C_0_0 - C_0_2);
//                                         const Scalar tmp45 = w45*(-C_2_4 - C_2_5 - C_2_6 - C_2_7);
//                                         const Scalar tmp46 = w36*(-C_0_5 - C_0_7);
//                                         const Scalar tmp47 = w40*(-C_1_6 - C_1_7);
//                                         const Scalar tmp48 = w33*(-C_0_4 - C_0_6);
//                                         const Scalar tmp49 = w34*(-C_1_4 - C_1_5);
//                                         const Scalar tmp50 = w38*(C_2_0 + C_2_1);
//                                         const Scalar tmp51 = w42*(-C_1_4 - C_1_5);
//                                         const Scalar tmp52 = w35*(C_2_2 + C_2_3);
//                                         const Scalar tmp53 = w37*(-C_1_6 - C_1_7);
//                                         const Scalar tmp54 = w39*(-C_0_1 - C_0_7);
//                                         const Scalar tmp55 = w45*(C_2_6 + C_2_7);
//                                         const Scalar tmp56 = w36*(-C_0_0 - C_0_6);
//                                         const Scalar tmp57 = w40*(-C_1_0 - C_1_1);
//                                         const Scalar tmp58 = w46*(C_2_4 + C_2_5);
//                                         const Scalar tmp59 = w34*(-C_1_2 - C_1_3);
//                                         const Scalar tmp60 = w38*(C_2_0 + C_2_1 + C_2_2 + C_2_3);
//                                         const Scalar tmp61 = w37*(C_1_0 + C_1_1 + C_1_4 + C_1_5);
//                                         const Scalar tmp62 = w39*(C_0_0 + C_0_2 + C_0_4 + C_0_6);
//                                         const Scalar tmp63 = w45*(C_2_4 + C_2_5 + C_2_6 + C_2_7);
//                                         const Scalar tmp64 = w36*(C_0_1 + C_0_3 + C_0_5 + C_0_7);
//                                         const Scalar tmp65 = w40*(C_1_2 + C_1_3 + C_1_6 + C_1_7);
//                                         const Scalar tmp66 = w41*(-C_0_5 - C_0_7);
//                                         const Scalar tmp67 = w39*(-C_0_4 - C_0_6);
//                                         const Scalar tmp68 = w36*(-C_0_1 - C_0_3);
//                                         const Scalar tmp69 = w33*(-C_0_0 - C_0_2);
//                                         const Scalar tmp70 = w38*(C_2_0 + C_2_3);
//                                         const Scalar tmp71 = w42*(C_1_2 + C_1_6);
//                                         const Scalar tmp72 = w41*(-C_0_2 - C_0_6);
//                                         const Scalar tmp73 = w37*(C_1_0 + C_1_4);
//                                         const Scalar tmp74 = w39*(-C_0_3 - C_0_7);
//                                         const Scalar tmp75 = w45*(C_2_4 + C_2_7);
//                                         const Scalar tmp76 = w36*(-C_0_0 - C_0_4);
//                                         const Scalar tmp77 = w40*(C_1_3 + C_1_7);
//                                         const Scalar tmp78 = w33*(-C_0_1 - C_0_5);
//                                         const Scalar tmp79 = w34*(C_1_1 + C_1_5);
//                                         const Scalar tmp80 = w39*(-C_0_1 - C_0_3 - C_0_5 - C_0_7);
//                                         const Scalar tmp81 = w36*(-C_0_0 - C_0_2 - C_0_4 - C_0_6);
//                                         const Scalar tmp82 = w38*(-C_2_4 - C_2_7);
//                                         const Scalar tmp83 = w42*(-C_1_1 - C_1_5);
//                                         const Scalar tmp84 = w41*(C_0_1 + C_0_5);
//                                         const Scalar tmp85 = w37*(-C_1_3 - C_1_7);
//                                         const Scalar tmp86 = w39*(C_0_0 + C_0_4);
//                                         const Scalar tmp87 = w45*(-C_2_0 - C_2_3);
//                                         const Scalar tmp88 = w36*(C_0_3 + C_0_7);
//                                         const Scalar tmp89 = w40*(-C_1_0 - C_1_4);
//                                         const Scalar tmp90 = w33*(C_0_2 + C_0_6);
//                                         const Scalar tmp91 = w34*(-C_1_2 - C_1_6);
//                                         const Scalar tmp92 = w38*(C_2_1 + C_2_2);
//                                         const Scalar tmp93 = w45*(C_2_5 + C_2_6);
//                                         const Scalar tmp94 = w37*(-C_1_2 - C_1_3 - C_1_6 - C_1_7);
//                                         const Scalar tmp95 = w40*(-C_1_0 - C_1_1 - C_1_4 - C_1_5);
//                                         const Scalar tmp96 = w42*(C_1_0 + C_1_1);
//                                         const Scalar tmp97 = w41*(C_0_0 + C_0_2);
//                                         const Scalar tmp98 = w37*(C_1_2 + C_1_3);
//                                         const Scalar tmp99 = w39*(C_0_1 + C_0_3);
//                                         const Scalar tmp100 = w36*(C_0_4 + C_0_6);
//                                         const Scalar tmp101 = w40*(C_1_4 + C_1_5);
//                                         const Scalar tmp102 = w33*(C_0_5 + C_0_7);
//                                         const Scalar tmp103 = w34*(C_1_6 + C_1_7);
//                                         const Scalar tmp104 = w38*(-C_2_2 - C_2_3);
//                                         const Scalar tmp105 = w35*(-C_2_0 - C_2_1);
//                                         const Scalar tmp106 = w41*(-C_0_3 - C_0_7);
//                                         const Scalar tmp107 = w37*(C_1_2 + C_1_3 + C_1_6 + C_1_7);
//                                         const Scalar tmp108 = w39*(-C_0_2 - C_0_6);
//                                         const Scalar tmp109 = w45*(-C_2_4 - C_2_5);
//                                         const Scalar tmp110 = w36*(-C_0_1 - C_0_5);
//                                         const Scalar tmp111 = w40*(C_1_0 + C_1_1 + C_1_4 + C_1_5);
//                                         const Scalar tmp112 = w33*(-C_0_0 - C_0_4);
//                                         const Scalar tmp113 = w46*(-C_2_6 - C_2_7);
//                                         const Scalar tmp114 = w42*(-C_1_0 - C_1_4);
//                                         const Scalar tmp115 = w41*(-C_0_0 - C_0_4);
//                                         const Scalar tmp116 = w37*(-C_1_2 - C_1_6);
//                                         const Scalar tmp117 = w39*(-C_0_1 - C_0_5);
//                                         const Scalar tmp118 = w36*(-C_0_2 - C_0_6);
//                                         const Scalar tmp119 = w40*(-C_1_1 - C_1_5);
//                                         const Scalar tmp120 = w33*(-C_0_3 - C_0_7);
//                                         const Scalar tmp121 = w34*(-C_1_3 - C_1_7);
//                                         const Scalar tmp122 = w38*(C_2_2 + C_2_3);
//                                         const Scalar tmp123 = w42*(C_1_6 + C_1_7);
//                                         const Scalar tmp124 = w35*(C_2_0 + C_2_1);
//                                         const Scalar tmp125 = w37*(C_1_4 + C_1_5);
//                                         const Scalar tmp126 = w39*(C_0_2 + C_0_4);
//                                         const Scalar tmp127 = w45*(C_2_4 + C_2_5);
//                                         const Scalar tmp128 = w36*(C_0_3 + C_0_5);
//                                         const Scalar tmp129 = w40*(C_1_2 + C_1_3);
//                                         const Scalar tmp130 = w46*(C_2_6 + C_2_7);
//                                         const Scalar tmp131 = w34*(C_1_0 + C_1_1);
//                                         const Scalar tmp132 = w38*(-C_2_1 - C_2_2);
//                                         const Scalar tmp133 = w37*(C_1_2 + C_1_7);
//                                         const Scalar tmp134 = w39*(C_0_1 + C_0_7);
//                                         const Scalar tmp135 = w36*(C_0_0 + C_0_6);
//                                         const Scalar tmp136 = w40*(C_1_0 + C_1_5);
//                                         const Scalar tmp137 = w45*(-C_2_5 - C_2_6);
//                                         const Scalar tmp138 = w38*(-C_2_4 - C_2_6);
//                                         const Scalar tmp139 = w35*(-C_2_5 - C_2_7);
//                                         const Scalar tmp140 = w41*(-C_0_0 - C_0_2);
//                                         const Scalar tmp141 = w37*(-C_1_3 - C_1_6);
//                                         const Scalar tmp142 = w39*(-C_0_1 - C_0_3);
//                                         const Scalar tmp143 = w45*(-C_2_1 - C_2_3);
//                                         const Scalar tmp144 = w36*(-C_0_4 - C_0_6);
//                                         const Scalar tmp145 = w40*(-C_1_1 - C_1_4);
//                                         const Scalar tmp146 = w33*(-C_0_5 - C_0_7);
//                                         const Scalar tmp147 = w46*(-C_2_0 - C_2_2);
//                                         const Scalar tmp148 = w39*(-C_0_3 - C_0_5);
//                                         const Scalar tmp149 = w36*(-C_0_2 - C_0_4);
//                                         const Scalar tmp150 = w38*(C_2_5 + C_2_6);
//                                         const Scalar tmp151 = w37*(-C_1_0 - C_1_5);
//                                         const Scalar tmp152 = w39*(-C_0_0 - C_0_6);
//                                         const Scalar tmp153 = w45*(C_2_1 + C_2_2);
//                                         const Scalar tmp154 = w36*(-C_0_1 - C_0_7);
//                                         const Scalar tmp155 = w40*(-C_1_2 - C_1_7);
//                                         const Scalar tmp156 = w41*(C_0_2 + C_0_6);
//                                         const Scalar tmp157 = w39*(C_0_3 + C_0_7);
//                                         const Scalar tmp158 = w36*(C_0_0 + C_0_4);
//                                         const Scalar tmp159 = w33*(C_0_1 + C_0_5);
//                                         const Scalar tmp160 = w38*(C_2_6 + C_2_7);
//                                         const Scalar tmp161 = w35*(C_2_4 + C_2_5);
//                                         const Scalar tmp162 = w45*(C_2_0 + C_2_1);
//                                         const Scalar tmp163 = w46*(C_2_2 + C_2_3);
//                                         const Scalar tmp164 = w38*(-C_2_0 - C_2_3);
//                                         const Scalar tmp165 = w37*(C_1_3 + C_1_6);
//                                         const Scalar tmp166 = w40*(C_1_1 + C_1_4);
//                                         const Scalar tmp167 = w45*(-C_2_4 - C_2_7);
//                                         const Scalar tmp168 = w39*(C_0_3 + C_0_5);
//                                         const Scalar tmp169 = w36*(C_0_2 + C_0_4);
//                                         const Scalar tmp170 = w38*(C_2_1 + C_2_3);
//                                         const Scalar tmp171 = w35*(C_2_0 + C_2_2);
//                                         const Scalar tmp172 = w41*(C_0_5 + C_0_7);
//                                         const Scalar tmp173 = w37*(C_1_1 + C_1_4);
//                                         const Scalar tmp174 = w39*(C_0_4 + C_0_6);
//                                         const Scalar tmp175 = w45*(C_2_4 + C_2_6);
//                                         const Scalar tmp176 = w36*(C_0_1 + C_0_3);
//                                         const Scalar tmp177 = w40*(C_1_3 + C_1_6);
//                                         const Scalar tmp178 = w33*(C_0_0 + C_0_2);
//                                         const Scalar tmp179 = w46*(C_2_5 + C_2_7);
//                                         const Scalar tmp180 = w38*(-C_2_1 - C_2_3);
//                                         const Scalar tmp181 = w42*(C_1_1 + C_1_5);
//                                         const Scalar tmp182 = w35*(-C_2_0 - C_2_2);
//                                         const Scalar tmp183 = w37*(C_1_3 + C_1_7);
//                                         const Scalar tmp184 = w39*(C_0_1 + C_0_3 + C_0_5 + C_0_7);
//                                         const Scalar tmp185 = w45*(-C_2_4 - C_2_6);
//                                         const Scalar tmp186 = w36*(C_0_0 + C_0_2 + C_0_4 + C_0_6);
//                                         const Scalar tmp187 = w40*(C_1_0 + C_1_4);
//                                         const Scalar tmp188 = w46*(-C_2_5 - C_2_7);
//                                         const Scalar tmp189 = w34*(C_1_2 + C_1_6);
//                                         const Scalar tmp190 = w38*(-C_2_0 - C_2_1);
//                                         const Scalar tmp191 = w35*(-C_2_2 - C_2_3);
//                                         const Scalar tmp192 = w41*(C_0_0 + C_0_4);
//                                         const Scalar tmp193 = w37*(-C_1_0 - C_1_1 - C_1_4 - C_1_5);
//                                         const Scalar tmp194 = w39*(C_0_1 + C_0_5);
//                                         const Scalar tmp195 = w45*(-C_2_6 - C_2_7);
//                                         const Scalar tmp196 = w36*(C_0_2 + C_0_6);
//                                         const Scalar tmp197 = w40*(-C_1_2 - C_1_3 - C_1_6 - C_1_7);
//                                         const Scalar tmp198 = w33*(C_0_3 + C_0_7);
//                                         const Scalar tmp199 = w46*(-C_2_4 - C_2_5);
//                                         const Scalar tmp200 = w38*(-C_2_6 - C_2_7);
//                                         const Scalar tmp201 = w42*(C_1_2 + C_1_3);
//                                         const Scalar tmp202 = w35*(-C_2_4 - C_2_5);
//                                         const Scalar tmp203 = w37*(C_1_0 + C_1_1);
//                                         const Scalar tmp204 = w45*(-C_2_0 - C_2_1);
//                                         const Scalar tmp205 = w40*(C_1_6 + C_1_7);
//                                         const Scalar tmp206 = w46*(-C_2_2 - C_2_3);
//                                         const Scalar tmp207 = w34*(C_1_4 + C_1_5);
//                                         const Scalar tmp208 = w37*(-C_1_1 - C_1_4);
//                                         const Scalar tmp209 = w39*(-C_0_2 - C_0_4);
//                                         const Scalar tmp210 = w36*(-C_0_3 - C_0_5);
//                                         const Scalar tmp211 = w40*(-C_1_3 - C_1_6);
//                                         const Scalar tmp212 = w38*(C_2_4 + C_2_7);
//                                         const Scalar tmp213 = w45*(C_2_0 + C_2_3);
//                                         const Scalar tmp214 = w41*(-C_0_1 - C_0_5);
//                                         const Scalar tmp215 = w39*(-C_0_0 - C_0_4);
//                                         const Scalar tmp216 = w36*(-C_0_3 - C_0_7);
//                                         const Scalar tmp217 = w33*(-C_0_2 - C_0_6);
//                                         const Scalar tmp218 = w42*(-C_1_3 - C_1_7);
//                                         const Scalar tmp219 = w37*(-C_1_1 - C_1_5);
//                                         const Scalar tmp220 = w40*(-C_1_2 - C_1_6);
//                                         const Scalar tmp221 = w34*(-C_1_0 - C_1_4);
//                                         const Scalar tmp222 = w39*(C_0_0 + C_0_6);
//                                         const Scalar tmp223 = w36*(C_0_1 + C_0_7);
//                                         const Scalar tmp224 = w38*(C_2_4 + C_2_5);
//                                         const Scalar tmp225 = w35*(C_2_6 + C_2_7);
//                                         const Scalar tmp226 = w45*(C_2_2 + C_2_3);
//                                         const Scalar tmp227 = w46*(C_2_0 + C_2_1);
//                                         const Scalar tmp228 = w38*(-C_2_0 - C_2_2);
//                                         const Scalar tmp229 = w42*(-C_1_2 - C_1_6);
//                                         const Scalar tmp230 = w35*(-C_2_1 - C_2_3);
//                                         const Scalar tmp231 = w37*(-C_1_0 - C_1_4);
//                                         const Scalar tmp232 = w39*(-C_0_0 - C_0_2 - C_0_4 - C_0_6);
//                                         const Scalar tmp233 = w45*(-C_2_5 - C_2_7);
//                                         const Scalar tmp234 = w36*(-C_0_1 - C_0_3 - C_0_5 - C_0_7);
//                                         const Scalar tmp235 = w40*(-C_1_3 - C_1_7);
//                                         const Scalar tmp236 = w46*(-C_2_4 - C_2_6);
//                                         const Scalar tmp237 = w34*(-C_1_1 - C_1_5);
//                                         const Scalar tmp238 = w42*(C_1_0 + C_1_4);
//                                         const Scalar tmp239 = w37*(C_1_2 + C_1_6);
//                                         const Scalar tmp240 = w40*(C_1_1 + C_1_5);
//                                         const Scalar tmp241 = w34*(C_1_3 + C_1_7);
//                                         const Scalar tmp242 = w38*(-C_2_4 - C_2_5);
//                                         const Scalar tmp243 = w42*(-C_1_0 - C_1_1);
//                                         const Scalar tmp244 = w35*(-C_2_6 - C_2_7);
//                                         const Scalar tmp245 = w37*(-C_1_2 - C_1_3);
//                                         const Scalar tmp246 = w45*(-C_2_2 - C_2_3);
//                                         const Scalar tmp247 = w40*(-C_1_4 - C_1_5);
//                                         const Scalar tmp248 = w46*(-C_2_0 - C_2_1);
//                                         const Scalar tmp249 = w34*(-C_1_6 - C_1_7);
//                                         const Scalar tmp250 = w42*(-C_1_6 - C_1_7);
//                                         const Scalar tmp251 = w37*(-C_1_4 - C_1_5);
//                                         const Scalar tmp252 = w40*(-C_1_2 - C_1_3);
//                                         const Scalar tmp253 = w34*(-C_1_0 - C_1_1);
//                                         const Scalar tmp254 = w38*(C_2_5 + C_2_7);
//                                         const Scalar tmp255 = w35*(C_2_4 + C_2_6);
//                                         const Scalar tmp256 = w45*(C_2_0 + C_2_2);
//                                         const Scalar tmp257 = w46*(C_2_1 + C_2_3);
//                                         const Scalar tmp258 = w38*(-C_2_4 - C_2_5 - C_2_6 - C_2_7);
//                                         const Scalar tmp259 = w45*(-C_2_0 - C_2_1 - C_2_2 - C_2_3);
//                                         const Scalar tmp260 = w38*(C_2_4 + C_2_6);
//                                         const Scalar tmp261 = w35*(C_2_5 + C_2_7);
//                                         const Scalar tmp262 = w45*(C_2_1 + C_2_3);
//                                         const Scalar tmp263 = w46*(C_2_0 + C_2_2);
//                                         EM_S[INDEX4(k,m,0,0,numEq,numComp,8)]+=-C_0_0*w50 - C_0_1*w41 - C_0_6*w33 - C_0_7*w49 + C_1_0*w47 - C_1_2*w42 - C_1_5*w34 + C_1_7*w48 - C_2_0*w43 - C_2_3*w35 - C_2_4*w46 - C_2_7*w44 + tmp132 + tmp137 + tmp208 + tmp209 + tmp210 + tmp211;
//                                         EM_S[INDEX4(k,m,0,1,numEq,numComp,8)]+=C_0_0*w50 + C_0_1*w41 + C_0_6*w33 + C_0_7*w49 + tmp126 + tmp128 + tmp242 + tmp243 + tmp244 + tmp245 + tmp246 + tmp247 + tmp248 + tmp249;
//                                         EM_S[INDEX4(k,m,0,2,numEq,numComp,8)]+=-C_1_0*w47 + C_1_2*w42 + C_1_5*w34 - C_1_7*w48 + tmp138 + tmp139 + tmp140 + tmp142 + tmp143 + tmp144 + tmp146 + tmp147 + tmp173 + tmp177;
//                                         EM_S[INDEX4(k,m,0,3,numEq,numComp,8)]+=tmp100 + tmp101 + tmp102 + tmp103 + tmp40 + tmp45 + tmp96 + tmp97 + tmp98 + tmp99;
//                                         EM_S[INDEX4(k,m,0,4,numEq,numComp,8)]+=C_2_0*w43 + C_2_3*w35 + C_2_4*w46 + C_2_7*w44 + tmp114 + tmp115 + tmp116 + tmp117 + tmp118 + tmp119 + tmp120 + tmp121 + tmp92 + tmp93;
//                                         EM_S[INDEX4(k,m,0,5,numEq,numComp,8)]+=tmp192 + tmp193 + tmp194 + tmp196 + tmp197 + tmp198 + tmp224 + tmp225 + tmp226 + tmp227;
//                                         EM_S[INDEX4(k,m,0,6,numEq,numComp,8)]+=tmp232 + tmp234 + tmp238 + tmp239 + tmp240 + tmp241 + tmp260 + tmp261 + tmp262 + tmp263;
//                                         EM_S[INDEX4(k,m,0,7,numEq,numComp,8)]+=tmp60 + tmp61 + tmp62 + tmp63 + tmp64 + tmp65;
//                                         EM_S[INDEX4(k,m,1,0,numEq,numComp,8)]+=-C_0_0*w41 - C_0_1*w50 - C_0_6*w49 - C_0_7*w33 + tmp148 + tmp149 + tmp242 + tmp243 + tmp244 + tmp245 + tmp246 + tmp247 + tmp248 + tmp249;
//                                         EM_S[INDEX4(k,m,1,1,numEq,numComp,8)]+=C_0_0*w41 + C_0_1*w50 + C_0_6*w49 + C_0_7*w33 + C_1_1*w47 - C_1_3*w42 - C_1_4*w34 + C_1_6*w48 - C_2_1*w43 - C_2_2*w35 - C_2_5*w46 - C_2_6*w44 + tmp151 + tmp155 + tmp164 + tmp167 + tmp168 + tmp169;
//                                         EM_S[INDEX4(k,m,1,2,numEq,numComp,8)]+=tmp101 + tmp103 + tmp40 + tmp42 + tmp44 + tmp45 + tmp46 + tmp48 + tmp96 + tmp98;
//                                         EM_S[INDEX4(k,m,1,3,numEq,numComp,8)]+=-C_1_1*w47 + C_1_3*w42 + C_1_4*w34 - C_1_6*w48 + tmp20 + tmp21 + tmp22 + tmp23 + tmp24 + tmp25 + tmp26 + tmp27 + tmp28 + tmp29;
//                                         EM_S[INDEX4(k,m,1,4,numEq,numComp,8)]+=tmp193 + tmp197 + tmp214 + tmp215 + tmp216 + tmp217 + tmp224 + tmp225 + tmp226 + tmp227;
//                                         EM_S[INDEX4(k,m,1,5,numEq,numComp,8)]+=C_2_1*w43 + C_2_2*w35 + C_2_5*w46 + C_2_6*w44 + tmp70 + tmp75 + tmp83 + tmp84 + tmp85 + tmp86 + tmp88 + tmp89 + tmp90 + tmp91;
//                                         EM_S[INDEX4(k,m,1,6,numEq,numComp,8)]+=tmp60 + tmp61 + tmp63 + tmp65 + tmp80 + tmp81;
//                                         EM_S[INDEX4(k,m,1,7,numEq,numComp,8)]+=tmp181 + tmp183 + tmp184 + tmp186 + tmp187 + tmp189 + tmp254 + tmp255 + tmp256 + tmp257;
//                                         EM_S[INDEX4(k,m,2,0,numEq,numComp,8)]+=-C_1_0*w42 + C_1_2*w47 + C_1_5*w48 - C_1_7*w34 + tmp138 + tmp139 + tmp140 + tmp141 + tmp142 + tmp143 + tmp144 + tmp145 + tmp146 + tmp147;
//                                         EM_S[INDEX4(k,m,2,1,numEq,numComp,8)]+=tmp100 + tmp102 + tmp40 + tmp41 + tmp43 + tmp45 + tmp47 + tmp49 + tmp97 + tmp99;
//                                         EM_S[INDEX4(k,m,2,2,numEq,numComp,8)]+=-C_0_2*w50 - C_0_3*w41 - C_0_4*w33 - C_0_5*w49 + C_1_0*w42 - C_1_2*w47 - C_1_5*w48 + C_1_7*w34 - C_2_1*w35 - C_2_2*w43 - C_2_5*w44 - C_2_6*w46 + tmp152 + tmp154 + tmp164 + tmp165 + tmp166 + tmp167;
//                                         EM_S[INDEX4(k,m,2,3,numEq,numComp,8)]+=C_0_2*w50 + C_0_3*w41 + C_0_4*w33 + C_0_5*w49 + tmp200 + tmp201 + tmp202 + tmp203 + tmp204 + tmp205 + tmp206 + tmp207 + tmp222 + tmp223;
//                                         EM_S[INDEX4(k,m,2,4,numEq,numComp,8)]+=tmp229 + tmp231 + tmp232 + tmp234 + tmp235 + tmp237 + tmp260 + tmp261 + tmp262 + tmp263;
//                                         EM_S[INDEX4(k,m,2,5,numEq,numComp,8)]+=tmp60 + tmp62 + tmp63 + tmp64 + tmp94 + tmp95;
//                                         EM_S[INDEX4(k,m,2,6,numEq,numComp,8)]+=C_2_1*w35 + C_2_2*w43 + C_2_5*w44 + C_2_6*w46 + tmp70 + tmp71 + tmp72 + tmp73 + tmp74 + tmp75 + tmp76 + tmp77 + tmp78 + tmp79;
//                                         EM_S[INDEX4(k,m,2,7,numEq,numComp,8)]+=tmp107 + tmp111 + tmp156 + tmp157 + tmp158 + tmp159 + tmp160 + tmp161 + tmp162 + tmp163;
//                                         EM_S[INDEX4(k,m,3,0,numEq,numComp,8)]+=tmp40 + tmp41 + tmp42 + tmp43 + tmp44 + tmp45 + tmp46 + tmp47 + tmp48 + tmp49;
//                                         EM_S[INDEX4(k,m,3,1,numEq,numComp,8)]+=-C_1_1*w42 + C_1_3*w47 + C_1_4*w48 - C_1_6*w34 + tmp20 + tmp21 + tmp22 + tmp24 + tmp25 + tmp26 + tmp28 + tmp29 + tmp33 + tmp37;
//                                         EM_S[INDEX4(k,m,3,2,numEq,numComp,8)]+=-C_0_2*w41 - C_0_3*w50 - C_0_4*w49 - C_0_5*w33 + tmp200 + tmp201 + tmp202 + tmp203 + tmp204 + tmp205 + tmp206 + tmp207 + tmp54 + tmp56;
//                                         EM_S[INDEX4(k,m,3,3,numEq,numComp,8)]+=C_0_2*w41 + C_0_3*w50 + C_0_4*w49 + C_0_5*w33 + C_1_1*w42 - C_1_3*w47 - C_1_4*w48 + C_1_6*w34 - C_2_0*w35 - C_2_3*w43 - C_2_4*w44 - C_2_7*w46 + tmp132 + tmp133 + tmp134 + tmp135 + tmp136 + tmp137;
//                                         EM_S[INDEX4(k,m,3,4,numEq,numComp,8)]+=tmp60 + tmp63 + tmp80 + tmp81 + tmp94 + tmp95;
//                                         EM_S[INDEX4(k,m,3,5,numEq,numComp,8)]+=tmp184 + tmp186 + tmp218 + tmp219 + tmp220 + tmp221 + tmp254 + tmp255 + tmp256 + tmp257;
//                                         EM_S[INDEX4(k,m,3,6,numEq,numComp,8)]+=tmp106 + tmp107 + tmp108 + tmp110 + tmp111 + tmp112 + tmp160 + tmp161 + tmp162 + tmp163;
//                                         EM_S[INDEX4(k,m,3,7,numEq,numComp,8)]+=C_2_0*w35 + C_2_3*w43 + C_2_4*w44 + C_2_7*w46 + tmp1 + tmp2 + tmp3 + tmp4 + tmp6 + tmp7 + tmp8 + tmp9 + tmp92 + tmp93;
//                                         EM_S[INDEX4(k,m,4,0,numEq,numComp,8)]+=-C_2_0*w46 - C_2_3*w44 - C_2_4*w43 - C_2_7*w35 + tmp0 + tmp114 + tmp115 + tmp116 + tmp117 + tmp118 + tmp119 + tmp120 + tmp121 + tmp5;
//                                         EM_S[INDEX4(k,m,4,1,numEq,numComp,8)]+=tmp190 + tmp191 + tmp192 + tmp193 + tmp194 + tmp195 + tmp196 + tmp197 + tmp198 + tmp199;
//                                         EM_S[INDEX4(k,m,4,2,numEq,numComp,8)]+=tmp228 + tmp230 + tmp232 + tmp233 + tmp234 + tmp236 + tmp238 + tmp239 + tmp240 + tmp241;
//                                         EM_S[INDEX4(k,m,4,3,numEq,numComp,8)]+=tmp258 + tmp259 + tmp61 + tmp62 + tmp64 + tmp65;
//                                         EM_S[INDEX4(k,m,4,4,numEq,numComp,8)]+=-C_0_2*w33 - C_0_3*w49 - C_0_4*w50 - C_0_5*w41 - C_1_1*w34 + C_1_3*w48 + C_1_4*w47 - C_1_6*w42 + C_2_0*w46 + C_2_3*w44 + C_2_4*w43 + C_2_7*w35 + tmp150 + tmp151 + tmp152 + tmp153 + tmp154 + tmp155;
//                                         EM_S[INDEX4(k,m,4,5,numEq,numComp,8)]+=C_0_2*w33 + C_0_3*w49 + C_0_4*w50 + C_0_5*w41 + tmp222 + tmp223 + tmp50 + tmp51 + tmp52 + tmp53 + tmp55 + tmp57 + tmp58 + tmp59;
//                                         EM_S[INDEX4(k,m,4,6,numEq,numComp,8)]+=C_1_1*w34 - C_1_3*w48 - C_1_4*w47 + C_1_6*w42 + tmp23 + tmp27 + tmp30 + tmp31 + tmp32 + tmp34 + tmp35 + tmp36 + tmp38 + tmp39;
//                                         EM_S[INDEX4(k,m,4,7,numEq,numComp,8)]+=tmp10 + tmp11 + tmp12 + tmp13 + tmp14 + tmp15 + tmp16 + tmp17 + tmp18 + tmp19;
//                                         EM_S[INDEX4(k,m,5,0,numEq,numComp,8)]+=tmp190 + tmp191 + tmp193 + tmp195 + tmp197 + tmp199 + tmp214 + tmp215 + tmp216 + tmp217;
//                                         EM_S[INDEX4(k,m,5,1,numEq,numComp,8)]+=-C_2_1*w46 - C_2_2*w44 - C_2_5*w43 - C_2_6*w35 + tmp82 + tmp83 + tmp84 + tmp85 + tmp86 + tmp87 + tmp88 + tmp89 + tmp90 + tmp91;
//                                         EM_S[INDEX4(k,m,5,2,numEq,numComp,8)]+=tmp258 + tmp259 + tmp61 + tmp65 + tmp80 + tmp81;
//                                         EM_S[INDEX4(k,m,5,3,numEq,numComp,8)]+=tmp180 + tmp181 + tmp182 + tmp183 + tmp184 + tmp185 + tmp186 + tmp187 + tmp188 + tmp189;
//                                         EM_S[INDEX4(k,m,5,4,numEq,numComp,8)]+=-C_0_2*w49 - C_0_3*w33 - C_0_4*w41 - C_0_5*w50 + tmp50 + tmp51 + tmp52 + tmp53 + tmp54 + tmp55 + tmp56 + tmp57 + tmp58 + tmp59;
//                                         EM_S[INDEX4(k,m,5,5,numEq,numComp,8)]+=C_0_2*w49 + C_0_3*w33 + C_0_4*w41 + C_0_5*w50 - C_1_0*w34 + C_1_2*w48 + C_1_5*w47 - C_1_7*w42 + C_2_1*w46 + C_2_2*w44 + C_2_5*w43 + C_2_6*w35 + tmp134 + tmp135 + tmp208 + tmp211 + tmp212 + tmp213;
//                                         EM_S[INDEX4(k,m,5,6,numEq,numComp,8)]+=tmp10 + tmp11 + tmp13 + tmp15 + tmp17 + tmp19 + tmp66 + tmp67 + tmp68 + tmp69;
//                                         EM_S[INDEX4(k,m,5,7,numEq,numComp,8)]+=C_1_0*w34 - C_1_2*w48 - C_1_5*w47 + C_1_7*w42 + tmp170 + tmp171 + tmp172 + tmp173 + tmp174 + tmp175 + tmp176 + tmp177 + tmp178 + tmp179;
//                                         EM_S[INDEX4(k,m,6,0,numEq,numComp,8)]+=tmp228 + tmp229 + tmp230 + tmp231 + tmp232 + tmp233 + tmp234 + tmp235 + tmp236 + tmp237;
//                                         EM_S[INDEX4(k,m,6,1,numEq,numComp,8)]+=tmp258 + tmp259 + tmp62 + tmp64 + tmp94 + tmp95;
//                                         EM_S[INDEX4(k,m,6,2,numEq,numComp,8)]+=-C_2_1*w44 - C_2_2*w46 - C_2_5*w35 - C_2_6*w43 + tmp71 + tmp72 + tmp73 + tmp74 + tmp76 + tmp77 + tmp78 + tmp79 + tmp82 + tmp87;
//                                         EM_S[INDEX4(k,m,6,3,numEq,numComp,8)]+=tmp104 + tmp105 + tmp107 + tmp109 + tmp111 + tmp113 + tmp156 + tmp157 + tmp158 + tmp159;
//                                         EM_S[INDEX4(k,m,6,4,numEq,numComp,8)]+=C_1_1*w48 - C_1_3*w34 - C_1_4*w42 + C_1_6*w47 + tmp30 + tmp31 + tmp32 + tmp33 + tmp34 + tmp35 + tmp36 + tmp37 + tmp38 + tmp39;
//                                         EM_S[INDEX4(k,m,6,5,numEq,numComp,8)]+=tmp10 + tmp12 + tmp14 + tmp15 + tmp16 + tmp18 + tmp250 + tmp251 + tmp252 + tmp253;
//                                         EM_S[INDEX4(k,m,6,6,numEq,numComp,8)]+=-C_0_0*w33 - C_0_1*w49 - C_0_6*w50 - C_0_7*w41 - C_1_1*w48 + C_1_3*w34 + C_1_4*w42 - C_1_6*w47 + C_2_1*w44 + C_2_2*w46 + C_2_5*w35 + C_2_6*w43 + tmp133 + tmp136 + tmp209 + tmp210 + tmp212 + tmp213;
//                                         EM_S[INDEX4(k,m,6,7,numEq,numComp,8)]+=C_0_0*w33 + C_0_1*w49 + C_0_6*w50 + C_0_7*w41 + tmp122 + tmp123 + tmp124 + tmp125 + tmp126 + tmp127 + tmp128 + tmp129 + tmp130 + tmp131;
//                                         EM_S[INDEX4(k,m,7,0,numEq,numComp,8)]+=tmp258 + tmp259 + tmp80 + tmp81 + tmp94 + tmp95;
//                                         EM_S[INDEX4(k,m,7,1,numEq,numComp,8)]+=tmp180 + tmp182 + tmp184 + tmp185 + tmp186 + tmp188 + tmp218 + tmp219 + tmp220 + tmp221;
//                                         EM_S[INDEX4(k,m,7,2,numEq,numComp,8)]+=tmp104 + tmp105 + tmp106 + tmp107 + tmp108 + tmp109 + tmp110 + tmp111 + tmp112 + tmp113;
//                                         EM_S[INDEX4(k,m,7,3,numEq,numComp,8)]+=-C_2_0*w44 - C_2_3*w46 - C_2_4*w35 - C_2_7*w43 + tmp0 + tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6 + tmp7 + tmp8 + tmp9;
//                                         EM_S[INDEX4(k,m,7,4,numEq,numComp,8)]+=tmp10 + tmp15 + tmp250 + tmp251 + tmp252 + tmp253 + tmp66 + tmp67 + tmp68 + tmp69;
//                                         EM_S[INDEX4(k,m,7,5,numEq,numComp,8)]+=C_1_0*w48 - C_1_2*w34 - C_1_5*w42 + C_1_7*w47 + tmp141 + tmp145 + tmp170 + tmp171 + tmp172 + tmp174 + tmp175 + tmp176 + tmp178 + tmp179;
//                                         EM_S[INDEX4(k,m,7,6,numEq,numComp,8)]+=-C_0_0*w49 - C_0_1*w33 - C_0_6*w41 - C_0_7*w50 + tmp122 + tmp123 + tmp124 + tmp125 + tmp127 + tmp129 + tmp130 + tmp131 + tmp148 + tmp149;
//                                         EM_S[INDEX4(k,m,7,7,numEq,numComp,8)]+=C_0_0*w49 + C_0_1*w33 + C_0_6*w41 + C_0_7*w50 - C_1_0*w48 + C_1_2*w34 + C_1_5*w42 - C_1_7*w47 + C_2_0*w44 + C_2_3*w46 + C_2_4*w35 + C_2_7*w43 + tmp150 + tmp153 + tmp165 + tmp166 + tmp168 + tmp169;
//                                     }
//                                 }
//                             } else { // constant data
//                                 for (index_t k=0; k<numEq; k++) {
//                                     for (index_t m=0; m<numComp; m++) {
//                                         const Scalar wC0 = C_p[INDEX3(k,m,0,numEq,numComp)]*w55;
//                                         const Scalar wC1 = C_p[INDEX3(k,m,1,numEq,numComp)]*w56;
//                                         const Scalar wC2 = C_p[INDEX3(k,m,2,numEq,numComp)]*w54;
//                                         EM_S[INDEX4(k,m,0,0,numEq,numComp,8)]+= 4.*wC0 + 4.*wC1 + 4.*wC2;
//                                         EM_S[INDEX4(k,m,0,1,numEq,numComp,8)]+=-4.*wC0 + 2.*wC1 + 2.*wC2;
//                                         EM_S[INDEX4(k,m,0,2,numEq,numComp,8)]+= 2.*wC0 - 4.*wC1 + 2.*wC2;
//                                         EM_S[INDEX4(k,m,0,3,numEq,numComp,8)]+=-2.*wC0 - 2.*wC1 +    wC2;
//                                         EM_S[INDEX4(k,m,0,4,numEq,numComp,8)]+= 2.*wC0 + 2.*wC1 - 4.*wC2;
//                                         EM_S[INDEX4(k,m,0,5,numEq,numComp,8)]+=-2.*wC0 +    wC1 - 2.*wC2;
//                                         EM_S[INDEX4(k,m,0,6,numEq,numComp,8)]+=    wC0 - 2.*wC1 - 2.*wC2;
//                                         EM_S[INDEX4(k,m,0,7,numEq,numComp,8)]+=   -wC0 -    wC1 -    wC2;
//                                         EM_S[INDEX4(k,m,1,0,numEq,numComp,8)]+= 4.*wC0 + 2.*wC1 + 2.*wC2;
//                                         EM_S[INDEX4(k,m,1,1,numEq,numComp,8)]+=-4.*wC0 + 4.*wC1 + 4.*wC2;
//                                         EM_S[INDEX4(k,m,1,2,numEq,numComp,8)]+= 2.*wC0 - 2.*wC1 +    wC2;
//                                         EM_S[INDEX4(k,m,1,3,numEq,numComp,8)]+=-2.*wC0 - 4.*wC1 + 2.*wC2;
//                                         EM_S[INDEX4(k,m,1,4,numEq,numComp,8)]+= 2.*wC0 +    wC1 - 2.*wC2;
//                                         EM_S[INDEX4(k,m,1,5,numEq,numComp,8)]+=-2.*wC0 + 2.*wC1 - 4.*wC2;
//                                         EM_S[INDEX4(k,m,1,6,numEq,numComp,8)]+=    wC0 -    wC1 -    wC2;
//                                         EM_S[INDEX4(k,m,1,7,numEq,numComp,8)]+=   -wC0 - 2.*wC1 - 2.*wC2;
//                                         EM_S[INDEX4(k,m,2,0,numEq,numComp,8)]+= 2.*wC0 + 4.*wC1 + 2.*wC2;
//                                         EM_S[INDEX4(k,m,2,1,numEq,numComp,8)]+=-2.*wC0 + 2.*wC1 +    wC2;
//                                         EM_S[INDEX4(k,m,2,2,numEq,numComp,8)]+= 4.*wC0 - 4.*wC1 + 4.*wC2;
//                                         EM_S[INDEX4(k,m,2,3,numEq,numComp,8)]+=-4.*wC0 - 2.*wC1 + 2.*wC2;
//                                         EM_S[INDEX4(k,m,2,4,numEq,numComp,8)]+=    wC0 + 2.*wC1 - 2.*wC2;
//                                         EM_S[INDEX4(k,m,2,5,numEq,numComp,8)]+=   -wC0 +    wC1 -    wC2;
//                                         EM_S[INDEX4(k,m,2,6,numEq,numComp,8)]+= 2.*wC0 - 2.*wC1 - 4.*wC2;
//                                         EM_S[INDEX4(k,m,2,7,numEq,numComp,8)]+=-2.*wC0 -    wC1 - 2.*wC2;
//                                         EM_S[INDEX4(k,m,3,0,numEq,numComp,8)]+= 2.*wC0 + 2.*wC1 +    wC2;
//                                         EM_S[INDEX4(k,m,3,1,numEq,numComp,8)]+=-2.*wC0 + 4.*wC1 + 2.*wC2;
//                                         EM_S[INDEX4(k,m,3,2,numEq,numComp,8)]+= 4.*wC0 - 2.*wC1 + 2.*wC2;
//                                         EM_S[INDEX4(k,m,3,3,numEq,numComp,8)]+=-4.*wC0 - 4.*wC1 + 4.*wC2;
//                                         EM_S[INDEX4(k,m,3,4,numEq,numComp,8)]+=    wC0 +    wC1 -    wC2;
//                                         EM_S[INDEX4(k,m,3,5,numEq,numComp,8)]+=   -wC0 + 2.*wC1 - 2.*wC2;
//                                         EM_S[INDEX4(k,m,3,6,numEq,numComp,8)]+= 2.*wC0 -    wC1 - 2.*wC2;
//                                         EM_S[INDEX4(k,m,3,7,numEq,numComp,8)]+=-2.*wC0 - 2.*wC1 - 4.*wC2;
//                                         EM_S[INDEX4(k,m,4,0,numEq,numComp,8)]+= 2.*wC0 + 2.*wC1 + 4.*wC2;
//                                         EM_S[INDEX4(k,m,4,1,numEq,numComp,8)]+=-2.*wC0 +    wC1 + 2.*wC2;
//                                         EM_S[INDEX4(k,m,4,2,numEq,numComp,8)]+=    wC0 - 2.*wC1 + 2.*wC2;
//                                         EM_S[INDEX4(k,m,4,3,numEq,numComp,8)]+=   -wC0 -    wC1 +    wC2;
//                                         EM_S[INDEX4(k,m,4,4,numEq,numComp,8)]+= 4.*wC0 + 4.*wC1 - 4.*wC2;
//                                         EM_S[INDEX4(k,m,4,5,numEq,numComp,8)]+=-4.*wC0 + 2.*wC1 - 2.*wC2;
//                                         EM_S[INDEX4(k,m,4,6,numEq,numComp,8)]+= 2.*wC0 - 4.*wC1 - 2.*wC2;
//                                         EM_S[INDEX4(k,m,4,7,numEq,numComp,8)]+=-2.*wC0 - 2.*wC1 -    wC2;
//                                         EM_S[INDEX4(k,m,5,0,numEq,numComp,8)]+= 2.*wC0 +    wC1 + 2.*wC2;
//                                         EM_S[INDEX4(k,m,5,1,numEq,numComp,8)]+=-2.*wC0 + 2.*wC1 + 4.*wC2;
//                                         EM_S[INDEX4(k,m,5,2,numEq,numComp,8)]+=    wC0 -    wC1 +    wC2;
//                                         EM_S[INDEX4(k,m,5,3,numEq,numComp,8)]+=   -wC0 - 2.*wC1 + 2.*wC2;
//                                         EM_S[INDEX4(k,m,5,4,numEq,numComp,8)]+= 4.*wC0 + 2.*wC1 - 2.*wC2;
//                                         EM_S[INDEX4(k,m,5,5,numEq,numComp,8)]+=-4.*wC0 + 4.*wC1 - 4.*wC2;
//                                         EM_S[INDEX4(k,m,5,6,numEq,numComp,8)]+= 2.*wC0 - 2.*wC1 -    wC2;
//                                         EM_S[INDEX4(k,m,5,7,numEq,numComp,8)]+=-2.*wC0 - 4.*wC1 - 2.*wC2;
//                                         EM_S[INDEX4(k,m,6,0,numEq,numComp,8)]+=    wC0 + 2.*wC1 + 2.*wC2;
//                                         EM_S[INDEX4(k,m,6,1,numEq,numComp,8)]+=   -wC0 +    wC1 +    wC2;
//                                         EM_S[INDEX4(k,m,6,2,numEq,numComp,8)]+= 2.*wC0 - 2.*wC1 + 4.*wC2;
//                                         EM_S[INDEX4(k,m,6,3,numEq,numComp,8)]+=-2.*wC0 -    wC1 + 2.*wC2;
//                                         EM_S[INDEX4(k,m,6,4,numEq,numComp,8)]+= 2.*wC0 + 4.*wC1 - 2.*wC2;
//                                         EM_S[INDEX4(k,m,6,5,numEq,numComp,8)]+=-2.*wC0 + 2.*wC1 -    wC2;
//                                         EM_S[INDEX4(k,m,6,6,numEq,numComp,8)]+= 4.*wC0 - 4.*wC1 - 4.*wC2;
//                                         EM_S[INDEX4(k,m,6,7,numEq,numComp,8)]+=-4.*wC0 - 2.*wC1 - 2.*wC2;
//                                         EM_S[INDEX4(k,m,7,0,numEq,numComp,8)]+=    wC0 +    wC1 +    wC2;
//                                         EM_S[INDEX4(k,m,7,1,numEq,numComp,8)]+=   -wC0 + 2.*wC1 + 2.*wC2;
//                                         EM_S[INDEX4(k,m,7,2,numEq,numComp,8)]+= 2.*wC0 -    wC1 + 2.*wC2;
//                                         EM_S[INDEX4(k,m,7,3,numEq,numComp,8)]+=-2.*wC0 - 2.*wC1 + 4.*wC2;
//                                         EM_S[INDEX4(k,m,7,4,numEq,numComp,8)]+= 2.*wC0 + 2.*wC1 -    wC2;
//                                         EM_S[INDEX4(k,m,7,5,numEq,numComp,8)]+=-2.*wC0 + 4.*wC1 - 2.*wC2;
//                                         EM_S[INDEX4(k,m,7,6,numEq,numComp,8)]+= 4.*wC0 - 2.*wC1 - 2.*wC2;
//                                         EM_S[INDEX4(k,m,7,7,numEq,numComp,8)]+=-4.*wC0 - 4.*wC1 - 4.*wC2;
//                                     }
//                                 }
//                             }
//                         }
//                         ///////////////
//                         // process D //
//                         ///////////////
//                         if (!D.isEmpty()) {
//                             const Scalar* D_p = D.getSampleDataRO(id, zero);
//                             if (D.actsExpanded()) {
//                                 for (index_t k=0; k<numEq; k++) {
//                                     for (index_t m=0; m<numComp; m++) {
//                                         const Scalar D_0 = D_p[INDEX3(k,m,0,numEq,numComp)];
//                                         const Scalar D_1 = D_p[INDEX3(k,m,1,numEq,numComp)];
//                                         const Scalar D_2 = D_p[INDEX3(k,m,2,numEq,numComp)];
//                                         const Scalar D_3 = D_p[INDEX3(k,m,3,numEq,numComp)];
//                                         const Scalar D_4 = D_p[INDEX3(k,m,4,numEq,numComp)];
//                                         const Scalar D_5 = D_p[INDEX3(k,m,5,numEq,numComp)];
//                                         const Scalar D_6 = D_p[INDEX3(k,m,6,numEq,numComp)];
//                                         const Scalar D_7 = D_p[INDEX3(k,m,7,numEq,numComp)];
//                                         const Scalar tmp0 = w59*(D_3 + D_7);
//                                         const Scalar tmp1 = w57*(D_0 + D_4);
//                                         const Scalar tmp2 = w58*(D_1 + D_2 + D_5 + D_6);
//                                         const Scalar tmp3 = w60*(D_0 + D_1 + D_2 + D_3);
//                                         const Scalar tmp4 = w61*(D_4 + D_5 + D_6 + D_7);
//                                         const Scalar tmp5 = w59*(D_1 + D_3);
//                                         const Scalar tmp6 = w57*(D_4 + D_6);
//                                         const Scalar tmp7 = w58*(D_0 + D_2 + D_5 + D_7);
//                                         const Scalar tmp8 = w59*(D_4 + D_6);
//                                         const Scalar tmp9 = w57*(D_1 + D_3);
//                                         const Scalar tmp10 = w60*(D_4 + D_5 + D_6 + D_7);
//                                         const Scalar tmp11 = w61*(D_0 + D_1 + D_2 + D_3);
//                                         const Scalar tmp12 = w59*(D_4 + D_5);
//                                         const Scalar tmp13 = w57*(D_2 + D_3);
//                                         const Scalar tmp14 = w58*(D_0 + D_1 + D_6 + D_7);
//                                         const Scalar tmp15 = w58*(D_0 + D_1 + D_2 + D_3 + D_4 + D_5 + D_6 + D_7);
//                                         const Scalar tmp16 = w59*(D_2 + D_6);
//                                         const Scalar tmp17 = w57*(D_1 + D_5);
//                                         const Scalar tmp18 = w58*(D_0 + D_3 + D_4 + D_7);
//                                         const Scalar tmp19 = w59*(D_1 + D_5);
//                                         const Scalar tmp20 = w57*(D_2 + D_6);
//                                         const Scalar tmp21 = w60*(D_0 + D_1 + D_4 + D_5);
//                                         const Scalar tmp22 = w61*(D_2 + D_3 + D_6 + D_7);
//                                         const Scalar tmp23 = w59*(D_0 + D_4);
//                                         const Scalar tmp24 = w57*(D_3 + D_7);
//                                         const Scalar tmp25 = w59*(D_6 + D_7);
//                                         const Scalar tmp26 = w57*(D_0 + D_1);
//                                         const Scalar tmp27 = w58*(D_2 + D_3 + D_4 + D_5);
//                                         const Scalar tmp28 = w60*(D_0 + D_5 + D_6);
//                                         const Scalar tmp29 = w61*(D_1 + D_2 + D_7);
//                                         const Scalar tmp30 = w59*(D_0 + D_2);
//                                         const Scalar tmp31 = w57*(D_5 + D_7);
//                                         const Scalar tmp32 = w58*(D_1 + D_3 + D_4 + D_6);
//                                         const Scalar tmp33 = w60*(D_1 + D_2 + D_7);
//                                         const Scalar tmp34 = w61*(D_0 + D_5 + D_6);
//                                         const Scalar tmp35 = w60*(D_1 + D_4 + D_7);
//                                         const Scalar tmp36 = w61*(D_0 + D_3 + D_6);
//                                         const Scalar tmp37 = w60*(D_1 + D_2 + D_4);
//                                         const Scalar tmp38 = w61*(D_3 + D_5 + D_6);
//                                         const Scalar tmp39 = w59*(D_5 + D_7);
//                                         const Scalar tmp40 = w57*(D_0 + D_2);
//                                         const Scalar tmp41 = w60*(D_0 + D_2 + D_4 + D_6);
//                                         const Scalar tmp42 = w61*(D_1 + D_3 + D_5 + D_7);
//                                         const Scalar tmp43 = w60*(D_2 + D_3 + D_6 + D_7);
//                                         const Scalar tmp44 = w61*(D_0 + D_1 + D_4 + D_5);
//                                         const Scalar tmp45 = w60*(D_2 + D_4 + D_7);
//                                         const Scalar tmp46 = w61*(D_0 + D_3 + D_5);
//                                         const Scalar tmp47 = w59*(D_2 + D_3);
//                                         const Scalar tmp48 = w57*(D_4 + D_5);
//                                         const Scalar tmp49 = w60*(D_3 + D_5 + D_6);
//                                         const Scalar tmp50 = w61*(D_1 + D_2 + D_4);
//                                         const Scalar tmp51 = w60*(D_0 + D_3 + D_5);
//                                         const Scalar tmp52 = w61*(D_2 + D_4 + D_7);
//                                         const Scalar tmp53 = w60*(D_0 + D_3 + D_6);
//                                         const Scalar tmp54 = w61*(D_1 + D_4 + D_7);
//                                         const Scalar tmp55 = w60*(D_1 + D_3 + D_5 + D_7);
//                                         const Scalar tmp56 = w61*(D_0 + D_2 + D_4 + D_6);
//                                         const Scalar tmp57 = w59*(D_0 + D_1);
//                                         const Scalar tmp58 = w57*(D_6 + D_7);
//                                         EM_S[INDEX4(k,m,0,0,numEq,numComp,8)]+=D_0*w62 + D_7*w63 + tmp49 + tmp50;
//                                         EM_S[INDEX4(k,m,0,1,numEq,numComp,8)]+=tmp27 + tmp57 + tmp58;
//                                         EM_S[INDEX4(k,m,0,2,numEq,numComp,8)]+=tmp30 + tmp31 + tmp32;
//                                         EM_S[INDEX4(k,m,0,3,numEq,numComp,8)]+=tmp10 + tmp11;
//                                         EM_S[INDEX4(k,m,0,4,numEq,numComp,8)]+=tmp2 + tmp23 + tmp24;
//                                         EM_S[INDEX4(k,m,0,5,numEq,numComp,8)]+=tmp43 + tmp44;
//                                         EM_S[INDEX4(k,m,0,6,numEq,numComp,8)]+=tmp55 + tmp56;
//                                         EM_S[INDEX4(k,m,0,7,numEq,numComp,8)]+=tmp15;
//                                         EM_S[INDEX4(k,m,1,0,numEq,numComp,8)]+=tmp27 + tmp57 + tmp58;
//                                         EM_S[INDEX4(k,m,1,1,numEq,numComp,8)]+=D_1*w62 + D_6*w63 + tmp45 + tmp46;
//                                         EM_S[INDEX4(k,m,1,2,numEq,numComp,8)]+=tmp10 + tmp11;
//                                         EM_S[INDEX4(k,m,1,3,numEq,numComp,8)]+=tmp5 + tmp6 + tmp7;
//                                         EM_S[INDEX4(k,m,1,4,numEq,numComp,8)]+=tmp43 + tmp44;
//                                         EM_S[INDEX4(k,m,1,5,numEq,numComp,8)]+=tmp18 + tmp19 + tmp20;
//                                         EM_S[INDEX4(k,m,1,6,numEq,numComp,8)]+=tmp15;
//                                         EM_S[INDEX4(k,m,1,7,numEq,numComp,8)]+=tmp41 + tmp42;
//                                         EM_S[INDEX4(k,m,2,0,numEq,numComp,8)]+=tmp30 + tmp31 + tmp32;
//                                         EM_S[INDEX4(k,m,2,1,numEq,numComp,8)]+=tmp10 + tmp11;
//                                         EM_S[INDEX4(k,m,2,2,numEq,numComp,8)]+=D_2*w62 + D_5*w63 + tmp35 + tmp36;
//                                         EM_S[INDEX4(k,m,2,3,numEq,numComp,8)]+=tmp14 + tmp47 + tmp48;
//                                         EM_S[INDEX4(k,m,2,4,numEq,numComp,8)]+=tmp55 + tmp56;
//                                         EM_S[INDEX4(k,m,2,5,numEq,numComp,8)]+=tmp15;
//                                         EM_S[INDEX4(k,m,2,6,numEq,numComp,8)]+=tmp16 + tmp17 + tmp18;
//                                         EM_S[INDEX4(k,m,2,7,numEq,numComp,8)]+=tmp21 + tmp22;
//                                         EM_S[INDEX4(k,m,3,0,numEq,numComp,8)]+=tmp10 + tmp11;
//                                         EM_S[INDEX4(k,m,3,1,numEq,numComp,8)]+=tmp5 + tmp6 + tmp7;
//                                         EM_S[INDEX4(k,m,3,2,numEq,numComp,8)]+=tmp14 + tmp47 + tmp48;
//                                         EM_S[INDEX4(k,m,3,3,numEq,numComp,8)]+=D_3*w62 + D_4*w63 + tmp28 + tmp29;
//                                         EM_S[INDEX4(k,m,3,4,numEq,numComp,8)]+=tmp15;
//                                         EM_S[INDEX4(k,m,3,5,numEq,numComp,8)]+=tmp41 + tmp42;
//                                         EM_S[INDEX4(k,m,3,6,numEq,numComp,8)]+=tmp21 + tmp22;
//                                         EM_S[INDEX4(k,m,3,7,numEq,numComp,8)]+=tmp0 + tmp1 + tmp2;
//                                         EM_S[INDEX4(k,m,4,0,numEq,numComp,8)]+=tmp2 + tmp23 + tmp24;
//                                         EM_S[INDEX4(k,m,4,1,numEq,numComp,8)]+=tmp43 + tmp44;
//                                         EM_S[INDEX4(k,m,4,2,numEq,numComp,8)]+=tmp55 + tmp56;
//                                         EM_S[INDEX4(k,m,4,3,numEq,numComp,8)]+=tmp15;
//                                         EM_S[INDEX4(k,m,4,4,numEq,numComp,8)]+=D_3*w63 + D_4*w62 + tmp33 + tmp34;
//                                         EM_S[INDEX4(k,m,4,5,numEq,numComp,8)]+=tmp12 + tmp13 + tmp14;
//                                         EM_S[INDEX4(k,m,4,6,numEq,numComp,8)]+=tmp7 + tmp8 + tmp9;
//                                         EM_S[INDEX4(k,m,4,7,numEq,numComp,8)]+=tmp3 + tmp4;
//                                         EM_S[INDEX4(k,m,5,0,numEq,numComp,8)]+=tmp43 + tmp44;
//                                         EM_S[INDEX4(k,m,5,1,numEq,numComp,8)]+=tmp18 + tmp19 + tmp20;
//                                         EM_S[INDEX4(k,m,5,2,numEq,numComp,8)]+=tmp15;
//                                         EM_S[INDEX4(k,m,5,3,numEq,numComp,8)]+=tmp41 + tmp42;
//                                         EM_S[INDEX4(k,m,5,4,numEq,numComp,8)]+=tmp12 + tmp13 + tmp14;
//                                         EM_S[INDEX4(k,m,5,5,numEq,numComp,8)]+=D_2*w63 + D_5*w62 + tmp53 + tmp54;
//                                         EM_S[INDEX4(k,m,5,6,numEq,numComp,8)]+=tmp3 + tmp4;
//                                         EM_S[INDEX4(k,m,5,7,numEq,numComp,8)]+=tmp32 + tmp39 + tmp40;
//                                         EM_S[INDEX4(k,m,6,0,numEq,numComp,8)]+=tmp55 + tmp56;
//                                         EM_S[INDEX4(k,m,6,1,numEq,numComp,8)]+=tmp15;
//                                         EM_S[INDEX4(k,m,6,2,numEq,numComp,8)]+=tmp16 + tmp17 + tmp18;
//                                         EM_S[INDEX4(k,m,6,3,numEq,numComp,8)]+=tmp21 + tmp22;
//                                         EM_S[INDEX4(k,m,6,4,numEq,numComp,8)]+=tmp7 + tmp8 + tmp9;
//                                         EM_S[INDEX4(k,m,6,5,numEq,numComp,8)]+=tmp3 + tmp4;
//                                         EM_S[INDEX4(k,m,6,6,numEq,numComp,8)]+=D_1*w63 + D_6*w62 + tmp51 + tmp52;
//                                         EM_S[INDEX4(k,m,6,7,numEq,numComp,8)]+=tmp25 + tmp26 + tmp27;
//                                         EM_S[INDEX4(k,m,7,0,numEq,numComp,8)]+=tmp15;
//                                         EM_S[INDEX4(k,m,7,1,numEq,numComp,8)]+=tmp41 + tmp42;
//                                         EM_S[INDEX4(k,m,7,2,numEq,numComp,8)]+=tmp21 + tmp22;
//                                         EM_S[INDEX4(k,m,7,3,numEq,numComp,8)]+=tmp0 + tmp1 + tmp2;
//                                         EM_S[INDEX4(k,m,7,4,numEq,numComp,8)]+=tmp3 + tmp4;
//                                         EM_S[INDEX4(k,m,7,5,numEq,numComp,8)]+=tmp32 + tmp39 + tmp40;
//                                         EM_S[INDEX4(k,m,7,6,numEq,numComp,8)]+=tmp25 + tmp26 + tmp27;
//                                         EM_S[INDEX4(k,m,7,7,numEq,numComp,8)]+=D_0*w63 + D_7*w62 + tmp37 + tmp38;
//                                     }
//                                  }
//                             } else { // constant data
//                                 for (index_t k=0; k<numEq; k++) {
//                                     for (index_t m=0; m<numComp; m++) {
//                                         const Scalar wD0 = 8.*D_p[INDEX2(k, m, numEq)]*w58;
//                                         EM_S[INDEX4(k,m,0,0,numEq,numComp,8)]+=8.*wD0;
//                                         EM_S[INDEX4(k,m,0,1,numEq,numComp,8)]+=4.*wD0;
//                                         EM_S[INDEX4(k,m,0,2,numEq,numComp,8)]+=4.*wD0;
//                                         EM_S[INDEX4(k,m,0,3,numEq,numComp,8)]+=2.*wD0;
//                                         EM_S[INDEX4(k,m,0,4,numEq,numComp,8)]+=4.*wD0;
//                                         EM_S[INDEX4(k,m,0,5,numEq,numComp,8)]+=2.*wD0;
//                                         EM_S[INDEX4(k,m,0,6,numEq,numComp,8)]+=2.*wD0;
//                                         EM_S[INDEX4(k,m,0,7,numEq,numComp,8)]+=   wD0;
//                                         EM_S[INDEX4(k,m,1,0,numEq,numComp,8)]+=4.*wD0;
//                                         EM_S[INDEX4(k,m,1,1,numEq,numComp,8)]+=8.*wD0;
//                                         EM_S[INDEX4(k,m,1,2,numEq,numComp,8)]+=2.*wD0;
//                                         EM_S[INDEX4(k,m,1,3,numEq,numComp,8)]+=4.*wD0;
//                                         EM_S[INDEX4(k,m,1,4,numEq,numComp,8)]+=2.*wD0;
//                                         EM_S[INDEX4(k,m,1,5,numEq,numComp,8)]+=4.*wD0;
//                                         EM_S[INDEX4(k,m,1,6,numEq,numComp,8)]+=   wD0;
//                                         EM_S[INDEX4(k,m,1,7,numEq,numComp,8)]+=2.*wD0;
//                                         EM_S[INDEX4(k,m,2,0,numEq,numComp,8)]+=4.*wD0;
//                                         EM_S[INDEX4(k,m,2,1,numEq,numComp,8)]+=2.*wD0;
//                                         EM_S[INDEX4(k,m,2,2,numEq,numComp,8)]+=8.*wD0;
//                                         EM_S[INDEX4(k,m,2,3,numEq,numComp,8)]+=4.*wD0;
//                                         EM_S[INDEX4(k,m,2,4,numEq,numComp,8)]+=2.*wD0;
//                                         EM_S[INDEX4(k,m,2,5,numEq,numComp,8)]+=   wD0;
//                                         EM_S[INDEX4(k,m,2,6,numEq,numComp,8)]+=4.*wD0;
//                                         EM_S[INDEX4(k,m,2,7,numEq,numComp,8)]+=2.*wD0;
//                                         EM_S[INDEX4(k,m,3,0,numEq,numComp,8)]+=2.*wD0;
//                                         EM_S[INDEX4(k,m,3,1,numEq,numComp,8)]+=4.*wD0;
//                                         EM_S[INDEX4(k,m,3,2,numEq,numComp,8)]+=4.*wD0;
//                                         EM_S[INDEX4(k,m,3,3,numEq,numComp,8)]+=8.*wD0;
//                                         EM_S[INDEX4(k,m,3,4,numEq,numComp,8)]+=   wD0;
//                                         EM_S[INDEX4(k,m,3,5,numEq,numComp,8)]+=2.*wD0;
//                                         EM_S[INDEX4(k,m,3,6,numEq,numComp,8)]+=2.*wD0;
//                                         EM_S[INDEX4(k,m,3,7,numEq,numComp,8)]+=4.*wD0;
//                                         EM_S[INDEX4(k,m,4,0,numEq,numComp,8)]+=4.*wD0;
//                                         EM_S[INDEX4(k,m,4,1,numEq,numComp,8)]+=2.*wD0;
//                                         EM_S[INDEX4(k,m,4,2,numEq,numComp,8)]+=2.*wD0;
//                                         EM_S[INDEX4(k,m,4,3,numEq,numComp,8)]+=   wD0;
//                                         EM_S[INDEX4(k,m,4,4,numEq,numComp,8)]+=8.*wD0;
//                                         EM_S[INDEX4(k,m,4,5,numEq,numComp,8)]+=4.*wD0;
//                                         EM_S[INDEX4(k,m,4,6,numEq,numComp,8)]+=4.*wD0;
//                                         EM_S[INDEX4(k,m,4,7,numEq,numComp,8)]+=2.*wD0;
//                                         EM_S[INDEX4(k,m,5,0,numEq,numComp,8)]+=2.*wD0;
//                                         EM_S[INDEX4(k,m,5,1,numEq,numComp,8)]+=4.*wD0;
//                                         EM_S[INDEX4(k,m,5,2,numEq,numComp,8)]+=   wD0;
//                                         EM_S[INDEX4(k,m,5,3,numEq,numComp,8)]+=2.*wD0;
//                                         EM_S[INDEX4(k,m,5,4,numEq,numComp,8)]+=4.*wD0;
//                                         EM_S[INDEX4(k,m,5,5,numEq,numComp,8)]+=8.*wD0;
//                                         EM_S[INDEX4(k,m,5,6,numEq,numComp,8)]+=2.*wD0;
//                                         EM_S[INDEX4(k,m,5,7,numEq,numComp,8)]+=4.*wD0;
//                                         EM_S[INDEX4(k,m,6,0,numEq,numComp,8)]+=2.*wD0;
//                                         EM_S[INDEX4(k,m,6,1,numEq,numComp,8)]+=   wD0;
//                                         EM_S[INDEX4(k,m,6,2,numEq,numComp,8)]+=4.*wD0;
//                                         EM_S[INDEX4(k,m,6,3,numEq,numComp,8)]+=2.*wD0;
//                                         EM_S[INDEX4(k,m,6,4,numEq,numComp,8)]+=4.*wD0;
//                                         EM_S[INDEX4(k,m,6,5,numEq,numComp,8)]+=2.*wD0;
//                                         EM_S[INDEX4(k,m,6,6,numEq,numComp,8)]+=8.*wD0;
//                                         EM_S[INDEX4(k,m,6,7,numEq,numComp,8)]+=4.*wD0;
//                                         EM_S[INDEX4(k,m,7,0,numEq,numComp,8)]+=   wD0;
//                                         EM_S[INDEX4(k,m,7,1,numEq,numComp,8)]+=2.*wD0;
//                                         EM_S[INDEX4(k,m,7,2,numEq,numComp,8)]+=2.*wD0;
//                                         EM_S[INDEX4(k,m,7,3,numEq,numComp,8)]+=4.*wD0;
//                                         EM_S[INDEX4(k,m,7,4,numEq,numComp,8)]+=2.*wD0;
//                                         EM_S[INDEX4(k,m,7,5,numEq,numComp,8)]+=4.*wD0;
//                                         EM_S[INDEX4(k,m,7,6,numEq,numComp,8)]+=4.*wD0;
//                                         EM_S[INDEX4(k,m,7,7,numEq,numComp,8)]+=8.*wD0;
//                                     }
//                                 }
//                             }
//                         }
//                         ///////////////
//                         // process X //
//                         ///////////////
//                         if (!X.isEmpty()) {
//                             const Scalar* X_p = X.getSampleDataRO(id, zero);
//                             if (X.actsExpanded()) {
//                                 for (index_t k=0; k<numEq; k++) {
//                                     const Scalar X_0_0 = X_p[INDEX3(k,0,0,numEq,3)];
//                                     const Scalar X_1_0 = X_p[INDEX3(k,1,0,numEq,3)];
//                                     const Scalar X_2_0 = X_p[INDEX3(k,2,0,numEq,3)];
//                                     const Scalar X_0_1 = X_p[INDEX3(k,0,1,numEq,3)];
//                                     const Scalar X_1_1 = X_p[INDEX3(k,1,1,numEq,3)];
//                                     const Scalar X_2_1 = X_p[INDEX3(k,2,1,numEq,3)];
//                                     const Scalar X_0_2 = X_p[INDEX3(k,0,2,numEq,3)];
//                                     const Scalar X_1_2 = X_p[INDEX3(k,1,2,numEq,3)];
//                                     const Scalar X_2_2 = X_p[INDEX3(k,2,2,numEq,3)];
//                                     const Scalar X_0_3 = X_p[INDEX3(k,0,3,numEq,3)];
//                                     const Scalar X_1_3 = X_p[INDEX3(k,1,3,numEq,3)];
//                                     const Scalar X_2_3 = X_p[INDEX3(k,2,3,numEq,3)];
//                                     const Scalar X_0_4 = X_p[INDEX3(k,0,4,numEq,3)];
//                                     const Scalar X_1_4 = X_p[INDEX3(k,1,4,numEq,3)];
//                                     const Scalar X_2_4 = X_p[INDEX3(k,2,4,numEq,3)];
//                                     const Scalar X_0_5 = X_p[INDEX3(k,0,5,numEq,3)];
//                                     const Scalar X_1_5 = X_p[INDEX3(k,1,5,numEq,3)];
//                                     const Scalar X_2_5 = X_p[INDEX3(k,2,5,numEq,3)];
//                                     const Scalar X_0_6 = X_p[INDEX3(k,0,6,numEq,3)];
//                                     const Scalar X_1_6 = X_p[INDEX3(k,1,6,numEq,3)];
//                                     const Scalar X_2_6 = X_p[INDEX3(k,2,6,numEq,3)];
//                                     const Scalar X_0_7 = X_p[INDEX3(k,0,7,numEq,3)];
//                                     const Scalar X_1_7 = X_p[INDEX3(k,1,7,numEq,3)];
//                                     const Scalar X_2_7 = X_p[INDEX3(k,2,7,numEq,3)];
//                                     const Scalar tmp0 = w72*(X_0_6 + X_0_7);
//                                     const Scalar tmp1 = w66*(X_2_0 + X_2_4);
//                                     const Scalar tmp2 = w64*(X_0_0 + X_0_1);
//                                     const Scalar tmp3 = w68*(X_2_1 + X_2_2 + X_2_5 + X_2_6);
//                                     const Scalar tmp4 = w65*(X_1_0 + X_1_2);
//                                     const Scalar tmp5 = w70*(X_2_3 + X_2_7);
//                                     const Scalar tmp6 = w67*(X_1_1 + X_1_3 + X_1_4 + X_1_6);
//                                     const Scalar tmp7 = w71*(X_1_5 + X_1_7);
//                                     const Scalar tmp8 = w69*(X_0_2 + X_0_3 + X_0_4 + X_0_5);
//                                     const Scalar tmp9 = w72*(-X_0_6 - X_0_7);
//                                     const Scalar tmp10 = w66*(X_2_1 + X_2_5);
//                                     const Scalar tmp11 = w64*(-X_0_0 - X_0_1);
//                                     const Scalar tmp12 = w68*(X_2_0 + X_2_3 + X_2_4 + X_2_7);
//                                     const Scalar tmp13 = w65*(X_1_1 + X_1_3);
//                                     const Scalar tmp14 = w70*(X_2_2 + X_2_6);
//                                     const Scalar tmp15 = w67*(X_1_0 + X_1_2 + X_1_5 + X_1_7);
//                                     const Scalar tmp16 = w71*(X_1_4 + X_1_6);
//                                     const Scalar tmp17 = w69*(-X_0_2 - X_0_3 - X_0_4 - X_0_5);
//                                     const Scalar tmp18 = w72*(X_0_4 + X_0_5);
//                                     const Scalar tmp19 = w66*(X_2_2 + X_2_6);
//                                     const Scalar tmp20 = w64*(X_0_2 + X_0_3);
//                                     const Scalar tmp21 = w65*(-X_1_0 - X_1_2);
//                                     const Scalar tmp22 = w70*(X_2_1 + X_2_5);
//                                     const Scalar tmp23 = w67*(-X_1_1 - X_1_3 - X_1_4 - X_1_6);
//                                     const Scalar tmp24 = w71*(-X_1_5 - X_1_7);
//                                     const Scalar tmp25 = w69*(X_0_0 + X_0_1 + X_0_6 + X_0_7);
//                                     const Scalar tmp26 = w72*(-X_0_4 - X_0_5);
//                                     const Scalar tmp27 = w66*(X_2_3 + X_2_7);
//                                     const Scalar tmp28 = w64*(-X_0_2 - X_0_3);
//                                     const Scalar tmp29 = w65*(-X_1_1 - X_1_3);
//                                     const Scalar tmp30 = w70*(X_2_0 + X_2_4);
//                                     const Scalar tmp31 = w67*(-X_1_0 - X_1_2 - X_1_5 - X_1_7);
//                                     const Scalar tmp32 = w71*(-X_1_4 - X_1_6);
//                                     const Scalar tmp33 = w69*(-X_0_0 - X_0_1 - X_0_6 - X_0_7);
//                                     const Scalar tmp34 = w72*(X_0_2 + X_0_3);
//                                     const Scalar tmp35 = w66*(-X_2_0 - X_2_4);
//                                     const Scalar tmp36 = w64*(X_0_4 + X_0_5);
//                                     const Scalar tmp37 = w68*(-X_2_1 - X_2_2 - X_2_5 - X_2_6);
//                                     const Scalar tmp38 = w65*(X_1_4 + X_1_6);
//                                     const Scalar tmp39 = w70*(-X_2_3 - X_2_7);
//                                     const Scalar tmp40 = w71*(X_1_1 + X_1_3);
//                                     const Scalar tmp41 = w72*(-X_0_2 - X_0_3);
//                                     const Scalar tmp42 = w66*(-X_2_1 - X_2_5);
//                                     const Scalar tmp43 = w64*(-X_0_4 - X_0_5);
//                                     const Scalar tmp44 = w68*(-X_2_0 - X_2_3 - X_2_4 - X_2_7);
//                                     const Scalar tmp45 = w65*(X_1_5 + X_1_7);
//                                     const Scalar tmp46 = w70*(-X_2_2 - X_2_6);
//                                     const Scalar tmp47 = w71*(X_1_0 + X_1_2);
//                                     const Scalar tmp48 = w72*(X_0_0 + X_0_1);
//                                     const Scalar tmp49 = w66*(-X_2_2 - X_2_6);
//                                     const Scalar tmp50 = w64*(X_0_6 + X_0_7);
//                                     const Scalar tmp51 = w65*(-X_1_4 - X_1_6);
//                                     const Scalar tmp52 = w70*(-X_2_1 - X_2_5);
//                                     const Scalar tmp53 = w71*(-X_1_1 - X_1_3);
//                                     const Scalar tmp54 = w72*(-X_0_0 - X_0_1);
//                                     const Scalar tmp55 = w66*(-X_2_3 - X_2_7);
//                                     const Scalar tmp56 = w64*(-X_0_6 - X_0_7);
//                                     const Scalar tmp57 = w65*(-X_1_5 - X_1_7);
//                                     const Scalar tmp58 = w70*(-X_2_0 - X_2_4);
//                                     const Scalar tmp59 = w71*(-X_1_0 - X_1_2);
//                                     EM_F[INDEX2(k,0,numEq)]+=tmp0 + tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6 + tmp7 + tmp8;
//                                     EM_F[INDEX2(k,1,numEq)]+=tmp10 + tmp11 + tmp12 + tmp13 + tmp14 + tmp15 + tmp16 + tmp17 + tmp9;
//                                     EM_F[INDEX2(k,2,numEq)]+=tmp12 + tmp18 + tmp19 + tmp20 + tmp21 + tmp22 + tmp23 + tmp24 + tmp25;
//                                     EM_F[INDEX2(k,3,numEq)]+=tmp26 + tmp27 + tmp28 + tmp29 + tmp3 + tmp30 + tmp31 + tmp32 + tmp33;
//                                     EM_F[INDEX2(k,4,numEq)]+=tmp15 + tmp25 + tmp34 + tmp35 + tmp36 + tmp37 + tmp38 + tmp39 + tmp40;
//                                     EM_F[INDEX2(k,5,numEq)]+=tmp33 + tmp41 + tmp42 + tmp43 + tmp44 + tmp45 + tmp46 + tmp47 + tmp6;
//                                     EM_F[INDEX2(k,6,numEq)]+=tmp31 + tmp44 + tmp48 + tmp49 + tmp50 + tmp51 + tmp52 + tmp53 + tmp8;
//                                     EM_F[INDEX2(k,7,numEq)]+=tmp17 + tmp23 + tmp37 + tmp54 + tmp55 + tmp56 + tmp57 + tmp58 + tmp59;
//                                 }
//                             } else { // constant data
//                                 for (index_t k=0; k<numEq; k++) {
//                                     const Scalar wX0 = 18.*X_p[INDEX2(k, 0, numEq)]*w55;
//                                     const Scalar wX1 = 18.*X_p[INDEX2(k, 1, numEq)]*w56;
//                                     const Scalar wX2 = 18.*X_p[INDEX2(k, 2, numEq)]*w54;
//                                     EM_F[INDEX2(k,0,numEq)]+= wX0 + wX1 + wX2;
//                                     EM_F[INDEX2(k,1,numEq)]+=-wX0 + wX1 + wX2;
//                                     EM_F[INDEX2(k,2,numEq)]+= wX0 - wX1 + wX2;
//                                     EM_F[INDEX2(k,3,numEq)]+=-wX0 - wX1 + wX2;
//                                     EM_F[INDEX2(k,4,numEq)]+= wX0 + wX1 - wX2;
//                                     EM_F[INDEX2(k,5,numEq)]+=-wX0 + wX1 - wX2;
//                                     EM_F[INDEX2(k,6,numEq)]+= wX0 - wX1 - wX2;
//                                     EM_F[INDEX2(k,7,numEq)]+=-wX0 - wX1 - wX2;
//                                 }
//                             }
//                         }
//                         ///////////////
//                         // process Y //
//                         ///////////////
//                         if (!Y.isEmpty()) {
//                             const Scalar* Y_p = Y.getSampleDataRO(id, zero);
//                             if (Y.actsExpanded()) {
//                                 for (index_t k=0; k<numEq; k++) {
//                                     const Scalar Y_0 = Y_p[INDEX2(k, 0, numEq)];
//                                     const Scalar Y_1 = Y_p[INDEX2(k, 1, numEq)];
//                                     const Scalar Y_2 = Y_p[INDEX2(k, 2, numEq)];
//                                     const Scalar Y_3 = Y_p[INDEX2(k, 3, numEq)];
//                                     const Scalar Y_4 = Y_p[INDEX2(k, 4, numEq)];
//                                     const Scalar Y_5 = Y_p[INDEX2(k, 5, numEq)];
//                                     const Scalar Y_6 = Y_p[INDEX2(k, 6, numEq)];
//                                     const Scalar Y_7 = Y_p[INDEX2(k, 7, numEq)];
//                                     const Scalar tmp0 = w76*(Y_3 + Y_5 + Y_6);
//                                     const Scalar tmp1 = w75*(Y_1 + Y_2 + Y_4);
//                                     const Scalar tmp2 = w76*(Y_2 + Y_4 + Y_7);
//                                     const Scalar tmp3 = w75*(Y_0 + Y_3 + Y_5);
//                                     const Scalar tmp4 = w76*(Y_1 + Y_4 + Y_7);
//                                     const Scalar tmp5 = w75*(Y_0 + Y_3 + Y_6);
//                                     const Scalar tmp6 = w76*(Y_0 + Y_5 + Y_6);
//                                     const Scalar tmp7 = w75*(Y_1 + Y_2 + Y_7);
//                                     const Scalar tmp8 = w76*(Y_1 + Y_2 + Y_7);
//                                     const Scalar tmp9 = w75*(Y_0 + Y_5 + Y_6);
//                                     const Scalar tmp10 = w76*(Y_0 + Y_3 + Y_6);
//                                     const Scalar tmp11 = w75*(Y_1 + Y_4 + Y_7);
//                                     const Scalar tmp12 = w76*(Y_0 + Y_3 + Y_5);
//                                     const Scalar tmp13 = w75*(Y_2 + Y_4 + Y_7);
//                                     const Scalar tmp14 = w76*(Y_1 + Y_2 + Y_4);
//                                     const Scalar tmp15 = w75*(Y_3 + Y_5 + Y_6);
//                                     EM_F[INDEX2(k,0,numEq)]+=Y_0*w74 + Y_7*w77 + tmp0 + tmp1;
//                                     EM_F[INDEX2(k,1,numEq)]+=Y_1*w74 + Y_6*w77 + tmp2 + tmp3;
//                                     EM_F[INDEX2(k,2,numEq)]+=Y_2*w74 + Y_5*w77 + tmp4 + tmp5;
//                                     EM_F[INDEX2(k,3,numEq)]+=Y_3*w74 + Y_4*w77 + tmp6 + tmp7;
//                                     EM_F[INDEX2(k,4,numEq)]+=Y_3*w77 + Y_4*w74 + tmp8 + tmp9;
//                                     EM_F[INDEX2(k,5,numEq)]+=Y_2*w77 + Y_5*w74 + tmp10 + tmp11;
//                                     EM_F[INDEX2(k,6,numEq)]+=Y_1*w77 + Y_6*w74 + tmp12 + tmp13;
//                                     EM_F[INDEX2(k,7,numEq)]+=Y_0*w77 + Y_7*w74 + tmp14 + tmp15;
//                                 }
//                             } else { // constant data
//                                 for (index_t k=0; k<numEq; k++) {
//                                     EM_F[INDEX2(k,0,numEq)]+=216.*Y_p[k]*w58;
//                                     EM_F[INDEX2(k,1,numEq)]+=216.*Y_p[k]*w58;
//                                     EM_F[INDEX2(k,2,numEq)]+=216.*Y_p[k]*w58;
//                                     EM_F[INDEX2(k,3,numEq)]+=216.*Y_p[k]*w58;
//                                     EM_F[INDEX2(k,4,numEq)]+=216.*Y_p[k]*w58;
//                                     EM_F[INDEX2(k,5,numEq)]+=216.*Y_p[k]*w58;
//                                     EM_F[INDEX2(k,6,numEq)]+=216.*Y_p[k]*w58;
//                                     EM_F[INDEX2(k,7,numEq)]+=216.*Y_p[k]*w58;
//                                 }
//                             }
//                         }

//                         // add to matrix (if add_EM_S) and RHS (if add_EM_F)
//                         const index_t firstNode=m_NN[0]*m_NN[1]*k2+m_NN[0]*k1+k0;
//                         domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F,
//                                 add_EM_S, add_EM_F, firstNode, numEq, numComp);
//                     } // end k0 loop
//                 } // end k1 loop
//             } // end k2 loop
//         } // end of colouring
//     } // end of parallel region
}

/****************************************************************************/
// PDE SYSTEM BOUNDARY
/****************************************************************************/

template<class Scalar>
void DefaultAssembler3D<Scalar>::assemblePDEBoundarySystem(
                                        AbstractSystemMatrix* mat, Data& rhs,
                                        const Data& d, const Data& y) const
{
    // TODO

//     dim_t numEq, numComp;
//     if (!mat)
//         numEq=numComp=(rhs.isEmpty() ? 1 : rhs.getDataPointSize());
//     else {
//         numEq=mat->getRowBlockSize();
//         numComp=mat->getColumnBlockSize();
//     }
//     const double SQRT3 = 1.73205080756887719318;
//     const double w12 = m_dx[0]*m_dx[1]/144;
//     const double w10 = w12*(-SQRT3 + 2);
//     const double w11 = w12*(SQRT3 + 2);
//     const double w13 = w12*(-4*SQRT3 + 7);
//     const double w14 = w12*(4*SQRT3 + 7);
//     const double w7 = m_dx[0]*m_dx[2]/144;
//     const double w5 = w7*(-SQRT3 + 2);
//     const double w6 = w7*(SQRT3 + 2);
//     const double w8 = w7*(-4*SQRT3 + 7);
//     const double w9 = w7*(4*SQRT3 + 7);
//     const double w2 = m_dx[1]*m_dx[2]/144;
//     const double w0 = w2*(-SQRT3 + 2);
//     const double w1 = w2*(SQRT3 + 2);
//     const double w3 = w2*(-4*SQRT3 + 7);
//     const double w4 = w2*(4*SQRT3 + 7);
//     const dim_t NE0 = m_NE[0];
//     const dim_t NE1 = m_NE[1];
//     const dim_t NE2 = m_NE[2];
//     const bool add_EM_S = !d.isEmpty();
//     const bool add_EM_F = !y.isEmpty();
//     const Scalar zero = static_cast<Scalar>(0);
//     rhs.requireWrite();

// #pragma omp parallel
//     {
//         vector<Scalar> EM_S(8*8*numEq*numComp, zero);
//         vector<Scalar> EM_F(8*numEq, zero);

//         if (domain->m_faceOffset[0] > -1) {
//             if (add_EM_S)
//                 fill(EM_S.begin(), EM_S.end(), zero);
//             if (add_EM_F)
//                 fill(EM_F.begin(), EM_F.end(), zero);

//             for (index_t k2_0=0; k2_0<2; k2_0++) { // colouring
// #pragma omp for
//                 for (index_t k2=k2_0; k2<NE2; k2+=2) {
//                     for (index_t k1=0; k1<NE1; ++k1) {
//                         const index_t e = INDEX2(k1,k2,NE1);
//                         ///////////////
//                         // process d //
//                         ///////////////
//                         if (add_EM_S) {
//                             const Scalar* d_p = d.getSampleDataRO(id, zero);
//                             if (d.actsExpanded()) {
//                                 for (index_t k=0; k<numEq; k++) {
//                                     for (index_t m=0; m<numComp; m++) {
//                                         const Scalar d_0 = d_p[INDEX3(k,m,0,numEq,numComp)];
//                                         const Scalar d_1 = d_p[INDEX3(k,m,1,numEq,numComp)];
//                                         const Scalar d_2 = d_p[INDEX3(k,m,2,numEq,numComp)];
//                                         const Scalar d_3 = d_p[INDEX3(k,m,3,numEq,numComp)];
//                                         const Scalar tmp0 = w0*(d_0 + d_1);
//                                         const Scalar tmp1 = w1*(d_2 + d_3);
//                                         const Scalar tmp2 = w0*(d_0 + d_2);
//                                         const Scalar tmp3 = w1*(d_1 + d_3);
//                                         const Scalar tmp4 = w0*(d_1 + d_3);
//                                         const Scalar tmp5 = w1*(d_0 + d_2);
//                                         const Scalar tmp6 = w0*(d_2 + d_3);
//                                         const Scalar tmp7 = w1*(d_0 + d_1);
//                                         const Scalar tmp8 = w2*(d_0 + d_3);
//                                         const Scalar tmp9 = w2*(d_1 + d_2);
//                                         const Scalar tmp10 = w2*(d_0 + d_1 + d_2 + d_3);
//                                         EM_S[INDEX4(k,m,0,0,numEq,numComp,8)] = d_0*w4 + d_3*w3 + tmp9;
//                                         EM_S[INDEX4(k,m,0,2,numEq,numComp,8)] = tmp6 + tmp7;
//                                         EM_S[INDEX4(k,m,0,4,numEq,numComp,8)] = tmp4 + tmp5;
//                                         EM_S[INDEX4(k,m,0,6,numEq,numComp,8)] = tmp10;
//                                         EM_S[INDEX4(k,m,2,0,numEq,numComp,8)] = tmp6 + tmp7;
//                                         EM_S[INDEX4(k,m,2,2,numEq,numComp,8)] = d_1*w4 + d_2*w3 + tmp8;
//                                         EM_S[INDEX4(k,m,2,4,numEq,numComp,8)] = tmp10;
//                                         EM_S[INDEX4(k,m,2,6,numEq,numComp,8)] = tmp2 + tmp3;
//                                         EM_S[INDEX4(k,m,4,0,numEq,numComp,8)] = tmp4 + tmp5;
//                                         EM_S[INDEX4(k,m,4,2,numEq,numComp,8)] = tmp10;
//                                         EM_S[INDEX4(k,m,4,4,numEq,numComp,8)] = d_1*w3 + d_2*w4 + tmp8;
//                                         EM_S[INDEX4(k,m,4,6,numEq,numComp,8)] = tmp0 + tmp1;
//                                         EM_S[INDEX4(k,m,6,0,numEq,numComp,8)] = tmp10;
//                                         EM_S[INDEX4(k,m,6,2,numEq,numComp,8)] = tmp2 + tmp3;
//                                         EM_S[INDEX4(k,m,6,4,numEq,numComp,8)] = tmp0 + tmp1;
//                                         EM_S[INDEX4(k,m,6,6,numEq,numComp,8)] = d_0*w3 + d_3*w4 + tmp9;
//                                     }
//                                  }
//                             } else { // constant data
//                                 for (index_t k=0; k<numEq; k++) {
//                                     for (index_t m=0; m<numComp; m++) {
//                                         const Scalar wd0 = 4.*d_p[INDEX2(k, m, numEq)]*w2;
//                                         EM_S[INDEX4(k,m,0,0,numEq,numComp,8)] = 4.*wd0;
//                                         EM_S[INDEX4(k,m,0,2,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,0,4,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,0,6,numEq,numComp,8)] =    wd0;
//                                         EM_S[INDEX4(k,m,2,0,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,2,2,numEq,numComp,8)] = 4.*wd0;
//                                         EM_S[INDEX4(k,m,2,4,numEq,numComp,8)] =    wd0;
//                                         EM_S[INDEX4(k,m,2,6,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,4,0,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,4,2,numEq,numComp,8)] =    wd0;
//                                         EM_S[INDEX4(k,m,4,4,numEq,numComp,8)] = 4.*wd0;
//                                         EM_S[INDEX4(k,m,4,6,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,6,0,numEq,numComp,8)] =    wd0;
//                                         EM_S[INDEX4(k,m,6,2,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,6,4,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,6,6,numEq,numComp,8)] = 4.*wd0;
//                                     }
//                                 }
//                             }
//                         }
//                         ///////////////
//                         // process y //
//                         ///////////////
//                         if (add_EM_F) {
//                             const Scalar* y_p = y.getSampleDataRO(id, zero);
//                             if (y.actsExpanded()) {
//                                 for (index_t k=0; k<numEq; k++) {
//                                     const Scalar y_0 = y_p[INDEX2(k, 0, numEq)];
//                                     const Scalar y_1 = y_p[INDEX2(k, 1, numEq)];
//                                     const Scalar y_2 = y_p[INDEX2(k, 2, numEq)];
//                                     const Scalar y_3 = y_p[INDEX2(k, 3, numEq)];
//                                     const Scalar tmp0 = 6.*w2*(y_1 + y_2);
//                                     const Scalar tmp1 = 6.*w2*(y_0 + y_3);
//                                     EM_F[INDEX2(k,0,numEq)] = tmp0 + 6.*w0*y_3 + 6.*w1*y_0;
//                                     EM_F[INDEX2(k,2,numEq)] = tmp1 + 6.*w0*y_2 + 6.*w1*y_1;
//                                     EM_F[INDEX2(k,4,numEq)] = tmp1 + 6.*w0*y_1 + 6.*w1*y_2;
//                                     EM_F[INDEX2(k,6,numEq)] = tmp0 + 6.*w0*y_0 + 6.*w1*y_3;
//                                 }
//                             } else { // constant data
//                                 for (index_t k=0; k<numEq; k++) {
//                                     EM_F[INDEX2(k,0,numEq)] = 36.*w2*y_p[k];
//                                     EM_F[INDEX2(k,2,numEq)] = 36.*w2*y_p[k];
//                                     EM_F[INDEX2(k,4,numEq)] = 36.*w2*y_p[k];
//                                     EM_F[INDEX2(k,6,numEq)] = 36.*w2*y_p[k];
//                                 }
//                             }
//                         }
//                         const index_t firstNode=m_NN[0]*m_NN[1]*k2+m_NN[0]*k1;
//                         domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F,
//                                 add_EM_S, add_EM_F, firstNode, numEq, numComp);
//                     } // k1 loop
//                 } // k2 loop
//             } // colouring
//         } // face 0

//         if (domain->m_faceOffset[1] > -1) {
//             if (add_EM_S)
//                 fill(EM_S.begin(), EM_S.end(), zero);
//             if (add_EM_F)
//                 fill(EM_F.begin(), EM_F.end(), zero);

//             for (index_t k2_0=0; k2_0<2; k2_0++) { // colouring
// #pragma omp for
//                 for (index_t k2=k2_0; k2<NE2; k2+=2) {
//                     for (index_t k1=0; k1<NE1; ++k1) {
//                         const index_t e = domain->m_faceOffset[1]+INDEX2(k1,k2,NE1);
//                         ///////////////
//                         // process d //
//                         ///////////////
//                         if (add_EM_S) {
//                             const Scalar* d_p = d.getSampleDataRO(id, zero);
//                             if (d.actsExpanded()) {
//                                 for (index_t k=0; k<numEq; k++) {
//                                     for (index_t m=0; m<numComp; m++) {
//                                         const Scalar d_0 = d_p[INDEX3(k,m,0,numEq,numComp)];
//                                         const Scalar d_1 = d_p[INDEX3(k,m,1,numEq,numComp)];
//                                         const Scalar d_2 = d_p[INDEX3(k,m,2,numEq,numComp)];
//                                         const Scalar d_3 = d_p[INDEX3(k,m,3,numEq,numComp)];
//                                         const Scalar tmp0 = w0*(d_0 + d_2);
//                                         const Scalar tmp1 = w1*(d_1 + d_3);
//                                         const Scalar tmp2 = w0*(d_2 + d_3);
//                                         const Scalar tmp3 = w1*(d_0 + d_1);
//                                         const Scalar tmp4 = w0*(d_1 + d_3);
//                                         const Scalar tmp5 = w1*(d_0 + d_2);
//                                         const Scalar tmp6 = w2*(d_0 + d_3);
//                                         const Scalar tmp7 = w2*(d_1 + d_2);
//                                         const Scalar tmp8 = w0*(d_0 + d_1);
//                                         const Scalar tmp9 = w1*(d_2 + d_3);
//                                         const Scalar tmp10 = w2*(d_0 + d_1 + d_2 + d_3);
//                                         EM_S[INDEX4(k,m,1,1,numEq,numComp,8)] = d_0*w4 + d_3*w3 + tmp7;
//                                         EM_S[INDEX4(k,m,1,3,numEq,numComp,8)] = tmp2 + tmp3;
//                                         EM_S[INDEX4(k,m,1,5,numEq,numComp,8)] = tmp4 + tmp5;
//                                         EM_S[INDEX4(k,m,1,7,numEq,numComp,8)] = tmp10;
//                                         EM_S[INDEX4(k,m,3,1,numEq,numComp,8)] = tmp2 + tmp3;
//                                         EM_S[INDEX4(k,m,3,3,numEq,numComp,8)] = d_1*w4 + d_2*w3 + tmp6;
//                                         EM_S[INDEX4(k,m,3,5,numEq,numComp,8)] = tmp10;
//                                         EM_S[INDEX4(k,m,3,7,numEq,numComp,8)] = tmp0 + tmp1;
//                                         EM_S[INDEX4(k,m,5,1,numEq,numComp,8)] = tmp4 + tmp5;
//                                         EM_S[INDEX4(k,m,5,3,numEq,numComp,8)] = tmp10;
//                                         EM_S[INDEX4(k,m,5,5,numEq,numComp,8)] = d_1*w3 + d_2*w4 + tmp6;
//                                         EM_S[INDEX4(k,m,5,7,numEq,numComp,8)] = tmp8 + tmp9;
//                                         EM_S[INDEX4(k,m,7,1,numEq,numComp,8)] = tmp10;
//                                         EM_S[INDEX4(k,m,7,3,numEq,numComp,8)] = tmp0 + tmp1;
//                                         EM_S[INDEX4(k,m,7,5,numEq,numComp,8)] = tmp8 + tmp9;
//                                         EM_S[INDEX4(k,m,7,7,numEq,numComp,8)] = d_0*w3 + d_3*w4 + tmp7;
//                                     }
//                                  }
//                             } else { // constant data
//                                 for (index_t k=0; k<numEq; k++) {
//                                     for (index_t m=0; m<numComp; m++) {
//                                         const Scalar wd0 = 4.*d_p[INDEX2(k, m, numEq)]*w2;
//                                         EM_S[INDEX4(k,m,1,1,numEq,numComp,8)] = 4.*wd0;
//                                         EM_S[INDEX4(k,m,1,3,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,1,5,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,1,7,numEq,numComp,8)] =    wd0;
//                                         EM_S[INDEX4(k,m,3,1,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,3,3,numEq,numComp,8)] = 4.*wd0;
//                                         EM_S[INDEX4(k,m,3,5,numEq,numComp,8)] =    wd0;
//                                         EM_S[INDEX4(k,m,3,7,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,5,1,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,5,3,numEq,numComp,8)] =    wd0;
//                                         EM_S[INDEX4(k,m,5,5,numEq,numComp,8)] = 4.*wd0;
//                                         EM_S[INDEX4(k,m,5,7,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,7,1,numEq,numComp,8)] =    wd0;
//                                         EM_S[INDEX4(k,m,7,3,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,7,5,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,7,7,numEq,numComp,8)] = 4.*wd0;
//                                     }
//                                 }
//                             }
//                         }
//                         ///////////////
//                         // process y //
//                         ///////////////
//                         if (add_EM_F) {
//                             const Scalar* y_p = y.getSampleDataRO(id, zero);
//                             if (y.actsExpanded()) {
//                                 for (index_t k=0; k<numEq; k++) {
//                                     const Scalar y_0 = y_p[INDEX2(k, 0, numEq)];
//                                     const Scalar y_1 = y_p[INDEX2(k, 1, numEq)];
//                                     const Scalar y_2 = y_p[INDEX2(k, 2, numEq)];
//                                     const Scalar y_3 = y_p[INDEX2(k, 3, numEq)];
//                                     const Scalar tmp0 = 6.*w2*(y_1 + y_2);
//                                     const Scalar tmp1 = 6.*w2*(y_0 + y_3);
//                                     EM_F[INDEX2(k,1,numEq)] = tmp0 + 6.*w0*y_3 + 6.*w1*y_0;
//                                     EM_F[INDEX2(k,3,numEq)] = tmp1 + 6.*w0*y_2 + 6.*w1*y_1;
//                                     EM_F[INDEX2(k,5,numEq)] = tmp1 + 6.*w0*y_1 + 6.*w1*y_2;
//                                     EM_F[INDEX2(k,7,numEq)] = tmp0 + 6.*w0*y_0 + 6.*w1*y_3;
//                                 }
//                             } else { // constant data
//                                 for (index_t k=0; k<numEq; k++) {
//                                     EM_F[INDEX2(k,1,numEq)] = 36.*w2*y_p[k];
//                                     EM_F[INDEX2(k,3,numEq)] = 36.*w2*y_p[k];
//                                     EM_F[INDEX2(k,5,numEq)] = 36.*w2*y_p[k];
//                                     EM_F[INDEX2(k,7,numEq)] = 36.*w2*y_p[k];
//                                 }
//                             }
//                         }
//                         const index_t firstNode=m_NN[0]*m_NN[1]*k2+m_NN[0]*(k1+1)-2;
//                         domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F,
//                                 add_EM_S, add_EM_F, firstNode, numEq, numComp);
//                     } // k1 loop
//                 } // k2 loop
//             } // colouring
//         } // face 1

//         if (domain->m_faceOffset[2] > -1) {
//             if (add_EM_S)
//                 fill(EM_S.begin(), EM_S.end(), zero);
//             if (add_EM_F)
//                 fill(EM_F.begin(), EM_F.end(), zero);

//             for (index_t k2_0=0; k2_0<2; k2_0++) { // colouring
// #pragma omp for
//                 for (index_t k2=k2_0; k2<NE2; k2+=2) {
//                     for (index_t k0=0; k0<NE0; ++k0) {
//                         const index_t e = domain->m_faceOffset[2]+INDEX2(k0,k2,NE0);
//                         ///////////////
//                         // process d //
//                         ///////////////
//                         if (add_EM_S) {
//                             const Scalar* d_p = d.getSampleDataRO(id, zero);
//                             if (d.actsExpanded()) {
//                                 for (index_t k=0; k<numEq; k++) {
//                                     for (index_t m=0; m<numComp; m++) {
//                                         const Scalar d_0 = d_p[INDEX3(k,m,0,numEq,numComp)];
//                                         const Scalar d_1 = d_p[INDEX3(k,m,1,numEq,numComp)];
//                                         const Scalar d_2 = d_p[INDEX3(k,m,2,numEq,numComp)];
//                                         const Scalar d_3 = d_p[INDEX3(k,m,3,numEq,numComp)];
//                                         const Scalar tmp0 = w5*(d_0 + d_1);
//                                         const Scalar tmp1 = w6*(d_2 + d_3);
//                                         const Scalar tmp2 = w5*(d_0 + d_2);
//                                         const Scalar tmp3 = w6*(d_1 + d_3);
//                                         const Scalar tmp4 = w5*(d_1 + d_3);
//                                         const Scalar tmp5 = w6*(d_0 + d_2);
//                                         const Scalar tmp6 = w7*(d_0 + d_3);
//                                         const Scalar tmp7 = w7*(d_0 + d_1 + d_2 + d_3);
//                                         const Scalar tmp8 = w7*(d_1 + d_2);
//                                         const Scalar tmp9 = w5*(d_2 + d_3);
//                                         const Scalar tmp10 = w6*(d_0 + d_1);
//                                         EM_S[INDEX4(k,m,0,0,numEq,numComp,8)] = d_0*w9 + d_3*w8 + tmp8;
//                                         EM_S[INDEX4(k,m,0,1,numEq,numComp,8)] = tmp10 + tmp9;
//                                         EM_S[INDEX4(k,m,0,4,numEq,numComp,8)] = tmp4 + tmp5;
//                                         EM_S[INDEX4(k,m,0,5,numEq,numComp,8)] = tmp7;
//                                         EM_S[INDEX4(k,m,1,0,numEq,numComp,8)] = tmp10 + tmp9;
//                                         EM_S[INDEX4(k,m,1,1,numEq,numComp,8)] = d_1*w9 + d_2*w8 + tmp6;
//                                         EM_S[INDEX4(k,m,1,4,numEq,numComp,8)] = tmp7;
//                                         EM_S[INDEX4(k,m,1,5,numEq,numComp,8)] = tmp2 + tmp3;
//                                         EM_S[INDEX4(k,m,4,0,numEq,numComp,8)] = tmp4 + tmp5;
//                                         EM_S[INDEX4(k,m,4,1,numEq,numComp,8)] = tmp7;
//                                         EM_S[INDEX4(k,m,4,4,numEq,numComp,8)] = d_1*w8 + d_2*w9 + tmp6;
//                                         EM_S[INDEX4(k,m,4,5,numEq,numComp,8)] = tmp0 + tmp1;
//                                         EM_S[INDEX4(k,m,5,0,numEq,numComp,8)] = tmp7;
//                                         EM_S[INDEX4(k,m,5,1,numEq,numComp,8)] = tmp2 + tmp3;
//                                         EM_S[INDEX4(k,m,5,4,numEq,numComp,8)] = tmp0 + tmp1;
//                                         EM_S[INDEX4(k,m,5,5,numEq,numComp,8)] = d_0*w8 + d_3*w9 + tmp8;
//                                     }
//                                  }
//                             } else { // constant data
//                                 for (index_t k=0; k<numEq; k++) {
//                                     for (index_t m=0; m<numComp; m++) {
//                                         const Scalar wd0 = 4.*d_p[INDEX2(k, m, numEq)]*w7;
//                                         EM_S[INDEX4(k,m,0,0,numEq,numComp,8)] = 4.*wd0;
//                                         EM_S[INDEX4(k,m,0,1,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,0,4,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,0,5,numEq,numComp,8)] =    wd0;
//                                         EM_S[INDEX4(k,m,1,0,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,1,1,numEq,numComp,8)] = 4.*wd0;
//                                         EM_S[INDEX4(k,m,1,4,numEq,numComp,8)] =    wd0;
//                                         EM_S[INDEX4(k,m,1,5,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,4,0,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,4,1,numEq,numComp,8)] =    wd0;
//                                         EM_S[INDEX4(k,m,4,4,numEq,numComp,8)] = 4.*wd0;
//                                         EM_S[INDEX4(k,m,4,5,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,5,0,numEq,numComp,8)] =    wd0;
//                                         EM_S[INDEX4(k,m,5,1,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,5,4,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,5,5,numEq,numComp,8)] = 4.*wd0;
//                                     }
//                                 }
//                             }
//                         }
//                         ///////////////
//                         // process y //
//                         ///////////////
//                         if (add_EM_F) {
//                             const Scalar* y_p = y.getSampleDataRO(id, zero);
//                             if (y.actsExpanded()) {
//                                 for (index_t k=0; k<numEq; k++) {
//                                     const Scalar y_0 = y_p[INDEX2(k, 0, numEq)];
//                                     const Scalar y_1 = y_p[INDEX2(k, 1, numEq)];
//                                     const Scalar y_2 = y_p[INDEX2(k, 2, numEq)];
//                                     const Scalar y_3 = y_p[INDEX2(k, 3, numEq)];
//                                     const Scalar tmp0 = 6.*w7*(y_1 + y_2);
//                                     const Scalar tmp1 = 6.*w7*(y_0 + y_3);
//                                     EM_F[INDEX2(k,0,numEq)] = tmp0 + 6.*w5*y_3 + 6.*w6*y_0;
//                                     EM_F[INDEX2(k,1,numEq)] = tmp1 + 6.*w5*y_2 + 6.*w6*y_1;
//                                     EM_F[INDEX2(k,4,numEq)] = tmp1 + 6.*w5*y_1 + 6.*w6*y_2;
//                                     EM_F[INDEX2(k,5,numEq)] = tmp0 + 6.*w5*y_0 + 6.*w6*y_3;
//                                 }
//                             } else { // constant data
//                                 for (index_t k=0; k<numEq; k++) {
//                                     EM_F[INDEX2(k,0,numEq)] = 36.*w7*y_p[k];
//                                     EM_F[INDEX2(k,1,numEq)] = 36.*w7*y_p[k];
//                                     EM_F[INDEX2(k,4,numEq)] = 36.*w7*y_p[k];
//                                     EM_F[INDEX2(k,5,numEq)] = 36.*w7*y_p[k];
//                                 }
//                             }
//                         }
//                         const index_t firstNode=m_NN[0]*m_NN[1]*k2+k0;
//                         domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F,
//                                 add_EM_S, add_EM_F, firstNode, numEq, numComp);
//                     } // k0 loop
//                 } // k2 loop
//             } // colouring
//         } // face 2

//         if (domain->m_faceOffset[3] > -1) {
//             if (add_EM_S)
//                 fill(EM_S.begin(), EM_S.end(), zero);
//             if (add_EM_F)
//                 fill(EM_F.begin(), EM_F.end(), zero);

//             for (index_t k2_0=0; k2_0<2; k2_0++) { // colouring
// #pragma omp for
//                 for (index_t k2=k2_0; k2<NE2; k2+=2) {
//                     for (index_t k0=0; k0<NE0; ++k0) {
//                         const index_t e = domain->m_faceOffset[3]+INDEX2(k0,k2,NE0);
//                         ///////////////
//                         // process d //
//                         ///////////////
//                         if (add_EM_S) {
//                             const Scalar* d_p = d.getSampleDataRO(id, zero);
//                             if (d.actsExpanded()) {
//                                 for (index_t k=0; k<numEq; k++) {
//                                     for (index_t m=0; m<numComp; m++) {
//                                         const Scalar d_0 = d_p[INDEX3(k,m,0,numEq,numComp)];
//                                         const Scalar d_1 = d_p[INDEX3(k,m,1,numEq,numComp)];
//                                         const Scalar d_2 = d_p[INDEX3(k,m,2,numEq,numComp)];
//                                         const Scalar d_3 = d_p[INDEX3(k,m,3,numEq,numComp)];
//                                         const Scalar tmp0 = w5*(d_0 + d_2);
//                                         const Scalar tmp1 = w6*(d_1 + d_3);
//                                         const Scalar tmp2 = w5*(d_1 + d_3);
//                                         const Scalar tmp3 = w6*(d_0 + d_2);
//                                         const Scalar tmp4 = w7*(d_0 + d_1 + d_2 + d_3);
//                                         const Scalar tmp5 = w5*(d_0 + d_1);
//                                         const Scalar tmp6 = w6*(d_2 + d_3);
//                                         const Scalar tmp7 = w7*(d_0 + d_3);
//                                         const Scalar tmp8 = w7*(d_1 + d_2);
//                                         const Scalar tmp9 = w5*(d_2 + d_3);
//                                         const Scalar tmp10 = w6*(d_0 + d_1);
//                                         EM_S[INDEX4(k,m,2,2,numEq,numComp,8)] = d_0*w9 + d_3*w8 + tmp8;
//                                         EM_S[INDEX4(k,m,2,3,numEq,numComp,8)] = tmp10 + tmp9;
//                                         EM_S[INDEX4(k,m,2,6,numEq,numComp,8)] = tmp2 + tmp3;
//                                         EM_S[INDEX4(k,m,2,7,numEq,numComp,8)] = tmp4;
//                                         EM_S[INDEX4(k,m,3,2,numEq,numComp,8)] = tmp10 + tmp9;
//                                         EM_S[INDEX4(k,m,3,3,numEq,numComp,8)] = d_1*w9 + d_2*w8 + tmp7;
//                                         EM_S[INDEX4(k,m,3,6,numEq,numComp,8)] = tmp4;
//                                         EM_S[INDEX4(k,m,3,7,numEq,numComp,8)] = tmp0 + tmp1;
//                                         EM_S[INDEX4(k,m,6,2,numEq,numComp,8)] = tmp2 + tmp3;
//                                         EM_S[INDEX4(k,m,6,3,numEq,numComp,8)] = tmp4;
//                                         EM_S[INDEX4(k,m,6,6,numEq,numComp,8)] = d_1*w8 + d_2*w9 + tmp7;
//                                         EM_S[INDEX4(k,m,6,7,numEq,numComp,8)] = tmp5 + tmp6;
//                                         EM_S[INDEX4(k,m,7,2,numEq,numComp,8)] = tmp4;
//                                         EM_S[INDEX4(k,m,7,3,numEq,numComp,8)] = tmp0 + tmp1;
//                                         EM_S[INDEX4(k,m,7,6,numEq,numComp,8)] = tmp5 + tmp6;
//                                         EM_S[INDEX4(k,m,7,7,numEq,numComp,8)] = d_0*w8 + d_3*w9 + tmp8;
//                                     }
//                                  }
//                             } else { // constant data
//                                 for (index_t k=0; k<numEq; k++) {
//                                     for (index_t m=0; m<numComp; m++) {
//                                         const Scalar wd0 = 4.*d_p[INDEX2(k, m, numEq)]*w7;
//                                         EM_S[INDEX4(k,m,2,2,numEq,numComp,8)] = 4.*wd0;
//                                         EM_S[INDEX4(k,m,2,3,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,2,6,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,2,7,numEq,numComp,8)] =   wd0;
//                                         EM_S[INDEX4(k,m,3,2,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,3,3,numEq,numComp,8)] = 4.*wd0;
//                                         EM_S[INDEX4(k,m,3,6,numEq,numComp,8)] =   wd0;
//                                         EM_S[INDEX4(k,m,3,7,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,6,2,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,6,3,numEq,numComp,8)] =   wd0;
//                                         EM_S[INDEX4(k,m,6,6,numEq,numComp,8)] = 4.*wd0;
//                                         EM_S[INDEX4(k,m,6,7,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,7,2,numEq,numComp,8)] =   wd0;
//                                         EM_S[INDEX4(k,m,7,3,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,7,6,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,7,7,numEq,numComp,8)] = 4.*wd0;
//                                     }
//                                 }
//                             }
//                         }
//                         ///////////////
//                         // process y //
//                         ///////////////
//                         if (add_EM_F) {
//                             const Scalar* y_p = y.getSampleDataRO(id, zero);
//                             if (y.actsExpanded()) {
//                                 for (index_t k=0; k<numEq; k++) {
//                                     const Scalar y_0 = y_p[INDEX2(k, 0, numEq)];
//                                     const Scalar y_1 = y_p[INDEX2(k, 1, numEq)];
//                                     const Scalar y_2 = y_p[INDEX2(k, 2, numEq)];
//                                     const Scalar y_3 = y_p[INDEX2(k, 3, numEq)];
//                                     const Scalar tmp0 = 6.*w7*(y_1 + y_2);
//                                     const Scalar tmp1 = 6.*w7*(y_0 + y_3);
//                                     EM_F[INDEX2(k,2,numEq)] = tmp0 + 6.*w5*y_3 + 6.*w6*y_0;
//                                     EM_F[INDEX2(k,3,numEq)] = tmp1 + 6.*w5*y_2 + 6.*w6*y_1;
//                                     EM_F[INDEX2(k,6,numEq)] = tmp1 + 6.*w5*y_1 + 6.*w6*y_2;
//                                     EM_F[INDEX2(k,7,numEq)] = tmp0 + 6.*w5*y_0 + 6.*w6*y_3;
//                                 }
//                             } else { // constant data
//                                 for (index_t k=0; k<numEq; k++) {
//                                     EM_F[INDEX2(k,2,numEq)] = 36.*w7*y_p[k];
//                                     EM_F[INDEX2(k,3,numEq)] = 36.*w7*y_p[k];
//                                     EM_F[INDEX2(k,6,numEq)] = 36.*w7*y_p[k];
//                                     EM_F[INDEX2(k,7,numEq)] = 36.*w7*y_p[k];
//                                 }
//                             }
//                         }
//                         const index_t firstNode=m_NN[0]*m_NN[1]*k2+m_NN[0]*(m_NN[1]-2)+k0;
//                         domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F,
//                                 add_EM_S, add_EM_F, firstNode, numEq, numComp);
//                     } // k0 loop
//                 } // k2 loop
//             } // colouring
//         } // face 3

//         if (domain->m_faceOffset[4] > -1) {
//             if (add_EM_S)
//                 fill(EM_S.begin(), EM_S.end(), zero);
//             if (add_EM_F)
//                 fill(EM_F.begin(), EM_F.end(), zero);

//             for (index_t k1_0=0; k1_0<2; k1_0++) { // colouring
// #pragma omp for
//                 for (index_t k1=k1_0; k1<NE1; k1+=2) {
//                     for (index_t k0=0; k0<NE0; ++k0) {
//                         const index_t e = domain->m_faceOffset[4]+INDEX2(k0,k1,NE0);
//                         ///////////////
//                         // process d //
//                         ///////////////
//                         if (add_EM_S) {
//                             const Scalar* d_p = d.getSampleDataRO(id, zero);
//                             if (d.actsExpanded()) {
//                                 for (index_t k=0; k<numEq; k++) {
//                                     for (index_t m=0; m<numComp; m++) {
//                                         const Scalar d_0 = d_p[INDEX3(k,m,0,numEq,numComp)];
//                                         const Scalar d_1 = d_p[INDEX3(k,m,1,numEq,numComp)];
//                                         const Scalar d_2 = d_p[INDEX3(k,m,2,numEq,numComp)];
//                                         const Scalar d_3 = d_p[INDEX3(k,m,3,numEq,numComp)];
//                                         const Scalar tmp0 = w10*(d_0 + d_2);
//                                         const Scalar tmp1 = w11*(d_1 + d_3);
//                                         const Scalar tmp2 = w12*(d_0 + d_1 + d_2 + d_3);
//                                         const Scalar tmp3 = w12*(d_1 + d_2);
//                                         const Scalar tmp4 = w10*(d_1 + d_3);
//                                         const Scalar tmp5 = w11*(d_0 + d_2);
//                                         const Scalar tmp6 = w12*(d_0 + d_3);
//                                         const Scalar tmp7 = w10*(d_0 + d_1);
//                                         const Scalar tmp8 = w11*(d_2 + d_3);
//                                         const Scalar tmp9 = w10*(d_2 + d_3);
//                                         const Scalar tmp10 = w11*(d_0 + d_1);
//                                         EM_S[INDEX4(k,m,0,0,numEq,numComp,8)] = d_0*w14 + d_3*w13 + tmp3;
//                                         EM_S[INDEX4(k,m,0,1,numEq,numComp,8)] = tmp10 + tmp9;
//                                         EM_S[INDEX4(k,m,0,2,numEq,numComp,8)] = tmp4 + tmp5;
//                                         EM_S[INDEX4(k,m,0,3,numEq,numComp,8)] = tmp2;
//                                         EM_S[INDEX4(k,m,1,0,numEq,numComp,8)] = tmp10 + tmp9;
//                                         EM_S[INDEX4(k,m,1,1,numEq,numComp,8)] = d_1*w14 + d_2*w13 + tmp6;
//                                         EM_S[INDEX4(k,m,1,2,numEq,numComp,8)] = tmp2;
//                                         EM_S[INDEX4(k,m,1,3,numEq,numComp,8)] = tmp0 + tmp1;
//                                         EM_S[INDEX4(k,m,2,0,numEq,numComp,8)] = tmp4 + tmp5;
//                                         EM_S[INDEX4(k,m,2,1,numEq,numComp,8)] = tmp2;
//                                         EM_S[INDEX4(k,m,2,2,numEq,numComp,8)] = d_1*w13 + d_2*w14 + tmp6;
//                                         EM_S[INDEX4(k,m,2,3,numEq,numComp,8)] = tmp7 + tmp8;
//                                         EM_S[INDEX4(k,m,3,0,numEq,numComp,8)] = tmp2;
//                                         EM_S[INDEX4(k,m,3,1,numEq,numComp,8)] = tmp0 + tmp1;
//                                         EM_S[INDEX4(k,m,3,2,numEq,numComp,8)] = tmp7 + tmp8;
//                                         EM_S[INDEX4(k,m,3,3,numEq,numComp,8)] = d_0*w13 + d_3*w14 + tmp3;
//                                     }
//                                  }
//                             } else { // constant data
//                                 for (index_t k=0; k<numEq; k++) {
//                                     for (index_t m=0; m<numComp; m++) {
//                                         const Scalar wd0 = 4.*d_p[INDEX2(k, m, numEq)]*w12;
//                                         EM_S[INDEX4(k,m,0,0,numEq,numComp,8)] = 4.*wd0;
//                                         EM_S[INDEX4(k,m,0,1,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,0,2,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,0,3,numEq,numComp,8)] =    wd0;
//                                         EM_S[INDEX4(k,m,1,0,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,1,1,numEq,numComp,8)] = 4.*wd0;
//                                         EM_S[INDEX4(k,m,1,2,numEq,numComp,8)] =    wd0;
//                                         EM_S[INDEX4(k,m,1,3,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,2,0,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,2,1,numEq,numComp,8)] =    wd0;
//                                         EM_S[INDEX4(k,m,2,2,numEq,numComp,8)] = 4.*wd0;
//                                         EM_S[INDEX4(k,m,2,3,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,3,0,numEq,numComp,8)] =    wd0;
//                                         EM_S[INDEX4(k,m,3,1,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,3,2,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,3,3,numEq,numComp,8)] = 4.*wd0;
//                                     }
//                                 }
//                             }
//                         }
//                         ///////////////
//                         // process y //
//                         ///////////////
//                         if (add_EM_F) {
//                             const Scalar* y_p = y.getSampleDataRO(id, zero);
//                             if (y.actsExpanded()) {
//                                 for (index_t k=0; k<numEq; k++) {
//                                     const Scalar y_0 = y_p[INDEX2(k, 0, numEq)];
//                                     const Scalar y_1 = y_p[INDEX2(k, 1, numEq)];
//                                     const Scalar y_2 = y_p[INDEX2(k, 2, numEq)];
//                                     const Scalar y_3 = y_p[INDEX2(k, 3, numEq)];
//                                     const Scalar tmp0 = 6.*w12*(y_1 + y_2);
//                                     const Scalar tmp1 = 6.*w12*(y_0 + y_3);
//                                     EM_F[INDEX2(k,0,numEq)] = tmp0 + 6.*w10*y_3 + 6.*w11*y_0;
//                                     EM_F[INDEX2(k,1,numEq)] = tmp1 + 6.*w10*y_2 + 6.*w11*y_1;
//                                     EM_F[INDEX2(k,2,numEq)] = tmp1 + 6.*w10*y_1 + 6.*w11*y_2;
//                                     EM_F[INDEX2(k,3,numEq)] = tmp0 + 6.*w10*y_0 + 6.*w11*y_3;
//                                 }
//                             } else { // constant data
//                                 for (index_t k=0; k<numEq; k++) {
//                                     EM_F[INDEX2(k,0,numEq)] = 36.*w12*y_p[k];
//                                     EM_F[INDEX2(k,1,numEq)] = 36.*w12*y_p[k];
//                                     EM_F[INDEX2(k,2,numEq)] = 36.*w12*y_p[k];
//                                     EM_F[INDEX2(k,3,numEq)] = 36.*w12*y_p[k];
//                                 }
//                             }
//                         }
//                         const index_t firstNode=m_NN[0]*k1+k0;
//                         domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F,
//                                 add_EM_S, add_EM_F, firstNode, numEq, numComp);
//                     } // k0 loop
//                 } // k1 loop
//             } // colouring
//         } // face 4

//         if (domain->m_faceOffset[5] > -1) {
//             if (add_EM_S)
//                 fill(EM_S.begin(), EM_S.end(), zero);
//             if (add_EM_F)
//                 fill(EM_F.begin(), EM_F.end(), zero);

//             for (index_t k1_0=0; k1_0<2; k1_0++) { // colouring
// #pragma omp for
//                 for (index_t k1=k1_0; k1<NE1; k1+=2) {
//                     for (index_t k0=0; k0<NE0; ++k0) {
//                         const index_t e = domain->m_faceOffset[5]+INDEX2(k0,k1,NE0);
//                         ///////////////
//                         // process d //
//                         ///////////////
//                         if (add_EM_S) {
//                             const Scalar* d_p = d.getSampleDataRO(id, zero);
//                             if (d.actsExpanded()) {
//                                 for (index_t k=0; k<numEq; k++) {
//                                     for (index_t m=0; m<numComp; m++) {
//                                         const Scalar d_0 = d_p[INDEX3(k,m,0,numEq,numComp)];
//                                         const Scalar d_1 = d_p[INDEX3(k,m,1,numEq,numComp)];
//                                         const Scalar d_2 = d_p[INDEX3(k,m,2,numEq,numComp)];
//                                         const Scalar d_3 = d_p[INDEX3(k,m,3,numEq,numComp)];
//                                         const Scalar tmp0 = w12*(d_0 + d_1 + d_2 + d_3);
//                                         const Scalar tmp1 = w10*(d_1 + d_3);
//                                         const Scalar tmp2 = w11*(d_0 + d_2);
//                                         const Scalar tmp3 = w10*(d_2 + d_3);
//                                         const Scalar tmp4 = w11*(d_0 + d_1);
//                                         const Scalar tmp5 = w10*(d_0 + d_1);
//                                         const Scalar tmp6 = w11*(d_2 + d_3);
//                                         const Scalar tmp7 = w12*(d_1 + d_2);
//                                         const Scalar tmp8 = w10*(d_0 + d_2);
//                                         const Scalar tmp9 = w11*(d_1 + d_3);
//                                         const Scalar tmp10 = w12*(d_0 + d_3);
//                                         EM_S[INDEX4(k,m,4,4,numEq,numComp,8)] = d_0*w14 + d_3*w13 + tmp7;
//                                         EM_S[INDEX4(k,m,5,4,numEq,numComp,8)] = tmp3 + tmp4;
//                                         EM_S[INDEX4(k,m,6,4,numEq,numComp,8)] = tmp1 + tmp2;
//                                         EM_S[INDEX4(k,m,7,4,numEq,numComp,8)] = tmp0;
//                                         EM_S[INDEX4(k,m,4,5,numEq,numComp,8)] = tmp3 + tmp4;
//                                         EM_S[INDEX4(k,m,5,5,numEq,numComp,8)] = d_1*w14 + d_2*w13 + tmp10;
//                                         EM_S[INDEX4(k,m,6,5,numEq,numComp,8)] = tmp0;
//                                         EM_S[INDEX4(k,m,7,5,numEq,numComp,8)] = tmp8 + tmp9;
//                                         EM_S[INDEX4(k,m,4,6,numEq,numComp,8)] = tmp1 + tmp2;
//                                         EM_S[INDEX4(k,m,5,6,numEq,numComp,8)] = tmp0;
//                                         EM_S[INDEX4(k,m,6,6,numEq,numComp,8)] = d_1*w13 + d_2*w14 + tmp10;
//                                         EM_S[INDEX4(k,m,7,6,numEq,numComp,8)] = tmp5 + tmp6;
//                                         EM_S[INDEX4(k,m,4,7,numEq,numComp,8)] = tmp0;
//                                         EM_S[INDEX4(k,m,5,7,numEq,numComp,8)] = tmp8 + tmp9;
//                                         EM_S[INDEX4(k,m,6,7,numEq,numComp,8)] = tmp5 + tmp6;
//                                         EM_S[INDEX4(k,m,7,7,numEq,numComp,8)] = d_0*w13 + d_3*w14 + tmp7;
//                                     }
//                                  }
//                             } else { // constant data
//                                 for (index_t k=0; k<numEq; k++) {
//                                     for (index_t m=0; m<numComp; m++) {
//                                         const Scalar wd0 = 4.*d_p[INDEX2(k, m, numEq)]*w12;
//                                         EM_S[INDEX4(k,m,4,4,numEq,numComp,8)] = 4.*wd0;
//                                         EM_S[INDEX4(k,m,5,4,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,6,4,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,7,4,numEq,numComp,8)] =    wd0;
//                                         EM_S[INDEX4(k,m,4,5,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,5,5,numEq,numComp,8)] = 4.*wd0;
//                                         EM_S[INDEX4(k,m,6,5,numEq,numComp,8)] =    wd0;
//                                         EM_S[INDEX4(k,m,7,5,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,4,6,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,5,6,numEq,numComp,8)] =    wd0;
//                                         EM_S[INDEX4(k,m,6,6,numEq,numComp,8)] = 4.*wd0;
//                                         EM_S[INDEX4(k,m,7,6,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,4,7,numEq,numComp,8)] =    wd0;
//                                         EM_S[INDEX4(k,m,5,7,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,6,7,numEq,numComp,8)] = 2.*wd0;
//                                         EM_S[INDEX4(k,m,7,7,numEq,numComp,8)] = 4.*wd0;
//                                     }
//                                 }
//                             }
//                         }
//                         ///////////////
//                         // process y //
//                         ///////////////
//                         if (add_EM_F) {
//                             const Scalar* y_p = y.getSampleDataRO(id, zero);
//                             if (y.actsExpanded()) {
//                                 for (index_t k=0; k<numEq; k++) {
//                                     const Scalar y_0 = y_p[INDEX2(k, 0, numEq)];
//                                     const Scalar y_1 = y_p[INDEX2(k, 1, numEq)];
//                                     const Scalar y_2 = y_p[INDEX2(k, 2, numEq)];
//                                     const Scalar y_3 = y_p[INDEX2(k, 3, numEq)];
//                                     const Scalar tmp0 = 6.*w12*(y_1 + y_2);
//                                     const Scalar tmp1 = 6.*w12*(y_0 + y_3);
//                                     EM_F[INDEX2(k,4,numEq)] = tmp0 + 6.*w10*y_3 + 6.*w11*y_0;
//                                     EM_F[INDEX2(k,5,numEq)] = tmp1 + 6.*w10*y_2 + 6.*w11*y_1;
//                                     EM_F[INDEX2(k,6,numEq)] = tmp1 + 6.*w10*y_1 + 6.*w11*y_2;
//                                     EM_F[INDEX2(k,7,numEq)] = tmp0 + 6.*w10*y_0 + 6.*w11*y_3;
//                                 }
//                             } else { // constant data
//                                 for (index_t k=0; k<numEq; k++) {
//                                     EM_F[INDEX2(k,4,numEq)] = 36.*w12*y_p[k];
//                                     EM_F[INDEX2(k,5,numEq)] = 36.*w12*y_p[k];
//                                     EM_F[INDEX2(k,6,numEq)] = 36.*w12*y_p[k];
//                                     EM_F[INDEX2(k,7,numEq)] = 36.*w12*y_p[k];
//                                 }
//                             }
//                         }
//                         const index_t firstNode=m_NN[0]*m_NN[1]*(m_NN[2]-2)+m_NN[0]*k1+k0;
//                         domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F,
//                                 add_EM_S, add_EM_F, firstNode, numEq, numComp);
//                     } // k0 loop
//                 } // k1 loop
//             } // colouring
//         } // face 5
//     } // end of parallel region
}

/****************************************************************************/
// PDE SYSTEM REDUCED
/****************************************************************************/

template<class Scalar>
void DefaultAssembler3D<Scalar>::assemblePDESystemReduced(
                                        AbstractSystemMatrix* mat, Data& rhs,
                                        const Data& A, const Data& B,
                                        const Data& C, const Data& D,
                                        const Data& X, const Data& Y) const
{
    // TODO

//     dim_t numEq, numComp;
//     if (!mat)
//         numEq=numComp=(rhs.isEmpty() ? 1 : rhs.getDataPointSize());
//     else {
//         numEq=mat->getRowBlockSize();
//         numComp=mat->getColumnBlockSize();
//     }

//     const double w0 = m_dx[0]/16;
//     const double w1 = m_dx[1]/16;
//     const double w2 = m_dx[2]/16;
//     const double w3 = m_dx[0]*m_dx[1]/32;
//     const double w4 = m_dx[0]*m_dx[2]/32;
//     const double w5 = m_dx[1]*m_dx[2]/32;
//     const double w6 = m_dx[0]*m_dx[1]/(16*m_dx[2]);
//     const double w7 = m_dx[0]*m_dx[2]/(16*m_dx[1]);
//     const double w8 = m_dx[1]*m_dx[2]/(16*m_dx[0]);
//     const double w9 = m_dx[0]*m_dx[1]*m_dx[2]/64;
//     const dim_t NE0 = m_NE[0];
//     const dim_t NE1 = m_NE[1];
//     const dim_t NE2 = m_NE[2];
//     const bool add_EM_S = (!A.isEmpty() || !B.isEmpty() || !C.isEmpty() || !D.isEmpty());
//     const bool add_EM_F = (!X.isEmpty() || !Y.isEmpty());
//     const Scalar zero = static_cast<Scalar>(0);
//     rhs.requireWrite();

// #pragma omp parallel
//     {
//         vector<Scalar> EM_S(8*8*numEq*numComp, zero);
//         vector<Scalar> EM_F(8*numEq, zero);

//         for (index_t k2_0=0; k2_0<2; k2_0++) { // colouring
// #pragma omp for
//             for (index_t k2=k2_0; k2<NE2; k2+=2) {
//                 for (index_t k1=0; k1<NE1; ++k1) {
//                     for (index_t k0=0; k0<NE0; ++k0)  {
//                         const index_t e = k0 + NE0*k1 + NE0*NE1*k2;
//                         if (add_EM_S)
//                             fill(EM_S.begin(), EM_S.end(), zero);
//                         if (add_EM_F)
//                             fill(EM_F.begin(), EM_F.end(), zero);

//                         ///////////////
//                         // process A //
//                         ///////////////
//                         if (!A.isEmpty()) {
//                             const Scalar* A_p = A.getSampleDataRO(id, zero);
//                             for (index_t k=0; k<numEq; k++) {
//                                 for (index_t m=0; m<numComp; m++) {
//                                     const Scalar Aw00 = A_p[INDEX4(k,0,m,0,numEq,3,numComp)]*w8;
//                                     const Scalar Aw10 = A_p[INDEX4(k,1,m,0,numEq,3,numComp)]*w2;
//                                     const Scalar Aw20 = A_p[INDEX4(k,2,m,0,numEq,3,numComp)]*w1;
//                                     const Scalar Aw01 = A_p[INDEX4(k,0,m,1,numEq,3,numComp)]*w2;
//                                     const Scalar Aw11 = A_p[INDEX4(k,1,m,1,numEq,3,numComp)]*w7;
//                                     const Scalar Aw21 = A_p[INDEX4(k,2,m,1,numEq,3,numComp)]*w0;
//                                     const Scalar Aw02 = A_p[INDEX4(k,0,m,2,numEq,3,numComp)]*w1;
//                                     const Scalar Aw12 = A_p[INDEX4(k,1,m,2,numEq,3,numComp)]*w0;
//                                     const Scalar Aw22 = A_p[INDEX4(k,2,m,2,numEq,3,numComp)]*w6;
//                                     EM_S[INDEX4(k,m,0,0,numEq,numComp,8)]+= Aw00 + Aw01 + Aw02 + Aw10 + Aw11 + Aw12 + Aw20 + Aw21 + Aw22;
//                                     EM_S[INDEX4(k,m,1,0,numEq,numComp,8)]+=-Aw00 - Aw01 - Aw02 + Aw10 + Aw11 + Aw12 + Aw20 + Aw21 + Aw22;
//                                     EM_S[INDEX4(k,m,2,0,numEq,numComp,8)]+= Aw00 + Aw01 + Aw02 - Aw10 - Aw11 - Aw12 + Aw20 + Aw21 + Aw22;
//                                     EM_S[INDEX4(k,m,3,0,numEq,numComp,8)]+=-Aw00 - Aw01 - Aw02 - Aw10 - Aw11 - Aw12 + Aw20 + Aw21 + Aw22;
//                                     EM_S[INDEX4(k,m,4,0,numEq,numComp,8)]+= Aw00 + Aw01 + Aw02 + Aw10 + Aw11 + Aw12 - Aw20 - Aw21 - Aw22;
//                                     EM_S[INDEX4(k,m,5,0,numEq,numComp,8)]+=-Aw00 - Aw01 - Aw02 + Aw10 + Aw11 + Aw12 - Aw20 - Aw21 - Aw22;
//                                     EM_S[INDEX4(k,m,6,0,numEq,numComp,8)]+= Aw00 + Aw01 + Aw02 - Aw10 - Aw11 - Aw12 - Aw20 - Aw21 - Aw22;
//                                     EM_S[INDEX4(k,m,7,0,numEq,numComp,8)]+=-Aw00 - Aw01 - Aw02 - Aw10 - Aw11 - Aw12 - Aw20 - Aw21 - Aw22;
//                                     EM_S[INDEX4(k,m,0,1,numEq,numComp,8)]+=-Aw00 + Aw01 + Aw02 - Aw10 + Aw11 + Aw12 - Aw20 + Aw21 + Aw22;
//                                     EM_S[INDEX4(k,m,1,1,numEq,numComp,8)]+= Aw00 - Aw01 - Aw02 - Aw10 + Aw11 + Aw12 - Aw20 + Aw21 + Aw22;
//                                     EM_S[INDEX4(k,m,2,1,numEq,numComp,8)]+=-Aw00 + Aw01 + Aw02 + Aw10 - Aw11 - Aw12 - Aw20 + Aw21 + Aw22;
//                                     EM_S[INDEX4(k,m,3,1,numEq,numComp,8)]+= Aw00 - Aw01 - Aw02 + Aw10 - Aw11 - Aw12 - Aw20 + Aw21 + Aw22;
//                                     EM_S[INDEX4(k,m,4,1,numEq,numComp,8)]+=-Aw00 + Aw01 + Aw02 - Aw10 + Aw11 + Aw12 + Aw20 - Aw21 - Aw22;
//                                     EM_S[INDEX4(k,m,5,1,numEq,numComp,8)]+= Aw00 - Aw01 - Aw02 - Aw10 + Aw11 + Aw12 + Aw20 - Aw21 - Aw22;
//                                     EM_S[INDEX4(k,m,6,1,numEq,numComp,8)]+=-Aw00 + Aw01 + Aw02 + Aw10 - Aw11 - Aw12 + Aw20 - Aw21 - Aw22;
//                                     EM_S[INDEX4(k,m,7,1,numEq,numComp,8)]+= Aw00 - Aw01 - Aw02 + Aw10 - Aw11 - Aw12 + Aw20 - Aw21 - Aw22;
//                                     EM_S[INDEX4(k,m,0,2,numEq,numComp,8)]+= Aw00 - Aw01 + Aw02 + Aw10 - Aw11 + Aw12 + Aw20 - Aw21 + Aw22;
//                                     EM_S[INDEX4(k,m,1,2,numEq,numComp,8)]+=-Aw00 + Aw01 - Aw02 + Aw10 - Aw11 + Aw12 + Aw20 - Aw21 + Aw22;
//                                     EM_S[INDEX4(k,m,2,2,numEq,numComp,8)]+= Aw00 - Aw01 + Aw02 - Aw10 + Aw11 - Aw12 + Aw20 - Aw21 + Aw22;
//                                     EM_S[INDEX4(k,m,3,2,numEq,numComp,8)]+=-Aw00 + Aw01 - Aw02 - Aw10 + Aw11 - Aw12 + Aw20 - Aw21 + Aw22;
//                                     EM_S[INDEX4(k,m,4,2,numEq,numComp,8)]+= Aw00 - Aw01 + Aw02 + Aw10 - Aw11 + Aw12 - Aw20 + Aw21 - Aw22;
//                                     EM_S[INDEX4(k,m,5,2,numEq,numComp,8)]+=-Aw00 + Aw01 - Aw02 + Aw10 - Aw11 + Aw12 - Aw20 + Aw21 - Aw22;
//                                     EM_S[INDEX4(k,m,6,2,numEq,numComp,8)]+= Aw00 - Aw01 + Aw02 - Aw10 + Aw11 - Aw12 - Aw20 + Aw21 - Aw22;
//                                     EM_S[INDEX4(k,m,7,2,numEq,numComp,8)]+=-Aw00 + Aw01 - Aw02 - Aw10 + Aw11 - Aw12 - Aw20 + Aw21 - Aw22;
//                                     EM_S[INDEX4(k,m,0,3,numEq,numComp,8)]+=-Aw00 - Aw01 + Aw02 - Aw10 - Aw11 + Aw12 - Aw20 - Aw21 + Aw22;
//                                     EM_S[INDEX4(k,m,1,3,numEq,numComp,8)]+= Aw00 + Aw01 - Aw02 - Aw10 - Aw11 + Aw12 - Aw20 - Aw21 + Aw22;
//                                     EM_S[INDEX4(k,m,2,3,numEq,numComp,8)]+=-Aw00 - Aw01 + Aw02 + Aw10 + Aw11 - Aw12 - Aw20 - Aw21 + Aw22;
//                                     EM_S[INDEX4(k,m,3,3,numEq,numComp,8)]+= Aw00 + Aw01 - Aw02 + Aw10 + Aw11 - Aw12 - Aw20 - Aw21 + Aw22;
//                                     EM_S[INDEX4(k,m,4,3,numEq,numComp,8)]+=-Aw00 - Aw01 + Aw02 - Aw10 - Aw11 + Aw12 + Aw20 + Aw21 - Aw22;
//                                     EM_S[INDEX4(k,m,5,3,numEq,numComp,8)]+= Aw00 + Aw01 - Aw02 - Aw10 - Aw11 + Aw12 + Aw20 + Aw21 - Aw22;
//                                     EM_S[INDEX4(k,m,6,3,numEq,numComp,8)]+=-Aw00 - Aw01 + Aw02 + Aw10 + Aw11 - Aw12 + Aw20 + Aw21 - Aw22;
//                                     EM_S[INDEX4(k,m,7,3,numEq,numComp,8)]+= Aw00 + Aw01 - Aw02 + Aw10 + Aw11 - Aw12 + Aw20 + Aw21 - Aw22;
//                                     EM_S[INDEX4(k,m,0,4,numEq,numComp,8)]+= Aw00 + Aw01 - Aw02 + Aw10 + Aw11 - Aw12 + Aw20 + Aw21 - Aw22;
//                                     EM_S[INDEX4(k,m,1,4,numEq,numComp,8)]+=-Aw00 - Aw01 + Aw02 + Aw10 + Aw11 - Aw12 + Aw20 + Aw21 - Aw22;
//                                     EM_S[INDEX4(k,m,2,4,numEq,numComp,8)]+= Aw00 + Aw01 - Aw02 - Aw10 - Aw11 + Aw12 + Aw20 + Aw21 - Aw22;
//                                     EM_S[INDEX4(k,m,3,4,numEq,numComp,8)]+=-Aw00 - Aw01 + Aw02 - Aw10 - Aw11 + Aw12 + Aw20 + Aw21 - Aw22;
//                                     EM_S[INDEX4(k,m,4,4,numEq,numComp,8)]+= Aw00 + Aw01 - Aw02 + Aw10 + Aw11 - Aw12 - Aw20 - Aw21 + Aw22;
//                                     EM_S[INDEX4(k,m,5,4,numEq,numComp,8)]+=-Aw00 - Aw01 + Aw02 + Aw10 + Aw11 - Aw12 - Aw20 - Aw21 + Aw22;
//                                     EM_S[INDEX4(k,m,6,4,numEq,numComp,8)]+= Aw00 + Aw01 - Aw02 - Aw10 - Aw11 + Aw12 - Aw20 - Aw21 + Aw22;
//                                     EM_S[INDEX4(k,m,7,4,numEq,numComp,8)]+=-Aw00 - Aw01 + Aw02 - Aw10 - Aw11 + Aw12 - Aw20 - Aw21 + Aw22;
//                                     EM_S[INDEX4(k,m,0,5,numEq,numComp,8)]+=-Aw00 + Aw01 - Aw02 - Aw10 + Aw11 - Aw12 - Aw20 + Aw21 - Aw22;
//                                     EM_S[INDEX4(k,m,1,5,numEq,numComp,8)]+= Aw00 - Aw01 + Aw02 - Aw10 + Aw11 - Aw12 - Aw20 + Aw21 - Aw22;
//                                     EM_S[INDEX4(k,m,2,5,numEq,numComp,8)]+=-Aw00 + Aw01 - Aw02 + Aw10 - Aw11 + Aw12 - Aw20 + Aw21 - Aw22;
//                                     EM_S[INDEX4(k,m,3,5,numEq,numComp,8)]+= Aw00 - Aw01 + Aw02 + Aw10 - Aw11 + Aw12 - Aw20 + Aw21 - Aw22;
//                                     EM_S[INDEX4(k,m,4,5,numEq,numComp,8)]+=-Aw00 + Aw01 - Aw02 - Aw10 + Aw11 - Aw12 + Aw20 - Aw21 + Aw22;
//                                     EM_S[INDEX4(k,m,5,5,numEq,numComp,8)]+= Aw00 - Aw01 + Aw02 - Aw10 + Aw11 - Aw12 + Aw20 - Aw21 + Aw22;
//                                     EM_S[INDEX4(k,m,6,5,numEq,numComp,8)]+=-Aw00 + Aw01 - Aw02 + Aw10 - Aw11 + Aw12 + Aw20 - Aw21 + Aw22;
//                                     EM_S[INDEX4(k,m,7,5,numEq,numComp,8)]+= Aw00 - Aw01 + Aw02 + Aw10 - Aw11 + Aw12 + Aw20 - Aw21 + Aw22;
//                                     EM_S[INDEX4(k,m,0,6,numEq,numComp,8)]+= Aw00 - Aw01 - Aw02 + Aw10 - Aw11 - Aw12 + Aw20 - Aw21 - Aw22;
//                                     EM_S[INDEX4(k,m,1,6,numEq,numComp,8)]+=-Aw00 + Aw01 + Aw02 + Aw10 - Aw11 - Aw12 + Aw20 - Aw21 - Aw22;
//                                     EM_S[INDEX4(k,m,2,6,numEq,numComp,8)]+= Aw00 - Aw01 - Aw02 - Aw10 + Aw11 + Aw12 + Aw20 - Aw21 - Aw22;
//                                     EM_S[INDEX4(k,m,3,6,numEq,numComp,8)]+=-Aw00 + Aw01 + Aw02 - Aw10 + Aw11 + Aw12 + Aw20 - Aw21 - Aw22;
//                                     EM_S[INDEX4(k,m,4,6,numEq,numComp,8)]+= Aw00 - Aw01 - Aw02 + Aw10 - Aw11 - Aw12 - Aw20 + Aw21 + Aw22;
//                                     EM_S[INDEX4(k,m,5,6,numEq,numComp,8)]+=-Aw00 + Aw01 + Aw02 + Aw10 - Aw11 - Aw12 - Aw20 + Aw21 + Aw22;
//                                     EM_S[INDEX4(k,m,6,6,numEq,numComp,8)]+= Aw00 - Aw01 - Aw02 - Aw10 + Aw11 + Aw12 - Aw20 + Aw21 + Aw22;
//                                     EM_S[INDEX4(k,m,7,6,numEq,numComp,8)]+=-Aw00 + Aw01 + Aw02 - Aw10 + Aw11 + Aw12 - Aw20 + Aw21 + Aw22;
//                                     EM_S[INDEX4(k,m,0,7,numEq,numComp,8)]+=-Aw00 - Aw01 - Aw02 - Aw10 - Aw11 - Aw12 - Aw20 - Aw21 - Aw22;
//                                     EM_S[INDEX4(k,m,1,7,numEq,numComp,8)]+= Aw00 + Aw01 + Aw02 - Aw10 - Aw11 - Aw12 - Aw20 - Aw21 - Aw22;
//                                     EM_S[INDEX4(k,m,2,7,numEq,numComp,8)]+=-Aw00 - Aw01 - Aw02 + Aw10 + Aw11 + Aw12 - Aw20 - Aw21 - Aw22;
//                                     EM_S[INDEX4(k,m,3,7,numEq,numComp,8)]+= Aw00 + Aw01 + Aw02 + Aw10 + Aw11 + Aw12 - Aw20 - Aw21 - Aw22;
//                                     EM_S[INDEX4(k,m,4,7,numEq,numComp,8)]+=-Aw00 - Aw01 - Aw02 - Aw10 - Aw11 - Aw12 + Aw20 + Aw21 + Aw22;
//                                     EM_S[INDEX4(k,m,5,7,numEq,numComp,8)]+= Aw00 + Aw01 + Aw02 - Aw10 - Aw11 - Aw12 + Aw20 + Aw21 + Aw22;
//                                     EM_S[INDEX4(k,m,6,7,numEq,numComp,8)]+=-Aw00 - Aw01 - Aw02 + Aw10 + Aw11 + Aw12 + Aw20 + Aw21 + Aw22;
//                                     EM_S[INDEX4(k,m,7,7,numEq,numComp,8)]+= Aw00 + Aw01 + Aw02 + Aw10 + Aw11 + Aw12 + Aw20 + Aw21 + Aw22;
//                                 }
//                             }
//                         }
//                         ///////////////
//                         // process B //
//                         ///////////////
//                         if (!B.isEmpty()) {
//                             const Scalar* B_p = B.getSampleDataRO(id, zero);
//                             for (index_t k=0; k<numEq; k++) {
//                                 for (index_t m=0; m<numComp; m++) {
//                                     const Scalar wB0 = B_p[INDEX3(k,0,m, numEq, 3)]*w5;
//                                     const Scalar wB1 = B_p[INDEX3(k,1,m, numEq, 3)]*w4;
//                                     const Scalar wB2 = B_p[INDEX3(k,2,m, numEq, 3)]*w3;
//                                     EM_S[INDEX4(k,m,0,0,numEq,numComp,8)]+=-wB0 - wB1 - wB2;
//                                     EM_S[INDEX4(k,m,0,1,numEq,numComp,8)]+=-wB0 - wB1 - wB2;
//                                     EM_S[INDEX4(k,m,0,2,numEq,numComp,8)]+=-wB0 - wB1 - wB2;
//                                     EM_S[INDEX4(k,m,0,3,numEq,numComp,8)]+=-wB0 - wB1 - wB2;
//                                     EM_S[INDEX4(k,m,0,4,numEq,numComp,8)]+=-wB0 - wB1 - wB2;
//                                     EM_S[INDEX4(k,m,0,5,numEq,numComp,8)]+=-wB0 - wB1 - wB2;
//                                     EM_S[INDEX4(k,m,0,6,numEq,numComp,8)]+=-wB0 - wB1 - wB2;
//                                     EM_S[INDEX4(k,m,0,7,numEq,numComp,8)]+=-wB0 - wB1 - wB2;
//                                     EM_S[INDEX4(k,m,1,0,numEq,numComp,8)]+= wB0 - wB1 - wB2;
//                                     EM_S[INDEX4(k,m,1,1,numEq,numComp,8)]+= wB0 - wB1 - wB2;
//                                     EM_S[INDEX4(k,m,1,2,numEq,numComp,8)]+= wB0 - wB1 - wB2;
//                                     EM_S[INDEX4(k,m,1,3,numEq,numComp,8)]+= wB0 - wB1 - wB2;
//                                     EM_S[INDEX4(k,m,1,4,numEq,numComp,8)]+= wB0 - wB1 - wB2;
//                                     EM_S[INDEX4(k,m,1,5,numEq,numComp,8)]+= wB0 - wB1 - wB2;
//                                     EM_S[INDEX4(k,m,1,6,numEq,numComp,8)]+= wB0 - wB1 - wB2;
//                                     EM_S[INDEX4(k,m,1,7,numEq,numComp,8)]+= wB0 - wB1 - wB2;
//                                     EM_S[INDEX4(k,m,2,0,numEq,numComp,8)]+=-wB0 + wB1 - wB2;
//                                     EM_S[INDEX4(k,m,2,1,numEq,numComp,8)]+=-wB0 + wB1 - wB2;
//                                     EM_S[INDEX4(k,m,2,2,numEq,numComp,8)]+=-wB0 + wB1 - wB2;
//                                     EM_S[INDEX4(k,m,2,3,numEq,numComp,8)]+=-wB0 + wB1 - wB2;
//                                     EM_S[INDEX4(k,m,2,4,numEq,numComp,8)]+=-wB0 + wB1 - wB2;
//                                     EM_S[INDEX4(k,m,2,5,numEq,numComp,8)]+=-wB0 + wB1 - wB2;
//                                     EM_S[INDEX4(k,m,2,6,numEq,numComp,8)]+=-wB0 + wB1 - wB2;
//                                     EM_S[INDEX4(k,m,2,7,numEq,numComp,8)]+=-wB0 + wB1 - wB2;
//                                     EM_S[INDEX4(k,m,3,0,numEq,numComp,8)]+= wB0 + wB1 - wB2;
//                                     EM_S[INDEX4(k,m,3,1,numEq,numComp,8)]+= wB0 + wB1 - wB2;
//                                     EM_S[INDEX4(k,m,3,2,numEq,numComp,8)]+= wB0 + wB1 - wB2;
//                                     EM_S[INDEX4(k,m,3,3,numEq,numComp,8)]+= wB0 + wB1 - wB2;
//                                     EM_S[INDEX4(k,m,3,4,numEq,numComp,8)]+= wB0 + wB1 - wB2;
//                                     EM_S[INDEX4(k,m,3,5,numEq,numComp,8)]+= wB0 + wB1 - wB2;
//                                     EM_S[INDEX4(k,m,3,6,numEq,numComp,8)]+= wB0 + wB1 - wB2;
//                                     EM_S[INDEX4(k,m,3,7,numEq,numComp,8)]+= wB0 + wB1 - wB2;
//                                     EM_S[INDEX4(k,m,4,0,numEq,numComp,8)]+=-wB0 - wB1 + wB2;
//                                     EM_S[INDEX4(k,m,4,1,numEq,numComp,8)]+=-wB0 - wB1 + wB2;
//                                     EM_S[INDEX4(k,m,4,2,numEq,numComp,8)]+=-wB0 - wB1 + wB2;
//                                     EM_S[INDEX4(k,m,4,3,numEq,numComp,8)]+=-wB0 - wB1 + wB2;
//                                     EM_S[INDEX4(k,m,4,4,numEq,numComp,8)]+=-wB0 - wB1 + wB2;
//                                     EM_S[INDEX4(k,m,4,5,numEq,numComp,8)]+=-wB0 - wB1 + wB2;
//                                     EM_S[INDEX4(k,m,4,6,numEq,numComp,8)]+=-wB0 - wB1 + wB2;
//                                     EM_S[INDEX4(k,m,4,7,numEq,numComp,8)]+=-wB0 - wB1 + wB2;
//                                     EM_S[INDEX4(k,m,5,0,numEq,numComp,8)]+= wB0 - wB1 + wB2;
//                                     EM_S[INDEX4(k,m,5,1,numEq,numComp,8)]+= wB0 - wB1 + wB2;
//                                     EM_S[INDEX4(k,m,5,2,numEq,numComp,8)]+= wB0 - wB1 + wB2;
//                                     EM_S[INDEX4(k,m,5,3,numEq,numComp,8)]+= wB0 - wB1 + wB2;
//                                     EM_S[INDEX4(k,m,5,4,numEq,numComp,8)]+= wB0 - wB1 + wB2;
//                                     EM_S[INDEX4(k,m,5,5,numEq,numComp,8)]+= wB0 - wB1 + wB2;
//                                     EM_S[INDEX4(k,m,5,6,numEq,numComp,8)]+= wB0 - wB1 + wB2;
//                                     EM_S[INDEX4(k,m,5,7,numEq,numComp,8)]+= wB0 - wB1 + wB2;
//                                     EM_S[INDEX4(k,m,6,0,numEq,numComp,8)]+=-wB0 + wB1 + wB2;
//                                     EM_S[INDEX4(k,m,6,1,numEq,numComp,8)]+=-wB0 + wB1 + wB2;
//                                     EM_S[INDEX4(k,m,6,2,numEq,numComp,8)]+=-wB0 + wB1 + wB2;
//                                     EM_S[INDEX4(k,m,6,3,numEq,numComp,8)]+=-wB0 + wB1 + wB2;
//                                     EM_S[INDEX4(k,m,6,4,numEq,numComp,8)]+=-wB0 + wB1 + wB2;
//                                     EM_S[INDEX4(k,m,6,5,numEq,numComp,8)]+=-wB0 + wB1 + wB2;
//                                     EM_S[INDEX4(k,m,6,6,numEq,numComp,8)]+=-wB0 + wB1 + wB2;
//                                     EM_S[INDEX4(k,m,6,7,numEq,numComp,8)]+=-wB0 + wB1 + wB2;
//                                     EM_S[INDEX4(k,m,7,0,numEq,numComp,8)]+= wB0 + wB1 + wB2;
//                                     EM_S[INDEX4(k,m,7,1,numEq,numComp,8)]+= wB0 + wB1 + wB2;
//                                     EM_S[INDEX4(k,m,7,2,numEq,numComp,8)]+= wB0 + wB1 + wB2;
//                                     EM_S[INDEX4(k,m,7,3,numEq,numComp,8)]+= wB0 + wB1 + wB2;
//                                     EM_S[INDEX4(k,m,7,4,numEq,numComp,8)]+= wB0 + wB1 + wB2;
//                                     EM_S[INDEX4(k,m,7,5,numEq,numComp,8)]+= wB0 + wB1 + wB2;
//                                     EM_S[INDEX4(k,m,7,6,numEq,numComp,8)]+= wB0 + wB1 + wB2;
//                                     EM_S[INDEX4(k,m,7,7,numEq,numComp,8)]+= wB0 + wB1 + wB2;
//                                 }
//                             }
//                         }
//                         ///////////////
//                         // process C //
//                         ///////////////
//                         if (!C.isEmpty()) {
//                             const Scalar* C_p = C.getSampleDataRO(id, zero);
//                             for (index_t k=0; k<numEq; k++) {
//                                 for (index_t m=0; m<numComp; m++) {
//                                     const Scalar wC0 = C_p[INDEX3(k, m, 0, numEq, numComp)]*w5;
//                                     const Scalar wC1 = C_p[INDEX3(k, m, 1, numEq, numComp)]*w4;
//                                     const Scalar wC2 = C_p[INDEX3(k, m, 2, numEq, numComp)]*w3;
//                                     EM_S[INDEX4(k,m,0,0,numEq,numComp,8)]+=-wC0 - wC1 - wC2;
//                                     EM_S[INDEX4(k,m,1,0,numEq,numComp,8)]+=-wC0 - wC1 - wC2;
//                                     EM_S[INDEX4(k,m,2,0,numEq,numComp,8)]+=-wC0 - wC1 - wC2;
//                                     EM_S[INDEX4(k,m,3,0,numEq,numComp,8)]+=-wC0 - wC1 - wC2;
//                                     EM_S[INDEX4(k,m,4,0,numEq,numComp,8)]+=-wC0 - wC1 - wC2;
//                                     EM_S[INDEX4(k,m,5,0,numEq,numComp,8)]+=-wC0 - wC1 - wC2;
//                                     EM_S[INDEX4(k,m,6,0,numEq,numComp,8)]+=-wC0 - wC1 - wC2;
//                                     EM_S[INDEX4(k,m,7,0,numEq,numComp,8)]+=-wC0 - wC1 - wC2;
//                                     EM_S[INDEX4(k,m,0,1,numEq,numComp,8)]+= wC0 - wC1 - wC2;
//                                     EM_S[INDEX4(k,m,1,1,numEq,numComp,8)]+= wC0 - wC1 - wC2;
//                                     EM_S[INDEX4(k,m,2,1,numEq,numComp,8)]+= wC0 - wC1 - wC2;
//                                     EM_S[INDEX4(k,m,3,1,numEq,numComp,8)]+= wC0 - wC1 - wC2;
//                                     EM_S[INDEX4(k,m,4,1,numEq,numComp,8)]+= wC0 - wC1 - wC2;
//                                     EM_S[INDEX4(k,m,5,1,numEq,numComp,8)]+= wC0 - wC1 - wC2;
//                                     EM_S[INDEX4(k,m,6,1,numEq,numComp,8)]+= wC0 - wC1 - wC2;
//                                     EM_S[INDEX4(k,m,7,1,numEq,numComp,8)]+= wC0 - wC1 - wC2;
//                                     EM_S[INDEX4(k,m,0,2,numEq,numComp,8)]+=-wC0 + wC1 - wC2;
//                                     EM_S[INDEX4(k,m,1,2,numEq,numComp,8)]+=-wC0 + wC1 - wC2;
//                                     EM_S[INDEX4(k,m,2,2,numEq,numComp,8)]+=-wC0 + wC1 - wC2;
//                                     EM_S[INDEX4(k,m,3,2,numEq,numComp,8)]+=-wC0 + wC1 - wC2;
//                                     EM_S[INDEX4(k,m,4,2,numEq,numComp,8)]+=-wC0 + wC1 - wC2;
//                                     EM_S[INDEX4(k,m,5,2,numEq,numComp,8)]+=-wC0 + wC1 - wC2;
//                                     EM_S[INDEX4(k,m,6,2,numEq,numComp,8)]+=-wC0 + wC1 - wC2;
//                                     EM_S[INDEX4(k,m,7,2,numEq,numComp,8)]+=-wC0 + wC1 - wC2;
//                                     EM_S[INDEX4(k,m,0,3,numEq,numComp,8)]+= wC0 + wC1 - wC2;
//                                     EM_S[INDEX4(k,m,1,3,numEq,numComp,8)]+= wC0 + wC1 - wC2;
//                                     EM_S[INDEX4(k,m,2,3,numEq,numComp,8)]+= wC0 + wC1 - wC2;
//                                     EM_S[INDEX4(k,m,3,3,numEq,numComp,8)]+= wC0 + wC1 - wC2;
//                                     EM_S[INDEX4(k,m,4,3,numEq,numComp,8)]+= wC0 + wC1 - wC2;
//                                     EM_S[INDEX4(k,m,5,3,numEq,numComp,8)]+= wC0 + wC1 - wC2;
//                                     EM_S[INDEX4(k,m,6,3,numEq,numComp,8)]+= wC0 + wC1 - wC2;
//                                     EM_S[INDEX4(k,m,7,3,numEq,numComp,8)]+= wC0 + wC1 - wC2;
//                                     EM_S[INDEX4(k,m,0,4,numEq,numComp,8)]+=-wC0 - wC1 + wC2;
//                                     EM_S[INDEX4(k,m,1,4,numEq,numComp,8)]+=-wC0 - wC1 + wC2;
//                                     EM_S[INDEX4(k,m,2,4,numEq,numComp,8)]+=-wC0 - wC1 + wC2;
//                                     EM_S[INDEX4(k,m,3,4,numEq,numComp,8)]+=-wC0 - wC1 + wC2;
//                                     EM_S[INDEX4(k,m,4,4,numEq,numComp,8)]+=-wC0 - wC1 + wC2;
//                                     EM_S[INDEX4(k,m,5,4,numEq,numComp,8)]+=-wC0 - wC1 + wC2;
//                                     EM_S[INDEX4(k,m,6,4,numEq,numComp,8)]+=-wC0 - wC1 + wC2;
//                                     EM_S[INDEX4(k,m,7,4,numEq,numComp,8)]+=-wC0 - wC1 + wC2;
//                                     EM_S[INDEX4(k,m,0,5,numEq,numComp,8)]+= wC0 - wC1 + wC2;
//                                     EM_S[INDEX4(k,m,1,5,numEq,numComp,8)]+= wC0 - wC1 + wC2;
//                                     EM_S[INDEX4(k,m,2,5,numEq,numComp,8)]+= wC0 - wC1 + wC2;
//                                     EM_S[INDEX4(k,m,3,5,numEq,numComp,8)]+= wC0 - wC1 + wC2;
//                                     EM_S[INDEX4(k,m,4,5,numEq,numComp,8)]+= wC0 - wC1 + wC2;
//                                     EM_S[INDEX4(k,m,5,5,numEq,numComp,8)]+= wC0 - wC1 + wC2;
//                                     EM_S[INDEX4(k,m,6,5,numEq,numComp,8)]+= wC0 - wC1 + wC2;
//                                     EM_S[INDEX4(k,m,7,5,numEq,numComp,8)]+= wC0 - wC1 + wC2;
//                                     EM_S[INDEX4(k,m,0,6,numEq,numComp,8)]+=-wC0 + wC1 + wC2;
//                                     EM_S[INDEX4(k,m,1,6,numEq,numComp,8)]+=-wC0 + wC1 + wC2;
//                                     EM_S[INDEX4(k,m,2,6,numEq,numComp,8)]+=-wC0 + wC1 + wC2;
//                                     EM_S[INDEX4(k,m,3,6,numEq,numComp,8)]+=-wC0 + wC1 + wC2;
//                                     EM_S[INDEX4(k,m,4,6,numEq,numComp,8)]+=-wC0 + wC1 + wC2;
//                                     EM_S[INDEX4(k,m,5,6,numEq,numComp,8)]+=-wC0 + wC1 + wC2;
//                                     EM_S[INDEX4(k,m,6,6,numEq,numComp,8)]+=-wC0 + wC1 + wC2;
//                                     EM_S[INDEX4(k,m,7,6,numEq,numComp,8)]+=-wC0 + wC1 + wC2;
//                                     EM_S[INDEX4(k,m,0,7,numEq,numComp,8)]+= wC0 + wC1 + wC2;
//                                     EM_S[INDEX4(k,m,1,7,numEq,numComp,8)]+= wC0 + wC1 + wC2;
//                                     EM_S[INDEX4(k,m,2,7,numEq,numComp,8)]+= wC0 + wC1 + wC2;
//                                     EM_S[INDEX4(k,m,3,7,numEq,numComp,8)]+= wC0 + wC1 + wC2;
//                                     EM_S[INDEX4(k,m,4,7,numEq,numComp,8)]+= wC0 + wC1 + wC2;
//                                     EM_S[INDEX4(k,m,5,7,numEq,numComp,8)]+= wC0 + wC1 + wC2;
//                                     EM_S[INDEX4(k,m,6,7,numEq,numComp,8)]+= wC0 + wC1 + wC2;
//                                     EM_S[INDEX4(k,m,7,7,numEq,numComp,8)]+= wC0 + wC1 + wC2;
//                                 }
//                             }
//                         }
//                         ///////////////
//                         // process D //
//                         ///////////////
//                         if (!D.isEmpty()) {
//                             const Scalar* D_p = D.getSampleDataRO(id, zero);
//                             for (index_t k=0; k<numEq; k++) {
//                                 for (index_t m=0; m<numComp; m++) {
//                                     const Scalar wD = D_p[INDEX2(k, m, numEq)]*w9;
//                                     EM_S[INDEX4(k,m,0,0,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,1,0,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,2,0,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,3,0,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,4,0,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,5,0,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,6,0,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,7,0,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,0,1,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,1,1,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,2,1,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,3,1,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,4,1,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,5,1,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,6,1,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,7,1,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,0,2,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,1,2,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,2,2,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,3,2,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,4,2,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,5,2,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,6,2,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,7,2,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,0,3,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,1,3,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,2,3,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,3,3,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,4,3,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,5,3,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,6,3,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,7,3,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,0,4,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,1,4,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,2,4,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,3,4,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,4,4,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,5,4,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,6,4,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,7,4,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,0,5,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,1,5,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,2,5,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,3,5,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,4,5,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,5,5,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,6,5,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,7,5,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,0,6,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,1,6,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,2,6,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,3,6,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,4,6,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,5,6,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,6,6,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,7,6,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,0,7,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,1,7,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,2,7,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,3,7,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,4,7,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,5,7,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,6,7,numEq,numComp,8)]+=wD;
//                                     EM_S[INDEX4(k,m,7,7,numEq,numComp,8)]+=wD;
//                                 }
//                             }
//                         }
//                         ///////////////
//                         // process X //
//                         ///////////////
//                         if (!X.isEmpty()) {
//                             const Scalar* X_p = X.getSampleDataRO(id, zero);
//                             for (index_t k=0; k<numEq; k++) {
//                                 const Scalar wX0 = 8.*X_p[INDEX2(k, 0, numEq)]*w5;
//                                 const Scalar wX1 = 8.*X_p[INDEX2(k, 1, numEq)]*w4;
//                                 const Scalar wX2 = 8.*X_p[INDEX2(k, 2, numEq)]*w3;
//                                 EM_F[INDEX2(k,0,numEq)]+=-wX0 - wX1 - wX2;
//                                 EM_F[INDEX2(k,1,numEq)]+= wX0 - wX1 - wX2;
//                                 EM_F[INDEX2(k,2,numEq)]+=-wX0 + wX1 - wX2;
//                                 EM_F[INDEX2(k,3,numEq)]+= wX0 + wX1 - wX2;
//                                 EM_F[INDEX2(k,4,numEq)]+=-wX0 - wX1 + wX2;
//                                 EM_F[INDEX2(k,5,numEq)]+= wX0 - wX1 + wX2;
//                                 EM_F[INDEX2(k,6,numEq)]+=-wX0 + wX1 + wX2;
//                                 EM_F[INDEX2(k,7,numEq)]+= wX0 + wX1 + wX2;
//                             }
//                         }
//                         ///////////////
//                         // process Y //
//                         ///////////////
//                         if (!Y.isEmpty()) {
//                             const Scalar* Y_p = Y.getSampleDataRO(id, zero);
//                             for (index_t k=0; k<numEq; k++) {
//                                 EM_F[INDEX2(k,0,numEq)]+=8.*Y_p[k]*w9;
//                                 EM_F[INDEX2(k,1,numEq)]+=8.*Y_p[k]*w9;
//                                 EM_F[INDEX2(k,2,numEq)]+=8.*Y_p[k]*w9;
//                                 EM_F[INDEX2(k,3,numEq)]+=8.*Y_p[k]*w9;
//                                 EM_F[INDEX2(k,4,numEq)]+=8.*Y_p[k]*w9;
//                                 EM_F[INDEX2(k,5,numEq)]+=8.*Y_p[k]*w9;
//                                 EM_F[INDEX2(k,6,numEq)]+=8.*Y_p[k]*w9;
//                                 EM_F[INDEX2(k,7,numEq)]+=8.*Y_p[k]*w9;
//                             }
//                         }

//                         // add to matrix (if add_EM_S) and RHS (if add_EM_F)
//                         const index_t firstNode=m_NN[0]*m_NN[1]*k2+m_NN[0]*k1+k0;
//                         domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F,
//                                 add_EM_S, add_EM_F, firstNode, numEq, numComp);
//                     } // end k0 loop
//                 } // end k1 loop
//             } // end k2 loop
//         } // end of colouring
//     } // end of parallel region
}

/****************************************************************************/
// PDE SYSTEM REDUCED BOUNDARY
/****************************************************************************/

template<class Scalar>
void DefaultAssembler3D<Scalar>::assemblePDEBoundarySystemReduced(
                                        AbstractSystemMatrix* mat, Data& rhs,
                                        const Data& d, const Data& y) const
{
    // TODO
//     dim_t numEq, numComp;
//     if (!mat)
//         numEq=numComp=(rhs.isEmpty() ? 1 : rhs.getDataPointSize());
//     else {
//         numEq=mat->getRowBlockSize();
//         numComp=mat->getColumnBlockSize();
//     }
//     const double w0 = m_dx[0]*m_dx[1]/16.;
//     const double w1 = m_dx[0]*m_dx[2]/16.;
//     const double w2 = m_dx[1]*m_dx[2]/16.;
//     const dim_t NE0 = m_NE[0];
//     const dim_t NE1 = m_NE[1];
//     const dim_t NE2 = m_NE[2];
//     const bool add_EM_S = !d.isEmpty();
//     const bool add_EM_F = !y.isEmpty();
//     const Scalar zero = static_cast<Scalar>(0);
//     rhs.requireWrite();

// #pragma omp parallel
//     {
//         vector<Scalar> EM_S(8*8*numEq*numComp, zero);
//         vector<Scalar> EM_F(8*numEq, zero);

//         if (domain->m_faceOffset[0] > -1) {
//             if (add_EM_S)
//                 fill(EM_S.begin(), EM_S.end(), zero);
//             if (add_EM_F)
//                 fill(EM_F.begin(), EM_F.end(), zero);

//             for (index_t k2_0=0; k2_0<2; k2_0++) { // colouring
// #pragma omp for
//                 for (index_t k2=k2_0; k2<NE2; k2+=2) {
//                     for (index_t k1=0; k1<NE1; ++k1) {
//                         const index_t e = INDEX2(k1,k2,NE1);
//                         ///////////////
//                         // process d //
//                         ///////////////
//                         if (add_EM_S) {
//                             const Scalar* d_p = d.getSampleDataRO(id, zero);
//                             for (index_t k=0; k<numEq; k++) {
//                                 for (index_t m=0; m<numComp; m++) {
//                                     const Scalar tmp0 = d_p[INDEX2(k, m, numEq)]*w2;
//                                     EM_S[INDEX4(k,m,0,0,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,2,0,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,4,0,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,6,0,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,0,2,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,2,2,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,4,2,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,6,2,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,0,4,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,2,4,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,4,4,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,6,4,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,0,6,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,2,6,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,4,6,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,6,6,numEq,numComp,8)]=tmp0;
//                                 }
//                             }
//                         }
//                         ///////////////
//                         // process y //
//                         ///////////////
//                         if (add_EM_F) {
//                             const Scalar* y_p = y.getSampleDataRO(id, zero);
//                             for (index_t k=0; k<numEq; k++) {
//                                 EM_F[INDEX2(k,0,numEq)] = 4.*w2*y_p[k];
//                                 EM_F[INDEX2(k,2,numEq)] = 4.*w2*y_p[k];
//                                 EM_F[INDEX2(k,4,numEq)] = 4.*w2*y_p[k];
//                                 EM_F[INDEX2(k,6,numEq)] = 4.*w2*y_p[k];
//                             }
//                         }
//                         const index_t firstNode=m_NN[0]*m_NN[1]*k2+m_NN[0]*k1;
//                         domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F,
//                                 add_EM_S, add_EM_F, firstNode, numEq, numComp);
//                     } // k1 loop
//                 } // k2 loop
//             } // colouring
//         } // face 0

//         if (domain->m_faceOffset[1] > -1) {
//             if (add_EM_S)
//                 fill(EM_S.begin(), EM_S.end(), zero);
//             if (add_EM_F)
//                 fill(EM_F.begin(), EM_F.end(), zero);

//             for (index_t k2_0=0; k2_0<2; k2_0++) { // colouring
// #pragma omp for
//                 for (index_t k2=k2_0; k2<NE2; k2+=2) {
//                     for (index_t k1=0; k1<NE1; ++k1) {
//                         const index_t e = domain->m_faceOffset[1]+INDEX2(k1,k2,NE1);
//                         ///////////////
//                         // process d //
//                         ///////////////
//                         if (add_EM_S) {
//                             const Scalar* d_p = d.getSampleDataRO(id, zero);
//                             for (index_t k=0; k<numEq; k++) {
//                                 for (index_t m=0; m<numComp; m++) {
//                                     const Scalar tmp0 = d_p[INDEX2(k, m, numEq)]*w2;
//                                     EM_S[INDEX4(k,m,1,1,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,3,1,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,5,1,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,7,1,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,1,3,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,3,3,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,5,3,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,7,3,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,1,5,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,3,5,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,5,5,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,7,5,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,1,7,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,3,7,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,5,7,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,7,7,numEq,numComp,8)]=tmp0;
//                                 }
//                             }
//                         }
//                         ///////////////
//                         // process y //
//                         ///////////////
//                         if (add_EM_F) {
//                             const Scalar* y_p = y.getSampleDataRO(id, zero);
//                             for (index_t k=0; k<numEq; k++) {
//                                 EM_F[INDEX2(k,1,numEq)] = 4.*w2*y_p[k];
//                                 EM_F[INDEX2(k,3,numEq)] = 4.*w2*y_p[k];
//                                 EM_F[INDEX2(k,5,numEq)] = 4.*w2*y_p[k];
//                                 EM_F[INDEX2(k,7,numEq)] = 4.*w2*y_p[k];
//                             }
//                         }
//                         const index_t firstNode=m_NN[0]*m_NN[1]*k2+m_NN[0]*(k1+1)-2;
//                         domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F,
//                                 add_EM_S, add_EM_F, firstNode, numEq, numComp);
//                     } // k1 loop
//                 } // k2 loop
//             } // colouring
//         } // face 1

//         if (domain->m_faceOffset[2] > -1) {
//             if (add_EM_S)
//                 fill(EM_S.begin(), EM_S.end(), zero);
//             if (add_EM_F)
//                 fill(EM_F.begin(), EM_F.end(), zero);

//             for (index_t k2_0=0; k2_0<2; k2_0++) { // colouring
// #pragma omp for
//                 for (index_t k2=k2_0; k2<NE2; k2+=2) {
//                     for (index_t k0=0; k0<NE0; ++k0) {
//                         const index_t e = domain->m_faceOffset[2]+INDEX2(k0,k2,NE0);
//                         ///////////////
//                         // process d //
//                         ///////////////
//                         if (add_EM_S) {
//                             const Scalar* d_p = d.getSampleDataRO(id, zero);
//                             for (index_t k=0; k<numEq; k++) {
//                                 for (index_t m=0; m<numComp; m++) {
//                                     const Scalar tmp0 = d_p[INDEX2(k, m, numEq)]*w1;
//                                     EM_S[INDEX4(k,m,0,0,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,1,0,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,4,0,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,5,0,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,0,1,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,1,1,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,4,1,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,5,1,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,0,4,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,1,4,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,4,4,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,5,4,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,0,5,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,1,5,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,4,5,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,5,5,numEq,numComp,8)]=tmp0;
//                                 }
//                             }
//                         }
//                         ///////////////
//                         // process y //
//                         ///////////////
//                         if (add_EM_F) {
//                             const Scalar* y_p = y.getSampleDataRO(id, zero);
//                             for (index_t k=0; k<numEq; k++) {
//                                 EM_F[INDEX2(k,0,numEq)] = 4.*w1*y_p[k];
//                                 EM_F[INDEX2(k,1,numEq)] = 4.*w1*y_p[k];
//                                 EM_F[INDEX2(k,4,numEq)] = 4.*w1*y_p[k];
//                                 EM_F[INDEX2(k,5,numEq)] = 4.*w1*y_p[k];
//                             }
//                         }
//                         const index_t firstNode=m_NN[0]*m_NN[1]*k2+k0;
//                         domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F,
//                                 add_EM_S, add_EM_F, firstNode, numEq, numComp);
//                     } // k0 loop
//                 } // k2 loop
//             } // colouring
//         } // face 2

//         if (domain->m_faceOffset[3] > -1) {
//             if (add_EM_S)
//                 fill(EM_S.begin(), EM_S.end(), zero);
//             if (add_EM_F)
//                 fill(EM_F.begin(), EM_F.end(), zero);

//             for (index_t k2_0=0; k2_0<2; k2_0++) { // colouring
// #pragma omp for
//                 for (index_t k2=k2_0; k2<NE2; k2+=2) {
//                     for (index_t k0=0; k0<NE0; ++k0) {
//                         const index_t e = domain->m_faceOffset[3]+INDEX2(k0,k2,NE0);
//                         ///////////////
//                         // process d //
//                         ///////////////
//                         if (add_EM_S) {
//                             const Scalar* d_p = d.getSampleDataRO(id, zero);
//                             for (index_t k=0; k<numEq; k++) {
//                                 for (index_t m=0; m<numComp; m++) {
//                                     const Scalar tmp0 = d_p[INDEX2(k, m, numEq)]*w1;
//                                     EM_S[INDEX4(k,m,2,2,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,3,2,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,6,2,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,7,2,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,2,3,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,3,3,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,6,3,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,7,3,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,2,6,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,3,6,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,6,6,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,7,6,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,2,7,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,3,7,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,6,7,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,7,7,numEq,numComp,8)]=tmp0;
//                                 }
//                             }
//                         }
//                         ///////////////
//                         // process y //
//                         ///////////////
//                         if (add_EM_F) {
//                             const Scalar* y_p = y.getSampleDataRO(id, zero);
//                             for (index_t k=0; k<numEq; k++) {
//                                 EM_F[INDEX2(k,2,numEq)] = 4.*w1*y_p[k];
//                                 EM_F[INDEX2(k,3,numEq)] = 4.*w1*y_p[k];
//                                 EM_F[INDEX2(k,6,numEq)] = 4.*w1*y_p[k];
//                                 EM_F[INDEX2(k,7,numEq)] = 4.*w1*y_p[k];
//                             }
//                         }
//                         const index_t firstNode=m_NN[0]*m_NN[1]*k2+m_NN[0]*(m_NN[1]-2)+k0;
//                         domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F,
//                                 add_EM_S, add_EM_F, firstNode, numEq, numComp);
//                     } // k0 loop
//                 } // k2 loop
//             } // colouring
//         } // face 3

//         if (domain->m_faceOffset[4] > -1) {
//             if (add_EM_S)
//                 fill(EM_S.begin(), EM_S.end(), zero);
//             if (add_EM_F)
//                 fill(EM_F.begin(), EM_F.end(), zero);

//             for (index_t k1_0=0; k1_0<2; k1_0++) { // colouring
// #pragma omp for
//                 for (index_t k1=k1_0; k1<NE1; k1+=2) {
//                     for (index_t k0=0; k0<NE0; ++k0) {
//                         const index_t e = domain->m_faceOffset[4]+INDEX2(k0,k1,NE0);
//                         ///////////////
//                         // process d //
//                         ///////////////
//                         if (add_EM_S) {
//                             const Scalar* d_p = d.getSampleDataRO(id, zero);
//                             for (index_t k=0; k<numEq; k++) {
//                                 for (index_t m=0; m<numComp; m++) {
//                                     const Scalar tmp0 = d_p[INDEX2(k, m, numEq)]*w0;
//                                     EM_S[INDEX4(k,m,0,0,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,1,0,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,2,0,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,3,0,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,0,1,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,1,1,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,2,1,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,3,1,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,0,2,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,1,2,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,2,2,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,3,2,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,0,3,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,1,3,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,2,3,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,3,3,numEq,numComp,8)]=tmp0;
//                                 }
//                             }
//                         }
//                         ///////////////
//                         // process y //
//                         ///////////////
//                         if (add_EM_F) {
//                             const Scalar* y_p = y.getSampleDataRO(id, zero);
//                             for (index_t k=0; k<numEq; k++) {
//                                 EM_F[INDEX2(k,0,numEq)] = 4.*w0*y_p[k];
//                                 EM_F[INDEX2(k,1,numEq)] = 4.*w0*y_p[k];
//                                 EM_F[INDEX2(k,2,numEq)] = 4.*w0*y_p[k];
//                                 EM_F[INDEX2(k,3,numEq)] = 4.*w0*y_p[k];
//                             }
//                         }
//                         const index_t firstNode=m_NN[0]*k1+k0;
//                         domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F,
//                                 add_EM_S, add_EM_F, firstNode, numEq, numComp);
//                     } // k0 loop
//                 } // k1 loop
//             } // colouring
//         } // face 4

//         if (domain->m_faceOffset[5] > -1) {
//             if (add_EM_S)
//                 fill(EM_S.begin(), EM_S.end(), zero);
//             if (add_EM_F)
//                 fill(EM_F.begin(), EM_F.end(), zero);

//             for (index_t k1_0=0; k1_0<2; k1_0++) { // colouring
// #pragma omp for
//                 for (index_t k1=k1_0; k1<NE1; k1+=2) {
//                     for (index_t k0=0; k0<NE0; ++k0) {
//                         const index_t e = domain->m_faceOffset[5]+INDEX2(k0,k1,NE0);
//                         ///////////////
//                         // process d //
//                         ///////////////
//                         if (add_EM_S) {
//                             const Scalar* d_p = d.getSampleDataRO(id, zero);
//                             for (index_t k=0; k<numEq; k++) {
//                                 for (index_t m=0; m<numComp; m++) {
//                                     const Scalar tmp0 = d_p[INDEX2(k, m, numEq)]*w0;
//                                     EM_S[INDEX4(k,m,4,4,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,5,4,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,6,4,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,7,4,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,4,5,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,5,5,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,6,5,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,7,5,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,4,6,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,5,6,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,6,6,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,7,6,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,4,7,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,5,7,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,6,7,numEq,numComp,8)]=tmp0;
//                                     EM_S[INDEX4(k,m,7,7,numEq,numComp,8)]=tmp0;
//                                 }
//                             }
//                         }
//                         ///////////////
//                         // process y //
//                         ///////////////
//                         if (add_EM_F) {
//                             const Scalar* y_p = y.getSampleDataRO(id, zero);
//                             for (index_t k=0; k<numEq; k++) {
//                                 EM_F[INDEX2(k,4,numEq)] = 4.*w0*y_p[k];
//                                 EM_F[INDEX2(k,5,numEq)] = 4.*w0*y_p[k];
//                                 EM_F[INDEX2(k,6,numEq)] = 4.*w0*y_p[k];
//                                 EM_F[INDEX2(k,7,numEq)] = 4.*w0*y_p[k];
//                             }
//                         }
//                         const index_t firstNode=m_NN[0]*m_NN[1]*(m_NN[2]-2)+m_NN[0]*k1+k0;
//                         domain->addToMatrixAndRHS(mat, rhs, EM_S, EM_F,
//                                 add_EM_S, add_EM_F, firstNode, numEq, numComp);
//                     } // k0 loop
//                 } // k1 loop
//             } // colouring
//         } // face 5
//     } // end of parallel region
}

// instantiate our two supported versions
template class DefaultAssembler3D<escript::DataTypes::real_t>;
template class DefaultAssembler3D<escript::DataTypes::cplx_t>;

} // namespace oxley

