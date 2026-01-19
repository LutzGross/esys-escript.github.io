/*****************************************************************************
*
* Copyright (c) 2003-2019 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
*
*****************************************************************************/

#include <cmath>
#include <iostream>
#include <random>
#include <vector>

#include <escript/Assert.h>

#include <oxley/RefinementAlgorithms.h>
#include <oxley/OtherAlgorithms.h>
#include <oxley/OxleyData.h>
#include <oxley/OxleyException.h>

#include <p4est_bits.h>
#include <p8est_bits.h>

// This file contains various callback functions to decide on refinement.
namespace oxley {

int refine_uniform(p4est_t * p4est, p4est_topidx_t tree, p4est_quadrant_t * quadrant)
{
    p4estData * forestData = (p4estData *) p4est->user_pointer;
    return quadrant->level <= forestData->max_levels_refinement;
}

int refine_uniform(p8est_t * p4est, p4est_topidx_t tree, p8est_quadrant_t * quadrant)
{
    p8estData * octantData = (p8estData *) p4est->user_pointer;
    return quadrant->level <= octantData->max_levels_refinement;
}

int refine_mare2dem(p4est_t * p4est, p4est_topidx_t tree, p4est_quadrant_t * quadrant)
{
    //TODO
//     p4estData * forestData = (p4estData *) p4est->user_pointer;
//     quadrantData * quadData = (quadrantData *) quadrant->p.user_data;
//     std::unordered_map<long,double> * current_solution = forestData->current_solution;
//     std::unordered_map<DoublePair,long,boost::hash<DoublePair>> * NodeIDs = forestData->NodeIDs;

//     // Get the solution value at the current node
//     p4est_qcoord_t xy[2] = {quadrant->x,quadrant->y};
//     long lni = NodeIDs->find(std::make_pair(xy[0],xy[1]))->second;
//     double quad_solution = current_solution->find(lni)->second;

// #ifdef OXLEY_ENABLE_DEBUG
//     std::cout << "refine_mare2dem: " << lni << " (" << xy[0] << "," << xy[1] << ")";
// #endif

//     // Get the Node IDs at the neighbouring nodes
//     double xyz[3] = {0};
//     long neighbour_nodeIDs[4] = {0};

//     double lx = forestData->m_dx[0][quadrant->level];
//     double ly = forestData->m_dx[1][quadrant->level];

//     p4est_qcoord_to_vertex(p4est->connectivity, quadData->treeid, quadrant->x,    quadrant->y+ly, xyz); //N    
//     neighbour_nodeIDs[0] = NodeIDs->find(std::make_pair(xyz[0],xyz[1]))->second;
//     p4est_qcoord_to_vertex(p4est->connectivity, quadData->treeid, quadrant->x,    quadrant->y-ly, xyz); //S
//     neighbour_nodeIDs[1] = NodeIDs->find(std::make_pair(xyz[0],xyz[1]))->second;
//     p4est_qcoord_to_vertex(p4est->connectivity, quadData->treeid, quadrant->x+lx, quadrant->y,    xyz); //E
//     neighbour_nodeIDs[2] = NodeIDs->find(std::make_pair(xyz[0],xyz[1]))->second;
//     p4est_qcoord_to_vertex(p4est->connectivity, quadData->treeid, quadrant->x-lx, quadrant->y,    xyz); //W
//     neighbour_nodeIDs[3] = NodeIDs->find(std::make_pair(xyz[0],xyz[1]))->second;
    
//     // Get the solution at the neighbouring four nodes
//     double average = 0;
//     average += forestData->m_origin[1] >= forestData->m_length[1] ? 0 : current_solution->find(neighbour_nodeIDs[0])->second;
//     average += forestData->m_origin[1] <= forestData->m_length[1] ? 0 : current_solution->find(neighbour_nodeIDs[1])->second;
//     average += forestData->m_origin[0] >= forestData->m_length[0] ? 0 : current_solution->find(neighbour_nodeIDs[2])->second;
//     average += forestData->m_origin[0] <= forestData->m_length[0] ? 0 : current_solution->find(neighbour_nodeIDs[3])->second;
//     average /= 4.0;

// #ifdef OXLEY_ENABLE_DEBUG
//     std::cout << std::endl;
// #endif

//     // Make a decision
//     return (std::abs(average - quad_solution) > MARE2DEM_TOL) && (quadrant->level < forestData->max_levels_refinement);
    return 0;
}

int refine_mare2dem(p8est_t * p8est, p8est_topidx_t tree, p8est_quadrant_t * quadrant)
{
    p8estData * forestData = (p8estData *) p8est->user_pointer;
    quadrantData * quadData = (quadrantData *) quadrant->p.user_data;
    std::unordered_map<long,double> * current_solution = forestData->current_solution;
    std::unordered_map<DoubleTuple,long,boost::hash<DoubleTuple>> * NodeIDs = forestData->NodeIDs;

    // Get the solution value at the current node
    p4est_qcoord_t xyz[3] = {quadrant->x,quadrant->y,quadrant->z};
    long lni = NodeIDs->find(std::make_tuple(xyz[0],xyz[1],xyz[2]))->second;
    double quad_solution = current_solution->find(lni)->second;

#ifdef OXLEY_ENABLE_DEBUG
    std::cout << "refine_mare2dem: " << lni << " (" << xyz[0] << "," << xyz[1] << ")";
#endif

    // Get the Node IDs at the neighbouring nodes
    double xy[3] = {0};
    long neighbour_nodeIDs[8] = {0};

    double lx = forestData->m_dx[0][quadrant->level];
    double ly = forestData->m_dx[1][quadrant->level];
    double lz = forestData->m_dx[2][quadrant->level];

    p8est_qcoord_to_vertex(p8est->connectivity, quadData->treeid, quadrant->x,    quadrant->y+ly, quadrant->z, xy); //N    
    neighbour_nodeIDs[0] = NodeIDs->find(std::make_tuple(xy[0],xy[1],xy[2]))->second;
    p8est_qcoord_to_vertex(p8est->connectivity, quadData->treeid, quadrant->x,    quadrant->y-ly, quadrant->z, xy); //S
    neighbour_nodeIDs[1] = NodeIDs->find(std::make_tuple(xy[0],xy[1],xy[2]))->second;
    p8est_qcoord_to_vertex(p8est->connectivity, quadData->treeid, quadrant->x+lx, quadrant->y,    quadrant->z, xy); //E
    neighbour_nodeIDs[2] = NodeIDs->find(std::make_tuple(xy[0],xy[1],xy[2]))->second;
    p8est_qcoord_to_vertex(p8est->connectivity, quadData->treeid, quadrant->x-lx, quadrant->y,    quadrant->z, xy); //W
    neighbour_nodeIDs[3] = NodeIDs->find(std::make_tuple(xy[0],xy[1],xy[2]))->second;
    p8est_qcoord_to_vertex(p8est->connectivity, quadData->treeid, quadrant->x,    quadrant->y+ly, quadrant->z+lz, xy); //N    
    neighbour_nodeIDs[4] = NodeIDs->find(std::make_tuple(xy[0],xy[1],xy[2]))->second;
    p8est_qcoord_to_vertex(p8est->connectivity, quadData->treeid, quadrant->x,    quadrant->y-ly, quadrant->z+lz, xy); //S
    neighbour_nodeIDs[5] = NodeIDs->find(std::make_tuple(xy[0],xy[1],xy[2]))->second;
    p8est_qcoord_to_vertex(p8est->connectivity, quadData->treeid, quadrant->x+lx, quadrant->y,    quadrant->z+lz, xy); //E
    neighbour_nodeIDs[6] = NodeIDs->find(std::make_tuple(xy[0],xy[1],xy[2]))->second;
    p8est_qcoord_to_vertex(p8est->connectivity, quadData->treeid, quadrant->x-lx, quadrant->y,    quadrant->z+lz, xy); //W
    neighbour_nodeIDs[7] = NodeIDs->find(std::make_tuple(xy[0],xy[1],xy[2]))->second;
    
    // Get the solution at the neighbouring four nodes
    double average = 0;
    average += forestData->m_origin[1] >= forestData->m_length[1] ? 0 : current_solution->find(neighbour_nodeIDs[0])->second;
    average += forestData->m_origin[1] <= forestData->m_length[1] ? 0 : current_solution->find(neighbour_nodeIDs[1])->second;
    average += forestData->m_origin[0] >= forestData->m_length[0] ? 0 : current_solution->find(neighbour_nodeIDs[2])->second;
    average += forestData->m_origin[0] <= forestData->m_length[0] ? 0 : current_solution->find(neighbour_nodeIDs[3])->second;
    //TODO fix below
    average += forestData->m_origin[1] >= forestData->m_length[1] ? 0 : current_solution->find(neighbour_nodeIDs[4])->second;
    average += forestData->m_origin[1] <= forestData->m_length[1] ? 0 : current_solution->find(neighbour_nodeIDs[5])->second;
    average += forestData->m_origin[0] >= forestData->m_length[0] ? 0 : current_solution->find(neighbour_nodeIDs[6])->second;
    average += forestData->m_origin[0] <= forestData->m_length[0] ? 0 : current_solution->find(neighbour_nodeIDs[7])->second;
    average /= 8.0;

#ifdef OXLEY_ENABLE_DEBUG
    std::cout << std::endl;
#endif

    // Make a decision
    return (std::abs(average - quad_solution) > MARE2DEM_TOL) && (quadrant->level < forestData->max_levels_refinement);
    // return 0;
}

// Boundaries
int refine_north(p4est_t * p4est, p4est_topidx_t tree, p4est_quadrant_t * quadrant)
{
    // pointers
    p4estData * forestData = (p4estData *) p4est->user_pointer;
    quadrantData * quadData = (quadrantData *) quadrant->p.user_data;
    double xy[3];
    p4est_qcoord_to_vertex(p4est->connectivity, quadData->treeid, quadrant->x, quadrant->y, xy);

    // refine everything to the right 
    double dx = forestData->refinement_depth;
    double domain_length = forestData->m_length[0];
    double y = domain_length - dx;

    // ne spatial coordinate
    p4est_qcoord_t l = P4EST_QUADRANT_LEN(quadrant->level);
    double xyE[3] = {-1};
    ESYS_ASSERT(quadData->treeid!=-1, "refine_north: invalid treeid");
    p4est_qcoord_to_vertex(p4est->connectivity, quadData->treeid, 
                                    quadrant->x+l, quadrant->y+l, xyE);

    float tol = 1e-8;
    bool do_refinement = (std::abs(xy[1] - y) > tol*l // to the north of the line
                        || (xyE[1] == domain_length)) // or on boundary
                        && (quadrant->level < forestData->max_levels_refinement); // above the limit

    return do_refinement;
}

int refine_south(p4est_t * p4est, p4est_topidx_t tree, p4est_quadrant_t * quadrant)
{
    // pointers
    p4estData * forestData = (p4estData *) p4est->user_pointer;
    quadrantData * quadData = (quadrantData *) quadrant->p.user_data;
    double xy[3];
    p4est_qcoord_to_vertex(p4est->connectivity, quadData->treeid, quadrant->x, quadrant->y, xy);

    // refine everything to the right 
    double dx = forestData->refinement_depth;
    double domain_length = forestData->m_length[0];
    double y = dx;

    // ne spatial coordinate
    p4est_qcoord_t l = P4EST_QUADRANT_LEN(quadrant->level);
    double xyE[3] = {-1};
    ESYS_ASSERT(quadData->treeid!=-1, "refine_south: invalid treeid");
    p4est_qcoord_to_vertex(p4est->connectivity, quadData->treeid, 
                                    quadrant->x+l, quadrant->y+l, xyE);

    float tol = 1e-8;
    bool do_refinement = ( (std::abs(xy[1] - y) < tol*l) // to the south of the line
                        || (xy[1] == 0)) // or on boundary
                        && (quadrant->level < forestData->max_levels_refinement); // above the limit

    return do_refinement;
}

int refine_east(p4est_t * p4est, p4est_topidx_t tree, p4est_quadrant_t * quadrant)
{
    // pointers
    p4estData * forestData = (p4estData *) p4est->user_pointer;
    quadrantData * quadData = (quadrantData *) quadrant->p.user_data;
    double xy[3];
    p4est_qcoord_to_vertex(p4est->connectivity, quadData->treeid, quadrant->x, quadrant->y, xy);

    // refine everything to the right 
    double dx = forestData->refinement_depth;
    double domain_length = forestData->m_length[0];
    double y = domain_length - dx;

    // ne spatial coordinate
    p4est_qcoord_t l = P4EST_QUADRANT_LEN(quadrant->level);
    double xyE[3] = {-1};
    ESYS_ASSERT(quadData->treeid!=-1, "refine_east: invalid treeid");
    p4est_qcoord_to_vertex(p4est->connectivity, quadData->treeid, 
                                    quadrant->x+l, quadrant->y+l, xyE);

    float tol = 1e-8;
    bool do_refinement = ( std::abs(xy[0] - y) > tol*l // to the right of the line
                        || (xyE[0] == domain_length)) // or on boundary
                        && (quadrant->level < forestData->max_levels_refinement); // above the limit

    return do_refinement;
}

int refine_west(p4est_t * p4est, p4est_topidx_t tree, p4est_quadrant_t * quadrant)
{
    // pointers
    p4estData * forestData = (p4estData *) p4est->user_pointer;
    quadrantData * quadData = (quadrantData *) quadrant->p.user_data;
    double xy[3];
    p4est_qcoord_to_vertex(p4est->connectivity, quadData->treeid, quadrant->x, quadrant->y, xy);

    // refine everything to the right 
    double dx = forestData->refinement_depth;
    double domain_length = forestData->m_length[0];
    double y = dx;

    // ne spatial coordinate
    p4est_qcoord_t l = P4EST_QUADRANT_LEN(quadrant->level);
    double xyE[3] = {-1};
    ESYS_ASSERT(quadData->treeid!=-1, "refine_west: invalid treeid");
    p4est_qcoord_to_vertex(p4est->connectivity, quadData->treeid, 
                                    quadrant->x+l, quadrant->y+l, xyE);

    float tol = 1e-8; 
    bool do_refinement = ( std::abs(xy[0] - y) < tol*l // to the west of the line
                        || (xy[0] == 0)) // or on boundary
                        && (quadrant->level < forestData->max_levels_refinement); // above the limit

    return do_refinement;
}

int refine_north(p8est_t * p8est, p8est_topidx_t tree, p8est_quadrant_t * quadrant)
{
    p8estData * forestData = (p8estData *) p8est->user_pointer;
    double dx = forestData->refinement_depth;
    double m_NX = forestData->m_NX[1];
    double domain_length = forestData->m_length[1];
    int steps = dx / m_NX;

    p4est_qcoord_t l = P4EST_QUADRANT_LEN(quadrant->level);
    quadrantData * quadData = (quadrantData *) quadrant->p.user_data;
    double xy[3];
    p8est_qcoord_to_vertex(p8est->connectivity, quadData->treeid, quadrant->x, quadrant->y, quadrant->z, xy);

    float tol = 1e-8;

    return (std::abs(xy[1] - (domain_length - m_NX)) >= tol*l)
        && (quadrant->level < (steps+1)) 
        && (quadrant->level < forestData->max_levels_refinement);
}

int refine_south(p8est_t * p8est, p8est_topidx_t tree, p8est_quadrant_t * quadrant)
{
    p8estData * forestData = (p8estData *) p8est->user_pointer;
    double dx = forestData->refinement_depth;
    double m_NX = forestData->m_NX[1];
    int steps = dx / m_NX;

    quadrantData * quadData = (quadrantData *) quadrant->p.user_data;
    double xy[3];
    p8est_qcoord_to_vertex(p8est->connectivity, quadData->treeid, quadrant->x, quadrant->y, quadrant->z, xy);

    p4est_qcoord_t l = P4EST_QUADRANT_LEN(quadrant->level);
    float tol = 1e-8;

    return  (std::abs(xy[1] - m_NX) <= tol*l)
        && (quadrant->level < (steps+1)) 
        && (quadrant->level < forestData->max_levels_refinement);
}

int refine_east(p8est_t * p8est, p8est_topidx_t tree, p8est_quadrant_t * quadrant)
{
    p8estData * forestData = (p8estData *) p8est->user_pointer;
    double dx = forestData->refinement_depth;
    double m_NX = forestData->m_NX[0];
    double domain_length = forestData->m_length[0];
    int steps = dx / m_NX;

    quadrantData * quadData = (quadrantData *) quadrant->p.user_data;
    double xy[3];
    p8est_qcoord_to_vertex(p8est->connectivity, quadData->treeid, quadrant->x, quadrant->y, quadrant->z, xy);

    float tol = 1e-8;
    p4est_qcoord_t l = P4EST_QUADRANT_LEN(quadrant->level);

    return (std::abs(xy[0] - (domain_length - m_NX)) >= tol*l)
        && (quadrant->level < (steps+1)) 
        && (quadrant->level < forestData->max_levels_refinement);
}

int refine_west(p8est_t * p8est, p8est_topidx_t tree, p8est_quadrant_t * quadrant)
{
    p8estData * forestData = (p8estData *) p8est->user_pointer;
    double dx = forestData->refinement_depth;
    double m_NX = forestData->m_NX[0];
    int steps = dx / m_NX;

    quadrantData * quadData = (quadrantData *) quadrant->p.user_data;
    double xy[3];
    p8est_qcoord_to_vertex(p8est->connectivity, quadData->treeid, quadrant->x, quadrant->y, quadrant->z, xy);

    float tol = 1e-8;
    p4est_qcoord_t l = P4EST_QUADRANT_LEN(quadrant->level);

    return (std::abs(xy[0] - m_NX) <= tol*l)
        && (quadrant->level < (steps+1)) 
        && (quadrant->level < forestData->max_levels_refinement);
}

int refine_top(p8est_t * p8est, p8est_topidx_t tree, p8est_quadrant_t * quadrant)
{
    p8estData * forestData = (p8estData *) p8est->user_pointer;
    double dx = forestData->refinement_depth;
    double m_NX = forestData->m_NX[2];
    double domain_length = forestData->m_length[2];
    int steps = dx / m_NX;

    quadrantData * quadData = (quadrantData *) quadrant->p.user_data;
    double xy[3];
    p8est_qcoord_to_vertex(p8est->connectivity, quadData->treeid, quadrant->x, quadrant->y, quadrant->z, xy);

    float tol = 1e-8;
    p4est_qcoord_t l = P4EST_QUADRANT_LEN(quadrant->level);

    return (std::abs(xy[2] - (domain_length - m_NX)) >= tol*l)
        && (quadrant->level < (steps+1)) 
        && (quadrant->level < forestData->max_levels_refinement);
}

int refine_bottom(p8est_t * p8est, p8est_topidx_t tree, p8est_quadrant_t * quadrant)
{
    p8estData * forestData = (p8estData *) p8est->user_pointer;
    double dx = forestData->refinement_depth;
    double m_NX = forestData->m_NX[2];
    int steps = dx / m_NX;

    quadrantData * quadData = (quadrantData *) quadrant->p.user_data;
    double xy[3];
    p8est_qcoord_to_vertex(p8est->connectivity, quadData->treeid, quadrant->x, quadrant->y, quadrant->z, xy);

    float tol = 1e-8;
    p4est_qcoord_t l = P4EST_QUADRANT_LEN(quadrant->level);

    return (std::abs(xy[2] - m_NX) <= tol*l) 
        && (quadrant->level < (steps+1)) 
        && (quadrant->level < forestData->max_levels_refinement);
}

// Returns true if the point (x,y) is inside the square whose corners are (x0,y0) and (x1,y1)
inline bool point_in_square(double x, double y, double x0, double y0, double x1, double y1)
{
    float tol = 1e-8;
    bool on_boundary = (std::abs(x-x0)>= tol) && (std::abs(x-x1)<= tol) && (std::abs(y-y0)>= tol) && (std::abs(y-y1)<= tol);
    bool in_square = (x > x0) && (x < x1) && (y > y0) && (y < y1);
    return on_boundary || in_square;
}

// Returns true if the point (x,y) is inside the square whose corners are (x0,y0) and (x1,y1)
inline bool point_in_box(double x, double y, double z, double x0, double y0, double z0, double x1, double y1, double z1)
{
    float tol = 1e-8;
    bool on_boundary = (std::abs(x-x0)>=tol) && (std::abs(x-x1)<=tol) && (std::abs(y-y0)>=tol) && (std::abs(y-y1)<=tol) && (std::abs(z-z0)>=tol) && (std::abs(z-z1)<=tol);
    bool in_octant = (x > x0) && (x < x1) && (y > y0) && (y < y1) && (z > z0) && (z < z1);
    return on_boundary || in_octant;
}

// Returns true if the line segment connecting (x1,y1) and (x2,y2) and the line segment connecting
// (x3,y3) and (x4,y4) intersect.
inline bool intersection(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4)
{
    double t = -(x1*(y3-y2)+x2*(y1-y3)+x3*(y2-y1))/(x1*(y4-y3)+x2*(y3-y4)+x4*(y2-y1)+x3*(y1-y2));
    double s =  (x1*(y4-y3)+x3*(y1-y4)+x4*(y3-y1))/(x1*(y4-y3)+x2*(y3-y4)+x4*(y2-y1)+x3*(y1-y2));
    return (s <= 1) && (s >= 0) && (t <= 1) && (t >= 0);
}

// Returns true if the line segment connecting (x1,y1,z2) and (x2,y2,z2) and the line segment connecting
// (x3,y3,z3) and (x4,y4,z4) intersect.
inline bool intersection3D(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4)
{
    double t = -(x1*(y3-y2)+x2*(y1-y3)+x3*(y2-y1))/(x1*(y4-y3)+x2*(y3-y4)+x4*(y2-y1)+x3*(y1-y2));
    double s =  (x1*(y4-y3)+x3*(y1-y4)+x4*(y3-y1))/(x1*(y4-y3)+x2*(y3-y4)+x4*(y2-y1)+x3*(y1-y2));
    return (s<=1) && (s>=0) && (t<=1) && (t>=0);
}

class point_info // for the sake of brevity in refine_region
{
public:
    double x,y;
};

int refine_region(p4est_t * p4est, p4est_topidx_t tree, p4est_quadrant_t * quadrant)
{
    p4estData * forestData = (p4estData *) p4est->user_pointer;
    quadrantData * quadData = (quadrantData *) quadrant->p.user_data;
    double xy[3];
    p4est_qcoord_to_vertex(p4est->connectivity, quadData->treeid, quadrant->x, quadrant->y, xy);
    
    // The corners of the refinement region
    double x[2]={forestData->refinement_boundaries[0],forestData->refinement_boundaries[1]};
    double y[2]={forestData->refinement_boundaries[2],forestData->refinement_boundaries[3]};

    bool x_in_region = (xy[0] >= x[0]) && (xy[0] <= x[1]);
    bool y_in_region = (xy[1] >= y[0]) && (xy[1] <= y[1]);

    bool refinement = x_in_region && y_in_region;

    return refinement && (quadrant->level < forestData->max_levels_refinement);
}

int refine_point(p4est_t * p4est, p4est_topidx_t tree, p4est_quadrant_t * quadrant)
{
    p4estData * forestData = (p4estData *) p4est->user_pointer;
    quadrantData * quadData = (quadrantData *) quadrant->p.user_data;
    double p[2] = {forestData->refinement_boundaries[0], forestData->refinement_boundaries[1]};
    double xy1[3];
    p4est_qcoord_to_vertex(p4est->connectivity, quadData->treeid, quadrant->x, quadrant->y, xy1);

    double xy2[3] = {0};
    p4est_qcoord_t l = P4EST_QUADRANT_LEN(quadrant->level);
    p4est_qcoord_to_vertex(p4est->connectivity, tree, quadrant->x+l, quadrant->y+l, xy2);

    float tol = 1e-8;

    // Check if the point is inside the quadrant
    bool on_boundary =   (std::abs(p[0]-xy1[0]) >= 1e-8) && (std::abs(p[0]-xy2[0]) <= 1e-8)
                      && (std::abs(p[1]-xy1[1]) >= 1e-8) && (std::abs(p[1]-xy2[1]) <= 1e-8);
    bool in_quadrant =   (p[0] >= xy1[0]) && (p[0] <= xy2[0])
                      && (p[1] >= xy1[1]) && (p[1] <= xy2[1]);

    return  (on_boundary || in_quadrant) &&
            (quadrant->level < forestData->max_levels_refinement);
}

int refine_circle(p4est_t * p4est, p4est_topidx_t tree, p4est_quadrant_t * quadrant)
{
    p4estData * forestData = (p4estData *) p4est->user_pointer;
    quadrantData * quadData = (quadrantData *) quadrant->p.user_data;
    double center[2] = {forestData->refinement_boundaries[0], forestData->refinement_boundaries[1]};
    double r = forestData->refinement_boundaries[2];
    double xy1[3];
    p4est_qcoord_to_vertex(p4est->connectivity, quadData->treeid, quadrant->x, quadrant->y, xy1);

    // Upper right point
    double p[3] = {0};
    p4est_qcoord_t l = P4EST_QUADRANT_LEN(quadrant->level);
    p4est_qcoord_to_vertex(p4est->connectivity, tree, quadrant->x+l, quadrant->y+l, p);

    // Center of the quadrant
    p[0] = 0.5*(xy1[0]+p[0]);
    p[1] = 0.5*(xy1[1]+p[1]);

    float tol = 1e-8;

    // Check if the point is inside the circle
    bool do_refinement = std::abs((center[0]-p[0])*(center[0]-p[0]) 
                           + (center[1]-p[1])*(center[1]-p[1]) - r*r) < tol;

    return  do_refinement &&
            (quadrant->level < forestData->max_levels_refinement);
}

int refine_mask(p4est_t * p4est, p4est_topidx_t tree, p4est_quadrant_t * quadrant)
{
    p4estData * forestData = (p4estData *) p4est->user_pointer;
    quadrantData * quadData = (quadrantData *) quadrant->p.user_data;
    escript::Data mask = forestData->mask;
    

    // get the mask value at this point
    long nodeid = quadData->nodeid;
    escript::DataTypes::real_t *dummy(0);
    // const escript::DataTypes::real_t * maskvalue = forestData->mask->getSampleDataRO(nodeid, *dummy);

    // check the value
    // bool do_refinement = (*maskvalue != 0);
    // return do_refinement;
    return false; //TODO
}

int refine_region(p8est_t * p8est, p8est_topidx_t tree, p8est_quadrant_t * quadrant)
{
    p8estData * forestData = (p8estData *) p8est->user_pointer;
    octantData * quadData = (octantData *) quadrant->p.user_data;
    double xy[3];
    p8est_qcoord_to_vertex(p8est->connectivity, quadData->treeid, quadrant->x, quadrant->y, quadrant->z, xy);
    
    // The corners of the refinement region
    double x[2]={forestData->refinement_boundaries[0],forestData->refinement_boundaries[1]};
    double y[2]={forestData->refinement_boundaries[2],forestData->refinement_boundaries[3]};
    double z[2]={forestData->refinement_boundaries[4],forestData->refinement_boundaries[5]};

    // The point being examined
    double l[3] = {forestData->m_dx[0][P4EST_MAXLEVEL-quadrant->level],
                   forestData->m_dx[1][P4EST_MAXLEVEL-quadrant->level],
                   forestData->m_dx[2][P4EST_MAXLEVEL-quadrant->level]};
    float tol = 1e-8;

    bool x_on_boundary = (std::abs(xy[0] - x[0]) >= tol) && (std::abs(xy[0] - x[1]) <= tol);
    bool y_on_boundary = (std::abs(xy[1] - y[0]) >= tol) && (std::abs(xy[1] - y[1]) <= tol);
    bool z_on_boundary = (std::abs(xy[2] - z[0]) >= tol) && (std::abs(xy[2] - z[1]) <= tol);
    bool on_boundary = x_on_boundary && y_on_boundary && z_on_boundary;
    bool x_in_region = (xy[0] >= x[0]) && (xy[0] <= x[1]);
    bool y_in_region = (xy[1] >= y[0]) && (xy[1] <= y[1]);
    bool z_in_region = (xy[2] >= z[0]) && (xy[2] <= z[1]);
    bool in_region = x_in_region && y_in_region && z_in_region;

    
    return (on_boundary || in_region) && (quadrant->level < forestData->max_levels_refinement);
}

int refine_point(p8est_t * p8est, p8est_topidx_t tree, p8est_quadrant_t * quadrant)
{
    p8estData * forestData = (p8estData *) p8est->user_pointer;
    octantData * octData = (octantData *) quadrant->p.user_data;
    double p[3] = {forestData->refinement_boundaries[0], 
                   forestData->refinement_boundaries[1],
                   forestData->refinement_boundaries[2]};

    double xy1[3] = {0};
    p8est_qcoord_to_vertex(p8est->connectivity, tree, quadrant->x, quadrant->y, quadrant->z, xy1);

    double xy2[3] = {0};
    p8est_qcoord_t l = P8EST_QUADRANT_LEN(quadrant->level);
    p8est_qcoord_to_vertex(p8est->connectivity, tree, quadrant->x+l, quadrant->y+l, quadrant->z+l, xy2);

    float tol = 1e-8;
    bool on_boundary =    (std::abs(p[0] - xy1[0]) >= tol) && (std::abs(p[0] - xy2[0]) <= tol)
                       && (std::abs(p[1] - xy1[1]) >= tol) && (std::abs(p[1] - xy2[1]) <= tol)
                       && (std::abs(p[2] - xy1[2]) >= tol) && (std::abs(p[2] - xy2[2]) <= tol);
    bool inside_octant =  (p[0] > xy1[0]) && (p[0] < xy2[0])
                       && (p[1] > xy1[1]) && (p[1] < xy2[1])
                       && (p[2] > xy1[2]) && (p[2] < xy2[2]);
    bool above_threshold = ((int) quadrant->level) < forestData->max_levels_refinement ;

    return  (on_boundary || inside_octant) && above_threshold;
}

int refine_sphere(p8est_t * p8est, p8est_topidx_t tree, p8est_quadrant_t * quadrant)
{
    p8estData * forestData = (p8estData *) p8est->user_pointer;
    quadrantData * quadData = (quadrantData *) quadrant->p.user_data;
    double center[3] = {forestData->refinement_boundaries[0], 
                        forestData->refinement_boundaries[1],
                        forestData->refinement_boundaries[2]};
    double r = forestData->refinement_boundaries[3];
    double xy1[3];
    p8est_qcoord_to_vertex(p8est->connectivity, quadData->treeid, quadrant->x, quadrant->y, quadrant->z, xy1);

    // Upper right point
    double p[3] = {0};
    p8est_qcoord_t l = P4EST_QUADRANT_LEN(quadrant->level);
    p8est_qcoord_to_vertex(p8est->connectivity, tree, quadrant->x+l, quadrant->y+l, quadrant->z+l, p);

    // Center of the quadrant
    p[0] = 0.5*(xy1[0]+p[0]);
    p[1] = 0.5*(xy1[1]+p[1]);
    p[2] = 0.5*(xy1[2]+p[2]);

    float tol = 1e-8;

    // Check if the point is inside the sphere
    bool do_refinement = std::abs((center[0]-p[0])*(center[0]-p[0]) + 
                             (center[1]-p[1])*(center[1]-p[1]) +
                             (center[2]-p[2])*(center[2]-p[2]) - r*r) < tol;

    return  do_refinement &&
            (quadrant->level < forestData->max_levels_refinement);
}

int refine_mask(p8est_t * p8est, p8est_topidx_t tree, p8est_quadrant_t * quadrant)
{
    p8estData * forestData = (p8estData *) p8est->user_pointer;
    quadrantData * quadData = (quadrantData *) quadrant->p.user_data;
    escript::Data mask = forestData->mask;

    // get the mask value at this point
    long nodeid = quadData->nodeid;
    escript::DataTypes::real_t *dummy(0);
    const escript::DataTypes::real_t * maskvalue = mask.getSampleDataRO(nodeid, *dummy);

    // check the value
    bool do_refinement = (*maskvalue != 0);
    return do_refinement;
}


void print_quad_debug_info(p4est_iter_volume_info_t * info, p4est_quadrant_t * quadrant)
{
    double xy[2] = {0.0,0.0};
    if(!p4est_quadrant_is_valid(quadrant))
        std::cout << "WARNING! Invalid quadrant: " << info->treeid << "." << info->quadid << std::endl;
    p4est_qcoord_to_vertex(info->p4est->connectivity, info->treeid, quadrant->x, quadrant->y, &xy[0]);
    std::cout << info->treeid << "." << info->quadid <<
        ": qcoords = " << quadrant->x << ", " << quadrant->y <<
        "\tquadlevel = " << quadrant->level << ", pad8=" << quadrant->pad8 << ", pad16=" << quadrant->pad16 <<
        "\t(x,y) = (" << xy[0] << "," << xy[1] << ")";

    quadrantData *quaddata = (quadrantData *) quadrant->p.user_data;
    std::cout << ",\tnodetag = " << quaddata->nodeTag
        << " locat: " << &quadrant
        << std::endl;
}

void print_quad_debug_info(p8est_iter_volume_info_t * info, p8est_quadrant_t * quadrant)
{
    double xyz[3] = {0.0,0.0,0.0};
    if(!p8est_quadrant_is_valid(quadrant))
        std::cout << "WARNING! Invalid quadrant" << info->treeid << "." << info->quadid << std::endl;
    p8est_qcoord_to_vertex(info->p4est->connectivity, info->treeid, quadrant->x, quadrant->y, quadrant->z, &xyz[0]);
    std::cout << info->treeid << "." << info->quadid << ": (x,y,z) = (" << xyz[0] << "," << xyz[1] << "," << xyz[2] << ")" << std::endl;
}

void gce_first_pass(p4est_iter_volume_info_t * info, void *tmp)
{
    if(!p4est_quadrant_is_valid(info->quad))
        throw OxleyException("Invalid quadrant");
}

void gce_second_pass(p4est_iter_volume_info_t * info, void *tmp)
{
    // Note: For the numbering scheme used by p4est, cf. Burstedde et al. (2011)

    // Get some pointers
    quadrantData * quaddata = (quadrantData *) info->quad->p.user_data;
    addSurfaceData * surfaceinfo = (addSurfaceData *) tmp;

    if(surfaceinfo->oldTag == -1)
    {
        // Spatial indices of x, y, z to pass on
        double x = info->quad->x;
        double y = info->quad->y;
        long n = surfaceinfo->x.size();

        // Work out the length of the quadrant
        double l = P4EST_QUADRANT_LEN(info->quad->level);

        // Check that the point is inside the domain defined by the function
        double xy[2][2];
        p4est_qcoord_to_vertex(info->p4est->connectivity, info->treeid, x, y, xy[0]);
        p4est_qcoord_to_vertex(info->p4est->connectivity, info->treeid, x+l, y, xy[1]);
        bool insideDomain = (xy[0][0] >= surfaceinfo->xmin) && (xy[1][0] <= surfaceinfo->xmax);

        if(insideDomain)
        {
            // This variable records whether a node is above, on (true) or below (false) the curve
            signed ab[4] = {0};
            double increment[4][2] = {{0,0},{l,0},{0,l},{l,l}};
#pragma omp parallel for
            for(int i = 0; i < 4; i++)
                ab[i] = aboveCurve(surfaceinfo->x, surfaceinfo->y,
                                    info->p4est->connectivity, info->treeid, n,
                                    x+increment[i][0], y+increment[i][1]);

            signed sum = 0;
            for(int i = 0; i < 4; i++)
                sum += ab[i];

            if(std::abs(sum) != 4)
                surfaceinfo->oldTag = quaddata->quadTag;
        }
    }
}

void gce_third_pass(p4est_iter_volume_info_t * info, void *tmp)
{
    // Get some pointers
    quadrantData * quaddata = (quadrantData *) info->quad->p.user_data;
    addSurfaceData * surfaceinfo = (addSurfaceData *) tmp;

    double x = info->quad->x;
    double y = info->quad->y;
    long n = surfaceinfo->x.size();

    // Check that the point is inside the domain defined by the function
    double l = P4EST_QUADRANT_LEN(info->quad->level);
    double xy[2][2];
    p4est_qcoord_to_vertex(info->p4est->connectivity, info->treeid, x, y, xy[0]);
    p4est_qcoord_to_vertex(info->p4est->connectivity, info->treeid, x+l, y, xy[1]);
    bool insideDomain = (xy[0][0] >= surfaceinfo->xmin) && (xy[1][0] <= surfaceinfo->xmax);

    bool needToUpdate = (quaddata->quadTag == surfaceinfo->oldTag) || (quaddata->quadTag == -1);

    if(insideDomain && needToUpdate)
    {
        bool above = aboveCurve(surfaceinfo->x, surfaceinfo->y,
            info->p4est->connectivity, info->treeid, n, x, y);
        if(above)
            quaddata->quadTag = surfaceinfo->newTag;
        else
            quaddata->quadTag = surfaceinfo->oldTag;
    }
}

int refine_gce(p4est_t * p4est, p4est_topidx_t tree, p4est_quadrant_t * quad)
{
    // Get some pointers
    addSurfaceData * surfaceinfo = (addSurfaceData *) p4est->user_pointer;
    quadrantData * quaddata = (quadrantData *) quad->p.user_data;

    // If we are in the region being defined
    if(quaddata->quadTag == surfaceinfo->oldTag || quaddata->quadTag == -1)
    {
        // The length of the coordinate vector
        long n = surfaceinfo->x.size();

        // Work out the length of the quadrant
        double l = P4EST_QUADRANT_LEN(quad->level);

        // Spatial indices of x, y, z to pass on
        double x = quad->x;
        double y = quad->y;

        // Check to see if we are inside the domain defined by z[x,y]
        double xy[2][2];
        p4est_qcoord_to_vertex(p4est->connectivity, tree, x, y, xy[0]);
        p4est_qcoord_to_vertex(p4est->connectivity, tree, x+l, y, xy[1]);
        bool insideDomain = (xy[0][0] >= surfaceinfo->xmin) && (xy[1][0] <= surfaceinfo->xmax);
        if(!insideDomain)
            return false;

        // This variable records whether a node is above, on (true) or below (false) the curve
        signed ab[4] = {0};
        double increment[4][2] = {{0,0},{l,0},{0,l},{l,l}};
#pragma omp parallel for
        for(int i = 0; i < 4; i++)
            ab[i] = aboveCurve(surfaceinfo->x, surfaceinfo->y, p4est->connectivity, tree, n, x+increment[i][0], y+increment[i][1]);

        signed sum = 0;
        for(int i = 0; i < 4; i++)
            sum += ab[i];

        return (std::abs(sum) == 4) ? false : true;
    }
    else
    {
        return false;
    }
}

// This is used when debugging to check to see if there are any invalid quads
void gce_first_pass(p8est_iter_volume_info_t * info, void *quad_data)
{
    if(!p8est_quadrant_is_valid(info->quad))
        throw OxleyException("Invalid quadrant");
}

// This works out the tag of the region we are about to subdivide
void gce_second_pass(p8est_iter_volume_info_t * info, void *tmp)
{
    // Note: For the numbering scheme used by p4est, cf. Burstedde et al. (2011)

    // Get some pointers
    octantData * quaddata = (octantData *) info->quad->p.user_data;
    addSurfaceData * surfaceinfo = (addSurfaceData *) tmp;

    if(surfaceinfo->oldTag == -1)
    {
        // Work out the length of the quadrant
        double l = P8EST_QUADRANT_LEN(info->quad->level);

        // Spatial indices of x, y, z to pass on
        double x = info->quad->x;
        double y = info->quad->y;
        double z = info->quad->z;

        long nx = surfaceinfo->x.size();
        long ny = surfaceinfo->y.size();

        // Check that the point is inside the domain defined by the function
        double xy[3][3];
        p8est_qcoord_to_vertex(info->p4est->connectivity, info->treeid, x, y, z, xy[0]);
        p8est_qcoord_to_vertex(info->p4est->connectivity, info->treeid, x+l, y, z, xy[1]);
        p8est_qcoord_to_vertex(info->p4est->connectivity, info->treeid, x, y+l, z, xy[2]);

        bool insideDomain = (xy[0][0] >= surfaceinfo->xmin) && (xy[1][0] <= surfaceinfo->xmax)
                         && (xy[1][1] >= surfaceinfo->ymin) && (xy[2][1] <= surfaceinfo->ymax);

        if(insideDomain)
        {
            // This variable records whether a node is above, on (true) or below (false) the curve
            signed ab[8] = {0};
            double increment[8][3] = {{0,0,0},{l,0,0},{0,l,0},{l,l,0},{0,0,l},{l,0,l},{0,l,l},{l,l,l}};
#pragma omp parallel for
            for(int i = 0; i < 8; i++)
                ab[i] = aboveSurface(surfaceinfo->x, surfaceinfo->y, surfaceinfo->z,
                    info->p4est->connectivity, info->treeid, nx, ny,
                    x+increment[i][0], y+increment[i][1], z+increment[i][2]);

            signed sum = 0;
            for(int i = 0; i < 8; i++)
                sum += ab[i];

            if(std::abs(sum) != 8)
                surfaceinfo->oldTag = quaddata->octantTag;
        }
    }
}

// This tags quads that are in the old region but are not being refined
void gce_third_pass(p8est_iter_volume_info_t * info, void *tmp)
{
    // Get some pointers
    octantData * quadData = (octantData *) info->quad->p.user_data;
    addSurfaceData * surfaceinfo = (addSurfaceData *) tmp;

    double x = info->quad->x;
    double y = info->quad->y;
    double z = info->quad->z;
    long nx = surfaceinfo->x.size();
    long ny = surfaceinfo->y.size();
    double l = P8EST_QUADRANT_LEN(info->quad->level);

    // Work out if we are inside the domain defined by z[x,y]
    double xy[3][3];
    p8est_qcoord_to_vertex(info->p4est->connectivity, info->treeid, x, y, z, xy[0]);
    p8est_qcoord_to_vertex(info->p4est->connectivity, info->treeid, x+l, y, z, xy[1]);
    p8est_qcoord_to_vertex(info->p4est->connectivity, info->treeid, x, y+l, z, xy[2]);
    bool insideDomain = (xy[0][0] >= surfaceinfo->xmin) && (xy[1][0] <= surfaceinfo->xmax)
                    && (xy[0][1] >= surfaceinfo->ymin) && (xy[2][1] <= surfaceinfo->ymax);

    bool needToUpdate = (quadData->octantTag == surfaceinfo->oldTag) || quadData->octantTag == -1;
    if(insideDomain && needToUpdate)
    {
        bool above = aboveSurface(surfaceinfo->x, surfaceinfo->y, surfaceinfo->z,
                        info->p4est->connectivity, info->treeid, nx, ny, x, y, z);
        if(above)
           quadData->octantTag = surfaceinfo->newTag;
        else
           quadData->octantTag = surfaceinfo->oldTag;
    }
}

// This decided whether or not to refine a quadrant
int refine_gce(p8est_t * p8est, p4est_topidx_t tree, p8est_quadrant_t * quad)
{
    // Get the tag info
    addSurfaceData * surfaceinfo = (addSurfaceData *) p8est->user_pointer;
    int oldTag = surfaceinfo->oldTag;
    octantData * quaddata = (octantData *) quad->p.user_data;

    if(quaddata->octantTag == oldTag || quaddata->octantTag == -1)
    {
        // Get the length of the quadrant
        double l = P8EST_QUADRANT_LEN(quad->level);

        // Spatial indices of x, y, z to pass on
        p4est_qcoord_t x = quad->x;
        p4est_qcoord_t y = quad->y;
        p4est_qcoord_t z = quad->z;

        long nx = surfaceinfo->x.size();
        long ny = surfaceinfo->y.size();

        // Check to see if we are inside the domain defined by z[x,y]
        double xy[3][3];
        p8est_qcoord_to_vertex(p8est->connectivity, tree, x, y, z, xy[0]);
        p8est_qcoord_to_vertex(p8est->connectivity, tree, x+l, y, z, xy[1]);
        p8est_qcoord_to_vertex(p8est->connectivity, tree, x, y+l, z, xy[2]);
        bool insideDomain = (xy[0][0] >= surfaceinfo->xmin) && (xy[1][0] <= surfaceinfo->xmax)
                        && (xy[0][1] >= surfaceinfo->ymin) && (xy[2][1] <= surfaceinfo->ymax);
        if(!insideDomain)
            return false;

        signed ab[8] = {0};
        double increment[8][3] = {{0,0,0},{l,0,0},{0,l,0},{l,l,0},{0,0,l},{l,0,l},{0,l,l},{l,l,l}};
#pragma omp parallel for
        for(int i = 0; i < 8; i++)
            ab[i] = aboveSurface(surfaceinfo->x, surfaceinfo->y, surfaceinfo->z,
                p8est->connectivity, tree, nx, ny,
                x+increment[i][0], y+increment[i][1], z+increment[i][2]);

        signed sum = 0;
        for(int i = 0; i < 8; i++)
            sum += ab[i];

        return (std::abs(sum) == 8) ? false : true;
    }
    else
    {
        return false;
    }
}

// void refine_copy_parent_quadrant(p4est_t * p4est, p4est_topidx_t tree,
//                                  int num_outgoing,
//                                  p4est_quadrant_t * outgoing[],
//                                  int num_incoming,
//                                  p4est_quadrant_t * incoming[])
// {
//     if(num_incoming == 4 && num_outgoing == 1)
//     {
//         // parent user data
//         // quadrantData *childData = (quadrantData *) incoming[0]->p.user_data;
//         quadrantData *parentData = (quadrantData *) outgoing[0]->p.user_data;

//         // Averaging
//         parentData->u = 0.0;
//         for(int i = 0; i < 4; i++)
//         {
//             quadrantData *childData = (quadrantData *) incoming[i]->p.user_data;
//             parentData->u+=childData->u;
//         }
//         parentData->u=0.25*parentData->u;

//         // Tags
//         quadrantData *childData = (quadrantData *) incoming[0]->p.user_data;
//         parentData->quadTag=childData->quadTag;

//         // Update the spatial coordinates
//         p4est_qcoord_to_vertex(p4est->connectivity, tree, outgoing[0]->x, outgoing[0]->y, &parentData->xy[0]);
//     }
//     else if(num_incoming == 1 && num_outgoing == 4)
//     {
//         quadrantData *parentData = (quadrantData *) incoming[0]->p.user_data;

//         // Loop over the four children
//         for(int i = 0; i < 4; i++){
//             quadrantData *childData = (quadrantData *) outgoing[i]->p.user_data;
//             childData->u=parentData->u;
//           	childData->quadTag=parentData->quadTag;

//             // Update the spatial coordinates
//             p4est_qcoord_to_vertex(p4est->connectivity, tree, outgoing[i]->x, outgoing[i]->y, &childData->xy[0]);
//         }
//     }
//     else
//     {
//         throw OxleyException("refine_copy_parent_quadrant: Unknown error.");
//     }
// }

void refine_copy_parent_octant(p8est_t * p8est, p4est_topidx_t tree,
                                 int num_outgoing, p8est_quadrant_t * outgoing[],
                                 int num_incoming, p8est_quadrant_t * incoming[])
{
    if(num_incoming == 8 && num_outgoing == 1)
    {
        // parent user data
        octantData *parentData = (octantData *) outgoing[0]->p.user_data;

        // Averaging
        parentData->u = 0.0;
        for(int i = 0; i < 8; i++)
        {
            octantData *childData = (octantData *) incoming[i]->p.user_data;
            parentData->u+=childData->u;
        }
        parentData->u=0.125*parentData->u;

        // Tags
        octantData *childData = (octantData *) incoming[0]->p.user_data;
        parentData->nodeTag=childData->nodeTag;
        parentData->octantTag=childData->octantTag;

        // metadata
        parentData->treeid=childData->treeid;
        parentData->owner=childData->owner;

        // Update the spatial coordinates
        p8est_qcoord_to_vertex(p8est->connectivity, tree, outgoing[0]->x, outgoing[0]->y, outgoing[0]->z, &parentData->xyz[0]);
    }
    else if(num_incoming == 1 && num_outgoing == 8)
    {
        octantData *parentData = (octantData *) incoming[0]->p.user_data;

        // Loop over the eight children
        for(int i = 0; i < 8; i++){
            octantData *childData = (octantData *) outgoing[i]->p.user_data;
            childData->u=parentData->u;
            childData->nodeTag=parentData->nodeTag;
            childData->octantTag=parentData->octantTag;
            childData->treeid=parentData->treeid;
            childData->owner=parentData->owner;

            // Update the spatial coordinates
            p8est_qcoord_to_vertex(p8est->connectivity, tree, outgoing[i]->x, outgoing[i]->y, outgoing[i]->z,&childData->xyz[0]);
        }
    }
    else
    {
        throw OxleyException("refine_copy_parent_quadrant: Unknown error.");
    }
}

void get_interpolateNodesOnElementWorker_data(p4est_iter_volume_info_t * info, void *fxx)
{
    //TODO

    // memcpy(&f_00[0], in.getSampleDataRO(INDEX2(k0,k1, m_NN[0]), sentinel), numComp*sizeof(S));
    // memcpy(&f_01[0], in.getSampleDataRO(INDEX2(k0,k1+1, m_NN[0]), sentinel), numComp*sizeof(S));
    // memcpy(&f_10[0], in.getSampleDataRO(INDEX2(k0+1,k1, m_NN[0]), sentinel), numComp*sizeof(S));
    // memcpy(&f_11[0], in.getSampleDataRO(INDEX2(k0+1,k1+1, m_NN[0]), sentinel), numComp*sizeof(S));

}

void get_interpolateNodesOnElementWorker_data(p8est_iter_volume_info_t * info, void *fxx)
{
    //TODO
}


void get_interpolateNodesOnFacesWorker_data(p4est_iter_volume_info_t * info, void *fxx)
{
    //TODO

//     node->x == 0 || node->x == P4EST_ROOT_LEN ||
//           node->y == 0 || node->y == P4EST_ROOT_LEN 

        //direction 0
//                     memcpy(&f_00[0], in.getSampleDataRO(INDEX2(0,k1, m_NN[0]), sentinel), numComp*sizeof(S));
//                     memcpy(&f_01[0], in.getSampleDataRO(INDEX2(0,k1+1, m_NN[0]), sentinel), numComp*sizeof(S));

    // face 1
    // memcpy(&f_10[0], in.getSampleDataRO(INDEX2(m_NN[0]-1,k1, m_NN[0]), sentinel), numComp*sizeof(S));
//                     memcpy(&f_11[0], in.getSampleDataRO(INDEX2(m_NN[0]-1,k1+1, m_NN[0]), sentinel), numComp*sizeof(S));


    //face 2

// memcpy(&f_00[0], in.getSampleDataRO(INDEX2(k0,0, m_NN[0]), sentinel), numComp*sizeof(S));
//                     memcpy(&f_10[0], in.getSampleDataRO(INDEX2(k0+1,0, m_NN[0]), sentinel), numComp*sizeof(S));


    //face 3
// memcpy(&f_01[0], in.getSampleDataRO(INDEX2(k0,m_NN[1]-1, m_NN[0]), sentinel), numComp*sizeof(S));
//                     memcpy(&f_11[0], in.getSampleDataRO(INDEX2(k0+1,m_NN[1]-1, m_NN[0]), sentinel), numComp*sizeof(S));


    /////////////////////////////////////////////////////////
    // shared
    /////////////////////////////////////////////////////////

    //face 0
    // memcpy(&f_00[0], in.getSampleDataRO(INDEX2(0,k1, m_NN[0]), sentinel), numComp*sizeof(S));
    // memcpy(&f_01[0], in.getSampleDataRO(INDEX2(0,k1+1, m_NN[0]), sentinel), numComp*sizeof(S));

    //face 1
    // memcpy(&f_10[0], in.getSampleDataRO(INDEX2(m_NN[0]-1,k1, m_NN[0]), sentinel), numComp*sizeof(S));
    //memcpy(&f_11[0], in.getSampleDataRO(INDEX2(m_NN[0]-1,k1+1, m_NN[0]), sentinel), numComp*sizeof(S));

    //face 2
// memcpy(&f_00[0], in.getSampleDataRO(INDEX2(k0,0, m_NN[0]), sentinel), numComp*sizeof(S));
//                     memcpy(&f_10[0], in.getSampleDataRO(INDEX2(k0+1,0, m_NN[0]), sentinel), numComp*sizeof(S));


    //face 3
// memcpy(&f_01[0], in.getSampleDataRO(INDEX2(k0,m_NN[1]-1, m_NN[0]), sentinel), numComp*sizeof(S));
//                     memcpy(&f_11[0], in.getSampleDataRO(INDEX2(k0+1,m_NN[1]-1, m_NN[0]), sentinel), numComp*sizeof(S));

}

// void get_interpolateNodesOnFacesWorker_data(p8est_iter_volume_info_t * info, void *fxx);
// {
//     //TODO
// }


void update_node_faceoffset(p4est_iter_volume_info_t * info, void *fxx)
{
    quadrantData * quaddata = (quadrantData *) info->quad->p.user_data;
    quaddata->m_faceOffset[0] = info->quad->x == 0 ? true : false;
    quaddata->m_faceOffset[1] = info->quad->x == P4EST_ROOT_LEN ? true : false;
    quaddata->m_faceOffset[2] = info->quad->y == 0 ? true : false;
    quaddata->m_faceOffset[3] = info->quad->y == P4EST_ROOT_LEN ? true : false;
}

void update_node_faceoffset(p8est_iter_volume_info_t * info, void *fxx)
{
    octantData * octdata = (octantData *) info->quad->p.user_data;
    octdata->m_faceOffset[0] = info->quad->x == 0 ? true : false;
    octdata->m_faceOffset[1] = info->quad->x == P8EST_ROOT_LEN ? true : false;
    octdata->m_faceOffset[2] = info->quad->y == 0 ? true : false;
    octdata->m_faceOffset[3] = info->quad->y == P8EST_ROOT_LEN ? true : false;
    octdata->m_faceOffset[4] = info->quad->z == 0 ? true : false;
    octdata->m_faceOffset[5] = info->quad->z == P8EST_ROOT_LEN ? true : false;
}

void update_RC(p4est_iter_face_info_t *info, void *user_data)
{
    //Get some pointers
    update_RC_data * data = (update_RC_data *) user_data;
    sc_array_t * sides = &(info->sides);

    p4est_iter_face_side_t * side = p4est_iter_fside_array_index_int(sides, 0);
    
    p4est_quadrant_t * quad;
    if(side->is_hanging)
        return;
    quad = side->is.full.quad;

    double xy0[3], xyA[3], xyB[3];
    
    // Do nothing if this isn't a lower quadrant
    p4est_qcoord_to_vertex(data->p4est->connectivity, side->treeid, quad->x, quad->y, xy0);
    quad_info tmp;
    tmp.x=xy0[0];
    tmp.y=xy0[1];
    tmp.level=quad->level;
    bool lower_quadrant=false;
    for(int i=0;i<data->pQuadInfo->size();i++)
    {
        if((tmp.x     == data->pQuadInfo[0][i].x)
        && (tmp.y     == data->pQuadInfo[0][i].y)
        && (tmp.level == data->pQuadInfo[0][i].level))
        {
            lower_quadrant=true;
            break;
        }
    }
    if(!lower_quadrant)
        return;
    
    // Calculate the length of the side
    p4est_qcoord_t l = P4EST_QUADRANT_LEN(quad->level);
    int fn = (int) side->face;
    long lx[4][2] = {{0,0},{l,l},{0,l},{0,l}};
    long ly[4][2] = {{0,l},{0,l},{0,0},{l,l}};

    p4est_qcoord_to_vertex(data->p4est->connectivity, side->treeid, quad->x+lx[fn][0], quad->y+ly[fn][0], xyA);
    long lni0 = data->pNodeIDs->find(std::make_pair(xyA[0],xyA[1]))->second;
    #ifdef OXLEY_ENABLE_DEBUG_UPDATE_RC_EXTRA
        std::cout << "(" << xyA[0] << ", " << xyA[1] << ")\t";
    #endif
    p4est_qcoord_to_vertex(data->p4est->connectivity, side->treeid, quad->x+lx[fn][1], quad->y+ly[fn][1], xyB);
    long lni1 = data->pNodeIDs->find(std::make_pair(xyB[0],xyB[1]))->second;
    #ifdef OXLEY_ENABLE_DEBUG_UPDATE_RC_EXTRA
        std::cout << "--\t(" << xyB[0] << ", " << xyB[1] << ")";
        std::cout << std::endl;
    #endif

    std::vector<long> * idx0 = &data->indices[0][lni0];
    std::vector<long> * idx1 = &data->indices[0][lni1];

    bool dup = false;
    for(int i = 1; i < idx0[0][0] + 1; i++)
        if(idx0[0][i] == lni1)
        {
            dup = true;
            break;
        }

    #ifdef OXLEY_ENABLE_DEBUG_UPDATE_RC_EXTRA
        std::cout << "level= " << (int) quad->level ;
        std::cout << "; hanging side " << (int) side->is_hanging;
        std::cout << "; face= " << fn;
        std::cout << "; (x,y)=(" << xy0[0] << ", " << xy0[1] << ")      \t";
        if(dup)
        {
            std::cout << " connection " << lni0 << "---" << lni1;
            std::cout << "\t(dupliate)" << std::endl;
        }
        else
        {
            std::cout << " connection " << lni0 << "---" << lni1;
            std::cout << "\t(not dupliate)" << std::endl;    
        }        
    #endif

    if(dup == false)
    {
#ifdef DOXLEY_ENABLE_DEBUG
        std::cout << "update_RC " << lni1 << ": (" << xy[0] << ", " << xy[1] << ")" << std::endl; // coordinates
#endif
        idx0[0][0]++;
        idx1[0][0]++;
        ESYS_ASSERT(idx0[0][0]<=4, "update_RC index out of bound");
        ESYS_ASSERT(idx1[0][0]<=4, "update_RC index out of bound");
        idx0[0][idx0[0][0]]=lni1;
        idx1[0][idx1[0][0]]=lni0;
    }
}

// ae tmp
void get_coords(p8est_connectivity_t * connectivity,
                        p4est_topidx_t treeid,
                        p4est_qcoord_t x, p4est_qcoord_t y, p4est_qcoord_t z,
                        double vxyz[3])
{
    const double       *vertices = connectivity->vertices;
    const p4est_topidx_t *vindices;
    // int                 xi, yi, zi;
    double              wx[2], wy[2], wz[2];

    vindices = connectivity->tree_to_vertex + 8 * treeid;

    vxyz[0] = vxyz[1] = vxyz[2] = 0.;

    double divisor = ((p4est_qcoord_t) 1 << P8EST_MAXLEVEL);

    wx[1] = (double) x / divisor;
    wx[0] = 1. - wx[1];

    wy[1] = (double) y / divisor;
    wy[0] = 1. - wy[1];

    wz[1] = (double) z / divisor;
    wz[0] = 1. - wz[1];

    // double w[6] = {x/divisor, (1-x)/divisor, y/divisor, (1-y)/divisor, z/divisor, (1-z)/divisor};

    int ii[8][3] ={{0,0,0},{0,0,1},{0,1,0},{0,1,1},{1,0,0},{1,0,1},{1,1,0},{1,1,1}};

    for(int i = 0; i < 8; i++) {
        double xfactor = wz[ii[i][0]] * wy[ii[i][1]] * wx[ii[i][2]];
        p4est_topidx_t vindex = 3 * (*vindices++);
        vxyz[0] += xfactor * vertices[vindex++];
        vxyz[1] += xfactor * vertices[vindex++];
        vxyz[2] += xfactor * vertices[vindex++];
    }
}

void update_RC(p8est_iter_edge_info *info, void *user_data)
{
    //Get some pointers
    update_RC_data_brick * data = (update_RC_data_brick *) user_data;
    sc_array_t * sides = &(info->sides);

    p8est_iter_edge_side_t * side = p8est_iter_eside_array_index_int(sides, 0); 
    if(side->is_hanging!='\000')
        return;
    p8est_quadrant_t * oct = side->is.full.quad;
    if(oct == nullptr) // oct is allocated to another MPI process by p4est but p4est still calls this function of the octant(?) 
        return;
    double xy0[3], xyA[3], xyB[3];
    p8est_qcoord_to_vertex(data->p8est->connectivity, side->treeid, oct->x, oct->y, oct->z, xy0);

    // Calculate the length of the side
    p8est_qcoord_t l = P8EST_QUADRANT_LEN(oct->level);
    int fn = (int) side->edge;

    //                           0     1     2     3     4     5     6     7     8     9     10    11       
    std::vector<std::vector<p8est_qcoord_t>> lx = {{0,l},{0,l},{0,l},{0,l},{0,0},{l,l},{0,0},{l,l},{0,0},{l,l},{0,0},{l,l}};
    std::vector<std::vector<p8est_qcoord_t>> ly = {{0,0},{l,l},{0,0},{l,l},{0,l},{0,l},{0,l},{0,l},{0,0},{0,0},{l,l},{l,l}};
    std::vector<std::vector<p8est_qcoord_t>> lz = {{0,0},{0,0},{l,l},{l,l},{0,0},{0,0},{l,l},{l,l},{0,l},{0,l},{0,l},{0,l}};

    // Get the neighbouring coordinates
    // p8est_qcoord_to_vertex(data->p8est->connectivity, side->treeid, oct->x+lx[fn][0], oct->y+ly[fn][0], oct->z+lz[fn][0], xyA);
    // long lni0 = data->pNodeIDs->find(std::make_tuple(xyA[0],xyA[1],xyA[2]))->second;
    // p8est_qcoord_to_vertex(data->p8est->connectivity, side->treeid, oct->x+lx[fn][1], oct->y+ly[fn][1], oct->z+lz[fn][1], xyB);
    // long lni1 = data->pNodeIDs->find(std::make_tuple(xyB[0],xyB[1],xyB[2]))->second;

    get_coords(data->p8est->connectivity, side->treeid, oct->x+lx[fn][0], oct->y+ly[fn][0], oct->z+lz[fn][0], xyA);
    long lni0 = data->pNodeIDs->find(std::make_tuple(xyA[0],xyA[1],xyA[2]))->second;
    get_coords(data->p8est->connectivity, side->treeid, oct->x+lx[fn][1], oct->y+ly[fn][1], oct->z+lz[fn][1], xyB);
    long lni1 = data->pNodeIDs->find(std::make_tuple(xyB[0],xyB[1],xyB[2]))->second;

    // ae tmp
    // std::cout << "lni01 = " << lni0 << ", " << lni1 << std::endl;







    IndexVector * idx0 = &data->indices[0][lni0];
    IndexVector * idx1 = &data->indices[0][lni1];

    // Check for duplicates
    bool dup = false;
    for(int i = 1; i < idx0[0][0] + 1; i++)
        if(idx0[0][i] == lni1)
        {
            dup = true;
            break;
        }

    // If this is a new coordinate, add it to the index
    if(dup == false)
    {
        idx0[0][0]++;
        idx1[0][0]++;
        ESYS_ASSERT(idx0[0][0]<7, "update_RC index out of bound");
        ESYS_ASSERT(idx1[0][0]<7, "update_RC index out of bound");
        idx0[0][idx0[0][0]]=lni1;
        idx1[0][idx1[0][0]]=lni0;
    }
}


void update_connections(p4est_iter_volume_info_t *info, void *user_data)
{
    //Get some pointers
    getConnections_data * data = (getConnections_data *) user_data;
    p4est_quadrant_t * quad = info->quad;
    
    // Coordinates
    p4est_qcoord_t length = P4EST_QUADRANT_LEN(quad->level);
    double xy[3];
    long lx[4] = {0,length,0,length};
    long ly[4] = {0,0,length,length};
    long lni[4] = {-1};
#pragma omp parallel for
    for(int i = 0; i < 4; i++)
    {
        p4est_qcoord_to_vertex(data->p4est->connectivity, info->treeid, quad->x+lx[i], quad->y+ly[i], xy);
        xy[0]+=data->m_origin[0];
        xy[1]+=data->m_origin[1];
        lni[i] = data->pNodeIDs->find(std::make_pair(xy[0],xy[1]))->second;
    }


#pragma omp parallel for
    for(int i = 0; i < 4; i++)
    {
        std::vector<escript::DataTypes::index_t> * temp = &data->indices[0][i];
        for(int j = 0; j < 4; j++)
        {
            bool dup = false;
            for(int k = 0; k < data->indices[0][i].size(); k++)
                if(temp[0][k] == lni[j])
                {
                    dup = true;
                    break;
                }
            if(dup == false)
                temp->push_back(lni[j]);
        }
    }

    // std::cout << "xy = " << xy[0] << ", " << xy[1] << std::endl; // coordinates
}

int refine_nodesToNodesFiner(p4est_t * p4est, p4est_topidx_t tree, p4est_quadrant_t * quadrant)
{
    quadrantData * quadData = (quadrantData *) quadrant->p.user_data;
    return quadData->needs_refinement;
}

void refine_copy_parent_quadrant_data(p4est_t * p4est, p4est_topidx_t tree,
                                 int num_outgoing,
                                 p4est_quadrant_t * outgoing[],
                                 int num_incoming,
                                 p4est_quadrant_t * incoming[])
{
    if(num_incoming == 4 && num_outgoing == 1)
    {
        quadrantData *childData = (quadrantData *) incoming[0]->p.user_data;
        if(childData->u_real==nullptr)
        {
            cplx_t new_value;
            for(int n = 0; n < 4; n++)
            {
                quadrantData *childData = (quadrantData *) incoming[n]->p.user_data;
                new_value+=*(childData->u_cplx);
            }

            quadrantData *parentData = (quadrantData *) outgoing[0]->p.user_data;
            new_value*=0.25;
            parentData->u_cplx=&new_value;    
        }
        else
        {
            real_t new_value;
            for(int n = 0; n < 4; n++)
            {
                quadrantData *childData = (quadrantData *) incoming[n]->p.user_data;
                new_value+=*(childData->u_real);
            }

            quadrantData *parentData = (quadrantData *) outgoing[0]->p.user_data;
            new_value*=0.25;
            parentData->u_real=&new_value;  
        }    
    }
    else if(num_incoming == 1 && num_outgoing == 4)
    {
        quadrantData *parentData = (quadrantData *) outgoing[0]->p.user_data;
        if(parentData->u_real==nullptr)
        {
            for(int n = 0; n < 4; n++)
            {
                quadrantData *childData = (quadrantData *) incoming[n]->p.user_data;
                cplx_t new_value = *(parentData->u_cplx);
                childData->u_cplx=&new_value;
            }
        }
        else
        {
            for(int n = 0; n < 4; n++)
            {
                quadrantData *childData = (quadrantData *) incoming[n]->p.user_data;
                real_t new_value = *(parentData->u_real);
                childData->u_real=&new_value;
            }
        }        
    }
    else
    {
        throw OxleyException("refine_copy_parent_quadrant_data: Unknown error.");
    }
}

void refine_copy_parent_element_data(p4est_t * p4est, p4est_topidx_t tree,
                                 int num_outgoing,
                                 p4est_quadrant_t * outgoing[],
                                 int num_incoming,
                                 p4est_quadrant_t * incoming[])
{
    if(num_incoming == 4 && num_outgoing == 1)
    {
        // cplx_t new_value;
        // for(int n = 0; n < 4; n++)
        // {
        //     quadrantData *childData = (quadrantData *) incoming[n]->p.user_data;
        //     new_value+=*(childData->u);
        // }

        // quadrantData *parentData = (quadrantData *) outgoing[0]->p.user_data;
        // new_value*=0.25;
        // parentData->u=&new_value;        
    }
    else if(num_incoming == 1 && num_outgoing == 4)
    {
        // quadrantData *parentData = (quadrantData *) outgoing[0]->p.user_data;
        // for(int n = 0; n < 4; n++)
        // {
        //     quadrantData *childData = (quadrantData *) incoming[n]->p.user_data;
        //     cplx_t new_value = *(parentData->u);
        //     childData->u=&new_value;
        // }
    }
    else
    {
        throw OxleyException("refine_copy_parent_element_data: Unknown error.");
    }
}


} // namespace oxley
