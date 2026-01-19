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

#include <escript/Data.h>
#include <escript/DataTypes.h>

#include <oxley/Brick.h>
#include <oxley/Oxley.h>
#include <oxley/OxleyDomain.h>
#include <oxley/Rectangle.h>

#include <unordered_map>
#include <utility>

#include <boost/functional/hash.hpp>
#include <boost/python/numpy.hpp>

#include "p4est/p4est_iterate.h"
#include "p4est/p8est_iterate.h"

#ifndef __OXLEY_DATA_H__
#define __OXLEY_DATA_H__

// Macroes for array indexing
#define INDEX2(_X1_,_X2_,_N1_) ((_X1_)+(_N1_)*(_X2_))
#define INDEX3(_X1_,_X2_,_X3_,_N1_,_N2_) ((_X1_)+(_N1_)*INDEX2(_X2_,_X3_,_N2_))
#define INDEX4(_X1_,_X2_,_X3_,_X4_,_N1_,_N2_,_N3_) ((_X1_)+(_N1_)*INDEX3(_X2_,_X3_,_X4_,_N2_,_N3_))
#define INDEX5(_X1_,_X2_,_X3_,_X4_,_X5_,_N1_,_N2_,_N3_,_N4_) ((_X1_)+(_N1_)*INDEX4(_X2_,_X3_,_X4_,_X5_,_N2_,_N3_,_N4_))
#define INDEX6(_X1_,_X2_,_X3_,_X4_,_X5_,_X6_,_N1_,_N2_,_N3_,_N4_,_N5_) ((_X1_)+(_N1_)*INDEX5(_X2_,_X3_,_X4_,_X5_,_X6_,_N2_,_N3_,_N4_,_N5_))

////////////////////////////////////////////////////////////////////////
// This file contains the data structures used by Rectangle and Brick
////////////////////////////////////////////////////////////////////////

// Forward declarations
struct addSurfaceData;

//This structure describes the information that is stored at each
//quadrant / octant in the p4est / p8est
struct quadrantData
{
	// The quadrant's tag
	long quadTag = 0;

	// Node tag
	double nodeTag = 0;

	// Spatial coordinates of the corner node that defines the quadrant
	double xy[2] = {0.0,0.0};

	// treeid index
	long treeid = -1;

	// Number of the MPI process that owns this quadrant
	int owner = -1;

	// faceOffset[i]=-1 if face i is not an external face, otherwise it is
    // the index of that face (where i: 0=left, 1=right, 2=bottom, 3=top)
    // escript::DataTypes::IndexVector m_faceOffset;
    bool m_faceOffset[4] = {false};

    // bottom left node id
    long nodeid = -1;

    // variables used during interpolation
    bool needs_refinement=false;
    long ids[4];
    const escript::DataTypes::real_t * u_real;
    const escript::DataTypes::cplx_t * u_cplx;
};

struct borderNodeInfo
{
	int nodeid=-1;
	int neighbours[8]={0};
	int8_t level;
	p4est_qcoord_t x;
	p4est_qcoord_t y;
	p4est_qcoord_t z; // Not used by Rectangle
	// p4est_quadrant_t * quad;
	p4est_topidx_t treeid=-1;
};

struct hangingNodeInfo // used by Rectangle
{
	p4est_qcoord_t x;
	p4est_qcoord_t y;
	int8_t level=-1;
	p4est_topidx_t treeid=-1;
	int8_t face_type=-1; 
	
	p4est_qcoord_t neighbour_x;
	p4est_qcoord_t neighbour_y;
	p4est_qcoord_t neighbour_l;
	p4est_topidx_t neighbour_tree;

	signed int position=-1; // position within the parent quadrant
	p4est_quadrant_t parent; // parent quadrant
	p4est_topidx_t parentTreeid;
};

struct hangingFaceInfo // used for hanging faces in Brick
{
	p4est_qcoord_t x;
	p4est_qcoord_t y;
	p4est_qcoord_t z;
	int8_t level;
	p4est_topidx_t treeid;
	int8_t face_type={-1}; 
	
	p4est_qcoord_t neighbour_x;
	p4est_qcoord_t neighbour_y;
	p4est_qcoord_t neighbour_z;
	p4est_qcoord_t neighbour_level;
	p4est_topidx_t neighbour_tree;
};

struct hangingEdgeInfo // used for hanging edges in Brick
{
	long nodeid=-1;

	p4est_qcoord_t x;
	p4est_qcoord_t y;
	p4est_qcoord_t z; 
	int8_t level;
	p4est_topidx_t treeid;
	int8_t edge_type={-1}; 
	
	p4est_qcoord_t neighbour_x;
	p4est_qcoord_t neighbour_y;
	p4est_qcoord_t neighbour_z;  
	p4est_qcoord_t neighbour_level;
	p4est_topidx_t neighbour_tree;
};

struct octantData
{
	double u = 0.0;

	// The octant's tag
	long octantTag = 0;

	// Node tag
	double nodeTag = 0;

	// Spatial coordinates of the corner node that defines the octant
	double xyz[3] = {0.0,0.0,0.0};

	// treeid index
	long treeid = -1;

	// Number of the MPI process that owns this quadrant
	int owner = -1;

	// faceOffset[i]=-1 if face i is not an external face, otherwise it is
    // the index of that face (where i: 0=left, 1=right, 2=bottom, 3=top)
    // escript::DataTypes::IndexVector m_faceOffset;
    bool m_faceOffset[6] = {false};

    long nodeid = -1;
};

//This structure describes the information that is stored with the p4est
class p4estData
{
public:
	// origin of domain
    double m_origin[2] = {0.0};

    // extent of the domain
    double m_lxy[2] = {0.0};

    // side lengths of domain
    double m_length[2] = {0.0};

    // grid spacings / cell sizes of domain for each level of refinement
    double m_dx[2][P4EST_MAXLEVEL+1] = {{0}};

    // initial grid spacing
    double m_NX[2] = {0};

    // number of face elements per edge (left, right, bottom, top)
    escript::DataTypes::dim_t m_faceCount[4] = {-1};

    // vector that maps each node to a DOF index (used for the coupler)
    escript::DataTypes::IndexVector m_dofMap;

    // periodic boundary conditions
    bool periodic[2] {false, false};

	// maximum levels of recursion to use during refinement
	int max_levels_refinement = 0;
	double refinement_depth=0.0;
	double refinement_boundaries[4]={0.0};
	escript::Data mask; // a mask

	// Pointer to the current solution and Node ID info
	// std::unordered_map<long,double> * current_solution;
	std::unordered_map<DoublePair,long,boost::hash<DoublePair>> * NodeIDs;

	void assign_info(addSurfaceData * tmp) {info=tmp;};

	addSurfaceData * borrow_info(){return info;};

	// used by the interpolation algorithms
	// const escript::Data source;
    // escript::Data target;
    // const oxley::OxleyDomainRect_ptr other;

	// This is here to temporarily store information
private:
	addSurfaceData * info;
};

// class interpolationData
// {
// public:
// 	const escript::Data source;
//     escript::Data target;
//     const oxley::OxleyDomainRect_ptr other;

//     // interpolationData(const escript::Data source, 
//     // 						escript::Data target, 
//     // 				  const oxley::OxleyDomainRect_ptr other);
// };

// interpolationData::interpolationData(const escript::Data s, 
//                                      escript::Data t, 
//                                      const oxley::OxleyDomainRect_ptr o)
// {
//     source = s;
//     target = t;
//     other = o;
// }


class p8estData
{
public:
	// origin of domain
    double m_origin[3] = {0.0,0.0,0.0};

    // extent of the domain
    double m_lxyz[3] = {0.0};

    // side lengths of domain
    double m_length[3] = {0.0,0.0,0.0};

    // grid spacings / cell sizes of domain for each level of refinement
    double m_dx[3][P4EST_MAXLEVEL+1] = {{0}};

    // number of spatial subdivisions
    int m_NX[3] = {0,0,0};

    // total number of elements in each dimension
    escript::DataTypes::dim_t m_gNE[3] = {0,0,0};

    // number of elements for this rank in each dimension including shared
    escript::DataTypes::dim_t m_NE[3] = {0,0,0};

    // periodic boundary conditions
    bool periodic[3] {false, false, false};

	// maximum levels of recursion to use during refinement
	int max_levels_refinement = 0;
	double refinement_depth=0.0;
	double refinement_boundaries[6]={0.0};
	escript::Data mask; // a pointer to a mask

	// Pointer to the current solution and Node ID info
	std::unordered_map<long,double> * current_solution;
	std::unordered_map<DoubleTuple,long,boost::hash<DoubleTuple>> * NodeIDs;

	void assign_info(addSurfaceData * tmp) {info=tmp;};

	addSurfaceData * borrow_info(){return info;};

	// used by the interpolation algorithms
	const escript::Data* source;
    escript::Data * target;
    const oxley::OxleyDomainBrick_ptr * other;

	// This is here to temporarily store information
private:
	addSurfaceData * info;
};

// This structure temporarily stores information used by the addSurface function
struct addSurfaceData {

	int oldTag = -1;
	int newTag = -1;
	std::vector<double> x;
	std::vector<double> y;
	std::vector<double> z;

	// The domain in which the function z[x,y] is defined
	double xmin,xmax,ymin,ymax;

};

struct quad_info {
	int level;
	//bottom left coordinate
	double x;
	double y;
};

struct oct_info {
	int level;
	//bottom left coordinate
	double x;
	double y;
	double z;
};


struct update_RC_data {

	std::unordered_map<DoublePair,long,boost::hash<DoublePair>> * pNodeIDs; 
	// std::unordered_map<long,bool> * phangingNodeIDs; 
	p4est_t * p4est;
	std::vector< std::vector<long> > * indices;
	double m_origin[2]={0};

	std::vector<quad_info> * pQuadInfo;
};

struct getConnections_data {

	const std::unordered_map<DoublePair,long,boost::hash<DoublePair>> * pNodeIDs; 
	p4est_t * p4est;
	std::vector< std::vector<escript::DataTypes::index_t> > * indices;
	double m_origin[2]={0};
};

struct update_RC_data_brick {

	std::unordered_map<DoubleTuple,long,boost::hash<DoubleTuple>> * pNodeIDs; 
	// std::unordered_map<long,bool> * phangingNodeIDs; 
	p8est_t * p8est;
	std::vector<oxley::IndexVector> * indices;
	double m_origin[3]={0};

	std::vector<oct_info> * pOctInfo;
};

struct getConnections_data_brick {

	const std::unordered_map<DoubleTuple,long,boost::hash<DoubleTuple>> * pNodeIDs; 
	p8est_t * p8est;
	std::vector< std::vector<escript::DataTypes::index_t> > * indices;
	double m_origin[3]={0};
};

// Tracks information used by the assembler
template<class Scalar>
struct assembly_data_d {

	std::vector<Scalar> * EM_S;
	std::vector<Scalar> * EM_F;
};

// bool operator==(const p4estData A, const p4estData B)
// {
	
// 	return (A.m_origin[0] == B.m_origin[0])
// 		&& (A.m_origin[1] == B.m_origin[1])
// 		&& (A.m_length[0] == B.m_length[0])
// 		&& (A.m_length[1] == B.m_length[1]);
// };

// bool operator==(const p8estData A, const p8estData B)
// {
	
// 	return (A.m_origin[0] == B.m_origin[0])
// 		&& (A.m_origin[1] == B.m_origin[1])
// 		&& (A.m_origin[2] == B.m_origin[2])
// 		&& (A.m_length[0] == B.m_length[0])
// 		&& (A.m_length[1] == B.m_length[1])
// 		&& (A.m_length[2] == B.m_length[2]);
// };

template <typename S> 
struct interpolateNodesOnElementsWorker_Data {

	S sentinel;
	int offset;
	double * fxx;

};

template <typename S> 
struct interpolateNodesOnFacesWorker_Data {

	S sentinel;
	int offset;
	double * fxx;
	int direction=-1;
	bool shared=false;

};

namespace oxley {

// Call back function that copies quadrant tags onto tagVector
void getQuadTagVector(p4est_iter_volume_info_t * info, void *tagVector);
void getQuadTagVector(p8est_iter_volume_info_t * info, void *tagVector);

// Call back function that copies coordinate info onto tagVector
void getXCoordVector(p4est_iter_volume_info_t * info, void *tagVector);
void getYCoordVector(p4est_iter_volume_info_t * info, void *tagVector);
void getXCoordVector(p8est_iter_volume_info_t * info, void *tagVector);
void getYCoordVector(p8est_iter_volume_info_t * info, void *tagVector);
void getZCoordVector(p8est_iter_volume_info_t * info, void *tagVector);

// Call back function that copies node information onto tagVector
void getNodeTagVector(p4est_iter_volume_info_t * info, void *tagVector);
void getNodeTagVector(p8est_iter_volume_info_t * info, void *tagVector);

} //namespace oxley

#endif
