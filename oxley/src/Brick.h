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

#ifndef __OXLEY_BRICK_H__
#define __OXLEY_BRICK_H__

#include <unordered_map>
#include <utility>

#include <boost/functional/hash.hpp>

#include <escript/Data.h>
#include <escript/EsysMPI.h>

#include <oxley/AbstractAssembler.h>
#include <oxley/Oxley.h>
#include <oxley/OxleyDomain.h>
#include <oxley/OxleyData.h>
#include <oxley/RefinementType.h>
#include <oxley/RefinementZone.h>

#include <oxley/tictoc.h>

#include "p4est/p8est_io.h"
#include "p4est/p8est.h"
#include "p4est/p8est_connectivity.h"
#include "p4est/p8est_lnodes.h"

#include <boost/python.hpp>
#ifdef ESYS_HAVE_BOOST_NUMPY
#include <boost/python/numpy.hpp>
#endif

using namespace boost::python;

namespace oxley {


#define NORMAL 1;
#define HANGINGEDGE 2;
#define HANGINGFACE 3;


/**
   \brief
   Brick is the 2-dimensional implementation of an Oxleydomain.
*/
class Brick: public OxleyDomain
{

    template<class Scalar> friend class DefaultAssembler3D;

public:

    /**
       \brief creates a rectangular mesh with n0 x n1 x n2 elements over the
              domain [x0,x1] x [y0,y1] x [z0,z1].
       \param
    */
    // Brick();

#if 0  // DEPRECATED: mpiInfo must be handed over from caller
    Brick(int order, dim_t n0, dim_t n1, dim_t n2, double x0, double y0, double z0,
      double x1, double y1, double z1, int d0, int d1, int d2,
      const std::vector<double>& points, const std::vector<int>& tags,
      const TagMap& tagnamestonums,
      int periodic0, int periodic1, int periodic2);
#endif

    /**
       \brief creates a brick mesh with custom MPI communicator
       \param jmpi MPI communicator info from caller (required)
    */
    Brick(escript::JMPI jmpi, int order, dim_t n0, dim_t n1, dim_t n2,
      double x0, double y0, double z0, double x1, double y1, double z1,
      int d0, int d1, int d2,
      const std::vector<double>& points, const std::vector<int>& tags,
      const TagMap& tagnamestonums,
      int periodic0, int periodic1, int periodic2);

    // DANGEROUS: If update is false then the mesh is not properly initialised
    Brick(oxley::Brick& B, int order, bool update);

    /**
       \brief creates a rectangular mesh from numpy arrays [x,y].
            Requires boost numpy
       \param
    */
#ifdef ESYS_HAVE_BOOST_NUMPY
    Brick(int order, dim_t n0, dim_t n1, dim_t n2, 
      boost::python::numpy::ndarray x, boost::python::numpy::ndarray y, boost::python::numpy::ndarray z);
#endif

    /**
       \brief
       Destructor.
    */
    ~Brick();

    /**
       \brief
       returns a description for this domain
    */
    virtual std::string getDescription() const;

    /**
       \brief
       dumps the mesh to a file with the given name
       \param filename The name of the output file
    */
    virtual void dump(const std::string& filename) const;

    /**
       \brief
       writes the current mesh to a file with the given name
       \param filename The name of the file to write to
    */
    virtual void write(const std::string& filename) const;

    /**
       \brief equality operator
    */
    virtual bool operator==(const escript::AbstractDomain& other) const;

    /**
       \brief
       interpolates data given on source onto target where source and target
       are given on different domains
    */
    virtual void interpolateAcross(escript::Data& target,
                                   const escript::Data& source) const;

    /**
       \brief
       returns true if this rank owns the sample id.
    */
    virtual bool ownSample(int fs_code, index_t id) const;

    /**
       \brief
       returns the number of data points summed across all MPI processes
    */
    virtual dim_t getNumDataPointsGlobal() const;

    /**
       \brief
       copies the surface normals at data points into out. The actual function
       space to be considered is defined by out. out has to be defined on this
       domain.
    */
    virtual void setToNormal(escript::Data& out) const;

    /**
       \brief
       copies the size of samples into out. The actual function space to be
       considered is defined by out. out has to be defined on this domain.
    */
    virtual void setToSize(escript::Data& out) const;

    /**
     * \brief
       Returns a Data object filled with random data passed through filter.
    */
    virtual escript::Data randomFill(const escript::DataTypes::ShapeType& shape,
       const escript::FunctionSpace& what, long seed, const boost::python::tuple& filter) const;

    /**
       \brief
       returns the array of reference numbers for a function space type
       \param fsType The function space type
    */
    const dim_t* borrowSampleReferenceIDs(int fsType) const;

    /**
       \brief
       writes the mesh to a VTK file
       \param filename The file name
    */
    virtual void writeToVTK(std::string filename, bool writeMesh) const;

    /**
       \brief
       writes the mesh to file
    */
    virtual void saveMesh(std::string filename) ;

    #ifdef ESYS_HAVE_TRILINOS
    /**
       \brief
       writes the mesh to file
    */
    virtual void loadMesh(std::string filename) ;

    /**
       \brief
       refines the mesh
       \param algorithmname The algorithm to use
    */
    virtual void refineMesh(std::string algorithmname);

    /**
       \brief
       refines the mesh near a boundary
       \param maxRecursion Max levels of recursion
       \param algorithmname The algorithm to use
    */
    virtual void refineBoundary(std::string boundary, double dx);

    /**
       \brief
       refines the mesh within the interior of a region bound by 
       x0, x1, y0, y1
       \param x0 boundary of the region
       \param x1 boundary of the region
       \param y0 boundary of the region
       \param y1 boundary of the region
    */
    virtual void refineRegion(double x0, double x1, double y0, double y1, double z0, double z1);

        /**
       \brief
       refines the mesh around the point x0, y0
       \param x0 spatial coordinate of point
       \param y0 spatial coordinate of point
    */
    virtual void refinePoint(double x0, double y0, double z0);

    /**
       \brief
       refines a circle on the mesh with center x0, y0 and radius r
       \param x0 spatial coordinate of center of the circle
       \param y0 spatial coordinate of center of the circle
       \param r radius of the circle
    */
    virtual void refineSphere(double x0, double y0, double z0, double r);
    #endif //ESYS_HAVE_TRILINOS

    /**
       \brief
         refines a region defined by a mask
    */
    virtual void refineMask(escript::Data mask);

    /**
       \brief
       sets the number of levels of refinement
    */
    virtual void setRefinementLevels(int refinementlevels)
    {
        forestData->max_levels_refinement = refinementlevels;
    };

    /**
       \brief
       sets the number of levels of refinement
    */
    virtual int getRefinementLevels() const
    {
        return forestData->max_levels_refinement;
    };

    /**
       \brief
       returns a Data object containing the coordinate information
    */
    virtual escript::Data getX() const;

    /**
       \brief
       returns the number of vertices (int)
    */
    int getNumVertices() const { return connectivity->num_vertices;};

    virtual dim_t findNode(const double *coords) const;

    /**
      \brief
      sets adaptive refinement on or off
    */
    void setAdaptiveRefinement(bool status) { adaptive_refinement = status ; } ;

    /**
       \brief equality operator
    */
    // virtual bool operator==(const escript::AbstractDomain& other) const;

    /**
       \brief inequality operator
    */
    // virtual bool operator!=(const escript::AbstractDomain& other) const {
    //     return !(operator==(other));
    // }

    // These functions are used internally
    p8est_t * borrow_p8est() const { return p8est;};

    p8estData * borrow_forestData() { return forestData;};

    p8est_connectivity_t * borrow_connectivity() const { return connectivity; };

    void * borrow_temp_data() { return temp_data; };

    void set_temp_data(addSurfaceData * x) { temp_data = x; };

    void clear_temp_data() { free(temp_data); };

    void print_debug_report(std::string);

    // /**
    //    \brief
    //    returns a Data object containing the coordinate information
    // */
    // int getNumVertices() const { return connectivity->num_vertices;};

    // virtual dim_t findNode(const double *coords) const;

    // virtual void nodesToDOF(escript::Data& out, const escript::Data& in) const;
    
    // virtual dim_t getDofOfNode(dim_t node) const;

    // Used by weipa
    const long getNodeId(double x, double y, double z);

    /**
       \brief
       returns a vector of rank numbers where vec[i]=n means that rank n
       'owns' element/face element i.
    */
    virtual RankVector getOwnerVector(int fsType) const;

    /**
       \brief
       returns the number of face elements in the order
       (left,right,bottom,top) on current MPI rank
    */
    virtual const dim_t* getNumFacesPerBoundary() const { return m_faceCount; }

    // This is not private as it is used by weipa
    // A p8est
    p8est_t * p8est;
    std::unordered_map<DoubleTuple,long,boost::hash<DoubleTuple>> NodeIDs; //global ids of the nodes

    /**
       \brief
       Returns the ID numbers of the neighbouring four nodes
    */
    void getNeighouringNodeIDs(int8_t level, p8est_qcoord_t x, p8est_qcoord_t y, p8est_qcoord_t z, 
                                             p8est_topidx_t treeid, long (&ids) [8]) const;

    /**
       \brief
       Returns true if the node is on the boundary
    */
    bool isBoundaryNode(p8est_quadrant_t * quad, int n, p8est_topidx_t treeid, p8est_qcoord_t length) const;

    /**
       \brief
       Returns true if the node is on the top or right boundaries
    */
    bool isUpperBoundaryNode(p8est_quadrant_t * quad, int n, p8est_topidx_t treeid, p8est_qcoord_t length) const;

    /**
       \brief
       Returns true if the node is on the bottom or left boundaries
    */
    bool isLowerBoundaryNode(p8est_quadrant_t * quad, int n, p8est_topidx_t treeid, p8est_qcoord_t length) const;

    /**
       \brief
       Returns true if the node is on the left boundary
    */
    bool isLeftBoundaryNode(p8est_quadrant_t * quad, int n, p8est_topidx_t treeid, p8est_qcoord_t length) const;

    /**
       \brief
       Returns true if the node is on the right boundary
    */
    bool isRightBoundaryNode(p8est_quadrant_t * quad, int n, p8est_topidx_t treeid, p8est_qcoord_t length) const;

    /**
       \brief
       Returns true if the node is on the bottom boundary
    */
    bool isBottomBoundaryNode(p8est_quadrant_t * quad, int n, p8est_topidx_t treeid, p8est_qcoord_t length) const;

    /**
       \brief
       Returns true if the node is on the top boundary
    */
    bool isTopBoundaryNode(p8est_quadrant_t * quad, int n, p8est_topidx_t treeid, p8est_qcoord_t length) const;

    /**
       \brief
       Returns true if the node is on the bottom boundary
    */
    bool isAboveBoundaryNode(p8est_quadrant_t * quad, int n, p8est_topidx_t treeid, p8est_qcoord_t length) const;

    /**
       \brief
       Returns true if the node is on the top boundary
    */
    bool isBelowBoundaryNode(p8est_quadrant_t * quad, int n, p8est_topidx_t treeid, p8est_qcoord_t length) const;

    /**
      \brief
      Applies a refinementzone
   */
    escript::Domain_ptr apply_refinementzone(RefinementZone R);

////////////////////////////////
private:

    // A ghost
    p8est_ghost_t * ghost;

    // The data structure in p8est
    p8estData * forestData;

    // This object records the connectivity of the p8est octants
    p8est_connectivity_t * connectivity;

    // This structure records the node numbering information
    p8est_lnodes * nodes;
    long * nodeIncrements;

    // Indices
    std::vector<IndexVector> * indices;

    // Pointer that records the location of a temporary data structure
    void * temp_data;

    // Brick needs to keep track of this information
    std::unordered_map<DoubleTuple,long,boost::hash<DoubleTuple>> treeIDs; //global ids of the hanging nodes
    std::vector<long> octantIDs; // IDs of the octants
    std::vector<oct_info> octantInfo;

    std::vector<borderNodeInfo> NodeIDsTop;
    std::vector<borderNodeInfo> NodeIDsBottom;
    std::vector<borderNodeInfo> NodeIDsLeft;
    std::vector<borderNodeInfo> NodeIDsRight;
    std::vector<borderNodeInfo> NodeIDsAbove;
    std::vector<borderNodeInfo> NodeIDsBelow;

    std::vector<DoubleTuple> NormalNodes;
    std::vector<long> HangingFaceNodes;
    std::vector<long> HangingEdgeNodes;
    std::vector<int> is_hanging; // element x is true if node id x is a hanging node

    std::vector<hangingFaceInfo> hanging_face_orientation;
    std::vector<hangingEdgeInfo> hanging_edge_orientation;

    std::vector<std::vector<long>> neighbours;

    // Row and column indices in CRS format
    IndexVector myRows;
    IndexVector myColumns;

    // vector that maps each node to a DOF index (used for the coupler)
    IndexVector m_dofMap;

        /// faceOffset[i]=-1 if face i is not an external face, otherwise it is
    /// the index of that face (where i: 0=left, 1=right, 2=bottom, 3=top)
    IndexVector m_faceOffset;

    // 
    IndexVector m_nodeId;

    // tolerance used when comparing doubletuples
    double tuple_tolerance=0.0;

    p8est_connectivity_t *
    new_brick_connectivity (int n0, int n1, int n2, int periodic_a, int periodic_b, int periodic_c,
                               double x0, double x1, double y0, double y1, double z0, double z1);

    virtual Assembler_ptr createAssembler(std::string type,
                                          const DataMap& options) const;

    virtual void updateMeshInformation();

    /**
      \brief
      Returns the ID of a quad from the ID of it's bottom left node
    */
    long getQuadID(long nodeid) const;

    template<typename Scalar>
    void assembleIntegrateImpl(std::vector<Scalar>& integrals, const escript::Data& arg) const;

////////////////////////////////
protected:
    
    /**
       \brief
       Returns the number of nodes
    */
    virtual dim_t getNumNodes() const;

    /**
       \brief
       Returns the number of nodes
    */
    virtual dim_t getNumHangingNodes() const;

    /**
       \brief
       Returns the number of elements
    */
    virtual dim_t getNumElements() const;

    /**
       \brief
       Returns the number of face elements
    */
    dim_t getNumFaceElements() const;

    /**
       \brief
       Returns the number of degrees of freedom
    */
    inline dim_t getNumDOF() const;

    /**
       \brief
       Returns true if the face is hanging
    */
    bool isHangingFace(p8est_lnodes_code_t face_code, int n) const;

    /**
       \brief
       Returns true if the node is hanging
    */
    bool isHangingNode(p8est_lnodes_code_t face_code, int n) const;

    /**
       \brief
       Returns true if the node is hanging
    */
    bool getHangingInfo(p8est_lnodes_code_t face_code, int hanging_faces[], int hanging_edges[]) const;

    /**
       \brief
       Returns true if the node is hanging
    */
    bool hasHanging(p8est_lnodes_code_t face_code) const;

    /**
       \brief
       Returns true if the node is hanging
    */
    bool onUpperBoundary(double x, double y, double z) const;

    /**
       \brief
       Returns true if the node is hanging
    */
    int getNumHangingNodes() { return num_hanging; };

    /**
       \brief
       Regenerates the Ghost information
    */
    void reset_ghost();

    /**
       \brief
       Updates NodeIncrements
    */
    void updateNodeIncrements();

    /**
       \brief
       Renumbers the nodes (called after refinement)
    */
    void renumberNodes();

    /**
       \brief
       Returns true if vector has duplicate point. 
    */
    bool gotPoint(DoubleTuple point, std::vector<DoubleTuple> vector);

    /**
       \brief
       Returns true if vector has duplicate point. 
    */
    bool hasDuplicate(DoubleTuple point, std::vector<DoubleTuple> vector, bool serial );

    /**
       \brief
       Returns true if vector has point. 
    */
    bool gotAlready(DoubleTuple point, std::vector<DoubleTuple> vector);

    /**
       \brief
       A faster version of p8est_qcoord_to_vertex
    */
    void p8est_qcoord_to_vertex_fast(p8est_connectivity_t * connectivity,
                        p4est_topidx_t treeid,
                        p4est_qcoord_t x, p4est_qcoord_t y, p4est_qcoord_t z,
                        double vxyz[3]);

    /**
       \brief
       Updates the NodeIDs of the hanging Nodes
    */
    void updateHangingNodeIDs();

    /**
       \brief
       Updates myRows and myColumns (used to store node connectivity information Yale Formay)
    */
    void updateRowsColumns();

    /**
       \brief
       Updates TreeIDs
    */
    void updateTreeIDs();

    /**
     * \brief
     * Updates the mesh after refinement
   */
   void updateMesh();

   /**
     * \brief
     * Updates the mesh after refinement
   */
   void updateMeshBackend();


   /**
    * \brief
    * Toggles automatic mesh updates in the refinement functions
    */
   void AutomaticMeshUpdateOnOff(bool autoMeshUpdates);
   bool autoMeshUpdates = true;

    /**
       \brief
       Updates the NodeID information
    */
    virtual void assembleCoordinates(escript::Data& arg) const;
    virtual void assembleGradient(escript::Data& out, const escript::Data& in) const;


    virtual void assembleIntegrate(std::vector<real_t>& integrals,
                                   const escript::Data& arg) const;
    virtual void assembleIntegrate(std::vector<cplx_t>& integrals,
                                   const escript::Data& arg) const;
    virtual std::vector<IndexVector> getConnections(bool includeShared=false) const;

#ifdef ESYS_HAVE_TRILINOS
    virtual esys_trilinos::TrilinosGraph_ptr getTrilinosGraph() const;
#endif
#ifdef ESYS_HAVE_PASO
    virtual paso::SystemMatrixPattern_ptr getPasoMatrixPattern(bool reducedRowOrder, bool reducedColOrder) const;
#endif

#ifdef ESYS_HAVE_PASO
    // the Paso System Matrix pattern
    mutable paso::SystemMatrixPattern_ptr m_pattern;
#endif

#ifdef ESYS_HAVE_TRILINOS
    /// Trilinos graph structure, cached for efficiency
    mutable esys_trilinos::TrilinosGraph_ptr m_graph;
#endif
    
    // INTERPOLATION


    virtual void interpolateNodesOnElements(escript::Data& out, const escript::Data& in, bool reduced) const;   
    virtual void interpolateNodesOnFaces(escript::Data& out, const escript::Data& in, bool reduced) const;
    virtual void nodesToDOF(escript::Data& out, const escript::Data& in) const;
    virtual dim_t getDofOfNode(dim_t node) const;
    virtual void populateSampleIds();
    
    // INTERPOLATION (FROM COARSE TO FINE MESHES AND VICE VERSA)

    /**
       \brief
       Checks that the given interpolation is possible, throws and OxleyException if it is not.
    */
    void validateInterpolationAcross(int fsType_source, const escript::AbstractDomain& domain, int fsType_target) const;
    void interpolateNodesToNodesFiner(const escript::Data& source, escript::Data& target, const Brick& other) const;
    void interpolateNodesToElementsFiner(const escript::Data& source, escript::Data& target, const Brick& other) const;
    void interpolateElementsToElementsCoarser(const escript::Data& source, escript::Data& target, const Brick& other) const;
    void interpolateElementsToElementsFiner(const escript::Data& source, escript::Data& target, const Brick& other) const;
    void interpolateReducedToElementsFiner(const escript::Data& source, escript::Data& target, const Brick& other) const;
    void interpolateReducedToReducedFiner(const escript::Data& source, escript::Data& target, const Brick& other) const;

    template <typename S>
    void interpolateNodesOnElementsWorker(escript::Data& out,
                                  const escript::Data& in, bool reduced, S sentinel) const;   
    template <typename S>
    void interpolateNodesOnFacesWorker(escript::Data& out,
                                         const escript::Data& in,
                                         bool reduced, S sentinel) const; 

    template<typename Scalar>
    void assembleGradientImpl(escript::Data& out,
                              const escript::Data& in) const;

    template<typename Scalar> void addToMatrixAndRHS(escript::AbstractSystemMatrix* S, escript::Data& F,
           const std::vector<Scalar>& EM_S, const std::vector<Scalar>& EM_F,
           bool addS, bool addF, borderNodeInfo quad, int nEq=1, int nComp=1) const;
    template<typename Scalar> void addToMatrixAndRHS(escript::AbstractSystemMatrix* S, escript::Data& F,
           const std::vector<Scalar>& EM_S, const std::vector<Scalar>& EM_F,
           bool addS, bool addF, index_t e, index_t t, int nEq=1, int nComp=1) const;

    // Updates m_faceOffset for each quadrant
    void updateFaceOffset();

    /**
       \brief
       Updates m_faceCount
    */
    void updateFaceElementCount();

    // vector with first node id on each rank
    IndexVector m_nodeDistribution;
    void updateNodeDistribution();
    IndexVector m_elementId;
    void updateElementIds();
    IndexVector m_faceId;


// #ifdef ESYS_HAVE_PASO
//     // the Paso System Matrix pattern
//     mutable paso::SystemMatrixPattern_ptr m_pattern;
// #endif

    IndexVector getNodeDistribution() const;

    void addPoints(const std::vector<double>& coords, const std::vector<int>& tags);

    // Initial number of nodes
    int m_NN[3] = {0};
    // Initial number of divisions
    long m_NE[3] = {0};
    // Initial spacing
    double m_NX[3] = {0};
    /// number of face elements per edge (left, right, bottom, top)
    dim_t m_faceCount[8];

    // The number of hanging nodes in the mesh
    int num_hanging;

    // Timer used to profile code
    // TicTocClock oxleytimer;

    // Modified versions of the p4est library functions that use openMP
    int p8est_connectivity_is_valid_fast(p8est_connectivity_t * conn); 

};

// typedef POINTER_WRAPPER_CLASS(Brick) OxleyDomainBrick_ptr;

} //end namespace


#endif //__OXLEY_BRICK_H__
