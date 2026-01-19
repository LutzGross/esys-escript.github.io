/*****************************************************************************
*
* Copyright (c) 2003-2019 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0f
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
*
*****************************************************************************/

#ifndef __OXLEY_RECTANGLE_H__
#define __OXLEY_RECTANGLE_H__

#include <unordered_map>
#include <utility>

#include <boost/functional/hash.hpp>

#include <escript/Data.h>
#include <escript/EsysMPI.h>

#include <oxley/AbstractAssembler.h>
#include <oxley/Oxley.h>
#include <oxley/OxleyData.h>
#include <oxley/OxleyDomain.h>
#include <oxley/RefinementZone.h>

#include <oxley/tictoc.h>

#include "p4est/p4est.h"
#include "p4est/p4est_connectivity.h"
#include "p4est/p4est_lnodes.h"

#include <boost/python.hpp>
#ifdef ESYS_HAVE_BOOST_NUMPY
#include <boost/python/numpy.hpp>
#endif

using namespace boost::python;

namespace oxley {

//forward declaration
class interpolateWorker_Data;

/**
   \brief
   Rectangle is the 2-dimensional implementation of a Oxleydomain.
*/
class Rectangle: public OxleyDomain
{

    template<class Scalar> friend class DefaultAssembler2D;

public:

    /**
       \brief creates a rectangular mesh with n0 x n1 elements over the
              rectangle [x0,x1] x [y0,y1].
       \param
    */
#if 0  // DEPRECATED: mpiInfo must be handed over from caller
    Rectangle(int order, dim_t n0, dim_t n1,
        double x0, double y0, double x1, double y1,
        int d0, int d1,
        const std::vector<double>& points, const std::vector<int>& tags,
        const TagMap& tagnamestonums,
        int periodic0, int periodic1);
#endif

    /**
       \brief creates a rectangular mesh with n0 x n1 elements over the
              rectangle [x0,x1] x [y0,y1] with a custom MPI communicator.
       \param jmpi MPI communicator info from caller (required)
    */
    Rectangle(escript::JMPI jmpi, int order, dim_t n0, dim_t n1,
        double x0, double y0, double x1, double y1,
        int d0, int d1,
        const std::vector<double>& points, const std::vector<int>& tags,
        const TagMap& tagnamestonums,
        int periodic0, int periodic1);

    /**
       \brief creates a rectangular mesh from numpy arrays [x,y].
            Requires boost numpy
       \param
    */
#ifdef ESYS_HAVE_BOOST_NUMPY
    Rectangle(int order, dim_t n0, dim_t n1, boost::python::numpy::ndarray x, boost::python::numpy::ndarray y);
#endif

    /**
       \brief creates a rectangular mesh from an existing Rectangle.
       \param
    */
    Rectangle(const oxley::Rectangle& rect, int order);

    /**
       \brief
       Destructor.
    */
    ~Rectangle();

    /**
       \brief
       returns a description for this domain
    */
    virtual std::string getDescription() const;

    /**
       \brief
       writes the current mesh to a file with the given name
       \param filename The name of the file to write to
    */
    virtual void write(const std::string& filename) const;

    /**
       \brief
       dumps the mesh to a file with the given name
       \param filename The name of the output file
    */
    virtual void dump(const std::string& filename) const;

    /**
       \brief equality operator
    */
    virtual bool operator==(const escript::AbstractDomain& other) const;

    /**
       \brief
       interpolates data given on source onto target where source and target
       are given on different domains
    */
    virtual void interpolateAcross(escript::Data& target, const escript::Data& source) const;

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

    escript::Data randomFillWorker(const escript::DataTypes::ShapeType& shape,
       long seed, const boost::python::tuple& filter) const;    

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
       \param writeMesh whether to write tag information
    */
    virtual void writeToVTK(std::string filename, bool writeMesh) const;

   #ifdef ESYS_HAVE_TRILINOS
    /**
       \brief
       writes the mesh to file
    */
    virtual void saveMesh(std::string filename);

    /**
       \brief
       writes the mesh to file
    */
    virtual void loadMesh(std::string filename);

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
    virtual void refineRegion(double x0, double x1, double y0, double y1);

    /**
       \brief
       refines the mesh around the point x0, y0
       \param x0 spatial coordinate of point
       \param y0 spatial coordinate of point
    */
    virtual void refinePoint(double x0, double y0);

    /**
       \brief
       refines a circle on the mesh with center x0, y0 and radius r
       \param x0 spatial coordinate of center of the circle
       \param y0 spatial coordinate of center of the circle
       \param r radius of the circle
    */
    virtual void refineCircle(double x0, double y0, double r);
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
        forestData.max_levels_refinement = refinementlevels;
    };

    /**
       \brief
       sets the number of levels of refinement
    */
    virtual int getRefinementLevels() const
    {
        return forestData.max_levels_refinement;
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
       \brief inequality operator
    */
    // virtual bool operator!=(const escript::AbstractDomain& other) const {
    //     return !(operator==(other));
    // }

    // These functions are used internally and are not exposed to python
    p4est_t * borrow_p4est() { return p4est;};
    p4estData * borrow_forestData() { return &forestData;};
    p4est_connectivity_t * borrow_connectivity()  { return connectivity; };
    void * borrow_temp_data() { return temp_data; };
    void set_temp_data(addSurfaceData * x) { temp_data = x; };
    void clear_temp_data() { free(temp_data); };
    void print_debug_report(std::string);

    // virtual void finaliseRhs(escript::Data& rhs);

    // Used by weipa
    const long getNodeId(double x, double y);

    /**
       \brief
       returns the number of face elements in the order
       (left,right,bottom,top) on current MPI rank
    */
    virtual const dim_t* getNumFacesPerBoundary() const { return m_faceCount; }

    /**
       \brief
       returns a vector of rank numbers where vec[i]=n means that rank n
       'owns' element/face element i.
    */
    virtual RankVector getOwnerVector(int fsType) const;

    // not private as it is needed by weipa
    // A p4est
    p4est_t * p4est;

    // Rectangle needs to keep track of this information
    std::unordered_map<DoublePair,long,boost::hash<DoublePair>> NodeIDs; //global ids of the nodes

        /**
       \brief
       Returns true if the node is on the left boundary
    */
    bool isLeftBoundaryNode(p4est_quadrant_t * quad, int n, p4est_topidx_t treeid, p4est_qcoord_t length) const;

    /**
       \brief
       Returns true if the node is on the right boundary
    */
    bool isRightBoundaryNode(p4est_quadrant_t * quad, int n, p4est_topidx_t treeid, p4est_qcoord_t length) const;

    /**
       \brief
       Returns true if the node is on the bottom boundary
    */
    bool isBottomBoundaryNode(p4est_quadrant_t * quad, int n, p4est_topidx_t treeid, p4est_qcoord_t length) const;

    /**
       \brief
       Returns true if the node is on the top boundary
    */
    bool isTopBoundaryNode(p4est_quadrant_t * quad, int n, p4est_topidx_t treeid, p4est_qcoord_t length) const;

    /**
       \brief
       Returns the ID numbers of the neighbouring four nodes
    */
    void getNeighouringNodeIDs(int8_t level, p4est_qcoord_t x, p4est_qcoord_t y, p4est_topidx_t treeid, long (&ids) [4]) const;

    /**
       \brief
       Returns true if the node is on the border and should be treated as hanging
    */
    bool checkHangingBorderNode(p4est_quadrant_t * quad, p4est_qcoord_t x, p4est_qcoord_t y, 
                                 p4est_topidx_t treeid, int n) const;

    /**
       \brief
       Returns the face code of the hanging border node
       (cf. p6est_lnodes.h lines 57-116)
    */
    int getHangingBorderNodeFacecode(p4est_quadrant_t * quad, int8_t level, p4est_qcoord_t x, p4est_qcoord_t y, 
                                 p4est_topidx_t treeid, int n, int boundary_code) const;

    /**
       \brief
       A version of p4est_qcoord_to_vertex that does not test for exceptions
       *** USE WITH CAUTION ***
    */
    void p4est_qcoord_to_vertex_mod (p4est_connectivity_t * connectivity,
                        p4est_topidx_t treeid,
                        p4est_qcoord_t x, p4est_qcoord_t y,
                        double vxyz[3]) const;

    /**
      \brief
      Applies a refinementzone
   */
    escript::Domain_ptr apply_refinementzone(RefinementZone R);

   /**
     * \brief
     * Updates the mesh after refinement
   */
   void updateMesh();

   /**
    * \brief
    * Toggles automatic mesh updates in the refinement functions
    */
   void AutomaticMeshUpdateOnOff(bool autoMeshUpdates);
   bool autoMeshUpdates = true;


    // Data used by the interpolation functions
    // const interpolationData interp_data;

private:
    // The data structure in p4est
    p4estData forestData;

    // Connectivity information
    p4est_connectivity_t * connectivity;

    // Node numbering
    p4est_lnodes_t * nodes;
    long nodeIncrements[MAXTREES] = {0};

    // Pointer that records the location of a temporary data structure
    void * temp_data;

    // std::unordered_map<long,bool> hangingNodeIDs; //global ids of the hanging nodes
    std::vector<bool> is_hanging; // element x is true if node id x is a hanging node
    // std::vector<std::vector<long>> is_hanging_face; // if face x-y is hanging then element x is y
    std::unordered_map<DoublePair,long,boost::hash<DoublePair>> treeIDs; //global ids of the hanging nodes
    std::vector<long> quadrantIDs; // IDs of the quadrants
    std::vector<quad_info> quadrantInfo;

    std::vector<borderNodeInfo> NodeIDsTop;
    std::vector<borderNodeInfo> NodeIDsBottom;
    std::vector<borderNodeInfo> NodeIDsLeft;
    std::vector<borderNodeInfo> NodeIDsRight;

    std::vector<hangingNodeInfo> hanging_face_orientation;
    
    // Row and column indices in CRS format
    IndexVector myRows;
    IndexVector myColumns;

    // vector that maps each node to a DOF index (used for the coupler)
    IndexVector m_dofMap;

    // 
    IndexVector m_nodeId;

    // This is a modified version of the p4est library function new_connectivity
p4est_connectivity_t *
new_rectangle_connectivity(int mi, int ni, int periodic_a, int periodic_b, 
    double x0, double x1, double y0, double y1);
  
    /**
       \brief
       creates and returns an assembler of the requested type.
    */
    virtual Assembler_ptr createAssembler(std::string type, const DataMap& options) const;

    virtual void updateMeshInformation();

    /**
      \brief
      Returns the ID of a quad from the ID of it's bottom left node
    */
    long getQuadID(long nodeid) const;

    template<typename Scalar>
    void assembleIntegrateImpl(std::vector<Scalar>& integrals, const escript::Data& arg) const;


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

    // virtual dim_t getNumDOFInAxis(unsigned axis) const;
    // virtual index_t getFirstInDim(unsigned axis) const;
    // virtual IndexVector getDiagonalIndices(bool upperOnly) const;

    /**
       \brief
       Returns true if the node is on the boundary
    */
    bool isBoundaryNode(p4est_quadrant_t * quad, int n, p4est_topidx_t treeid, p4est_qcoord_t length) const;

    /**
       \brief
       Returns true if the node is on the top or right boundaries
    */
    bool isUpperBoundaryNode(p4est_quadrant_t * quad, int n, p4est_topidx_t treeid, p4est_qcoord_t length) const;

    /**
       \brief
       Returns true if the node is on the bottom or left boundaries
    */
    bool isLowerBoundaryNode(p4est_quadrant_t * quad, int n, p4est_topidx_t treeid, p4est_qcoord_t length) const;

    /**
       \brief
       Returns true if the face is hanging
    */
    bool isHangingFace(p4est_lnodes_code_t face_code, int n) const;

    /**
       \brief
       Returns true if the node is hanging
    */
    bool isHangingNode(p4est_lnodes_code_t face_code, int n) const;

    /**
       \brief
       Returns true if the node is hanging
    */
    bool getHangingNodes(p4est_lnodes_code_t face_code, int hanging[]) const;

    /**
       \brief
       Returns true if the node is hanging
    */
    int getNumHangingNodes() { return num_hanging; };

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
       \brief
       Renumbers the hanging nodes (called after refinement)
       (currently not used anywhere in the code)
    */
    // void renumberHangingNodes();

    /**
       \brief
       Updates the NodeID information
    */
    virtual void assembleCoordinates(escript::Data& arg) const;

    virtual void assembleGradient(escript::Data& out, const escript::Data& in) const;

    /**
       \brief
       Updates NodeIncrements
    */
    void updateQuadrantIDinformation();

    
    virtual void assembleIntegrate(std::vector<real_t>& integrals, const escript::Data& arg) const;
    virtual void assembleIntegrate(std::vector<cplx_t>& integrals, const escript::Data& arg) const;
    virtual std::vector<IndexVector> getConnections(bool includeShared=false) const;
#ifdef ESYS_HAVE_TRILINOS
    virtual esys_trilinos::TrilinosGraph_ptr getTrilinosGraph() const;
#endif
#ifdef ESYS_HAVE_PASO
    virtual paso::SystemMatrixPattern_ptr getPasoMatrixPattern(bool reducedRowOrder, bool reducedColOrder) const;
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
    void interpolateNodesToNodesFiner(const escript::Data& source, escript::Data& target, const Rectangle& other) const;
    template <typename S> 
    void interpolateNodesToNodesWorker(const escript::Data& source, escript::Data& target, const Rectangle& other) const;
    void interpolateNodesToNodesCoarser(const escript::Data& source, escript::Data& target, const Rectangle& other)  const;
    void interpolateNodesToElementsFiner(const escript::Data& source, escript::Data& target, const Rectangle& other) const;
    void interpolateElementsToElementsCoarser(const escript::Data& source, escript::Data& target, const Rectangle& other) const;
    void interpolateElementsToElementsFiner(const escript::Data& source, escript::Data& target, const Rectangle& other) const;
    template <typename S> 
    void interpolateElementsToElementsWorker(const escript::Data& source, escript::Data& target, const Rectangle& other) const;
    void interpolateReducedToElementsFiner(const escript::Data& source, escript::Data& target, const Rectangle& other) const;
    void interpolateReducedToReducedFiner(const escript::Data& source, escript::Data& target, const Rectangle& other) const;
   
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


#ifdef ESYS_HAVE_PASO
    // the Paso System Matrix pattern
    mutable paso::SystemMatrixPattern_ptr m_pattern;
#endif

#ifdef ESYS_HAVE_TRILINOS
    /// Trilinos graph structure, cached for efficiency
    mutable esys_trilinos::TrilinosGraph_ptr m_graph;
#endif

    IndexVector getNodeDistribution() const;

    void addPoints(const std::vector<double>& coords, const std::vector<int>& tags);

    // Initial number of nodes
    int m_NN[2] = {0};
    // Initial number of divisions
    long m_NE[2] = {0};
    // Initial spacing
    double m_NX[2] = {0};
    /// number of face elements per edge (left, right, bottom, top)
    dim_t m_faceCount[4];

    /// faceOffset[i]=-1 if face i is not an external face, otherwise it is
    /// the index of that face (where i: 0=left, 1=right, 2=bottom, 3=top)
    IndexVector m_faceOffset;

    // The number of hanging nodes in the mesh
    int num_hanging;

    TicTocClock oxleytimer;


};

// typedef POINTER_WRAPPER_CLASS(Rectangle) OxleyDomainRect_ptr;

// class interpolateWorker_Data {
// public:
//    const escript::Data * source;
//    const escript::Data * target;
//    const oxley::Rectangle * other;

//    interpolateWorker_Data(const escript::Data * source, 
//                           const escript::Data * target, 
//                           const oxley::Rectangle * other);
//    ~interpolateWorker_Data();

// };



} //end namespace


#endif //__OXLEY_RECTANGLE_H__
