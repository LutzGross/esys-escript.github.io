/*****************************************************************************
*
* Copyright (c) 2003-2019 by The University of Queensland
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

#ifndef __OXLEY_BRICK_H__
#define __OXLEY_BRICK_H__

#include <escript/EsysMPI.h>
#include <escript/SubWorld.h>

#include <oxley/OxleyDomain.h>
#include <oxley/OxleyData.h>

#include <p8est.h>
#include <p8est_connectivity.h>
#include <p8est_lnodes.h>

#include <boost/python.hpp>

using namespace boost::python;

namespace oxley {

/**
   \brief
   Brick is the 2-dimensional implementation of an Oxleydomain.
*/
class Brick: public OxleyDomain
{
public:

    /**
       \brief creates a rectangular mesh with n0 x n1 x n2 elements over the
              rectangle [x0,x1] x [y0,y1] x [z0,z1].
       \param
    */
    // Brick();

    Brick(int order, dim_t n0, dim_t n1, dim_t n2, double x0, double y0, double z0,
      double x1, double y1, double z1, int d0, int d1, int d2,
      int periodic0, int periodic1, int periodic2);

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
       \brief
       interpolates data given on source onto target where source and target
       are given on different domains
    */
    virtual void interpolateAcross(escript::Data& target,
                                   const escript::Data& source) const;

    /**
       \brief
       determines whether interpolation from source to target is possible
    */
    virtual bool probeInterpolationAcross(int, const escript::AbstractDomain&,
            int) const;

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
       returns a pointer to the pXest
    */
    virtual p8est_t* borrow_p8est() { return p8est; };

    /**
       \brief
       refines the mesh
       \param maxRecursion Max levels of recursion
       \param algorithmname The algorithm to use
    */
    virtual void refineMesh(int maxRecursion, std::string algorithmname);

    /**
       \brief
       sets the number of levels of refinement
    */
    virtual void setRefinementLevels(int refinementlevels)
    {
        forestData->max_levels_refinement = refinementlevels;
    };

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

    /**
       \brief
       returns a Data object containing the coordinate information
    */
    int getNumVertices() const { return connectivity->num_vertices;};

    virtual dim_t findNode(const double *coords) const;

    virtual void nodesToDOF(escript::Data& out, const escript::Data& in) const;
    
    virtual dim_t getDofOfNode(dim_t node) const;

    // These functions are used internally
    p8est_t * borrow_p4est() const { return p8est;};

    p8estData * borrow_forestData() { return forestData;};

    p8est_connectivity_t * borrow_connectivity() const { return connectivity; };

    void * borrow_temp_data() { return temp_data; };

    void set_temp_data(addSurfaceData * x) { temp_data = x; };

    void clear_temp_data() { free(temp_data); };

    void print_debug_report(std::string);

////////////////////////////////
protected:

#ifdef ESYS_HAVE_TRILINOS
    virtual esys_trilinos::const_TrilinosGraph_ptr getTrilinosGraph() const;
#endif
#ifdef ESYS_HAVE_PASO
    virtual paso::SystemMatrixPattern_ptr getPasoMatrixPattern(bool reducedRowOrder, bool reducedColOrder) const;
#endif

    virtual dim_t getNumNodes() const;
    virtual dim_t getNumElements() const;
    virtual dim_t getNumFaceElements() const;
    virtual dim_t getNumDOF() const;
    virtual index_t getFirstInDim(unsigned axis) const;

    virtual void assembleCoordinates(escript::Data& arg) const;
    virtual std::vector<IndexVector> getConnections(bool includeShared=false) const;

#ifdef ESYS_HAVE_PASO
    // the Paso System Matrix pattern
    mutable paso::SystemMatrixPattern_ptr m_pattern;
#endif

#ifdef ESYS_HAVE_TRILINOS
    /// Trilinos graph structure, cached for efficiency
    mutable esys_trilinos::const_TrilinosGraph_ptr m_graph;
#endif
    
    virtual void interpolateNodesOnElements(escript::Data& out,
                                  const escript::Data& in, bool reduced) const;
    virtual void interpolateNodesOnFaces(escript::Data& out,
                                         const escript::Data& in,
                                         bool reduced) const;

    template <typename S>
    void interpolateNodesOnElementsWorker(escript::Data& out,
                                  const escript::Data& in, bool reduced, S sentinel) const;
    template <typename S>         
    void interpolateNodesOnFacesWorker(escript::Data& out,
                                         const escript::Data& in,
                                         bool reduced, S sentinel) const;

    virtual void assembleGradient(escript::Data& out,
                                  const escript::Data& in) const;

    template<typename Scalar>
    void assembleGradientImpl(escript::Data& out,
                              const escript::Data& in) const;

////////////////////////////////
private:

    // A p8est
    p8est_t * p8est;

    // The data structure in p8est
    p8estData * forestData;

    // This object records the connectivity of the p8est quadrants
    p8est_connectivity_t * connectivity;

    // This structure records the node numbering information
    p8est_lnodes * nodes;

    // This ghost is needed to initialise the node numbering structure p4est_lnodes
    p8est_ghost_t * nodes_ghost;

    // Pointer that records the location of a temporary data structure
    void * temp_data;

};

typedef POINTER_WRAPPER_CLASS(Brick) OxleyDomainBrick_ptr;

} //end namespace


#endif //__OXLEY_BRICK_H__
