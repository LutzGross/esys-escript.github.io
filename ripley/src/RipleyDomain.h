
/*******************************************************
*
* Copyright (c) 2003-2011 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/

#ifndef __RIPLEY_DOMAIN_H__
#define __RIPLEY_DOMAIN_H__

#include <ripley/NodeFile.h>
#include <ripley/ElementFile.h>
#include <ripley/RipleyError.h>
#include <ripley/SystemMatrixAdapter.h>
#include <ripley/TransportProblemAdapter.h>

extern "C" {
//#include "ripley/Assemble.h"
#include "paso/SystemMatrixPattern.h"
#include "paso/Transport.h"
}
#include <escript/AbstractContinuousDomain.h>
#include <escript/FunctionSpace.h>
#include <escript/FunctionSpaceFactory.h>

#include <boost/python/dict.hpp>
#include <boost/python/extract.hpp>

namespace ripley {

/**
   \brief
   RipleyDomain implements the AbstractContinuousDomain interface
   for the Ripley library.

   Description:
   RipleyDomain implements the AbstractContinuousDomain interface
   for the Ripley library.
*/

class RipleyDomain : public escript::AbstractContinuousDomain
{
public:

    /**
       \brief
       Constructor with name, dimensionality and mpi info structure.
    */
    RIPLEY_DLL_API
    RipleyDomain(const std::string& name, dim_t numDim, Esys_MPIInfo *mpiInfo);

    /**
       \brief
       Destructor for RipleyDomain. As specified in the constructor
       this calls Ripley_Mesh_free for the pointer given to the
       constructor.
    */
    RIPLEY_DLL_API
    ~RipleyDomain();

    /**
       \brief
       returns the number of processors used for this domain
    */
    RIPLEY_DLL_API
    virtual int getMPISize() const;

    /**
       \brief
       returns the MPI rank of this processor
    */
    RIPLEY_DLL_API
    virtual int getMPIRank() const;

    /**
       \brief
       if compiled for MPI then executes an MPI_Barrier, else does nothing
    */
    RIPLEY_DLL_API
    virtual void MPIBarrier() const;

    /**
       \brief
       returns true if on MPI processor 0, else false
    */
    RIPLEY_DLL_API
    virtual bool onMasterProcessor() const;

    /**
       \brief
       returns the MPI communicator
    */
    RIPLEY_DLL_API
  #ifdef ESYS_MPI
    MPI_Comm
  #else
    unsigned int
  #endif
    getMPIComm() const;

    /**
       \brief
       writes the current mesh to a file with the given name.
       \param fileName The name of the file to write to.
    */
    RIPLEY_DLL_API
    void write(const std::string& fileName) const;

    /**
       \brief
       \param full whether to print all data, including coordinates etc.
    */
    RIPLEY_DLL_API
    void Print_Mesh_Info(const bool full=false) const;

    /**
       \brief
       dumps the mesh to a file with the given name.
       \param fileName The name of the output file
    */
    RIPLEY_DLL_API
    void dump(const std::string& fileName) const;

    /**
       \brief
       returns the tag key for the given sample number.
       \param functionSpaceType The function space type.
       \param sampleNo The sample number.
    */
    RIPLEY_DLL_API
    int getTagFromSampleNo(int functionSpaceType, int sampleNo) const;

    /**
       \brief
       returns the reference number of the given sample number.
       \param functionSpaceType The function space type.
    */
    RIPLEY_DLL_API
    const int* borrowSampleReferenceIDs(int functionSpaceType) const;

    /**
       \brief
       returns true if the given integer is a valid function space type
       for this domain.
    */
    RIPLEY_DLL_API
    virtual bool isValidFunctionSpaceType(int functionSpaceType) const;

    /**
       \brief
       returns a description for this domain
    */
    RIPLEY_DLL_API
    virtual std::string getDescription() const;

    /**
       \brief
       returns a description for the given function space type code
    */
    RIPLEY_DLL_API
    virtual std::string functionSpaceTypeAsString(int functionSpaceType) const;

    /**
       \brief
       builds the table of function space type names
    */
    RIPLEY_DLL_API
    void setFunctionSpaceTypeNames();

    /**
       \brief
       returns a continuous FunctionSpace code
    */
    RIPLEY_DLL_API
    virtual int getContinuousFunctionCode() const;

    /**
       \brief
       returns a continuous on reduced order nodes FunctionSpace code
    */
    RIPLEY_DLL_API
    virtual int getReducedContinuousFunctionCode() const;

    /**
       \brief
       returns a function FunctionSpace code
    */
    RIPLEY_DLL_API
    virtual int getFunctionCode() const;

    /**
       \brief
       returns a function with reduced integration order FunctionSpace code
    */
    RIPLEY_DLL_API
    virtual int getReducedFunctionCode() const;

    /**
       \brief
       returns a function on boundary FunctionSpace code
    */
    RIPLEY_DLL_API
    virtual int getFunctionOnBoundaryCode() const;

    /**
       \brief
       returns a function on boundary with reduced integration order
       FunctionSpace code
    */
    RIPLEY_DLL_API
    virtual int getReducedFunctionOnBoundaryCode() const;

    /**
       \brief
       return a FunctionOnContactZero code
    */
    RIPLEY_DLL_API
    virtual int getFunctionOnContactZeroCode() const;

    /**
       \brief
       returns a FunctionOnContactZero code with reduced integration order
    */
    RIPLEY_DLL_API
    virtual int getReducedFunctionOnContactZeroCode() const;

    /**
       \brief
       returns a FunctionOnContactOne code
    */
    RIPLEY_DLL_API
    virtual int getFunctionOnContactOneCode() const;

    /**
       \brief
       returns a FunctionOnContactOne code with reduced integration order
    */
    RIPLEY_DLL_API
    virtual int getReducedFunctionOnContactOneCode() const;

    /**
       \brief
       returns a Solution FunctionSpace code
    */
    RIPLEY_DLL_API
    virtual int getSolutionCode() const;

    /**
       \brief
       returns a ReducedSolution FunctionSpace code
    */
    RIPLEY_DLL_API
    virtual int getReducedSolutionCode() const;

    /**
       \brief
       returns a DiracDeltaFunctions FunctionSpace code
    */
    RIPLEY_DLL_API
    virtual int getDiracDeltaFunctionsCode() const;

    /**
       \brief
    */
    typedef std::map<int, std::string> FunctionSpaceNamesMapType;

    /**
       \brief
    */
    RIPLEY_DLL_API
    virtual int getDim() const;

    /**
       \brief
        Returns a status indicator of the domain. The status identifier should be
        unique over the lifetime of the object but may be updated if changes to
        the domain happen, e.g. modifications to its geometry.
    */
    RIPLEY_DLL_API
    virtual StatusType getStatus() const;

    /**
       \brief
       returns the number of data points summed across all MPI processes
    */
    RIPLEY_DLL_API
    virtual int getNumDataPointsGlobal() const;

    /**
       \brief
       returns the number of data points per sample, and the number of samples
       as a pair.
       \param functionSpaceCode The function space type
    */
    RIPLEY_DLL_API
    virtual std::pair<int,int> getDataShape(int functionSpaceCode) const;

    /**
       \brief
       copies the location of data points into arg. The domain of arg has to
       match this.
    */
    RIPLEY_DLL_API
    virtual void setToX(escript::Data& arg) const;

    /**
       \brief
       sets a map from a clear tag name to a tag key
       \param name tag name
       \param tag tag key
    */
    RIPLEY_DLL_API
    virtual void setTagMap(const std::string& name, int tag);

    /**
       \brief
       returns the tag key for tag name
       \param name tag name
    */
    RIPLEY_DLL_API
    virtual int getTag(const std::string& name) const;

    /**
       \brief
       returns true if name is a defined tag name
       \param name tag name to be checked
    */
    RIPLEY_DLL_API
    virtual bool isValidTagName(const std::string& name) const;

    /**
       \brief
       returns all tag names in a single string separated by commas
    */
    RIPLEY_DLL_API
    virtual std::string showTagNames() const;

    /**
       \brief
       assigns new location to the domain
    */
    RIPLEY_DLL_API
    virtual void setNewX(const escript::Data& arg);

    /**
       \brief
       interpolates data given on source onto target where source and target
       have to be given on the same domain
    */
    RIPLEY_DLL_API
    virtual void interpolateOnDomain(escript::Data& target,const escript::Data& source) const;

    /**
       \brief
    */
    RIPLEY_DLL_API
    virtual bool probeInterpolationOnDomain(int functionSpaceType_source,int functionSpaceType_target) const;

    /**
       \brief
       given a vector of FunctionSpace typecodes, passes back a code which all
       can be interpolated to
       \return true if result is valid, false if not
    */
    RIPLEY_DLL_API
    bool
    commonFunctionSpace(const std::vector<int>& fs, int& resultcode) const;

    /**
       \brief
       interpolates data given on source onto target where source and target
       are given on different domains
    */
    RIPLEY_DLL_API
    virtual void interpolateACross(escript::Data& target, const escript::Data& source) const;

    /**
       \brief
       determines whether interpolation from source to target is possible.
    */
    RIPLEY_DLL_API
    virtual bool probeInterpolationACross(int functionSpaceType_source,const escript::AbstractDomain& targetDomain, int functionSpaceType_target) const;

    /**
       \brief
       copies the surface normals at data points into out. The actual function
       space to be considered is defined by out. out has to be defined on this.
    */
    RIPLEY_DLL_API
    virtual void setToNormal(escript::Data& out) const;

    /**
       \brief
       copies the size of samples into out. The actual function space to be
       considered is defined by out. out has to be defined on this.
    */
    RIPLEY_DLL_API
    virtual void setToSize(escript::Data& out) const;

    /**
       \brief
       copies the gradient of arg into grad. The actual function space to be
       considered for the gradient is defined by grad. arg and grad have to be
       defined on this.
    */
    RIPLEY_DLL_API
    virtual void setToGradient(escript::Data& grad,const escript::Data& arg) const;

    /**
       \brief
       copies the integrals of the function defined by arg into integrals.
       arg has to be defined on this.
    */
    RIPLEY_DLL_API
    virtual void setToIntegrals(std::vector<double>& integrals,const escript::Data& arg) const;

    /**
       \brief
       returns the identifier of the matrix type to be used for the global
       stiffness matrix when a particular solver, package, perconditioner,
       and symmetric matrix is used.
       \param solver
       \param preconditioner
       \param package
       \param symmetry
    */
    RIPLEY_DLL_API
    virtual int getSystemMatrixTypeId(const int solver, const int preconditioner, const int package, const bool symmetry) const;

    /**
       \brief
       returns the identifier of the transport problem type to be used when a
       particular solver, perconditioner, package and symmetric matrix is used.
       \param solver
       \param preconditioner
       \param package
       \param symmetry
    */
    RIPLEY_DLL_API
    virtual int getTransportTypeId(const int solver, const int preconditioner, const int package, const bool symmetry) const;

    /**
       \brief
       returns true if data on this domain and a function space of type
       functionSpaceCode has to be considered as cell centered data.
    */
    RIPLEY_DLL_API
    virtual bool isCellOriented(int functionSpaceCode) const;

    /**
       \brief
    */
    RIPLEY_DLL_API
    virtual bool ownSample(int fs_code, index_t id) const;

    /**
       \brief
       adds a PDE onto the stiffness matrix mat and a rhs
    */
    RIPLEY_DLL_API
    virtual void addPDEToSystem(
                       escript::AbstractSystemMatrix& mat, escript::Data& rhs,
                       const escript::Data& A, const escript::Data& B, const escript::Data& C,
                       const escript::Data& D, const escript::Data& X, const escript::Data& Y,
                       const escript::Data& d, const escript::Data& y,
               const escript::Data& d_contact, const escript::Data& y_contact,
                       const escript::Data& d_dirac, const escript::Data& y_dirac) const;


    /**
       \brief
       adds a PDE onto the lumped stiffness matrix mat
    */
    RIPLEY_DLL_API
    virtual void addPDEToLumpedSystem(
                       escript::Data& mat,
                       const escript::Data& D,
                       const escript::Data& d,
                       const escript::Data& d_dirac,
                       const bool useHRZ) const;

    /**
       \brief
       adds a PDE onto the stiffness matrix mat and a rhs
    */
    RIPLEY_DLL_API
    virtual void addPDEToRHS(escript::Data& rhs,
                       const escript::Data& X, const escript::Data& Y,
                       const escript::Data& y, const escript::Data& y_contact,
                       const escript::Data& y_dirac) const;

    /**
       \brief
       adds a PDE onto a transport problem
    */
    RIPLEY_DLL_API
    virtual void addPDEToTransportProblem(
                       escript::AbstractTransportProblem& tp, escript::Data& source,
                       const escript::Data& M,
                       const escript::Data& A, const escript::Data& B, const escript::Data& C,const  escript::Data& D,
                       const  escript::Data& X,const  escript::Data& Y,
                       const escript::Data& d, const escript::Data& y,
                       const escript::Data& d_contact,const escript::Data& y_contact,
                       const escript::Data& d_dirac,const escript::Data& y_dirac) const;


    /**
       \brief
       creates a SystemMatrixAdapter stiffness matrix and initializes it with zeros
    */
    RIPLEY_DLL_API
    escript::ASM_ptr newSystemMatrix(
                        const int row_blocksize,
                        const escript::FunctionSpace& row_functionspace,
                        const int column_blocksize,
                        const escript::FunctionSpace& column_functionspace,
                        const int type) const;

    /**
     \brief
      creates a TransportProblemAdapter
    */
    RIPLEY_DLL_API
    escript::ATP_ptr newTransportProblem(
                        const bool useBackwardEuler,
                        const int blocksize,
                        const escript::FunctionSpace& functionspace,
                        const int type) const;

    /**
       \brief
       returns locations in the FEM nodes
    */
    RIPLEY_DLL_API
    virtual escript::Data getX() const;

    /**
       \brief
       returns boundary normals at the quadrature point on the face elements
    */
    RIPLEY_DLL_API
    virtual escript::Data getNormal() const;

    /**
       \brief returns the element size
    */
    RIPLEY_DLL_API
    virtual escript::Data getSize() const;

    /**
       \brief equality operator
    */
    RIPLEY_DLL_API
    virtual bool operator==(const escript::AbstractDomain& other) const;

    /**
       \brief unequality operator
    */
    RIPLEY_DLL_API
    virtual bool operator!=(const escript::AbstractDomain& other) const;

    /**
       \brief
       assigns new tag newTag to all samples of functionspace with a positive
       value of mask for any its sample point
    */
    RIPLEY_DLL_API
    virtual void setTags(const int functionSpaceType, const int newTag, const escript::Data& mask) const;

    /**
       \brief
       returns the number of tags in use
    */
    RIPLEY_DLL_API
    virtual int getNumberOfTagsInUse(int functionSpaceCode) const;

    /**
       \brief
       returns a pointer to the list of tags in use
    */
    RIPLEY_DLL_API
    virtual const int* borrowListOfTagsInUse(int functionSpaceCode) const;

    /**
       \brief
       checks if this domain allows tags for the specified functionSpaceCode
    */
    RIPLEY_DLL_API
    virtual
    bool canTag(int functionSpaceCode) const;

    /**
       \brief
       returns the approximation order used for a function space
    */
    RIPLEY_DLL_API
    virtual
    int getApproximationOrder(const int functionSpaceCode) const { return 1; }

    /**
       \brief
       returns true if this domain supports contact elements, false otherwise
    */
    RIPLEY_DLL_API
    bool supportsContactElements() const { return false; }

    /**
       \brief
       prepares the mesh for further calculations
    */
    void prepare(bool optimize);

    /**
       \brief
       creates node/DOF mappings
    */
    void createMappings(const IndexVector &dofDistribution,
                        const IndexVector &nodeDistribution);

    NodeFile_ptr getNodes() const { return m_nodes; }
    ElementFile_ptr getElements() const { return m_elements; }
    ElementFile_ptr getFaceElements() const { return m_faceElements; }
    ElementFile_ptr getPoints() const { return m_points; }

    void setNodes(NodeFile_ptr file) { m_nodes = file; }
    void setElements(ElementFile_ptr file) { m_elements = file; }
    void setFaceElements(ElementFile_ptr file) { m_faceElements = file; }
    void setPoints(ElementFile_ptr file) { m_points = file; }

protected:
    /**
       \brief
       Redistributes elements to minimize communication during assemblage.
    */
    void optimizeElementOrdering(void);

    /**
       \brief
       Optimizes the labeling of the DOFs on each processor.
    */
    void optimizeDOFLabeling(const IndexVector &distribution);

    /**
       \brief
       Resets the tags in the node file and all element files.
    */
    void updateTagsInUse(void);

    /**
       \brief
       Redistributes the Nodes and Elements including overlap according to
       the DOF_distribution. It will create an element coloring but will not
       create any mappings.
    */
    void distributeByRankOfDOF(const IndexVector &dofDistribution);

    /**
       \brief
       Tries to reduce the coloring for all element files.
    */
    void createColoring(const IndexVector &node_localDOF_map);

    /**
       \brief
       At input the element nodes refer to the numbering defined by the global
       Id assigned to the nodes in the NodeFile. It is also not ensured that
       all nodes referred by an element are actually available on the process.
       At the output, a local node labeling is used and all nodes are
       available. In particular the numbering of the element nodes is between
       0 and NodeFile->numNodes. The method does not create a distribution of
       the degrees of freedom.
    */
    void resolveNodeIds(void);

    /**
       \brief
       Assigns new node reference numbers to elements.
       If k is the old node, the new node is newNode[k-offset].
    */
    void relabelElementNodes(const IndexVector &newNode, index_t offset);

    /**
      \brief
      optimizes the distribution of DOFs across processors using ParMETIS.
      On return a new distribution is given and the globalDOF are relabeled
      accordingly but the mesh has not been redistributed yet.
    */
    void optimizeDOFDistribution(RankVector &distribution);

    /**
       \brief
       Marks the used nodes with offset.
    */
    void markNodes(IndexVector &mask, index_t offset);

    /**
       \brief
       Returns a reference to the matrix pattern.
    */
    Paso_SystemMatrixPattern *getPattern(bool reduce_row_order, bool reduce_col_order) const;

    /**
       \brief
       Creates a new matrix pattern.
    */
    Paso_SystemMatrixPattern *makePattern(bool reduce_row_order, bool reduce_col_order) const;

private:
    std::string m_name;
    Esys_MPIInfo *m_mpiInfo;
    NodeFile_ptr m_nodes;
    ElementFile_ptr m_elements;
    ElementFile_ptr m_faceElements;
    ElementFile_ptr m_points;
    TagMap m_tagMap;
    Paso_SystemMatrixPattern *m_fullFullPattern;
    Paso_SystemMatrixPattern *m_fullReducedPattern;
    Paso_SystemMatrixPattern *m_reducedFullPattern;
    Paso_SystemMatrixPattern *m_reducedReducedPattern;

    static FunctionSpaceNamesMapType m_functionSpaceTypeNames;
};

} // end of namespace ripley

#endif // __RIPLEY_DOMAIN_H__

