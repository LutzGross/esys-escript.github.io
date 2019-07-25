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



#include <oxley/Oxley.h>
#ifndef __OXLEY_EXCEPTION_H__
#include <oxley/OxleyException.h>

#include <escript/EsysMPI.h>
#include <escript/AbstractContinuousDomain.h>

// #include <p4est/sc_mpi.h>


namespace oxley {

/* 
This class is the parent of Oxley Rectangle and Brick
*/
class OxleyDomain : public escript::AbstractContinuousDomain
{
public:
    /**
       \brief
       Constructor 
    */
    OxleyDomain(dim_t dim, int order);

    /**
       \brief
       Destructor
    */
    ~OxleyDomain();

    /**
       \brief
       returns a description for this domain
    */
    virtual std::string getDescription() const;
  

    /**
     \brief
     returns a reference to the MPI information wrapper for this domain
    */
    virtual escript::JMPI getMPI() const { return m_mpiInfo; }

    /**
       \brief
       returns the number of processors used for this domain
    */
    virtual int getMPISize() const { return m_mpiInfo->size; }

    /**
       \brief
       returns the MPI rank of this processor
    */
    virtual int getMPIRank() const { return m_mpiInfo->rank; } 

    /**
       \brief
       if compiled for MPI then executes an MPI_Barrier, else does nothing
    */
    virtual void MPIBarrier() const {
#ifdef ESYS_MPI
        MPI_Barrier(m_mpiInfo->comm);
#endif
    }

    /**
       \brief
       returns true if on MPI processor 0, else false
    */
    virtual bool onMasterProcessor() const { return getMPIRank()==0; }

    /**
       \brief
       returns the MPI communicator
    */
    MPI_Comm getMPIComm() const { return m_mpiInfo->comm; }

    /**
       \brief
       returns a description for the given function space type code
    */
    virtual std::string functionSpaceTypeAsString(int fsType) const;

    /**
       \brief
       returns the number of spatial dimensions of the domain
    */
    virtual int getDim() const { return m_numDim; }

    /**
       \brief equality operator
    */
    virtual bool operator==(const escript::AbstractDomain& other) const;

    /**
       \brief inequality operator
    */
    virtual bool operator!=(const escript::AbstractDomain& other) const {
        return !(operator==(other));
    }

    /**
       \brief
       writes the current mesh to a file with the given name
       \param filename The name of the file to write to
    */
    virtual void write(const std::string& filename) const = 0;

    /**
       \brief
       dumps the mesh to a file with the given name
       \param filename The name of the output file
    */
    virtual void dump(const std::string& filename) const = 0;


    /**
       \brief
       returns the tag key for the given sample number
       \param fsType The function space type
       \param sampleNo The sample number
    */
    int getTagFromSampleNo(int fsType, dim_t sampleNo) const;

    /**
       \brief
       returns true if this rank owns the sample id on given function space
    */
    virtual bool ownSample(int fsType, index_t id) const = 0;

    /**
       \brief
       sets a map from a clear tag name to a tag key
       \param name tag name
       \param tag tag key
    */
    virtual void setTagMap(const std::string& name, int tag) {
        m_tagMap[name] = tag;
    }

    /**
       \brief
       returns the tag key for tag name
       \param name tag name
    */
    virtual int getTag(const std::string& name) const {
        if (m_tagMap.find(name) != m_tagMap.end()) {
            return m_tagMap.find(name)->second;
        } else {
            throw OxleyException("getTag: invalid tag name");
        }
    }

    /**
       \brief
       assigns new tag newTag to all samples of given function space with a
       positive value of mask for any of its sample points
    */
    virtual void setTags(int fsType, int newTag, const escript::Data& mask) const;

    /**
       \brief
       returns all tag names in a single string separated by commas
    */
    virtual std::string showTagNames() const;


    /**
       \brief
       returns the array of reference numbers for a function space type
       \param fsType The function space type
    */
    const index_t* borrowSampleReferenceIDs(int fsType) const = 0;

    /**
       \brief
       copies the surface normals at data points into out. The actual function
       space to be considered is defined by out. out has to be defined on this
       domain.
    */
    virtual void setToNormal(escript::Data& out) const = 0;

    /**
       \brief
       copies the size of samples into out. The actual function space to be
       considered is defined by out. out has to be defined on this domain.
    */
    virtual void setToSize(escript::Data& out) const = 0;

    /**
       \brief
       interpolates data given on source onto target where source and target
       have to be given on the same domain
    */
    virtual void interpolateOnDomain(escript::Data& target,
            const escript::Data& source) const;

    /**
       \brief
       returns true if data on fsType_source can be interpolated onto
       fsType_target, false otherwise
    */
    virtual bool probeInterpolationOnDomain(int fsType_source,
            int fsType_target) const;

    /**
       \brief Preferred direction of interpolation.
       If you really need to test for a particular direction, then use
       probeInterpolation.
       \return 0 for not possible, 1 for possible and preferred, -1 other
               direction preferred (does not mean this direction is possible)
    */
    virtual signed char preferredInterpolationOnDomain(int fsType_source,
                                                       int fsType_target) const;

    /**
       \brief
       given a vector of FunctionSpace type codes, passes back a code which all
       can be interpolated to
       \return true if result is valid, false if not
    */
    bool
    commonFunctionSpace(const std::vector<int>& fs, int& resultcode) const;

    /**
       \brief
       interpolates data given on source onto target where source and target
       are given on different domains
    */
    virtual void interpolateAcross(escript::Data& target,
                                   const escript::Data& source) const = 0;

    /**
       \brief
       determines whether interpolation from source to target is possible
    */
    virtual bool probeInterpolationAcross(int, const escript::AbstractDomain&,
            int) const = 0;

    /**
       \brief
       returns locations in the SEM nodes
    */
    virtual escript::Data getX() const;

    /**
       \brief
       returns boundary normals at the quadrature point on the face elements
    */
    virtual escript::Data getNormal() const;

    /**
       \brief returns the element size
    */
    virtual escript::Data getSize() const;

    /**
       \brief
       copies the location of data points into arg. The domain of arg has to
       match this domain.
    */
    virtual void setToX(escript::Data& arg) const;


    /**
       \brief
       copies the gradient of 'in' into 'out'. The actual function space to be
       considered for the gradient is defined by 'in'. Both arguments have to
       be defined on this domain.
    */
    virtual void setToGradient(escript::Data& out,
            const escript::Data& in) const;

    /**
       \brief
       returns true if data on this domain and given function space type has
       to be considered as cell centered data
    */
    virtual bool isCellOriented(int fsType) const;

    /**
       \brief
       returns the number of tags in use for a function space type
    */
    virtual int getNumberOfTagsInUse(int fsType) const;

    /**
       \brief
       returns a pointer to the list of tags in use for a function space type
    */
    virtual const int* borrowListOfTagsInUse(int fsType) const;

    /**
       \brief
       checks if this domain allows tags for the specified function space type
    */
    virtual bool canTag(int fsType) const;

    /**
       \brief
       returns the approximation order used for a function space
    */
    virtual int getApproximationOrder(int fsType) const { return 1; }

    /**
       \brief
       returns true if this domain supports contact elements, false otherwise
    */
    virtual bool supportsContactElements() const { return false; }

    /**
       \brief
       writes the mesh to file
    */
    virtual void writeToVTK(std::string filename) const;

protected:
    int m_numDim;
    escript::JMPI m_mpiInfo;
    TagMap m_tagMap;

};

} // end of namespace oxley


#endif //__OXLEY_EXCEPTION_H__