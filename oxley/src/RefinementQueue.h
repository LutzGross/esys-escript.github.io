
#ifndef _OXLEY_REFINEMENTQUEUE
#define _OXLEY_REFINEMENTQUEUE

#include <iostream>
#include <string>
#include <vector>

#include <escript/AbstractContinuousDomain.h>
#include <escript/DataTypes.h>
#include <escript/Pointers.h>

#include <oxley/RefinementType.h>


namespace oxley {

class RefinementQueue;
class RefinementQueue2D;
class RefinementQueue3D;

typedef POINTER_WRAPPER_CLASS(RefinementQueue)   RefinementQueue_Ptr;
typedef POINTER_WRAPPER_CLASS(RefinementQueue2D) RefinementQueue2D_Ptr;
typedef POINTER_WRAPPER_CLASS(RefinementQueue3D) RefinementQueue3D_Ptr;
typedef POINTER_WRAPPER_CLASS(const RefinementQueue)   const_RefinementQueue_Ptr;
typedef POINTER_WRAPPER_CLASS(const RefinementQueue2D) const_RefinementQueue2D_Ptr;
typedef POINTER_WRAPPER_CLASS(const RefinementQueue3D) const_RefinementQueue3D_Ptr;

/**
    \brief
    Abstract class RefinementQueue.

    A queue collects refinement operations and applies them to an oxley
    domain, returning a NEW domain: the source is left untouched. That is why
    the operations live here rather than on the domain, where they used to
    mutate it in place - a domain and the Data defined on it went out of step
    the moment anything refined, with nothing in the type system to say so.

    The refinement a domain is BORN with is separate, and stays a constructor
    argument (Rectangle/Brick take refine_level, an int or a per-block list).
*/
class RefinementQueue
{
public:

	/**
       \brief
       Constructor
    */
	RefinementQueue();

	/**
       \brief
       Destructor
    */
	~RefinementQueue();

	/**
       \brief
       Add to queue
    */
	virtual void addToQueue(RefinementType R);

    /**
       \brief
       Returns the length of the queue
    */
    virtual int getNumberOfOperations()
    {
        return queue.size();
    };

    /**
       \brief
       Returns the n^th refinement
    */
    virtual RefinementType getRefinement(int n);

	/**
       \brief
       A queue of refinements 
    */
	std::vector<RefinementType> queue;

   /**
       \brief
       Returns the length of the queue
    */
    virtual void setRefinementLevel(int n);

    int refinement_levels;

    /**
       \brief
       Prints the current queue to console
    */
    virtual void print();

    /**
       \brief
       Removes the nt^th item from the queue
    */
    virtual void deleteFromQueue(int n);

private:

protected:
};


class RefinementQueue2D : public RefinementQueue
{
public:
   RefinementQueue2D();
   ~RefinementQueue2D();
	/**
       \brief
       RefinementAlgorithms
    */
    /// refines every element of the domain, the old refineMesh("uniform")
    void refineUniform(int level);
    void refinePoint(float x0, float y0, int level);
    void refineRegion(float x0, float y0, float x1, float y1, int level);
    void refineCircle(float x0, float y0, float r, int level);
    void refineBorder(Border b, float dx, int level);
    /// as above, naming the border: top, bottom, left, right (or north,
    /// south, west, east). The Border enum is not exposed to python, so this
    /// is the form callers actually use.
    void refineBorder(std::string border, float dx, int level);
    void refineMask(escript::Data d, int level);

    /**
       \brief
       Applies the queued refinements to a domain and returns the result.

       The domain passed in is NOT modified: everything is done on a copy, so
       the caller keeps a usable handle on the coarser mesh and on any Data
       defined over it.

       \param domain the domain to refine; must be 2D (a Rectangle)
       \return the refined domain
    */
    escript::Domain_ptr apply(escript::Domain_ptr domain);

    /**
       \brief
       Prints the current queue to console
    */
    void print();

    /**
       \brief
       Removes the nt^th item from the queue
    */
    void deleteFromQueue(int n);

};

class RefinementQueue3D : public RefinementQueue
{
public:
   RefinementQueue3D();
   ~RefinementQueue3D();
	/**
       \brief
       RefinementAlgorithms
    */
    /// refines every element of the domain, the old refineMesh("uniform")
    void refineUniform(int level);
    void refinePoint(float x0, float y0, float z0, int level);
    void refineRegion(float x0, float y0, float z0, float x1, float y1, float z1, int level);
    void refineSphere(float x0, float y0, float z0, float r, int level);
    void refineBorder(Border b, float dx, int level);
    /// as above, naming the border: top, bottom, left, right, front, back.
    /// The Border enum is not exposed to python, so this is the form
    /// callers actually use.
    void refineBorder(std::string border, float dx, int level);
    void refineMask(escript::Data d, int level);

    /**
       \brief
       Applies the queued refinements to a domain and returns the result.

       The domain passed in is NOT modified: everything is done on a copy, so
       the caller keeps a usable handle on the coarser mesh and on any Data
       defined over it.

       \param domain the domain to refine; must be 3D (a Brick)
       \return the refined domain
    */
    escript::Domain_ptr apply(escript::Domain_ptr domain);

    /**
       \brief
       Prints the current queue to console
    */
    void print();

    /**
       \brief
       Removes the nt^th item from the queue
    */
    void deleteFromQueue(int n);
};

} //namespace oxley


#endif //_OXLEY_REFINEMENTQUEUE
