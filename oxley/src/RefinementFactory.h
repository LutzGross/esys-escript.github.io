
#ifndef _OXLEY_REFINEMENTFACTORY
#define _OXLEY_REFINEMENTFACTORY

#include <iostream>
#include <string>
#include <vector>

#include <escript/AbstractContinuousDomain.h>
#include <escript/DataTypes.h>
#include <escript/Pointers.h>

#include <oxley/RefinementType.h>


namespace oxley {

class RefinementFactory;
class RefinementFactory2D;
class RefinementFactory3D;

typedef POINTER_WRAPPER_CLASS(RefinementFactory)   RefinementFactory_Ptr;
typedef POINTER_WRAPPER_CLASS(RefinementFactory2D) RefinementFactory2D_Ptr;
typedef POINTER_WRAPPER_CLASS(RefinementFactory3D) RefinementFactory3D_Ptr;
typedef POINTER_WRAPPER_CLASS(const RefinementFactory)   const_RefinementFactory_Ptr;
typedef POINTER_WRAPPER_CLASS(const RefinementFactory2D) const_RefinementFactory2D_Ptr;
typedef POINTER_WRAPPER_CLASS(const RefinementFactory3D) const_RefinementFactory3D_Ptr;

/**
    \brief
    Abstract class RefinementFactory.

    A factory collects refinement operations and applies them to an oxley
    domain, returning a NEW domain: the source is left untouched. That is why
    the operations live here rather than on the domain, where they used to
    mutate it in place - a domain and the Data defined on it went out of step
    the moment anything refined, with nothing in the type system to say so.

    The refinement a domain is BORN with is separate, and stays a constructor
    argument (Rectangle/Brick take refine_level, an int or a per-block list).
*/
class RefinementFactory
{
public:

	/**
       \brief
       Constructor
    */
	RefinementFactory();

	/**
       \brief
       Destructor
    */
	~RefinementFactory();

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


class RefinementFactory2D : public RefinementFactory
{
public:
   RefinementFactory2D();
   ~RefinementFactory2D();
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

class RefinementFactory3D : public RefinementFactory
{
public:
   RefinementFactory3D();
   ~RefinementFactory3D();
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


#endif //_OXLEY_REFINEMENTFACTORY
