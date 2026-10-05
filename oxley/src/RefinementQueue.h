
#ifndef _OXLEY_REFINEMENTQUEUE
#define _OXLEY_REFINEMENTQUEUE

#include <array>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include <escript/AbstractContinuousDomain.h>
#include <escript/Data.h>
#include <escript/DataTypes.h>
#include <escript/Pointers.h>

#include <oxley/RefinementType.h>


namespace oxley {

class RefinementQueue;
class RefinementQueue2D;
class RefinementQueue3D;
// the domains a queue grows a new forest of. Only referred to here; including
// their headers walks into the OxleyData.h <-> Rectangle.h include cycle, which
// is also why each apply() is defined in the domain's own translation unit.
class Rectangle;
class Brick;

typedef POINTER_WRAPPER_CLASS(RefinementQueue)   RefinementQueue_Ptr;
typedef POINTER_WRAPPER_CLASS(RefinementQueue2D) RefinementQueue2D_Ptr;
typedef POINTER_WRAPPER_CLASS(RefinementQueue3D) RefinementQueue3D_Ptr;
typedef POINTER_WRAPPER_CLASS(const RefinementQueue)   const_RefinementQueue_Ptr;
typedef POINTER_WRAPPER_CLASS(const RefinementQueue2D) const_RefinementQueue2D_Ptr;
typedef POINTER_WRAPPER_CLASS(const RefinementQueue3D) const_RefinementQueue3D_Ptr;

/// the masks handed to apply(), by the tag refineMask queued them under
typedef std::map<std::string, escript::Data> MaskMap;

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
    /// throws unless masks holds exactly the tags refineMask queued; who
    /// names the caller in the message
    void checkMasks(const MaskMap& masks, const std::string& who);
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
    /// queues a refinement where the mask named tag is positive. The queue is
    /// a template that is not tied to a domain, so it holds only the name: the
    /// mask itself is handed to apply().
    void refineMask(std::string tag, int level);

    /**
       \brief
       Applies the queued refinements to a domain and returns the result.

       The domain passed in is NOT modified: everything is done on a copy, so
       the caller keeps a usable handle on the coarser mesh and on any Data
       defined over it.

       Collective: every rank of the domain's communicator must call it, with
       the same queue.

       \param domain the domain to refine; must be 2D (a Rectangle)
       \param masks the mask for each tag refineMask queued, a scalar Data
              on domain. Every tag must have one, and every mask a tag.
       \return the refined domain
    */
    escript::Domain_ptr apply(escript::Domain_ptr domain,
                              const MaskMap& masks = MaskMap());

    /**
       \brief
       Prints the current queue to console
    */
    void print();

private:
    // The forest-growing steps, one per queued operation. These used to be
    // refine* methods on Rectangle; they moved here because a domain must not
    // refine itself - every Data object over it would then refer to a mesh that
    // is gone. They act on the NEW domain apply() has just built and nobody
    // else holds yet, which is the only forest it is safe to change.
    //
    // Defined in Rectangle.cpp, next to apply(), for the same include-cycle
    // reason apply() is.
    void growByAlgorithm(Rectangle& dom, std::string algorithmname);
    void growAtBorder(Rectangle& dom, std::string boundaryname, double dx);
    void growInRegion(Rectangle& dom, double x0, double x1, double y0, double y1);
    void growAtPoint(Rectangle& dom, double x0, double y0);
    void growAtCircle(Rectangle& dom, double x0, double y0, double r);
    void growFromMask(Rectangle& dom, const std::set<std::array<long,4>>& masked);
    /// the leaves of source where mask is positive, gathered from all ranks
    void collectMasked(Rectangle& source, const escript::Data& mask,
                       std::set<std::array<long,4>>& masked);

public:

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
    /// see RefinementQueue2D::refineMask
    void refineMask(std::string tag, int level);

    /**
       \brief
       Applies the queued refinements to a domain and returns the result.

       The domain passed in is NOT modified: everything is done on a copy, so
       the caller keeps a usable handle on the coarser mesh and on any Data
       defined over it.

       Collective: every rank of the domain's communicator must call it, with
       the same queue.

       \param domain the domain to refine; must be 3D (a Brick)
       \param masks the mask for each tag refineMask queued, a scalar Data
              on domain. Every tag must have one, and every mask a tag.
       \return the refined domain
    */
    escript::Domain_ptr apply(escript::Domain_ptr domain,
                              const MaskMap& masks = MaskMap());

    /**
       \brief
       Prints the current queue to console
    */
    void print();

private:
    // The forest-growing steps, one per queued operation. These used to be
    // refine* methods on Brick; they moved here because a domain must not
    // refine itself - every Data object over it would then refer to a mesh that
    // is gone. They act on the NEW domain apply() has just built and nobody
    // else holds yet, which is the only forest it is safe to change.
    //
    // Defined in Brick.cpp, next to apply(), for the same include-cycle reason
    // apply() is.
    void growByAlgorithm(Brick& dom, std::string algorithmname);
    void growAtBorder(Brick& dom, std::string boundaryname, double dx);
    void growInRegion(Brick& dom, double x0, double x1, double y0,
                      double y1, double z0, double z1);
    void growAtPoint(Brick& dom, double x0, double y0, double z0);
    void growAtSphere(Brick& dom, double x0, double y0, double z0, double r);
    void growFromMask(Brick& dom, const std::set<std::array<long,5>>& masked);
    /// the leaves of source where mask is positive, gathered from all ranks
    void collectMasked(Brick& source, const escript::Data& mask,
                       std::set<std::array<long,5>>& masked);

public:

    /**
       \brief
       Removes the nt^th item from the queue
    */
    void deleteFromQueue(int n);
};

} //namespace oxley


#endif //_OXLEY_REFINEMENTQUEUE
