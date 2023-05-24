
#ifndef _OXLEY_REFINEMENTZONE
#define _OXLEY_REFINEMENTZONE

#include <iostream>
#include <vector>

#include <escript/Pointers.h>

#include <oxley/RefinementType.h>


namespace oxley {

class RefinementZone;
class RefinementZone2D;
class RefinementZone3D;

typedef POINTER_WRAPPER_CLASS(RefinementZone)   RefinementZone_Ptr;
typedef POINTER_WRAPPER_CLASS(RefinementZone2D) RefinementZone2D_Ptr;
typedef POINTER_WRAPPER_CLASS(RefinementZone3D) RefinementZone3D_Ptr;
typedef POINTER_WRAPPER_CLASS(const RefinementZone)   const_RefinementZone_Ptr;
typedef POINTER_WRAPPER_CLASS(const RefinementZone2D) const_RefinementZone2D_Ptr;
typedef POINTER_WRAPPER_CLASS(const RefinementZone3D) const_RefinementZone3D_Ptr;

/**
    \brief
    Abstract class RefinementZone 
*/
class RefinementZone
{
public:

	/**
       \brief
       Constructor
    */
	RefinementZone();

	/**
       \brief
       Destructor
    */
	~RefinementZone();

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

private:

protected:
};


class RefinementZone2D : public RefinementZone
{
public:
   RefinementZone2D();
   ~RefinementZone2D();
	/**
       \brief
       RefinementAlgorithms
    */
    void refinePoint(float x0, float y0);
    void refineRegion(float x0, float y0, float x1, float y1);
    void refineCircle(float x0, float y0, float r);
    void refineBorder(Border b, float dx);
};

class RefinementZone3D : public RefinementZone
{
public:
   RefinementZone3D();
   ~RefinementZone3D();
	/**
       \brief
       RefinementAlgorithms
    */
    void refinePoint(float x0, float y0, float z0);
    void refineRegion(float x0, float y0, float z0, float x1, float y1, float z1);
    void refineSphere(float x0, float y0, float z0, float r);
    void refineBorder(Border b, float dx);
};

} //namespace oxley


#endif //_OXLEY_REFINEMENTZONE
