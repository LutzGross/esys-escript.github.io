
#ifndef _OXLEY_REFINEMENTZONE
#define _OXLEY_REFINEMENTZONE

#include <vector>

#include <escript/Pointers.h>

#include <oxley/RefinementType.h>


namespace oxley {

class RefinementZone;
class RefinementZone2D;
class RefinementZone3D;

typedef boost::shared_ptr<RefinementZone>   RefinementZone_Ptr;
typedef boost::shared_ptr<RefinementZone2D> RefinementZone2D_Ptr;
typedef boost::shared_ptr<RefinementZone3D> RefinementZone3D_Ptr;

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
	virtual void AddToQueue(RefinementType R);

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

private:

protected:
};


class RefinementZone2D : public RefinementZone
{
public:
	/**
       \brief
       RefinementAlgorithms
    */
    void refinePoint(double x0, double y0);
    void refineRegion(double x0, double y0, double x1, double y1);
    void refineCircle(double x0, double y0, double r);
    void refineBorder(Border b, double dx);
};

class RefinementZone3D : public RefinementZone
{
public:
	/**
       \brief
       RefinementAlgorithms
    */
    void refinePoint(double x0, double y0, double z0);
    void refineRegion(double x0, double y0, double z0, double x1, double y1, double z1);
    void refineSphere(double x0, double y0, double z0, double r);
    void refineBorder(Border b, double dx);
};

} //namespace oxley


#endif //_OXLEY_REFINEMENTZONE
