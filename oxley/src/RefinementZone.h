
#include <vector>

#include <oxley/RefinementType.h>


/**
    \brief
    Abstract class RefinementZone 
*/
class RefinementZone
{
public:

	/**
       \brief
       Add to queue
    */
	void AddToQueue(RefinementType R);

private:

	/**
       \brief
       A queue of refinements 
    */
	std::vector<RefinementType> queue;

protected:
};


class RefinementZone2D : public RefinementZone
{
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
	/**
       \brief
       RefinementAlgorithms
    */
    void refinePoint(double x0, double y0, double z0);
    void refineRegion(double x0, double y0, double z0, double x1, double y1, double z1);
    void refineSphere(double x0, double y0, double z0, double r);
    void refineBorder(Border b, double dx);
};