
#include <oxley/RefinementZone.h>
#include <oxley/OxleyException.h>

namespace oxley {

void RefinementZone::AddToQueue(RefinementType R)
{	
	queue.push_back(R);
};

RefinementZone::RefinementZone()
{

};

RefinementZone::~RefinementZone()
{

};

RefinementType RefinementZone::getRefinement(int n)
{
    if(n <= getNumberOfOperations())
        return queue[n];
    else
        throw OxleyException("Number is greater than queue length");
};

void RefinementZone2D::refinePoint(double x0, double y0)
{
	Point2DRefinement refine(x0,y0);
	Point2DRefinement * pRefine = &refine;
	RefinementType * alg = new RefinementType;
	alg = static_cast<RefinementType*>(pRefine);
	AddToQueue(*alg);
};

void RefinementZone2D::refineRegion(double x0, double y0, double x1, double y1)
{
	Region2DRefinement refine(x0,y0,x1,y1);
	Region2DRefinement * pRefine = &refine;
	RefinementType * alg = new RefinementType;
	alg = static_cast<RefinementType*>(pRefine);
	AddToQueue(*alg);
};

void RefinementZone2D::refineCircle(double x0, double y0, double r)
{
	CircleRefinement refine(x0,y0,r);
	CircleRefinement * pRefine = &refine;
	RefinementType * alg = new RefinementType;
	alg = static_cast<RefinementType*>(pRefine);
	AddToQueue(*alg);
};

void RefinementZone2D::refineBorder(Border b, double dx)
{
	Border2DRefinement refine(b,dx);
	Border2DRefinement * pRefine = &refine;
	RefinementType * alg = new RefinementType;
	alg = static_cast<RefinementType*>(pRefine);
	AddToQueue(*alg);
};

void RefinementZone3D::refinePoint(double x0, double y0, double z0)
{
	Point3DRefinement refine(x0,y0,z0);
	Point3DRefinement * pRefine = &refine;
	RefinementType * alg = new RefinementType;
	alg = static_cast<RefinementType*>(pRefine);
	AddToQueue(*alg);
};

void RefinementZone3D::refineRegion(double x0, double y0, double z0, double x1, double y1, double z1)
{
	Region3DRefinement refine(x0,y0,z0,x1,y1,z1);
	Region3DRefinement * pRefine = &refine;
	RefinementType * alg = new RefinementType;
	alg = static_cast<RefinementType*>(pRefine);
	AddToQueue(*alg);
};

void RefinementZone3D::refineSphere(double x0, double y0, double z0, double r)
{
	SphereRefinement refine(x0,y0,z0,r);
	SphereRefinement * pRefine = &refine;
	RefinementType * alg = new RefinementType;
	alg = static_cast<RefinementType*>(pRefine);
	AddToQueue(*alg);
};

void RefinementZone3D::refineBorder(Border b, double dx)
{
	Border3DRefinement refine(b,dx);
	Border3DRefinement * pRefine = &refine;
	RefinementType * alg = new RefinementType;
	alg = static_cast<RefinementType*>(pRefine);
	AddToQueue(*alg);
};



} //namespace oxley