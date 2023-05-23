
#include <oxley/RefinementZone.h>
#include <oxley/OxleyException.h>

namespace oxley {

void RefinementZone::addToQueue(RefinementType R)
{	
	queue.push_back(R);
	R.levels=refinement_levels;
};

RefinementZone::RefinementZone()
{
	refinement_levels=1;
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

void RefinementZone::setRefinementLevel(int n)
{
	if(n >= 0)
    	refinement_levels=n;
    else
    	throw OxleyException("The levels of refinement must be equal to or greater than zero.");
};

RefinementZone2D::RefinementZone2D()
{

};

RefinementZone2D::~RefinementZone2D()
{

};

void RefinementZone2D::refinePoint(float x0, float y0)
{
	Point2DRefinement refine(x0,y0);
	Point2DRefinement * pRefine = &refine;
	RefinementType * alg = new RefinementType;
	alg = static_cast<RefinementType*>(pRefine);
	addToQueue(*alg);
};

void RefinementZone2D::refineRegion(float x0, float y0, float x1, float y1)
{
	Region2DRefinement refine(x0,y0,x1,y1);
	Region2DRefinement * pRefine = &refine;
	RefinementType * alg = new RefinementType;
	alg = static_cast<RefinementType*>(pRefine);
	addToQueue(*alg);
};

void RefinementZone2D::refineCircle(float x0, float y0, float r)
{
	CircleRefinement refine(x0,y0,r);
	CircleRefinement * pRefine = &refine;
	RefinementType * alg = new RefinementType;
	alg = static_cast<RefinementType*>(pRefine);
	addToQueue(*alg);
};

void RefinementZone2D::refineBorder(Border b, float dx)
{
	Border2DRefinement refine(b,dx);
	Border2DRefinement * pRefine = &refine;
	RefinementType * alg = new RefinementType;
	alg = static_cast<RefinementType*>(pRefine);
	addToQueue(*alg);
};

RefinementZone3D::RefinementZone3D()
{

};

RefinementZone3D::~RefinementZone3D()
{

};

void RefinementZone3D::refinePoint(float x0, float y0, float z0)
{
	Point3DRefinement refine(x0,y0,z0);
	Point3DRefinement * pRefine = &refine;
	RefinementType * alg = new RefinementType;
	alg = static_cast<RefinementType*>(pRefine);
	addToQueue(*alg);
};

void RefinementZone3D::refineRegion(float x0, float y0, float z0, float x1, float y1, float z1)
{
	Region3DRefinement refine(x0,y0,z0,x1,y1,z1);
	Region3DRefinement * pRefine = &refine;
	RefinementType * alg = new RefinementType;
	alg = static_cast<RefinementType*>(pRefine);
	addToQueue(*alg);
};

void RefinementZone3D::refineSphere(float x0, float y0, float z0, float r)
{
	SphereRefinement refine(x0,y0,z0,r);
	SphereRefinement * pRefine = &refine;
	RefinementType * alg = new RefinementType;
	alg = static_cast<RefinementType*>(pRefine);
	addToQueue(*alg);
};

void RefinementZone3D::refineBorder(Border b, float dx)
{
	Border3DRefinement refine(b,dx);
	Border3DRefinement * pRefine = &refine;
	RefinementType * alg = new RefinementType;
	alg = static_cast<RefinementType*>(pRefine);
	addToQueue(*alg);
};



} //namespace oxley