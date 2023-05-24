
#include "oxley/RefinementZone.h"
#include "oxley/OxleyException.h"

namespace oxley {

RefinementZone::RefinementZone()
{
	refinement_levels=0;
};

RefinementZone::~RefinementZone()
{

};

void RefinementZone::addToQueue(RefinementType R)
{	
	R.levels=refinement_levels;
	queue.push_back(R);
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
	RefinementType refine;
	refine.Point2DRefinement(x0,y0,refinement_levels);
	addToQueue(refine);
};

void RefinementZone2D::refineRegion(float x0, float y0, float x1, float y1)
{
	RefinementType refine;
	refine.Region2DRefinement(x0,y0,x1,y1,refinement_levels);
	addToQueue(refine);
};

void RefinementZone2D::refineCircle(float x0, float y0, float r)
{
	RefinementType refine;
	refine.CircleRefinement(x0,y0,r,refinement_levels);
	addToQueue(refine);
};

void RefinementZone2D::refineBorder(Border b, float dx)
{
	RefinementType refine;
	refine.Border2DRefinement(b,dx,refinement_levels);
	addToQueue(refine);
};

RefinementZone3D::RefinementZone3D()
{

};

RefinementZone3D::~RefinementZone3D()
{

};

void RefinementZone3D::refinePoint(float x0, float y0, float z0)
{
	RefinementType refine;
	refine.Point3DRefinement(x0,y0,z0,refinement_levels);
	addToQueue(refine);
};

void RefinementZone3D::refineRegion(float x0, float y0, float z0, float x1, float y1, float z1)
{
	RefinementType refine;
	refine.Region3DRefinement(x0,y0,z0,x1,y1,z1,refinement_levels);
	addToQueue(refine);
};

void RefinementZone3D::refineSphere(float x0, float y0, float z0, float r)
{
	RefinementType refine;
	refine.SphereRefinement(x0,y0,z0,r,refinement_levels);
	addToQueue(refine);
};

void RefinementZone3D::refineBorder(Border b, float dx)
{
	RefinementType refine;
	refine.Border3DRefinement(b,dx,refinement_levels);
	addToQueue(refine);
};



} //namespace oxley