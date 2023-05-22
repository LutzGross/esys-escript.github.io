
#include <oxley/OxleyException.h>
#include <oxley/RefinementType.h>


Point2DRefinement::Point2DRefinement(double x, double y)
{
	flavour=POINT2D;
	x0=x;
	y0=y;
};

Point2DRefinement::~Point2DRefinement()
{

};

Region2DRefinement::Region2DRefinement(double x00, double y00, double x11, double y11)
{
	flavour=REGION2D;
	x0=x00;
	y0=y00;
	x1=x11;
	y1=y11;
};

Region2DRefinement::~Region2DRefinement()
{

};

Point3DRefinement::Point3DRefinement(double x, double y, double z)
{
	flavour=POINT3D;
	x0=x;
	y0=y;
	z0=z;
};

Point3DRefinement::~Point3DRefinement()
{

};

Region3DRefinement::Region3DRefinement(double x00, double y00, double z00, 
										double x11, double y11, double z11)
{
	flavour=REGION3D;
	x0=x00;
	y0=y00;
	z0=z00;
	x1=x11;
	y1=y11;
	z1=z11;
};

Region3DRefinement::~Region3DRefinement()
{

};

CircleRefinement::CircleRefinement(double x, double y, double r0)
{
	flavour=CIRCLE;
	x0=x;
	y0=y;
	r=r0;
};

CircleRefinement::~CircleRefinement()
{

};

SphereRefinement::SphereRefinement(double x, double y, double z, double r0)
{
	flavour=SPHERE;
	x0=x;
	y0=y;
	z0=z;
	r=r0;
};

SphereRefinement::~SphereRefinement()
{

};


Border2DRefinement::Border2DRefinement(Border border, double dx)
{
	if(border == TOP || border == BOTTOM)
		throw oxley::OxleyException("Invalid border.");

	b = border;
	depth = dx;
};

Border2DRefinement::~Border2DRefinement()
{

};


Border3DRefinement::Border3DRefinement(Border border, double dx)
{
	b = border;
	depth = dx;
};

Border3DRefinement::~Border3DRefinement()
{

};
