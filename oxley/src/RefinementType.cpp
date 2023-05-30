


#include <oxley/OxleyException.h>
#include <oxley/RefinementType.h>

namespace oxley {

RefinementType::RefinementType()
{
	x0	 = 0.0;
	y0	 = 0.0;
	z0	 = 0.0;
	x1	 = 0.0;
	y1	 = 0.0;
	z1	 = 0.0;
	r 	 = 0.0;
	dims 	 = 1;
	b 	 = NORTH;
	depth = 0.0;
	flavour = POINT2D;
	levels = 1;
	data = new escript::Data();
};

RefinementType::~RefinementType()
{
	// delete [] data;
};

void RefinementType::Point2DRefinement(double x, double y, int r_levels)
{
	levels=r_levels;
	flavour=POINT2D;
	x0=x;
	y0=y;
};

void RefinementType::Point2DRefinement()
{
	levels=-1;
	flavour=POINT2D;
	x0=-1;
	y0=-1;
};

void RefinementType::Region2DRefinement(double x00, double y00, double x11, double y11, int r_levels)
{
	levels=r_levels;
	flavour=REGION2D;
	x0=x00;
	y0=y00;
	x1=x11;
	y1=y11;
};

void RefinementType::Region2DRefinement()
{
	levels=-1;
	flavour=REGION2D;
	x0=-1;
	y0=-1;
	x1=-1;
	y1=-1;
};

void RefinementType::Point3DRefinement(double x, double y, double z, int r_levels)
{
	levels=r_levels;
	flavour=POINT3D;
	x0=x;
	y0=y;
	z0=z;
};

void RefinementType::Point3DRefinement()
{
	levels=-1;
	flavour=POINT3D;
	x0=-1;
	y0=-1;
	z0=-1;
};

void RefinementType::Region3DRefinement(double x00, double y00, double z00, 
										double x11, double y11, double z11, int r_levels)
{
	levels=r_levels;
	flavour=REGION3D;
	x0=x00;
	y0=y00;
	z0=z00;
	x1=x11;
	y1=y11;
	z1=z11;
};

void RefinementType::Region3DRefinement()
{
	levels=-1;
	flavour=REGION3D;
	x0=-1;
	y0=-1;
	z0=-1;
	x1=-1;
	y1=-1;
	z1=-1;
};

void RefinementType::CircleRefinement(double x, double y, double r0, int r_levels)
{
	levels=r_levels;
	flavour=CIRCLE;
	x0=x;
	y0=y;
	r=r0;
};

void RefinementType::CircleRefinement()
{
	levels=-1;
	flavour=CIRCLE;
	x0=-1;
	y0=-1;
	r=-1;
};

void RefinementType::SphereRefinement(double x, double y, double z, double r0, int r_levels)
{
	levels=r_levels;
	flavour=SPHERE;
	x0=x;
	y0=y;
	z0=z;
	r=r0;
};

void RefinementType::SphereRefinement()
{
	levels=-1;
	flavour=SPHERE;
	x0=-1;
	y0=-1;
	z0=-1;
	r=-1;
};

void RefinementType::Border2DRefinement(Border border, double dx, int r_levels)
{
	levels=r_levels;
	flavour=BOUNDARY;

	if(border == TOP || border == BOTTOM)
		throw oxley::OxleyException("Invalid border.");

	b = border;
	depth = dx;
};

void RefinementType::Border2DRefinement()
{
	levels=-1;
	flavour=BOUNDARY;
	b = TOP;
	depth = -1;
};

void RefinementType::Border3DRefinement(Border border, double dx, int r_levels)
{
	levels=r_levels;
	flavour=BOUNDARY;

	b = border;
	depth = dx;
};

void RefinementType::Border3DRefinement()
{
	levels=-1;
	flavour=BOUNDARY;

	b = TOP;
	depth = -1;
};

void RefinementType::Mask2DRefinement(escript::Data * d, int r_levels)
{
	levels=r_levels;
	flavour=MASK2D;
	data->copy(*d);
	// data = new escript::Data(d->copySelf());
};

void RefinementType::Mask3DRefinement(escript::Data * d, int r_levels)
{
	levels=r_levels;
	flavour=MASK3D;
	data->copy(*d);
	// data = new escript::Data(d->copySelf());
};

} //namespace oxley