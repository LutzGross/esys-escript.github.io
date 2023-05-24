
#ifndef _OXLEY_REFINEMENTTYPE
#define _OXLEY_REFINEMENTTYPE

namespace oxley {

enum RefinementAlgorithm { POINT2D, POINT3D, REGION2D, REGION3D, MASK, CIRCLE, SPHERE, BOUNDARY };
enum Border { NORTH, SOUTH, EAST, WEST, TOP, BOTTOM };

/**
	\brief
	Abstract refinement type class
*/
class RefinementType
{

public:

	double x0, y0, z0;
	double x1, y1, z1;
	double r;
	int dims;
	Border b;
	double depth;
	RefinementAlgorithm flavour;
	int levels;

	RefinementType();
	~RefinementType();

	void Point2DRefinement(double x, double y, int levels);
	void Point2DRefinement();

	void Region2DRefinement(double x00, double y00, double x11, double y11, int levels);
	void Region2DRefinement();
	void Point3DRefinement(double x, double y, double z, int levels);
	void Point3DRefinement();
	void Region3DRefinement(double x00, double y00, double z00, 
							double x11, double y11, double z11, int levels);
	void Region3DRefinement();
	void CircleRefinement(double x, double y, double r0, int levels);
	void CircleRefinement();
	void SphereRefinement(double x, double y, double z, double r0, int levels);
	void SphereRefinement();

	void Border2DRefinement(Border border, double dx, int levels);
	void Border2DRefinement();
	void Border3DRefinement(Border border, double dx, int levels);
	void Border3DRefinement();

};

} //namespace oxley

#endif //_OXLEY_REFINEMENTTYPE