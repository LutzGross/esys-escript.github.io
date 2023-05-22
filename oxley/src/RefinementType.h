

enum RefinementAlgorithm { POINT, REGION2D, REGION3D, MASK, CIRCLE, SPHERE };
enum Border { NORTH, SOUTH, EAST, WEST, TOP, BOTTOM };

/**
	\brief
	Abstract refinement type class
*/
class RefinementType
{

public:

	RefinementAlgorithm flavour;

};

class PointRefinement : public RefinementType
{
public:
	double x0, y0;
	
	/**
       \brief
       Constructor
    */
	PointRefinement(double x, double y);

	/**
       \brief
       Destructor
    */
	~PointRefinement();
};

class Region2DRefinement : public RefinementType
{
public:
	double x0, y0;
	double x1, y1;

	/**
       \brief
       Constructor
    */
	Region2DRefinement(double x0, double y0, double x1, double y1);

	/**
       \brief
       Destructor
    */
	~Region2DRefinement();
};

class Region3DRefinement : public RefinementType
{
public:
	int x0, y0, z0;
	int x1, y1, z1;

	/**
       \brief
       Constructor
    */
	Region3DRefinement(double x0, double y0, double z0, double x1, double y1, double z1);

	/**
       \brief
       Destructor
    */
	~Region3DRefinement();
};

class CircleRefinement : public RefinementType
{
public:
	int x0, y0;
	int r;

	/**
       \brief
       Constructor
    */
	CircleRefinement(double x0, double y0, double r);

	/**
       \brief
       Destructor
    */
	~CircleRefinement();
};

class SphereRefinement : public RefinementType
{
public:
	int x0, y0, z0;
	int r;

	/**
       \brief
       Constructor
    */
	SphereRefinement(double x0, double y0, double z0, double r);

	/**
       \brief
       Destructor
    */
	~SphereRefinement();
};

class Border2DRefinement : public RefinementType
{
public:
	Border b;
	double depth;

	/**
       \brief
       Constructor
    */
	Border2DRefinement(Border b, double depth);

	/**
       \brief
       Destructor
    */
	~Border2DRefinement();
};

class Border3DRefinement : public RefinementType
{
public:
	Border b;
	double depth;

	/**
       \brief
       Constructor
    */
	Border3DRefinement(Border b, double depth);

	/**
       \brief
       Destructor
    */
	~Border3DRefinement();
};