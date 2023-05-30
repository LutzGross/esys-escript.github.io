
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

void RefinementZone::deleteFromQueue(int n)
{
	throw OxleyException("Unknown error.");
}

void RefinementZone::print()
{
   	throw OxleyException("Unknown error.");
};

RefinementZone2D::RefinementZone2D()
{

};

RefinementZone2D::~RefinementZone2D()
{

};

void RefinementZone2D::refinePoint(float x0, float y0, int level)
{
    if(level == -1)
        level=refinement_levels;
	RefinementType refine;
	refine.Point2DRefinement(x0,y0,level);
	addToQueue(refine);
};

void RefinementZone2D::refineRegion(float x0, float y0, float x1, float y1, int level)
{
    if(level == -1)
        level=refinement_levels;
	RefinementType refine;
	refine.Region2DRefinement(x0,y0,x1,y1,level);
	addToQueue(refine);
};

void RefinementZone2D::refineCircle(float x0, float y0, float r, int level)
{
    if(level == -1)
        level=refinement_levels;
	RefinementType refine;
	refine.CircleRefinement(x0,y0,r,level);
	addToQueue(refine);
};

void RefinementZone2D::refineBorder(Border b, float dx, int level)
{
    if(level == -1)
        level=refinement_levels;
	RefinementType refine;
	refine.Border2DRefinement(b,dx,level);
	addToQueue(refine);
};

void RefinementZone2D::refineMask(escript::Data d)
{
	RefinementType refine;
	refine.Mask2DRefinement(&d,refinement_levels);
	addToQueue(refine);	
};

void RefinementZone2D::print()
{
	for(int i = 0; i < queue.size(); i++)
	{
		std::cout << i << ": ";
		RefinementType Refinement = queue[i];
		int l = Refinement.levels;
		switch(Refinement.flavour)
		{
			case POINT2D:
            {
                double x=Refinement.x0;
                double y=Refinement.y0;
                std::cout << "Point (" << x << ", " << y << "), level=" << l << std::endl;
                break;
            }
            case REGION2D:
            {
                double x0=Refinement.x0;
                double y0=Refinement.y0;
                double x1=Refinement.x1;
                double y1=Refinement.y1;
                std::cout << "Region (" << x0 << ", " << y0 << ") ("
                						<< x1 << ", " << y1 << "), level=" << l << std::endl;
                break;
            }
            case CIRCLE:
            {
                double x0=Refinement.x0;
                double y0=Refinement.y0;
                double r0=Refinement.r;
                std::cout << "Circle at (" << x0 << ", " << y0 << ") with r = " << r0 << ", level=" << l << std::endl;
                break;
            }
            case BOUNDARY:
            {
                double dx=Refinement.depth;
                std::cout << "Boundary ";
                switch(Refinement.b)
                {
                    case NORTH:
                    {
                        std::cout << " (Top)";
                        break;
                    }
                    case SOUTH:
                    {
                        std::cout << " (Bottom)";
                        break;
                    }
                    case WEST:
                    {
                        std::cout << " (Left)";
                        break;
                    }
                    case EAST:
                    {
                        std::cout << " (Right)";
                        break;
                    }
                    case TOP:
                    case BOTTOM:
                    default:
                    {
                        throw OxleyException("Invalid border direction.");
                    }
                }
                std::cout << " to depth " << dx << std::endl;
                break;
            }
        	case MASK2D:
    		{
    			std::cout << "A mask refinement " << std::endl;
    			break;
    		}
        	case POINT3D:
        	case REGION3D:
        	case MASK3D:
        	case SPHERE:
        	default:
    		{
    			throw OxleyException("Unknown RefinementType")	;
    		}
        }
	}
};

void RefinementZone2D::deleteFromQueue(int n)
{
	if(n <= queue.size())
		queue.erase(queue.begin()+n-1);
	else
		throw OxleyException("n must be smaller than the length of the queue");
}

RefinementZone3D::RefinementZone3D()
{

};

RefinementZone3D::~RefinementZone3D()
{

};

void RefinementZone3D::refinePoint(float x0, float y0, float z0, int level)
{
    if(level == -1)
        level=refinement_levels;
	RefinementType refine;
	refine.Point3DRefinement(x0,y0,z0,refinement_levels);
	addToQueue(refine);
};

void RefinementZone3D::refineRegion(float x0, float y0, float z0, float x1, float y1, float z1, int level)
{
    if(level == -1)
        level=refinement_levels;
	RefinementType refine;
	refine.Region3DRefinement(x0,y0,z0,x1,y1,z1,level);
	addToQueue(refine);
};

void RefinementZone3D::refineSphere(float x0, float y0, float z0, float r, int level)
{
    if(level == -1)
        level=refinement_levels;
	RefinementType refine;
	refine.SphereRefinement(x0,y0,z0,r,level);
	addToQueue(refine);
};

void RefinementZone3D::refineBorder(Border b, float dx, int level)
{
    if(level == -1)
        level=refinement_levels;
	RefinementType refine;
	refine.Border3DRefinement(b,dx,level);
	addToQueue(refine);
};

void RefinementZone3D::refineMask(escript::Data d)
{
	RefinementType refine;
	refine.Mask3DRefinement(&d, refinement_levels);
	addToQueue(refine);	
};

void RefinementZone3D::print()
{
	for(int i = 0; i < queue.size(); i++)
	{
		std::cout << i << ": ";
		RefinementType Refinement = queue[i];
		int l = Refinement.levels;
		switch(Refinement.flavour)
		{
			case POINT3D:
            {
                double x=Refinement.x0;
                double y=Refinement.y0;
                double z=Refinement.z0;
                std::cout << "Point (" << x << ", " << y << ", " << z << "), level=" << l << std::endl;
                break;
            }
            case REGION3D:
            {
                double x0=Refinement.x0;
                double y0=Refinement.y0;
                double z0=Refinement.z0;
                double x1=Refinement.x1;
                double y1=Refinement.y1;
                double z1=Refinement.z1;
                std::cout << "Region (" << x0 << ", " << y0 << ", " << z0 << ") ("
                						<< x1 << ", " << y1 << ", " << z1 << "), level=" << l << std::endl;
                break;
            }
            case SPHERE:
            {
                double x0=Refinement.x0;
                double y0=Refinement.y0;
                double z0=Refinement.z0;
                double r0=Refinement.r;
                std::cout << "Sphere at (" << x0 << ", " << y0 << ", " << z0 << ") with r = " 
                						   << r0 << ", level=" << l << std::endl;
                break;
            }
            case BOUNDARY:
            {
                double dx=Refinement.depth;
                std::cout << "Boundary ";
                switch(Refinement.b)
                {
                    case NORTH:
                    {
                        std::cout << " (Top)";
                        break;
                    }
                    case SOUTH:
                    {
                        std::cout << " (Bottom)";
                        break;
                    }
                    case WEST:
                    {
                        std::cout << " (Left)";
                        break;
                    }
                    case EAST:
                    {
                        std::cout << " (Right)";
                        break;
                    }
                    case TOP:
                	{
                		std::cout << " (Top)";
                		break;
                	}
                    case BOTTOM:
                    {
                		std::cout << " (Bottom)";
                		break;
                	}
                    default:
                    {
                        throw OxleyException("Invalid border direction.");
                    }
                }
                std::cout << " to depth " << dx << std::endl;
                break;
            }
        	case MASK3D:
    		{
    			std::cout << "A mask refinement " << std::endl;
    			break;
    		}
        	case POINT2D:
        	case REGION2D:
        	case MASK2D:
        	case CIRCLE:
        	default:
        		throw OxleyException("Unknown RefinementType");
        }
	}
};

void RefinementZone3D::deleteFromQueue(int n)
{
	if(n <= queue.size())
		queue.erase(queue.begin()+n-1);
	else
		throw OxleyException("n must be smaller than the length of the queue");
}

} //namespace oxley