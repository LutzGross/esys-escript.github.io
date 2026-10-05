
#include "oxley/RefinementQueue.h"
#include "oxley/OxleyException.h"

#include <cctype>
#include <string>

// apply() is NOT defined here. It needs the concrete domain classes, and
// including Rectangle.h/Brick.h from this file walks into the OxleyData.h <->
// Rectangle.h include cycle. Each apply() therefore lives in the translation
// unit of the domain it builds: RefinementQueue2D::apply in Rectangle.cpp,
// RefinementQueue3D::apply in Brick.cpp.

namespace oxley {

namespace {
/// a mask tag is handed to apply() as a keyword argument, apply(domain,
/// tag=mask), so it has to be a python identifier, and cannot be "domain"
void checkMaskTag(const std::string& tag)
{
    bool ok = !tag.empty() && (std::isalpha((unsigned char) tag[0]) || tag[0] == '_');
    for(size_t i = 1; ok && i < tag.size(); i++)
        ok = std::isalnum((unsigned char) tag[i]) || tag[i] == '_';
    if(!ok)
        throw OxleyException("refineMask: the tag '" + tag + "' is not a valid "
                "name. It is passed to apply() as a keyword argument, so it "
                "must be a python identifier.");
    if(tag == "domain")
        throw OxleyException("refineMask: 'domain' cannot be a tag, it is "
                "the name of apply()'s first argument.");
}
} // anonymous namespace

void RefinementQueue::checkMasks(const MaskMap& masks, const std::string& who)
{
    std::set<std::string> tags;
    for(size_t i = 0; i < queue.size(); i++)
        if(queue[i].flavour == MASK2D || queue[i].flavour == MASK3D)
            tags.insert(queue[i].tag);
    for(const std::string& tag : tags)
        if(masks.find(tag) == masks.end())
            throw OxleyException(who + ": no mask given for the tag '" + tag
                    + "'. Pass it as apply(domain, " + tag + "=mask).");
    // an unknown name is most likely a misspelt tag, which would otherwise
    // silently refine nothing
    for(MaskMap::const_iterator m = masks.begin(); m != masks.end(); ++m)
        if(tags.count(m->first) == 0)
            throw OxleyException(who + ": '" + m->first + "' is not the tag "
                    "of any mask refinement in the queue.");
}

RefinementQueue::RefinementQueue()
{
	refinement_levels=0;
}

RefinementQueue::~RefinementQueue()
{

}

void RefinementQueue::addToQueue(RefinementType R)
{	
	queue.push_back(R);
}

RefinementType RefinementQueue::getRefinement(int n)
{
    if(n <= getNumberOfOperations())
        return queue[n];
    else
        throw OxleyException("Number is greater than queue length");
}

void RefinementQueue::setRefinementLevel(int n)
{
	if(n >= 0)
    	refinement_levels=n;
    else
    	throw OxleyException("The levels of refinement must be equal to or greater than zero.");
}

void RefinementQueue::deleteFromQueue(int n)
{
	throw OxleyException("Unknown error.");
}

void RefinementQueue::print()
{
   	throw OxleyException("Unknown error.");
}

RefinementQueue2D::RefinementQueue2D()
{

}

RefinementQueue2D::~RefinementQueue2D()
{

}

void RefinementQueue2D::refineUniform(int level)
{
    if(level == -1)
        level=refinement_levels;
	RefinementType refine;
	refine.UniformRefinement(level);
	addToQueue(refine);
}

void RefinementQueue2D::refinePoint(float x0, float y0, int level)
{
    if(level == -1)
        level=refinement_levels;
	RefinementType refine;
	refine.Point2DRefinement(x0,y0,level);
	addToQueue(refine);
}

void RefinementQueue2D::refineRegion(float x0, float y0, float x1, float y1, int level)
{
    if(level == -1)
        level=refinement_levels;
	RefinementType refine;
	refine.Region2DRefinement(x0,y0,x1,y1,level);
	addToQueue(refine);
}

void RefinementQueue2D::refineCircle(float x0, float y0, float r, int level)
{
    if(level == -1)
        level=refinement_levels;
	RefinementType refine;
	refine.CircleRefinement(x0,y0,r,level);
	addToQueue(refine);
}

void RefinementQueue2D::refineBorder(Border b, float dx, int level)
{
    if(level == -1)
        level=refinement_levels;
	RefinementType refine;
	refine.Border2DRefinement(b,dx,level);
	addToQueue(refine);
}

void RefinementQueue2D::refineMask(std::string tag, int level)
{
    checkMaskTag(tag);
    if(level == -1)
        level=refinement_levels;
	RefinementType refine;
	refine.Mask2DRefinement(tag,level);
	addToQueue(refine);
}

void RefinementQueue2D::print()
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
    			std::cout << "Mask '" << Refinement.tag << "', level=" << l << std::endl;
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
}

void RefinementQueue2D::deleteFromQueue(int n)
{
	if(n <= queue.size())
		queue.erase(queue.begin()+n-1);
	else
		throw OxleyException("n must be smaller than the length of the queue");
}

RefinementQueue3D::RefinementQueue3D()
{

}

RefinementQueue3D::~RefinementQueue3D()
{

}

void RefinementQueue3D::refineUniform(int level)
{
    if(level == -1)
        level=refinement_levels;
	RefinementType refine;
	refine.UniformRefinement(level);
	addToQueue(refine);
}

void RefinementQueue3D::refinePoint(float x0, float y0, float z0, int level)
{
    if(level == -1)
        level=refinement_levels;
	RefinementType refine;
	refine.Point3DRefinement(x0,y0,z0,level);
	addToQueue(refine);
}

void RefinementQueue3D::refineRegion(float x0, float y0, float z0, float x1, float y1, float z1, int level)
{
    if(level == -1)
        level=refinement_levels;
	RefinementType refine;
	refine.Region3DRefinement(x0,y0,z0,x1,y1,z1,level);
	addToQueue(refine);
}

void RefinementQueue3D::refineSphere(float x0, float y0, float z0, float r, int level)
{
    if(level == -1)
        level=refinement_levels;
	RefinementType refine;
	refine.SphereRefinement(x0,y0,z0,r,level);
	addToQueue(refine);
}

void RefinementQueue3D::refineBorder(Border b, float dx, int level)
{
    if(level == -1)
        level=refinement_levels;
	RefinementType refine;
	refine.Border3DRefinement(b,dx,level);
	addToQueue(refine);
}

void RefinementQueue3D::refineMask(std::string tag, int level)
{
    checkMaskTag(tag);
    if(level == -1)
        level=refinement_levels;
	RefinementType refine;
	refine.Mask3DRefinement(tag,level);
	addToQueue(refine);
}

void RefinementQueue3D::print()
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
                        std::cout << " (North)";
                        break;
                    }
                    case SOUTH:
                    {
                        std::cout << " (South)";
                        break;
                    }
                    case WEST:
                    {
                        std::cout << " (West)";
                        break;
                    }
                    case EAST:
                    {
                        std::cout << " (East)";
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
    			std::cout << "Mask '" << Refinement.tag << "', level=" << l << std::endl;
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
}

void RefinementQueue3D::deleteFromQueue(int n)
{
	if(n <= queue.size())
		queue.erase(queue.begin()+n-1);
	else
		throw OxleyException("n must be smaller than the length of the queue");
}


namespace {
/// maps a border name onto the enum, accepting both the geographic and the
/// geometric spelling, which the old domain-level refineBoundary took too
Border borderFromName(std::string name, int dim)
{
    for(size_t i = 0; i < name.size(); ++i)
        name[i] = std::tolower(name[i]);
    if(name == "left" || name == "west")     return WEST;    // x minimal
    if(name == "right" || name == "east")    return EAST;    // x maximal
    if(dim == 2)
    {
        if(name == "top" || name == "north")     return NORTH;   // y maximal
        if(name == "bottom" || name == "south")  return SOUTH;   // y minimal
    }
    else
    {
        // in 3D top and bottom are the faces normal to z, north and south
        // (back and front) those normal to y
        if(name == "north" || name == "back")    return NORTH;   // y maximal
        if(name == "south" || name == "front")   return SOUTH;   // y minimal
        if(name == "top")                        return TOP;     // z maximal
        if(name == "bottom")                     return BOTTOM;  // z minimal
    }
    throw OxleyException("refineBorder: unknown border '" + name + "'.");
}
} // anonymous namespace

void RefinementQueue2D::refineBorder(std::string border, float dx, int level)
{
    refineBorder(borderFromName(border, 2), dx, level);
}

void RefinementQueue3D::refineBorder(std::string border, float dx, int level)
{
    refineBorder(borderFromName(border, 3), dx, level);
}

} //namespace oxley
