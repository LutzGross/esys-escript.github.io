
#include "OctTree.h"
#include <cmath>


using namespace buckley;


// Note: This code assumes x,y,z are positive
OctTree::OctTree(double x, double y, double z)
:p(fabs(x)/2, fabs(y)/2, fabs(z)/2, fabs(x), fabs(y), fabs(z),0)
{
    side[0]=fabs(x);
    side[1]=fabs(y);
    side[2]=fabs(z);
}

OctTree::~OctTree(){}

// guarantees that all leaves are at depth >=d
void OctTree::allSplit(unsigned int d)
{
    p.allSplit(d);  
}

void OctTree::allCollapse(unsigned d)
{
   p.collapseAll(d);  
  
}

void OctTree::collapsePoint(double x, double y, double z, unsigned int d)
{
    p.collapsePoint(x,y,z,d);  
  
}

void OctTree::splitPoint(double x, double y, double z, unsigned int desdepth)
{
    if ((x<0) || (y<0) || (z<0) ||
        (x>side[0]) || (y>side[1]) || (z>side[2]))
    {
        return;
    }
    p.splitPoint(x, y, z, desdepth);
}

namespace 
{

// I haven't decided whether I want to have a generally mutable interface to nodes or not
// the use of const_cast is ok coz it's me doing it
void zorderlabel(const OctCell& c, void* v)
{
    const_cast<OctCell&>(c).id=(*reinterpret_cast<unsigned*>(v))++;  
}

void countleaves(const OctCell& c, void* v)
{
    (*reinterpret_cast<int*>(v))++;  
}
  
}

void OctTree::assignIDs()
{
    unsigned id=0;
    p.doLeafWalk(zorderlabel, &id);
}


void OctTree::walkLeaves(cellfunct c, void* v)
{
    p.doLeafWalk(c, v);  
}

unsigned  OctTree::leafCount()
{
    unsigned c;
    p.doLeafWalk(countleaves, &c);
    return c;
}

