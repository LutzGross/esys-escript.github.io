
#include "OctTree.h"
#include <cmath>


using namespace refine;


// Note: This code assumes x,y,z are positive
OctTree::OctTree(double x, double y, double z)
:p(fabs(x)/2, fabs(y)/2, fabs(z)/2, fabs(x), fabs(y), fabs(z))
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

void OctTree::splitPoint(double x, double y, double z)
{
    if ((x<0) || (y<0) || (z<0) ||
        (x>side[0]) || (y>side[1]) || (z>side[2]))
    {
        return;
    }
    p.splitPoint(x,y,z);
}


void OctTree::walkLeaves(cellfunct c, void* v)
{
    p.doLeafWalk(c, v);  
}

