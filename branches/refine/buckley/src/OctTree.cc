
#include "OctTree.h"
#include <cmath>
#include <cstring>  // for memset
#include "FaceConsts.h"

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

void clearIDs(const OctCell& c, void* v)
{
    if (c.leaf)
    {
        memset(const_cast<OctCell&>(c).leafinfo->pmap,0,sizeof(unkid)*8);
    }
}

// set the unknown ids in pmap
// the void argument is a pair <childnum, nextid>
// where childnum is the position of this cell in its parents array and nextid is the next available id to give to a corner
void setUnkids(const OctCell& c, int kidnum, void* v)
{
    bool bigger[6];	// record if the neighbour cell on that face is bigger/shallower than us
    for (unsigned i=0;i<6;++i)
    {
	bigger[i]=(c.leafinfo->next[i])?
		    (c.leafinfo->next[i]->owner->depth < c.depth) :
		    false;
    }
    // we need to assign an id to each corner or mark it as hanging
    for (unsigned i=0;i<8;++i)
    {
	if (!c.leafinfo->pmap[i])	// no id yet
	{
	    // now we need to workout whether the node hangs
	    if ((canhang[kidnum][i]) && (bigger[facestouch[i][0]] || bigger[facestouch[i][1]] || bigger[facestouch[i][2]]))
	    {
		// the node hangs
		c.leafinfo->pmap[i]=HANG_NODE;
	    }
	    else
	    {
	        // the node does not hang
		c.leafinfo->pmap[i]=(*reinterpret_cast<unkid*>(v))++;	      
	    }
	}
    }
}


}

void OctTree::assignIDs()
{
    unsigned id=0;
    p.doLeafWalk(zorderlabel, &id);
    p.doLeafWalk(clearIDs, 0);
    unkid i=HANG_NODE+1;
    p.doLeafWalkWithKids(setUnkids, 0, &i);
}


void OctTree::walkLeaves(cellfunct c, void* v)
{
    p.doLeafWalk(c, v);  
}



void OctTree::walkWithKids(cellfunct2 c, void* v)
{
    p.doLeafWalkWithKids(c, 0, v);  	// we need to pick a childnum to start but a level 0 leaf is a trivial
}					// tree anyway

unsigned  OctTree::leafCount()
{
    unsigned c;
    p.doLeafWalk(countleaves, &c);
    return c;
}

void OctTree::debug()
{
   p.debug(false);
}

