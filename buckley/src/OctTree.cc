
#include "OctTree.h"
#include <cmath>
#include <cstring>  // for memset
#include "FaceConsts.h"

using namespace buckley;


// Note: This code assumes x,y,z are positive
OctTree::OctTree(double x, double y, double z)
:p(fabs(x)/2, fabs(y)/2, fabs(z)/2, fabs(x), fabs(y), fabs(z),0), modified(true)
{
    side[0]=fabs(x);
    side[1]=fabs(y);
    side[2]=fabs(z);
    leafcount=1;
}

OctTree::~OctTree(){}

// guarantees that all leaves are at depth >=d
void OctTree::allSplit(unsigned int d)
{
    modified=true;
    p.allSplit(d, &leafcount);  
}

void OctTree::allCollapse(unsigned d)
{
    modified=true;  
   p.collapseAll(d, &leafcount);    
}

void OctTree::collapsePoint(double x, double y, double z, unsigned int d)
{
    modified=true;  
    p.collapsePoint(x,y,z,d, &leafcount);  
  
}

void OctTree::splitPoint(double x, double y, double z, unsigned int desdepth)
{  
    if ((x<0) || (y<0) || (z<0) ||
        (x>side[0]) || (y>side[1]) || (z>side[2]))
    {
        return;
    }
    modified=true;    
    p.splitPoint(x, y, z, desdepth, &leafcount);
}

namespace 
{

void zorderlabel(const OctCell& c, void* v)
{
    c.id=(*reinterpret_cast<unsigned*>(v))++;  
}

void countleaves(const OctCell& c, void* v)
{
    (*reinterpret_cast<int*>(v))++;  
}

void clearIDs(const OctCell& c, void* v)
{
    if (c.leaf)
    {
        memset(c.leafinfo->pmap,0,sizeof(unkid)*8);
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


typedef struct
{
   const OctCell** ca;
   unkid id;  
} copier;

void copytoarr(const OctCell& c, void* v)
{
    copier* cs=reinterpret_cast<copier*>(v);
    cs->ca[cs->id++]=&c;
}


}


const OctCell** OctTree::process(unkid& numunk) const
{
    assignIDs();
    // nasty hack process for now
    // ultimately, this should be done by starting each thread midway through the traversal
    const OctCell** temp=new const OctCell*[leafcount];
    copier c;
    c.ca=temp;
    c.id=0;
    p.doLeafWalk_const(copytoarr, &c);
    numunk=c.id;
    return temp;
}

void OctTree::assignIDs() const
{
    unsigned id=0;
    p.doLeafWalk_const(zorderlabel, &id);
    p.doLeafWalk_const(clearIDs, 0);
    unkid i=HANG_NODE+1;
    p.doLeafWalkWithKids_const(setUnkids, 0, &i);
}


void OctTree::walkLeaves(const_cellfn c, void* v) const
{
    p.doLeafWalk_const(c, v);  
}

void OctTree::walkLeaves(cellfn c, void* v)
{
    p.doLeafWalk(c, v);  
}



void OctTree::walkWithKids(const_cellfn2 c, void* v)
{
    p.doLeafWalkWithKids_const(c, 0, v);  	// we need to pick a childnum to start but a level 0 leaf is a trivial
}					// tree anyway

unsigned  OctTree::leafCountByWalk() const
{
    unsigned c;
    p.doLeafWalk_const(countleaves, &c);
    return c;
}

unsigned  OctTree::leafCount() const
{
    return leafcount;
}

void OctTree::debug()
{
   p.debug(false);
}

