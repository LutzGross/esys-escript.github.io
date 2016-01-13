
#include <map>
#include "OctTree.h"
#include <cmath>
#include <cstring>  // for memset
#include "FaceConsts.h"


#include <iostream>

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


void OctTree::getBB(double bb[6]) const
{
    bb[0]=p.centre[0]-p.sides[0]/2;
    bb[1]=p.centre[1]-p.sides[1]/2;
    bb[2]=p.centre[2]-p.sides[2]/2;
    
    bb[3]=p.centre[0]+p.sides[0]/2;
    bb[4]=p.centre[1]+p.sides[1]/2;
    bb[5]=p.centre[2]+p.sides[2]/2;
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
    memset(c.leafinfo->pmap,0,sizeof(unkid)*8);
    unsigned& max_depth=*reinterpret_cast<unsigned*>(v);
    if (max_depth<c.depth)
    {
        max_depth=c.depth;
    }
}

    
// later we may need to look at using the neighbour relationships to identify shared nodes or even use a hashmap
// but for now we'll just do this


typedef struct 
{
   unsigned x;
   unsigned y;
   unsigned z;
  
} triple;

inline bool operator<(const triple& t1, const triple& t2)
{
//    return (t1.x+t1.y+t1.z) < (t2.x+t2.y+t2.z);
    if (t1.x< t2.x) return true;
    if (t1.x > t2.x) return false;
    if (t1.y< t2.y) return true;
    if (t1.y > t2.y) return false;
    if (t1.z< t2.z) return true;
    if (t1.z > t2.z) return false;
    return false;
}

typedef std::map<triple, unkid> unkmap;

    
typedef struct
{
    unkmap pointmap;
    unsigned max_depth;  
    double factor;
    double xmin;  
    double ymin;
    double zmin;
    unkid nextid;
} unkstruct;
    
    
    

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
	    else   // the node does not hang
	    {
	        unkstruct& astruct=*reinterpret_cast<unkstruct*>(v);
		triple t;
		double x,y,z;
		// need to compute the spatial coords for this child
		// I suspect this will be horribly inefficient
		c.childCoords(i, x, y, z);
		t.x=unsigned((x-astruct.xmin)/astruct.factor);
		t.y=unsigned((y-astruct.ymin)/astruct.factor);
		t.z=unsigned((z-astruct.zmin)/astruct.factor);	
	        unkmap::iterator it=astruct.pointmap.find(t);
		if (it!=astruct.pointmap.end())
		{
		    c.leafinfo->pmap[i]=it->second;		  
		}
		else
		{
		    astruct.pointmap[t]=astruct.nextid;
		    c.leafinfo->pmap[i]=astruct.nextid++;
		}
	    }
	}
    }
}


typedef struct
{
   const OctCell** ca;
   unkid id;
   std::vector<const OctCell*>* face_cells;
} copier;


// This is not threadsafe!!!
void copytoarr(const OctCell& c, void* v)
{
    copier* cs=reinterpret_cast<copier*>(v);
    cs->ca[cs->id++]=&c;
    for (int i=0;i<6;++i)
    {
        if (c.leafinfo->next[i]==0)	// this face is external
	{
	    cs->face_cells[i].push_back(&c);    
	}
    }
}


}


const OctCell** OctTree::process(unkid& numunk, std::vector<const OctCell*> face_cells[6]) const
{
    unkid maxunk=assignIDs();
    // nasty hack process for now
    // ultimately, this should be done by starting each thread midway through the traversal
    const OctCell** temp=new const OctCell*[leafcount];
    copier c;
    c.ca=temp;
    c.id=0;
    c.face_cells=face_cells;
    p.doLeafWalk_const(copytoarr, &c);
    numunk=maxunk-(HANG_NODE+1);
    return temp;
}


unkid OctTree::assignIDs() const
{
  
    unkid id=0;
    p.doLeafWalk_const(zorderlabel, &id);
    unsigned maxdepth=0;
    p.doLeafWalk_const(clearIDs, &maxdepth);	// will also work out the max depth for us
    unkstruct astruct;  
    astruct.nextid=HANG_NODE+1;
    astruct.max_depth=maxdepth;
    astruct.xmin=p.centre[0]-side[0];
    astruct.ymin=p.centre[1]-side[1];
    astruct.zmin=p.centre[2]-side[2];
    astruct.factor=pow(2, -(int)maxdepth);
    p.doLeafWalkWithKids_const(setUnkids, 0, &astruct);
    return astruct.nextid;
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

