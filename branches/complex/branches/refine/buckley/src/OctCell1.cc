#include <iostream>
#include "OctCell.h"

using namespace buckley;

using namespace std;

namespace 
{
   const double eps=1e-1;  
  
}

OctCell::OctCell(double cx, double cy, double cz, double wx, double wy, double wz, OctCell* progenitor)
{
    leaf=true;
    centre[0]=cx;
    centre[1]=cy;
    centre[2]=cz;
    sides[0]=wx;
    sides[1]=wy;
    sides[2]=wz;
    parent=progenitor;
    if (parent)
    {
        depth=parent->depth+1;
    }
    else
    {
        depth=0;
    }
    id=0;
}

OctCell::~OctCell()
{
    if (!leaf)
    {
        for (unsigned int i=0;i<8;++i)
        {
            delete kids[i];
        }
    }
}

// If we are allocating from a pool we need to keep a record of it somewhere
// maybe keep track (static?) of holding object with operators
// or define OctCell::new ?
void OctCell::split()
{
    if (leaf)
    {
        const double dx=sides[0]/2;
	const double dy=sides[1]/2;
	const double dz=sides[2]/2;
	const double cx=centre[0];
	const double cy=centre[1];
	const double cz=centre[2];
	kids[0]=new OctCell(cx-dx/2, cy-dy/2, cz-dz/2, dx, dy, dz, this);
	kids[1]=new OctCell(cx+dx/2, cy-dy/2, cz-dz/2, dx, dy, dz, this);
	kids[2]=new OctCell(cx+dx/2, cy+dy/2, cz-dz/2, dx, dy, dz, this);
	kids[3]=new OctCell(cx-dx/2, cy+dy/2, cz-dz/2, dx, dy, dz, this);
	kids[4]=new OctCell(cx-dx/2, cy-dy/2, cz+dz/2, dx, dy, dz, this);
	kids[5]=new OctCell(cx+dx/2, cy-dy/2, cz+dz/2, dx, dy, dz, this);
	kids[6]=new OctCell(cx+dx/2, cy+dy/2, cz+dz/2, dx, dy, dz, this);
	kids[7]=new OctCell(cx-dx/2, cy+dy/2, cz+dz/2, dx, dy, dz, this);
	leaf=false;
    }  
}

// could use the replacement delete operator to notify other parts of the system that a value has been removed.
void OctCell::merge()
{
    if (!leaf)
    {
        for (unsigned int i=0;i<8;i++)
	{
            delete kids[i];	  
	}
	leaf=true;
    }  
}

void OctCell::allSplit(unsigned int depth)
{
    if (depth>0)  
    {
        if (leaf)
	{
	    split();    
	}
	if (depth==1)
	{
	   return; 
	}
        for (unsigned int i=0;i<8;i++)
	{
	    kids[i]->split();
	    kids[i]->allSplit(depth-1);
	}
    }
}


// After a cell has been split, this method can be called to ensure that the tree is still "not too unbalanced"
void OctCell::outwardRefine(unsigned desireddepth)
{
    double minside=(sides[0]>sides[1])?sides[0]:sides[1];
    minside=(minside>sides[2])?sides[2]:minside;
    minside/=5;		// this will work provided that all other divisions have been safe 
    
//    cerr << "Outward buckley: " << centre[0] << " " << centre[1] << " " << centre[2] << "\n";
//    cerr << "-\n";
    upSplitPoint(centre[0]-sides[0]/2-minside, centre[1], centre[2], desireddepth);
//    cerr << "-\n";
    upSplitPoint(centre[0]+sides[0]/2+minside, centre[1], centre[2], desireddepth);
//    cerr << "-\n";
    upSplitPoint(centre[0], centre[1]-sides[1]/2-minside, centre[2], desireddepth);
//    cerr << "-\n";
    upSplitPoint(centre[0], centre[1]+sides[1]/2+minside, centre[2], desireddepth);
//    cerr << "-\n";
    upSplitPoint(centre[0], centre[1], centre[2]-sides[2]/2-minside, desireddepth);
//    cerr << "-\n";
    upSplitPoint(centre[0], centre[1], centre[2]+sides[2]/2+minside, desireddepth);
//    cerr << "----\n";
}

// Works up  the tree until it finds 
void OctCell::upSplitPoint(double x, double y, double z, unsigned d)
{
    if ((x<centre[0]-sides[0]/2) || (x>centre[0]+sides[0]/2) || (y<centre[1]-sides[1]/2) ||
        (y>centre[1]+sides[1]/2) || (z<centre[2]-sides[2]/2) || (z>centre[2]+sides[2]/2))
    {
        // It is outside this box so go up
	if (parent!=0)
	{
	    parent->upSplitPoint(x, y, z, d);
	}
        return;
    }
    else
    {
//        cerr << "Up to " << centre[0] << " " << centre[1] << " " << centre[2] << "\n";
      
        // it's in this cell somewhere so start moving down again
        splitPoint(x, y, z, d);       
    } 
}

// Works up  the tree until it finds 
void OctCell::upCollPoint(double x, double y, double z, unsigned d)
{
    if ((x<centre[0]-sides[0]/2) || (x>centre[0]+sides[0]/2) || (y<centre[1]-sides[1]/2) ||
        (y>centre[1]+sides[1]/2) || (z<centre[2]-sides[2]/2) || (z>centre[2]+sides[2]/2))
    {
        // It is outside this box so go up
	if (parent!=0)
	{
	    parent->upCollPoint(x, y, z, d);
	}
        return;
    }
    else
    {
//        cerr << "Up to " << centre[0] << " " << centre[1] << " " << centre[2] << "\n";
      
        // it's in this cell somewhere so start moving down again
        collapsePoint(x, y, z, d);       
    } 
}


// return the leaf which contains this point.
// Assumes that the point lies somewhere in this cell
OctCell* OctCell::findLeaf(double x, double y, double z)
{
    if (leaf)
    {
        return this;
    }
    int pz=(z>centre[2])?4:0;
    int py=(y>centre[1]);
    int px=(x>centre[0]);
    int child=pz+2*py+(px^py);	      
    return kids[child]->findLeaf(x,y,z);
}


// divide the leaf containing the point (if required)
// This needs to be more efficient but will do for now
void OctCell::splitPoint(double x, double y, double z, unsigned d)
{
    // find cell
    // if it doesn't need refining, then bail
    // buckley neighbours to be at least d-1
    // buckley this cell
    OctCell* start=this;
    do
    {
	OctCell* l=start->findLeaf(x, y, z);
	if (l->depth>=d)
	{
	    return;  
	}
	l->outwardRefine(l->depth);
	l->split();
	start=l;
    } while (start->depth+1 <d);
}


// collapses the children of this node
void OctCell::collapse()
{
    if (!leaf)
    {
        // collapse any grandkids
        for (unsigned i=0;i<8;++i)
	{
            kids[i]->collapse();	// so we know we are only dealing with leaves
	}
        for (unsigned i=0;i<8;++i)
	{
	    delete kids[i];
	}
	leaf=true;
    }  
}

void OctCell::collapseAll(unsigned int desdepth)
{
    if (leaf)
    {
        return;  
    }
    for (unsigned i=0;i<8;++i)
    {
	kids[i]->collapseAll(desdepth);
    }
    if (depth>=desdepth)
    {
        collapse();
    }
}


// The leaf which contains (x,y,z) should be at depth at most d
void OctCell::collapsePoint(double x, double y, double z, unsigned d)
{
    // find cell
    // if it doesn't need refining, then bail
    // buckley neighbours to be at least d-1
    // buckley this cell
    OctCell* start=this;
    do
    {
	OctCell* l=start->findLeaf(x, y, z);
	if (l->depth<=d)
	{
	    return;  
	}
	l=l->parent;
	l->outwardCollapse(l->depth+1);		// since we have gone up a level
	l->collapse();
	start=l;
    } while (start->depth+1 >d);
  
}


// This is horribly inefficient!!!
// After a cell has been, thisdegrem method can be called to ensure that the tree is still "not too unbalanced"
void OctCell::outwardCollapse(unsigned desireddepth)
{
    double minside=(sides[0]>sides[1])?sides[0]:sides[1];
    minside=(minside>sides[2])?sides[2]:minside;
    minside/=5;		// this will work provided that all other divisions have been safe 
    
//    cerr << "Outward collapse: " << centre[0] << " " << centre[1] << " " << centre[2] << "\n";
//    cerr << "-" << (centre[0]-sides[0]/2-minside) << ' ' <<  centre[1] << ' ' <<  centre[2]<<"\n";
    upCollPoint(centre[0]-sides[0]/2-minside, centre[1]-minside, centre[2]-minside, desireddepth);
    
    upCollPoint(centre[0]-sides[0]/2-minside, centre[1]-minside, centre[2]+minside, desireddepth);
    
    upCollPoint(centre[0]-sides[0]/2-minside, centre[1]+minside, centre[2]-minside, desireddepth);
    
    upCollPoint(centre[0]-sides[0]/2-minside, centre[1]+minside, centre[2]+minside, desireddepth);
    
    
//    cerr << "-\n";
    upCollPoint(centre[0]+sides[0]/2+minside, centre[1]-minside, centre[2]-minside, desireddepth);
    
    upCollPoint(centre[0]+sides[0]/2+minside, centre[1]-minside, centre[2]+minside, desireddepth);

    upCollPoint(centre[0]+sides[0]/2+minside, centre[1]+minside, centre[2]-minside, desireddepth);

    upCollPoint(centre[0]+sides[0]/2+minside, centre[1]+minside, centre[2]+minside, desireddepth);
    
    
//    cerr << "-\n";
    upCollPoint(centre[0]-minside, centre[1]-sides[1]/2-minside, centre[2]-minside, desireddepth);
    
    upCollPoint(centre[0]-minside, centre[1]-sides[1]/2-minside, centre[2]+minside, desireddepth);
    
    upCollPoint(centre[0]+minside, centre[1]-sides[1]/2-minside, centre[2]-minside, desireddepth);
    
    upCollPoint(centre[0]+minside, centre[1]-sides[1]/2-minside, centre[2]+minside, desireddepth);
    
    
    
    
    
//    cerr << "-\n";
    upCollPoint(centre[0]-minside, centre[1]+sides[1]/2+minside, centre[2]-minside, desireddepth);
    
    upCollPoint(centre[0]-minside, centre[1]+sides[1]/2+minside, centre[2]+minside, desireddepth);

    upCollPoint(centre[0]+minside, centre[1]+sides[1]/2+minside, centre[2]-minside, desireddepth);
    
    upCollPoint(centre[0]+minside, centre[1]+sides[1]/2+minside, centre[2]+minside, desireddepth);
    
    
//    cerr << "-\n";
    upCollPoint(centre[0]-minside, centre[1]-minside, centre[2]-sides[2]/2-minside, desireddepth);
    
    upCollPoint(centre[0]-minside, centre[1]+minside, centre[2]-sides[2]/2-minside, desireddepth);
    
    upCollPoint(centre[0]+minside, centre[1]-minside, centre[2]-sides[2]/2-minside, desireddepth);
    
    upCollPoint(centre[0]+minside, centre[1]+minside, centre[2]-sides[2]/2-minside, desireddepth);    
    
//    cerr << "-\n";
    upCollPoint(centre[0]-minside, centre[1]-minside, centre[2]+sides[2]/2+minside, desireddepth);

    upCollPoint(centre[0]-minside, centre[1]+minside, centre[2]+sides[2]/2+minside, desireddepth);

    upCollPoint(centre[0]+minside, centre[1]-minside, centre[2]+sides[2]/2+minside, desireddepth);

    upCollPoint(centre[0]+minside, centre[1]+minside, centre[2]+sides[2]/2+minside, desireddepth);
    
//    cerr << "----\n";
}

void OctCell::doLeafWalk(cellfunct c, void* v)
{
    if (leaf)
    {
        c(*this, v);
    }
    else
    {
        for (unsigned int i=0;i<8;++i)
        {
            kids[i]->doLeafWalk(c,v);
        }
    }
}




