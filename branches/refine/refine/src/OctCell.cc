#include <iostream>
#include "OctCell.h"

using namespace refine;

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

// find the leaf containing the point
// If the point lies on a boundary, then the code will pick a side
// This code assumes that the specified point must lie somewhere in this cell
//
// This code does not enforce smooth stepping down (ie we could have 1 cell neighbouring >2 cells)
OctCell* OctCell::splitPoint(double x, double y, double z, unsigned int desireddepth)
{
cerr << "SP " << x << " " << y << " " << z << "(" << desireddepth << ") at " << centre[0] << " " << centre[1] << " " << centre[2] << " ("
      << depth << ") " << endl;
    if (depth>=desireddepth)
    {
       return 0; 
    }
    if (leaf)
    {
        if (depth+1==desireddepth)	// note that this level is not checked for outward refinement
	{				// I'm assuming the final caller is doing that
	    cerr << "USplit at " << centre[0] << " " << centre[1] << " " << centre[2] << endl;
	    split();
	    return this;
	}
	else
	{
	    cerr << "Splitting Point " << centre[0] << " " << centre[1] << " " << centre[2] << endl;
	    // this situation means someone is trying to refine a cell multiple levels in one call
	    // not really what we had in mind so ...
	    split();
	    cerr << "My depth is " << depth << " and I created a cell with depth " << (depth+1) << endl;
	    outwardRefine(depth-1);
	    // now we need to do this again
	    
	    return splitPoint(x, y, z, desireddepth);
	}
    }
    else
    {
        int pz=(z>centre[2])?4:0;
        int py=(y>centre[1])?2:0;
        int px=(x>centre[0]);
        return kids[pz+py+px]->splitPoint(x, y, z, desireddepth);
    }
}

// After a cell has been split, this method can be called to ensure that the tree is still "not too unbalanced"
void OctCell::outwardRefine(unsigned desireddepth)
{
    double minside=(sides[0]>sides[1])?sides[0]:sides[1];
    minside=(minside>sides[2])?sides[2]:minside;
    minside/=5;		// this will work provided that all other divisions have been safe 
    
    cerr << "Outward refine: " << centre[0] << " " << centre[1] << " " << centre[2] << "\n";
    cerr << "-\n";
    upSplitPoint(centre[0]-sides[0]/2-minside, centre[1], centre[2], desireddepth);
    cerr << "-\n";
    upSplitPoint(centre[0]+sides[0]/2+minside, centre[1], centre[2], desireddepth);
    cerr << "-\n";
    upSplitPoint(centre[0], centre[1]-sides[1]/2-minside, centre[2], desireddepth);
    cerr << "-\n";
    upSplitPoint(centre[0], centre[1]+sides[1]/2+minside, centre[2], desireddepth);
    cerr << "-\n";
    upSplitPoint(centre[0], centre[1], centre[2]-sides[2]/2-minside, desireddepth);
    cerr << "-\n";
    upSplitPoint(centre[0], centre[1], centre[2]+sides[2]/2+minside, desireddepth);
    cerr << "----\n";
}


// Does the same as above but ensures the tree remains "not too unbalanced"
OctCell* OctCell::safeSplitPoint(double x, double y, double z, unsigned int desireddepth)
{
    if (!desireddepth)
    {
        return 0;  
    }
    OctCell* s=splitPoint(x,y,z, desireddepth);  
    if (!s)
    {
        return 0;  
    }
    s->outwardRefine(desireddepth-1);
    return this;
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
        cerr << "Up to " << centre[0] << " " << centre[1] << " " << centre[2] << "\n";
      
        // it's in this cell somewhere so start moving down again
        splitPoint(x, y, z, d);       
    } 
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




