
#include "OctCell.h"

using namespace refine;

OctCell::OctCell(double cx, double cy, double cz, double wx, double wy, double wz)
{
    leaf=true;
    centre[0]=cx;
    centre[1]=cy;
    centre[2]=cz;
    sides[0]=wx;
    sides[1]=wy;
    sides[2]=wz;
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
	kids[0]=new OctCell(cx-dx/2, cy-dy/2, cz-dz/2, dx, dy, dz);
	kids[1]=new OctCell(cx+dx/2, cy-dy/2, cz-dz/2, dx, dy, dz);
	kids[2]=new OctCell(cx+dx/2, cy+dy/2, cz-dz/2, dx, dy, dz);
	kids[3]=new OctCell(cx-dx/2, cy+dy/2, cz-dz/2, dx, dy, dz);
	kids[4]=new OctCell(cx-dx/2, cy-dy/2, cz+dz/2, dx, dy, dz);
	kids[5]=new OctCell(cx+dx/2, cy-dy/2, cz+dz/2, dx, dy, dz);
	kids[6]=new OctCell(cx+dx/2, cy+dy/2, cz+dz/2, dx, dy, dz);
	kids[7]=new OctCell(cx-dx/2, cy+dy/2, cz+dz/2, dx, dy, dz);
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
void OctCell::splitPoint(double x, double y, double z)
{
    if (leaf)
    {
        split();    
    }
    else
    {
        int pz=(z>centre[2])?4:0;
        int py=(y>centre[1])?2:0;
        int px=(x>centre[0]);
        kids[pz+py+px]->splitPoint(x, y, z);
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

