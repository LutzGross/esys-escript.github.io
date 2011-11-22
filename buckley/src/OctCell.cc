#include <iostream>
#include "OctCell.h"
#include "LeafInfo.h"

using namespace buckley;

using namespace std;

namespace 
{
   const double eps=1e-1;    
}


#include "FaceConsts.h"

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
    leafinfo=new LeafInfo(this);
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
    else
    {
        delete leafinfo;
    }
}

// perhaps this should be in the .h file to be inlined?
void OctCell::childCoords(unsigned  k, double& x, double& y, double& z) const
{
    x=centre[0];
    y=centre[1];
    z=centre[2];
    switch(k)
    {
      case 0: x-=sides[0]/2;  y-=sides[1]/2;  z-=sides[2]/2; break;  
      case 1: x+=sides[0]/2;  y-=sides[1]/2;  z-=sides[2]/2; break;  
      case 2: x+=sides[0]/2;  y+=sides[1]/2;  z-=sides[2]/2; break;  
      case 3: x-=sides[0]/2;  y+=sides[1]/2;  z-=sides[2]/2; break;  
      case 4: x-=sides[0]/2;  y-=sides[1]/2;  z+=sides[2]/2; break;  
      case 5: x+=sides[0]/2;  y-=sides[1]/2;  z+=sides[2]/2; break;  
      case 6: x+=sides[0]/2;  y+=sides[1]/2;  z+=sides[2]/2; break;  
      default: x-=sides[0]/2;  y+=sides[1]/2;  z+=sides[2]/2; break;  
    };
}

// get coords of quadrature points
void OctCell::quadCoords(unsigned k, double& x, double& y, double& z) const
{
    x=centre[0];
    y=centre[1];
    z=centre[2];
    switch(k)
    {
      case 0: x-=sides[0]/4;  y-=sides[1]/4;  z-=sides[2]/4; break;  
      case 1: x+=sides[0]/4;  y-=sides[1]/4;  z-=sides[2]/4; break;  
      case 2: x+=sides[0]/4;  y+=sides[1]/4;  z-=sides[2]/4; break;  
      case 3: x-=sides[0]/4;  y+=sides[1]/4;  z-=sides[2]/4; break;  
      case 4: x-=sides[0]/4;  y-=sides[1]/4;  z+=sides[2]/4; break;  
      case 5: x+=sides[0]/4;  y-=sides[1]/4;  z+=sides[2]/4; break;  
      case 6: x+=sides[0]/4;  y+=sides[1]/4;  z+=sides[2]/4; break;  
      default: x-=sides[0]/4;  y+=sides[1]/4;  z+=sides[2]/4; break;  
    };  
  
  
}



// If we are allocating from a pool we need to keep a record of it somewhere
// maybe keep track (static?) of holding object with operators
// or define OctCell::new ?
void OctCell::split(size_t* leafc)
{
    if (leaf)
    {
	for (int f=0;f<6;f++)
	{
	    if (leafinfo->next[f] && (depth>leafinfo->next[f]->owner->depth))
	    {
		leafinfo->next[f]->owner->split(leafc);
	    }  
	  
	}
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
	leafinfo->split(kids);
	delete leafinfo;
	leafinfo=0;
	(*leafc)+=7;
    }  
}

void OctCell::allSplit(unsigned int depth, size_t* leafc)
{
    if (depth>0)  
    {
        if (leaf)
	{
	    split(leafc);    
	}
	if (depth==1)
	{
	   return; 
	}
        for (unsigned int i=0;i<8;i++)
	{
	    kids[i]->split(leafc);
	    kids[i]->allSplit(depth-1, leafc);
	}
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
void OctCell::splitPoint(double x, double y, double z, unsigned d, size_t* leafc)
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
	l->split(leafc);
	start=l;
    } while (start->depth+1 <d);
}

// collapse the children of this cell into it
void OctCell::collapse(size_t* leafc)
{
    if (!leaf)
    {
        // collapse any grandkids
        for (unsigned i=0;i<8;++i)
	{
	    if (! kids[i]->leaf)
	    {
                kids[i]->collapse(leafc);	// so we know we are only dealing with leaves
	    }
	}
	// now we need to ensure that collapsing this would not break any neighbours
/*	// we do this by looking at faces from opposite corners

	for (int j=1;j<7;j+=2)
	{

	if (kids[0]->leafinfo->next[j])
	{
	    if (depth+1 <= kids[0]->leafinfo->next[j]->owner->depth)
	    {
	     	if (depth+1<kids[0]->leafinfo->next[j]->owner->depth)
		{
		    kids[0]->leafinfo->next[j]->owner->parent->collapse();	  
		}
		else	// entering else means that there is at least one leaf on the face which is 
		{	// at the correct level but that doesn't mean they all are
		    for (unsigned s=0;s<4;++s)
		    {
		        if (! kids[0]->leafinfo->next[j]->owner->parent->kids[faces[oppdir(j)][s]]->leaf)
			{
			    kids[0]->leafinfo->next[j]->owner->parent->collapse();
			}
		    }
		}
	    }
	}
	
	}
	
	for (int j=0;j<6;j+=2)
	{

	if (kids[6]->leafinfo->next[j])
	{
	    if (depth+1 <= kids[6]->leafinfo->next[j]->owner->depth)
	    {
	     	if (depth+1<kids[6]->leafinfo->next[j]->owner->depth)
		{
		    kids[6]->leafinfo->next[j]->owner->parent->collapse();	  
		}
		else	// entering else means that there is at least one leaf on the face which is 
		{	// at the correct level but that doesn't mean they all are
		    for (unsigned s=0;s<4;++s)
		    {
		        if (! kids[6]->leafinfo->next[j]->owner->parent->kids[faces[oppdir(j)][s]]->leaf)
			{
			    kids[6]->leafinfo->next[j]->owner->parent->collapse();
			}
		    }
		}
	    }
	}
	
	}*/	
	
/*	
	if (kids[0]->leafinfo->next[1] && (depth+1<kids[0]->leafinfo->next[1]->owner->depth))
	{
	    kids[0]->leafinfo->next[1]->owner->parent->collapse();	  
	}
	if (kids[0]->leafinfo->next[3] && (depth+1<kids[0]->leafinfo->next[3]->owner->depth))
	{
	    kids[0]->leafinfo->next[3]->owner->parent->collapse();	  
	}
	if (kids[0]->leafinfo->next[5] && (depth+1<kids[0]->leafinfo->next[5]->owner->depth))
	{
	    kids[0]->leafinfo->next[5]->owner->parent->collapse();	  
	}*/
/*
	if (kids[6]->leafinfo->next[0] && (depth+1<kids[6]->leafinfo->next[0]->owner->depth))
	{
	    kids[6]->leafinfo->next[0]->owner->parent->collapse();	  
	}
	if (kids[6]->leafinfo->next[2] && (depth+1<kids[6]->leafinfo->next[2]->owner->depth))
	{
	    kids[6]->leafinfo->next[2]->owner->parent->collapse();	  
	}
	if (kids[6]->leafinfo->next[4] && (depth+1<kids[6]->leafinfo->next[4]->owner->depth))
	{
	    kids[6]->leafinfo->next[4]->owner->parent->collapse();	  
	}*/
// }


    unsigned newdepth=depth;
    unsigned leafdepth=depth+1;

    for (int k=0;k<6;++k)
    {
	int kid=faces[k][0];
	LeafInfo* neigh=kids[kid]->leafinfo->next[k];
	
	if (!neigh || (neigh->owner->depth < leafdepth))	// no neigbour or leaf is already shallower than our leaf
	{
	    continue;
	}
	// now we need to check to see if there are any subdivided cells on this face
	OctCell* nparent=neigh->owner->parent;
	if (nparent->depth>newdepth)	// this is in case we ended up pointing at a divided cell on this face
	{
	    nparent=nparent->parent;  
	}
	for (unsigned i=0;i<4;++i)
	{
	    if (! nparent->kids[faces[oppdir(k)][i]]->leaf)	// if there are any non-leaves on the neighbour face
	    {
		nparent->kids[faces[oppdir(k)][i]]->collapse(leafc);
	    }      
	}
    }
	kids[0]->leafinfo->merge();	// pick a child and collapse it
	leaf=true;
        for (unsigned i=0;i<8;++i)
	{
	    delete kids[i];
	}
	leaf=true;
	(*leafc)-=7;
    }
}

void OctCell::collapseAll(unsigned int desdepth, size_t* leafc)
{
    if (leaf)
    {
        return;  
    }
    for (unsigned i=0;i<8;++i)
    {
	kids[i]->collapseAll(desdepth, leafc);
    }
    if (depth>=desdepth)
    {
        collapse(leafc);
    }
}


// The leaf which contains (x,y,z) should be at depth at most d
void OctCell::collapsePoint(double x, double y, double z, unsigned d, size_t* leafc)
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
	l->collapse(leafc);
	start=l;
    } while (start->depth+1 >d);
  
}


void OctCell::doLeafWalk_const(const_cellfn c, void* v) const
{
    if (leaf)
    {
        c(*this, v);
    }
    else
    {
        for (unsigned int i=0;i<8;++i)
        {
            kids[i]->doLeafWalk_const(c,v);
        }
    }  
}

void OctCell::doLeafWalk(cellfn c, void* v)
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

void OctCell::doLeafWalkWithKids_const(const_cellfn2 c, int k, void* v) const
{
    if (leaf)
    {
        c(*this, k, v);
    }
    else
    {
        for (unsigned int i=0;i<8;++i)
        {
            kids[i]->doLeafWalkWithKids_const(c,i, v);
        }
    }
  
  
  
}


void OctCell::debug(bool fromroot)
{
    linkCheck(fromroot);  
}


void OctCell::linkCheck(bool fromroot)
{
    if (fromroot && parent)
    {
       parent->debug(false);
    }
    else
    {
      if (leaf)
      {
// 	  if (leafinfo->owner!=this)
// 	  {
// 	      cout << this << " fails check\n";
// 	  }
//cerr << centre[0] << ", " << centre[1] << ", " << centre[2] << "\n";
	  for (unsigned i=0;i<6;++i)
	  {
	      if (leafinfo->next[i])
	      {
//cerr << i << ',';		
		  if (leafinfo->next[i]->owner->depth>10)
		  {
		      cout << this << " fails check ...\n";
		  }
	      }
	  }
//cerr << endl;	  
	  
      }
      else
      {
	  for (unsigned int i=0;i<8;++i)
	  {
	      kids[i]->debug(false);
	  }        
	
      }
    }
  
  
}

bool OctCell::whohas(LeafInfo* li, bool fromroot)
{
    if (fromroot && parent)
    {
       return parent->whohas(li, true);
    }
    else
    {
      bool has=false;
      if (leaf)
      {
	  if (leafinfo==li)
	  {
	      cerr << this << " has it as info\n";
	      has=true;
	  }
	  else
	  {
// 	  if (leafinfo->owner!=this)
// 	  {
// 	      cout << this << " fails check\n";
// 	  }
//cerr << centre[0] << ", " << centre[1] << ", " << centre[2] << "\n";
	  for (unsigned i=0;i<6;++i)
	  {
	      if (leafinfo->next[i]==li)
	      {
		  has=true;
		      cerr << this << " has it as " << i << " face. (d=";
		      cerr << depth << ") ";
		      cerr << centre[0] << ", " << centre[1] << ", " << centre[2] << endl;
	      }
	  }
	  
	  }
//cerr << endl;	  
	  
      }
      else
      {
	  for (unsigned int i=0;i<8;++i)
	  {
	      has|=kids[i]->whohas(li, false);
	  }        
	
      }
      return has;
    }
}

namespace
{
void countleaves(const OctCell& c, void* v)
{
    (*reinterpret_cast<int*>(v))++;  
}  


void writePoints(const OctCell& c, void* v)
{

    ostream& os=*(reinterpret_cast<ostream*>(v));
    static int i=1;
    os << i++ << " " << (c.centre[0]-c.sides[0]/2) << ' ' << (c.centre[1]-c.sides[1]/2) << ' ' << (c.centre[2]-c.sides[2]/2) << '\n';
    os << i++ << " " << (c.centre[0]+c.sides[0]/2) << ' ' << (c.centre[1]-c.sides[1]/2) << ' ' << (c.centre[2]-c.sides[2]/2) << '\n';
    os << i++ << " " << (c.centre[0]+c.sides[0]/2) << ' ' << (c.centre[1]+c.sides[1]/2) << ' ' << (c.centre[2]-c.sides[2]/2) << '\n';
    os << i++ << " " << (c.centre[0]-c.sides[0]/2) << ' ' << (c.centre[1]+c.sides[1]/2) << ' ' << (c.centre[2]-c.sides[2]/2) << '\n';

    os << i++ << " " << (c.centre[0]-c.sides[0]/2) << ' ' << (c.centre[1]-c.sides[1]/2) << ' ' << (c.centre[2]+c.sides[2]/2) << '\n';
    os << i++ << " " << (c.centre[0]+c.sides[0]/2) << ' ' << (c.centre[1]-c.sides[1]/2) << ' ' << (c.centre[2]+c.sides[2]/2) << '\n';
    os << i++ << " " << (c.centre[0]+c.sides[0]/2) << ' ' << (c.centre[1]+c.sides[1]/2) << ' ' << (c.centre[2]+c.sides[2]/2) << '\n';
    os << i++ << " " << (c.centre[0]-c.sides[0]/2) << ' ' << (c.centre[1]+c.sides[1]/2) << ' ' << (c.centre[2]+c.sides[2]/2) << '\n';
}

void writeEl(const OctCell& c, void* v)
{
    static int e=1;
    static int i=1;
    ostream& os=reinterpret_cast<ostream&>(v);
    os << e << " 5 0 ";		// element number
    os << i << ' ' << (i+1) << ' ' << (i+2) << ' ' << (i+3) << ' ';
    os << (i+4) << ' ' << (i+5) << ' ' << (i+6) << ' ' << (i+7) << endl;
    e++;
    i+=8;
}

  
}

void OctCell::gmshDump()
{
    if (parent)
    {
	parent->gmshDump();
    }
    else
    {
      ostream& os=cerr;
	  os << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n";

	  int c=0;
	  doLeafWalk_const(countleaves, &c);
	  os << c*8 << endl;
	  int i=1;
	  doLeafWalk_const(writePoints, &os);
	  os << "$EndNodes\n";
	  os << "$Elements\n";
	  os << c << endl;
	  doLeafWalk_const(writeEl, &os);
	  i=1;
	  os << "$EndElements\n";  
    }
}