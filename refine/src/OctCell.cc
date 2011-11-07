#include <iostream>
#include "OctCell.h"
#include "LeafInfo.h"

using namespace refine;

using namespace std;

namespace 
{
   const double eps=1e-1;    
}



namespace
{

bool gdebug=false;  
  
  
// the direction is:    
// +x -x +y -y +z -z
//  0  1  2  3  4  5
  
// return the opposite direction  (down for up etc)  
inline unsigned short oppdir(unsigned short d)
{
//    return (d%2)?(d-1):(d+1);  
    return d ^ 1;
}

// lists the kids on each face(direction)
static unsigned faces[6][4] = { {1, 2, 6, 5},  // +x
				{3, 0, 4, 7},  // -x
				{3, 2, 6, 7},  // +y
				{0, 1, 5, 4},  // -y
				{4, 5, 6, 7},  // +z
				{3, 2, 1, 0} };// -z


// Which child number touches me in the specified direction
unsigned short neighbour(unsigned short mypos, unsigned short dir)
{
		// direction		0  1  2  3  4  5  
    static unsigned short narr[8][6]={ {1, 1, 3, 3, 4, 4  }, 	// me at pos 0 
				       {0, 0, 2, 2, 5, 5  },	// me at pos 1
				       {3, 3, 1, 1, 6, 6  },	// me at pos 2
				       {2, 2, 0, 0, 7, 7  },	// me at pos 3
				       {5, 5, 7, 7, 0, 0  },	// me at pos 4
				       {4, 4, 6, 6, 1, 1  },	// me at pos 5,
				       {7, 7, 5, 5, 2, 2  },	// me at pos 6
				       {6, 6, 4, 4, 3, 3 }};	// me at pos 7
    return narr[mypos][dir];
}
  
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





// If we are allocating from a pool we need to keep a record of it somewhere
// maybe keep track (static?) of holding object with operators
// or define OctCell::new ?
void OctCell::split()
{
    if (leaf)
    {
      
for (int f=0;f<6;f++)
{
    if (leafinfo->next[f] && (depth>leafinfo->next[f]->owner->depth))
    {
//        cout << this << " depth is " << depth << " " << leafinfo->next[f]->owner << " is " << leafinfo->next[f]->owner->depth << endl;
	leafinfo->next[f]->owner->split();
    }  
  
}

      
/*      
// sanity check on division
for (int f=0;f<6;f++)
{
    if (leafinfo->next[f] && (depth+1<leafinfo->next[f]->owner->depth))
    {
        cerr << "My depth is " << depth << " theirs is " << leafinfo->next[f]->owner->depth << endl;
    }  
  
}*/
      
      
/*debug(true);      
for (int i=0;i<=depth;++i) {cerr << ' ';}
debug(true);
cerr << "Pre " << centre[0] << " " << centre[1] << " " << centre[2] << " d=" << depth<< endl;      
debug(true);*/
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
//debug(true);	
	leafinfo->split(kids);
//debug(true);
/*
if (whohas(leafinfo))
{
cerr << "BADDDDDD\n";  
//gmshDump();  
  
  
}*/

	delete leafinfo;		// someone is still using this leafinfo   arrrgggg
//debug(true);	
	
	leafinfo=0;
// cerr << "---\n";	
// debug(true);
// for (int i=0;i<=depth;++i) {cerr << ' ';}
// cerr << "Post " << centre[0] << " " << centre[1] << " " << centre[2] << endl;	
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
//linkCheck(true);	  
	    kids[i]->split();
//linkCheck(true);
// cerr << "-\n";
	    kids[i]->allSplit(depth-1);
//linkCheck(true);	    
	}
//cerr << "All split at depth " << depth << endl;	
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
        cerr << "Up to " << centre[0] << " " << centre[1] << " " << centre[2] << "\n";
      
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
    // refine neighbours to be at least d-1
    // refine this cell
    OctCell* start=this;
    do
    {
	OctCell* l=start->findLeaf(x, y, z);
	if (l->depth>=d)
	{
	    return;  
	}
	l->split();
	start=l;
    } while (start->depth+1 <d);
}

// collapse the children of this cell into it
void OctCell::collapse()
{
    if (!leaf)
    {
cerr << centre[0] << ", " << centre[1] << ", " << centre[2] << "     ->\n";      
        // collapse any grandkids
        for (unsigned i=0;i<8;++i)
	{
	    if (! kids[i]->leaf)
	    {
                kids[i]->collapse();	// so we know we are only dealing with leaves
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
/*    if (neigh->owner->depth > leafdepth)	// if we collapse then they will be 2 levels down - fix this first
    {
        cerr << "Neighbour's leaf is too deep - collapsing\n";
	neigh->owner->parent->collapse();
    }*/
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
	    nparent->kids[faces[oppdir(k)][i]]->collapse();
	}      
    }
}



// for (int k=0;k<6;++k)
// {
//     int kid=faces[k][0];
//     LeafInfo* neigh=kids[kid]->leafinfo->next[k];
//     
//     if (!neigh || (neigh->owner->depth<depth+1))
//     {
//         continue;
//     }
//     OctCell* nparent=neigh->owner->parent;
//     if (nparent->depth>depth)
//     {
//         nparent=nparent->parent;
//     }
//     for (int i=0;i<4;++i)
//     {
//        if (!kids[faces[oppdir(k)][i]]->leaf)
//        {
// 	   kids[faces[oppdir(k)][i]]->collapse();	 
//        }
//     }
// }



/*


for (int k=0;k<6;++k)
{
    int kid=faces[k][0];
    LeafInfo* neigh=kids[kid]->leafinfo->next[k];
    
    if (!neigh || (neigh->owner->depth<depth+1))
    {
        continue;
    }
    // find the node at the same level as us 
    OctCell* par=neigh->owner->parent;    



    if (parent->depth!=depth)
    {
       par=par->parent; 
    }
    if (par==parent)
    {
        continue;
    }    
    for (int i=0;i<4;++i)	// ensure that all direct children on the face closest to us are leaves
    {
        if (! par->kids[faces[oppdir(k)][i]]->leaf)
	{
            par->kids[faces[oppdir(k)][i]]->collapse();      
	}    
    }
}



// debugging
for (int k=0;k<6;++k)
{
    int kid=faces[k][0];
    LeafInfo* neigh=kids[kid]->leafinfo->next[k];
    
    if (!neigh || (neigh->owner->depth<depth))
    {
        continue;
    }
    
    // now we have a neigbour we can check its faces
    // first get the parent
    OctCell* par=neigh->owner->parent;
    if (parent==par)    {continue;}
    if (neigh->owner->depth>depth+1)
    {
cerr << "Too deep\n";      
        neigh->owner->parent->collapse();  
    }
    else
    {
      for (int i=0;i<4;++i)
      {
	  if (par->kids[faces[oppdir(k)][i]]->depth>depth+1)
	  {
	      par->kids[faces[oppdir(k)][i]]->collapse();
cerr << this << "  On face " << k << " collapsed " << i << endl;	      
	  }
	
      }
    }
}


// debugging
for (int k=0;k<6;++k)
{
    int kid=faces[k][0];
    LeafInfo* neigh=kids[kid]->leafinfo->next[k];
    
    if (!neigh || (neigh->owner->depth<depth))
    {
        continue;
    }
    
    // now we have a neigbour we can check its faces
    // first get the parent
    OctCell* par=neigh->owner->parent;
    if (parent==par)    {continue;}
    if (neigh->owner->depth>depth+1)
    {
        neigh->owner->parent->collapse();  
    }
    else
    {
      for (int i=0;i<4;++i)
      {
	  if (! par->kids[faces[oppdir(k)][i]]->leaf)
	  {
	      par->kids[faces[oppdir(k)][i]]->collapse();
cerr << this << "  On face " << k << " collapsed " << i << endl;	      
	  }
	
      }
    }
}




// debugging
for (int k=0;k<6;++k)
{
    int kid=faces[k][0];
    LeafInfo* neigh=kids[kid]->leafinfo->next[k];
    
    if (!neigh || (neigh->owner->depth<depth))
    {
        continue;
    }
    
    // now we have a neigbour we can check its faces
    // first get the parent
    OctCell* par=neigh->owner->parent;
    if (parent==par)    {continue;}
    if (neigh->owner->depth>depth+1)
    {
        neigh->owner->parent->collapse();  
    }
    else
    {
      for (int i=0;i<4;++i)
      {
	  if (! par->kids[faces[oppdir(k)][i]]->leaf)
	  {
cerr << this << "  On face " << k << " needed " << i << endl;	    
	  }
	
      }
    }
}*/

LeafInfo* t=kids[0]->leafinfo;


	kids[0]->leafinfo->merge();	// pick a child and collapse it
	leaf=true;
whohas(t);  
	
        for (unsigned i=0;i<8;++i)
	{
	    delete kids[i];
	}
	leaf=true;	
      
cerr << centre[0] << ", " << centre[1] << ", " << centre[2] << "     <-\n";      
      
    }
  
}

// void OctCell::collapse()
// {
//     if (!leaf)
//     {
//         // collapse any grandkids
//         for (unsigned i=0;i<8;++i)
// 	{
//             kids[i]->collapse();	// so we know we are only dealing with leaves
// 	}
//         for (unsigned i=0;i<8;++i)
// 	{
// 	    delete kids[i];
// 	}
// 	leaf=true;
//     }  
// }

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
    // refine neighbours to be at least d-1
    // refine this cell
    OctCell* start=this;
    do
    {
	OctCell* l=start->findLeaf(x, y, z);
	if (l->depth<=d)
	{
	    return;  
	}
	l=l->parent;
//	l->outwardCollapse(l->depth+1);		// since we have gone up a level
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
    
    cerr << "Outward collapse: " << centre[0] << " " << centre[1] << " " << centre[2] << "\n";
    cerr << "-" << (centre[0]-sides[0]/2-minside) << ' ' <<  centre[1] << ' ' <<  centre[2]<<"\n";
    upCollPoint(centre[0]-sides[0]/2-minside, centre[1]-minside, centre[2]-minside, desireddepth);
    
    upCollPoint(centre[0]-sides[0]/2-minside, centre[1]-minside, centre[2]+minside, desireddepth);
    
    upCollPoint(centre[0]-sides[0]/2-minside, centre[1]+minside, centre[2]-minside, desireddepth);
    
    upCollPoint(centre[0]-sides[0]/2-minside, centre[1]+minside, centre[2]+minside, desireddepth);
    
    
    cerr << "-\n";
    upCollPoint(centre[0]+sides[0]/2+minside, centre[1]-minside, centre[2]-minside, desireddepth);
    
    upCollPoint(centre[0]+sides[0]/2+minside, centre[1]-minside, centre[2]+minside, desireddepth);

    upCollPoint(centre[0]+sides[0]/2+minside, centre[1]+minside, centre[2]-minside, desireddepth);

    upCollPoint(centre[0]+sides[0]/2+minside, centre[1]+minside, centre[2]+minside, desireddepth);
    
    
    cerr << "-\n";
    upCollPoint(centre[0]-minside, centre[1]-sides[1]/2-minside, centre[2]-minside, desireddepth);
    
    upCollPoint(centre[0]-minside, centre[1]-sides[1]/2-minside, centre[2]+minside, desireddepth);
    
    upCollPoint(centre[0]+minside, centre[1]-sides[1]/2-minside, centre[2]-minside, desireddepth);
    
    upCollPoint(centre[0]+minside, centre[1]-sides[1]/2-minside, centre[2]+minside, desireddepth);
    
    
    
    
    
    cerr << "-\n";
    upCollPoint(centre[0]-minside, centre[1]+sides[1]/2+minside, centre[2]-minside, desireddepth);
    
    upCollPoint(centre[0]-minside, centre[1]+sides[1]/2+minside, centre[2]+minside, desireddepth);

    upCollPoint(centre[0]+minside, centre[1]+sides[1]/2+minside, centre[2]-minside, desireddepth);
    
    upCollPoint(centre[0]+minside, centre[1]+sides[1]/2+minside, centre[2]+minside, desireddepth);
    
    
    cerr << "-\n";
    upCollPoint(centre[0]-minside, centre[1]-minside, centre[2]-sides[2]/2-minside, desireddepth);
    
    upCollPoint(centre[0]-minside, centre[1]+minside, centre[2]-sides[2]/2-minside, desireddepth);
    
    upCollPoint(centre[0]+minside, centre[1]-minside, centre[2]-sides[2]/2-minside, desireddepth);
    
    upCollPoint(centre[0]+minside, centre[1]+minside, centre[2]-sides[2]/2-minside, desireddepth);    
    
    cerr << "-\n";
    upCollPoint(centre[0]-minside, centre[1]-minside, centre[2]+sides[2]/2+minside, desireddepth);

    upCollPoint(centre[0]-minside, centre[1]+minside, centre[2]+sides[2]/2+minside, desireddepth);

    upCollPoint(centre[0]+minside, centre[1]-minside, centre[2]+sides[2]/2+minside, desireddepth);

    upCollPoint(centre[0]+minside, centre[1]+minside, centre[2]+sides[2]/2+minside, desireddepth);
    
    cerr << "----\n";
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
    cout << e << " 5 0 ";		// element number
    cout << i << ' ' << (i+1) << ' ' << (i+2) << ' ' << (i+3) << ' ';
    cout << (i+4) << ' ' << (i+5) << ' ' << (i+6) << ' ' << (i+7) << endl;
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
	  doLeafWalk(countleaves, &c);
	  os << c*8 << endl;
	  int i=1;
	  doLeafWalk(writePoints, &os);
	  os << "$EndNodes\n";
	  os << "$Elements\n";
	  os << c << endl;
	  doLeafWalk(writeEl, &os);
	  i=1;
	  os << "$EndElements\n";  
    }
}