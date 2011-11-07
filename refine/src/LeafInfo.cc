#include <cstring>	// for memset
#include <cmath>
#include <iostream>
#include "LeafInfo.h"

using namespace refine;
using namespace std;

LeafInfo::LeafInfo(refine::OctCell* c)
{
    owner=c;
    memset(next, 0, sizeof(OctCell*)*6);
    memset(pmap, 0, sizeof(unkid)*8);
}


// The assumption here is that a LeafInfo will always be managed by an OctCell which will deal with any rolling deletes

LeafInfo::~LeafInfo()
{
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

// The order of the next pointers is:
// +x -x +y -y +z -z
//  0  1  2  3  4  5

// The owner of this leafinfo has been split into the cells given in the kids array
// Patch all the relevant pointers.
// kids are created as:   3 2    7 6
//                        0 1    4 5
void LeafInfo::split(OctCell* kids[8])
{

  
// if (fabs(owner->centre[0]-0.75)+fabs(owner->centre[1]-0.25)+fabs(owner->centre[2]-0.25)<0.01)
// {
//   cerr << "******\n";
//   
//   // do a sanity check on all neighbours
//   cerr << "Owners depth=" << owner->depth << endl;
//   
//   for (int i=0;i<6;++i)
//   {
//      if (next[i])
//      {
// 	cerr << i << " d=" << next[i]->owner->depth << endl;     
//      }
//      else
//      {
//         cerr << i << " null.\n";
//      }
//   }
//   
//   
// //gdebug=true;  
// }
// else
// {
// //gdebug=false;  
// }
	// now we form the intra split links
	kids[0]->leafinfo->next[0]=kids[1]->leafinfo;
	kids[0]->leafinfo->next[2]=kids[3]->leafinfo;
	kids[0]->leafinfo->next[4]=kids[4]->leafinfo;
	
	kids[4]->leafinfo->next[0]=kids[5]->leafinfo;
	kids[4]->leafinfo->next[2]=kids[7]->leafinfo;
	kids[4]->leafinfo->next[5]=kids[0]->leafinfo;
	
	kids[1]->leafinfo->next[1]=kids[0]->leafinfo;
	kids[1]->leafinfo->next[2]=kids[2]->leafinfo;
	kids[1]->leafinfo->next[4]=kids[5]->leafinfo;
	
	kids[5]->leafinfo->next[1]=kids[4]->leafinfo;
	kids[5]->leafinfo->next[2]=kids[6]->leafinfo;
	kids[5]->leafinfo->next[5]=kids[1]->leafinfo;
	
	kids[2]->leafinfo->next[1]=kids[3]->leafinfo;
	kids[2]->leafinfo->next[3]=kids[1]->leafinfo;
	kids[2]->leafinfo->next[4]=kids[6]->leafinfo;
	
	kids[6]->leafinfo->next[1]=kids[7]->leafinfo;
	kids[6]->leafinfo->next[3]=kids[5]->leafinfo;
	kids[6]->leafinfo->next[5]=kids[2]->leafinfo;
	
	kids[3]->leafinfo->next[0]=kids[2]->leafinfo;
	kids[3]->leafinfo->next[3]=kids[0]->leafinfo;
	kids[3]->leafinfo->next[4]=kids[7]->leafinfo;
	
	kids[7]->leafinfo->next[0]=kids[6]->leafinfo;
	kids[7]->leafinfo->next[3]=kids[4]->leafinfo;
	kids[7]->leafinfo->next[5]=kids[3]->leafinfo;      
  
	
	for (unsigned short f=0;f<6;++f)
	{
if (gdebug) {cout << "face " << f << endl;}	  
	    if (next[f])	// if we have neighbours in that direction
	    {
	        if (next[f]->owner->depth > owner->depth)	// if this neighbour is already split further than us
		{
		    OctCell* p=next[f]->owner->parent;
		    // walk each of the new children and link them to their neighbours on the same level
		    // this is just an unrolled for loop
		    kids[faces[f][0]]->leafinfo->next[f] = p->kids[neighbour(faces[f][0], f)]->leafinfo;
		    kids[faces[f][1]]->leafinfo->next[f] = p->kids[neighbour(faces[f][1], f)]->leafinfo;
		    kids[faces[f][2]]->leafinfo->next[f] = p->kids[neighbour(faces[f][2], f)]->leafinfo;
		    kids[faces[f][3]]->leafinfo->next[f] = p->kids[neighbour(faces[f][3], f)]->leafinfo;

if (gdebug) 
{
  cout << "neigbour is " << next[f] << endl;
  cout << faces[f][1] << " " << neighbour(1, f) << endl;
  cout << neighbour(faces[f][1], f) << endl;
/*  cout << ' ' << faces[f][0] << " -> " << neighbour(faces[f][0], f) << " " << oppdir(f) << endl;
  cout << ' ' << faces[f][1] << " -> " << neighbour(faces[f][1], f) << " " << oppdir(f) << endl;
  cout << ' ' << faces[f][2] << " -> " << neighbour(faces[f][2], f) << " " << oppdir(f) << endl;
  cout << ' ' << faces[f][3] << " -> " << neighbour(faces[f][3], f) << " " << oppdir(f) << endl;*/
}		    
		    
		    // now the incoming links from neighbours on this face
		    p->kids[neighbour(faces[f][0], f)]->leafinfo->next[oppdir(f)]=kids[faces[f][0]]->leafinfo;
		    p->kids[neighbour(faces[f][1], f)]->leafinfo->next[oppdir(f)]=kids[faces[f][1]]->leafinfo;
		    p->kids[neighbour(faces[f][2], f)]->leafinfo->next[oppdir(f)]=kids[faces[f][2]]->leafinfo;
		    p->kids[neighbour(faces[f][3], f)]->leafinfo->next[oppdir(f)]=kids[faces[f][3]]->leafinfo;

		  
		}
		else	// before we split, the neighbour has the same depth as us
		{
		    OctCell* p=next[f]->owner;	// all links point to this single leaf
		    
		    kids[faces[f][0]]->leafinfo->next[f] = p->leafinfo;
		    kids[faces[f][1]]->leafinfo->next[f] = p->leafinfo;
		    kids[faces[f][2]]->leafinfo->next[f] = p->leafinfo;
		    kids[faces[f][3]]->leafinfo->next[f] = p->leafinfo;
		    
		    p->leafinfo->next[oppdir(f)]=kids[faces[f][0]]->leafinfo;	// any child on the face will do
		  
		  
		}
	      
	    }
	  
	  
	}
}

// collapse this leaf (and all its siblings) into a single leaf
// Warning: When this operation is complete, there will be no outside links (from other leaves) pointing to any
// of the siblings but memory management is still the responsibility of the owner
void LeafInfo::merge()
{
    OctCell* pcell=owner->parent;
    if (pcell==0)
    {
       return;
    }
// cout << "Merging into parent " << pcell << endl;    
    LeafInfo* li=new LeafInfo(pcell);
    // fill in faces for new leafinfo
    // pick a cell on each face and access its neighbour field
    for (unsigned short f=0;f<6;++f)
    {
        // pick a child on this face
        unsigned short s=faces[f][0];
	
	LeafInfo* n=pcell->kids[s]->leafinfo->next[f];
	if (n)		// if we have neighbours in that direction
	{
	    if (n->owner->depth < owner->depth)		// if the neighbour is already shallower than us (prior to merge)
	    {
	        li->next[f]=n;
		// now to patch the return link
		n->next[oppdir(f)]=li;
	    }
	    else		// neighbour is on the same level as us
	    {
	        li->next[f]=n;
	        // need to patch all the return links on this face 
//		n->owner->parent->kids[oppdir(f)][0].leafinfo->next[oppdir(f)]=li;
		n->owner->parent->kids[faces[oppdir(f)][0]]->leafinfo->next[oppdir(f)]=li;
		n->owner->parent->kids[faces[oppdir(f)][1]]->leafinfo->next[oppdir(f)]=li;	
		n->owner->parent->kids[faces[oppdir(f)][2]]->leafinfo->next[oppdir(f)]=li;	
		n->owner->parent->kids[faces[oppdir(f)][3]]->leafinfo->next[oppdir(f)]=li;			

	    }	  
	}

    }
    pcell->leafinfo=li;
//    return li;
}