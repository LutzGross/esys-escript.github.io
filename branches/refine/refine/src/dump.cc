#include <iostream>
#include "OctTree.h"

using namespace std;
using namespace refine;

void f(const OctCell& c, void* v)
{

  
    static int i=1;
    cout << i++ << " " << (c.centre[0]-c.sides[0]/2) << ' ' << (c.centre[1]-c.sides[1]/2) << ' ' << (c.centre[2]-c.sides[2]/2) << '\n';
    cout << i++ << " " << (c.centre[0]+c.sides[0]/2) << ' ' << (c.centre[1]-c.sides[1]/2) << ' ' << (c.centre[2]-c.sides[2]/2) << '\n';
    cout << i++ << " " << (c.centre[0]+c.sides[0]/2) << ' ' << (c.centre[1]+c.sides[1]/2) << ' ' << (c.centre[2]-c.sides[2]/2) << '\n';
    cout << i++ << " " << (c.centre[0]-c.sides[0]/2) << ' ' << (c.centre[1]+c.sides[1]/2) << ' ' << (c.centre[2]-c.sides[2]/2) << '\n';

    cout << i++ << " " << (c.centre[0]-c.sides[0]/2) << ' ' << (c.centre[1]-c.sides[1]/2) << ' ' << (c.centre[2]+c.sides[2]/2) << '\n';
    cout << i++ << " " << (c.centre[0]+c.sides[0]/2) << ' ' << (c.centre[1]-c.sides[1]/2) << ' ' << (c.centre[2]+c.sides[2]/2) << '\n';
    cout << i++ << " " << (c.centre[0]+c.sides[0]/2) << ' ' << (c.centre[1]+c.sides[1]/2) << ' ' << (c.centre[2]+c.sides[2]/2) << '\n';
    cout << i++ << " " << (c.centre[0]-c.sides[0]/2) << ' ' << (c.centre[1]+c.sides[1]/2) << ' ' << (c.centre[2]+c.sides[2]/2) << '\n';
}

void cen_pts(const OctCell& c, void* v)
{
    int& i=*(reinterpret_cast<int*>(v));  
    cout << i++ << ' ' << c.centre[0] << ' ' << c.centre[1] << ' ' << c.centre[2] << endl;
}

void cen_elts(const OctCell& c, void* v)
{
    int& i=*(reinterpret_cast<int*>(v));  
    cout << i << " 15 3 0 "<<i <<" 0 " << i << endl;		// 15 is the code for a single point element
    i++;
}

void g(const OctCell& c, void* v)
{
    static int e=1;
    static int i=1;
    
    cout << e << " 5 0 ";		// element number
    cout << i << ' ' << (i+1) << ' ' << (i+2) << ' ' << (i+3) << ' ';
    cout << (i+4) << ' ' << (i+5) << ' ' << (i+6) << ' ' << (i+7) << endl;
    e++;
    i+=8;
}

void countleaves(const OctCell& c, void* v)
{
    (*reinterpret_cast<int*>(v))++;  
}

void dumpGrid(OctTree& ot)
{
      cout << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n";
    int c=0;
    ot.walkLeaves(countleaves, &c);
    cout << c*8 << endl;
    int i=1;
    ot.walkLeaves(f,0);
    cout << "$EndNodes\n";
    cout << "$Elements\n";
    cout << c << endl;
    ot.walkLeaves(g, 0);
    i=1;
    cout << "$EndElements\n";  
  
}

void dumpCen(OctTree& ot)
{
    int c=0;
    ot.walkLeaves(countleaves, &c);
    cout << c << endl;
    int i=1;
    ot.walkLeaves(cen_pts, &i); 
    cout << "$EndNodes\n";
    cout << "$Elements\n";
    cout << c << endl;
    i=1;
    ot.walkLeaves(cen_elts, &i);
    cout << "$EndElements\n";  
  
}


void printLInfo(const OctCell& c, void* v)
{
    if (c.leaf)
    {
	cout << c.leafinfo << " belongs to " << &c << " " << c.centre[0] << "," << c.centre[1] << "," << c.centre[2] << "," << endl;
    }  
}

void neigh(const OctCell& c, void* v)
{
    if (c.leaf)
    {
        cout << &c << endl;
        for (int i=0;i<6;++i)
	{
	    cout << "  " << i << ": ";
	    if (c.leafinfo->next[i])
	    {
	        cout << c.leafinfo->next[i];
		if (c.leafinfo->next[i]->owner)
		{
		  cout  << " at depth " << c.leafinfo->next[i]->owner->depth <<endl;	  
		}
		else
		{
		  cout << " NULL owner\n";
		}
	    }
	    else
	    {
	        cout << "null\n";
	    }
	}
    }
}


int main()
{
    int c=0;
//    cout << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n";
    OctTree ot(1,1,1);
    ot.allSplit(7);

    
//    ot.allCollapse(7);
//    ot.collapsePoint(0.120, 0.120, 0.120,1);

//      ot.collapsePoint(0.5, 0.5, 0.5, 1);

//      ot.collapsePoint(0.51, 0.51, 0.51, 1);
//      ot.collapsePoint(0.49, 0.49, 0.49, 2);

//return 0;   
//      ot.debug();
      

      
      ot.collapsePoint(0.49, 0.49, 0.49, 2);


//    dumpGrid(ot);      
//return 0;      
      
//     ot.walkLeaves(printLInfo, 0);
//     cout << "----------\n";
//     ot.walkLeaves(neigh, 0);
//     cout << "\n\n";  

   
      
      
      ot.collapsePoint(0.1, 0.49, 0.49, 3);
      
      
//    ot.walkLeaves(printLInfo, 0);
/*    cout << "----------\n";
    ot.walkLeaves(neigh, 0);
    cout << "\n\n";      */
      
      
      ot.collapsePoint(0.3, 0.49, 0.49, 3);
      ot.collapsePoint(0.7, 0.49, 0.49, 3);
      ot.collapsePoint(0.99, 0.49, 0.49,3);
      
      ot.collapsePoint(0.1, 0.99, 0.99,3);
      ot.collapsePoint(0.49, 0.99, 0.99,3);
      ot.collapsePoint(0.99, 0.99, 0.99,3);
      

//    ot.splitPoint(0.512, 0.126,0.99, 10);
    
//    ot.splitPoint(0.2, 0.8,0.01, 30);
     
     ot.assignIDs();
    dumpGrid(ot);
    //dumpCen(ot);

}
