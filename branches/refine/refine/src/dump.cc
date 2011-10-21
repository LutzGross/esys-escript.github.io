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


int main()
{
    int c=0;
    cout << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n";
    OctTree ot(1,1,1);
//    ot.allSplit(1);
    cerr << "{{{{{{{{{\n";
//    ot.splitPoint(0.374, 0.375,0.375, 3);
    ot.splitPoint(0.126, 0.126,0.125/2, 4);

    dumpGrid(ot);
    //dumpCen(ot);

}
