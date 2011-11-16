
#ifndef OCTCELL_H_2011
#define OCTCELL_H_2011
namespace buckley
{

class OctCell;  
class LeafInfo;
  
typedef void (*cellfunct)(const OctCell&, void*);  
typedef void (*cellfunct2)(const OctCell&, int k, void*);
typedef unsigned int unkid;			// Type representing the id of an unknown/point



class OctCell
{
public:
      OctCell(double x1, double y1, double z1, double x2, double y2, double z2, OctCell* par);
      ~OctCell();
      void split();		// split this cell into 8 children
      void collapse();		// remove all kids and make this a leaf
      void collapseAll(unsigned desdepth);
      void collapsePoint(double x, double y, double z, unsigned d);
      void allSplit(unsigned int depth);
      void splitPoint(double x, double y, double z, unsigned desdepth);
      void doLeafWalk(cellfunct c, void* v);
      void doLeafWalkWithKids(cellfunct2 c, int k, void* v);
      OctCell* findLeaf(double x, double y, double z);
      
      
      void linkCheck(bool fromroot);
      
      
      void debug(bool fromroot);
      bool whohas(LeafInfo* li, bool fromroot=true);
      
      void gmshDump();
//private:
      void upSplitPoint(double x, double y, double z, unsigned d);
      void upCollPoint(double x, double y, double z, unsigned d);      
      void outwardRefine(unsigned desireddepth);
      void outwardCollapse(unsigned desireddepth);
      
      bool leaf;
      double centre[3];
      double sides[3];		// dimensions in x,y,z
      OctCell* kids[8];
      unsigned int depth;
      OctCell* parent;
      unsigned int id;
      LeafInfo* leafinfo;
};

}

#include "LeafInfo.h"

#endif
