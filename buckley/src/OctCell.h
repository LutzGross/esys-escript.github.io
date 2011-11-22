
#ifndef OCTCELL_H_2011
#define OCTCELL_H_2011

#include <cstddef>	// to get size_t

namespace buckley
{

class OctCell;  
class LeafInfo;

typedef void (*const_cellfn)(const OctCell&, void*); 

typedef void (*cellfn)(OctCell&, void*);  
typedef void (*cellfn2)(const OctCell&, int k, void*);
typedef void (*const_cellfn2)(const OctCell&, int k, void*);

typedef unsigned int unkid;			// Type representing the id of an unknown/point

class OctCell
{
public:
      OctCell(double x1, double y1, double z1, double x2, double y2, double z2, OctCell* par);
      ~OctCell();
      void split(size_t* leafc);		// split this cell into 8 children
      void collapse(size_t* leafc);		// remove all kids and make this a leaf
      void collapseAll(unsigned desdepth, size_t* leafc);
      void collapsePoint(double x, double y, double z, unsigned desdepth, size_t* leafc);
      void allSplit(unsigned int depth, size_t* leafc);
      void splitPoint(double x, double y, double z, unsigned desdepth, size_t* leafc);
      void doLeafWalk_const(const_cellfn c, void* v) const;
      void doLeafWalk(cellfn c, void* v);      
      void doLeafWalkWithKids_const(const_cellfn2 c, int k, void* v) const;
      OctCell* findLeaf(double x, double y, double z);
      
      void childCoords(unsigned k, double& x, double& y, double& z) const;
      void quadCoords(unsigned k, double& x, double& y, double& z) const;
      
      
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
      mutable unsigned int id;
      LeafInfo* leafinfo;
};

}

#include "LeafInfo.h"

#endif
