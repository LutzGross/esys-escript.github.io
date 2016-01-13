

#include "OctCell.h"

namespace buckey
{

  
  

class OctTree
{
public:
    OctTree(double x, double y, double z);		// specifies the dimensions of bounding box
    ~OctTree();
    void allSplit(unsigned int d);
    void allCollapse(unsigned int d);
    void collapsePoint(double x, double y, double z, unsigned int d)    ;
    void splitPoint(double x, double y, double z, unsigned desdepth);
    void walkLeaves(cellfunct c, void* v);
    void assignIDs();
    unsigned leafCount();	// inefficient, don't use too much
private:
    double side[3];
    OctCell p;
};


}