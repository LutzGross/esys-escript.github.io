

#include "OctCell.h"

namespace buckley
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
    void walkWithKids(cellfunct2 c, void* v);
    void assignIDs();
    unsigned leafCount();	// inefficient, don't use too much
    
    void debug();
private:
    double side[3];
    OctCell p;
};


}