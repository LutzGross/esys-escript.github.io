

#include "OctCell.h"

namespace refine
{

  
  

class OctTree
{
public:
    OctTree(double x, double y, double z);		// specifies the dimensions of bounding box
    ~OctTree();
    void allSplit(unsigned int d);
    void splitPoint(double x, double y, double z);
    void walkLeaves(cellfunct c, void* v);
private:
    double side[3];
    OctCell p;
};


}