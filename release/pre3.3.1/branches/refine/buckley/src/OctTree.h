

#include "OctCell.h"
#include <vector>

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
    void walkLeaves(const_cellfn c, void* v) const;
    void walkLeaves(cellfn c, void* v);
    void walkWithKids(const_cellfn2 c, void* v);
    unkid assignIDs() const;
    unsigned leafCountByWalk() const;	// inefficient, don't use too much
    unsigned leafCount() const;
    
    const OctCell** process(unkid& numunk,std::vector<const OctCell*> face_cells[6]) const;	// tree mutations are complete for now
    
    
    void getBB(double bb[6]) const;	// get the bounding box of the tree
    void debug();
private:

    double side[3];
    OctCell p;
    bool modified;
    mutable size_t leafcount;
    friend class OctCell;
};


}