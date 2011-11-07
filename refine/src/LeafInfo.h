
#ifndef LEAFINFO_H_2011
#define LEAFINFO_H_2011

#include "OctCell.h"

namespace refine
{

class LeafInfo
{
public:
    LeafInfo(OctCell* c);
    ~LeafInfo();
    void split(OctCell* kids[8]);
    void merge();
    OctCell* owner;
    LeafInfo* next[6];
    unkid pmap[8];
    
    int z[111];		// padding to reduce odds of another struct landing on top of this
    
};

}

#endif