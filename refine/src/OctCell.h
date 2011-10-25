
namespace refine
{

class OctCell;  
  
typedef void (*cellfunct)(const OctCell&, void*);  



class OctCell
{
public:
      OctCell(double x1, double y1, double z1, double x2, double y2, double z2, OctCell* par);
      ~OctCell();
      void split();		// split this cell into 8 children
      void collapse();		// remove all kids and make this a leaf
      void allSplit(unsigned int depth);
      void splitPoint(double x, double y, double z, unsigned desdepth);
      void merge();
      void doLeafWalk(cellfunct c, void* v);  
      void outwardRefine(unsigned desireddepth);
      OctCell* findLeaf(double x, double y, double z);
//private:
      void upSplitPoint(double x, double y, double z, unsigned d);
      bool leaf;
      double centre[3];
      double sides[3];		// dimensions in x,y,z
      OctCell* kids[8];
      unsigned int depth;
      OctCell* parent;
};

}