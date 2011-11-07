all:	dump.cc OctTree.cc OctCell.cc LeafInfo.h LeafInfo.cc
	g++ -g dump.cc OctTree.cc OctCell.cc LeafInfo.cc

dom:	domtest.cc OctTree.cc OctTree.h OctCell.cc OctCell.h
	g++ -g domtest.cc OctTree.cc OctTree.h OctCell.cc OctCell.h

test2:	test2.cc OctTree.cc OctTree.h OctCell.cc OctCell.h LeafInfo.h LeafInfo.cc
	g++ -g test2.cc OctTree.cc OctCell.cc LeafInfo.cc   -o test2 

