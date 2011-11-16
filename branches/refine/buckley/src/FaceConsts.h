#ifndef FACECONSTS_H
#define FACECONSTS_H

namespace
{

// This is to be included in an anonymous namespace in a .cc file

  
// the direction is:    
// +x -x +y -y +z -z
//  0  1  2  3  4  5
  
// return the opposite direction  (down for up etc)  
inline unsigned short oppdir(unsigned short d)
{
//    return (d%2)?(d-1):(d+1);  
    return d ^ 1;
}

// lists the kids on each face(direction)
static unsigned faces[6][4] = { {1, 2, 6, 5},  // +x
				{3, 0, 4, 7},  // -x
				{3, 2, 6, 7},  // +y
				{0, 1, 5, 4},  // -y
				{4, 5, 6, 7},  // +z
				{3, 2, 1, 0} };// -z


// Which child number touches me in the specified direction
unsigned short neighbour(unsigned short mypos, unsigned short dir)
{
		// direction		0  1  2  3  4  5  
    static unsigned short narr[8][6]={ {1, 1, 3, 3, 4, 4  }, 	// me at pos 0 
				       {0, 0, 2, 2, 5, 5  },	// me at pos 1
				       {3, 3, 1, 1, 6, 6  },	// me at pos 2
				       {2, 2, 0, 0, 7, 7  },	// me at pos 3
				       {5, 5, 7, 7, 0, 0  },	// me at pos 4
				       {4, 4, 6, 6, 1, 1  },	// me at pos 5,
				       {7, 7, 5, 5, 2, 2  },	// me at pos 6
				       {6, 6, 4, 4, 3, 3 }};	// me at pos 7
    return narr[mypos][dir];
}


unsigned short facestouch[8][3] = {{1, 3, 5},	 // which faces touch the specified corner
				   {0, 3, 5},
				   {0, 2, 5},  
				   {1, 2, 5},
				   {1, 3, 4},
				   {0, 3, 4},
				   {0, 2, 4},
				   {1, 2, 4} };

// canhang[c][n] is true if a leaf which is the c'th child of its parent and its n'th node could hang
bool canhang[8][8] = {{false, true, true, true, true, true, false, true},	// 0
		      {true, false, true, true, true, true, true, false},	// 1
		      {true, true, false, true, false, true, true, true},	// 2
		      {true, true, true, false, true, false, true, true},	// 3
		      {true, true, false, true, false, true, true, true},	// 4
		      {true, true, true, false, true, false, true, true},	// 5
		      {false, true, true, true, true, true, false, true},	// 6
		      {true, false, true, true, true, true, true, false}	// 7
		      };
}	// end of namespace

// This really should be typed as unkid --- will fix that once I've untangled what goes where
const unsigned int HANG_NODE=1;

#endif