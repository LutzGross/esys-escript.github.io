/*****************************************************************************
*
* Copyright (c) 2014 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/


#include <vector>


/* This file contains two main classes for dealing with a large 3D region which has been divided 
 into a 3D Grid of Blocks (usually to be distributed). Each block is divided into 27 subblocks.
 The first and last subblocks in each dimension are cubes.
 
  class Block mangages a single block. It has methods for copying between a flat array (storing
 all the values in the 3D block) and buffers storing individual subblocks.
 These buffers can be used to transfer individual subblocks to other blocks (using some external 
 means).
   
  class BlockGrid deals with the position of a given block in relation to the rest of the blocks
 in the grid. (In an MPI setting, there would be one block per rank.) 
 It also required transfers of subblocks into and out of the block in order to make the whole
 (global) array consistant.
   
  Each block has a region "inset" wide in from each edge which is shared with nieghbouring blocks.
 Where these regions overlap with another block, the block closest to block 0,0,0 are used.
 Or more precisely, values move left->right then lowy->highy and finally lowz -> highz.

 Please don't mix external calls into this file, it may be useful to separate it for 
 debugging purposes.

 Types required:
      neighbourID_t - Stores the label of a neighbouring block.
		      In an MPI setting, this will be the type used refer to ranks
      coord_t       - Stores a position of a block in the grid of blocks (it could be
		      within one dimension or overall in a flat structure.
		      It is not (necessarily) the same as neighbourID_t because
		      coord_t should be _unsigned_ and there is no guarantee
		      that neighbourID_t will even be an integral type.
*/

#ifdef ESYS_MPI
#include <esysUtils/Esys_MPI.h>
typedef Esys_MPI_rank neighbourID_t;	// This should be the MPI_rank type
typedef unsigned coord_t;		// if we ever get more than 2^32 ranks, we have other problems


#else

typedef int neighbourID_t;	// This should be the MPI_rank type
typedef unsigned coord_t;		// if we ever get more than 2^32 ranks, we have other problems

#endif


typedef std::pair<neighbourID_t, int> neighpair;
typedef std::vector<neighpair> neighbourvector;



typedef struct 
{
public:
  neighbourID_t sourceID;	// ranks involved in communication
  neighbourID_t destID;
  int tag;
  char srcbuffid;		// number of buffer to use for coms
  char destbuffid;
} message;

typedef std::vector<message> messvec;



class BlockGrid
{
public:  
    BlockGrid(coord_t maxx, coord_t maxy, coord_t maxz);

    neighbourID_t getNID(coord_t x, coord_t y, coord_t z) const;


// generate all incoming com messages for this block.
// for each subblock (27 of them), there may be an x, y, z direction to search in  
void generateInNeighbours(coord_t blockx, coord_t blocky, coord_t blockz, messvec& v);


// generate all outgoing com messages for this block
void generateOutNeighbours(coord_t blockx, coord_t blocky, coord_t blockz, messvec& v);
private:
   coord_t xmax;
   coord_t ymax;
   coord_t zmax;   

};

/* Do not ask about buffers for sub-block 1,1,1 (also known as #13)
 They do not exist, such buffers woudl be:
   1) big
   2) unnecessary since the centre sub-block is not sent anywhere

Note that this class does not deal with data transfer between blocks
Sub-blocks are copied to and from buffers. Other code is required to
actually move the data.
*/
class Block
{
public: 
  
    // s? specifiy the [local] size (in points) of each dimension
    Block(size_t sx, size_t sy, size_t sz, size_t inset, size_t xmidlen, size_t ymidlen, size_t zmidlen);
    
    ~Block();    
    
    // Out buffers are loaded with the contents of the flat array and are
    // to be sent to other blocks
    double* getOutBuffer(char subx, char suby, char subz);
    double* getOutBuffer(char bid);

    // In buffers are populated from external communications
    // and copied back to teh flat array
    double* getInBuffer(char subx, char suby, char subz);
    double* getInBuffer(char bid);
    
    // return number of doubles in the given block
    size_t getBuffSize(char subx, char suby, char subz);
    size_t getBuffSize(char bid);
      
    // where does the subblock specified start in a source array
    size_t startOffset(char subx, char suby, char subz);

    // debug only
    void displayBlock(char subx, char suby, char subz, bool out);
    
    // Copy a 3d region from a flat array into a buffer
    void copyToBuffer(char buffid, double* src);

    // Copy a 3d region from a buffer into a flat array
    void copyFromBuffer(char buffid, double* dest);    
    
    
    void copyAllToBuffer(double* src);
    
    void copyUsedFromBuffer(double* dest);

    void setUsed(char buffid);    
  
private:    
  
    // determines the dimensions of each subblock
    void populateDimsTable(); 
    void populateOffsetTable(size_t inset, size_t xmidlen, size_t ymidlen, size_t zmidlen);
    void createBuffArrays(double* startaddress, double* buffptr[27], size_t inset, size_t xmidlen, size_t ymidlen, size_t zmidlen);  
  
  
  
    double* inbuff;
    double* outbuff;
    size_t buffoffsets[27];	// offsets of the various blocks within the buffer arrays
    size_t flatoffsets[27];	// starting point of each block within a flat array
    bool used[27];
    size_t dims[27][3];	// dimension of each subblock 
    size_t sx;
    size_t sy;
    size_t sz;
    size_t inset;
    size_t xmidlen;
    size_t ymidlen;
    size_t zmidlen;    
    double* inbuffptr[27];
    double* outbuffptr[27];


    
};

// Returns the MPI message tag to use for a transfer between the two subblocks
int getTag(char sourcex, char sourcey, char sourcez, char targetx, char targety, char targetz);

// computes the tag based on the destination and the direction it comes from
// the booleans indicate whether a negative shift in that direction is required
int getTag(char destx, char desty, char destz, bool deltax, bool deltay, bool deltaz);


// the booleans indicate whether a negative shift in that direction is required
char getSrcBuffID(char destx, char desty, char destz, bool deltax, bool deltay, bool deltaz);


