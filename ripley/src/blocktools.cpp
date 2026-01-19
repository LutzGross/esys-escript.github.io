
/*****************************************************************************
*
* Copyright (c) 2014-2018 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#include "blocktools.h"

#include <cstring>	// for memset
#include <iostream>	// for the debug method

using namespace std;

BlockGrid::BlockGrid(coord_t x, coord_t y, coord_t z)
 : xmax(x), ymax(y), zmax(z)
{}

neighbourID_t BlockGrid::getNID(coord_t x, coord_t y, coord_t z) const
{
    return x+y*(xmax+1)+z*(xmax+1)*(ymax+1);    
}

// generate all incoming com messages for this block.
// for each subblock (27 of them), there may be an x, y, z direction to search in  
void generateInNeighbours(coord_t blockx, coord_t blocky, coord_t blockz, messvec& v);


// generate all outgoing com messages for this block
void generateOutNeighbours(coord_t blockx, coord_t blocky, coord_t blockz, messvec& v);

// generate all incoming com messages for this block.
// for each subblock (27 of them), there may be an x, y, z direction to search in  
void BlockGrid::generateInNeighbours(coord_t blockx, coord_t blocky, coord_t blockz, messvec& v)
{
    neighbourID_t myid=getNID(blockx, blocky, blockz);
    unsigned char deltax=0;
    unsigned char deltay=0;
    unsigned char deltaz=blockz?1:0;	// do not look to a lower layer if we are on the bottom layer of _blocks_.
    for (unsigned char z=0;z<3;++z)
    {
      deltay=blocky?1:0;		// do not look at a lower row of blocks if there isn't one
      for (unsigned char y=0;y<3;++y)
      {
	deltax=blockx?1:0;
	for (unsigned char x=0;x<3;++x)
	{
	  // now we will have a set of possible directions to look in (delta?==0 means don't attempt to
	  // change this component	  
	  if (deltax+deltay+deltaz)	// if we have an import to do
	  {
	      coord_t srcx=blockx-deltax;
	      coord_t srcy=blocky-deltay;
	      coord_t srcz=blockz-deltaz;
	      message m;
	      m.sourceID=getNID(srcx, srcy, srcz);
	      m.destID=myid;
	      m.tag=getTag(x,y,z,deltax==1, deltay==1, deltaz==1);
	      m.srcbuffid=getSrcBuffID(x,y,z,deltax==1, deltay==1, deltaz==1);
	      m.destbuffid=x+y*3+z*9;
	      v.push_back(m);
	  }
	  deltax=0;	// we are no longer on the left hand edge
	}
	deltay=0;	// only y=0 imports
      }
      deltaz=0;	// since only the bottom sublayer looks to the row below
    }  
}


// generate all outgoing com messages for this block
void BlockGrid::generateOutNeighbours(coord_t blockx, coord_t blocky, coord_t blockz, messvec& v)
{
    messvec vv;
    neighbourID_t myid=getNID(blockx, blocky, blockz);
    for (unsigned char z=0;z<2;++z)
    {
        if (z && (blockz==zmax))	// don't look up if there is no up
	{
	    break;
	}
        for (unsigned char y=0;y<2;++y)
	{
	    if (y && (blocky==ymax))
	    {
	        break;
	    }
	    for (unsigned char x=0;x<2;++x)
	    {
		if (x && (blockx==xmax))
		{
		    break;
		}
		if (x+y+z>0)	// don't look at your own neighbours
		{
		    generateInNeighbours(blockx+x, blocky+y, blockz+z, vv);
		}
	    }
	}
    }
    // Now we reverse the direction of the coms (since we want Out not In
    for (int i=0;i<vv.size();++i)
    {
	if (vv[i].sourceID==myid)
	{
	    v.push_back(vv[i]);
	}    
    }  
}

  
double* Block::getOutBuffer(unsigned char subx, unsigned char suby, unsigned char subz)
{
    unsigned char bid=subx+suby*3+subz*3*3;	// (bid is "blockid")
    if (bid==13)	// there is no buffer for block 13
    {
	return 0;	// don't ask for this buffer because refusal may offend
    }
    return outbuffptr[bid];	
}

double* Block::getInBuffer(unsigned char subx, unsigned char suby, unsigned char subz)
{
    unsigned char bid=subx+suby*3+subz*3*3;	// (bid is "blockid")
    if (bid==13)	// there is no buffer for block 13
    {
	return 0;	// don't ask for this buffer because refusal may offend
    }
    return inbuffptr[bid];
}

size_t Block::getBuffSize(unsigned char subx, unsigned char suby, unsigned char subz)
{
    unsigned char bid=subx+suby*3+subz*3*3;	// (bid is "blockid")
    if (bid==13)	// there is no buffer for block 13
    {
	return 0;	
    }
    return dims[bid][0]*dims[bid][1]*dims[bid][2]*dpsize;	
}

double* Block::getOutBuffer(unsigned char bid)
{
    if (bid==13)	// there is no buffer for block 13
    {
	return 0;	// don't ask for this buffer because refusal may offend
    }
    return outbuffptr[bid];	
}

double* Block::getInBuffer(unsigned char bid)
{
    if (bid==13)	// there is no buffer for block 13
    {
	return 0;	// don't ask for this buffer because refusal may offend
    }
    return inbuffptr[bid];
}

size_t Block::getBuffSize(unsigned char bid)
{
    if (bid==13)	// there is no buffer for block 13
    {
	return 0;	
    }
    return dims[bid][0]*dims[bid][1]*dims[bid][2]*dpsize;	
}



Block::~Block()
{
    delete[] inbuff;
    delete[] outbuff;
}

void Block::populateDimsTable()
{
    for (int i=0;i<27;++i)
    {
	for (int j=0;j<3;++j)
	{
	    dims[i][j]=inset;
	}
    }
    for (int i=1;i<27;i+=3)
    {
	dims[i][0]=xmidlen;
    }
    for (int l=0;l<3;++l)
    {
	dims[3+l*9][1]=ymidlen;
	dims[4+l*9][1]=ymidlen;
	dims[5+l*9][1]=ymidlen;
    }
    for (int i=9;i<18;++i)
    {
	dims[i][2]=zmidlen;
    }
}


// gives the relative start location within a flat 
void Block::populateOffsetTable(size_t inset, size_t xmidlen, size_t ymidlen, size_t zmidlen)
{
  
    size_t cur=0;
    for (int i=0;i<27;++i)
    {
	flatoffsets[i]=cur;
	cur+=dims[i][0]*dims[i][1]*dims[i][2]*dpsize;
    }
    for (int i=0;i<13;++i)
    {
	buffoffsets[i]=flatoffsets[i];
    }
    buffoffsets[13]=0;
    for (int i=14;i<27;++i)
    {
	buffoffsets[i]=flatoffsets[i]-flatoffsets[14]+flatoffsets[13];
    }

}


void Block::createBuffArrays(double* startaddress, double* buffptr[27], size_t inset, size_t xmidlen, size_t ymidlen, size_t zmidlen)
{
    buffptr[0]=startaddress;
    for (int i=0;i<27;++i)
    {
	buffptr[i]=startaddress+buffoffsets[i];
    }
    buffptr[13]=0;	// since the middle should never be sent anywhere
}

void Block::setUsed(unsigned char buffid)
{
    used[buffid]=true;
}

void Block::copyAllToBuffer(double* src)
{
    for (unsigned char i=0;i<13;++i)
    {
	copyToBuffer(i, src);
      
    }
    for (unsigned char i=14;i<27;++i)
    {
	copyToBuffer(i, src);
    }
}
    
void Block::copyUsedFromBuffer(double* dest)
{
    for (unsigned char i=0;i<27;++i)
    {
	if (used[i])
	{
	    copyFromBuffer(i, dest);
	}
    }
}

// s? specifiy the size (in points) of each dimension
// maxb? gives the largest block number in each dimension in the overall grid (number from zero)
Block::Block(size_t sx, size_t sy, size_t sz, size_t inset, size_t xmidlen, 
	     size_t ymidlen, size_t zmidlen, unsigned int dpp) : dpsize(dpp)
{
    this->sx=sx;
    this->sy=sy;
    this->sz=sz;
    this->inset=inset;
    this->xmidlen=xmidlen;
    this->ymidlen=ymidlen;
    this->zmidlen=zmidlen;
    populateDimsTable();	
    
    size_t totalbuff=0;
    for (int i=0;i<27;++i)
    {
	used[i]=false;
	totalbuff+=dims[i][0]*dims[i][1]*dims[i][2];
    }
    totalbuff-=dims[13][0]*dims[13][1]*dims[13][2];
	
    totalbuff*=dpsize;		// to account for non-scalars
    inbuff=new double[totalbuff];
    outbuff=new double[totalbuff];
    memset(inbuff,0,totalbuff*sizeof(double));
    memset(outbuff,0,totalbuff*sizeof(double));
    populateOffsetTable(inset, xmidlen, ymidlen, zmidlen);
    createBuffArrays(inbuff, inbuffptr, inset, xmidlen, ymidlen, zmidlen);
    createBuffArrays(outbuff, outbuffptr, inset, xmidlen, ymidlen, zmidlen);

}  

// where does the subblock specified start in a source array
size_t Block::startOffset(unsigned char subx, unsigned char suby, unsigned char subz)
{
    size_t off=0;
    off+=((subx==0) ? 0 :((subx==1)?inset:(inset+xmidlen) ));
    // now the y component
    size_t ystep=((suby==0) ? 0 :((suby==1)?inset:(inset+ymidlen) ));
    size_t zstep=((subz==0) ? 0 :((subz==1)?inset:(inset+zmidlen) ));
    off+=ystep*(2*inset+xmidlen);
    off+=zstep*(2*inset+xmidlen)*(2*inset+ymidlen);
    return off*dpsize;
}


// This is only intended for debugging, so I have not made it fast
void Block::displayBlock(unsigned char subx, unsigned char suby, unsigned char subz, bool out)
{
    
    unsigned char bid=subx+suby*3+subz*3*3;	// (bid is "blockid")	
    double* b=out?outbuffptr[bid]:inbuffptr[bid];
    for (int z=0;z<dims[bid][2];++z)
    {
	std::cout << std::endl << "Row " << z << std::endl;
      
	for (int y=0;y<dims[bid][1];++y)
	{
	    for (int x=0;x<dims[bid][0];++x)
	    {
		if (dpsize==1)
		{
		    std::cout << b[x+y*dims[bid][0]+z*dims[bid][1]*dims[bid][0]] << ' ';
		}
		else
		{
		    std::cout << '(';
		    for (int i=0;i<dpsize;++i)
		    {
			std::cout << b[(x+y*dims[bid][0]+z*dims[bid][1]*dims[bid][0])*dpsize+i] << ", ";
		    }
		    std::cout << ") ";
		}	
	    }
	    std::cout << std::endl;
	}
    }
}    


// Copy a 3d region from a flat array into a buffer
void Block::copyToBuffer(unsigned char bid, double* src)
{
    if (bid==13)	// there is no buffer for block 13
    {
	return;
    }
    // where does the strided content start?
    double* start=src+startOffset(bid%3, bid%9/3, bid/9);
    double* dest=outbuffptr[bid];

    size_t zlim=dims[bid][2];	// how big is the block
    size_t ylim=dims[bid][1];
    size_t xlim=dims[bid][0];
    size_t totaly=(2*inset+ymidlen);
    for (size_t z=0;z<zlim;++z)
    {
	for (size_t y=0;y<ylim;++y)
	{ 
	    memcpy(dest, start, xlim*sizeof(double)*dpsize);
	    dest+=xlim*dpsize;
	    start+=(2*inset+xmidlen)*dpsize;		
	}
	// we are at the end of the slab so we need to jump to the next level up
	start+=((totaly-ylim)*(2*inset+xmidlen))*dpsize;
    }
}


// Copy a 3d region from a buffer into a flat array
void Block::copyFromBuffer(unsigned char bid, double* dest)
{
    if (bid==13)	// there is no buffer for block 13
    {
	return;
    }      
    double* start=dest+startOffset(bid%3, bid%9/3, bid/9);
    double* src=inbuffptr[bid];
    size_t zlim=dims[bid][2];	// how big is the block
    size_t ylim=dims[bid][1];
    size_t xlim=dims[bid][0];
    size_t totaly=(2*inset+ymidlen);
    for (size_t z=0;z<zlim;++z)
    {
	for (size_t y=0;y<ylim;++y)
	{
	    memcpy(start, src, xlim*sizeof(double)*dpsize);
	    src+=xlim*dpsize;
	    start+=(2*inset+xmidlen)*dpsize;
	}
	// we are at the end of the slab so we need to jump to the next level up
	start+=(totaly-ylim)*(2*inset+xmidlen)*dpsize;
    }	
}

// Returns the MPI message tag to use for a transfer between the two subblocks
int getTag(unsigned char sourcex, unsigned char sourcey, unsigned char sourcez, unsigned char targetx, unsigned char targety, unsigned char targetz)
{
    return sourcex*100000+sourcey*10000+sourcez*1000+targetx*100+targety*10+targetz;
}

// computes the tag based on the destination and the direction it comes from
// the booleans indicate whether a negative shift in that direction is required
int getTag(unsigned char destx, unsigned char desty, unsigned char destz, bool deltax, bool deltay, bool deltaz)
{
    unsigned char sourcex=deltax?2:destx;
    unsigned char sourcey=deltay?2:desty;
    unsigned char sourcez=deltaz?2:destz;
  
    return sourcex*100000+sourcey*10000+sourcez*1000+destx*100+desty*10+destz;  
}


// the booleans indicate whether a negative shift in that direction is required
unsigned char getSrcBuffID(unsigned char destx, unsigned char desty, unsigned char destz, bool deltax, bool deltay, bool deltaz)
{
    unsigned char sourcex=deltax?2:destx;
    unsigned char sourcey=deltay?2:desty;
    unsigned char sourcez=deltaz?2:destz;
  
    return sourcex+sourcey*3+sourcez*9;
}

