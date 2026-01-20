
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
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


BlockGrid2::BlockGrid2(coord_t x, coord_t y)
 : xmax(x), ymax(y)
{}

neighbourID_t BlockGrid2::getNID(coord_t x, coord_t y) const
{
    return x+y*(xmax+1);    
}

// generate all incoming com messages for this block.
// for each subblock (27 of them), there may be an x, y, z direction to search in  
void BlockGrid2::generateInNeighbours(coord_t blockx, coord_t blocky, messvec& v)
{
    neighbourID_t myid=getNID(blockx, blocky);
    unsigned char deltax=0;
    unsigned char deltay=0;
    deltay=blocky?1:0;		// do not look at a lower row of blocks if there isn't one
    for (unsigned char y=0;y<3;++y)
    {
      deltax=blockx?1:0;
      for (unsigned char x=0;x<3;++x)
      {
	// now we will have a set of possible directions to look in (delta?==0 means don't attempt to
	// change this component	  
	if (deltax+deltay)	// if we have an import to do
	{
	    coord_t srcx=blockx-deltax;
	    coord_t srcy=blocky-deltay;
	    message m;
	    m.sourceID=getNID(srcx, srcy);
	    m.destID=myid;
	    m.tag=getTag2(x,y,deltax==1, deltay==1);
	    m.srcbuffid=getSrcBuffID2(x,y,deltax==1, deltay==1);
	    m.destbuffid=x+y*3;
	    v.push_back(m);
	}
	deltax=0;	// we are no longer on the left hand edge
      }
      deltay=0;	// only y=0 imports
    } 
}


// generate all outgoing com messages for this block
void BlockGrid2::generateOutNeighbours(coord_t blockx, coord_t blocky, messvec& v)
{
    messvec vv;
    neighbourID_t myid=getNID(blockx, blocky);
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
	    if (x+y>0)	// don't look at your own neighbours
	    {
		generateInNeighbours(blockx+x, blocky+y, vv);
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

  
double* Block2::getOutBuffer(unsigned char subx, unsigned char suby)
{
    unsigned char bid=subx+suby*3;	// (bid is "blockid")
    if (bid==4)	// there is no buffer for block 4
    {
	return 0;	// don't ask for this buffer because refusal may offend
    }
    return outbuffptr[bid];	
}

double* Block2::getInBuffer(unsigned char subx, unsigned char suby)
{
    unsigned char bid=subx+suby*3;	// (bid is "blockid")
    if (bid==4)	// there is no buffer for block 4
    {
	return 0;	// don't ask for this buffer because refusal may offend
    }
    return inbuffptr[bid];
}

size_t Block2::getBuffSize(unsigned char subx, unsigned char suby)
{
    unsigned char bid=subx+suby*3;	// (bid is "blockid")
    if (bid==4)	// there is no buffer for block 4
    {
	return 0;	
    }
    return dims[bid][0]*dims[bid][1]*dpsize;	
}

double* Block2::getOutBuffer(unsigned char bid)
{
    if (bid==4)	// there is no buffer for block 4
    {
	return 0;	// don't ask for this buffer because refusal may offend
    }
    return outbuffptr[bid];	
}

double* Block2::getInBuffer(unsigned char bid)
{
    if (bid==4)	// there is no buffer for block 4
    {
	return 0;	// don't ask for this buffer because refusal may offend
    }
    return inbuffptr[bid];
}

size_t Block2::getBuffSize(unsigned char bid)
{
    if (bid==4)	// there is no buffer for block 4
    {
	return 0;	
    }
    return dims[bid][0]*dims[bid][1]*dpsize;	
}



Block2::~Block2()
{
    delete[] inbuff;
    delete[] outbuff;
}

void Block2::populateDimsTable()
{
    // start off by setting everything to inset*inset
    for (int i=0;i<9;++i)
    {
	for (int j=0;j<2;++j)
	{
	    dims[i][j]=inset;
	}
    }
    // fix middle column
    dims[1][0]=dims[7][0]=dims[4][0]=xmidlen;
    // fix E and W middle
    dims[3][1]=dims[4][1]=dims[5][1]=ymidlen;
}


// gives the relative start location within a flat 
void Block2::populateOffsetTable(size_t inset, size_t xmidlen, size_t ymidlen)
{
  
    size_t cur=0;
    for (int i=0;i<9;++i)
    {
	flatoffsets[i]=cur;
	cur+=dims[i][0]*dims[i][1]*dpsize;
    }
    for (int i=0;i<4;++i)
    {
	buffoffsets[i]=flatoffsets[i];
    }
    buffoffsets[4]=0;
    for (int i=5;i<9;++i)
    {
	buffoffsets[i]=flatoffsets[i]-flatoffsets[5]+flatoffsets[4];
    }

}


void Block2::createBuffArrays(double* startaddress, double* buffptr[27], size_t inset, size_t xmidlen, size_t ymidlen)
{
    buffptr[0]=startaddress;
    for (int i=0;i<9;++i)
    {
	buffptr[i]=startaddress+buffoffsets[i];
    }
    buffptr[4]=0;	// since the middle should never be sent anywhere
}

void Block2::setUsed(unsigned char buffid)
{
    used[buffid]=true;
}

void Block2::copyAllToBuffer(double* src)
{
    for (unsigned char i=0;i<4;++i)
    {
	copyToBuffer(i, src);
      
    }
    for (unsigned char i=5;i<9;++i)
    {
	copyToBuffer(i, src);
    }
}
    
void Block2::copyUsedFromBuffer(double* dest)
{
    for (unsigned char i=0;i<9;++i)
    {
	if (used[i])
	{
	    copyFromBuffer(i, dest);
	}
    }
}

// s? specifiy the size (in points) of each dimension
// maxb? gives the largest block number in each dimension in the overall grid (number from zero)
Block2::Block2(size_t sx, size_t sy, size_t inset, size_t xmidlen, 
	     size_t ymidlen, unsigned int dpp) : dpsize(dpp)
{
    this->sx=sx;
    this->sy=sy;
    this->inset=inset;
    this->xmidlen=xmidlen;
    this->ymidlen=ymidlen;
    populateDimsTable();	
    
    size_t totalbuff=0;
    for (int i=0;i<9;++i)
    {
	used[i]=false;
	totalbuff+=dims[i][0]*dims[i][1];
    }
    totalbuff-=dims[4][0]*dims[4][1];
	
    totalbuff*=dpsize;		// to account for non-scalars
    inbuff=new double[totalbuff];
    outbuff=new double[totalbuff];
    memset(inbuff,0,totalbuff*sizeof(double));
    memset(outbuff,0,totalbuff*sizeof(double));
    populateOffsetTable(inset, xmidlen, ymidlen);
    createBuffArrays(inbuff, inbuffptr, inset, xmidlen, ymidlen);
    createBuffArrays(outbuff, outbuffptr, inset, xmidlen, ymidlen);

}  

// where does the subblock specified start in a source array
size_t Block2::startOffset(unsigned char subx, unsigned char suby)
{
    size_t off=0;
    off+=((subx==0) ? 0 :((subx==1)?inset:(inset+xmidlen) ));
    // now the y component
    size_t ystep=((suby==0) ? 0 :((suby==1)?inset:(inset+ymidlen) ));
    off+=ystep*(2*inset+xmidlen);
    return off*dpsize;
}


// This is only intended for debugging, so I have not made it fast
void Block2::displayBlock(unsigned char subx, unsigned char suby, bool out)
{
    
    unsigned char bid=subx+suby*3;	// (bid is "blockid")	
    double* b=out?outbuffptr[bid]:inbuffptr[bid];
    for (int y=0;y<dims[bid][1];++y)
    {
	for (int x=0;x<dims[bid][0];++x)
	{
	    if (dpsize==1)
	    {
		std::cout << b[x+y*dims[bid][0]] << ' ';
	    }
	    else
	    {
		std::cout << '(';
		for (int i=0;i<dpsize;++i)
		{
		    std::cout << b[(x+y*dims[bid][0])*dpsize+i] << ", ";
		}
		std::cout << ") ";
	    }	
	}
	std::cout << std::endl;
    }
}    


// Copy a 3d region from a flat array into a buffer
void Block2::copyToBuffer(unsigned char bid, double* src)
{
    if (bid==4)	// there is no buffer for block 4
    {
	return;
    }
    // where does the strided content start?
    double* start=src+startOffset(bid%3, bid/3);
    double* dest=outbuffptr[bid];
    size_t ylim=dims[bid][1];	// how big is the block
    size_t xlim=dims[bid][0];
    
// // testing for mismatch in the copy
// for (int i=0;i<ylim*xlim*dpsize;++i)
// {
// dest[i]=-42;    
// }
    
    for (size_t y=0;y<ylim;++y)
    {     
	memcpy(dest, start, xlim*sizeof(double)*dpsize);
	dest+=xlim*dpsize;
	start+=(2*inset+xmidlen)*dpsize;		
    }
}


// Copy a 3d region from a buffer into a flat array
void Block2::copyFromBuffer(unsigned char bid, double* dest)
{
    if (bid==4)	// there is no buffer for block 4
    {
	return;
    }      
    double* start=dest+startOffset(bid%3, bid/3);
    double* src=inbuffptr[bid];
    size_t ylim=dims[bid][1];	// how big is the block
    size_t xlim=dims[bid][0];
    for (size_t y=0;y<ylim;++y)
    {
	memcpy(start, src, xlim*sizeof(double)*dpsize);
	src+=xlim*dpsize;
	start+=(2*inset+xmidlen)*dpsize;
    }
}

// Returns the MPI message tag to use for a transfer between the two subblocks
int getTag2(unsigned char sourcex, unsigned char sourcey, unsigned char targetx, unsigned char targety)
{
    return sourcex*10000+sourcey*1000+targetx*100+targety*10;
}

// computes the tag based on the destination and the direction it comes from
// the booleans indicate whether a negative shift in that direction is required
int getTag2(unsigned char destx, unsigned char desty, bool deltax, bool deltay)
{
    unsigned char sourcex=deltax?2:destx;
    unsigned char sourcey=deltay?2:desty;
    return sourcex*10000+sourcey*1000+destx*100+desty*10;  
}


// the booleans indicate whether a negative shift in that direction is required
unsigned char getSrcBuffID2(unsigned char destx, unsigned char desty, bool deltax, bool deltay)
{
    unsigned char sourcex=deltax?2:destx;
    unsigned char sourcey=deltay?2:desty;
  
    return sourcex+sourcey*3;
}

