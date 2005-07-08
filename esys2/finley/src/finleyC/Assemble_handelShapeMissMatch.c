/* $Id$ */

/**************************************************************/

/*    assemblage routines: */

/* these routines handel cases where there ia a miss match in the number of nodes and the number of shape funtion: */

/* _in: V[N,numNodes0,numNodes1] is reduced into V[N,numShapes0,numShapes1] */
/* _out: V[N,numShapes0,numShapes1] is expanded into V[N,numNodes0,numNodes1] */
/* in both cases the leading dimenisons are unchanged */


/**************************************************************/

/*   Copyrights by ACcESS Australia, 2003,2004 */
/*   author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Common.h"
#include "Assemble.h"


/**************************************************************/

/* this takes the mean value: */

void Finley_Assemble_handelShapeMissMatch_Mean_in(dim_t N, dim_t numNodes0,dim_t numNodes1, double* V,dim_t numShapes0, dim_t numShapes1) {
  dim_t l0=MIN(numShapes0,numNodes0-numShapes0);
  dim_t l1=MIN(numShapes1,numNodes1-numShapes1);
  dim_t i,k0,k1;
  if (MAX(l0,l1)>0) {
    if (numNodes1==numShapes1) {
      for (k0=0;k0<l0;k0++) {
	for (i=0;i<N;i++) V[INDEX2(i,k0,N)]= ( V[INDEX2(i,k0,N)]+V[INDEX2(i,k0+numShapes0,N)])/2.;
      }
    } else {
      for (k0=0;k0<l0;k0++) {
	for (k1=0;k1<l1;k1++) {
	  for (i=0;i<N;i++) {
	    V[INDEX3(i,k0,k1,N,numNodes0)]=
	      ( V[INDEX3(i,k0           ,k1           ,N,numNodes0)]+
		V[INDEX3(i,k0+numShapes0,k1           ,N,numNodes0)]+
		V[INDEX3(i,k0           ,k1+numShapes1,N,numNodes0)]+
		V[INDEX3(i,k0+numShapes0,k1+numShapes1,N,numNodes0)])/4.;
	  }
	}
      }
    }
  }
}

void Finley_Assemble_handelShapeMissMatch_Mean_out(dim_t N, dim_t numNodes0,dim_t numNodes1, double* V,dim_t numShapes0, dim_t numShapes1) {
  double RTMP;
  dim_t i,k0,k1;
  dim_t l0=MIN(numShapes0,numNodes0-numShapes0);
  dim_t l1=MIN(numShapes1,numNodes1-numShapes1);
  if (MAX(l0,l1)>0) {
    if (numNodes1==numShapes1) {
      for (k0=0;k0<l0;k0++) {
	for (i=0;i<N;i++) {
	  RTMP=V[INDEX2(i,k0,N)];
	  V[INDEX2(i,k0           ,N)]=0.;
	  V[INDEX2(i,k0+numShapes0,N)]=0.;
   
	  V[INDEX2(i,k0           ,N)]+=RTMP/2.;
	  V[INDEX2(i,k0+numShapes0,N)]+=RTMP/2.;
	}
      }
    } else {
      for (k0=0;k0<l0;k0++) {
	for (k1=0;k1<l1;k1++) {
	  for (i=0;i<N;i++) {
	    RTMP=V[INDEX3(i,k0,k1,N,numNodes0)];
	    V[INDEX3(i,k0           ,k1           ,N,numNodes0)]=0.;
	    V[INDEX3(i,k0+numShapes0,k1           ,N,numNodes0)]=0.;
	    V[INDEX3(i,k0           ,k1+numShapes1,N,numNodes0)]=0.;
	    V[INDEX3(i,k0+numShapes0,k1+numShapes1,N,numNodes0)]=0.;
   
	    V[INDEX3(i,k0           ,k1           ,N,numNodes0)]+=RTMP/4.;
	    V[INDEX3(i,k0+numShapes0,k1           ,N,numNodes0)]+=RTMP/4.;
	    V[INDEX3(i,k0           ,k1+numShapes1,N,numNodes0)]+=RTMP/4.;
	    V[INDEX3(i,k0+numShapes0,k1+numShapes1,N,numNodes0)]+=RTMP/4.;
	  }
	}
      }
    }
  }
}

void Finley_Assemble_handelShapeMissMatch_Step_in(dim_t N, dim_t numNodes0,dim_t numNodes1, double* V,dim_t numShapes0, dim_t numShapes1) {
  dim_t l0=MIN(numShapes0,numNodes0-numShapes0);
  dim_t l1=MIN(numShapes1,numNodes1-numShapes1);
  dim_t i,k0,k1;
  if (MAX(l0,l1)>0) {
    if (numNodes1==numShapes1) {
      for (k0=0;k0<l0;k0++) {
	for (i=0;i<N;i++) {
	  V[INDEX2(i,k0,N)]= -V[INDEX2(i,k0,N)]+V[INDEX2(i,k0+numShapes0,N)];
	}
      }
    } else {
      for (k0=0;k0<l0;k0++) {
	for (k1=0;k1<l1;k1++) {
	  for (i=0;i<N;i++) {
	    V[INDEX3(i,k0,k1,N,numNodes0)]=
	      V[INDEX3(i,k0           ,k1           ,N,numNodes0)]+
	      V[INDEX3(i,k0+numShapes0,k1           ,N,numNodes0)]-
	      V[INDEX3(i,k0           ,k1+numShapes1,N,numNodes0)]-
	      V[INDEX3(i,k0+numShapes0,k1+numShapes1,N,numNodes0)];
	  }
	}
      }
    }
  }
}

void Finley_Assemble_handelShapeMissMatch_Step_out(dim_t N, dim_t numNodes0,dim_t numNodes1, double* V,dim_t numShapes0, dim_t numShapes1) {
  double RTMP;
  dim_t i,k0,k1;
  dim_t l0=MIN(numShapes0,numNodes0-numShapes0);
  dim_t l1=MIN(numShapes1,numNodes1-numShapes1);
  if (MAX(l0,l1)>0) {
    if (numNodes1==numShapes1) {
      for (k0=0;k0<l0;k0++) {
	for (i=0;i<N;i++) {
	  RTMP=V[INDEX2(i,k0,N)];
	  V[INDEX2(i,k0           ,N)]=0.;
	  V[INDEX2(i,k0+numShapes0,N)]=0.;
   
	  V[INDEX2(i,k0           ,N)]-=RTMP;
	  V[INDEX2(i,k0+numShapes0,N)]+=RTMP;
	}
      }
    } else {
      for (k0=0;k0<l0;k0++) {
	for (k1=0;k1<l1;k1++) {
	  for (i=0;i<N;i++) {
	    RTMP=V[INDEX3(i,k0,k1,N,numNodes0)];
	    V[INDEX3(i,k0           ,k1           ,N,numNodes0)]=0.;
	    V[INDEX3(i,k0+numShapes0,k1           ,N,numNodes0)]=0.;
	    V[INDEX3(i,k0           ,k1+numShapes1,N,numNodes0)]=0.;
	    V[INDEX3(i,k0+numShapes0,k1+numShapes1,N,numNodes0)]=0.;
   
	    V[INDEX3(i,k0           ,k1           ,N,numNodes0)]+=RTMP;
	    V[INDEX3(i,k0+numShapes0,k1           ,N,numNodes0)]-=RTMP;
	    V[INDEX3(i,k0           ,k1+numShapes1,N,numNodes0)]-=RTMP;
	    V[INDEX3(i,k0+numShapes0,k1+numShapes1,N,numNodes0)]+=RTMP;
	  }
	}
      }
    }
  }
}

/*
 * $Log$
 * Revision 1.2  2005/07/08 04:07:47  jgs
 * Merge of development branch back to main trunk on 2005-07-08
 *
 * Revision 1.1.1.1.2.1  2005/06/29 02:34:48  gross
 * some changes towards 64 integers in finley
 *
 * Revision 1.1.1.1  2004/10/26 06:53:57  jgs
 * initial import of project esys2
 *
 * Revision 1.1  2004/07/02 04:21:13  gross
 * Finley C code has been included
 *
 *
 */
