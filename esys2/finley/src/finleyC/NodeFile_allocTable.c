/* $Id$ */
/**************************************************************/

/*   Finley: Mesh: NodeFile */

/**************************************************************/

/*   Copyrights by ACcESS Australia 2003 */
/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Finley.h"
#include "NodeFile.h"

/**************************************************************/

/*  allocates the node table within an node file to hold numNodes of nodes. The LinearTo mapping, if it exists, */
/*  is deallocated. use Finley_Mesh_setLinearMesh to create a new one. */

void Finley_NodeFile_allocTable(Finley_NodeFile* in ,int numNodes) {
  maybelong *Id2=NULL, *Tag2=NULL, *degreeOfFreedom2=NULL, *reducedDegreeOfFreedom2=NULL, *toReduced2=NULL;
  double *Coordinates2=NULL;
  maybelong n,i;
  
  /*  allocate memory: */
  
  Id2=(maybelong*) MEMALLOC(numNodes*sizeof(maybelong));
  Coordinates2=(double*) MEMALLOC(numNodes*in->numDim*sizeof(double));
  Tag2=(maybelong*) MEMALLOC(numNodes*sizeof(maybelong));
  degreeOfFreedom2=(maybelong*) MEMALLOC(numNodes*sizeof(maybelong));
  reducedDegreeOfFreedom2=(maybelong*) MEMALLOC(numNodes*sizeof(maybelong));
  toReduced2=(maybelong*) MEMALLOC(numNodes*sizeof(maybelong));
  
  /*  if fine, deallocate the old table and replace by new: */
  
  if (Finley_checkPtr(Id2) || Finley_checkPtr(Coordinates2) || Finley_checkPtr(Tag2) ||
      Finley_checkPtr(degreeOfFreedom2) || Finley_checkPtr(reducedDegreeOfFreedom2) || Finley_checkPtr(toReduced2) ) {
    MEMFREE(Id2);
    MEMFREE(Coordinates2);
    MEMFREE(Tag2);
    MEMFREE(degreeOfFreedom2);
    MEMFREE(reducedDegreeOfFreedom2);
    MEMFREE(toReduced2);
  } else { 
    Finley_NodeFile_deallocTable(in);
    in->Id=Id2;
    in->Coordinates=Coordinates2;
    in->Tag=Tag2;
    in->degreeOfFreedom=degreeOfFreedom2;
    in->reducedDegreeOfFreedom=reducedDegreeOfFreedom2;
    in->toReduced=toReduced2;
    in->numNodes=numNodes;
    in->numDegreesOfFreedom=numNodes;
    in->reducedNumDegreesOfFreedom=numNodes;
    in->reducedNumNodes=numNodes;
    /* this initialization makes sure that data are located on the right processor */
    #pragma omp parallel for private(n,i) schedule(static)
    for (n=0;n<numNodes;n++) {
       for (i=0;i<in->numDim;i++) in->Coordinates[INDEX2(i,n,in->numDim)]=0.;
       in->Id[n]=-1;
       in->Tag[n]=-1;
       in->degreeOfFreedom[n]=n;
       in->reducedDegreeOfFreedom[n]=n;
       in->toReduced[n]=n;
    }
  }
  return;
}

/*  deallocates the node table within an node file: */

void Finley_NodeFile_deallocTable(Finley_NodeFile* in) {
  if (in!=NULL) {
    MEMFREE(in->Id);
    MEMFREE(in->Coordinates);
    MEMFREE(in->Tag);
    MEMFREE(in->degreeOfFreedom);
    MEMFREE(in->reducedDegreeOfFreedom);
    MEMFREE(in->toReduced);
    in->numNodes=0;
    in->numDegreesOfFreedom=0;
    in->reducedNumDegreesOfFreedom=0;
    in->reducedNumNodes=0;
  }
}
/* 
* $Log$
* Revision 1.1  2004/10/26 06:53:57  jgs
* Initial revision
*
* Revision 1.1.1.1  2004/06/24 04:00:40  johng
* Initial version of eys using boost-python.
*
*
*/
