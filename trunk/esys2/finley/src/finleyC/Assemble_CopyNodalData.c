/* $Id$ */

/**************************************************************/

/*    assemblage routines: copies data between different types nodal representation   */

/**************************************************************/

/*   Copyrights by ACcESS Australia, 2003,2004 */
/*   author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "escript/Data/DataC.h"
#include "Finley.h"
#include "Util.h"
#include "Assemble.h"
#include "NodeFile.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/******************************************************************************************************/


void Finley_Assemble_CopyNodalData(Finley_NodeFile* nodes,escriptDataC* out,escriptDataC* in) {
    if (nodes==NULL) return;
    int n,i;
    int numComps=getDataPointSize(out);
    int in_data_type=getFunctionSpaceType(in);
    int out_data_type=getFunctionSpaceType(out);

    /* check out and in */
    if (numComps!=getDataPointSize(in)) {
       Finley_ErrorCode=TYPE_ERROR;
       sprintf(Finley_ErrorMsg,"number of components of input and output Data do not match.");
    } else if (!isExpanded(out)) {
       Finley_ErrorCode=TYPE_ERROR;
       sprintf(Finley_ErrorMsg,"expanded Data object is expected for output data.");
    }

    if (in_data_type == FINLEY_NODES) {
        if (! numSamplesEqual(in,1,nodes->numNodes)) {
               Finley_ErrorCode=TYPE_ERROR;
               sprintf(Finley_ErrorMsg,"illegal number of samples of input Data object");
       }
    } else if (in_data_type == FINLEY_DEGREES_OF_FREEDOM) {
        if (! numSamplesEqual(in,1,nodes->numDegreesOfFreedom)) {
               Finley_ErrorCode=TYPE_ERROR;
               sprintf(Finley_ErrorMsg,"illegal number of samples of input Data object");
       }
    } else if (in_data_type == FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
        if (! numSamplesEqual(in,1,nodes->reducedNumDegreesOfFreedom)) {
               Finley_ErrorCode=TYPE_ERROR;
               sprintf(Finley_ErrorMsg,"illegal number of samples of input Data object");
       }
    } else {
       Finley_ErrorCode=TYPE_ERROR;
       sprintf(Finley_ErrorMsg,"illegal function space type for target object");
    }
    
    if (out_data_type == FINLEY_NODES) {
        if (! numSamplesEqual(out,1,nodes->numNodes)) {
               Finley_ErrorCode=TYPE_ERROR;
               sprintf(Finley_ErrorMsg,"illegal number of samples of output Data object");
       }
    } else if (out_data_type == FINLEY_DEGREES_OF_FREEDOM) {
        if (! numSamplesEqual(out,1,nodes->numDegreesOfFreedom)) {
               Finley_ErrorCode=TYPE_ERROR;
               sprintf(Finley_ErrorMsg,"illegal number of samples of output Data object");
       }
    } else if (out_data_type == FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
        if (! numSamplesEqual(out,1,nodes->reducedNumDegreesOfFreedom)) {
               Finley_ErrorCode=TYPE_ERROR;
               sprintf(Finley_ErrorMsg,"illegal number of samples of output Data object");
       }
    } else {
       Finley_ErrorCode=TYPE_ERROR;
       sprintf(Finley_ErrorMsg,"illegal function space type for source object");
    }

    /* now we can start */

    if (Finley_ErrorCode==NO_ERROR) {
        if (in_data_type == FINLEY_NODES) {
           if  (out_data_type == FINLEY_NODES) {
              #pragma omp parallel for private(n) schedule(static)
              for (n=0;n<nodes->numNodes;n++) 
                   Finley_copyDouble(numComps,getSampleData(in,n),getSampleData(out,n));
           } else if (out_data_type == FINLEY_DEGREES_OF_FREEDOM) {
              #pragma omp parallel for private(n) schedule(static)
              for (n=0;n<nodes->numNodes;n++) 
                   Finley_copyDouble(numComps,getSampleData(in,n),getSampleData(out,nodes->degreeOfFreedom[n]));
           } else if (out_data_type == FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
              #pragma omp parallel for private(n,i) schedule(static)
              for (i=0;i<nodes->numNodes;i++) {
                   n=nodes->reducedDegreeOfFreedom[i];
                   if (n>=0) Finley_copyDouble(numComps,getSampleData(in,i),getSampleData(out,n)); 
              }
           }
        } else if (in_data_type == FINLEY_DEGREES_OF_FREEDOM) {
           if  (out_data_type == FINLEY_NODES) {
              #pragma omp parallel for private(n) schedule(static)
              for (n=0;n<nodes->numNodes;n++) 
                   Finley_copyDouble(numComps,getSampleData(in,nodes->degreeOfFreedom[n]),getSampleData(out,n));
           } else if (out_data_type == FINLEY_DEGREES_OF_FREEDOM) {
              #pragma omp parallel for private(n) schedule(static)
              for (n=0;n<nodes->numDegreesOfFreedom;n++) 
                    Finley_copyDouble(numComps,getSampleData(in,n),getSampleData(out,n));
           } else if (out_data_type == FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
              #pragma omp parallel for private(n,i) schedule(static)
              for (i=0;i<nodes->numNodes;i++) {
                  n=nodes->reducedDegreeOfFreedom[i];
                  if (n>=0) Finley_copyDouble(numComps,getSampleData(in,nodes->degreeOfFreedom[i]),getSampleData(out,n));
              }
           }
       } else if (in_data_type == FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
           if (out_data_type == FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
              #pragma omp parallel for private(n) schedule(static)
              for (n=0;n<nodes->reducedNumDegreesOfFreedom;n++) 
                    Finley_copyDouble(numComps,getSampleData(in,n),getSampleData(out,n));
           } else {
             Finley_ErrorCode=TYPE_ERROR;
             sprintf(Finley_ErrorMsg,"cannot copy from data on reduced degrees of freedom");
           }
      }
   }
   return;
}
/*
 * $Log$
 * Revision 1.1  2004/10/26 06:53:56  jgs
 * Initial revision
 *
 * Revision 1.2  2004/07/21 05:00:54  gross
 * name changes in DataC
 *
 * Revision 1.1  2004/07/02 04:21:13  gross
 * Finley C code has been included
 *
 *
 */
