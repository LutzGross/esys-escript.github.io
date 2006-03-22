/*
 ************************************************************
 *          Copyright 2006 by ACcESS MNRF                   *
 *                                                          *
 *              http://www.access.edu.au                    *
 *       Primary Business: Queensland, Australia            *
 *  Licensed under the Open Software License version 3.0    *
 *     http://www.opensource.org/licenses/osl-3.0.php       *
 *                                                          *
 ************************************************************
*/
/**************************************************************/

/*    assemblage routines: copies data between different types nodal representation   */

/**************************************************************/

/*   author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "Util.h"
#include "Assemble.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/******************************************************************************************************/


void Finley_Assemble_CopyNodalData(Finley_NodeFile* nodes,escriptDataC* out,escriptDataC* in) {
    if (nodes==NULL) return;
    dim_t n,i;
    dim_t numComps=getDataPointSize(out);
    type_t in_data_type=getFunctionSpaceType(in);
    type_t out_data_type=getFunctionSpaceType(out);
    Finley_resetError();

    /* check out and in */
    if (numComps!=getDataPointSize(in)) {
       Finley_setError(TYPE_ERROR,"__FILE__: number of components of input and output Data do not match.");
    } else if (!isExpanded(out)) {
       Finley_setError(TYPE_ERROR,"__FILE__: expanded Data object is expected for output data.");
    }

    if (in_data_type == FINLEY_NODES) {
        if (! numSamplesEqual(in,1,nodes->numNodes)) {
               Finley_setError(TYPE_ERROR,"__FILE__: illegal number of samples of input Data object");
       }
    } else if (in_data_type == FINLEY_DEGREES_OF_FREEDOM) {
        if (! numSamplesEqual(in,1,nodes->numDegreesOfFreedom)) {
               Finley_setError(TYPE_ERROR,"__FILE__: illegal number of samples of input Data object");
       }
    } else if (in_data_type == FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
        if (! numSamplesEqual(in,1,nodes->reducedNumDegreesOfFreedom)) {
               Finley_setError(TYPE_ERROR,"__FILE__: illegal number of samples of input Data object");
       }
    } else {
       Finley_setError(TYPE_ERROR,"__FILE__: illegal function space type for target object");
    }
    
    if (out_data_type == FINLEY_NODES) {
        if (! numSamplesEqual(out,1,nodes->numNodes)) {
               Finley_setError(TYPE_ERROR,"__FILE__: illegal number of samples of output Data object");
       }
    } else if (out_data_type == FINLEY_DEGREES_OF_FREEDOM) {
        if (! numSamplesEqual(out,1,nodes->numDegreesOfFreedom)) {
               Finley_setError(TYPE_ERROR,"__FILE__: illegal number of samples of output Data object");
       }
    } else if (out_data_type == FINLEY_REDUCED_DEGREES_OF_FREEDOM) {
        if (! numSamplesEqual(out,1,nodes->reducedNumDegreesOfFreedom)) {
               Finley_setError(TYPE_ERROR,"__FILE__: illegal number of samples of output Data object");
       }
    } else {
       Finley_setError(TYPE_ERROR,"__FILE__: illegal function space type for source object");
    }

    /* now we can start */

    if (Finley_noError()) {
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
           } else if (out_data_type == FINLEY_DEGREES_OF_FREEDOM && numSamplesEqual(out,1,nodes->reducedNumDegreesOfFreedom)) {
              #pragma omp parallel for private(n) schedule(static)
              for (i=0;i<nodes->numNodes;i++) {
                  n=nodes->reducedDegreeOfFreedom[i];
                  if (n>=0) Finley_copyDouble(numComps,getSampleData(in,n),getSampleData(out,nodes->degreeOfFreedom[i]));
              }
           } else if (out_data_type == FINLEY_NODES && numSamplesEqual(out,1,nodes->reducedNumDegreesOfFreedom)) {
              #pragma omp parallel for private(n) schedule(static)
              for (i=0;i<nodes->numNodes;i++) {
                  n=nodes->reducedDegreeOfFreedom[i];
                  if (n>=0) Finley_copyDouble(numComps,getSampleData(in,n),getSampleData(out,i));
              }
           } else {
             Finley_setError(TYPE_ERROR,"__FILE__: cannot copy from data on reduced degrees of freedom");
           }
        }
   }
   return;
}
/*
 * $Log$
 * Revision 1.4  2005/09/15 03:44:21  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-15
 *
 * Revision 1.3  2005/08/12 01:45:42  jgs
 * erge of development branch dev-02 back to main trunk on 2005-08-12
 *
 * Revision 1.2.2.3  2005/09/07 06:26:17  gross
 * the solver from finley are put into the standalone package paso now
 *
 * Revision 1.2.2.2  2005/08/09 02:23:12  gross
 * print statement removed
 *
 * Revision 1.2.2.1  2005/08/03 09:55:33  gross
 * ContactTest is passing now./mk install!
 *
 * Revision 1.2  2005/07/08 04:07:45  jgs
 * Merge of development branch back to main trunk on 2005-07-08
 *
 * Revision 1.1.1.1.2.1  2005/06/29 02:34:46  gross
 * some changes towards 64 integers in finley
 *
 * Revision 1.1.1.1  2004/10/26 06:53:56  jgs
 * initial import of project esys2
 *
 * Revision 1.2  2004/07/21 05:00:54  gross
 * name changes in DataC
 *
 * Revision 1.1  2004/07/02 04:21:13  gross
 * Finley C code has been included
 *
 *
 */
