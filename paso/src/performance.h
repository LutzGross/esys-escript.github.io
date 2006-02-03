/* $Id$ */

/**************************************************************/

/* Paso: perfomance monitor interface using papi*/

/**************************************************************/

/* Copyrights by ACcESS Australia 2003/04 */
/* Author: gross@access.edu.au */

/**************************************************************/

#ifndef INC_PASO_PERFORMANCE
#define INC_PASO_PERFORMANCE

#ifdef PAPI

#include <papi.h>

Paso_Performance* Performance_init();
void Performance_start_solver(Paso_Performance&);
void Performance_start_preconditioner_init(Paso_Performance&);
void Performance_start_preconditioner(Paso_Performance&);
void Performance_start_mvm(Paso_Performance&);
void Performance_end_preconditioner(Paso_Performance&);
void Performance_end_preconditioner_init(Paso_Performance&);
void Performance_end_mvm(Paso_Performance&);
void Performance_end_solver(Paso_Performance&);
void Performance_finalized(Paso_Performance&);

struct Paso_Performance {


};
typedef Paso_Performance struct Paso_Performance;

#endif


void Performance_init(void*) {
   #ifdef PAPI
      /* Initialize the PAPI library */
      retval = PAPI_library_init(PAPI_VER_CURRENT);
      if (retval != PAPI_VER_CURRENT && retval > 0) {
        fprintf(stderr,"PAPI library version mismatch!\n");
        exit(1);
      }
      if (retval < 0) {
        fprintf(stderr, “Initialization error!\n”);
        exit(1);
      }
   #endif
   NUM_EVENTS=5;
   events={PAPI_FP_OPS,PAPI_L1_DCM,PAPI_L2_DCM,PAPI_L3_DCM};

   out->all
   out->solver
   out->preconditioner_init
   out->preconditioner
   out->mvm
}

Performance_start(Paso_Performance& pp) {
    pp->



}
