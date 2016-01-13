/* $Id$ */

/**************************************************************/

/*   Paso: system matrix pattern                            */

/**************************************************************/

/*   Copyrights by ACcESS Australia 2004,2005 */
/*   Author: gross@access.edu.au */

/**************************************************************/

#ifndef INC_PASO_SYSTEMMATRIXPATTERN
#define INC_PASO_SYSTEMMATRIXPATTERN

#ifdef MSVC
#ifdef PASO_EXPORTS
#define PASO_DLL __declspec(dllexport)
#else
#define PASO_DLL __declspec(dllimport)
#endif
#else
#define PASO_DLL
#endif


#include "Common.h"

#ifdef __cplusplus
extern "C" {
#endif

/**************************************************************/

/* Typically INDEX_OFFSET and PTR_OFFSET are 1 if the matrix is
   handed over to a library written in FORTRAN. */

#define INDEX_OFFSET 0
#define PTR_OFFSET 0

typedef struct Paso_SystemMatrixPattern {
  dim_t n_ptr;
  dim_t n_index;
  index_t* ptr;
  index_t* index;
  dim_t len;
  dim_t reference_counter;

} Paso_SystemMatrixPattern;


/*  interfaces: */

PASO_DLL Paso_SystemMatrixPattern* Paso_SystemMatrixPattern_alloc(dim_t, index_t*,index_t*);
PASO_DLL Paso_SystemMatrixPattern* Paso_SystemMatrixPattern_reference(Paso_SystemMatrixPattern*);
PASO_DLL void Paso_SystemMatrixPattern_dealloc(Paso_SystemMatrixPattern*);
PASO_DLL int Paso_comparIndex(const void *,const void *);
PASO_DLL Paso_SystemMatrixPattern* Paso_SystemMatrixPattern_getSubpattern(Paso_SystemMatrixPattern*,dim_t,index_t*,index_t*);
PASO_DLL void Paso_SystemMatrixPattern_mis(Paso_SystemMatrixPattern* pattern_p, index_t* mis_marker);
PASO_DLL Paso_SystemMatrixPattern* Paso_SystemMatrixPattern_unrollBlocks(Paso_SystemMatrixPattern*, dim_t,dim_t);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* #ifndef INC_PASO_SYSTEMPATTERN */

/*
 * $Log$
 * Revision 1.2  2005/09/15 03:44:38  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-15
 *
 * Revision 1.1.2.1  2005/09/05 06:29:47  gross
 * These files have been extracted from finley to define a stand alone libray for iterative
 * linear solvers on the ALTIX. main entry through Paso_solve. this version compiles but
 * has not been tested yet.
 *
 *
 */
