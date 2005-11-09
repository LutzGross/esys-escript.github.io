/* $Id$ */

/**************************************************************/

/*   Paso: SystemMatrix and SystemVector */

/**************************************************************/

/*   Copyrights by ACcESS Australia 2003,2004,2005 */
/*   Author: gross@access.edu.au */

/**************************************************************/

#ifndef INC_PASO_OPTIONS
#define INC_PASO_OPTIONS

/* solver options */

#define PASO_DEFAULT 0
#define PASO_DIRECT 1
#define PASO_CHOLEVSKY 2
#define PASO_PCG 3
#define PASO_CR 4
#define PASO_CGS 5
#define PASO_BICGSTAB 6
#define PASO_SSOR 7
#define PASO_ILU0 8
#define PASO_ILUT 9
#define PASO_JACOBI 10
#define PASO_GMRES 11
#define PASO_PRES20 12
#define PASO_LUMPING 13
#define PASO_SCSL 14
#define PASO_MKL 15
#define PASO_UMFPACK 16
#define PASO_NO_REORDERING 17
#define PASO_MINIMUM_FILL_IN 18
#define PASO_NESTED_DISSECTION 19
#define PASO_ITERATIVE 20
#define PASO_PASO 21


typedef struct {
    index_t method;
    index_t package;
    bool_t symmetric;
    double tolerance;

    bool_t verbose;
    bool_t reordering;
    double final_residual;
    index_t preconditioner;
    dim_t iter_max;
    dim_t iter;
    double drop_tolerance;
    double drop_storage;
    dim_t truncation;
    dim_t restart;


} Paso_Options;

/*  interfaces: */

void Paso_Options_setDefaults(Paso_Options* in);
index_t Paso_Options_getPackage(index_t solver,index_t package, bool_t symmetry);
index_t Paso_Options_getSolver(index_t solver,index_t package, bool_t symmetry);


#endif /* #ifndef INC_PASO_OPTIONS */



/*
 * $Log$
 * Revision 1.2  2005/09/15 03:44:38  jgs
 * Merge of development branch dev-02 back to main trunk on 2005-09-15
 *
 * Revision 1.1.2.2  2005/09/07 00:59:08  gross
 * some inconsistent renaming fixed to make the linking work.
 *
 * Revision 1.1.2.1  2005/09/05 06:29:46  gross
 * These files have been extracted from finley to define a stand alone libray for iterative
 * linear solvers on the ALTIX. main entry through Paso_solve. this version compiles but
 * has not been tested yet.
 *
 *
 */
