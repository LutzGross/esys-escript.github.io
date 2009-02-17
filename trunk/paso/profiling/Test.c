#include <stdio.h>
#include "paso/Common.h"
#include "paso/Solver.h"
#include "paso/SystemMatrix.h"
#include "Paso_tests.h"
#include <math.h>

#define PI (3.141592653589793)

double Lsup(double* x, int n) {
    double max=0;
    int i;
    
    for (i=0;i<n;i++) {
        max=MAX(ABS(x[i]),max);
    }
    
    return max;
}

int main (int argc, char *argv[]) {
    Paso_SystemMatrix *A = NULL;
    double *b,*x,*x_ref;
    dim_t i,n,level=1;
    double *error;
    double Lsuperror;
   
    Paso_Options options;
    Paso_Options_setDefaults(&options);

   
    if (argc<2) {
        fprintf(stderr,"Please enter the filename\n");
        return -1;    
    }

    if (argc>2) {
      level=atoi(argv[2]);
    }
    
    A=MEMALLOC(1,Paso_SystemMatrix);
   
    A=Paso_SystemMatrix_loadMM_toCSR(argv[1]);

    if (A==NULL) {
      fprintf(stderr,"CSR Matrix not Loaded\n");
      return 0;
    }
    n=Paso_SystemMatrix_getTotalNumRows(A);
    b=MEMALLOC(n,double);
    x=MEMALLOC(n,double);
    x_ref=MEMALLOC(n,double);
    error=MEMALLOC(n,double);
    
    if(argc==4) {
        Paso_RHS_loadMM_toCSR(argv[3],b,n);
        Paso_test_run(A,b,level);
    }
    else {
        for(i=0;i<n;i++) {
          x_ref[i]=cos(i*PI/n);
        }
        /*  raw scaled vector update operation: out = alpha * A * in + beta * out */
        Paso_SystemMatrix_MatrixVector_CSR_OFFSET0(DBLE(1),A,x_ref,DBLE(0),b);

        if (level==3) {
            options.method=PASO_PCG;
            options.verbose=FALSE;
            options.preconditioner=PASO_JACOBI;
            

    
            Paso_solve(A,x,b,&options);

            
            for(i=0;i<n;i++) {
              error[i]=x[i]-x_ref[i];
            }
            Lsuperror=Lsup(error,n)/Lsup(x_ref,n);
            fprintf(stdout,"Lsup error Jacobi %e\n",Lsuperror);
            
            A->solver=NULL;
            options.method=PASO_PCG;
            options.preconditioner=PASO_AMG;
     
            Paso_solve(A,x,b,&options);
            
            for(i=0;i<n;i++) {
              error[i]=x[i]-x_ref[i];
            }
            Lsuperror=Lsup(error,n)/Lsup(x_ref,n);
            fprintf(stdout,"Lsup error AMG %e\n",Lsuperror);
        }
        else {
            Paso_test_run(A,b,level);            
        }
 
    }


    
    MEMFREE(b);
    MEMFREE(x);
    MEMFREE(x_ref);
    MEMFREE(error);
    Paso_SystemMatrix_free(A);
return 1;
}


