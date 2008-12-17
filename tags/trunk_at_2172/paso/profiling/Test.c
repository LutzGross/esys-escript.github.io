#include <stdio.h>
#include "paso/Common.h"
#include "paso/Solver.h"
#include "paso/SystemMatrix.h"
#include "Paso_tests.h"

int main (int argc, char *argv[]) {
    Paso_SystemMatrix *A = NULL;
    double *b,*x;
    dim_t i,n,level=1;
   
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
    
    if(argc==4) {
        Paso_RHS_loadMM_toCSR(argv[3],b,n);
    }
    else {
        for(i=0;i<n;i++) {
          b[i]=1;
        }    
    }
    
    
    Paso_test_run(A,b,level);
    MEMFREE(b);
    MEMFREE(x);
    Paso_SystemMatrix_free(A);
return 1;
}                                                                                                                                                                                                     
