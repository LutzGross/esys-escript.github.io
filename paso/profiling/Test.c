#include <stdio.h>
#include "paso/Common.h"
#include "paso/Solver.h"
#include "paso/SystemMatrix.h"
#include "Paso_tests.h"

int main (int argc, char *argv[]) {
    Paso_SystemMatrix *A = NULL;
    double *b,*x;
    dim_t i,n;
   
    if (argc<2) {
        fprintf(stderr,"Please enter the filename\n");
        return -1;    
    }
    
    A=MEMALLOC(1,Paso_SystemMatrix);
   
    A=Paso_SystemMatrix_loadMM_toCSR(argv[1]);
    if (A==NULL) {
      printf("CSR Matrix not Loaded\n");
      return 0;
    }
    n=Paso_SystemMatrix_getTotalNumRows(A);
    b=MEMALLOC(n,double);
    x=MEMALLOC(n,double);
    for(i=0;i<n;i++)
    {
     b[i]=1;
    }
    
    Paso_test_run(A,b,1);
    
    MEMFREE(b);
    MEMFREE(x);
    Paso_SystemMatrix_free(A);
return 1;
}                                                                                                                                                                                                     
