#include <stdio.h>
#include <unistd.h>
#include "paso/Common.h"
#include "paso/Solver.h"
#include "paso/SystemMatrix.h"
#include "Paso_tests.h"
#include "getopt.h"
#include <math.h>

#define PI (3.141592653589793)

/*
 Usage: PasoTests -f filename [-s solver] [-p preconditioner] [-l level] [-r rhs matrix] [-c coarsening method] [-t threshold parameter for AMG coarsening] 
        filename - matrix to be loaded in CSR Matrix-Market format 
        solver   - PCG, GMRES, PRES20, TFQMR and MINRES
        preconditioner - ILU0, RILU, JACOBI, GS and AMG
        level    - options are 1,2 and 3
                   0 - default option just solves with default of specified parameters 
                   1 - test all solvers with default preconditioner
                   2 - test all preconditioners with default solver
                   3 - compare solution obtained by using AMG and Jacobi precondioners
        coarsening method - YS, RS, AGG and STD.
        threshold parameter for AMG - this is the threshold value used in AMG in courenening process. Default is 0.25. For YS and AGG, please, use 0.05.            
*/

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
    dim_t i,n,level=0;
    double *error;
    double Lsuperror;
   
    int c;
    char *filename,*solver,*prec,*rhs;
    extern char *optarg;
    extern int optopt;

    Paso_Options options;
    Paso_Options_setDefaults(&options);

    
    
    options.verbose=TRUE;

    while ((c = getopt(argc, argv, "s:p:f:r:l:t:c:h")) != -1) {
      switch(c) {
        case 's':
            solver=optarg;
            if (strcmp(solver,"PCG")==0)
                options.method=PASO_PCG;
            else if (strcmp(solver,"GMRES")==0)
                options.method=PASO_GMRES;
            else if (strcmp(solver,"PRES20")==0)    
                options.method=PASO_PRES20;
            else if (strcmp(solver,"BICGSTAB")==0)
                options.method=PASO_BICGSTAB;
            else if (strcmp(solver,"TFQMR")==0)
                options.method=PASO_TFQMR;
            else if (strcmp(solver,"MINRES")==0)
                options.method=PASO_MINRES;
            else if (strcmp(solver,"DIRECT")==0) {
                options.package=PASO_MKL;
                options.verbose=1;
            }
        break;
        case 'p':
            prec=optarg;
            if (strcmp(prec,"JACOBI")==0) 
                options.preconditioner=PASO_JACOBI;
            else if (strcmp(prec,"RILU")==0) 
                options.preconditioner=PASO_RILU;
            else if (strcmp(prec,"ILU0")==0)
                options.preconditioner=PASO_ILU0;
            else if (strcmp(prec,"GS")==0)
                options.preconditioner=PASO_GS;
            else if (strcmp(prec,"AMG")==0) {
                options.preconditioner=PASO_AMG;
            }
            else if (strcmp(prec,"AMLI")==0) {
                options.preconditioner=PASO_AMG;
            }
        break;
        case 'f':
            filename = optarg;
            A=MEMALLOC(1,Paso_SystemMatrix);
            A=Paso_SystemMatrix_loadMM_toCSR(filename);
            n=Paso_SystemMatrix_getTotalNumRows(A);
            b=MEMALLOC(n,double);
            x=MEMALLOC(n,double);
            x_ref=MEMALLOC(n,double);
            error=MEMALLOC(n,double);
            for(i=0;i<n;i++) {
             x_ref[i]=cos(i*PI/n);
            }
            Paso_SystemMatrix_MatrixVector_CSR_OFFSET0(DBLE(1),A,x_ref,DBLE(0),b);
            break;
        case 'r':
            rhs=optarg;
            if (A==NULL) {
             fprintf(stderr,"System matrix is not loaded yet.\n");
             break;
            }
            n=Paso_SystemMatrix_getTotalNumRows(A);
            Paso_RHS_loadMM_toCSR(rhs,b,n);
            break;
        case 'l':
            level=atoi(optarg);
            break;
        case 't':
            options.coarsening_threshold=atof(optarg);
        case 'c':
            prec=optarg;
            if (strcmp(prec,"RS")==0) 
                options.coarsening_method=PASO_RUGE_STUEBEN_COARSENING;
            else if (strcmp(prec,"YS")==0) 
                options.coarsening_method=PASO_YAIR_SHAPIRA_COARSENING;
            else if (strcmp(prec,"AGG")==0)
                options.coarsening_method=PASO_AGGREGATION_COARSENING;
            else if (strcmp(prec,"STD")==0)
                options.coarsening_method=PASO_STANDARD_COARSENING;
        break;
        case '?':
            printf("unknown arg %c\n", optopt);
            break;
        case 'h':
            printf("Usage: PasoTests -f filename [-s solver] [-p preconditioner] [-l level] [-r rhs vector] [-c coarsening method] [-t threshold parameter for AMG coarsening] \n");
            printf("\t filename - matrix to be loaded in CSR Matrix-Market format\n");
            printf("\t solver   - PCG, GMRES, PRES20, TFQMR and MINRES\n");
            printf("\t preconditioner - ILU0, RILU, JACOBI, GS and AMG\n");
            printf("\t level    - options are 1,2 and 3\n");
            printf("\t\t 0 - default option just solves with default of specified parameters\n");
            printf("\t\t 1 - test all solvers with default preconditioner\n");            
            printf("\t\t 2 - test all preconditioners with default solver\n");
            printf("\t\t 3 - compare solution obtained by using AMG and Jacobi precondioners\n");            
            printf("\trhs vector - right hand side vector in CSR Matrix Market format.\n");
            printf("\tcoarsening method - YS, RS, AGG and STD.\n");
            printf("\tthreshold parameter for AMG - this is the threshold value used in AMG in courenening process. Default is 0.25. For YS and AGG, please, use 0.05.\n");            
            break;
        }
    }
    
    if (A==NULL) {
      /*fprintf(stderr,"CSR Matrix not Loaded\n");*/
      return 0;
    }
   
    
    if (level==0) {
        Paso_solve(A,x,b,&options);
    }
    else if (level==3) {
        options.method=PASO_PCG;
        options.verbose=TRUE;
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
    else if (level==4) {
        options.method=PASO_PCG;
        options.verbose=TRUE;
        options.preconditioner=PASO_AMG;
        
        Paso_Solver_setPreconditioner(A,&options);
        Paso_Solver_solvePreconditioner(A,x,b);
 
 
        for(i=0;i<n;i++) {
          error[i]=x[i]-x_ref[i];
        }
        Lsuperror=Lsup(error,n)/Lsup(x_ref,n);
        fprintf(stdout,"Lsup error AMG as a solver %e\n\n",Lsuperror);
        
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

   if (A!=NULL) {
    MEMFREE(b);
    MEMFREE(x);
    MEMFREE(x_ref);
    MEMFREE(error);
    Paso_SystemMatrix_free(A);
   }
    
return 1;
}


