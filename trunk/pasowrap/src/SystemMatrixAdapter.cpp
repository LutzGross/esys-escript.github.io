
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

#include "SystemMatrixAdapter.h" 
#include <escript/SolverOptions.h>

using namespace std;

namespace paso {

struct null_deleter
{
   void operator()(void const *ptr) const
      {
      }   
};

PASOWRAP_DLL_API
SystemMatrixAdapter::SystemMatrixAdapter()
{
   throw PasoException("Error - Illegal to generate default SystemMatrixAdapter.");
}

PASOWRAP_DLL_API
SystemMatrixAdapter::SystemMatrixAdapter(Paso_SystemMatrix* system_matrix,
                                         const int row_blocksize,
                                         const escript::FunctionSpace& row_functionspace,
                                         const int column_blocksize,
                                         const escript::FunctionSpace& column_functionspace):
AbstractSystemMatrix(row_blocksize,row_functionspace,column_blocksize,column_functionspace)
{
   m_system_matrix.reset(system_matrix,null_deleter());
}

PASOWRAP_DLL_API
SystemMatrixAdapter::~SystemMatrixAdapter()
{ 
   if (m_system_matrix.unique()) {
      Paso_SystemMatrix* mat=m_system_matrix.get();
      Paso_SystemMatrix_free(mat);
   }
}

PASOWRAP_DLL_API
Paso_SystemMatrix* SystemMatrixAdapter::getPaso_SystemMatrix() const 
{
   return m_system_matrix.get();
}

PASOWRAP_DLL_API
void SystemMatrixAdapter::ypAx(escript::Data& y,escript::Data& x) const 
{
   Paso_SystemMatrix* mat=getPaso_SystemMatrix();

   if ( x.getDataPointSize()  != getColumnBlockSize()) {
      throw PasoException("matrix vector product : column block size does not match the number of components in input.");
   } else if (y.getDataPointSize() != getRowBlockSize()) {
      throw PasoException("matrix vector product : row block size does not match the number of components in output.");
   } else if ( x.getFunctionSpace()  != getColumnFunctionSpace()) {
      throw PasoException("matrix vector product : column function space and function space of input don't match.");
   } else if (y.getFunctionSpace() != getRowFunctionSpace()) {
      throw PasoException("matrix vector product : row function space and function space of output don't match.");
   }
   x.expand();
   y.expand();
   x.requireWrite();
   y.requireWrite();
   double* x_dp=x.getSampleDataRW(0);
   double* y_dp=y.getSampleDataRW(0);
   Paso_SystemMatrix_MatrixVector(1., mat,x_dp, 1.,y_dp);
   checkPasoError();
}

PASOWRAP_DLL_API
int SystemMatrixAdapter::mapOptionToPaso(const int option)  {

   switch (option) {
       case  escript::ESCRIPT_DEFAULT:
          return PASO_DEFAULT;
       case  escript::ESCRIPT_DIRECT:
          return PASO_DIRECT;
       case  escript::ESCRIPT_CHOLEVSKY:
          return PASO_CHOLEVSKY;
       case  escript::ESCRIPT_PCG:
          return PASO_PCG;
       case  escript::ESCRIPT_CR:
          return PASO_CR;
       case  escript::ESCRIPT_CGS:
          return PASO_CGS;
       case  escript::ESCRIPT_BICGSTAB:
          return PASO_BICGSTAB;
       case  escript::ESCRIPT_ILU0:
          return PASO_ILU0;
       case  escript::ESCRIPT_ILUT:
          return PASO_ILUT;
       case  escript::ESCRIPT_JACOBI:
          return PASO_JACOBI;
       case  escript::ESCRIPT_GMRES:
          return PASO_GMRES;
       case  escript::ESCRIPT_PRES20:
          return PASO_PRES20;
       case  escript::ESCRIPT_LUMPING:
          return PASO_LUMPING;
       case  escript::ESCRIPT_NO_REORDERING:
          return PASO_NO_REORDERING;
       case  escript::ESCRIPT_MINIMUM_FILL_IN:
          return PASO_MINIMUM_FILL_IN;
       case  escript::ESCRIPT_NESTED_DISSECTION:
          return PASO_NESTED_DISSECTION;
       case  escript::ESCRIPT_MKL:
          return PASO_MKL;
       case  escript::ESCRIPT_UMFPACK:
          return PASO_UMFPACK;
       case  escript::ESCRIPT_ITERATIVE:
          return PASO_ITERATIVE;
       case  escript::ESCRIPT_PASO:
          return PASO_PASO;
       case  escript::ESCRIPT_AMG:
          return PASO_AMG;
       case  escript::ESCRIPT_AMLI:
          return PASO_AMLI;
       case  escript::ESCRIPT_REC_ILU:
          return PASO_REC_ILU;
       case  escript::ESCRIPT_TRILINOS:
          return PASO_TRILINOS;
       case  escript::ESCRIPT_NONLINEAR_GMRES:
          return PASO_NONLINEAR_GMRES;
       case  escript::ESCRIPT_TFQMR :
          return PASO_TFQMR;
       case  escript::ESCRIPT_MINRES:
          return PASO_MINRES;
       case  escript::ESCRIPT_GAUSS_SEIDEL:
          return PASO_GAUSS_SEIDEL;
       case  escript::ESCRIPT_RILU:
          return PASO_RILU;
       case  escript::ESCRIPT_DEFAULT_REORDERING:
          return PASO_DEFAULT_REORDERING;
       case  escript::ESCRIPT_SUPER_LU:
          return PASO_SUPER_LU;
       case  escript::ESCRIPT_PASTIX:
          return PASO_PASTIX;
       case  escript::ESCRIPT_YAIR_SHAPIRA_COARSENING:
          return PASO_YAIR_SHAPIRA_COARSENING;
       case  escript::ESCRIPT_RUGE_STUEBEN_COARSENING:
          return PASO_RUGE_STUEBEN_COARSENING;
       case  escript::ESCRIPT_STANDARD_COARSENING:
          return PASO_STANDARD_COARSENING;   
       case  escript::ESCRIPT_AGGREGATION_COARSENING:
          return PASO_AGGREGATION_COARSENING;
       case  escript::ESCRIPT_NO_PRECONDITIONER:
          return PASO_NO_PRECONDITIONER;
       case escript::ESCRIPT_CLASSIC_INTERPOLATION_WITH_FF_COUPLING:
         return PASO_CLASSIC_INTERPOLATION_WITH_FF_COUPLING;
       case escript::ESCRIPT_CLASSIC_INTERPOLATION:
         return PASO_CLASSIC_INTERPOLATION;
       case escript::ESCRIPT_DIRECT_INTERPOLATION:
         return PASO_DIRECT_INTERPOLATION;
       case escript::ESCRIPT_BOOMERAMG:
         return PASO_BOOMERAMG;
       case escript::ESCRIPT_CIJP_FIXED_RANDOM_COARSENING:
         return PASO_CIJP_FIXED_RANDOM_COARSENING;
       case escript::ESCRIPT_CIJP_COARSENING:
         return PASO_CIJP_COARSENING;
       case escript::ESCRIPT_FALGOUT_COARSENING:
         return PASO_FALGOUT_COARSENING;
       case escript::ESCRIPT_PMIS_COARSENING:
         return PASO_PMIS_COARSENING;
       case escript::ESCRIPT_HMIS_COARSENING:
         return PASO_HMIS_COARSENING;
       case  escript::ESCRIPT_LINEAR_CRANK_NICOLSON:
         return PASO_LINEAR_CRANK_NICOLSON;
       case  escript::ESCRIPT_CRANK_NICOLSON:
         return PASO_CRANK_NICOLSON;
       case  escript::ESCRIPT_BACKWARD_EULER:         
         return PASO_BACKWARD_EULER;
       default:
           stringstream temp;
           temp << "Error - Cannot map option value "<< option << " onto Paso";
           throw PasoException(temp.str());
    }
}

PASOWRAP_DLL_API
int SystemMatrixAdapter::getSystemMatrixTypeId(const int solver, const int preconditioner,
        const int package, const bool symmetry, Esys_MPIInfo* mpiInfo)
{
    int out=Paso_SystemMatrix_getSystemMatrixTypeId(mapOptionToPaso(solver),
            mapOptionToPaso(preconditioner), mapOptionToPaso(package),
            symmetry?1:0, mpiInfo);
    checkPasoError();
    return out;
}

PASOWRAP_DLL_API
void SystemMatrixAdapter::Print_Matrix_Info(const bool full=false) const
{
    Paso_SystemMatrix* mat=m_system_matrix.get();
    int first_row_index  = mat->row_distribution->first_component[mat->mpi_info->rank];
    int last_row_index   = mat->row_distribution->first_component[mat->mpi_info->rank+1]-1;
    int first_col_index  = mat->col_distribution->first_component[mat->mpi_info->rank];
    int last_col_index   = mat->col_distribution->first_component[mat->mpi_info->rank+1]-1;

    fprintf(stdout, "Print_Matrix_Info running on CPU %d of %d\n", mat->mpi_info->rank, mat->mpi_info->size);

    switch (mat->type) {
        case MATRIX_FORMAT_DEFAULT:      
            fprintf(stdout, "\tMatrix type MATRIX_FORMAT_DEFAULT\n"); 
            break;
        case MATRIX_FORMAT_CSC:
            fprintf(stdout, "\tMatrix type MATRIX_FORMAT_CSC\n"); 
            break;
        case MATRIX_FORMAT_BLK1:
            fprintf(stdout, "\tMatrix type MATRIX_FORMAT_BLK1\n");
            break;
        case MATRIX_FORMAT_OFFSET1:
            fprintf(stdout, "\tMatrix type MATRIX_FORMAT_OFFSET1\n"); 
            break;
        case MATRIX_FORMAT_TRILINOS_CRS:
            fprintf(stdout, "\tMatrix type MATRIX_FORMAT_TRILINOS_CRS\n"); 
            break;
        default:
            fprintf(stdout, "\tMatrix type unknown\n"); 
            break;
   }

   fprintf(stdout, "\trow indices run from %d to %d\n", first_row_index, last_row_index);
   fprintf(stdout, "\tcol indices run from %d to %d\n", first_col_index, last_col_index);
   fprintf(stdout, "\tmainBlock numRows %d\n", mat->mainBlock->numRows);
   fprintf(stdout, "\tmainBlock numCols %d\n", mat->mainBlock->numCols);
   fprintf(stdout, "\tmainBlock pattern numOutput %d\n", mat->mainBlock->pattern->numOutput);
   fprintf(stdout, "\tcol_coupleBlock numRows %d\n", mat->col_coupleBlock->numRows);
   fprintf(stdout, "\tcol_coupleBlock numCols %d\n", mat->col_coupleBlock->numCols);
   fprintf(stdout, "\tcol_coupleBlock pattern numOutput %d\n", mat->col_coupleBlock->pattern->numOutput);
   fprintf(stdout, "\trow_coupleBlock numRows %d\n", mat->row_coupleBlock->numRows);
   fprintf(stdout, "\trow_coupleBlock numCols %d\n", mat->row_coupleBlock->numCols);
   fprintf(stdout, "\trow_coupleBlock pattern numOutput %d\n", mat->row_coupleBlock->pattern->numOutput);
   fprintf(stdout, "\trow_block_size %d\n", mat->row_block_size);
   fprintf(stdout, "\tcol_block_size %d\n", mat->col_block_size);
   fprintf(stdout, "\tblock_size %d\n", mat->block_size);
   fprintf(stdout, "\tlogical_row_block_size %d\n", mat->logical_row_block_size);
   fprintf(stdout, "\tlogical_col_block_size %d\n", mat->logical_col_block_size);

}

PASOWRAP_DLL_API
void SystemMatrixAdapter::setToSolution(escript::Data& out,escript::Data& in, boost::python::object& options) const
{
   Paso_SystemMatrix* mat=getPaso_SystemMatrix();
   Paso_Options paso_options;
   options.attr("resetDiagnostics")();
   escriptToPasoOptions(&paso_options,options);
   if ( out.getDataPointSize()  != getColumnBlockSize()) {
      throw PasoException("solve : column block size does not match the number of components of solution.");
   } else if ( in.getDataPointSize() != getRowBlockSize()) {
      throw PasoException("solve : row block size does not match the number of components of  right hand side.");
   } else if ( out.getFunctionSpace()  != getColumnFunctionSpace()) {
      throw PasoException("solve : column function space and function space of solution don't match.");
   } else if (in.getFunctionSpace() != getRowFunctionSpace()) {
      throw PasoException("solve : row function space and function space of right hand side don't match.");
   }
   out.expand();
   in.expand();
   double* out_dp=out.getSampleDataRW(0);        
   double* in_dp=in.getSampleDataRW(0);                
   Paso_solve(mat,out_dp,in_dp,&paso_options);
   pasoToEscriptOptions(&paso_options,options);
   checkPasoError();

}

PASOWRAP_DLL_API
void SystemMatrixAdapter::nullifyRowsAndCols(escript::Data& row_q,escript::Data& col_q, const double mdv) const
{
   Paso_SystemMatrix* mat = getPaso_SystemMatrix();
   if ( col_q.getDataPointSize()  != getColumnBlockSize()) {
      throw PasoException("nullifyRowsAndCols : column block size does not match the number of components of column mask.");
   } else if ( row_q.getDataPointSize() != getRowBlockSize()) {
      throw PasoException("nullifyRowsAndCols : row block size does not match the number of components of row mask.");
   } else if ( col_q.getFunctionSpace()  != getColumnFunctionSpace()) {
      throw PasoException("nullifyRowsAndCols : column function space and function space of column mask don't match.");
   } else if (row_q.getFunctionSpace() != getRowFunctionSpace()) {
      throw PasoException("nullifyRowsAndCols : row function space and function space of row mask don't match.");
   }
   row_q.expand();
   col_q.expand();
   row_q.requireWrite();
   col_q.requireWrite();
   double* row_q_dp=row_q.getSampleDataRW(0);
   double* col_q_dp=col_q.getSampleDataRW(0);
   Paso_SystemMatrix_nullifyRowsAndCols(mat,row_q_dp,col_q_dp, mdv);
   checkPasoError();
}

PASOWRAP_DLL_API
void SystemMatrixAdapter::saveMM(const std::string& fileName) const
{
   if( fileName.size() == 0 )
   {
      throw PasoException("Null file name!");
   }

   char *fName = TMPMEMALLOC(fileName.size()+1,char);
        
   strcpy(fName,fileName.c_str());
   Paso_SystemMatrix* mat = getPaso_SystemMatrix();
   Paso_SystemMatrix_saveMM(mat,fName);
   checkPasoError();
   TMPMEMFREE(fName);

}

PASOWRAP_DLL_API
void SystemMatrixAdapter::saveHB(const std::string& fileName) const
{
   if( fileName.size() == 0 )
   {
      throw PasoException("Null file name!");
   }

   char *fName = TMPMEMALLOC(fileName.size()+1,char);

   strcpy(fName,fileName.c_str());
   Paso_SystemMatrix* mat = getPaso_SystemMatrix();
   Paso_SystemMatrix_saveHB(mat,fName);
   checkPasoError();
   TMPMEMFREE(fName);

}

PASOWRAP_DLL_API
void SystemMatrixAdapter::resetValues() const
{
   Paso_SystemMatrix* mat = getPaso_SystemMatrix();
   Paso_SystemMatrix_setValues(mat,0.);
   Paso_solve_free(mat);
   checkPasoError();
}

PASOWRAP_DLL_API
void SystemMatrixAdapter::pasoToEscriptOptions(const Paso_Options* paso_options,boost::python::object& options) 
{
#define SET(__key__,__val__,__type__) options.attr("_updateDiagnostics")(__key__,(__type__)paso_options->__val__)
   SET("num_iter", num_iter, int);
   SET("num_level", num_level, int);
   SET("num_inner_iter", num_inner_iter, int);
   SET("time", time, double);
   SET("set_up_time", set_up_time, double);
   SET("net_time", net_time, double);
   SET("residual_norm", residual_norm, double);
   SET("converged",converged, bool);
   SET("time_step_backtracking_used", time_step_backtracking_used,bool);
   SET("coarse_level_sparsity",coarse_level_sparsity,double);
   SET("num_coarse_unknowns",num_coarse_unknowns,int);
#undef SET
}

PASOWRAP_DLL_API
void SystemMatrixAdapter::escriptToPasoOptions(Paso_Options* paso_options, const boost::python::object& options) 
{

    escript::SolverBuddy sb = boost::python::extract<escript::SolverBuddy>(options);

    Paso_Options_setDefaults(paso_options);
    paso_options->method = mapOptionToPaso(sb.getSolverMethod());
    paso_options->package = mapOptionToPaso(sb.getPackage());
    paso_options->verbose = sb.isVerbose();
    paso_options->symmetric = sb.isSymmetric();
    paso_options->tolerance = sb.getTolerance();
    paso_options->absolute_tolerance = sb.getAbsoluteTolerance();
    paso_options->inner_tolerance = sb.getInnerTolerance();
    paso_options->adapt_inner_tolerance = sb.adaptInnerTolerance();
    paso_options->reordering = mapOptionToPaso(sb.getReordering());
    paso_options->preconditioner = mapOptionToPaso(sb.getPreconditioner());
    paso_options->iter_max = sb.getIterMax();
    paso_options->inner_iter_max = sb.getInnerIterMax();
    paso_options->drop_tolerance = sb.getDropTolerance();
    paso_options->drop_storage = sb.getDropStorage();
    paso_options->truncation = sb.getTruncation();
    paso_options->restart = sb._getRestartForC();
    paso_options->sweeps = sb.getNumSweeps();
    paso_options->pre_sweeps = sb.getNumPreSweeps();
    paso_options->post_sweeps = sb.getNumPostSweeps();
    paso_options->level_max = sb.getLevelMax();
    paso_options->min_coarse_matrix_size = sb.getMinCoarseMatrixSize();
    paso_options->coarsening_threshold = sb.getCoarseningThreshold();
    paso_options->accept_failed_convergence = sb.acceptConvergenceFailure();
    paso_options->coarsening_method = mapOptionToPaso(sb.getCoarsening());
    paso_options->smoother = mapOptionToPaso(sb.getSmoother());
    paso_options->relaxation_factor = sb.getRelaxationFactor();
    paso_options->use_local_preconditioner = sb.useLocalPreconditioner();
    paso_options->min_coarse_sparsity = sb.getMinCoarseMatrixSparsity();
    paso_options->refinements = sb.getNumRefinements();
    paso_options->coarse_matrix_refinements = sb.getNumCoarseMatrixRefinements();
    paso_options->usePanel = sb.usePanel();
    paso_options->interpolation_method = sb.getAMGInterpolation();
    paso_options->diagonal_dominance_threshold = sb.getDiagonalDominanceThreshold();
    paso_options->ode_solver = sb.getODESolver();
}
 

}  // end of namespace
