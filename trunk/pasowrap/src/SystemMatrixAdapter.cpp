
/*****************************************************************************
*
* Copyright (c) 2003-2015 by University of Queensland
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

#define ESNEEDPYTHON
#include "esysUtils/first.h"


#include "SystemMatrixAdapter.h" 
#include <escript/SolverOptions.h>
#include <paso/Solver.h>

using namespace std;

namespace paso {

SystemMatrixAdapter::SystemMatrixAdapter()
{
   throw PasoException("Error - Illegal to generate default SystemMatrixAdapter.");
}

SystemMatrixAdapter::SystemMatrixAdapter(SystemMatrix_ptr system_matrix,
        int row_blocksize, const escript::FunctionSpace& row_fs,
        int column_blocksize, const escript::FunctionSpace& column_fs) :
    AbstractSystemMatrix(row_blocksize,row_fs,column_blocksize,column_fs),
    m_system_matrix(system_matrix)
{
}

SystemMatrixAdapter::~SystemMatrixAdapter()
{ 
}

SystemMatrix_ptr SystemMatrixAdapter::getPaso_SystemMatrix() const 
{
   return m_system_matrix;
}

void SystemMatrixAdapter::ypAx(escript::Data& y, escript::Data& x) const 
{
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
   SystemMatrix_MatrixVector(1., m_system_matrix, x_dp, 1., y_dp);
   checkPasoError();
}

int SystemMatrixAdapter::mapOptionToPaso(int option)
{
    switch (option) {
        case escript::SO_DEFAULT:
            return PASO_DEFAULT;

        case escript::SO_PACKAGE_MKL:
            return PASO_MKL;
        case escript::SO_PACKAGE_PASO:
            return PASO_PASO;
        case escript::SO_PACKAGE_PASTIX:
            return PASO_PASTIX;
        case escript::SO_PACKAGE_SUPER_LU:
            return PASO_SUPER_LU;
        case escript::SO_PACKAGE_TRILINOS:
            return PASO_TRILINOS;
        case escript::SO_PACKAGE_UMFPACK:
            return PASO_UMFPACK;

        case escript::SO_METHOD_BICGSTAB:
            return PASO_BICGSTAB;
        case escript::SO_METHOD_CGS:
            return PASO_CGS;
        case escript::SO_METHOD_CHOLEVSKY:
            return PASO_CHOLEVSKY;
        case escript::SO_METHOD_CR:
            return PASO_CR;
        case escript::SO_METHOD_DIRECT:
            return PASO_DIRECT;
        case escript::SO_METHOD_GMRES:
            return PASO_GMRES;
        case escript::SO_METHOD_ITERATIVE:
            return PASO_ITERATIVE;
        case escript::SO_METHOD_MINRES:
            return PASO_MINRES;
        case escript::SO_METHOD_NONLINEAR_GMRES:
            return PASO_NONLINEAR_GMRES;
        case escript::SO_METHOD_PCG:
            return PASO_PCG;
        case escript::SO_METHOD_PRES20:
            return PASO_PRES20;
        case escript::SO_METHOD_TFQMR:
            return PASO_TFQMR;

        case escript::SO_PRECONDITIONER_AMG:
            return PASO_AMG;
        case escript::SO_PRECONDITIONER_AMLI:
            return PASO_AMLI;
        case escript::SO_PRECONDITIONER_BOOMERAMG:
            return PASO_BOOMERAMG;
        case escript::SO_PRECONDITIONER_GAUSS_SEIDEL:
            return PASO_GAUSS_SEIDEL;
        case escript::SO_PRECONDITIONER_ILU0:
            return PASO_ILU0;
        case escript::SO_PRECONDITIONER_ILUT:
            return PASO_ILUT;
        case escript::SO_PRECONDITIONER_JACOBI:
            return PASO_JACOBI;
        case escript::SO_PRECONDITIONER_NONE:
            return PASO_NO_PRECONDITIONER;
        case escript::SO_PRECONDITIONER_REC_ILU:
            return PASO_REC_ILU;
        case escript::SO_PRECONDITIONER_RILU:
            return PASO_RILU;

        case escript::SO_ODESOLVER_BACKWARD_EULER:         
            return PASO_BACKWARD_EULER;
        case escript::SO_ODESOLVER_CRANK_NICOLSON:
            return PASO_CRANK_NICOLSON;
        case escript::SO_ODESOLVER_LINEAR_CRANK_NICOLSON:
            return PASO_LINEAR_CRANK_NICOLSON;

        case escript::SO_INTERPOLATION_CLASSIC:
            return PASO_CLASSIC_INTERPOLATION;
        case escript::SO_INTERPOLATION_CLASSIC_WITH_FF_COUPLING:
            return PASO_CLASSIC_INTERPOLATION_WITH_FF_COUPLING;
        case escript::SO_INTERPOLATION_DIRECT:
            return PASO_DIRECT_INTERPOLATION;

        case escript::SO_COARSENING_AGGREGATION:
            return PASO_AGGREGATION_COARSENING;
        case escript::SO_COARSENING_CIJP:
            return PASO_CIJP_COARSENING;
        case escript::SO_COARSENING_CIJP_FIXED_RANDOM:
            return PASO_CIJP_FIXED_RANDOM_COARSENING;
        case escript::SO_COARSENING_FALGOUT:
            return PASO_FALGOUT_COARSENING;
        case escript::SO_COARSENING_HMIS:
            return PASO_HMIS_COARSENING;
        case escript::SO_COARSENING_PMIS:
            return PASO_PMIS_COARSENING;
        case escript::SO_COARSENING_RUGE_STUEBEN:
            return PASO_RUGE_STUEBEN_COARSENING;
        case escript::SO_COARSENING_STANDARD:
            return PASO_STANDARD_COARSENING;   
        case escript::SO_COARSENING_YAIR_SHAPIRA:
            return PASO_YAIR_SHAPIRA_COARSENING;

        case escript::SO_REORDERING_DEFAULT:
            return PASO_DEFAULT_REORDERING;
        case escript::SO_REORDERING_MINIMUM_FILL_IN:
            return PASO_MINIMUM_FILL_IN;
        case escript::SO_REORDERING_NESTED_DISSECTION:
            return PASO_NESTED_DISSECTION;
        case escript::SO_REORDERING_NONE:
            return PASO_NO_REORDERING;

        default:
            stringstream temp;
            temp << "Error - Cannot map option value "<< option << " onto Paso";
            throw PasoException(temp.str());
    }
}

int SystemMatrixAdapter::getSystemMatrixTypeId(int solver, int preconditioner,
        int package, bool symmetry, const esysUtils::JMPI& mpiInfo)
{
    int out=SystemMatrix::getSystemMatrixTypeId(mapOptionToPaso(solver),
            mapOptionToPaso(preconditioner), mapOptionToPaso(package),
            symmetry?1:0, mpiInfo);
    checkPasoError();
    return out;
}

void SystemMatrixAdapter::Print_Matrix_Info(bool full=false) const
{
    int first_row_index = m_system_matrix->row_distribution->first_component[m_system_matrix->mpi_info->rank];
    int last_row_index  = m_system_matrix->row_distribution->first_component[m_system_matrix->mpi_info->rank+1]-1;
    int first_col_index = m_system_matrix->col_distribution->first_component[m_system_matrix->mpi_info->rank];
    int last_col_index  = m_system_matrix->col_distribution->first_component[m_system_matrix->mpi_info->rank+1]-1;

    std::cout << "Print_Matrix_Info running on CPU "
        << m_system_matrix->mpi_info->rank << " of "
        << m_system_matrix->mpi_info->size << std::endl;

    switch (m_system_matrix->type) {
        case MATRIX_FORMAT_DEFAULT:      
            std::cout << "\tMatrix type MATRIX_FORMAT_DEFAULT" << std::endl;
            break;
        case MATRIX_FORMAT_CSC:
            std::cout << "\tMatrix type MATRIX_FORMAT_CSC" << std::endl;
            break;
        case MATRIX_FORMAT_BLK1:
            std::cout << "\tMatrix type MATRIX_FORMAT_BLK1" << std::endl;
            break;
        case MATRIX_FORMAT_OFFSET1:
            std::cout << "\tMatrix type MATRIX_FORMAT_OFFSET1" << std::endl;
            break;
        case MATRIX_FORMAT_TRILINOS_CRS:
            std::cout << "\tMatrix type MATRIX_FORMAT_TRILINOS_CRS" << std::endl;
            break;
        default:
            std::cout << "\tMatrix type unknown" << std::endl;
            break;
    }

    std::cout << "\trow indices run from " << first_row_index << " to "
              << last_row_index << std::endl;
    std::cout << "\tcol indices run from " << first_col_index << " to "
              << last_col_index << std::endl;
    std::cout << "\tmainBlock numRows " << m_system_matrix->mainBlock->numRows
              << std::endl;
    std::cout << "\tmainBlock numCols " << m_system_matrix->mainBlock->numCols
              << std::endl;
    std::cout << "\tmainBlock pattern numOutput "
              << m_system_matrix->mainBlock->pattern->numOutput << std::endl;
    std::cout << "\tcol_coupleBlock numRows "
              << m_system_matrix->col_coupleBlock->numRows << std::endl;
    std::cout << "\tcol_coupleBlock numCols "
              << m_system_matrix->col_coupleBlock->numCols << std::endl;
    std::cout << "\tcol_coupleBlock pattern numOutput "
              << m_system_matrix->col_coupleBlock->pattern->numOutput
              << std::endl;
    std::cout << "\trow_coupleBlock numRows "
              << m_system_matrix->row_coupleBlock->numRows << std::endl;
    std::cout << "\trow_coupleBlock numCols "
              << m_system_matrix->row_coupleBlock->numCols << std::endl;
    std::cout << "\trow_coupleBlock pattern numOutput "
              << m_system_matrix->row_coupleBlock->pattern->numOutput
              << std::endl;
    std::cout << "\trow_block_size " << m_system_matrix->row_block_size
              << std::endl;
    std::cout << "\tcol_block_size " << m_system_matrix->col_block_size
              << std::endl;
    std::cout << "\tblock_size " << m_system_matrix->block_size << std::endl;
    std::cout << "\tlogical_row_block_size "
              << m_system_matrix->logical_row_block_size << std::endl;
    std::cout << "\tlogical_col_block_size "
              << m_system_matrix->logical_col_block_size << std::endl;
}

void SystemMatrixAdapter::setToSolution(escript::Data& out, escript::Data& in,
                                        boost::python::object& options) const
{
   Options paso_options;
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
   paso::solve(m_system_matrix, out_dp, in_dp, &paso_options);
   pasoToEscriptOptions(&paso_options,options);
   checkPasoError();
}

void SystemMatrixAdapter::nullifyRowsAndCols(escript::Data& row_q,
                                             escript::Data& col_q, double mdv)
{
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
   m_system_matrix->nullifyRowsAndCols(row_q_dp, col_q_dp, mdv);
   checkPasoError();
}

void SystemMatrixAdapter::saveMM(const std::string& filename) const
{
   m_system_matrix->saveMM(filename.c_str());
   checkPasoError();
}

void SystemMatrixAdapter::saveHB(const std::string& filename) const
{
   m_system_matrix->saveHB(filename.c_str());
   checkPasoError();
}

void SystemMatrixAdapter::resetValues()
{
   m_system_matrix->setValues(0.);
   solve_free(m_system_matrix.get());
   checkPasoError();
}

void SystemMatrixAdapter::pasoToEscriptOptions(const Options* paso_options,
                                               boost::python::object& options) 
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

void SystemMatrixAdapter::escriptToPasoOptions(Options* paso_options,
                                         const boost::python::object& options) 
{
    escript::SolverBuddy sb = boost::python::extract<escript::SolverBuddy>(options);

    paso_options->setDefaults();
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
    paso_options->ode_solver = mapOptionToPaso(sb.getODESolver());
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
}
 

}  // end of namespace

