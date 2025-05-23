
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF (Piro_ENABLE_NOX)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    ThyraSolverTpetra
    EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_Ifpack2
    SOURCES
    Main_ThyraSolver_Tpetra.cpp
    MockModelEval_A_Tpetra.cpp
    MockModelEval_A_Tpetra.hpp
    MatrixBased_LOWS.cpp
    MatrixBased_LOWS.hpp
    NUM_MPI_PROCS 1-4
    ARGS -v
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADDED_EXE_TARGET_NAME_OUT ThyraSolver_NAME
    ADDED_TESTS_NAMES_OUT ThyraSolverTpetra_NAME
  )

  IF (Piro_ENABLE_Tempus)
    #IF (TPL_ENABLE_MPI)
      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        AnalysisDriverTempus_MSD
        EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_ROL
        EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_Teko
        EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_Ifpack2
        EXCLUDE_IF_NOT_TRUE Ifpack2_ENABLE_Amesos2
        SOURCES
        Main_AnalysisDriver_Tempus.cpp
        MockModelEval_A_Tpetra.cpp
        MockModelEval_A_Tpetra.hpp
        MockModelEval_B_Tpetra.cpp
        MockModelEval_B_Tpetra.hpp
        MockModelEval_B_Tpetra_2_parameters.cpp
        MockModelEval_B_Tpetra_2_parameters.hpp
        MatrixBased_LOWS.cpp
        MatrixBased_LOWS.hpp
        MassSpringDamperModel.cpp
        MassSpringDamperModel.hpp
        NUM_MPI_PROCS 1
        ARGS input_Analysis_ROL_ReducedSpace_Transient_MSD.xml
        PASS_REGULAR_EXPRESSION "TEST PASSED"
        ADDED_TESTS_NAMES_OUT AnalysisDriverTempus_MSD_NAME 
      )
    #ENDIF()
    #TRIBITS_ADD_EXECUTABLE_AND_TEST(
    #  AnalysisDriverTempus
    #  EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_ROL
    #  EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_Teko
    #  EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_Ifpack2
    #  EXCLUDE_IF_NOT_TRUE Ifpack2_ENABLE_Amesos2
    #  SOURCES
    #  Main_AnalysisDriver_Tempus.cpp
    #  MockModelEval_A_Tpetra.cpp
    #  MockModelEval_A_Tpetra.hpp
    #  MockModelEval_B_Tpetra.cpp
    #  MockModelEval_B_Tpetra.hpp
    #  MockModelEval_B_Tpetra_2_parameters.cpp
    #  MockModelEval_B_Tpetra_2_parameters.hpp
    #  MatrixBased_LOWS.cpp
    #  MatrixBased_LOWS.hpp
    #  MassSpringDamperModel.cpp
    #  MassSpringDamperModel.hpp
    #  NUM_MPI_PROCS 1-4
    #  ARGS -v
    #  PASS_REGULAR_EXPRESSION "TEST PASSED"
    #  ADDED_TESTS_NAMES_OUT AnalysisDriverTempus_NAME 
    #)
  ENDIF()

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    AnalysisDriverTpetra
    EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_ROL
    EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_Teko
    EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_Ifpack2
    EXCLUDE_IF_NOT_TRUE Ifpack2_ENABLE_Amesos2
    SOURCES
    Main_AnalysisDriver_Tpetra.cpp
    MockModelEval_A_Tpetra.cpp
    MockModelEval_A_Tpetra.hpp
    MockModelEval_B_Tpetra.cpp
    MockModelEval_B_Tpetra.hpp
    MockModelEval_B_Tpetra_2_parameters.cpp
    MockModelEval_B_Tpetra_2_parameters.hpp
    MatrixBased_LOWS.cpp
    MatrixBased_LOWS.hpp
    NUM_MPI_PROCS 1-4
    ARGS -v
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADDED_TESTS_NAMES_OUT AnalysisDriverTpetra_NAME 
  )
  
  IF (ThyraSolverTpetra_NAME)
    SET(ThyraSolverTpetra_EXENAME ThyraSolverTpetra)
  ELSE()
    SET(ThyraSolverTpetra_EXENAME)
  ENDIF()

  IF (AnalysisDriverTpetra_NAME)
    SET(AnalysisDriverTpetra_EXENAME AnalysisDriverTpetra)
  ELSE()
    SET(AnalysisDriverTpetra_EXENAME)
  ENDIF()

ENDIF(Piro_ENABLE_NOX)

IF (Piro_ENABLE_Tempus)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    TempusSolver_SensitivitySinCosUnitTests
    SOURCES
    Piro_TempusSolver_SensitivitySinCosUnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
    Piro_Test_ThyraSupport.hpp
    SinCosModel.cpp
    SinCosModel.hpp
    NUM_MPI_PROCS 1-4
    STANDARD_PASS_OUTPUT
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    TempusSolver_AdjointSensitivitySinCosUnitTests
    SOURCES
    Piro_TempusSolver_AdjointSensitivitySinCosUnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
    Piro_Test_ThyraSupport.hpp
    SinCosModel.cpp
    SinCosModel.hpp
    NUM_MPI_PROCS 1-4
    STANDARD_PASS_OUTPUT
  )

ENDIF (Piro_ENABLE_Tempus)

TRIBITS_COPY_FILES_TO_BINARY_DIR(copyTestInputFiles
  DEST_FILES   
        input_Solve_NOX_ForwardSensitivities_Thyra.xml
        input_Solve_NOX_AdjointSensitivities_ImplicitAdjointME_Thyra.xml
        input_Solve_NOX_AdjointSensitivities_ExplicitAdjointME_Thyra.xml
        input_Analysis_ROL_ReducedSpace_LineSearch.xml
        input_Analysis_ROL_ReducedSpace_LineSearch_AdjointSensitivities_CheckGradients.xml
        input_Analysis_ROL_ReducedSpace_LineSearch_HessianBasedDotProduct.xml
        input_Analysis_ROL_ReducedSpace_Transient.xml
        input_Analysis_ROL_ReducedSpace_Transient_2_parameters.xml
        input_Analysis_ROL_ReducedSpace_Transient_MSD.xml
        input_Analysis_ROL_ReducedSpace_TrustRegion_HessianBasedDotProduct.xml
        input_Analysis_ROL_ReducedSpace_TrustRegion_HessianBasedDotProduct_2_parameters.xml
        input_Analysis_ROL_ReducedSpace_TrustRegion_BoundConstrained_NOXSolver.xml
        input_Analysis_ROL_ReducedSpace_TrustRegion_BoundConstrained_ExplicitAdjointME_NOXSolver.xml
        input_Analysis_ROL_FullSpace_AugmentedLagrangian_BoundConstrained.xml
        input_Tempus_BackwardEuler_SinCos.xml 
  SOURCE_DIR   ${PACKAGE_SOURCE_DIR}/test
  SOURCE_PREFIX "_"
  EXEDEPS ${ThyraSolverTpetra_EXENAME} ${AnalysisDriverTpetra_EXENAME}
)


IF (PIRO_HAVE_EPETRA_STACK)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    UnitTests
    SOURCES
    Piro_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
    MockModelEval_A.cpp
    MockModelEval_A.hpp
    MockModelEval_C.cpp
    MockModelEval_C.hpp
    MockModelEval_D.cpp
    MockModelEval_D.hpp
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Epetra_MatrixFreeOperator_UnitTests
    SOURCES
    Piro_Epetra_MatrixFreeOperator_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
    Piro_Test_EpetraSupport.hpp
    MockModelEval_A.hpp
    MockModelEval_A.cpp
    NUM_MPI_PROCS 1-4
    STANDARD_PASS_OUTPUT
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    MatrixFreeDecorator_UnitTests
    SOURCES
    Piro_MatrixFreeDecorator_UnitTests.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
    Piro_Test_ThyraSupport.hpp
    MockModelEval_A.hpp
    MockModelEval_A.cpp
    NUM_MPI_PROCS 1-4
    STANDARD_PASS_OUTPUT
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    EvalModel
    SOURCES
    Main_EvalModel.cpp
    MockModelEval_A.cpp
    MockModelEval_A.hpp
    NUM_MPI_PROCS 1-4
    ARGS -v
    PASS_REGULAR_EXPRESSION "TEST PASSED"
  )

  IF (Piro_ENABLE_NOX)

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      ThyraSolver
      SOURCES
      Main_ThyraSolver.cpp
      MockModelEval_A.cpp
      MockModelEval_A.hpp
      NUM_MPI_PROCS 1-4
      ARGS -v
      PASS_REGULAR_EXPRESSION "TEST PASSED"
      ADDED_EXE_TARGET_NAME_OUT ThyraSolver_NAME
    )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      SecondOrderIntegrator
      SOURCES
      Main_SecondOrderIntegrator.cpp
      MockModelEval_B.cpp
      MockModelEval_B.hpp
      ObserveSolution_Epetra.hpp
      NUM_MPI_PROCS 1
      ARGS -v
      PASS_REGULAR_EXPRESSION "TEST PASSED"
    )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      NOXSolver_UnitTests
      SOURCES
      Piro_NOXSolver_UnitTests.cpp
      ${TEUCHOS_STD_UNIT_TEST_MAIN}
      MockModelEval_A.cpp
      MockModelEval_A.hpp
      Piro_Test_ThyraSupport.hpp
      Piro_Test_WeakenedModelEvaluator.hpp
      Piro_Test_MockObserver.hpp
      NUM_MPI_PROCS 1-4
      STANDARD_PASS_OUTPUT
    )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      LOCASolver_UnitTests
      SOURCES
      Piro_LOCASolver_UnitTests.cpp
      ${TEUCHOS_STD_UNIT_TEST_MAIN}
      MockModelEval_A.cpp
      MockModelEval_A.hpp
      Piro_Test_ThyraSupport.hpp
      NUM_MPI_PROCS 1-4
      STANDARD_PASS_OUTPUT
    )

    IF (ThyraSolver_NAME)
     SET(ThyraSolver_EXENAME ThyraSolver)
    ELSE()
     SET(ThyraSolver_EXENAME)
    ENDIF()

  ENDIF(Piro_ENABLE_NOX)

  IF (Piro_ENABLE_Tempus)
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      TempusSolver_UnitTests
      SOURCES
      Piro_TempusSolver_UnitTests.cpp
      ${TEUCHOS_STD_UNIT_TEST_MAIN}
      Piro_Test_ThyraSupport.hpp
      MockModelEval_A.hpp
      MockModelEval_A.cpp
      NUM_MPI_PROCS 1-4
      STANDARD_PASS_OUTPUT
    )
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      TempusSolverForwardOnly_UnitTests
      SOURCES
      Piro_TempusSolverForwardOnly_UnitTests.cpp
      ${TEUCHOS_STD_UNIT_TEST_MAIN}
      Piro_Test_ThyraSupport.hpp
      MockModelEval_A.hpp
      MockModelEval_A.cpp
      NUM_MPI_PROCS 1-4
      STANDARD_PASS_OUTPUT
    )
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      TempusSolver_SensitivityUnitTests
      SOURCES
      Piro_TempusSolver_SensitivityUnitTests.cpp
      ${TEUCHOS_STD_UNIT_TEST_MAIN}
      Piro_Test_ThyraSupport.hpp
      MockModelEval_A.hpp
      MockModelEval_A.cpp
      NUM_MPI_PROCS 1
      STANDARD_PASS_OUTPUT
    )
    TRIBITS_ADD_TEST(
      TempusSolver_SensitivityUnitTests
      NUM_MPI_PROCS 4
      STANDARD_PASS_OUTPUT
    )
    
  ENDIF (Piro_ENABLE_Tempus)

  TRIBITS_COPY_FILES_TO_BINARY_DIR(copyEpetraTestInputFiles
    DEST_FILES   input_Solve_NOX_1.xml
                 input_Solve_NOX_2.xml
                 input_Solve_NOX_3.xml                 
		         input_Solve_NOX_Adjoint.xml
                 input_Solve_LOCA_1.xml
                 input_Solve_VV.xml
                 input_Solve_TR.xml
                 input_Solve_NB.xml
	             input_SGSolve.xml
	             input_problem1.xml
	             input_problem2.xml
	             input_coupled.xml
	             input_problem1_sg.xml
	             input_problem2_sg.xml
	             input_coupled_sg.xml
    SOURCE_DIR   ${PACKAGE_SOURCE_DIR}/test
    SOURCE_PREFIX "_"
    EXEDEPS ${ThyraSolver_EXENAME} UnitTests NOXSolver_UnitTests
  )

ENDIF(PIRO_HAVE_EPETRA_STACK)
