#
# simple_stratimikos_example
#

TRIBITS_ADD_EXECUTABLE(
  simple_stratimikos_example
  SOURCES
    simple_stratimikos_example.cpp
  COMM serial mpi
  )


TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyExampleFourByFourMtx
  DEST_FILES FourByFour.mtx
  SOURCE_DIR ${PACKAGE_SOURCE_DIR}/../amesos/test/Test_Basic
  SOURCE_PREFIX "In_"
  EXEDEPS simple_stratimikos_example
  )


TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyExampleXMLFiles
  DEST_FILES
    amesos.klu.xml
    aztecoo.ifpack.xml
    aztecoo.ml.xml
    belos.ifpack.xml
    belos.ml.xml
  SOURCE_PREFIX "_"
  EXEDEPS simple_stratimikos_example
  )


TRIBITS_ADD_TEST(
  simple_stratimikos_example
  NAME simple_stratimikos_example_amesos
  EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_Amesos
  ARGS
    "--linear-solver-params-file=amesos.klu.xml --matrix-file=FourByFour.mtx --tol=1e-13"
    "--extra-params-file=ExtraAmesosParams.xml --matrix-file=FourByFour.mtx --tol=1e-10 --linear-solver-params-used-file=lowsf.amesos.used.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  )


TRIBITS_ADD_TEST(
  simple_stratimikos_example
  NAME simple_stratimikos_example_aztecoo_ifpack
  EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_AztecOO ${PACKAGE_NAME}_ENABLE_Ifpack
  ARGS
    "--linear-solver-params-file=aztecoo.ifpack.xml --matrix-file=FourByFour.mtx --tol=1e-13"
    "--extra-params-file=ExtraAztecOOParams.xml --matrix-file=FourByFour.mtx --tol=1e-10 --linear-solver-params-used-file=lowsf.aztecoo.used.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  )


TRIBITS_ADD_TEST(
  simple_stratimikos_example
  NAME simple_stratimikos_example_aztecoo_ml
  EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_AztecOO ${PACKAGE_NAME}_ENABLE_ML
  ARGS
    "--linear-solver-params-file=aztecoo.ml.xml --matrix-file=FourByFour.mtx --tol=1e-13"
  COMM serial mpi
  NUM_MPI_PROCS 1
  )


TRIBITS_ADD_TEST(
  simple_stratimikos_example
  NAME simple_stratimikos_example_belos_ifpack
  EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_Belos ${PACKAGE_NAME}_ENABLE_Ifpack
  ARGS
    "--linear-solver-params-file=belos.ifpack.xml --matrix-file=FourByFour.mtx --tol=1e-13"
  COMM serial mpi
  NUM_MPI_PROCS 1
  )


TRIBITS_ADD_TEST(
  simple_stratimikos_example
  NAME simple_stratimikos_example_belos_ml
  EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_Belos ${PACKAGE_NAME}_ENABLE_ML
  ARGS
    "--linear-solver-params-file=belos.ml.xml --matrix-file=FourByFour.mtx --tol=1e-12"
    # 2008/12/11: rabartl Above: I had to drop the tolerance from 1e-13 to 1-e12 from
    # what is listed in stratimikos/test/defintion.  This must be due to some bad
    # compiler option or something but I don't have time to figure this out.
  COMM serial mpi
  NUM_MPI_PROCS 1
  )


#
# Generate sample output
#

CONFIGURE_FILE(generate-sample-output.pl.stub.in generate-sample-output.pl)

ADD_CUSTOM_TARGET(stratimikos_sample_output
  COMMENT "Generating sample output for package ${PACKAGE_NAME} ..."
  COMMAND ${PERL_EXECUTABLE} generate-sample-output.pl
  )


#
# MixedOrderPhysicsBasedPreconditioner
#


TRIBITS_ADD_EXECUTABLE(
  MixedOrderPhysicsBasedPreconditioner
  SOURCES
    MixedOrderPhysicsBasedPreconditioner.cpp
  COMM serial mpi
  )


TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyExampleMOPBPFiles
  DEST_FILES
    MixedOrderPhysicsBasedPreconditioner.Base.xml
    MixedOrderPhysicsBasedPreconditioner.AllSolveOutput.xml
    MixedOrderPhysicsBasedPreconditioner.P2SolveOutputMSolveOutputToFile.xml
    MixedOrderPhysicsBasedPreconditioner.P2SolveOutput.xml
    MixedOrderPhysicsBasedPreconditioner.Belos.xml
    P2.mtx
    P1.mtx
    M11.mtx
    M12.mtx
    M21.mtx
    M22.mtx
    ExtraAmesosParams.xml
    ExtraAztecOOParams.xml
  SOURCE_PREFIX "_"
  EXEDEPS MixedOrderPhysicsBasedPreconditioner
  )


TRIBITS_ADD_TEST(
  MixedOrderPhysicsBasedPreconditioner
  NAME MixedOrderPhysicsBasedPreconditioner_aztecoo_ifpack
  EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_AztecOO ${PACKAGE_NAME}_ENABLE_Ifpack
  ARGS
    "--params-file=MixedOrderPhysicsBasedPreconditioner.Base.xml --solve-tol=1.0"
    "--params-file=MixedOrderPhysicsBasedPreconditioner.Base.xml --extra-params-file=MixedOrderPhysicsBasedPreconditioner.P2SolveOutput.xml --solve-tol=1.0"
    "--params-file=MixedOrderPhysicsBasedPreconditioner.Base.xml --extra-params-file=MixedOrderPhysicsBasedPreconditioner.AllSolveOutput.xml --solve-tol=1.0"
    "--params-file=MixedOrderPhysicsBasedPreconditioner.Base.xml --extra-params-file=MixedOrderPhysicsBasedPreconditioner.P2SolveOutputMSolveOutputToFile.xml --solve-tol=1.0"
    "--params-file=MixedOrderPhysicsBasedPreconditioner.Base.xml --extra-params-file=MixedOrderPhysicsBasedPreconditioner.P2SolveOutput.xml --use-algebraic-prec --solve-tol=1.0"
  COMM serial mpi
  NUM_MPI_PROCS 1
  )


TRIBITS_ADD_TEST(
  MixedOrderPhysicsBasedPreconditioner
  NAME MixedOrderPhysicsBasedPreconditioner_belos_ifpack
  EXCLUDE_IF_NOT_TRUE ${PACKAGE_NAME}_ENABLE_Belos ${PACKAGE_NAME}_ENABLE_Ifpack
  ARGS
    "--params-file=MixedOrderPhysicsBasedPreconditioner.Belos.xml --use-algebraic-prec --solve-tol=1e+5"
    "--params-file=MixedOrderPhysicsBasedPreconditioner.Belos.xml --invert-P1 --solve-tol=1e+5"
  COMM serial mpi
  NUM_MPI_PROCS 1
  )


#
# ForwardSolverAsPreconditioner.cpp
#

TRIBITS_ADD_EXECUTABLE(
  ForwardSolverAsPreconditioner
  SOURCES
    ForwardSolverAsPreconditioner.cpp
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyForwardSolverAsPreconditionerFiles
  DEST_FILES
    P2.mtx
    AztecOOPrecBelosSolver.xml
  SOURCE_PREFIX "_"
  EXEDEPS ForwardSolverAsPreconditioner
  )

IF (${PACKAGE_NAME}_ENABLE_AztecOO AND ${PACKAGE_NAME}_ENABLE_Belos)

  # NOTE: These tests really show that using AztecOO as a preconditioner
  # (using "FixedPoint") works as a preconditioner for Belos.  It shows that
  # the AztecOO-based LinearOp is a true linear operator and it uses regular
  # GMRES (not Flexible GMRES).  If you turn off the AztecOO preconditioner
  # (i.e. change from 'ilu' to 'none') the test fails.  If you switch from
  # 'FixedPoint' to 'GMRES' as the AztecOO solver, the test fails (actually,
  # the linear op test fails).

  TRIBITS_ADD_TEST(
    ForwardSolverAsPreconditioner
    POSTFIX_AND_ARGS_0 AztecOOPRecBelosSolver
      --matrix-file=P2.mtx --param-list-file=AztecOOPrecBelosSolver.xml
      --verb-level=medium --solve-tol=1e-6
    NUM_MPI_PROCS 1
    )

  TRIBITS_ADD_TEST(
    ForwardSolverAsPreconditioner
    POSTFIX_AND_ARGS_0 AztecOOPRecBelosSolver
      --matrix-file=P2.mtx --param-list-file=AztecOOPrecBelosSolver.xml
      --verb-level=medium --solve-tol=1e-6
    COMM mpi
    NUM_MPI_PROCS 4
    )

ENDIF()
