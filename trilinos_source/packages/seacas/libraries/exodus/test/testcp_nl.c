/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 * 
 * See packages/seacas/LICENSE for details
 */
/*****************************************************************************
 *
 * testcp_lrg - copy file test.exo created by testwt, normal_model to large_model
 *
 * author - Sandia National Laboratories
 *          Larry A. Schoof - Original
 *
 *
 * environment - UNIX
 *
 * entry conditions -
 *   input parameters:
 *
 * exit conditions -
 *
 * revision history -
 *
 *****************************************************************************/

#include "exodusII.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
  int exoid, exoid1, error, idum;
  int CPU_word_size, IO_word_size;

  float version;

  char *cdum = 0;

  ex_opts(EX_VERBOSE | EX_ABORT);

  /* open EXODUS II files */

  CPU_word_size = 0; /* sizeof(float) */
  IO_word_size  = 0; /* use size in file */

  exoid = ex_open("testcp.exo",   /* filename path */
                  EX_READ,        /* access mode = READ */
                  &CPU_word_size, /* CPU word size */
                  &IO_word_size,  /* IO word size */
                  &version);      /* ExodusII library version */

  printf("\nafter ex_open\n");
  if (exoid < 0) {
    exit(1);
  }

  printf("test.exo is an EXODUSII file; version %4.2f\n", version);
  printf("         CPU word size %1d\n", CPU_word_size);
  printf("         I/O word size %1d\n", IO_word_size);
  ex_inquire(exoid, EX_INQ_API_VERS, &idum, &version, cdum);
  printf("EXODUSII API; version %4.2f\n", version);

  CPU_word_size = 8; /* this really shouldn't matter for
                        the copy but tests the conversion
                        routines */
  IO_word_size = 4;

  exoid1 = ex_create("testcp_nl.exo",             /* filename */
                     EX_CLOBBER | EX_LARGE_MODEL, /* OK to overwrite, large */
                     &CPU_word_size,              /* CPU float word size in bytes */
                     &IO_word_size);              /* I/O float word size in bytes */

  printf("\nafter ex_create, exoid = %3d\n", exoid1);
  if (exoid1 < 0) {
    exit(1);
  }

  printf("         CPU word size %1d\n", CPU_word_size);
  printf("         I/O word size %1d\n", IO_word_size);

  error = ex_copy(exoid, exoid1);
  printf("\nafter ex_copy, error = %3d\n", error);

  error = ex_close(exoid);
  printf("\nafter ex_close, error = %3d\n", error);

  error = ex_close(exoid1);
  printf("\nafter ex_close, error = %3d\n", error);
  return 0;
}
