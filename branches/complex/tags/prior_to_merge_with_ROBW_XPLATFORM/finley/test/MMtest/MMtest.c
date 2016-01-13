/*
 ************************************************************
 *          Copyright 2006 by ACcESS MNRF                   *
 *                                                          *
 *              http://www.access.edu.au                    *
 *       Primary Business: Queensland, Australia            *
 *  Licensed under the Open Software License version 3.0    *
 *     http://www.opensource.org/licenses/osl-3.0.php       *
 *                                                          *
 ************************************************************
*/


#include "finley/finleyC/System.h"
#include "finley/finleyC/Finley.h"

#include <stdio.h>

int main( int argc, char *argv[] )
{
	Finley_SystemMatrix *fsm = NULL;

	if( argc != 3 )
	{
		printf( "usage: %s infile.mm outfile.mm\n", argv[0] );
		return -1;
	}

	fsm = Finley_SystemMatrix_loadMM_toCSR( argv[1] );
	if( Finley_ErrorCode != NO_ERROR )
	{
		printf( "Error:: %s\n", Finley_ErrorMsg );
		return -1;
	}

	Finley_SystemMatrix_saveMM( fsm, argv[2] );

	return 0;
}
