#include "finley/finleyC/System.h"
#include "finley/finleyC/Finley.h"

#include <stdio.h>

int main( int argc, char *argv[] )
{
	Finley_SystemMatrix *fsm = NULL;

	if( argc < 3 )
	{
		printf( "usage: %s infile.mm outfile.hb [extra_arg]\n", argv[0] );
		return -1;
	}

	if( argc == 4 )
		fsm = Finley_SystemMatrix_loadMM_toCSC( argv[1] );
	else
		fsm = Finley_SystemMatrix_loadMM_toCSR( argv[1] );

	if( Finley_ErrorCode != NO_ERROR )
	{
		printf( "Error:: %s\n", Finley_ErrorMsg );
		return -1;
	}

	Finley_SystemMatrix_saveHB( fsm, argv[2] );
	Finley_SystemMatrix_saveMM( fsm, "savedMM" );

	return 0;
}
