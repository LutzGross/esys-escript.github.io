/* $Id$ */

#ifndef INC_SETFINLEYERROR
#define INC_SETFINLEYERROR

/**************************************************************/

/* setFinleyError */

/**************************************************************/

/*   Copyrights by ACcESS Australia, 2004 */
/*   Version: $Id$ */

/**************************************************************/

/**
   \brief
   This function provides a safe way of setting the finley error code
   and error message. Will not exceed the error message buffer.
*/
void setFinleyError(enum Finley_ErrorCodeType errorCode, 
		      char* errorMess,int lenMess);

#endif

/*
 * $Log$
 * Revision 1.1  2004/10/26 06:53:57  jgs
 * Initial revision
 *
 * Revision 1.1.1.1  2004/06/24 04:00:40  johng
 * Initial version of eys using boost-python.
 *
 *
 */
