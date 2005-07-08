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
		      char* errorMess,size_t lenMess);

#endif

/*
 * $Log$
 * Revision 1.2  2005/07/08 04:07:59  jgs
 * Merge of development branch back to main trunk on 2005-07-08
 *
 * Revision 1.1.1.1.2.1  2005/06/29 02:34:58  gross
 * some changes towards 64 integers in finley
 *
 * Revision 1.1.1.1  2004/10/26 06:53:57  jgs
 * initial import of project esys2
 *
 * Revision 1.1.1.1  2004/06/24 04:00:40  johng
 * Initial version of eys using boost-python.
 *
 *
 */
