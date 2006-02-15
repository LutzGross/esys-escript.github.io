/*
 ******************************************************************************
 *                                                                            *
 *       COPYRIGHT  ACcESS 2003,2004,2005,2006 -  All Rights Reserved              *
 *                                                                            *
 * This software is the property of ACcESS. No part of this code              *
 * may be copied in any form or by any means without the expressed written    *
 * consent of ACcESS.  Copying, use or modification of this software          *
 * by any unauthorised person is illegal unless that person has a software    *
 * license agreement with ACcESS.                                             *
 *                                                                            *
 ******************************************************************************
*/

/**************************************************************/

/*   Finley: Mesh: sets new coordinates for elements */

/**************************************************************/

/*   Author: gross@access.edu.au */
/*   Version: $Id$ */

/**************************************************************/

#include "ElementFile.h"

/**************************************************************/


void Finley_ElementFile_setCoordinates(Finley_ElementFile* self,escriptDataC* newX) {
  self->volume_is_valid=FALSE;   
  self->DSDV_is_valid=FALSE;    
  self->DSLinearDV_is_valid=FALSE; 
  self->X_is_valid=FALSE;         
}
/*
* $Log$
*/
