/*=============================================================================

 ******************************************************************************
 *                                                                            *
 *       COPYRIGHT  ACcESS 2004 -  All Rights Reserved                        *
 *                                                                            *
 * This software is the property of ACcESS.  No part of this code             *
 * may be copied in any form or by any means without the expressed written    *
 * consent of ACcESS.  Copying, use or modification of this software          *
 * by any unauthorised person is illegal unless that person has a software    *
 * license agreement with ACcESS.                                             *
 *                                                                            *
 ****************************************************************************** 

******************************************************************************/

#ifndef ESYSTYPES_H
#define ESYSTYPES_H
/*
 * Use the integer types defined in the 1999 ISO C Standard
 * To specify a suitable Esys integer type
 */
#include <stdint.h>

#if ESYS_INT_BITS==64
typedef int64_t EsysIntType;
#else
typedef int32_t EsysIntType;
#endif

/*
 * A primative test to ensure the array index type is at least as large
 * as requested. Could put in another test if it is larger.
 * An obscure compile error will result if the array index type isn't large
 * enough 
 *
 */
static char EsysIntType_Too_Small[sizeof(EsysIntType)*8-ESYS_INT_BITS];


#endif
