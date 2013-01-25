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
 * A primitive test to ensure the array index type is at least as large
 * as requested. Could put in another test if it is larger.
 * An obscure compile error will result if the array index type isn't large
 * enough 
 */
static char EsysIntType_Too_Small[sizeof(EsysIntType)*8-ESYS_INT_BITS];

#endif
