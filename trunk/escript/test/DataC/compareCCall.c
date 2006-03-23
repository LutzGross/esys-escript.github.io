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
#include "compareCCall.h"

/*
 * Return true if all of the calls match the expected values.
 * Return false if results don't match.
 */
int compareCCall(struct escriptDataC* data, int typeResult) 
{
  int result=0;
  result+=(getFunctionSpaceType(data)==typeResult);
  return result;
}
