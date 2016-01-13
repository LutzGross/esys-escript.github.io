/* 
 *****************************************************************************
 *                                                                           *
 *       COPYRIGHT  ACcESS  -  All Rights Reserved                           *
 *                                                                           *
 * This software is the property of ACcESS. No part of this code             *
 * may be copied in any form or by any means without the expressed written   *
 * consent of ACcESS.  Copying, use or modification of this software         *
 * by any unauthorised person is illegal unless that person has a software   *
 * license agreement with ACcESS.                                            *
 *                                                                           *
 *****************************************************************************
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
