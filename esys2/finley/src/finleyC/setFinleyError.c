#include "Finley.h"
#include "setFinleyError.h"

enum Finley_ErrorCodeType Finley_ErrorCode=NO_ERROR;
char Finley_ErrorMsg[LenErrorMsg_MAX]={'\0'};

void setFinleyError(enum Finley_ErrorCodeType errorCode, char* errorMess,
		    int lenMess) 
{
  Finley_ErrorCode=errorCode;
  if (lenMess>=LenErrorMsg_MAX) {
    strncpy(Finley_ErrorMsg,errorMess,LenErrorMsg_MAX);
    /*
     * make sure the message is null terminated
     */
    Finley_ErrorMsg[LenErrorMsg_MAX-1]='\0';
  } else {
    strncpy(Finley_ErrorMsg,errorMess,lenMess);
  }
}


