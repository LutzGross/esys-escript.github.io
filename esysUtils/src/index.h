
/*****************************************************************************
*
* Copyright (c) 2003-2012 by University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/


#ifndef INC_ESYS_INDEX
#define INC_ESYS_INDEX

/************************************************************************************/

/*    Macros for array indexing       */

/************************************************************************************/

/************************************************************************************/

/*   some useful functions: */

#include <limits.h>


#define FALSE 0
#define TRUE 1
#define UNKNOWN -1
#define DBLE(_x_) (double)(_x_)
#define INDEX1(_X1_) (_X1_)
#define INDEX2(_X1_,_X2_,_N1_) ((_X1_)+(_N1_)*(_X2_))
#define INDEX3(_X1_,_X2_,_X3_,_N1_,_N2_) ((_X1_)+(_N1_)*INDEX2(_X2_,_X3_,_N2_))
#define INDEX4(_X1_,_X2_,_X3_,_X4_,_N1_,_N2_,_N3_) ((_X1_)+(_N1_)*INDEX3(_X2_,_X3_,_X4_,_N2_,_N3_))
#define INDEX5(_X1_,_X2_,_X3_,_X4_,_X5_,_N1_,_N2_,_N3_,_N4_) ((_X1_)+(_N1_)*INDEX4(_X2_,_X3_,_X4_,_X5_,_N2_,_N3_,_N4_))
#define INDEX6(_X1_,_X2_,_X3_,_X4_,_X5_,_X6_,_N1_,_N2_,_N3_,_N4_,_N5_) ((_X1_)+(_N1_)*INDEX5(_X2_,_X3_,_X4_,_X5_,_X6_,_N2_,_N3_,_N4_,_N5_))

#define MAX(_arg1_,_arg2_) ((_arg1_)>(_arg2_) ?  (_arg1_) : (_arg2_))
#define MAX3(_arg1_,_arg2_,_arg3_) MAX(_arg1_,MAX(_arg2_,_arg3_))
#define MIN(_arg1_,_arg2_) ((_arg1_)>(_arg2_) ?  (_arg2_) : (_arg1_)) 
#define MIN3(_arg1_,_arg2_,_arg3_) MIN(_arg1_,MIN(_arg2_,_arg3_))
#define ABS(_arg_) MAX((_arg_),-(_arg_))
#define SIGN(_arg_) ((_arg_)>0 ?  1  : ((_arg_)<0 ? -1 : 0 ))
#define SAMESIGN(_arg1_, _arg2_) ( ( ( (_arg1_)>=0 ) && ( (_arg2_)>=0 ) ) || ( ((_arg1_)<=0 ) && ( (_arg2_)<=0 ) ) )
#define SWAP(_a0_,_a1_,_type_) { \
                                _type_ s; \
                                s=(_a0_); \
                                _a0_= (_a1_); \
                                _a1_=s; \
                               }
#define XNOR(_a0_,_a1_) ( ( (_a0_) && (_a1_) ) || ( !(_a0_) && !(_a1_) ) )

#define INDEX_T_MAX INT_MAX
#define INDEX_T_MIN -INT_MAX

#endif 
