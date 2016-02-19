
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

#ifndef __ESYS_INDEX_H__
#define __ESYS_INDEX_H__

/*    Macros for array indexing       */

/****************************************************************************/

/****************************************************************************/

#define FALSE 0
#define TRUE 1
#define UNKNOWN -1
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
#define SAMESIGN(_arg1_, _arg2_) ( ( ( (_arg1_)>=0 ) && ( (_arg2_)>=0 ) ) || ( ((_arg1_)<=0 ) && ( (_arg2_)<=0 ) ) )
#define XNOR(_a0_,_a1_) ( ( (_a0_) && (_a1_) ) || ( !(_a0_) && !(_a1_) ) )

#endif // __ESYS_INDEX_H__

