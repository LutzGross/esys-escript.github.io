/* $Id$ */
/*
 ******************************************************************************
 *                                                                            *
 *       COPYRIGHT  ACcESS 2005 -  All Rights Reserved                        *
 *                                                                            *
 * This software is the property of ACcESS. No part of this code              *
 * may be copied in any form or by any means without the expressed written    *
 * consent of ACcESS.  Copying, use or modification of this software          *
 * by any unauthorised person is illegal unless that person has a software    *
 * license agreement with ACcESS.                                             *
 *                                                                            *
 ******************************************************************************
*/

#include "bruce/Bruce/BruceFactory.h"
#include "bruce/Bruce/BruceException.h"

#include <iostream>
#include <sstream>
#include <boost/python/extract.hpp>

using namespace std;
using namespace escript;

namespace bruce {

AbstractContinuousDomain* brick(int n0,int n1,int n2,
                                double l0,double l1,double l2)
{
    int numElements[]={n0,n1,n2};
    double length[]={l0,l1,l2};

    AbstractContinuousDomain* temp=new Bruce();
    return temp;
}

AbstractContinuousDomain* rectangle(int n0,int n1,
                                    double l0, double l1)
{
    int numElements[]={n0,n1};
    double length[]={l0,l1};

    AbstractContinuousDomain* temp=new Bruce();
    return temp;
}

}  // end of namespace
