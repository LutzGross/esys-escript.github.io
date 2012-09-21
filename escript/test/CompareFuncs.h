
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


#include "escript/DataC.h"

/**
   \brief
   Compare the result obtained with the C interface to the given value.
*/
int compareTypeCode(struct escriptDataC* data, int typeResult);
int compareNumSamples(struct escriptDataC* data, int numDataPointsPerSample, int numSamples);
int compareIsExpanded(struct escriptDataC* data, int expanded);
int compareIsEmpty(struct escriptDataC* data, int empty);
/*int comparePointShape(struct escriptDataC* data, int rank, int* dimensions);*/
/*int compareSampleDataWrite(struct escriptDataC* data, int sampleNo, double* sampleData);*/
/*int compareSampleDataRead(struct escriptDataC* data, int tag, double* sampleData);*/
