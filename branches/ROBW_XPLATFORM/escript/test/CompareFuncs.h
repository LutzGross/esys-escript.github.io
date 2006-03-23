#include "escriptcpp/DataC.h"

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
