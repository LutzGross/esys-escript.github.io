
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


#if !defined  escript_DataC_20040611_H
#define escript_DataC_20040611_H
#include "system_dep.h"

/**
   \brief
   Provide a wrapper around a Data object so Data may be accessed from C.

   Description:
   Provide a wrapper around a Data object so Data may be accessed from C.

*/
struct escriptDataC {
  void* m_dataPtr;
};

typedef struct escriptDataC escriptDataC; 

/**
   \brief
   Return the function space type code.
   \param data Input - C wrapper for Data.
*/
ESCRIPT_DLL_API int getFunctionSpaceType(escriptDataC* data);

/**
   \brief
   sets the int variable _fs_ the function space type of _data_ if the data are not empty.
   \param _fs_ Input/Output - variable to be updated.
   \param _data_ Input - C wrapper for Data.
*/
#define updateFunctionSpaceType(_fs_,_data_) _fs_=(isEmpty(_data_) ? _fs_ : getFunctionSpaceType(_data_))
                                                                                     
/**
   \brief
   is true if the function space type of _data_ is equal to _fs_ or is empty
   \param _fs_ Input - function space type to checked against
   \param _data_ Input - C wrapper for Data.
*/
#define functionSpaceTypeEqual(_fs_,_data_) ( (isEmpty(_data_) || _fs_==getFunctionSpaceType(_data_)) ) ? 1 : 0

/**
   \brief
   Returns the true if the data are empty or data is NULL.
   \param data Input - C wrapper for Data.
*/
ESCRIPT_DLL_API int isEmpty(escriptDataC* data);

/**
   \brief
   Return true if the input shape matches the data point shape for data
   \param data Input - C wrapper for Data.
   \param rank Input - number of dimensions.
   \param dimensions Input - 
*/
ESCRIPT_DLL_API int isDataPointShapeEqual(escriptDataC* data, int rank, int* dimensions);
/**
   \brief
   Return true if the number of data points per sample and the number 
   of samples equal the input values. In the case that data is empty or NULL,
   true is returned.
   \param data Input - C wrapper for Data.
   \param numDataPointsPerSample Input - number of data points per sample
   \param numSamples Input - number of samples
*/
ESCRIPT_DLL_API int numSamplesEqual(escriptDataC* data, int numDataPointsPerSample,
		    int numSamples);

/**
   \brief
   Returns the number of data points per sample
   \param data Input - C wrapper for Data.
*/
ESCRIPT_DLL_API int getNumDataPointsPerSample(escriptDataC* data);

/**
   \brief
   Returns the rank of the point data for the data. 
   \param data Input - C wrapper for Data.
*/
ESCRIPT_DLL_API int getDataPointRank(escriptDataC* data);

/**
   \brief
   Returns the value of the i-th component of the shape of the point data.
   \param data Input - C wrapper for Data.
   \param i Input - index of shape component.
*/
ESCRIPT_DLL_API int getDataPointShape(escriptDataC* data,int i);

/**
   \brief
   Return the number of doubles needed for each data point.
   \param data Input - C wrapper for Data.
*/
ESCRIPT_DLL_API int getDataPointSize(escriptDataC* data);

/*
   \brief
   Return the number of doubles stored for the Data object.
   Argument data may be NULL, in which case 0 is returnd.
   \param data Input - C wrapper for Data.

This function has been removed because it does not make sense for LazyData
*/
/*ESCRIPT_DLL_API int getLength(escriptDataC* data);*/

/**
   \brief
   Return true if data can be treated as expanded.
   
   Argument data may be NULL, in which case false is returnd.
   \param data Input - C wrapper for Data.
   \return true if data is expanded or the data is lazy but would resolve to expanded. False otherwise.
*/
ESCRIPT_DLL_API int isExpanded(escriptDataC* data);

/**
   \brief
   Return a pointer to the data for the given sample number.
   if data is empty NULL is returned.
   data may be NULL, in which case NULL is returnd.
  \param data Input - C wrapper for Data.
  \param sampleNo Input - The sample number.

*/
ESCRIPT_DLL_API double __const * getSampleDataRO(escriptDataC* data, int sampleNo);
/* Placement of __const might be important. See .cpp */


ESCRIPT_DLL_API double* getSampleDataRW(escriptDataC* data, int sampleNo);


/**
   \brief
   Return a pointer to the data for the given sample number.
   Fast version of getSampledataRO: does no error checking.
  \param data Input - C wrapper for Data.
  \param sampleNo Input - The sample number.

*/
ESCRIPT_DLL_API double __const* getSampleDataROFast(escriptDataC* data, int sampleNo);

/**
   \brief
   Return a pointer to the data for the given sample number.
   Fast version of getSampledataRW: does no error checking.
  \param data Input - C wrapper for Data.
  \param sampleNo Input - The sample number.
*/
ESCRIPT_DLL_API double* getSampleDataRWFast(escriptDataC* data, int sampleNo);


/**
   \brief
   Return getSampleDataRWFast(escriptDataC* data, 0) if there are samples.
   if not, returns NULL.
   \warning This function calls requireWrite if there are samples so do not use in parallel sections.
   \warning Please do not use this in new code.
  \param data Input - C wrapper for Data.
*/
ESCRIPT_DLL_API double* getDataRW(escriptDataC* data);


/**
   \brief Ensure that this object is ready for writing.
   It will be resolved and copied if it is currently shared.
   Use only in single threaded sections of code.
   Do not create new Data objects based on this one between this call and 
   writing to the object.
*/
ESCRIPT_DLL_API void requireWrite(escriptDataC* data);

#endif
