/* $Id$ */
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

#if !defined  escript_DataC_20040611_H
#define escript_DataC_20040611_H

#define ESCRIPT_MAX_DATA_RANK 4
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
int getFunctionSpaceType(escriptDataC* data);

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
int isEmpty(escriptDataC* data);

/**
   \brief
   Return true if the input shape matches the data point shape for data
   \param data Input - C wrapper for Data.
   \param rank Input - number of dimensions.
   \param dimensions Input - 
*/
int isDataPointShapeEqual(escriptDataC* data, int rank, int* dimensions);
/**
   \brief
   Return true if the number of data points per sample and the number 
   of samples equal the input values. In the case that data is empty or NULL,
   true is returned.
   \param data Input - C wrapper for Data.
   \param numDataPointsPerSample Input - number of data points per sample
   \param numSamples Input - number of samples
*/
int numSamplesEqual(escriptDataC* data, int numDataPointsPerSample,
		    int numSamples);

/**
   \brief
   Returns the rank of the point data for the data. 
   \param data Input - C wrapper for Data.
*/
int getDataPointRank(escriptDataC* data);

/**
   \brief
   Returns the value of the i-th component of the shape of the point data.
   \param data Input - C wrapper for Data.
   \param i Input - index of shape component.
*/
int getDataPointShape(escriptDataC* data,int i);

/**
   \brief
   Return the number of doubles needed for each data point.
   \param data Input - C wrapper for Data.
*/
int getDataPointSize(escriptDataC* data);
/**
   \brief
   Return the number of doubles stored for the Data object.
   Argument data may be NULL, in which case 0 is returnd.
   \param data Input - C wrapper for Data.
*/
int getLength(escriptDataC* data);
/**
   \brief
   Return true if data is expanded.
   Argument data may be NULL, in which case false is returnd.
   \param data Input - C wrapper for Data.
*/
int isExpanded(escriptDataC* data);
/**
   \brief
   Return a pointer to the data for the given sample number.
   if data is empty NULL is returned.
   data may be NULL, in which case NULL is returnd.
  \param data Input - C wrapper for Data.
  \param sampleNo Input - The sample number.
*/
double* getSampleData(escriptDataC* data, int sampleNo);

#endif
