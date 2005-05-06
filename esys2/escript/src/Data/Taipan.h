/* 
 ******************************************************************************
 *                                                                            *
 *       COPYRIGHT  ACcESS 2004 -  All Rights Reserved                        *
 *                                                                            *
 * This software is the property of ACcESS. No part of this code              *
 * may be copied in any form or by any means without the expressed written    *
 * consent of ACcESS.  Copying, use or modification of this software          *
 * by any unauthorised person is illegal unless that person has a software    *
 * license agreement with ACcESS.                                             *
 *                                                                            *
 ******************************************************************************
*/

#if !defined escript_Taipan_20050427_H
#define escript_Taipan_20050427_H

namespace escript {

/**
   \brief
   Taipan array manager, C++ version.
   Based on TaipanMemManager C module by Lutz Gross.

   Description:
   Taipan: data-array manager.

   The Taipan data-array manager holds a set of (dim x N) arrays distributed across a number of threads.
   If a (dim x N) array is requested via the Taipan allocator, the pool of managed arrays is searched for
   a free array of this size on the current number of threads. If none is available, a new one is allocated
   and added to the pool of managed arrays.

   When a managed array is deallocated, the array is marked as free but not returned to the system as long
   as at least one array with dimension N is in use. Otherwise all arrays with second dimsenion N are
   deallocated as it is assumed that arrays with second dimension N will not be used anymore.

   Template Parameters:
   For templates describe any conditions that the parameters used in the
   template must satisfy
*/

class Taipan {

 public:

  /**
     \brief
     Default constructor for Taipan data-array manager.

     Description:
     Default constructor for Taipan data-array manager.

     Preconditions:
     Describe any preconditions

     Throws:
     Describe any exceptions thrown
  */
  Taipan();

  /**
     \brief
     Default destructor for Taipan data-array manager.

     Description:
     Default destructor for Taipan data-array manager.

     Preconditions:
     Describe any preconditions

     Throws:
     Describe any exceptions thrown
  */
  ~Taipan();

  /**
     \brief
     Taipan data-array allocator.

     The parameter "dim" defines the contiguous "blocksize" within the array.
     Where the array is allocated accross multiple threads, it will be split
     on block boundaries only. N defines the number of "blocks" in the array.

     Where escript Data objects are concerned, N corresponds to the number of samples,
     and dim correesponds to the number of data-points per sample multiplied by the
     number of entries in a data-point.
  */
  double*
  new_array(int dim, int N);

  /**
     \brief
     Taipan data-array deallocator.

     The parameter "array" should be an array object that was returned by Taipan::new_array.
  */
  void
  delete_array(double* array);

  /**
     \brief
     Calculate the total number of arrays currently under management.
  */
  int
  num_arrays();

  /**
     \brief
     Calculate the total number of arrays of N blocks currently under management.
  */
  int
  num_arrays(int N);

  /**
     \brief
     Calculate the total number of free arrays of N blocks currently under management.
  */
  int
  num_free(int N);

  /**
     \brief
     Return the total number of array elements currently under management.
  */
  long
  num_elements();
 
 protected:

 private:

  typedef struct Taipan_MemTable {
    double* array;
    int dim;
    int N;
    int numThreads;
    bool free;
    struct Taipan_MemTable* next;
  } Taipan_MemTable;

  Taipan_MemTable* memTable_Root;

  long totalElements;

};

} // end of namespace

#endif
