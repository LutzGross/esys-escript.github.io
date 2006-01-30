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

#if !defined escript_DataProf_20050620_H
#define escript_DataProf_20050620_H

#include <iostream>
#include <sstream>
#include <string>

using namespace std;

namespace escript {

typedef struct profDataEntry {
  int interpolate;
  int grad;
  int integrate;
  int where;
  int unary;
  int binary;
  int reduction1;
  int reduction2;
  int slicing;
} profDataEntry;

/**
   \brief
   Class used for tracking profiling data within the escript Data class.

   Description:
   Give a detailed description of DataProf

   Template Parameters:
   For templates describe any conditions that the parameters used in the
   template must satisfy
*/
class DataProf {

 public:

  /**
     \brief
     Default constructor for DataProf

     Description:
     Default constructor for DataProf

     Preconditions:
     Describe any preconditions

     Throws:
     Describe any exceptions thrown
  */
  DataProf();

  /**
     \brief
     Default destructor for DataProf

     Description:
     Default destructor for DataProf

     Preconditions:
     Describe any preconditions

     Throws:
     Describe any exceptions thrown
  */
  ~DataProf();

  /**
     \brief
     Return a new table entry for tracking profiling data on
     an escript Data object.

     Description:
     Return a new table entry for tracking profiling data on
     an escript Data object.

     Preconditions:
     Describe any preconditions

     Throws:
     Describe any exceptions thrown
  */
  profDataEntry*
  newData();

  /**
     \brief
     Dump the contents of a profData entry as a formatted string.

     Description:
     Dump the contents of a profData entry as a formatted string.

     Preconditions:
     Describe any preconditions

     Throws:
     Describe any exceptions thrown
  */
  string
  dumpProf(profDataEntry* entry);

  /**
     \brief
     Compress the contents of profData table and dump as a formatted string.

     Description:
     Compress the contents of profData table and dump as a formatted string.

     Preconditions:
     Describe any preconditions

     Throws:
     Describe any exceptions thrown
  */
  string
  compProf();

 protected:

 private:

  typedef struct profDataTable {
    struct profDataEntry* data;
    struct profDataTable* next;
  } profDataTable;

  profDataTable* profDataTable_Root;

  int totalDataObjects;

};

} // end of namespace

#endif
