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

#include "escript/Data/Taipan.h"

#include <iostream>
#include <cassert>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

namespace escript {

Taipan::Taipan() :
  memTable_Root(0),
  totalElements(0)
{
}

Taipan::~Taipan() {

  Taipan_MemTable *tab;
  Taipan_MemTable *tab_next;

  // deallocate all managed arrays and the memTable
  tab = memTable_Root;
  while (tab != 0) {
    tab_next = tab->next;
    delete[] tab->array;
    delete tab;
    tab = tab_next;
  }

  // clear the MemTable root node
  memTable_Root = 0;

  // clear the totalElements counter
  totalElements = -1;
}

double*
Taipan::new_array(int dim, int N) {

  assert(totalElements >= 0);

  int len = 0;
                                                                                                                                       
  #ifdef _OPENMP
  int numThreads = omp_get_num_threads();
  #else
  int numThreads = 1;
  #endif

  Taipan_MemTable *tab;
  Taipan_MemTable *new_tab;
  Taipan_MemTable *tab_prev;

//  numThreads = omp_get_max_threads();

  // is a suitable array already available?
  if (memTable_Root != 0) {
    tab = memTable_Root;
    while (tab != 0) {
      if (tab->dim == dim &&
          tab->N == N &&
          tab->free &&
          tab->numThreads == numThreads) {
        tab->free = false;
        return tab->array;
      }
      tab_prev = tab;
      tab = tab->next;
    }
  }

  //
  // otherwise a new array must be allocated

  // create the corresponding memTable entry
  len = dim * N;
  new_tab = new Taipan_MemTable;
  new_tab->dim = dim;
  new_tab->N = N;
  new_tab->numThreads = numThreads;
  new_tab->free = false;
  new_tab->next = 0;
  if (memTable_Root == 0) {
    memTable_Root = new_tab;
  } else {
    tab_prev->next = new_tab;
  }

  // allocate and initialise the new array
  new_tab->array = new double[len];
  int i,j;
  #pragma omp parallel for private(i,j) schedule(static)
  for (i=0; i<N; i++) {
    for (j=0; j<dim; j++) {
      new_tab->array[j+dim*i]=0.0;
    }
  }
  totalElements += len;

  return new_tab->array;
}

void
Taipan::delete_array(double* array) {

  assert(totalElements >= 0);

  int N;
  int len = 0;
  bool found = false;

  Taipan_MemTable *tab;
  Taipan_MemTable *tab_next;
  Taipan_MemTable *tab_prev = 0;

  if (memTable_Root != 0) {

    // find the table entry for this array and mark it as free
    tab = memTable_Root;
    while (tab != 0) {
      if (tab->array == array) {
        N = tab->N;
        tab->free = true;
        found = true;
        break;
      }
      tab = tab->next;
    }
    if (!found) {
      // this wasn't an array under management, so quit now
      return;
    }

    // are there any N block arrays still in use?
    tab = memTable_Root;
    while (tab != 0) {
      if (tab->N==N && !tab->free) {
        return;
      }
      tab = tab->next;
    }

    // if not, all N block arrays are deallocated
    tab = memTable_Root;
    while (tab != 0) {
      tab_next = tab->next;
      if (tab->N == N) {
        delete[] tab->array;
        len += tab->dim * N;
        if (tab_prev != 0) {
          tab_prev->next = tab->next;
        } else {
          memTable_Root = tab->next;
        }
        delete tab;
      } else {
        tab_prev = tab;
      }
      tab = tab_next;
    }

    totalElements -= len;

  } else {
    // what to do if no arrays under management?
  }
}

int
Taipan::num_arrays() {

  assert(totalElements >= 0);

  int num_arrays = 0;

  Taipan_MemTable *tab;

  // count all managed arrays in the memTable
  tab = memTable_Root;
  while (tab != 0) {
    num_arrays++;
    tab = tab->next;
  }

  return num_arrays;
}

int
Taipan::num_arrays(int N) {

  assert(totalElements >= 0);

  int num_arrays = 0;

  Taipan_MemTable *tab;

  // count all managed arrays of N blocks in the memTable
  tab = memTable_Root;
  while (tab != 0) {
    if (tab->N == N) {
      num_arrays++;
    }
    tab = tab->next;
  }

  return num_arrays;
}

int
Taipan::num_free(int N) {

  assert(totalElements >= 0);

  int num_free = 0;

  Taipan_MemTable *tab;

  // count all free managed arrays of N blocks in the memTable
  tab = memTable_Root;
  while (tab != 0) {
    if (tab->N == N) {
      if (tab->free) {
        num_free++;
      }
    }
    tab = tab->next;
  }

  return num_free;
}

long
Taipan::num_elements() {
  assert(totalElements >= 0);
  return totalElements;
}

}  // end of namespace
