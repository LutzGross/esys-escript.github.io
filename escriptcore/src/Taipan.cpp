
/*****************************************************************************
*
* Copyright (c) 2003-2020 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/


#include "Taipan.h"

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
  // create and initialise a new StatTable
  statTable = new Taipan_StatTable;
  clear_stats();
}

Taipan::~Taipan() {

  long len=0;
  Taipan_MemTable *tab;
  Taipan_MemTable *tab_next;

  // dump memory usage statistics
  dump_stats();

  // deallocate StatTable object
  delete statTable;

  // deallocate all managed arrays and the memTable
  tab = memTable_Root;
  while (tab != 0) {
    tab_next = tab->next;
    len = tab->dim * tab->N;
    totalElements -= len;
    delete[] tab->array;
    delete tab;
    tab = tab_next;
  }

  assert(totalElements == 0);

  // clear the MemTable root node
  memTable_Root = 0;

  // reset totalElements counter
  totalElements = -1;
}

void
Taipan::release_unused_arrays()
{
  long len=0;
  Taipan_MemTable *tab;
  Taipan_MemTable *tab_next, *tab_prev=0;
  tab = memTable_Root;
  while (tab != 0) {
      tab_next = tab->next;
      if (tab->free) {
        delete[] tab->array;
        len += tab->dim * tab->N;
        if (tab_prev != 0) {
          tab_prev->next = tab->next;
        } else {
          memTable_Root = tab->next;
        }
        delete tab;
        // increment count of arrays deallocated
        statTable->deallocations++;
      } else {
        tab_prev = tab;
      }
      tab = tab_next;
  }
  totalElements -= len;
  statTable->deallocated_elements += len;
  cout << static_cast<double>(len*sizeof(double))/1048576 << " Mbytes unused memory has been released." << endl;
}


double*
Taipan::new_array(size_type dim, size_type N) {

  assert(totalElements >= 0);

  size_type len = 0;
  #ifdef _OPENMP
  int numThreads = omp_get_num_threads();
  #else
  int numThreads = 1;
  #endif

  Taipan_MemTable *tab;
  Taipan_MemTable *new_tab;
  Taipan_MemTable *tab_prev=0;

  // increment count of alloc operations called
  statTable->requests++;

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

  try
  {
     // allocate and initialise the new array
     new_tab->array = new double[len];
  }
  catch (...)
  {
     cerr << "Memory manager failed to create array of size " << len << " doubles" << endl;
     throw;
  }
  size_type i,j;
  if (N==1) {
    for (j=0; j<dim; j++) 
      new_tab->array[j]=0.0;
  } else if (N>1) {
    #pragma omp parallel for private(i,j) schedule(static)
    for (i=0; i<N; i++) {
      for (j=0; j<dim; j++) 
        new_tab->array[j+dim*i]=0.0;
    }
  }

  totalElements += len;

  // update maximum table size
  statTable->max_tab_size = (statTable->max_tab_size < totalElements) ? totalElements : statTable->max_tab_size;

  // increment count of arrays allocated
  statTable->allocations++;

  // increment count of elements allocated
  statTable->allocated_elements += len;

  return new_tab->array;
}

void
Taipan::delete_array(double* array) {

  assert(totalElements >= 0);

  size_type N;
  size_type len = 0;
  bool found = false;

  Taipan_MemTable *tab;
  Taipan_MemTable *tab_next;
  Taipan_MemTable *tab_prev = 0;

  // increment count of free operations called
  statTable->frees++;

  if (array == 0) {
    // have been given an empty array, so quit now
    return;
  }

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

    if (N<=1) {
      // we never deallocate arrays with N<=1, so quit now
      return;
    }

    // are there any N block arrays still in use?
    tab = memTable_Root;
    while (tab != 0) {
      if (tab->N==N && !tab->free) 
        return;
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
        // increment count of arrays deallocated
        statTable->deallocations++;
      } else {
        tab_prev = tab;
      }
      tab = tab_next;
    }
   
    totalElements -= len;

    // increment count of elements deallocated
    statTable->deallocated_elements += len;

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
Taipan::num_arrays(size_type N) {

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
Taipan::num_free(size_type N) {

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

void
Taipan::dump_stats() {

  assert(totalElements >= 0);
#ifdef TAIPAN_STATS
  double elMb=statTable->allocated_elements*8.0/1048576;
  double deelMb=statTable->deallocated_elements*8.0/1048576;
  double tszMb=statTable->max_tab_size*8.0/1048576;

  cout << "======= escript Mem Stats ===========================" << endl;
  cout << "Total Num requests:             " << statTable->requests << endl;
  cout << "Total Num releases:             " << statTable->frees << endl;
  cout << "Total Num allocated arrays:     " << statTable->allocations << endl;
  cout << "Total Num deallocated arrays:   " << statTable->deallocations << endl;
  cout << "Total Num allocated elements:   " << statTable->allocated_elements << " (" << elMb << " Mb)" << endl;
  cout << "Total Num deallocated elements: " << statTable->deallocated_elements << " (" << deelMb << " Mb)" << endl;
  cout << "Maximum memory buffer size:     " << statTable->max_tab_size << " (" << tszMb << " Mb)" << endl;
  cout << "Curr Num arrays:                " << num_arrays() << endl;
  cout << "Curr Num elements in buffer:    " << num_elements() << endl;
  cout << "==================================================" << endl;
#endif
}

void
Taipan::clear_stats() {

  assert(totalElements >= 0);

  statTable->requests=0;
  statTable->frees=0;
  statTable->allocations=0;
  statTable->deallocations=0;
  statTable->allocated_elements=0;
  statTable->deallocated_elements=0;
  statTable->max_tab_size=0;
}

}  // end of namespace
