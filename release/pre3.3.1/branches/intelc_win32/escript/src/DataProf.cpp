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

#include "DataProf.h"

#include <iostream>
#include <sstream>

using namespace std;

namespace escript {

DataProf::DataProf() {
  profDataTable_Root = 0;
  totalDataObjects = 0;
}

DataProf::~DataProf() {

  profDataTable* tempTable = profDataTable_Root;
  profDataTable* lastTable;

  cout << compProf() << endl;

  while (tempTable != 0) {
    lastTable = tempTable;
    tempTable = tempTable->next;
    delete lastTable->data;
    delete lastTable;
  }

  totalDataObjects = -1;

}

profDataEntry*
DataProf::newData(){

  profDataTable* tempTable;
  profDataEntry* tempEntry;

  // create the new table and entry
  tempTable = new profDataTable;
  tempEntry = new profDataEntry;

  // add the new table to the beginning of the list
  tempTable->next = profDataTable_Root;
  profDataTable_Root = tempTable;

  // add the new entry to the new table
  tempTable->data = tempEntry;

  // initialise the new entry
  tempEntry->interpolate = 0;
  tempEntry->grad = 0;
  tempEntry->integrate = 0;
  tempEntry->where = 0;
  tempEntry->unary = 0;
  tempEntry->binary = 0;
  tempEntry->reduction1 = 0;
  tempEntry->reduction2 = 0;
  tempEntry->slicing = 0;

  // increment the Data objects counter
  totalDataObjects++;

  // return the pointer to the new entry
  return tempEntry;

}

std::string
DataProf::dumpProf(profDataEntry* entry) {

  stringstream temp_str;

  temp_str << "=============================\n";
  temp_str << "interpolate: " << entry->interpolate << "\n";
  temp_str << "grad       : " << entry->grad << "\n";
  temp_str << "integrate  : " << entry->integrate << "\n";
  temp_str << "where      : " << entry->where << "\n";
  temp_str << "unary      : " << entry->unary << "\n";
  temp_str << "binary     : " << entry->binary << "\n";
  temp_str << "reduction1 : " << entry->reduction1 << "\n";
  temp_str << "reduction2 : " << entry->reduction2 << "\n";
  temp_str << "slicing    : " << entry->slicing << "\n";
  temp_str << "=============================\n ";
  temp_str << endl;

  return temp_str.str();

}

string
DataProf::compProf() {

  int comp_interpolate = 0;
  int comp_grad = 0;
  int comp_integrate = 0;
  int comp_where = 0;
  int comp_unary = 0;
  int comp_binary = 0;
  int comp_reduction1 = 0;
  int comp_reduction2 = 0;
  int comp_slicing = 0;

  profDataTable* tempTable = profDataTable_Root;

  while (tempTable != 0) {
    comp_interpolate += tempTable->data->interpolate;
    comp_grad += tempTable->data->grad;
    comp_integrate += tempTable->data->integrate;
    comp_where += tempTable->data->where;
    comp_unary += tempTable->data->unary;
    comp_binary += tempTable->data->binary;
    comp_reduction1 += tempTable->data->reduction1;
    comp_reduction2 += tempTable->data->reduction2;
    comp_slicing += tempTable->data->slicing;
    tempTable = tempTable->next;
  }

  stringstream temp_str;

  temp_str << "========== Op Stats ===================\n";
  temp_str << "Total objects: " << totalDataObjects << "\n";
  temp_str << "interpolate  : " << comp_interpolate << "\n";
  temp_str << "grad         : " << comp_grad << "\n";
  temp_str << "integrate    : " << comp_integrate << "\n";
  temp_str << "where        : " << comp_where << "\n";
  temp_str << "unary        : " << comp_unary << "\n";
  temp_str << "binary       : " << comp_binary << "\n";
  temp_str << "reduction1   : " << comp_reduction1 << "\n";
  temp_str << "reduction2   : " << comp_reduction2 << "\n";
  temp_str << "slicing      : " << comp_slicing << "\n";
  temp_str << "======================================= ";
  temp_str << endl;

  return temp_str.str();

}

}  // end of namespace
