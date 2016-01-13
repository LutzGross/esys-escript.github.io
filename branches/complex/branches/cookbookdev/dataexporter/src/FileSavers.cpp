
/*******************************************************
*
* Copyright (c) 2003-2009 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/

#include <escriptexport/FileSavers.h>
#include <escriptexport/EscriptDataset.h>
#include <escript/Data.h>

#include <boost/python/extract.hpp>
#include <iostream>

using boost::python::dict;
using escript::Domain_ptr;
using std::string;

namespace escriptexport
{

// extracts Data objects and their names from src into vars and names
void unpackDict(const dict& src, DataVec& vars, StringVec& names)
{
    int numData = boost::python::extract<int>(src.attr("__len__")());
    if (numData > 0) {
        boost::python::list keys = src.keys();

        for (int i=0; i<numData; i++) {
            string varName = boost::python::extract<string>(keys[i]);
            names.push_back(varName);
            escript::Data varData = boost::python::extract<escript::Data>(
                    src[keys[i]]);
            vars.push_back(varData);
        }
    }

}

void saveSilo(const string& filename, int cycle, double time,
              Domain_ptr domain, const dict& datavars)
{
    DataVec vars;
    StringVec varNames;

    unpackDict(datavars, vars, varNames);

    EscriptDataset_ptr dataset = EscriptDataset_ptr(new EscriptDataset());
    if (!dataset->initFromEscript(domain, vars, varNames))
        throw escript::DataException("saveSilo: Error initialising dataset. No Silo support?!");

    dataset->setCycleAndTime(cycle, time);
    dataset->saveSilo(filename);
}

} // namespace escriptexport

