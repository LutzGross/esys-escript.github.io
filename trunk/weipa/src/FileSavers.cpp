
/*******************************************************
*
* Copyright (c) 2003-2010 by University of Queensland
* Earth Systems Science Computational Center (ESSCC)
* http://www.uq.edu.au/esscc
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
*******************************************************/

#include <weipa/FileSavers.h>
#include <weipa/EscriptDataset.h>
#include <escript/Data.h>

#include <boost/python/extract.hpp>
#include <iostream>

using boost::python::dict;
using escript::Domain_ptr;
using std::string;

namespace weipa
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

    EscriptDataset_ptr dataset;
#ifdef PASO_MPI
    MPI_Comm comm = domain->getMPIComm();
    dataset.reset(new EscriptDataset(comm));
#else
    dataset.reset(new EscriptDataset());
#endif

    if (!dataset->initFromEscript(domain.get(), vars, varNames))
        throw escript::DataException("saveSilo: Error initialising dataset. No Silo support?!");

    dataset->setCycleAndTime(cycle, time);
    dataset->saveSilo(filename);
}

void saveVTK(const string& filename, int cycle, double time, Domain_ptr domain,
             const boost::python::dict& datavars, const string& metadata,
             const string& metadata_schema)
{
    DataVec vars;
    StringVec varNames;

    unpackDict(datavars, vars, varNames);

    EscriptDataset_ptr dataset;
#ifdef PASO_MPI
    MPI_Comm comm = domain->getMPIComm();
    dataset.reset(new EscriptDataset(comm));
#else
    dataset.reset(new EscriptDataset());
#endif

    if (!dataset->initFromEscript(domain.get(), vars, varNames))
        throw escript::DataException("saveVTK: Error initialising dataset.");

    dataset->setCycleAndTime(cycle, time);
    dataset->saveVTK(filename);
}

} // namespace weipa

