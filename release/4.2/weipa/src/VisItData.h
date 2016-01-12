
/*****************************************************************************
*
* Copyright (c) 2003-2016 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Open Software License version 3.0
* http://www.opensource.org/licenses/osl-3.0.php
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

#ifndef __WEIPA_VISITDATA_H__
#define __WEIPA_VISITDATA_H__

#include <escript/AbstractDomain.h>
#include <escript/Data.h>
#include <weipa/EscriptDataset.h>

#include <VisItInterfaceTypes_V2.h>

namespace weipa {

class VisItData {

public:
    VisItData() : runFlag(false) {}

    void publishData(EscriptDataset_ptr ds) { dataset=ds; }
    void setCommandNames(std::vector<std::string> names) { cmdNames=names; }
    void setSimulationStatus(bool running) { runFlag=running; }

    visit_handle getDomainList();
    visit_handle getMesh(const char* name);
    visit_handle getSimMetaData();
    visit_handle getVariable(const char* name);

private:
    void addExpressionMetadata(visit_handle smd, const std::string& name,
                               const std::string& def, int type);
    void addMeshMetadata(visit_handle smd, const std::string& name,
                         int dim, int numDoms);
    void addVariableMetadata(visit_handle smd, const std::string& name,
                             const std::string& meshName, int centering,
                             int rank);

    bool runFlag;
    EscriptDataset_ptr dataset;
    std::vector<std::string> cmdNames;
    std::map<std::string, DataVar_ptr> variables;
};

typedef boost::shared_ptr<VisItData> VisItData_ptr;


} // namespace weipa

#endif // __WEIPA_VISITDATA_H__

