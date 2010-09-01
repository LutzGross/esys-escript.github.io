
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

#ifndef __VISITDATA_H__
#define __VISITDATA_H__

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

#endif // __VISITDATA_H__

