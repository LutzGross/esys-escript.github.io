
/*****************************************************************************
*
* Copyright (c) 2014 by University of Queensland
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

#ifndef escript_SubWorld_H
#define escript_SubWorld_H

#include "esysUtils/Esys_MPI.h"
#include "AbstractDomain.h"
#include "Reducer.h"

namespace escript
{
  
/** class to hold a collection of MPI ranks and a communicator linking them
*/
class SubWorld : public boost::enable_shared_from_this<SubWorld>
{
public:
    SubWorld(esysUtils::JMPI& comm);
    ~SubWorld();
    void setDomain(Domain_ptr d);
    Domain_ptr getDomain();
    esysUtils::JMPI& getMPI();
    void addJob(boost::python::object j);
    char runJobs(std::string& errmsg);
    void clearJobs();
    void clearImportExports();
    void addVariable(std::string&, Reducer_ptr& red, bool manualimport);
    void removeVariable(std::string& name);  
    
    
    bool localTransport(std::vector<char>& vb, std::string& errmsg);
    
private:    
    esysUtils::JMPI mpiinfo;
    escript::Domain_ptr domain;
    std::vector<boost::python::object> jobvec;
    
typedef std::map<std::string, Reducer_ptr> str2reduce;  
typedef std::map<std::string, bool> str2bool;
    str2reduce reducemap;
    str2bool importmap;

};

typedef boost::shared_ptr<SubWorld> SubWorld_ptr;



}
#endif
