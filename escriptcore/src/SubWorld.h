
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

namespace escript
{
  
/** class to hold a collection of MPI ranks and a communicator linking them
*/
class SubWorld : public boost::enable_shared_from_this<SubWorld>
{
public:
    SubWorld(MPI_Comm comm);
    ~SubWorld();
    void setDomain(Domain_ptr d);
    Domain_ptr getDomain();
    MPI_Comm getComm();
    void addJob(boost::python::object j);
    void runJobs();
private:    
    MPI_Comm communicator;
    escript::Domain_ptr domain;
    std::vector<boost::python::object> jobvec;
};

typedef boost::shared_ptr<SubWorld> SubWorld_ptr;



}
#endif
