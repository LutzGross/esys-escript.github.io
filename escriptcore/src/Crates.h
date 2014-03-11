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

#ifndef ESCRIPT_CRATES_H
#define ESCRIPT_CRATES_H

#include <string>
#include <boost/shared_ptr.hpp>

namespace escript
{

class AbstractCrate
{
public:
    AbstractCrate(const std::string& name, const std::string& opname);
    // If we allow custom python functions as the operator, add another constructor chain
    
    virtual ~AbstractCrate();
    const std::string& getLabel();
protected:
    
    // true if "opname" is a legal operator for this crate
    virtual bool OpLegal(const std::string& opname)=0;	// pure virtual to force Crate to be abstract
    const std::string label;
    const std::string op;

};

typedef boost::shared_ptr<AbstractCrate> crate_ptr;

class DataCrate : public AbstractCrate
{
public:  
    DataCrate(const std::string& name, const std::string& opname);
    ~DataCrate();

protected:
    bool OpLegal(const std::string& opname);

};

class PickleCrate : public AbstractCrate
{
public:
    PickleCrate(const std::string& name, const std::string& opname);
    ~PickleCrate();
    
protected:
    bool OpLegal(const std::string& opname);

};




}

#endif