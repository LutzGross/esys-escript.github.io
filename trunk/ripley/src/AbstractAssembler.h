
/*****************************************************************************
*
* Copyright (c) 2003-2014 by University of Queensland
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
#ifndef __ESCRIPT_ABSTRACTASSEMBLER_H__
#define __ESCRIPT_ABSTRACTASSEMBLER_H__

#include <map>
#include <escript/Data.h>
#include <escript/Pointers.h>
#include <paso/SystemMatrix.h>

namespace ripley {
/* returns the data associated with the string key or an empty data object
   if the map does not contain the given key */
escript::Data unpackData(std::string target,
        std::map<std::string, escript::Data> mapping);

class AbstractAssembler;

typedef POINTER_WRAPPER_CLASS(AbstractAssembler) Assembler_ptr;
typedef POINTER_WRAPPER_CLASS(const AbstractAssembler) const_Assembler_ptr;

class AbstractAssembler : public REFCOUNT_BASE_CLASS(AbstractAssembler) {
public:
    virtual ~AbstractAssembler() {};
   
    /* The new interface for assemblers */
    virtual void assemblePDESingle(paso::SystemMatrix_ptr mat,
                    escript::Data& rhs,
                    std::map<std::string, escript::Data> coefs) const = 0;
    virtual void assemblePDEBoundarySingle(paso::SystemMatrix_ptr mat,
                    escript::Data& rhs,
                    std::map<std::string, escript::Data> coefs) const = 0;
    virtual void assemblePDESingleReduced(paso::SystemMatrix_ptr mat,
                    escript::Data& rhs,
                    std::map<std::string, escript::Data> coefs) const = 0;
    virtual void assemblePDEBoundarySingleReduced(paso::SystemMatrix_ptr mat,
                    escript::Data& rhs,
                    std::map<std::string, escript::Data> coefs) const = 0;
    virtual void assemblePDESystem(paso::SystemMatrix_ptr mat,
                    escript::Data& rhs,
                    std::map<std::string, escript::Data> coefs) const = 0;
    virtual void assemblePDEBoundarySystem(paso::SystemMatrix_ptr mat,
                    escript::Data& rhs,
                    std::map<std::string, escript::Data> coefs) const = 0;
    virtual void assemblePDESystemReduced(paso::SystemMatrix_ptr mat,
                    escript::Data& rhs,
                    std::map<std::string, escript::Data> coefs) const = 0;
    virtual void assemblePDEBoundarySystemReduced(paso::SystemMatrix_ptr mat,
                    escript::Data& rhs,
                    std::map<std::string, escript::Data> coefs) const = 0;

    virtual void collateFunctionSpaceTypes(std::vector<int>& fsTypes, 
            std::map<std::string, escript::Data> coefs) const = 0;
};

} // namespace escript


#endif // __ESCRIPT_ABSTRACTASSEMBLER_H__

