
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
* Development since 2012 by School of Earth Sciences
*
*****************************************************************************/
#ifndef __ABSTRACT_ASSEMBLER_H__
#define __ABSTRACT_ASSEMBLER_H__

#include <map>
#include <escript/Data.h>
#include <ripley/Ripley.h>
#include <ripley/RipleyException.h>
#include <paso/SystemMatrix.h>

namespace ripley {

escript::Data unpackData(std::string target,
        std::map<std::string, escript::Data> mapping);

class RipleyDomain;
/* returns the data associated with the string key or an empty data object
   if the map does not contain the given key */
escript::Data unpackData(std::string, std::map<std::string, escript::Data>);

class AbstractAssembler {
public:
    virtual ~AbstractAssembler() {};
    /* The default RipleyDomain assemblers, with original signatures */
    
    /// assembles a single PDE into the system matrix 'mat' and the right hand
    /// side 'rhs'
    virtual void assemblePDESingle(Paso_SystemMatrix* mat, escript::Data& rhs,
            const escript::Data& A, const escript::Data& B,
            const escript::Data& C, const escript::Data& D,
            const escript::Data& X, const escript::Data& Y) const {
        throw RipleyException("This assembler does not support "
                "old style signatures");
    }

    /// assembles boundary conditions of a single PDE into the system matrix
    /// 'mat' and the right hand side 'rhs'
    virtual void assemblePDEBoundarySingle(Paso_SystemMatrix* mat,
            escript::Data& rhs, const escript::Data& d,
            const escript::Data& y) const {
        throw RipleyException("This assembler does not support "
                "old style signatures");
    }

    /// assembles a single PDE with reduced order into the system matrix 'mat'
    /// and the right hand side 'rhs'
    virtual void assemblePDESingleReduced(Paso_SystemMatrix* mat,
            escript::Data& rhs, const escript::Data& A, const escript::Data& B,
            const escript::Data& C, const escript::Data& D,
            const escript::Data& X, const escript::Data& Y) const {
        throw RipleyException("This assembler does not support "
                "old style signatures");
    }
    
    /// assembles boundary conditions of a single PDE with reduced order into
    /// the system matrix 'mat' and the right hand side 'rhs'
    virtual void assemblePDEBoundarySingleReduced(Paso_SystemMatrix* mat,
            escript::Data& rhs, const escript::Data& d,
            const escript::Data& y) const {
        throw RipleyException("This assembler does not support "
                "old style signatures");
    }
    
    /// assembles a system of PDEs into the system matrix 'mat' and the right
    /// hand side 'rhs'
    virtual void assemblePDESystem(Paso_SystemMatrix* mat, escript::Data& rhs,
            const escript::Data& A, const escript::Data& B,
            const escript::Data& C, const escript::Data& D,
            const escript::Data& X, const escript::Data& Y) {
        throw RipleyException("This assembler does not support "
                "old style signatures");
    }
    
    /// assembles boundary conditions of a system of PDEs into the system
    /// matrix 'mat' and the right hand side 'rhs'
    virtual void assemblePDEBoundarySystem(Paso_SystemMatrix* mat,
            escript::Data& rhs, const escript::Data& d,
            const escript::Data& y) const {
        throw RipleyException("This assembler does not support "
                "old style signatures");
    }

    /// assembles a system of PDEs with reduced order into the system matrix
    /// 'mat' and the right hand side 'rhs'
    virtual void assemblePDESystemReduced(Paso_SystemMatrix* mat,
            escript::Data& rhs, const escript::Data& A, const escript::Data& B,
            const escript::Data& C, const escript::Data& D,
            const escript::Data& X, const escript::Data& Y) {
        throw RipleyException("This assembler does not support "
                "old style signatures");
    }

    /// assembles boundary conditions of a system of PDEs with reduced order
    /// into the system matrix 'mat' and the right hand side 'rhs'
    virtual void assemblePDEBoundarySystemReduced(Paso_SystemMatrix* mat,
            escript::Data& rhs, const escript::Data& d,
            const escript::Data& y) const {
        throw RipleyException("This assembler does not support "
                "old style signatures");
    }
    
    /* The new interface for assemblers */
    virtual void assemblePDESingle(Paso_SystemMatrix* mat, escript::Data& rhs,
        std::map<std::string, escript::Data> coefs) const = 0;
    virtual void assemblePDEBoundarySingle(Paso_SystemMatrix* mat,
            escript::Data& rhs,
            std::map<std::string, escript::Data> coefs) const = 0;
    virtual void assemblePDESingleReduced(Paso_SystemMatrix* mat,
            escript::Data& rhs,
            std::map<std::string, escript::Data> coefs) const = 0;
    virtual void assemblePDEBoundarySingleReduced(Paso_SystemMatrix* mat,
            escript::Data& rhs,
            std::map<std::string, escript::Data> coefs) const = 0;
    virtual void assemblePDESystem(Paso_SystemMatrix* mat, escript::Data& rhs,
            std::map<std::string, escript::Data> coefs) const = 0;
    virtual void assemblePDEBoundarySystem(Paso_SystemMatrix* mat,
            escript::Data& rhs,
            std::map<std::string, escript::Data> coefs) const = 0;
    virtual void assemblePDESystemReduced(Paso_SystemMatrix* mat,
            escript::Data& rhs,
            std::map<std::string, escript::Data> coefs) const = 0;
    virtual void assemblePDEBoundarySystemReduced(Paso_SystemMatrix* mat,
            escript::Data& rhs,
            std::map<std::string, escript::Data> coefs) const = 0;

};

}


#endif
