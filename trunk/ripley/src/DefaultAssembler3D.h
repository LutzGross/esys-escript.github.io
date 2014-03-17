
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
#ifndef __DEFAULT_ASSEMBLER3D_H__
#define __DEFAULT_ASSEMBLER3D_H__

#include <map>
#include <escript/Data.h>
#include <ripley/Ripley.h>
#include <ripley/RipleyException.h>
#include <ripley/AbstractAssembler.h>
#include <ripley/Brick.h>

namespace ripley {


class DefaultAssembler3D : public AbstractAssembler {
public:
    DefaultAssembler3D(Brick *dom, double *m_dx, dim_t *m_NX, dim_t *m_NE,
            dim_t *m_NN) : AbstractAssembler() {
        domain = dom;
        this->m_dx = m_dx;
        this->m_NX = m_NX;
        this->m_NE = m_NE;
        this->m_NN = m_NN;
    }
    ~DefaultAssembler3D(){};
    /* The default RipleyDomain assemblers, with original signatures */
    
    /// assembles a single PDE into the system matrix 'mat' and the right hand
    /// side 'rhs'
    virtual void assemblePDESingle(Paso_SystemMatrix* mat, escript::Data& rhs,
            const escript::Data& A, const escript::Data& B,
            const escript::Data& C, const escript::Data& D,
            const escript::Data& X, const escript::Data& Y) const;

    /// assembles boundary conditions of a single PDE into the system matrix
    /// 'mat' and the right hand side 'rhs'
    virtual void assemblePDEBoundarySingle(Paso_SystemMatrix* mat,
            escript::Data& rhs, const escript::Data& d,
            const escript::Data& y) const;

    /// assembles a single PDE with reduced order into the system matrix 'mat'
    /// and the right hand side 'rhs'
    virtual void assemblePDESingleReduced(Paso_SystemMatrix* mat,
            escript::Data& rhs, const escript::Data& A, const escript::Data& B,
            const escript::Data& C, const escript::Data& D,
            const escript::Data& X, const escript::Data& Y) const;

    /// assembles boundary conditions of a single PDE with reduced order into
    /// the system matrix 'mat' and the right hand side 'rhs'
    virtual void assemblePDEBoundarySingleReduced(Paso_SystemMatrix* mat,
            escript::Data& rhs, const escript::Data& d,
            const escript::Data& y) const;

    /// assembles a system of PDEs into the system matrix 'mat' and the right
    /// hand side 'rhs'
    virtual void assemblePDESystem(Paso_SystemMatrix* mat, escript::Data& rhs,
            const escript::Data& A, const escript::Data& B,
            const escript::Data& C, const escript::Data& D,
            const escript::Data& X, const escript::Data& Y) const;

    /// assembles boundary conditions of a system of PDEs into the system
    /// matrix 'mat' and the right hand side 'rhs'
    virtual void assemblePDEBoundarySystem(Paso_SystemMatrix* mat,
            escript::Data& rhs, const escript::Data& d,
            const escript::Data& y) const;

    /// assembles a system of PDEs with reduced order into the system matrix
    /// 'mat' and the right hand side 'rhs'
    virtual void assemblePDESystemReduced(Paso_SystemMatrix* mat,
            escript::Data& rhs, const escript::Data& A, const escript::Data& B,
            const escript::Data& C, const escript::Data& D,
            const escript::Data& X, const escript::Data& Y) const;

    /// assembles boundary conditions of a system of PDEs with reduced order
    /// into the system matrix 'mat' and the right hand side 'rhs'
    virtual void assemblePDEBoundarySystemReduced(Paso_SystemMatrix* mat,
            escript::Data& rhs, const escript::Data& d,
            const escript::Data& y) const;
    
    /* The new interface for assemblers */
    virtual void assemblePDESingle(Paso_SystemMatrix* mat, escript::Data& rhs,
        std::map<std::string, escript::Data> coefs) const;
    virtual void assemblePDEBoundarySingle(Paso_SystemMatrix* mat,
            escript::Data& rhs,
            std::map<std::string, escript::Data> coefs) const;
    virtual void assemblePDESingleReduced(Paso_SystemMatrix* mat,
            escript::Data& rhs,
            std::map<std::string, escript::Data> coefs) const;
    virtual void assemblePDEBoundarySingleReduced(Paso_SystemMatrix* mat,
            escript::Data& rhs,
            std::map<std::string, escript::Data> coefs) const;
    virtual void assemblePDESystem(Paso_SystemMatrix* mat, escript::Data& rhs,
            std::map<std::string, escript::Data> coefs) const;
    virtual void assemblePDEBoundarySystem(Paso_SystemMatrix* mat,
            escript::Data& rhs,
            std::map<std::string, escript::Data> coefs) const;
    virtual void assemblePDESystemReduced(Paso_SystemMatrix* mat,
            escript::Data& rhs,
            std::map<std::string, escript::Data> coefs) const;
    virtual void assemblePDEBoundarySystemReduced(Paso_SystemMatrix* mat,
            escript::Data& rhs,
            std::map<std::string, escript::Data> coefs) const;
protected:
    const Brick *domain;
    const double *m_dx;
    const dim_t *m_NX;
    const dim_t *m_NE;
    const dim_t *m_NN;
};

}


#endif
