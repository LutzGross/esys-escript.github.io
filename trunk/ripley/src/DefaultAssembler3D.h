
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
#ifndef __RIPLEY_DEFAULTASSEMBLER3D_H__
#define __RIPLEY_DEFAULTASSEMBLER3D_H__

#include <map>
#include <escript/Data.h>
#include <ripley/Ripley.h>
#include <ripley/RipleyException.h>
#include <escript/AbstractAssembler.h>
#include <ripley/Brick.h>

namespace ripley {


class DefaultAssembler3D : public escript::AbstractAssembler {
public:
    DefaultAssembler3D(escript::const_Domain_ptr dom, const double *m_dx, const dim_t *m_NX, 
            const dim_t *m_NE, const dim_t *m_NN) : escript::AbstractAssembler() {
        domain = boost::static_pointer_cast<const Brick>(dom);
        this->m_dx = m_dx;
        this->m_NX = m_NX;
        this->m_NE = m_NE;
        this->m_NN = m_NN;
    }
    ~DefaultAssembler3D(){};
    /* The default RipleyDomain assemblers, with original signatures */
    
    /// assembles a single PDE into the system matrix 'mat' and the right hand
    /// side 'rhs'
    virtual void assemblePDESingle(paso::SystemMatrix_ptr mat,
            escript::Data& rhs,
            const escript::Data& A, const escript::Data& B,
            const escript::Data& C, const escript::Data& D,
            const escript::Data& X, const escript::Data& Y) const;

    /// assembles boundary conditions of a single PDE into the system matrix
    /// 'mat' and the right hand side 'rhs'
    virtual void assemblePDEBoundarySingle(paso::SystemMatrix_ptr mat,
            escript::Data& rhs, const escript::Data& d,
            const escript::Data& y) const;

    /// assembles a single PDE with reduced order into the system matrix 'mat'
    /// and the right hand side 'rhs'
    virtual void assemblePDESingleReduced(paso::SystemMatrix_ptr mat,
            escript::Data& rhs, const escript::Data& A, const escript::Data& B,
            const escript::Data& C, const escript::Data& D,
            const escript::Data& X, const escript::Data& Y) const;

    /// assembles boundary conditions of a single PDE with reduced order into
    /// the system matrix 'mat' and the right hand side 'rhs'
    virtual void assemblePDEBoundarySingleReduced(paso::SystemMatrix_ptr mat,
            escript::Data& rhs, const escript::Data& d,
            const escript::Data& y) const;

    /// assembles a system of PDEs into the system matrix 'mat' and the right
    /// hand side 'rhs'
    virtual void assemblePDESystem(paso::SystemMatrix_ptr mat, escript::Data& rhs,
            const escript::Data& A, const escript::Data& B,
            const escript::Data& C, const escript::Data& D,
            const escript::Data& X, const escript::Data& Y) const;

    /// assembles boundary conditions of a system of PDEs into the system
    /// matrix 'mat' and the right hand side 'rhs'
    virtual void assemblePDEBoundarySystem(paso::SystemMatrix_ptr mat,
            escript::Data& rhs, const escript::Data& d,
            const escript::Data& y) const;

    /// assembles a system of PDEs with reduced order into the system matrix
    /// 'mat' and the right hand side 'rhs'
    virtual void assemblePDESystemReduced(paso::SystemMatrix_ptr mat,
            escript::Data& rhs, const escript::Data& A, const escript::Data& B,
            const escript::Data& C, const escript::Data& D,
            const escript::Data& X, const escript::Data& Y) const;

    /// assembles boundary conditions of a system of PDEs with reduced order
    /// into the system matrix 'mat' and the right hand side 'rhs'
    virtual void assemblePDEBoundarySystemReduced(paso::SystemMatrix_ptr mat,
            escript::Data& rhs, const escript::Data& d,
            const escript::Data& y) const;
    
    /* The new interface for assemblers */
    virtual void assemblePDESingle(paso::SystemMatrix_ptr mat, escript::Data& rhs,
        std::map<std::string, escript::Data> coefs) const;
    virtual void assemblePDEBoundarySingle(paso::SystemMatrix_ptr mat,
            escript::Data& rhs,
            std::map<std::string, escript::Data> coefs) const;
    virtual void assemblePDESingleReduced(paso::SystemMatrix_ptr mat,
            escript::Data& rhs,
            std::map<std::string, escript::Data> coefs) const;
    virtual void assemblePDEBoundarySingleReduced(paso::SystemMatrix_ptr mat,
            escript::Data& rhs,
            std::map<std::string, escript::Data> coefs) const;
    virtual void assemblePDESystem(paso::SystemMatrix_ptr mat, escript::Data& rhs,
            std::map<std::string, escript::Data> coefs) const;
    virtual void assemblePDEBoundarySystem(paso::SystemMatrix_ptr mat,
            escript::Data& rhs,
            std::map<std::string, escript::Data> coefs) const;
    virtual void assemblePDESystemReduced(paso::SystemMatrix_ptr mat,
            escript::Data& rhs,
            std::map<std::string, escript::Data> coefs) const;
    virtual void assemblePDEBoundarySystemReduced(paso::SystemMatrix_ptr mat,
            escript::Data& rhs,
            std::map<std::string, escript::Data> coefs) const;
            
    void collateFunctionSpaceTypes(std::vector<int>& fsTypes, 
            std::map<std::string, escript::Data> coefs) const;
protected:
    boost::shared_ptr<const Brick> domain;
    const double *m_dx;
    const dim_t *m_NX;
    const dim_t *m_NE;
    const dim_t *m_NN;
};

} // namespace ripley

#endif // __RIPLEY_DEFAULTASSEMBLER3D_H__

