
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
#ifndef __SPECKLEY_DEFAULTASSEMBLER2D_H__
#define __SPECKLEY_DEFAULTASSEMBLER2D_H__

#include <speckley/Rectangle.h>

namespace speckley {


class DefaultAssembler2D : public AbstractAssembler
{
public:
    DefaultAssembler2D(escript::const_Domain_ptr dom, const double *dx,
                       const dim_t *NE, const dim_t *NN)
        : AbstractAssembler(),
        m_dx(dx),
        m_NE(NE),
        m_NN(NN)
    {
        domain = boost::static_pointer_cast<const Rectangle>(dom);
    }

    ~DefaultAssembler2D() {}

    /* The default SpeckleyDomain assemblers, with original signatures */
    
    /// assembles a single PDE into the system matrix 'mat' and the right hand
    /// side 'rhs'
    virtual void assemblePDESingle(escript::AbstractSystemMatrix* mat, escript::Data& rhs,
            const escript::Data& A, const escript::Data& B,
            const escript::Data& C, const escript::Data& D,
            const escript::Data& X, const escript::Data& Y) const;

    /// assembles boundary conditions of a single PDE into the system matrix
    /// 'mat' and the right hand side 'rhs'
    virtual void assemblePDEBoundarySingle(escript::AbstractSystemMatrix* mat,
            escript::Data& rhs, const escript::Data& d,
            const escript::Data& y) const;

    /// assembles a single PDE with reduced order into the system matrix 'mat'
    /// and the right hand side 'rhs'
    virtual void assemblePDESingleReduced(escript::AbstractSystemMatrix* mat,
            escript::Data& rhs, const escript::Data& A, const escript::Data& B,
            const escript::Data& C, const escript::Data& D,
            const escript::Data& X, const escript::Data& Y) const;

    /// assembles boundary conditions of a single PDE with reduced order into
    /// the system matrix 'mat' and the right hand side 'rhs'
    virtual void assemblePDEBoundarySingleReduced(escript::AbstractSystemMatrix* mat,
            escript::Data& rhs, const escript::Data& d,
            const escript::Data& y) const;

    /// assembles a system of PDEs into the system matrix 'mat' and the right
    /// hand side 'rhs'
    virtual void assemblePDESystem(escript::AbstractSystemMatrix* mat, escript::Data& rhs,
            const escript::Data& A, const escript::Data& B,
            const escript::Data& C, const escript::Data& D,
            const escript::Data& X, const escript::Data& Y) const;

    /// assembles boundary conditions of a system of PDEs into the system
    /// matrix 'mat' and the right hand side 'rhs'
    virtual void assemblePDEBoundarySystem(escript::AbstractSystemMatrix* mat,
            escript::Data& rhs, const escript::Data& d,
            const escript::Data& y) const;

    /// assembles a system of PDEs with reduced order into the system matrix
    /// 'mat' and the right hand side 'rhs'
    virtual void assemblePDESystemReduced(escript::AbstractSystemMatrix* mat,
            escript::Data& rhs, const escript::Data& A, const escript::Data& B,
            const escript::Data& C, const escript::Data& D,
            const escript::Data& X, const escript::Data& Y) const;

    /// assembles boundary conditions of a system of PDEs with reduced order
    /// into the system matrix 'mat' and the right hand side 'rhs'
    virtual void assemblePDEBoundarySystemReduced(escript::AbstractSystemMatrix* mat,
            escript::Data& rhs, const escript::Data& d,
            const escript::Data& y) const;
    
    /* The new interface for assemblers */

    virtual void assemblePDESingle(escript::AbstractSystemMatrix* mat,
                                   escript::Data& rhs,
                                   const DataMap& coefs) const;
    virtual void assemblePDEBoundarySingle(escript::AbstractSystemMatrix* mat,
                                           escript::Data& rhs,
                                           const DataMap& coefs) const;
    virtual void assemblePDESingleReduced(escript::AbstractSystemMatrix* mat,
                                          escript::Data& rhs,
                                          const DataMap& coefs) const;
    virtual void assemblePDEBoundarySingleReduced(
                                          escript::AbstractSystemMatrix* mat,
                                          escript::Data& rhs,
                                          const DataMap& coefs) const;
    virtual void assemblePDESystem(escript::AbstractSystemMatrix* mat,
                                   escript::Data& rhs,
                                   const DataMap& coefs) const;
    virtual void assemblePDEBoundarySystem(escript::AbstractSystemMatrix* mat,
                                           escript::Data& rhs,
                                           const DataMap& coefs) const;
    virtual void assemblePDESystemReduced(escript::AbstractSystemMatrix* mat,
                                          escript::Data& rhs,
                                          const DataMap& coefs) const;
    virtual void assemblePDEBoundarySystemReduced(
                                          escript::AbstractSystemMatrix* mat,
                                          escript::Data& rhs,
                                          const DataMap& coefs) const;
            
    virtual void collateFunctionSpaceTypes(std::vector<int>& fsTypes,
                                           const DataMap& coefs) const;

protected:
    boost::shared_ptr<const Rectangle> domain;
    const double *m_dx;
    const dim_t *m_NE;
    const dim_t *m_NN;
};

} // namespace speckley

#endif // __SPECKLEY_DEFAULTASSEMBLER2D_H__

