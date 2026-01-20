
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/
#ifndef __SPECKLEY_WAVE_ASSEMBLER_2D_H__
#define __SPECKLEY_WAVE_ASSEMBLER_2D_H__

#include <speckley/Rectangle.h>

namespace speckley {


class WaveAssembler2D : public AbstractAssembler
{
public:
    WaveAssembler2D(escript::const_Domain_ptr dom, const double *dx,
                       const dim_t *NE, const dim_t *NN, const DataMap& c)
        : AbstractAssembler(),
        m_dx(dx),
        m_NE(NE),
        m_NN(NN)
    {
        domain = REFCOUNTNS::static_pointer_cast<const Rectangle>(dom);
        isHTI = isVTI = false;
        DataMap::const_iterator a = c.find("c12"), b = c.find("c23");
        if (c.find("c11") == c.end()
                    || c.find("c13") == c.end() || c.find("c33") == c.end()
                    || c.find("c44") == c.end() || c.find("c66") == c.end()
                    || (a == c.end() && b == c.end()))
            throw SpeckleyException("required constants missing for WaveAssembler");

        if (a != c.end() && b != c.end()) {
            throw SpeckleyException("WaveAssembler3D() doesn't support general form waves");
        } else if (a == c.end()) {
            c23 = b->second;
            isHTI = true;
            if (c23.getFunctionSpace().getTypeCode() != ReducedElements) {
                throw SpeckleyException("C tensor elements must be reduced");
            }
            if (c23.isEmpty()) {
                throw SpeckleyException("C tensor elements must not be empty");
            }
        } else if (b == c.end()) {
            c12 = a->second;
            isVTI = true;
            if (c12.getFunctionSpace().getTypeCode() != ReducedElements) {
                throw SpeckleyException("C tensor elements must be reduced");
            }
            if (c12.isEmpty()) {
                throw SpeckleyException("C tensor elements must not be empty");
            }
        } // final else case taken care of with the missing constants above
        c11 = c.find("c11")->second;
        c13 = c.find("c13")->second;
        c33 = c.find("c33")->second;
        c44 = c.find("c44")->second;
        c66 = c.find("c66")->second;
        if (c11.getFunctionSpace().getTypeCode() != ReducedElements
                || c13.getFunctionSpace().getTypeCode() != ReducedElements
                || c33.getFunctionSpace().getTypeCode() != ReducedElements
                || c44.getFunctionSpace().getTypeCode() != ReducedElements 
                || c66.getFunctionSpace().getTypeCode() != ReducedElements) {
            throw SpeckleyException("C tensor elements must be reduced");
        }
        if (c11.isEmpty()
                || c13.isEmpty()
                || c33.isEmpty()
                || c44.isEmpty() 
                || c66.isEmpty()) {
            throw SpeckleyException("C tensor elements must not be empty");
        }
    }

    ~WaveAssembler2D() {}

    /* The default SpeckleyDomain assemblers, with original signatures */
    
    /// assembles a single PDE into the system matrix 'mat' and the right hand
    /// side 'rhs'
    virtual void assemblePDESingle(escript::AbstractSystemMatrix* mat, escript::Data& rhs,
            const escript::Data& A, const escript::Data& B,
            const escript::Data& C, const escript::Data& D,
            const escript::Data& du, const escript::Data& Y) const;

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
            const escript::Data& du, const escript::Data& Y) const;

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
            const escript::Data& du, const escript::Data& Y) const;

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
            const escript::Data& du, const escript::Data& Y) const;

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
    POINTER_WRAPPER_CLASS(const Rectangle) domain;
    const double *m_dx;
    const dim_t *m_NE;
    const dim_t *m_NN;
    bool isHTI, isVTI;
    escript::Data c11, c12, c13, c23, c33, c44, c66;
};

} // namespace speckley

#endif // __SPECKLEY_DEFAULTASSEMBLER2D_H__

