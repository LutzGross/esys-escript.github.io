
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
#ifndef __RIPLEY_LAMEASSEMBLER3D_H__
#define __RIPLEY_LAMEASSEMBLER3D_H__

#include <ripley/Brick.h>

namespace ripley {


class LameAssembler3D : public AbstractAssembler
{
public:
    LameAssembler3D(escript::const_Domain_ptr dom, const double *dx,
                    const dim_t *NE, const dim_t *NN)
        : AbstractAssembler(),
        m_dx(dx),
        m_NE(NE),
        m_NN(NN)
    {
        domain = REFCOUNTNS::static_pointer_cast<const Brick>(dom);
    }
    ~LameAssembler3D(){};
    
    /* The new interface for assemblers */

    virtual void assemblePDESingle(escript::AbstractSystemMatrix* mat,
                              escript::Data& rhs, const DataMap& coefs) const;
    virtual void assemblePDEBoundarySingle(escript::AbstractSystemMatrix* mat,
                              escript::Data& rhs, const DataMap& coefs) const;
    virtual void assemblePDESingleReduced(escript::AbstractSystemMatrix* mat,
                              escript::Data& rhs, const DataMap& coefs) const;
    virtual void assemblePDEBoundarySingleReduced(
                              escript::AbstractSystemMatrix* mat,
                              escript::Data& rhs, const DataMap& coefs) const;
    virtual void assemblePDESystem(escript::AbstractSystemMatrix* mat,
                              escript::Data& rhs, const DataMap& coefs) const;
    virtual void assemblePDEBoundarySystem(escript::AbstractSystemMatrix* mat,
                              escript::Data& rhs, const DataMap& coefs) const;
    virtual void assemblePDESystemReduced(escript::AbstractSystemMatrix* mat,
                              escript::Data& rhs, const DataMap& coefs) const;
    virtual void assemblePDEBoundarySystemReduced(
                              escript::AbstractSystemMatrix* mat,
                              escript::Data& rhs, const DataMap& coefs) const;
            
    void collateFunctionSpaceTypes(std::vector<int>& fsTypes,
                                   const DataMap& coefs) const;

protected:
    POINTER_WRAPPER_CLASS(const Brick) domain;
    const double *m_dx;
    const dim_t *m_NE;
    const dim_t *m_NN;
};

} // namespace ripley


#endif // __RIPLEY_LAMEASSEMBLER3D_H__

