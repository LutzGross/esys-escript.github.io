
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
#ifndef __RIPLEY_LAMEASSEMBLER2D_H__
#define __RIPLEY_LAMEASSEMBLER2D_H__

#include <ripley/Rectangle.h>

namespace ripley {


class LameAssembler2D : public AbstractAssembler
{
public:
    LameAssembler2D(escript::const_Domain_ptr dom, const double *dx,
            const dim_t *NX, const dim_t *NE, const dim_t *NN)
        : AbstractAssembler(),
        m_dx(dx),
        m_NX(NX),
        m_NE(NE),
        m_NN(NN)
    {
        domain = boost::static_pointer_cast<const Rectangle>(dom);
    }
    ~LameAssembler2D() {}
    
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
    boost::shared_ptr<const Rectangle> domain;
    const double *m_dx;
    const dim_t *m_NX;
    const dim_t *m_NE;
    const dim_t *m_NN;
};

} // namespace ripley

#endif // __RIPLEY_LAMEASSEMBLER2D_H__

