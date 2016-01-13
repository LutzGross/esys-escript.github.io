
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

#include <map>
#include <escript/Data.h>
#include <ripley/Ripley.h>
#include <ripley/RipleyException.h>
#include <ripley/AbstractAssembler.h>
#include <ripley/Rectangle.h>

namespace ripley {


class LameAssembler2D : public AbstractAssembler {
public:
    LameAssembler2D(escript::const_Domain_ptr dom, const double *m_dx,
            const dim_t *m_NX, const dim_t *m_NE,
            const dim_t *m_NN) : AbstractAssembler() {
        domain = boost::static_pointer_cast<const Rectangle>(dom);
        this->m_dx = m_dx;
        this->m_NX = m_NX;
        this->m_NE = m_NE;
        this->m_NN = m_NN;
    }
    ~LameAssembler2D(){};
    
    /* The new interface for assemblers */
    virtual void assemblePDESingle(SystemMatrix* mat,
            escript::Data& rhs,
            std::map<std::string, escript::Data> coefs) const;
    virtual void assemblePDEBoundarySingle(SystemMatrix* mat,
            escript::Data& rhs,
            std::map<std::string, escript::Data> coefs) const;
    virtual void assemblePDESingleReduced(SystemMatrix* mat,
            escript::Data& rhs,
            std::map<std::string, escript::Data> coefs) const;
    virtual void assemblePDEBoundarySingleReduced(SystemMatrix* mat,
            escript::Data& rhs,
            std::map<std::string, escript::Data> coefs) const;
    virtual void assemblePDESystem(SystemMatrix* mat, escript::Data& rhs,
            std::map<std::string, escript::Data> coefs) const;
    virtual void assemblePDEBoundarySystem(SystemMatrix* mat,
            escript::Data& rhs,
            std::map<std::string, escript::Data> coefs) const;
    virtual void assemblePDESystemReduced(SystemMatrix* mat,
            escript::Data& rhs,
            std::map<std::string, escript::Data> coefs) const;
    virtual void assemblePDEBoundarySystemReduced(SystemMatrix* mat,
            escript::Data& rhs,
            std::map<std::string, escript::Data> coefs) const;
            
    void collateFunctionSpaceTypes(std::vector<int>& fsTypes, 
            std::map<std::string, escript::Data> coefs) const;
protected:
    boost::shared_ptr<const Rectangle> domain;
    const double *m_dx;
    const dim_t *m_NX;
    const dim_t *m_NE;
    const dim_t *m_NN;
};

} // namespace ripley

#endif // __RIPLEY_LAMEASSEMBLER2D_H__

