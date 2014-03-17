
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
#ifndef __LAME_ASSEMBLER3D_H__
#define __LAME_ASSEMBLER3D_H__

#include <map>
#include <escript/Data.h>
#include <ripley/Ripley.h>
#include <ripley/RipleyException.h>
#include <ripley/AbstractAssembler.h>
#include <ripley/Brick.h>

namespace ripley {


class LameAssembler3D : public AbstractAssembler {
public:
    LameAssembler3D(Brick *dom, double *m_dx, dim_t *m_NX, dim_t *m_NE,
            dim_t *m_NN) : AbstractAssembler() {
        domain = dom;
        this->m_dx = m_dx;
        this->m_NX = m_NX;
        this->m_NE = m_NE;
        this->m_NN = m_NN;
    }
    ~LameAssembler3D(){};
    
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
