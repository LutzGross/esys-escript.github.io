
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
#ifndef __WAVE_ASSEMBLER2D_H__
#define __WAVE_ASSEMBLER2D_H__

#include <map>
#include <escript/Data.h>
#include <ripley/Ripley.h>
#include <ripley/RipleyException.h>
#include <ripley/AbstractAssembler.h>
#include <ripley/Rectangle.h>

namespace ripley {


class WaveAssembler2D : public AbstractAssembler {
public:
    WaveAssembler2D(Rectangle *dom, double *m_dx, dim_t *m_NX, dim_t *m_NE,
                dim_t *m_NN, std::map<std::string, escript::Data> c);
    ~WaveAssembler2D(){};

    /* The only assembly function we care about right now*/
    void assemblePDESystem(Paso_SystemMatrix* mat, escript::Data& rhs,
            std::map<std::string, escript::Data> coefs) const;


    void assemblePDESingle(Paso_SystemMatrix* mat, escript::Data& rhs,
        std::map<std::string, escript::Data> coefs) const {throw RipleyException("This assembly not supported by this assembler");}
    void assemblePDEBoundarySingle(Paso_SystemMatrix* mat, escript::Data& rhs,
            std::map<std::string, escript::Data> coefs) const {throw RipleyException("This assembly not supported by this assembler");}
    void assemblePDESingleReduced(Paso_SystemMatrix* mat, escript::Data& rhs,
            std::map<std::string, escript::Data> coefs) const {throw RipleyException("This assembly not supported by this assembler");}
    void assemblePDEBoundarySingleReduced(Paso_SystemMatrix* mat,
            escript::Data& rhs,
            std::map<std::string, escript::Data> coefs) const {throw RipleyException("This assembly not supported by this assembler");}
    void assemblePDEBoundarySystem(Paso_SystemMatrix* mat, escript::Data& rhs,
            std::map<std::string, escript::Data> coefs) const {throw RipleyException("This assembly not supported by this assembler");}
    void assemblePDESystemReduced(Paso_SystemMatrix* mat, escript::Data& rhs,
            std::map<std::string, escript::Data> coefs) const {throw RipleyException("This assembly not supported by this assembler");}
    void assemblePDEBoundarySystemReduced(Paso_SystemMatrix* mat,
            escript::Data& rhs,
            std::map<std::string, escript::Data> coefs) const {throw RipleyException("This assembly not supported by this assembler");}

private:
    std::map<std::string, escript::Data> c;
    Rectangle *domain;
    double *m_dx;
    dim_t *m_NX;
    dim_t *m_NE;
    dim_t *m_NN;
    escript::Data c11, c12, c13, c33, c44, c66;
};

}


#endif
