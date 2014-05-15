
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
#ifndef __RIPLEY_WAVEASSEMBLER2D_H__
#define __RIPLEY_WAVEASSEMBLER2D_H__

#include <map>
#include <escript/Data.h>
#include <ripley/Ripley.h>
#include <ripley/RipleyException.h>
#include <ripley/AbstractAssembler.h>
#include <ripley/Rectangle.h>

namespace ripley {


class WaveAssembler2D : public AbstractAssembler {
public:
    WaveAssembler2D(escript::const_Domain_ptr dom, const double *m_dx, const dim_t *m_NX, 
            const dim_t *m_NE, const dim_t *m_NN, std::map<std::string, escript::Data> c);
    ~WaveAssembler2D(){};

    /* The only assembly function we care about right now*/
    void assemblePDESystem(paso::SystemMatrix_ptr mat, escript::Data& rhs,
            std::map<std::string, escript::Data> coefs) const;


    void assemblePDESingle(paso::SystemMatrix_ptr mat, escript::Data& rhs,
        std::map<std::string, escript::Data> coefs) const {throw RipleyException("This assembly not supported by this assembler");}
    void assemblePDEBoundarySingle(paso::SystemMatrix_ptr mat, escript::Data& rhs,
            std::map<std::string, escript::Data> coefs) const {throw RipleyException("This assembly not supported by this assembler");}
    void assemblePDESingleReduced(paso::SystemMatrix_ptr mat, escript::Data& rhs,
            std::map<std::string, escript::Data> coefs) const {throw RipleyException("This assembly not supported by this assembler");}
    void assemblePDEBoundarySingleReduced(paso::SystemMatrix_ptr mat,
            escript::Data& rhs,
            std::map<std::string, escript::Data> coefs) const {throw RipleyException("This assembly not supported by this assembler");}
    void assemblePDEBoundarySystem(paso::SystemMatrix_ptr mat, escript::Data& rhs,
            std::map<std::string, escript::Data> coefs) const {throw RipleyException("This assembly not supported by this assembler");}
    void assemblePDESystemReduced(paso::SystemMatrix_ptr mat, escript::Data& rhs,
            std::map<std::string, escript::Data> coefs) const {throw RipleyException("This assembly not supported by this assembler");}
    void assemblePDEBoundarySystemReduced(paso::SystemMatrix_ptr mat,
            escript::Data& rhs,
            std::map<std::string, escript::Data> coefs) const {throw RipleyException("This assembly not supported by this assembler");}

    void collateFunctionSpaceTypes(std::vector<int>& fsTypes, 
            std::map<std::string, escript::Data> coefs) const;

private:
    std::map<std::string, escript::Data> c;
    boost::shared_ptr<const Rectangle> domain;
    const double *m_dx;
    const dim_t *m_NX;
    const dim_t *m_NE;
    const dim_t *m_NN;
    escript::Data c11, c12, c13, c23, c33, c44, c66;
    bool isVTI, isHTI;
};

} // namespace ripley

#endif // __RIPLEY_WAVEASSEMBLER2D_H__

