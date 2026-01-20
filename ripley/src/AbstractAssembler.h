
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

#ifndef __RIPLEY_ABSTRACTASSEMBLER_H__
#define __RIPLEY_ABSTRACTASSEMBLER_H__

#include <escript/AbstractSystemMatrix.h>
#include <escript/Data.h>
#include <escript/Pointers.h>
#include <ripley/domainhelpers.h>

namespace ripley {

class AbstractAssembler;

typedef POINTER_WRAPPER_CLASS(AbstractAssembler) Assembler_ptr;
typedef POINTER_WRAPPER_CLASS(const AbstractAssembler) const_Assembler_ptr;

class AbstractAssembler : public REFCOUNT_BASE_CLASS(AbstractAssembler)
{
public:
    virtual ~AbstractAssembler() {}

    /* The new interface for assemblers */
    virtual void assemblePDESingle(escript::AbstractSystemMatrix* mat,
                    escript::Data& rhs, const DataMap& coefs) const = 0;
    virtual void assemblePDEBoundarySingle(escript::AbstractSystemMatrix* mat,
                    escript::Data& rhs, const DataMap& coefs) const = 0;
    virtual void assemblePDESingleReduced(escript::AbstractSystemMatrix* mat,
                    escript::Data& rhs, const DataMap& coefs) const = 0;
    virtual void assemblePDEBoundarySingleReduced(
                    escript::AbstractSystemMatrix* mat, escript::Data& rhs,
                    const DataMap& coefs) const = 0;
    virtual void assemblePDESystem(escript::AbstractSystemMatrix* mat,
                    escript::Data& rhs, const DataMap& coefs) const = 0;
    virtual void assemblePDEBoundarySystem(escript::AbstractSystemMatrix* mat,
                    escript::Data& rhs, const DataMap& coefs) const = 0;
    virtual void assemblePDESystemReduced(escript::AbstractSystemMatrix* mat,
                    escript::Data& rhs, const DataMap& coefs) const = 0;
    virtual void assemblePDEBoundarySystemReduced(
                    escript::AbstractSystemMatrix* mat, escript::Data& rhs,
                    const DataMap& coefs) const = 0;

    virtual void collateFunctionSpaceTypes(std::vector<int>& fsTypes,
                                           const DataMap& coefs) const = 0;

};

} // namespace ripley


#endif // __RIPLEY_ABSTRACTASSEMBLER_H__

