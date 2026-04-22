
/*****************************************************************************
*
* Copyright (c) 2003-2026 by the esys.escript Group
* https://github.com/LutzGross/esys-escript.github.io
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* See CREDITS file for contributors and development history
**
*****************************************************************************/

#include "SolverFramework.h"
#include "EsysException.h"

#ifdef ESYS_HAVE_TRILINOS
#include <Tpetra_Core.hpp>
#endif

#include <stdexcept>

namespace escript {

// ---- static members --------------------------------------------------------

int SolverFramework::s_trilinasInstances = 0;
std::mutex SolverFramework::s_trilinosMutex;

// ---- construction / destruction --------------------------------------------

SolverFramework::SolverFramework(int type)
    : m_type(type)
{
#ifdef ESYS_HAVE_TRILINOS
    if (type == TRILINOS) {
        std::lock_guard<std::mutex> lock(s_trilinosMutex);
        if (s_trilinasInstances == 0 && !Tpetra::isInitialized()) {
            int argc = 0;
            char** argv = nullptr;
            Tpetra::initialize(&argc, &argv);
        }
        ++s_trilinasInstances;
        return;
    }
#endif
    if (type == TRILINOS) {
        throw EsysException("SolverFramework: Trilinos requested but this "
                               "build was not compiled with Trilinos support.");
    }
#ifndef ESYS_HAVE_PASO
    if (type == PASO) {
        throw EsysException("SolverFramework: Paso requested but this "
                               "build was not compiled with Paso support.");
    }
#endif
}

SolverFramework::~SolverFramework()
{
#ifdef ESYS_HAVE_TRILINOS
    if (m_type == TRILINOS) {
        std::lock_guard<std::mutex> lock(s_trilinosMutex);
        --s_trilinasInstances;
        // Trilinos/Kokkos finalization is intentionally deferred to process
        // exit via Tpetra's own atexit mechanism.  Calling Tpetra::finalize()
        // here risks a double-free with long-lived static Belos objects
        // (see esys.trilinos/include/BelosMultiVecTraits_Tpetra.hpp patch and
        // Trilinos issue #11976).  This will be revisited once that upstream
        // bug is resolved.
    }
#endif
}

// ---- public methods --------------------------------------------------------

bool SolverFramework::isTrilinos() const
{
    return m_type == TRILINOS;
}

bool SolverFramework::isPaso() const
{
    return m_type == PASO;
}

std::string SolverFramework::getName() const
{
    switch (m_type) {
        case TRILINOS: return "TRILINOS";
        case PASO:     return "PASO";
        default:                  return "UNKNOWN";
    }
}

bool SolverFramework::isAvailable(int type)
{
    switch (type) {
        case TRILINOS:
#ifdef ESYS_HAVE_TRILINOS
            return true;
#else
            return false;
#endif
        case PASO:
#ifdef ESYS_HAVE_PASO
            return true;
#else
            return false;
#endif
        default:
            return false;
    }
}

std::shared_ptr<SolverFramework> SolverFramework::trilinos()
{
#ifdef ESYS_HAVE_TRILINOS
    return std::make_shared<SolverFramework>(TRILINOS);
#else
    throw EsysException("SolverFramework::trilinos: this build was not compiled "
                        "with Trilinos support.");
#endif
}

std::shared_ptr<SolverFramework> SolverFramework::paso()
{
#ifdef ESYS_HAVE_PASO
    return std::make_shared<SolverFramework>(PASO);
#else
    throw EsysException("SolverFramework::paso: this build was not compiled "
                        "with Paso support.");
#endif
}

std::shared_ptr<SolverFramework> SolverFramework::getDefault()
{
    // Created once; lives for the process lifetime so that Kokkos is never
    // finalized mid-run even if all user-visible domain objects are deleted.
    static std::shared_ptr<SolverFramework> defaultPkg;
    static std::once_flag flag;
    std::call_once(flag, []() {
#ifdef ESYS_HAVE_TRILINOS
        defaultPkg = std::make_shared<SolverFramework>(TRILINOS);
#elif defined(ESYS_HAVE_PASO)
        defaultPkg = std::make_shared<SolverFramework>(PASO);
#endif
    });
    return defaultPkg;
}

} // namespace escript
