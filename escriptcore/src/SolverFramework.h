
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

#ifndef __ESCRIPT_SOLVERFRAMEWORK_H__
#define __ESCRIPT_SOLVERFRAMEWORK_H__

#include "system_dep.h"
#include <memory>
#include <string>
#include <mutex>

namespace escript {

/**
   \brief
   Represents a solver backend (PASO or TRILINOS) and manages its runtime
   lifecycle.

   A SolverFramework object is created once and attached to a domain at
   construction time. The domain's framework determines which linear algebra
   library is used for all solves on that domain and is immutable after
   the domain has been set up.

   For Trilinos, the first SolverFramework(TRILINOS) instance initialises
   Kokkos/Tpetra; all subsequent instances share that initialisation.
   Finalisation is deferred to process exit (via the existing Trilinos
   atexit mechanism) to avoid ordering issues with long-lived static objects.

   Typical usage (Python side):
   \code
   from esys.escript import SolverFramework
   pkg = SolverFramework.trilinos()
   domain = Rectangle(20, 20, 1, framework=pkg)
   \endcode
*/
class ESCRIPT_DLL_API SolverFramework
{
public:
    /// Type constant for the Paso backend.
    static constexpr int PASO = 0;
    /// Type constant for the Trilinos backend.
    static constexpr int TRILINOS = 1;

    /// Constructs a SolverFramework for the given backend type.
    /// \param type  SolverFramework::PASO or SolverFramework::TRILINOS
    explicit SolverFramework(int type);

    ~SolverFramework();

    // Not copyable — one object, one lifecycle stake
    SolverFramework(const SolverFramework&) = delete;
    SolverFramework& operator=(const SolverFramework&) = delete;

    /// Returns the backend type identifier (SolverFramework::PASO or
    /// SolverFramework::TRILINOS).
    int getType() const { return m_type; }

    /// Returns true if this framework is the Trilinos backend.
    bool isTrilinos() const;

    /// Returns true if this framework is the Paso backend.
    bool isPaso() const;

    /// Returns a human-readable name for this framework ("PASO", "TRILINOS").
    std::string getName() const;

    /// Returns true if the given framework type was compiled into this build.
    /// \param type  SolverFramework::PASO or SolverFramework::TRILINOS
    static bool isAvailable(int type);

    /// Returns a SolverFramework for the Trilinos backend.
    /// Throws if Trilinos was not compiled into this build.
    static std::shared_ptr<SolverFramework> trilinos();

    /// Returns a SolverFramework for the Paso backend.
    /// Throws if Paso was not compiled into this build.
    static std::shared_ptr<SolverFramework> paso();

    /// Returns the process-wide default SolverFramework.
    /// The default is TRILINOS when available, PASO otherwise.
    /// The returned object is created once and lives for the process lifetime,
    /// so it is safe to call getDefault() from multiple domains without
    /// worrying about sequential init/finalize cycles.
    static std::shared_ptr<SolverFramework> getDefault();

private:
    int m_type;

    // Tracks the number of live Trilinos SolverFramework instances so we can
    // call Tpetra::initialize() exactly once.
    static int s_trilinasInstances;
    static std::mutex s_trilinosMutex;
};

} // namespace escript

#endif // __ESCRIPT_SOLVERFRAMEWORK_H__
