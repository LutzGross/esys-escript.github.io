// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
// 
// See packages/seacas/LICENSE for details

#ifndef IOSS_Ionit_Initializer_h
#define IOSS_Ionit_Initializer_h

namespace Ioss {
  /** \brief A special namespace for a class used to initialize Ioss.
   */
  namespace Init {
    /** \brief Initialization of the Ioss library.
     *
     *  An object of this type must be created before using any other
     *  functions or methods in the Ioss library.
     */
    class Initializer
    {
    public:
      Initializer();
      ~Initializer();
      static Initializer &initialize_ioss();
      // Copy constructor
      // Assignment operator
    };
  } // namespace Init
} // namespace Ioss
#endif
