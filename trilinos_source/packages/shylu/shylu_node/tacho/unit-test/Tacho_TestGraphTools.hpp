// clang-format off
/* =====================================================================================
Copyright 2022 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains
certain rights in this software.

SCR#:2790.0

This file is part of Tacho. Tacho is open source software: you can redistribute it
and/or modify it under the terms of BSD 2-Clause License
(https://opensource.org/licenses/BSD-2-Clause). A copy of the licese is also
provided under the main directory

Questions? Kyungjoo Kim at <kyukim@sandia.gov,https://github.com/kyungjoo-kim>

Sandia National Laboratories, Albuquerque, NM, USA
===================================================================================== */
// clang-format on
#ifndef __TACHO_TEST_GRAPH_TOOLS_HPP__
#define __TACHO_TEST_GRAPH_TOOLS_HPP__

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <Kokkos_Timer.hpp>

#include "Tacho_Util.hpp"
#include "Tacho_CrsMatrixBase.hpp"
#include "Tacho_MatrixMarket.hpp"

#include "Tacho_Graph.hpp"

#if defined(TACHO_HAVE_METIS)
#include "Tacho_GraphTools_Metis.hpp"
#endif

namespace Test {

  TEST( Graph, constructor ) {  
    
    ///
    /// host space crs matrix base
    ///
    const ordinal_type 
      m = 4,
      n = 4,
      nnz = 16;

    crs_matrix_base_type_host Ah(m, n, nnz);

    ordinal_type cnt = 0;
    for (ordinal_type i=0;i<m;++i) {
      Ah.RowPtrBegin(i) = cnt;
      for (ordinal_type j=0;j<n;++j,++cnt) {
        Ah.Col(cnt) = j;
        Ah.Value(cnt) = i*n+j;
      }
      Ah.RowPtrEnd(i) = cnt;
    }
  
    ///
    /// construction of graph
    ///
    Graph G(Ah);
  
    auto rptr = G.RowPtr();
    auto cidx = G.ColIdx();

    for (ordinal_type i=0;i<G.NumRows();++i) {
      const size_type jbeg = rptr(i), jend = rptr(i+1);
      for (size_type j=jbeg;j<jend;++j) {
        EXPECT_TRUE(cidx(j) != i);
      }
    }
  }

#if defined(TACHO_HAVE_METIS)
  TEST( Graph, metis ) {
    
    std::string filename = MM_TEST_FILE;
    crs_matrix_base_type_host Ah;
    MatrixMarket<value_type>::read(filename, Ah);

    Graph G(Ah);
    GraphTools_Metis M(G);

    M.reorder();

    ///
    /// perm and invperm should be properly setup 
    ///
    const auto perm = M.PermVector();
    const auto peri = M.InvPermVector();

    const ordinal_type m = G.NumRows();
    for (ordinal_type i=0;i<m;++i) {
      EXPECT_EQ(i, peri(perm(i)));
    }
  }
#endif
}

#endif
