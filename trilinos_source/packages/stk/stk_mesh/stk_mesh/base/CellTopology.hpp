// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#ifndef stk_mesh_fem_CellTopology_hpp
#define stk_mesh_fem_CellTopology_hpp

#ifdef HAVE_SHARDS_DEBUG
#define STK_MESH_FEM_CHECK_REQUIRE( S )  S
#else
#define STK_MESH_FEM_CHECK_REQUIRE( S ) /* empty */
#endif

#include <Shards_CellTopologyTraits.hpp>
#include <Shards_CellTopology.hpp>
#include <Shards_CellTopologyData.h>
#include <Shards_BasicTopologies.hpp>

namespace stk {
namespace mesh {

typedef shards::CellTopology CellTopology;

#ifndef STK_HIDE_DEPRECATED_CODE // Delete after 2019-07-18
template< typename id_type >
STK_DEPRECATED
int findPermutation( const CellTopology top ,
                     const id_type * const expected_node ,
                     const id_type * const actual_node )
{
  return shards::findPermutation( *top.getCellTopologyData() , expected_node , actual_node );
}
#endif

/** \} */

} // namespace mesh
} // namespace stk

#undef STK_MESH_FEM_CHECK_REQUIRE

#endif // stk_mesh_fem_CellTopology_hpp
