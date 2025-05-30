// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "Ioss_CodeTypes.h"           // for IntVector
#include "Ioss_ElementTopology.h"     // for ElementTopology
#include <Ioss_ElementVariableType.h> // for ElementVariableType
#include <Ioss_Tet40.h>
#include <cassert> // for assert

//------------------------------------------------------------------------
// Define a variable type for storage of this elements connectivity
namespace Ioss {
  const char *Tet40::name = "tetra40";
  class St_Tet40 : public ElementVariableType
  {
  public:
    static void factory() { static St_Tet40 registerThis; }

  protected:
    St_Tet40() : ElementVariableType(Ioss::Tet40::name, 40) {}
  };
} // namespace Ioss

// ========================================================================
namespace {
  struct Constants
  {
    static const int nnode     = 40;
    static const int nedge     = 6;
    static const int nedgenode = 4;
    static const int nface     = 4;
    static const int nfacenode = 13;
    static const int nfaceedge = 3;
    static int       edge_node_order[nedge][nedgenode];
    static int       face_node_order[nface][nfacenode];
    static int       face_edge_order[nface][nfaceedge];
    static int       nodes_per_face[nface + 1];
    static int       edges_per_face[nface + 1];
  };

  // Edge numbers are zero-based [0..number_edges)
  int Constants::edge_node_order[nedge][nedgenode] = // [edge][edge_node]
      {{0, 1, 4, 5}, {1, 2, 6, 7}, {2, 0, 8, 9}, {0, 3, 10, 13}, {1, 3, 11, 14}, {2, 3, 12, 15}};

  // Face numbers are zero-based [0..number_faces)
  int Constants::face_node_order[nface][nfacenode] = // [face][face_node]
      {{0, 1, 3, 4, 5, 11, 14, 13, 10, 20, 21, 31, 30},
       {1, 2, 3, 6, 7, 12, 15, 14, 11, 22, 23, 33, 32},
       {0, 3, 2, 10, 13, 15, 12, 8, 9, 24, 34, 35, 25},
       {0, 2, 1, 9, 8, 7, 6, 5, 4, 19, 18, 17, 16}};

  int Constants::face_edge_order[nface][nfaceedge] = // [face][face_edge]
      {{0, 4, 3}, {1, 5, 4}, {3, 5, 2}, {2, 1, 0}};

  // face 0 returns number of nodes for all faces if homogeneous
  //        returns -1 if faces have differing topology
  int Constants::nodes_per_face[nface + 1] = {nfacenode, nfacenode, nfacenode, nfacenode,
                                              nfacenode};

  // face 0 returns number of edges for all faces if homogeneous
  //        returns -1 if faces have differing topology
  int Constants::edges_per_face[nface + 1] = {nfaceedge, nfaceedge, nfaceedge, nfaceedge,
                                              nfaceedge};
} // namespace

void Ioss::Tet40::factory()
{
  static Ioss::Tet40 registerThis;
  Ioss::St_Tet40::factory();
}

Ioss::Tet40::Tet40() : Ioss::ElementTopology(Ioss::Tet40::name, "Tetrahedron_40")
{
  Ioss::ElementTopology::alias(Ioss::Tet40::name, "tet40");
  Ioss::ElementTopology::alias(Ioss::Tet40::name, "Solid_Tet_40_3D");
}

int Ioss::Tet40::parametric_dimension() const { return 3; }
int Ioss::Tet40::spatial_dimension() const { return 3; }
int Ioss::Tet40::order() const { return 3; }

int Ioss::Tet40::number_corner_nodes() const { return 4; }
int Ioss::Tet40::number_nodes() const { return Constants::nnode; }
int Ioss::Tet40::number_edges() const { return Constants::nedge; }
int Ioss::Tet40::number_faces() const { return Constants::nface; }

int Ioss::Tet40::number_nodes_edge(int /* edge */) const { return Constants::nedgenode; }

int Ioss::Tet40::number_nodes_face(int face) const
{
  // face is 1-based.  0 passed in for all faces.
  assert(face >= 0 && face <= number_faces());
  return Constants::nodes_per_face[face];
}

int Ioss::Tet40::number_edges_face(int face) const
{
  // face is 1-based.  0 passed in for all faces.
  assert(face >= 0 && face <= number_faces());
  return Constants::edges_per_face[face];
}

Ioss::IntVector Ioss::Tet40::edge_connectivity(int edge_number) const
{
  assert(edge_number > 0 && edge_number <= Constants::nedge);
  Ioss::IntVector connectivity(Constants::nedgenode);

  for (int i = 0; i < Constants::nedgenode; i++) {
    connectivity[i] = Constants::edge_node_order[edge_number - 1][i];
  }

  return connectivity;
}

Ioss::IntVector Ioss::Tet40::face_connectivity(int face_number) const
{
  assert(face_number > 0 && face_number <= number_faces());
  Ioss::IntVector connectivity(Constants::nodes_per_face[face_number]);

  for (int i = 0; i < Constants::nodes_per_face[face_number]; i++) {
    connectivity[i] = Constants::face_node_order[face_number - 1][i];
  }

  return connectivity;
}

Ioss::IntVector Ioss::Tet40::element_connectivity() const
{
  Ioss::IntVector connectivity(number_nodes());
  for (int i = 0; i < number_nodes(); i++) {
    connectivity[i] = i;
  }
  return connectivity;
}

Ioss::ElementTopology *Ioss::Tet40::face_type(int face_number) const
{
  // face_number == 0 returns topology for all faces if
  // all faces are the same topology; otherwise, returns nullptr
  // face_number is 1-based.

  assert(face_number >= 0 && face_number <= number_faces());
  IOSS_ASSERT_USED(face_number);
  return Ioss::ElementTopology::factory("tri13");
}

Ioss::ElementTopology *Ioss::Tet40::edge_type(int edge_number) const
{
  // edge_number == 0 returns topology for all edges if
  // all edges are the same topology; otherwise, returns nullptr
  // edge_number is 1-based.

  assert(edge_number >= 0 && edge_number <= number_edges());
  IOSS_ASSERT_USED(edge_number);
  return Ioss::ElementTopology::factory("edge4");
}

Ioss::IntVector Ioss::Tet40::face_edge_connectivity(int face_number) const
{
  assert(face_number > 0 && face_number <= Constants::nface);

  int             nface_edge = number_edges_face(face_number);
  Ioss::IntVector fcon(nface_edge);

  for (int i = 0; i < nface_edge; i++) {
    fcon[i] = Constants::face_edge_order[face_number - 1][i];
  }

  return fcon;
}
