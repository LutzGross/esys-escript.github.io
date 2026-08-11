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
*
*****************************************************************************/

#ifndef __OXLEY_FINLEYCONVERTER_H__
#define __OXLEY_FINLEYCONVERTER_H__

#include <oxley/OxleyDomain.h>

#include <escript/AbstractDomain.h>

namespace oxley {

/**
   \brief
   Builds a finley domain describing the same mesh as an oxley forest.

   The forest supplies the geometry and topology; finley owns the degrees of
   freedom, the parallel overlap and everything a PDE needs. Each rank hands
   over the part of the forest it holds, with nodes named by their lnodes global
   id, and finley resolves the references and distributes the mesh, so no
   communication happens here beyond the conformity check.

   Every octant is split into simplices - Tri3 in 2D, Tet4 in 3D - so that one
   element family serves both conforming and adaptive forests. A quad face is
   always cut by the diagonal through its lowest global node id, which both
   sides of the face compute independently and agree on, so the split needs no
   communication.

   Currently restricted to CONFORMING forests: the split of a hanging face
   needs global ids for the hanging positions, which lnodes does not give them.
   Converting a forest with hanging nodes throws.

   \param dom the forest to convert
   \param order integration order (1 or 2)
   \param reducedOrder reduced integration order (1 or 2)
   \param optimize whether to let finley repartition with ParMETIS. Pass false
                   to keep the partition the forest already has.
   \param simplices when false, emit one Rec4/Hex8 per octant instead of
                    splitting. Only for telling a fault in the finley handover
                    apart from a fault in the split - it is not a production
                    path, since a forest would then change element family the
                    moment refinement introduces a hanging node.
*/
/**
   \brief
   Copies a ContinuousFunction from an oxley forest onto its finley export.

   The two meshes name their shared nodes identically - the export hands finley
   the ids in MeshAccess::nodeFinleyId - so this is a lookup by global id, not
   an interpolation: a value is carried across unchanged.

   The export has nodes the forest does not: the positions materialised at a
   2:1 seam, which are no nodes of the oxley domain and so carry no value of
   their own. They take the average of their masters, which is what makes the
   result continuous - the same constraint the forest leaves implicit.

   Not a method on the domain and not routed through escript's interpolate():
   Data::interpolate dispatches on the SOURCE domain, so finley -> oxley would
   have to be implemented inside finley, and finley must not know about oxley
   (the build enforces that direction).

   \param source a Data on ContinuousFunction (or Solution) of an oxley domain
   \param target the finley domain built from that forest by toFinley()
   \return the same field on ContinuousFunction of the finley domain
*/
escript::Data toFinleyData(const escript::Data& source,
                           escript::Domain_ptr target);

/**
   \brief
   Copies a ContinuousFunction from a finley export back onto its oxley forest.

   The exact inverse of toFinleyData for the nodes the two meshes share. The
   values at the materialised seam positions are dropped: they are not nodes of
   the forest, so there is nowhere to put them.

   \param source a Data on ContinuousFunction (or Solution) of the export
   \param target the oxley domain the export was built from
   \return the same field on ContinuousFunction of the oxley domain
*/
escript::Data fromFinleyData(const escript::Data& source,
                             escript::Domain_ptr target);

/**
   \brief
   Copies a ReducedFunction - one value per element - onto the finley export.

   An octant becomes 2 to 6 triangles, so the value is REPLICATED onto each of
   them. The simplices carry ids that say which octant they came from, so no
   map has to be stored and none has to survive finley redistributing the mesh.

   \param source a Data on ReducedFunction of an oxley domain
   \param target the finley domain toFinley() built from that forest
*/
escript::Data toFinleyReducedData(const escript::Data& source,
                                  escript::Domain_ptr target);

/**
   \brief
   Brings a ReducedFunction back from the finley export onto the forest.

   The reverse of the replication: an octant's value is the AREA-WEIGHTED mean
   of its simplices, so a field constant over the octant returns unchanged no
   matter how the octant was split. The weights are recomputed on the oxley
   side from the split, which is a deterministic function of the hanging
   configuration - the finley side never has to send areas.

   This is the transfer the adaptive loop needs: the error indicator is
   computed on the finley mesh and comes home per octant.

   \param source a Data on ReducedFunction of the export
   \param target the oxley domain the export was built from
*/
escript::Data fromFinleyReducedData(const escript::Data& source,
                                    escript::Domain_ptr target);

/**
   \brief
   Copies a FunctionOnBoundary (or its reduced form) onto the finley export.

   The 2D split never subdivides a boundary edge - a hanging node is the
   midpoint of a face that HAS a finer neighbour, so it is interior - which
   makes the boundary faces correspond one to one. Measured: the two meshes put
   the same physical quadrature points on the boundary, in both spaces.

   They do not agree on the ORDER, though, so this is a permutation rather than
   a copy. Each side sorts a face's points by coordinate, which gives both the
   same order without either having to send coordinates.

   \param source a Data on FunctionOnBoundary of an oxley domain
   \param target the finley domain toFinley() built from that forest
*/
escript::Data toFinleyBoundaryData(const escript::Data& source,
                                   escript::Domain_ptr target);

/**
   \brief
   Brings a FunctionOnBoundary (or its reduced form) back onto the forest.

   The inverse of toFinleyBoundaryData; see there for why it is a permutation.

   \param source a Data on FunctionOnBoundary of the export
   \param target the oxley domain the export was built from
*/
escript::Data fromFinleyBoundaryData(const escript::Data& source,
                                     escript::Domain_ptr target);

/**
   \brief
   Copies a Function - values at the quadrature points - onto the export.

   The one transfer that is not a rearrangement: the two meshes use different
   rules. An octant carries the 2x2 Gauss points; a Tri3 carries its three edge
   midpoints, neither a subset of the other.

   So this evaluates rather than moves. The four values on an octant are
   unisolvent for a bilinear function, and that function is evaluated at the
   triangles' points - exact for any field that is bilinear on the octant,
   which includes every linear one, and an approximation otherwise.

   \param source a Data on Function of an oxley domain
   \param target the finley domain toFinley() built from that forest
*/
escript::Data toFinleyFunctionData(const escript::Data& source,
                                   escript::Domain_ptr target);

/**
   \brief
   Brings a Function back from the export onto the forest.

   The three midpoints of a triangle are unisolvent for a linear function, so
   each triangle's values rebuild one, and the octant's Gauss points are read
   off whichever triangle contains them. Exact for a field that is linear on
   the split - which is what the export's own P1 space produces - and an
   approximation otherwise.

   \param source a Data on Function of the export
   \param target the oxley domain the export was built from
*/
escript::Data fromFinleyFunctionData(const escript::Data& source,
                                     escript::Domain_ptr target);

escript::Domain_ptr toFinley(const OxleyDomain& dom, int order = -1,
                             int reducedOrder = -1, bool optimize = false,
                             bool simplices = true);

/**
   \brief
   Diagnostic: reports whether a degree-2 p4est_lnodes numbers the hanging
   positions of this forest.

   A hanging position is the centre of a coarse octant face or the midpoint of a
   coarse octant edge. At degree 1 it is not a node at all, which is why the
   corresponding element_nodes slot holds a far master instead. At degree 2 that
   position IS an independent node of the coarse neighbour, so the slot may
   carry its global id - if so, the degree-2 numbering can be used directly and
   no separate numbering of hanging positions is needed.

   Test: every corner slot of every octant is mapped to its geometric position,
   and the id <-> position pairing is checked in both directions. A global id at
   two different positions means the slot holds something other than the node at
   that corner (a far master); one position under two different global ids means
   the position is numbered twice, which would leave a crack between the elements
   that disagree. Both counts must be zero for the numbering to be usable, and
   then [2], [6] and the number of distinct corner positions all coincide.

   \return a vector of counts:
           [0] octants examined
           [1] corner slots examined
           [2] distinct global ids seen at corner slots
           [3] global ids appearing at MORE THAN ONE position (must be 0)
           [4] octants whose face_code reports hanging faces
           [5] total local nodes in the degree-2 lnodes
           [6] distinct corner positions seen
           [7] positions carrying MORE THAN ONE global id (must be 0)
*/
std::vector<long> lnodesDegree2Report(const OxleyDomain& dom);

/**
   \brief
   One entry per element_nodes slot of a degree-2 lnodes: the global id the slot
   holds, and where the slot itself sits. Raw material for the diagnostic above -
   what the report reduces to counts, this leaves open to inspection from Python.
*/
struct SlotRecord
{
    long element;      ///< local element (leaf) index
    int slot;          ///< slot within the element, lexicographic (z)yx order
    long gid;          ///< global node id the slot holds
    double x, y, z;    ///< position OF THE SLOT (not of the node it holds)
    int faceCode;      ///< the element's lnodes face_code (0 = nothing hangs)
    bool corner;       ///< is this slot one of the element's corners?
};

/**
   \brief
   Dumps every degree-2 lnodes slot of the local forest. Diagnostic only.
*/
std::vector<SlotRecord> lnodesDegree2Slots(const OxleyDomain& dom);

} // namespace oxley

#endif // __OXLEY_FINLEYCONVERTER_H__
