    # Oxley Reimplementation — Design

Status: **draft for discussion**. Supersedes the current `oxley/` domain/DOF/refinement
layer. Grounded in the reuse audit (see "Reuse map") and the API sketches
`oxley-designs/simple3D.py` and `oxley-designs/tools.py`.

Tags below: **[DECIDED]** = agreed direction, **[OPEN]** = needs a decision before/while implementing.

---

## 1. Motivation

Oxley wraps p4est/p8est to give escript adaptive (octree) FEM domains. The current
implementation is prototype-quality and has four *architectural* problems that
incremental cleanup cannot reach:

1. **Mutable domains.** `domain.refine(...)` mutates in place. escript `Data` binds to a
   domain's function spaces and caches sample/DOF counts, so in-place refinement silently
   invalidates every existing `Data` on that domain:
   `x=domain.getX(); domain.refine(); x+domain.getX()` → size mismatch → crash.
2. **Node identity by floating-point coordinate hashing.** Global nodes are keyed by
   `unordered_map<tuple<double,double,double>, long>` (~100+ call sites). Fragile (FP
   identity), non-scalable, and the root cause of the broken parallel path.
3. **Hanging nodes eliminated by a global sparse triple product.** The stiffness matrix is
   assembled over *all* lnodes including hanging ones, then condensed via `IZᵀ·A·IZ`
   (`finaliseAworker` → `IztAIz`). Inefficient, and the eager `IZ` construction
   (`initZ`/`initIZ` in the constructor) is the MPI abort site (Tpetra `CrsGraph`).
4. **Construction conflates blocks with elements.** `n0/n1/n2` are documented as "elements"
   but build p4est *trees* (`num_trees = n0·n1[·n2]`), each holding a single leaf — so the
   macro-mesh *is* the element mesh and p4est's "start coarse, adapt deep" advantage is lost.
   The API (a ripley clone) also exposes no way to set a structured base resolution or to
   drive refinement from geometry/features.

The result is a domain that works serially for simple uniform-ish cases but is "not useful
for practical applications." This document specifies a reimplementation of the domain/DOF/
refinement **middle layer** on top of the (sound) p4est layer and assembler kernels.

## 2. Goals & non-goals

**Goals**
- Immutable domains; refinement produces a **new** domain (functional).
- Refinement driven by **geometry/feature tasks** that also **tag** the regions they create,
  so material properties and BCs attach naturally (the "practical applications" gap).
- Correct **hanging-node** handling in **2D and 3D**, done at the element level.
- Correct **parallel** (MPI) node/DOF distribution.
- **Block-structured construction**: a coarse grid of *blocks* (p4est trees) with an optional
  per-block initial refinement level and background tag — a structured coarse model of the
  domain — that the feature-driven `Refiner` then adapts. Dimension-agnostic; explicit
  `origin`.
- Reuse the p4est layer and the conforming assembler kernels.

**Non-goals (for the first cut)**
- Higher element/integration order (stays fixed 2-point Gauss, integration order 3).
- De-refinement/coarsening of an existing domain (refinement only makes finer domains).
- Anisotropic refinement (p4est refines isotropically; base mesh may be anisotropic).

## 3. User-facing API

### 3.1 Construction — block-structured  [DECIDED]

A domain is built from a coarse grid of **blocks**. Each block is one p4est **tree**
(macro-cell); its leaves are the elements. This exposes p4est's real two-level structure
(blocks → refined leaves) instead of the current one-tree-per-element collapse (§1.4).

Dimension-agnostic primary constructor (2D or 3D from the length of `numBlocks`):

```python
Block(numBlocks=(20, 20, 10),        # number of BLOCKS (p4est trees) per axis
      length=(1000., 1000., 800.),   # physical extent per axis
      origin=(-500., -500., -800.),  # coordinate of the lower corner
      refine_level=0,                # uniform refinement of every block (see below)
      block_tags=None,               # background region tag per block (see below)
      comm=None,                     # mpi4py communicator (default COMM_WORLD)
      framework=None)                # SolverFramework (Paso/Trilinos)
```

`Rectangle`/`Brick` remain dimension-specific conveniences:

```python
Rectangle(numBlocks0, numBlocks1, length0=1., length1=1., origin0=0., origin1=0.,
          refine_level=0, block_tags=None, comm=None, framework=None)
Brick(numBlocks0, numBlocks1, numBlocks2, length0=1., length1=1., length2=1.,
      origin0=0., origin1=0., origin2=0., refine_level=0, block_tags=None, ...)
```

**Resolution model.** `refine_level = L` uniformly refines every block `L` times, so:
- elements per axis within a block = `2^L`;
- base element size along axis *i* = `length_i / (numBlocks_i · 2^L)`.

The `Refiner` (§3.2) then adds *adaptive* refinement above this uniform base.

**Per-block level and tag — scalar-or-array (broadcast) semantics.** [DECIDED; arrays are a
later phase — see §7]
- `refine_level`: an `int` (uniform) **or** an array of shape `numBlocks` giving a per-block
  initial level (structured pre-refinement — e.g. finer near the surface).
- `block_tags`: `None` **or** an array of shape `numBlocks` of tag **names** (strings, mapped
  through the escript tag map; integer ids also accepted) — a **background region tag** per
  block.

This makes `Block` a *structured coarse model → adaptive FEM domain* converter: a gridded
earth model (a material-id array, a resolution array) drops straight into `block_tags` /
`refine_level`, and the `Refiner` resolves embedded features on top.

**Consequence of per-block levels.** Differing levels across neighbouring blocks make the
mesh non-conforming *at block boundaries from the start* → hanging nodes and 2:1 balance
apply to the base mesh, not just after adaptive refinement, and blocks differing by more than
one level trigger balancing cascades (§4.7). The hanging-node machinery (§4.3) must therefore
be correct for the base mesh too.

Notes:
- **Blocks are the macro-mesh, not the MPI partition.** p4est partitions *leaves* (elements)
  across ranks by space-filling curve (§4.5), so a rank owns a slice of leaves that may span
  several blocks or part of one — hence "blocks", not "subdomains".
- `origin` + `length` replace the current `l0`-as-2-tuple overload.
- Legacy `d0/d1/d2`, `periodic*`, `order` are **not** accepted (already removed; the shim
  warns). Decomposition is p4est's; periodicity is a separate future feature.

### 3.2 Refinement — functional, feature-driven, tagging  [DECIDED shape, some OPEN details]

```python
refiner = Refiner(levels_max=5)                              # additional levels above base
refiner.add(Sphere(center=(0,0,-400), radius=200,
                   tagname="Anomaly", resolution=5.))        # refine volume, tag it
refiner.add(PlaneInterface(origin=(0,0,400), normal=(0,0,1),
                           tagname="Deep", resolution=10.))  # refine surface, tag below
domain2 = refiner(domain1)     # NEW domain; domain1 (and its Data) untouched
```

- `Refiner(levels_max)` accumulates tasks; `refiner(domain)` (== `refiner.refine(domain)`)
  returns a new immutable domain.
- **Task model** (`RefinementTask`): each task has
  - `resolution` — target element **size**; refine while element size `>` resolution
    (fixes the inverted `<=` in the sketch);
  - `tagname` — region tag applied to elements it selects (optional);
  - `isinterface` — `False`: refine elements fully inside the region; `True`: refine
    elements the surface passes through;
  - `check(x) -> bool` (vectorized over coordinate arrays) — inside/below test.
- Built-in tasks: `Sphere(center, radius, ...)`, `PlaneInterface(origin, normal, ...)`.
  Extensible: `Box`, `Cylinder`, `PolySurface`, `FieldThreshold(data, value)` (for
  solution-driven adaptivity), etc.
- **Levels are relative to the base** [DECIDED]: `levels_max` bounds the number of
  *additional* refinement levels above the uniform base (§3.1), so the absolute p4est level
  cap is `max(refine_level) + levels_max`. A task's `resolution` is an absolute element size;
  refinement stops when `element_size ≤ resolution` or the cap is hit.
- **Two-layer tagging** [DECIDED]:
  1. *Block tags* (§3.1) give each element a coarse structured **background** region.
  2. *Task tags* label embedded features and **override** the block (background) tag on the
     elements a task applies to. Children inherit the parent element's tag unless a task
     re-tags them.
  Element tags become the escript tag map, usable as `kappa.setTaggedValue("Anomaly", ...)`.
- **[OPEN]** Whether *interface* tasks re-tag (volume tasks clearly tag their interior; the
  "below/inside" side for an interface is ambiguous). Task-vs-task overlap precedence
  (add-order vs explicit priority) also TBD.
- **[OPEN]** Boundary/face tags (e.g. `"bottom"`). Base-mesh outer faces get canonical tags
  (`left/right/bottom/top/front/back`); how interface tasks interact with face tags TBD.

### 3.3 Function spaces & data  [DECIDED semantics]

On the returned domain the usual escript function spaces behave as:

| Function space | Points per element | Meaning |
|---|---|---|
| `ContinuousFunction`, `ReducedContinuousFunction` | nodes | all lnodes (incl. hanging, values *derived*) |
| `Function` | 4 (2D) / 8 (3D) | 2-pt Gauss |
| `ReducedFunction` | 1 | element centre |
| `FunctionOnBoundary` | 4 (3D faces) / 2 (2D edges) | boundary Gauss |
| `ReducedFunctionOnBoundary` | 1 | boundary element centre |
| `Solution`, `ReducedSolution` | DOFs | **non-hanging** nodes; labelling depends on framework |

**Invariant [DECIDED]:** *DOFs are the non-hanging nodes.* Hanging-node values are never
independent unknowns — they are derived from their masters via the constraint table (§4.3).

Cross-domain interpolation (`interpolate(x_on_domain2, ContinuousFunction(domain1))`) uses
`interpolateAcross`; combining `Data` from two different domains without interpolation raises
a clean "different domains" error (domain identity via `p4est_checksum`, already correct).

### 3.4 PDEs

Unchanged escript usage; the domain supplies matrices/vectors through the assembler:
```python
mypde = LinearSinglePDE(domain2)
kappa = Scalar(1., Function(domain2)); kappa.setTaggedValue("Anomaly", 10)
mypde.setValue(A=kappa*kronecker(domain2), y=..., q=whereZero(Solution(domain2).getX()[2]))
u = mypde.getSolution()
```

## 4. Core architecture

Three layers; the reimplementation replaces the middle.

```
  User API (Python)        Block/Rectangle/Brick, Refiner, tasks, tagging     [NEW]
  ------------------------------------------------------------------------
  Domain / DOF / Assembly  lnodes numbering, hanging table, condensation,     [REWRITE]
                           MPI ownership, function spaces, interpolation
  ------------------------------------------------------------------------
  p4est layer + kernels    connectivity, forest/ghost, 2:1 balance,           [REUSE]
                           2-pt Gauss element kernels, hanging detection
```

### 4.1 Global numbering via `p4est_lnodes`  [DECIDED]

Replace coordinate-hash node identity with **`p4est_lnodes` (degree 1)**, which provides:
- a global, topological node numbering with an **owned/shared (ghost) partition** per rank;
- `element_nodes` (the node ids of each element's corners);
- `face_code` per element — already used for hanging detection.

Node ids come straight from lnodes; no floating-point keys. DOF ids are derived by
compacting out hanging nodes (§4.3).

### 4.2 Hanging detection  [REUSE]

Keep the existing `face_code` decoders (`isHangingNode`/`getHangingInfo`, 2D and 3D). They
already turn an element's `face_code` into which local nodes/edges/faces hang.

### 4.3 Hanging → master constraint table  [DECIDED]  ← keystone

A single table, built in one element sweep from `lnodes.element_nodes` + `face_code`:

```
constraints:  hanging_node_gid  →  [ (master_gid, weight), ... ]
```

Standard interpolation weights: 2D edge-midpoint = ½·(2 endpoints); 3D edge-midpoint =
½·(2), face-centre = ¼·(4 corners). With lnodes numbering the masters already have global
ids, so the table is cheap and exact.

**DOF map**: `node_gid → dof_gid` assigns a DOF to each *non-hanging* node and marks hanging
nodes as "constrained" (pointing into `constraints`). DOF count = non-hanging node count.

This one table serves **both**:
- **Assembly** (§4.4) — element-level condensation.
- **Interpolation / nodal data** (§4.6) — hanging value = Σ weight·master value.

Retires `Z`/`IZ`, `initZ`/`initIZ`, `makeZ`/`makeIZ`, `IztAIz`, and makes the empty 3D
`assemblePDEHanging` stub unnecessary.

### 4.4 Element-level condensation assembly  [DECIDED]

Per element (reusing the existing 2-pt Gauss kernels to form the local matrix `A_e` and load
`b_e` over the element's lnodes):
1. Build the small element constraint operator `C_e` from `constraints` (identity rows for
   real local nodes; interpolation rows for hanging local nodes).
2. Condense locally: contributions of hanging local DOFs are distributed onto their masters
   (`A_e' = C_eᵀ A_e C_e`, `b_e' = C_eᵀ b_e`) — i.e. deal.II-style
   `distribute_local_to_global`.
3. Scatter **only real DOFs** into the global matrix/RHS.

No oversized global matrix; no global triple product.

### 4.5 Parallel (MPI)  [DECIDED approach]

- p4est distributes elements across ranks by space-filling curve (`p4est_new_ext` on the
  domain communicator + `p4est_partition`); no `d0×d1` block grid.
- lnodes gives each rank its **owned** nodes and its **ghost** (shared) nodes.
- Real DOFs are owned per rank → Trilinos/Paso **row map** = owned real DOFs.
- A hanging node whose masters are owned by another rank contributes to those masters as
  **ghost columns** — the standard Trilinos column-map/import pattern. No bespoke `IZ`.
- `ownSample` for Elements/FaceElements must be implemented (currently throws).

### 4.6 Interpolation & nodal data  [DECIDED]

- Reading nodal `Data` (e.g. `ContinuousFunction(domain).getX()`): fill real nodes directly;
  fill hanging nodes from `constraints` (Σ weight·master). Makes continuous fields consistent
  across non-conforming faces.
- `Solution → ContinuousFunction`: scatter DOFs to real nodes, then apply `constraints` to
  hanging nodes.
- `interpolateAcross` (domain→domain): complete the existing implementation for the missing
  function-space pairs; exploit the parent/child relationship when `domain2` is a refinement
  of `domain1` (interpolation can be exact along shared coarse elements).

### 4.7 Refinement engine  [DECIDED shape; performance OPEN]

`refiner.refine(domain)`:
1. Copy the base forest (leaves `domain` untouched).
2. Register task tags in the copy's tag map.
3. Iterate up to `levels_max` passes. Each pass, using a p4est refine callback:
   - a leaf is a **candidate** if its size `>` some task's `resolution`;
   - evaluate the relevant task `check()` at the candidate's (would-be child) coordinates;
   - **volume task**: mark for refinement if all children inside; tag children/parent by side;
   - **interface task**: mark for refinement if children straddle the surface;
   - stop when no leaf was refined in a pass.
4. `p4est_balance` (2:1) after refinement, capped at `max(refine_level) + levels_max`. The
   base mesh is likewise balanced at construction when per-block `refine_level`s differ
   (§3.1). **[OPEN]** balance each pass vs once at the end — balancing cascades extra
   refinement that can exceed a task's `resolution` locally.
5. `p4est_partition`; rebuild lnodes, constraints, DOF map, tags.
6. Construct and return the new domain.

**[OPEN] Where `check()` runs / performance.** Tasks are Python objects; calling `check()`
per element per pass from a C++ callback is too slow at scale. Proposed first cut: the pass
gathers candidate coordinates into a numpy array and calls `check()` **vectorized** (tasks
operate on arrays), then feeds the boolean result back to the p4est refine decision. Later
optimization: C++-native criteria for the common shapes.

### 4.8 Output — VTK / Silo  [DECIDED]

Output uses escript's **standard `weipa` path** (`saveVTK`, `saveSilo`) — identical to
finley/ripley/speckley — via the weipa Oxley adapter (`weipa/src/Oxley{Domain,Elements,
Nodes}`). Each leaf element is written as one unstructured cell (quad/hex). Silo needs the
`silo` build feature; VTK (`.vtu`) needs none. There is **no build cycle**: weipa depends on
the domain libraries, not vice-versa (`escriptcore → domains → weipa`).

**Decouple the adapter from oxley internals.** The current weipa Oxley adapter re-walks the
p4est forest itself (`p4est->trees`, `p4est_quadrant_array_index`, …) and re-derives node ids
through the domain's coordinate-hash map (`NodeIDs.find(make_pair(x,y))`,
`getNeighouringNodeIDs`). That makes the node numbering **co-owned by oxley and weipa**: any
change to oxley's numbering breaks the adapter, and the two must move in lockstep — a design
loop (not a build loop). Ripley/finley avoid this by having their adapters read the domain's
public node/element tables rather than re-deriving topology; oxley's adapter is the outlier.

Fix, as part of **Milestone A4**: oxley exposes a **public mesh-access interface** —
per-leaf element→node-id connectivity, node global ids + coordinates, the node→DOF map, and
tags — and the weipa Oxley adapter is rewritten to consume **only** that interface, with no
knowledge of p4est or the numbering scheme. The lnodes numbering (§4.1) then lives entirely
inside oxley, and the numbering rewrite stays contained (does not ripple into weipa beyond a
one-time adapter rewrite against the new interface). Also drop the adapter's leftover debug
`std::cout`s.

The bespoke `writeToVTK`/`saveMesh` (p4est's native `p4est_vtk`) is **demoted to an optional
developer/debug** dump of the forest (refinement level, MPI partition), not user-facing
output.

**Milestone B (non-conforming meshes):** unstructured VTK/Silo writes each leaf as its own
cell, so hanging-node meshes export directly. Nodal fields must write **constrained
(interpolated) values at hanging nodes** (§4.6) so continuous fields render without cracks —
no change to the file format, just correct hanging values.

## 5. Reuse map

| Component | Verdict |
|---|---|
| p4est/p8est connectivity builders, forest/ghost creation, 2:1 balance | **Reuse** |
| Domain identity `operator==` via `p4est_checksum` | **Reuse** |
| Conforming assembler kernels (2-pt Gauss, per-level `dx`), 2D & 3D `assemblePDESingle` | **Reuse (kernels)** |
| Hanging-node **detection** (`face_code` decoders) | **Reuse** |
| Global node/DOF numbering (coordinate hashing) | **Replace → lnodes** |
| Hanging elimination (`Z`/`IZ`, `IztAIz`, `initZ`/`initIZ`) | **Remove → element condensation** |
| 3D `assemblePDEHanging` (empty stub) | **Remove (moot)** |
| MPI ownership/distribution, `ownSample` | **Rewrite** |
| `interpolateAcross`, tagging | **Salvage & complete** |
| Domain/refinement API (in-place mutators) | **Replace → immutable functional** |

## 6. Function-space / DOF summary

- **Nodes**: lnodes (degree 1), global topological numbering, owned/ghost partition.
- **DOFs**: non-hanging nodes; `node→dof` map compacts out hanging nodes.
- **Hanging nodes**: derived via `constraints` (never independent unknowns).
- **Elements/boundary**: 2-pt Gauss (unchanged); `ReducedFunction` = centre.

## 7. Phased implementation plan

Sequenced into two milestones. **Milestone A first**: a *uniform, conforming* Block domain
that reproduces ripley/finley functionality — no adaptive refinement. A uniform
`refine_level` yields a fully conforming mesh (**no hanging nodes**), so the hanging-node
constraint table, the element-level condensation (degenerates to a plain scatter), the 3D
hanging assembler, the `Refiner`, and `interpolateAcross` are all **deferred to Milestone B**.
Milestone A still delivers the core rewrite — lnodes numbering and the MPI DOF distribution —
because the current code is broken there even for uniform meshes.

### Milestone A — conforming uniform domain (ripley/finley parity)  ← current focus

A1. **Construction.** `Block(numBlocks, length, origin, refine_level, comm, framework)` with a
    single **scalar, uniform** `refine_level` (conforming); `Rectangle`/`Brick` delegating.
    (`block_tags` / per-block arrays deferred to B.)
A2. **Mesh-access interface + numbering.** lnodes-based node numbering; expose a **public
    mesh-access interface** — per-leaf element→node connectivity, node global ids +
    coordinates, tags (and the node→DOF map, trivial while conforming). Replaces the
    coordinate-hash scheme; this is the interface both output and assembly consume (§4.8).
A3. **VTK / Silo output.**  ← *next after Block, per the output-early strategy.* Rewrite the
    weipa Oxley adapter on the A2 interface (no p4est/`NodeIDs` internals) so `saveVTK` /
    `saveSilo` work (§4.8); demote native `p4est_vtk` to a debug dump. This is the first
    **visual validation** of construction — block grid, `refine_level`, `origin`, tags — and
    a fast feedback loop for A1/A2 before any assembly exists.
A4. **Assembly.** Conforming element assembly reusing the 2-pt Gauss kernels; scatter real
    DOFs directly. Drop `Z`/`IZ`/`IztAIz`. Function spaces + `getDataShape`.
A5. **Nodal data & interpolation.** `getX`, within-domain `interpolate`,
    `Solution↔ContinuousFunction`, outer-face canonical tags, dirac points.
A6. **MPI.** Owned/ghost DOF maps + `ownSample`; parallel assembly. **Fixes the current
    `CrsGraph` abort.** Verify parallel solve matches serial.

*Acceptance:* the existing oxley Python test suite (which mirrors the ripley/finley tests —
`run_escriptOnOxley`, `run_linearPDEsOnOxley`, `run_utilOnOxley`, solver tests, …) passes
**serial and parallel**, and a reference Poisson / elasticity solve matches ripley.

### Milestone B — adaptivity

B1. **Hanging-node constraint table** (from lnodes `face_code`) + element-level condensation,
    2D then 3D. Verify a Poisson solve on a hanging-node mesh vs. analytic/finley reference.
B2. **Per-block `refine_level` / `block_tags` arrays** (introduce base-mesh hanging + the
    two-layer tagging model).
B3. **Refinement engine + tasks.** `Refiner`, `RefinementTask`, `Sphere`, `PlaneInterface`,
    functional `domain2 = refiner(domain1)`, feature tagging.
B4. **Cross-domain interpolation.** Complete `interpolateAcross`; solution-driven adaptivity
    (`FieldThreshold` task).

## 8. Risks / open questions

- **[OPEN]** Balance vs `resolution`/`levels_max` interaction (§4.7 step 4).
- **[OPEN]** Refinement `check()` performance & where it runs (§4.7).
- **[OPEN]** Interface-task re-tagging, and task-vs-task overlap precedence (§3.2); boundary/
  face tagging (§3.2).
- **[RISK]** Correctness of 3D face+edge hanging constraints (the current 3D path never
  implemented this) — needs careful verification in step 2/6. Per-block base levels (§3.1)
  exercise these constraints on the base mesh, so they must be right early.
- **[RISK]** Parallel constraints where masters are off-rank — must be covered by ghost
  columns and tested early (step 6).
- **[OPEN]** Solution-driven adaptivity (`FieldThreshold` task) — fits the same engine but
  needs an error/indicator field API; deferred.
