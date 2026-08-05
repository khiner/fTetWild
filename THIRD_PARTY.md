# Third-party code

fTetWild is MPL 2.0 (`LICENSE.MPL2`). It has no third-party dependencies at build time beyond
CLI11 (binary only) and Catch2 (tests only); everything else the algorithm needs has been brought
in-tree. The files below came from other projects and stay under their own terms.

They live in `src/` alongside first-party code rather than in an `external/` directory. This file
is what keeps their provenance findable. Each one also repeats its licence in its own header, and
that per-file notice is the actual grant: libigl, for instance, ships both `LICENSE.GPL` and
`LICENSE.MPL2`, and only the per-file headers say which applies.

## Files

| File | Upstream | Licence |
|---|---|---|
| `src/MeshCleanup.hpp`, `src/MeshCleanup.cpp` | libigl — `igl/{sort,sortrows,round,unique_rows,remove_duplicate_vertices,vertex_components,orientable_patches,bfs_orient}` | MPL 2.0 |
| `src/writeOBJ.h` | libigl — `igl/writeOBJ.h` | MPL 2.0 |
| `src/FastWindingNumber.h`, `src/FastWindingNumber.cpp` | [Side Effects Software](https://github.com/sideeffects/WindingNumber), via libigl | MIT |
| `src/predicates.c` | [Jonathan Richard Shewchuk](https://www.cs.cmu.edu/~quake/robust.html), CMU | Public domain |
| `src/AABBWrapper.h`, `src/AABBWrapper.cpp` (the `MeshFacetsAABBWithEps` tree only) | [geogram](https://github.com/BrunoLevy/geogram) — Bruno Levy, Inria | BSD 3-clause |
| `src/geo/*` | [geogram](https://github.com/BrunoLevy/geogram) — Bruno Levy, Inria | BSD 3-clause (`src/geo/LICENSE.geogram`) |
| `src/triangle_triangle_intersection.cpp` | Guigue & Devillers, *Journal of Graphics Tools* 8(1), 2003 | **none stated** |
| `src/getRSS.c` | David Robert Nadeau, [NadeauSoftware.com](http://NadeauSoftware.com/) | CC BY 3.0 |
| `src/MeshIO.cpp` (the Gmsh 2.2 writer only) | [PyMesh](https://github.com/PyMesh/PyMesh) — Qingnan Zhou | **none — see below** |

`src/Predicates.hpp` and `src/Predicates.cpp` are fTetWild's own wrapper over Shewchuk's
`predicates.c`, not vendored code.

## Obligations

- **MPL 2.0** is file-level copyleft, and it follows the code through edits and rewrites. This costs
  nothing today since fTetWild is MPL 2.0 itself, but it does mean these files can never be
  relicensed permissively. §3.3 is what permits the BSD/MIT/public-domain files alongside them in
  one Larger Work.
- **BSD 3-clause** and **MIT** require the copyright notice, condition list and disclaimer to travel
  with the source. BSD 3-clause also forbids using "Inria" or the "ALICE Project-Team" name to
  promote anything derived from it.
- **CC BY 3.0** requires attribution. Creative Commons discourages CC licences for software and
  this one is not OSI-approved, so `getRSS.c` is the untidiest entry here — though the obligation
  itself is only attribution.
- **Public domain** (`predicates.c`) carries no obligation. The attribution is kept because it
  documents the algorithm.

## Unresolved

Two files carry no grant of rights at all. Both predate this tree — upstream fTetWild ships them the
same way — and neither is a new exposure, but neither is settled either.

- **The Gmsh writer in `MeshIO.cpp` (PyMesh).** It carries a bare copyright line, and the PyMesh
  repository has no `LICENSE` file; the GitHub API reports its licence as `null`. This is the weaker
  of the two. It is also the easy one to fix: Gmsh MSH 2.2 is a published format, and a writer
  written from the spec would owe nothing to Zhou's code.
- **`triangle_triangle_intersection.cpp` (Guigue–Devillers).** Companion code to a published JGT
  paper, with a citation but no licence. Weaker concern in practice — it has been redistributed for
  two decades — and rewriting it is the wrong move regardless: the routine is numerically delicate,
  and a from-the-paper reimplementation would risk the byte-identical corpus match for no real gain.

Rewriting does not clear a licence in general. A rewrite of licensed code is a derivative work and
the obligations follow it; only a genuinely independent reimplementation from a specification
escapes. That is why the Gmsh note above is about writing from the published spec, not about editing the
existing code until it looks different.
