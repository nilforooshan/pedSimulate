# R package pedSimulate: Pedigree, genetic merit and phenotype simulation

## Version: 0.0.3

## Version: 0.0.4

- Applying mortality to each generation separately rather than on selection candidates, which may come from more than a single generation.

## Version: 0.1.1

- Added function `assortative` for assortative/disassortative mating.
- Added `"-P"` and `"-PA"` to arguments `fsel` and `msel` to select in opposite directions of `"P"` (phenotype) and `"PA"` (parent average).

## Version: 0.1.2

- Debugged the bug introduced in the previous version in the `simulatePed` function.

## Version: 1.0.0

- Major revision
- Functions `simulatePed` and `assortative` were combined into `simulatePed`, and arguments `fullsib` and `parentprogeny` were removed.
- Added function `appendPed` for simulating new generations from an existing pedigree and appending to it.
- Added function `fs_mate_finder` for finding fullsib matings in the pedigree.
- Added function `hs_mate_finder` for finding halfsib matings in the pedigree.
- Added function `pp_mate_finder` for finding parent-progeny matings in the pedigree.

## Version: 1.0.1

- Made a small debug.

## Version: 1.1.0

- Added agrumnts `f.order` and `m.order` to `simulatePed` and `appendPed` functions.

## Version: 1.1.1

- Debugged for sampling from a sample size of 1.

## Version: 1.2.0

- Added function `simulateGen` for simulating genotypes.

## Version: 1.2.1

- Fixed a small bug in `hs_mate_finder`.

## Version: 1.3.0

- A major debug: `%in% (-overlap.s:0)+i-1` was performing differently from `%in% ((-overlap.s:0)+i-1)`, similarly `%in% (-overlap.d:0)+i-1` and `%in% ((-overlap.d:0)+i-1)`.
- Added function `appendGen` for simulating genotypes for an appended pedigree to an existing pedigree with genotypes.

## Version: 1.3.1

- A small (same) debug in `simulatePed` and `appendPed`.

## Version: 1.3.2

- Added `seed` argument to functions for reproducible outputs.

## Version: 1.4.0

- Debugged function `appendPed` for making no updates to the existing pedigree (the pedigree upon which new generations are appended).
- "Bijma, P. & Rutten, M. (2002)" gone missing on the web. References were added.

## Version: 1.4.1

- Fixed the naming issue `pedSimulte-package` in the documentation.

## Version: 1.4.3

- `man/pedSimulte-package.Rd` was missing in the previous version.
