# Changelog
All notable changes to this project will be documented in this file.

## [1.2.0] - 2026-08-08
### Added
- Export the bacterial cgMLST MST as an edge list (`*_MST.tsv`) and as a sample-level Newick tree (`*_MST.nwk`), rooted at the weighted graph centre determined with `networkx`.
- Show the MST as a third tree panel (`cgMLST MST`) in the bacterial Microreact project, next to the phylogenetic tree and the time tree.
### Fixed
- Keep sequence types separated by zero allelic differences connected in the MST; `scipy` had been discarding those edges and splitting the tree into a forest.
### Changed
- Name the Microreact tree panels `Phylogenetic tree` and `Time tree` instead of `Tree`.
- Load `calculate_allelic_distance_and_plot_MST.py` from `bin/` in the MST tests instead of keeping a second copy under `tests/MST_bacteria/`.

## [1.1.0] - 2026-07-28
### Fixed
- Limit MST distance-calculation workers to the number of resolved cgMLST profiles.

## [1.0.2] - 2026-07-27
### Changed
- Skip local CPU-count validation when the `slurm` execution profile is selected.

## [1.0.1] - 2026-07-27
### Changed
- Skip local Docker image validation when the `slurm` execution profile is selected.

## [1.0.0] - 2025-12-22
### Changed
- Mark this as the first production-ready version of the program.
- Updated documentation in the `doc` folder and the `README.md` file 
### Added
- Introduce the `CHANGELOG.md` file.