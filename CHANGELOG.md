# Changelog
All notable changes to this project will be documented in this file.

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