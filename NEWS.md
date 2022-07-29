# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.6.1] - 2022-07-29

### Fixed

- Handling of periodic boundaries, when two periodic boundaries intersect in a common vertex or edge.

## [0.6.0] - 2022-02-11
### Changed
- Substantial change in the dependencies to allow parallel computations.

### Added
- Support for parallel computations via `GridapDistributed`.

## [0.5.0] - 2021-11-24
### Changed
- Restrict to Gridap 0.17 as required for surface meshes.

## [0.4.5] - 2021-11-24
### Fixed
- Some tests failing for Gridap 0.16.

## [0.4.4] - 2021-11-22

### Added
- Support for surfaces meshes.

## [0.4.3] - 2021-10-28

### Added
- Support for Gridap v0.17.

## [0.4.2] - 2021-06-08

### Added
- Support for Gridap v0.16.

## [0.4.1] - 2021-01-22

### Added
- Support for periodic boundary conditions.
- Automatic installation via `BinaryBuilder.jl`.

## Previous

A changelog was not maintained for previous versions.
