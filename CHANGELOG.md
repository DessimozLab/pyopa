# Changelog

All notable changes to pyopa are documented in this file. Versions that
were tagged but never published on PyPI are marked as such; their
changes are released with the following version.

## 0.9.0 - 2026-10-30

### Added
- `MutipleAlEnv.estimate_many_pam` estimates the PAM distances of all
  pairs of a set of aligned sequences in parallel with OpenMP and
  returns them in condensed pair order, as `scipy.spatial.distance.pdist`
  ([#28](https://github.com/DessimozLab/pyopa/pull/28)).
- `AlignmentProfile.align_double` takes an explicit `threshold` argument
  ([#29](https://github.com/DessimozLab/pyopa/pull/29)).
- Wheels for Linux aarch64 and for Python 3.14 and 3.15
  ([#27](https://github.com/DessimozLab/pyopa/pull/27)).
- Wheels for Windows x86_64 and arm64 (Python 3.11 or newer on arm64, as
  numpy has no Windows arm64 wheels for Python 3.10); Windows is tested in
  CI ([#33](https://github.com/DessimozLab/pyopa/pull/33)).

### Changed
- Requires Python 3.10 or newer; no more wheels for Python 3.6 to 3.9.
  Package metadata moved from `setup.py` to `pyproject.toml`
  ([#27](https://github.com/DessimozLab/pyopa/pull/27)).
- `MutipleAlEnv.estimate_pam` and `estimate_many_pam` are about 30 times
  faster: they evaluate `exp(pam * logPAM1)` from an eigendecomposition
  of logPAM1, only for the pairs of symbols in the alignment. Scores and
  PAM distances are unchanged, variances agree to about 1e-9 relative.
  If a custom logPAM1 is not reversible, the original implementation is
  used. logPAM1 is now read only when the `MutipleAlEnv` is created
  ([#30](https://github.com/DessimozLab/pyopa/pull/30)).
- Wheels are built against numpy 2.x and work with numpy 1.x and 2.x
  ([#27](https://github.com/DessimozLab/pyopa/pull/27)).
- The extension is built with OpenMP where the compiler supports it, and
  falls back to a serial build otherwise; set `PYOPA_REQUIRE_OPENMP=1` to
  make the build fail instead. On macOS, libomp is taken from
  `OPENMP_PREFIX` or Homebrew; the macOS wheels bundle a libomp built for
  macOS 11.0 ([#28](https://github.com/DessimozLab/pyopa/pull/28)).
- `-march=native` is no longer used, so wheels no longer depend on the
  CPU of the build machine ([#28](https://github.com/DessimozLab/pyopa/pull/28)).
- `create_scaled_matrices` is about 6 times faster, which halves the time
  of `generate_all_env`
  ([#29](https://github.com/DessimozLab/pyopa/pull/29)).
- `align_double(..., calculate_ranges=True)` stops the reverse alignment
  at the start of the alignment, which makes `align_strings` about 13 %
  faster ([#29](https://github.com/DessimozLab/pyopa/pull/29)).
- CI and wheel builds use uv and cibuildwheel 4
  ([#27](https://github.com/DessimozLab/pyopa/pull/27)).

### Fixed
- `EstimatePam` uses much less stack, so it can run in threads with small
  stacks such as on musl libc
  ([#28](https://github.com/DessimozLab/pyopa/pull/28)).
- `MutipleAlEnv.free` can be called repeatedly and refuses to free the
  matrices while `estimate_many_pam` is running; using freed matrices
  raises a `RuntimeError` instead of crashing
  ([#28](https://github.com/DessimozLab/pyopa/pull/28)).
- `align_double` truncated pointers on 64-bit Windows, where `long` is
  32 bits ([#33](https://github.com/DessimozLab/pyopa/pull/33)).
- Compiler warnings about implicit conversions with MSVC ([#33](https://github.com/DessimozLab/pyopa/pull/33)).

## 0.8.6 - 2025-04-24

### Fixed
- Publishing to PyPI.

## 0.8.5 - 2025-04-24 (not published on PyPI)

### Added
- Support and wheels for Python 3.12 and 3.13.

### Changed
- Builds with Cython 3.
- The classifiers list Python 3.8 to 3.13 instead of Python 2.7 and 3.4
  to 3.7 (wheels are still built for Python 3.6 and 3.7). CI no longer
  tests Python 3.7.
- Removed the Travis CI configuration and the old Docker based wheel
  build script.

## 0.8.4 - 2023-07-28

### Changed
- Linux aarch64 wheels are no longer built, as their build did not
  finish in CI.

## 0.8.3 - 2023-07-28 (not published on PyPI)

### Added
- Wheels for Python 3.11.
- Wheels for Linux aarch64.

### Changed
- Pinned the build requirements (setuptools, wheel, Cython < 3.0) in
  `pyproject.toml`.

## 0.8.2 - 2022-04-25

### Fixed
- Publishing to PyPI.

## 0.8.1 - 2022-04-25 (not published on PyPI)

### Added
- Support for ARM64, including Apple Silicon, by translating the SSE
  intrinsics with SIMDe.
- Wheels for macOS x86_64 and arm64, Linux manylinux2014 and musllinux,
  for Python 3.6 to 3.10, built and tested with cibuildwheel on GitHub
  Actions.

### Changed
- CI moved from Travis CI to GitHub Actions.

### Fixed
- The SIMDe headers are included in the sdist.
- A compiler warning in `EstimatePam`.

## 0.8.0 - 2018-11-14

### Added
- Wheels for Linux and macOS for Python 2.7 and 3.4 to 3.7.
- `pyproject.toml` declaring the build dependencies.
- The text of the MPL 2.0 license.

### Changed
- Restructured as the `pyopa` package with a Cython backend
  (`pyopa.backend`) and the default matrices included in the package.
- Avoids large stack allocations in the C code.
