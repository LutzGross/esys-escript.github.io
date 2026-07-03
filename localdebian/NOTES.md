# localdebian — upstream's mirror of the Debian packaging

This tree mirrors the packaging that lives with the Debian Science team at
https://salsa.debian.org/science-team/python-escript and is uploaded to the
archive. Upstream cannot upload directly; changes here are handed to an
Uploader (Alastair McKinstry) to build/lintian-check and sponsor.

## Provenance

Reconciled against the archive source package **python-escript 6.1.1-1**
(experimental), i.e. the `python-escript_6.1.1-1.debian.tar.xz` from
deb.debian.org, so this is a faithful copy of what is actually in Debian —
*not* the earlier standalone `sid_*` modernization draft (that approach was
superseded by McKinstry's upload and has been dropped).

Key traits of the 6.1.1-1 packaging:
- debhelper-compat 14; Standards-Version 4.7.4; `3.0 (quilt)`.
- Two arch binaries `python3-escript` (OpenMP+MPI) and `python3-escript-mpi`,
  plus `python-escript-doc` (Arch: all).
- I/O keeps NetCDF (libnetcdf-c++4-dev); MPI via `mpi-default-dev`; mpi4py.
- The real builds are driven by `scons/templates/debian_nompi_options.py` and
  `debian_options.py` (see d/rules `override_dh_auto_build`). The
  `sid_options.py` / `sid_options_mpi.py` files here are vestigial (only a
  commented-out d/rules line references them) and are kept as-shipped.
- doc-base files are parked under `hide/` (not installed as-is).

## Local delta on top of 6.1.1-1 (see d/changelog `6.1.1-2 UNRELEASED`)

- **Architecture restricted to 64-bit** for both `python3-escript` and
  `python3-escript-mpi`. esys-escript is 64-bit-only and never built on the
  32-bit ports: 6.1.1-1 FTBFS on armhf and i386 at the scons compiler sanity
  check (`conf.CheckCXX` in `site_scons/dependencies.py`), and 5.6-10 failed
  there too. All 64-bit release arches (amd64, arm64, ppc64el, s390x) plus the
  64-bit ports built and installed fine.

This delta is the change proposed to McKinstry for the next revision; see
`../debian-mckinstry-handoff.md`.