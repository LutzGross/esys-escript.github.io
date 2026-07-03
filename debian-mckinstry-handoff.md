Subject: esys-escript 6.1.1 — Debian packaging ready for sbuild + upload

Hi Alastair,

esys-escript 6.1.1 is released upstream, and I've prepared the Debian
packaging update (5.6-10 → 6.1.1-1). As an Uploader on python-escript,
would you be able to build/lintian-check and sponsor the upload?

Upstream tarball (uscan-trackable via the updated d/watch):
  https://github.com/LutzGross/esys-escript.github.io/archive/refs/tags/6.1.1.tar.gz

The packaging delta is attached as a single patch (debian-6.1.1-update.patch),
verified to apply cleanly on top of the current debian/latest (5.6-10). Summary
of the changes:

  * New upstream release 6.1.1.
  * I/O: enable HDF5 alongside NetCDF (escript 6.x supports both); add
    libhdf5-dev to Build-Depends.
  * MPI build now enables mpi4py: python3-mpi4py added to Build-Depends and to
    python3-escript-mpi Depends; sid_options_mpi.py sets mpi4py=True.
  * Drop the esys.downunder module + its inversion guide (removed upstream in
    6.x): d/control, d/python-escript-doc.install, remove doc-base.4.
  * python-escript-doc reconciled to the 6.x doc set: ship the Sphinx Python
    API HTML; drop the removed cookbook/install PDFs and their doc-base entries.
  * Refer to the project as "esys-escript" in descriptions/doc-base titles.
  * Standards-Version 4.7.2.
  * d/patches refreshed against 6.1: 8 obsolete patches dropped (upstreamed or
    targeting removed files), 4 rebased; the full series applies cleanly to the
    6.1.1 source. (Build/lintian validation in a chroot still TODO — I don't
    have a Debian build host.)
  * d/watch fixed to match dotted upstream tags (6.1.1).

Follow-up (buildd, 6.1.1-1 in experimental): it built and installed on all
64-bit release arches (amd64, arm64, ppc64el, s390x) plus the 64-bit ports,
but FTBFS on armhf and i386. Both fail identically at the very first scons
configure step (conf.CheckCXX() -> "Cannot run C++ compiler 'g++'"); esys-escript
is 64-bit-only and never built on 32-bit (5.6-10 also failed there). Please
restrict Architecture to 64-bit for python3-escript and python3-escript-mpi,
e.g.:

  Architecture: amd64 arm64 ppc64el s390x riscv64 loong64 ppc64 sparc64

so buildd stops attempting the 32-bit ports.

Happy to push this to a branch on salsa for review if that's easier than the
patch. Thanks very much!

Lutz
