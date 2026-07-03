# Debian quilt patch refresh: 5.6 → 6.1.0

Assessment of the 16 active patches in the salsa `debian/patches/series`, tested
against a clean `git archive 6.1.0` export. Status legend:

- **DROP** — obsolete (fix upstreamed, or target file gone/unused).
- **REFRESH** — still applies, only context offset; `quilt push; quilt refresh`.
- **REWORK** — still wanted but no longer applies; needs rebasing (and a build).
- **DONE** — already rebased here.

| Patch | Verdict | Evidence / action |
|---|---|---|
| `32bitboostextract` | **DROP** | Patches `scons/templates/sid_options.py`, which in 6.1 is just `from scons.templates.debian_options import *` and is **not used** by the Debian build (it uses `debian/sid_options.py`). |
| `enable-silo.patch` | **DROP** | Same — silo is already configured in `debian/sid_options.py`. |
| `mpi-multiarch.patch` | **DROP** | Same — multiarch MPI is already handled in `debian/sid_options.py` (reads `mpi-default-dev/debian_defaults`). |
| `make.patch` | **DROP** | No top-level `Makefile` in 6.1 (pure scons). |
| `g++10-fix.patch` | **DROP** | Already upstream: 6.1 `paso/src/SystemMatrix_copyRemoteCoupleBlock.cpp:267` uses `new SparseMatrix<real_t>(...)`. |
| `boost_numpy.patch` | **DROP** | Upstream rewrote the boost-numpy probe (`site_scons/dependencies.py` has its own `have_boost_numpy` logic). Verify the upstream probe finds Debian's `libboost-numpy-dev`. |
| `fixmathjax` | **REFRESH** | Applies with offset to `doc/sphinx_api/conf.py`. |
| `openmpi-version.patch` | **REFRESH** | Applies with offset to `SConstruct`. |
| `exception.patch` | **REFRESH** | Applies with offset to `weipa/src/EscriptDataset.cpp`. Confirm the NetCDF4/g++ workaround is still required. |
| `ignore-flags.patch` | **REFRESH** | Applies with offset to `site_scons/extractdebbuild.py`. |
| `py3.13.patch` | **DONE** | Regenerated against 6.1 (raw-string fix, now **3** occurrences in `site_scons/site_init.py`); verified `patch -p1` applies cleanly. |
| `py3.patch` | **DROP** | Obsolete: 6.1 is py3-native. `extractdebbuild.py` already decodes; `call_python_config` was rewritten to decode internally and return str; the gmsh `which` block it touched was removed. |
| `which.patch` | **DONE** | Rebased to the four `run-escript.in` `which`->`command -v` hunks (verified). The `dependencies.py` hunks were dropped: that gmsh `which` block was removed upstream. |
| `fix-dpkg-buildflags-on-c-c++.patch` | **DONE** | Rebased against 6.1 `extractdebbuild.py` (the strict `mycflags!=mycxxflags` check is unchanged upstream); verified. |
| `tex.patch` | **DROP** | Obsolete in 6.1. The cookbook hunks' files are gone; in `doc/user/linearPDE.tex` all the broken `\Refe{}` refs it patched were rewritten away except `\Refe{TRILINOS}`, which now resolves (`\label{TRILINOS}` is defined in `doc/user/trilinos.tex`). The FTBFS it worked around is fixed upstream. |
| `use-t1-encoding.patch` | **DONE** | Reworked to the single `doc/user/user.tex` hunk (install/cookbook/inversion guides removed in 6.x); verified `patch -p1` applies cleanly. |

## Applied here

- Dropped the 8 obsolete patch files and removed them from `series`
  (32bitboostextract, enable-silo, mpi-multiarch, make, g++10-fix, boost_numpy,
  tex, py3).
- Rebased and verified `py3.13.patch`, `use-t1-encoding.patch`, `which.patch`,
  `fix-dpkg-buildflags-on-c-c++.patch` against the 6.1.0 source.
- Rewrote `series` (8 patches, grouped: refresh-able vs rebased).

## Remaining (needs the maintainer's build loop)

All patches are now settled. The **REFRESH** group (4: fixmathjax,
openmpi-version, exception, ignore-flags) applies with offset only (`quilt
refresh` to clean up line numbers). The reworked group (py3.13, use-t1-encoding,
which, fix-dpkg-buildflags) is rebased and verified. **A full-stack apply of the
whole series against the 6.1.0 source succeeds** (`quilt push -a` equivalent).
Dropped as obsolete: 32bitboostextract, enable-silo, mpi-multiarch, make,
g++10-fix, boost_numpy, tex, **py3**. Still TODO: a real `gbp buildpackage` /
`sbuild` + `lintian` to confirm the package builds and the docs/sphinx_api emit
where the .install expects.
