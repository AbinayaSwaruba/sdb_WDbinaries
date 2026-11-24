# MESA donor–accretor binary evolution

This project prepares a hot subdwarf donor (sdB) and a white dwarf (WD)
accretor in separate MESA runs, then uses the two saved stellar models to
initialize and evolve a binary. The repository also includes a modified copy
of MESA's binary mass-transfer source file, `binary_mdot.f90`.

## Repository layout

| Path | Purpose | Model files |
| --- | --- | --- |
| `star1_base/` | Builds the sdB donor. | Produces `sdb.mod`. |
| `star2_base/` | Builds the WD accretor. | Produces `wd.mod`. |
| `binary_base/` | Loads the saved stellar models and evolves the binary. | Reads `sdb.mod` and `wd.mod`. |
| `binary_mdot.f90` | Modified version of MESA's binary mass-transfer source file. | Fortran source used at build time. |

## Before running

Use a MESA installation compatible with the inlists and the modified
`binary_mdot.f90`. Verify the MESA release against the version used to develop
the patch. The release used for the binary calculations should be recorded
alongside the results.

Check each work directory's inlists and run scripts for site-specific paths.
Use the model-save and model-load paths in the active inlists when locating
or copying the `.mod` files. The commands below assume the standard MESA
work-directory scripts `mk` and `rn` are present. Use the repository's
corresponding build and run scripts if they have different names.

## 1. Build the donor model

From the repository root, work in `star1_base/`:

```bash
cd star1_base
./mk && ./rn
```

This stage produces `sdb.mod`, the saved donor model. Check the run's
termination output and verify that the saved model is from the intended
evolutionary stage before using it in the binary. If several donor models
are being prepared, keep separate copies and record the parameters used
for each one.

## 2. Build the accretor model

From the repository root, work in `star2_base/`:

```bash
cd star2_base
./mk && ./rn
```

This stage produces `wd.mod`, the saved white dwarf model. Check that the
saved model corresponds to the intended accretor mass and structure. The
donor and accretor preparation runs can be performed independently.

## 3. Build with the modified mass-transfer source

`binary_mdot.f90` is a modified copy of the upstream MESA file of the same
name. To reproduce a binary calculation that uses it, compare the copy
against the corresponding file in the intended MESA release and inspect
the patch before installing it. For example, once the upstream file has
been located, use `diff -u /path/to/upstream/binary_mdot.f90 binary_mdot.f90`
from the repository root.

Apply the change in a dedicated MESA checkout or installation. If the
binary work directory links against a prebuilt MESA installation, rebuild
the affected MESA code and then rebuild the binary work directory against
that installation. Check that the build uses the modified source and that
the executable was rebuilt after applying the patch. Record both the
upstream MESA release and the repository commit containing the patch.

## 4. Evolve the binary

Place the two saved models at the paths expected by the active binary
inlists. For example, if the preparation runs save the models directly in
their work directories and the binary inlists expect both files in
`binary_base/`, run the following from the repository root:

```bash
cp -i star1_base/sdb.mod binary_base/
cp -i star2_base/wd.mod binary_base/
```

Confirm that the binary inlists load `sdb.mod` for star 1 and `wd.mod` for
star 2. Check the initial orbital settings, mass-transfer controls, and
stopping criteria for the run you intend to perform. Then build and run:

```bash
cd binary_base
./mk && ./rn
```

At startup, verify the loaded component masses and the initial orbital
period or separation against the intended setup. During analysis, use the
configured history and profile outputs to inspect the evolution of the
donor and accretor masses, the orbital period or separation, and the
mass-transfer rate. Check the termination message to establish where the
calculation stopped.

## Reproducing a run

Record the repository commit, the MESA release, the inlists or parameter
changes used for all three stages, and the exact pair of saved models
passed to the binary run. Keep separate output locations for separate
systems or parameter choices so their histories and profiles can be
traced back to the corresponding inputs.

If a run fails to load a model, first verify the filename and path from
the binary run's working directory. If the models load but the expected
mass-transfer behavior is missing, check the binary controls and confirm
that the executable links against the MESA build containing the patch.
