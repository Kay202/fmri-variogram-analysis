# Replication package

Code for *What a Variogram on a Voxel Lattice Can and Cannot Determine:
Non-Parametric Estimation, an Interpolation Signature, and the Limits of
Cluster-Extent Calibration*.

## What is distributed

The code, and nothing else. The GitHub repository and the Zenodo archive hold
the same thing: every script behind every table and figure in the paper, and
`MANIFEST.csv` naming the script and the inputs behind each one.

No data is distributed, in either place. `data/` is where the pipeline writes
its intermediates, and it is produced by the stages below rather than shipped
with them. A fresh download therefore rebuilds nothing on its own: the chain has
to be run from the public dataset, and its first stage needs MATLAB with SPM,
about 132 GB of disk, and the stimulus-presentation logs, which are not in the
OpenNeuro release and have to be requested from the authors of the source study.

What this package is for is reading the method and re-running it, not a
one-command rebuild of the figures.

## Running the output stage

Once `data/` exists, every table and figure in the manuscript comes from:

```
Rscript RunOutputs.R
```

Sixteen LaTeX table fragments land in `outputs/tables`, seventeen figures in
`outputs/figures`. It takes a few seconds and needs R only. The script counts
what appeared and names anything missing.

## Run order

Two drivers do all of this for you: `RunAnalysis.R` for steps 4 to 6, then
`RunOutputs.R` for steps 7 and 8. Run the scripts directly only if you want to
run one stage on its own.

Steps 1 to 3 are the imaging stage, and since no data is distributed they are
where everything starts. They need MATLAB with SPM, the fMRIPrep derivatives,
the stimulus-presentation logs, and about 132 GB of free disk; see "What
reproduces" below. Everything from step 4 on runs off the `data/` they produce.

```
1. codes/GLM.m                first level, both sessions; writes the residual store
2. codes/Variogram.m          the empirical variogram
   codes/VariogramDepth.m     variograms by depth shell
   codes/DepthProfile.m       the depth profile
3. codes/GeneratorCheck.m     optional; produces no manuscript table

4. codes/Simulation.R run     the simulation study; hours, checkpointed
5. codes/Analysis.R           shell models, interpolation signature, depth
                              gradient, reference covariance, identifiability
                              sweep
6. codes/Inference.m          the inference experiment, four stages

7. codes/Tables.R             every LaTeX table fragment
8. codes/Figures.R            every figure
```

The three scripts in step 2 read the residual store and are independent of each
other, so they can run in any order once step 1 is done.

Step 5 must come before step 6: the reference covariance is fitted around the
shell models, and the inference experiment reads both. Steps 7 and 8 can run in
either order.

Step 4 takes a mode argument and **must be given `run`**; with no argument it
performs a timing benchmark and writes nothing.

Step 6 resumes by itself. Each of its four stages writes one CSV and is skipped
when that CSV exists, so an interrupted run continues where it stopped. Delete a
CSV to force that stage to re-run.

## Layout

```
codes/                  the scripts you run
  helpers/              sourced or added to the path, not run directly
data/                   not distributed; written by the stages above
  out_matlab/           the per-run lag files the analysis starts from
  out_R/                one fit the analysis reads
  intermediates/        what Tables.R and Figures.R read
outputs/
  tables/  figures/     written at run time
RunAnalysis.R           driver for steps 4 to 6
RunOutputs.R            driver for steps 7 and 8
MANIFEST.csv            artefact -> script -> inputs
```

No script needs a path edited. Each one resolves its own location and finds
`data/` and `outputs/` from there. The one exception is the imaging stage, which
needs three locations on your machine; see "Configuring paths".

## What reproduces, and what does not

| stage | from | to | needs | time |
|---|---|---|---|---|
| imaging | OpenNeuro dataset | per-run lag files | MATLAB, SPM, the data | days |
| analysis | lag files | intermediate tables | R, and MATLAB for one stage | hours |
| outputs | intermediate tables | every table and figure | R | seconds |

The imaging stage cannot be included. The first-level residuals alone are 263 GB
across the 21 participants, and rebuilding them needs SPM, the imaging data, and
the stimulus-presentation logs, which are not in the OpenNeuro release. Its four
scripts are in `codes/` so the chain is legible, and `data/out_matlab` is where
they write. The dataset is at <https://openneuro.org/datasets/ds005891>.

The later two stages run as soon as `data/` exists, the output stage in seconds.
Neither runs from a fresh download, because `data/` is not distributed.

Six tables in the manuscript carry their numbers inline rather than being
generated: the Gaussian fit diagnostics, the octant table, the method
comparison, the first-level sensitivity table, and the two finite-sample tables.
`RunOutputs.R` does not rebuild them. The CSVs behind them are written to
`data/intermediates`.

## Configuring paths

The analysis and output stages need no configuration.

The imaging stage needs three locations that depend on your machine. Set them as
environment variables, or edit the defaults at the top of
`codes/helpers/vg_paths.m`:

| variable | what |
|---|---|
| `VG_SPM` | SPM installation (SPM 25.01.02 was used) |
| `VG_FMRIPREP` | fMRIPrep derivatives, the directory holding `sub-*/` |
| `VG_LOGS` | Presentation logs, the directory holding `Sub*/ExpInfo_run*.mat` |

`VG_SPMFIRST` is needed only by scripts that read the original first-level model
rather than the Model B residuals. `VG_ROOT` overrides the project root.

Anything not set, or set to a directory that does not exist, raises an error
naming the variable rather than failing later on a missing file.

## Requirements

R 4.5.1, with `data.table`, `splines`, `nnls`, `minpack.lm` and `readxl`.

MATLAB R2025a, with the Image Processing Toolbox for `bwconncomp`.

If SPM is on the saved MATLAB path it shadows `isfile`, which `readtable` calls
internally, and the failure surfaces as exit 255 with no catchable error. Every
MATLAB script therefore begins with `restoredefaultpath; rehash toolboxcache;`.
If SPM is on your saved path you will still see two startup warnings that
`isfile` and `isfolder` have the same names as built-ins; they are emitted before
the first line runs, and the `restoredefaultpath` clears the shadowing.

## License

Copyright (c) 2026 Godfred Arhin and Emmanuel Kubuafor.

Creative Commons Attribution 4.0 International. The full terms are in `LICENSE`,
which is the legal code verbatim so that GitHub and Zenodo both read it as
CC BY 4.0, and the same license is recorded on the archive minted from this
repository. Use the code as you like, with credit.
