# AnalysisNotes

Each Mu2e analysis note here sits in its own folder, together with everything needed to redo it: the data it used, the exact software versions, a step-by-step processing recipe, the analysis code, and the results.

## Notes

| note | topic | data |
|---|---|---|
| [kpp-cosmics-2026-10](kpp-cosmics-2026-10/) | KPP cosmic run: timing, calorimeter and CRV positions, Michel decays (reproduces DocDB 58468) | runs 123680, 123681 |
| [kpp-trk-2026-10](kpp-trk-2026-10/) | KPP `kpp_trk` runs: processing with PassN (three stopgaps, tracker addressing check), tracker-CRV timing and positions | runs 124984, 124986, 124989 |
| [crv-gap-2026-10](crv-gap-2026-10/) | CRV module-gap inefficiency from tracks extrapolated to the EX sector, data vs MC (reproduces DocDB 57978) | run 124155; MDC2025au extracted-cosmic MC |

Run 1A and Run 1B notes will follow.

## Layout of a note

```
<topic>-<yyyy-mm>/
  README.md     the note: goal, data, software, recipe, results, caveats
  scripts/      processing scripts (reco, ntupling); take the muse workdir as an argument
  fcl/          job configurations the recipe adds on top of the release
  *.py          analysis code
  plots/        figures and summary.json with the fitted numbers
```

Conventions:
- **Pin everything.** Record the git sha of every repository in the muse workdir, the dataset names, and the run numbers.
- **No hard-coded personal paths.** Scripts take the workdir, the inputs and the output directory as arguments.
- **No ROOT or art files in the repo** (see `.gitignore`). Write down where the outputs live instead.
- **Numbers in a README come from `summary.json`**, which the analysis code writes, not from copying values by hand.

## Where to work

Do not work in your home area (`$HOME`, nashome). It is small and shared, and a muse build alone does not fit there.

| what | where | typical size |
|---|---|---|
| muse workdir (clones + build) | `/exp/mu2e/app/users/$USER/...` | about 4 GB per workdir |
| outputs (art, ntuples, logs) | `/exp/mu2e/data/users/$USER/...` | 1-2 GB per run processed |
| plots, summary.json | the note folder in this repo | under 1 MB |

`/tmp` is not for outputs either. Raw files are read in place from `/pnfs` and never copied.

## Common setup

```bash
source /cvmfs/mu2e.opensciencegrid.org/setupmu2e-art.sh                # mu2einit
cd <muse workdir> && muse setup                                         # reco and ntupling
source /cvmfs/mu2e.opensciencegrid.org/bin/pyenv.sh ana                 # uproot, awkward, scipy, matplotlib
```

Use `muse setup` and `pyenv ana` in separate shells. From a script, source `pyenv.sh` directly and not through a pipe, or the environment is lost.
