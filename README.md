# AnalysisNotes

Each Mu2e analysis note here sits in its own folder, together with everything needed to redo it: the data it used, the exact software versions, a step-by-step processing recipe, the analysis code, and the results.

## Notes

| note | topic | data |
|---|---|---|
| [kpp-cosmics-2026-10](kpp-cosmics-2026-10/) | KPP cosmic run: timing, calorimeter and CRV positions, Michel decays (reproduces DocDB 58468) | runs 123680, 123681 |

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

## Common setup

```bash
source /cvmfs/mu2e.opensciencegrid.org/setupmu2e-art.sh                # mu2einit
cd <muse workdir> && muse setup                                         # reco and ntupling
source /cvmfs/mu2e.opensciencegrid.org/bin/pyenv.sh ana                 # uproot, awkward, scipy, matplotlib
```

Use `muse setup` and `pyenv ana` in separate shells. From a script, source `pyenv.sh` directly and not through a pipe, or the environment is lost.
