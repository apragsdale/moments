# Godambe / LRT calibration experiment: notes (updated 2026-11-10)

Where we stopped: the experiment script is written and a one-dataset pilot ran end to
end (2026-11-09). The full run will be on UW–Madison CHTC (HTCondor, Apptainer image
with moments from devel), one job per dataset (scenario + seed). The HTCondor files are
written (2026-11-10). Still open before running: the SFS M = 0 change (see "Next
steps").

## Goal

Check, with msprime simulations run through the full workflow, that

- confidence intervals from FIM, GIM (linear and log scale) and bootstrap refits have
  the right coverage, for LD statistics and for the SFS, on the same simulated data;
- the adjusted likelihood ratio test for migration has the right false positive rate
  under the null, and see how much power it has.

## Files in this directory (all untracked)

| Path | What it is |
| --- | --- |
| `godambe_calibration.py` | The experiment. `run --scenario S --seed N` simulates one dataset and analyzes it with LD and SFS; `simulate` / `analyze --method {ld,sfs}` do the stages separately; `summarize` pools all results. The docstring has usage. |
| `submit_replicates.sub` | HTCondor submit file: one job per row of the queue file (set to `queue_test.csv`; switch the `queue` line for the full run), 1 CPU, returns `data/{S}_{N}.pkl` and `results/{ld,sfs}_{S}_{N}.pkl`. |
| `run_replicate.sh` | Job script: sets one thread, runs `run --scenario $1 --seed $2 --regions $3 --boot-refits $4 --procs 1`. |
| `queue_inputs.csv` | Full run, `scenario, seed, regions, boot_refits`: mig 1–50, weak 1–50, null 1–100, all with 500 regions and 500 refits. |
| `queue_test.csv` | Cluster test: seeds 9001–9002 for each scenario, 40 regions, 20 refits (~1–2 h per job). Seeds differ from the full run, so outputs don't collide; delete them before `summarize`. |
| `NOTES.md` | This file. |

The semi-parametric prototypes (`calib*.py`) and the pilot outputs were deleted on
2026-11-10; their results are in the write-up and in the pilot section below.

Related, outside this directory:
- `../../notes_ld-godambe-review-and-calibration-2026-11-09.md`: the write-up of the review
  and the semi-parametric study (also at https://claude.ai/artifact/SGedFtxeod2MZASHBAsrzS).

## Experiment design (as in the script now)

- **Model:** IM. Ancestral Ne = 10,000; deme0 2,000; deme1 20,000; split 1,500
  generations ago; symmetric migration m. In genetic units (nu0, nu1, T, M, Ne) with
  M = 2 Ne m: (0.2, 2, 0.075, M, 10,000).
- **Scenarios / datasets:** `mig` M = 2 (m = 1e-4), 50 datasets; `weak` M = 0.2
  (m = 1e-5), 50 datasets; `null` M = 0, 100 datasets. Each dataset: 500 regions × 2 Mb,
  40 diploids per deme, u = 1.25e-8, r = 1e-8 (lowered from 1.5e-8 by the user).
- **Data per dataset:** each region simulated with msprime, written to VCF, parsed by
  `moments.LD.Parsing.compute_ld_statistics` (r bins as in the docs, 0 to 1e-3) and
  `moments.Spectrum.from_vcf` (polarized; msprime REF = ancestral). 500 bootstrap sets
  from seeded region resampling, shared by LD and SFS (LD sets built with
  `means_from_region_data`, as `get_bootstrap_sets` does internally; SFS sets are sums
  of the same regions' spectra).
- **Analysis per dataset and method:** fit nested (M = 0) and full model (full fit
  started from the nested optimum with M0 in {0.01, 0.1, 0.5} and from the truth; best
  kept), D = 2 ΔLL; LRT adjustments: none, current `LRT_adjust`, submatrix with centred
  J, inverse blocks with centred J (H, J, cU over all parameters, by finite differences:
  central, and second-order forward at M = 0, so M is never negative). For `mig` and
  `weak`: FIM, GIM, GIM (log), and 500 bootstrap refits (bootstrap SE and percentile
  CIs).
- **SFS:** projected to 40 haploids per deme for inference (full size 80 is much
  slower). Theta by multinomial scaling, reported as Ne = theta / (4 u L_total).
  SFS `GIM (log)` uses `moments.Godambe` with log=True, which still has the step-scaling
  problem (step = eps * log p); its intervals are for comparison only.

## Pilot (2026-11-09)

Environment: Python 3.13 venv with moments from branch `godambe_log_step`, msprime
1.4.5, scikit-allel 1.3.13, h5py, pandas, demes. (`bootstrap_data` needs at least as
many regions as statistics per bin, 15 for two populations, so the smallest pilot is 16
regions.)

| Stage | Pilot | Per unit | Full design (200 datasets × 500 regions) |
| --- | --- | --- | --- |
| simulate + parse | 16 regions: 224 s wall, 24 CPU-min on 8 cores | ~89 CPU-s per 2 Mb region | ~2,500 CPU-hours (~12 CPU-h per dataset) |
| analyze, LD | 53 s wall, 97 CPU-s (8 refits) | ~15 CPU-min per dataset with 100 refits (rough) | ~25 CPU-hours |
| analyze, SFS | 77 s wall, 159 CPU-s (8 refits) | ~25 CPU-min per dataset with 100 refits (rough) | ~40 CPU-hours |

Pilot timings: simulation and parsing dominate. Analysis timings are rough (fits and refits were not
timed separately). Disk: ~5.5 MB per dataset file (~1.1 GB for 200 datasets).

**Estimates for the current settings (2026-11-10)**, one dataset per job on one core,
500 bootstrap sets and 500 refits (rough: refit cost assumed equal to the pilot's
average optimization, ~7.5 CPU-s LD and ~12 CPU-s SFS; cluster cores may be slower).
The lower u and r should leave the parsing cost about the same: fewer SNPs, but wider
r bins in bp.

| Per dataset | `mig` / `weak` | `null` |
| --- | --- | --- |
| simulate + parse | ~12.5 h | ~12.5 h |
| analyze, LD | ~1 h | ~1 min (no refits) |
| analyze, SFS | ~1.7 h | ~1.5 min |

Totals: `mig` and `weak` ~760 CPU-h each, `null` ~1,250 CPU-h (almost all simulation),
~2,800 CPU-h overall. Requested wall time: ~24 h simulate, 2 h LD, 3–4 h SFS. Set
`OMP_NUM_THREADS=1` (and MKL/OpenBLAS equivalents). Disk: ~30 MB per dataset with 500
bootstrap sets (mostly SFS bootstraps), ~6 GB total.

The null datasets are what test the LRT false positive rate (the main open question);
100 gives a standard error of ~2.2 points at the 5% level. Option to get more from them:
evaluate the false positive rate on the first 100, 250 and all 500 regions.

## CHTC test run (2026-11-10)

6 jobs (2 per scenario, 40 regions, 20 refits, 1 CPU): all returned their files to
`data/` and `results/` (`preserve_relative_paths` works). CPU use ~0.98, memory ~1.1 GB,
disk ~620 MB (mostly the container image). Times on CHTC are ~2.5x the laptop pilot:

- simulate + parse: 189–243 s per region (mean 222), so ~31 h (up to 34 h) on one core
  for 500 regions;
- LD analysis: ~150 s without refits, ~20 s per refit (500 refits: ~2.8 h);
- SFS analysis: ~260 s without refits, ~39 s per refit (500 refits: ~5.5 h).

One core per job would be 31–39 h, so jobs now use `request_cpus = 4` (4 GB), passed to
`--procs`: ~10 h per mig/weak job, ~8 h per null job. Total ~7,000 CPU-hours; with 200
refits instead of 500, ~6,500.

## Finding from the pilot: SFS likelihood jumps at M = 0

`Spectrum.integrate` uses `Integration_nomig.integrate_nomig` when all migration rates
are exactly 0, a different integrator from the one used for any M > 0. On the pilot
SFS (n = 40 per deme) at the nested fit:

    ll(M = 0)     = -4393.7040
    ll(M = 1e-12) = -4393.7291   (jump of 0.025)
    ll(M = 1e-6)  = -4393.7287   (smooth above 0: slope ~444)

- Finite differences that start at M = 0 cross the jump. In the pilot, H_MM came out
  negative for steps up to 3e-3. Starting at M = 1e-12, H_MM = 1166.6 for every step
  from 1e-4 to 1e-2.
- Library impact: `moments.Godambe.LRT_adjust` (SFS) differentiates from M = 0 with
  eps = 0.01, which puts an error of about 0.025 / 1e-4 = 250 into H_MM ≈ 1167
  (~20%). D is also shifted by ~0.05 (nested fit uses the other integrator).
  Probably also affects Demes SFS inference with zero migration. Worth its own issue,
  or adding to the LRT issue.
- The LD likelihood has no such jump (H_MM stable from step 1e-4 to 1e-2).

## Next steps

1. **Script changes (proposed, not made):**
   - SFS: evaluate the no-migration model at M = 1e-12 instead of exactly 0 (in the
     nested fit and as the base point for the finite differences), so all SFS
     evaluations use the same integrator. Keep calling `LRT_adjust` at M = 0 as users
     do, so its results show the effect of the jump.
   - DONE 2026-11-10: seed per job. Each dataset is (scenario, seed); all random
     streams come from `SeedSequence([scenario index, seed, ...])` (regions:
     `[.., 0, region]`, bootstrap resampling: `[.., 1]`), so jobs are reproducible and
     never share seeds. Seeds are explicit in `queue_inputs.csv` rather than
     `$(Process)`, so resubmitting or adding datasets can't renumber them.
2. Decide whether to keep the SFS projection at 40 haploids per deme.
3. Cluster environment: the user is building an Apptainer image with moments from
   devel (needs the log=True step fix, commit 7bc7488 on `godambe_log_step`, merged
   first), msprime, scikit-allel, h5py, pandas, demes.
4. On CHTC: `mkdir -p log data results`, submit the test (`queue_test.csv`, 6 jobs)
   and check that all three files come back per job, plus run time and memory in the
   condor log. Then delete the test outputs, switch the `queue` line to
   `queue_inputs.csv`, and submit all 200. Jobs are ~15 h (mig/weak) and ~12.5 h (null) on one core;
   check that this fits CHTC's default job length, or split simulate and analyze into
   two jobs. Then `python godambe_calibration.py summarize` on the returned files.

## Related open items

- Issues to open or already drafted: J centring (also makes `score_stat` < 1 always);
  LRT too conservative (method fix: inverse blocks + centred J; Varin, Reid & Firth
  2011 §2.3; Chandler & Bate 2007 as an alternative); log=True step size (fixed on
  `godambe_log_step`).
- To come back to: unit tests that `moments.Demes.Inference` (optimize/optimize_LD,
  uncerts/uncerts_LD) matches `Inference`/`Godambe` for split_mig after converting
  physical <-> genetic units (fits to ~1e-4, SEs via the delta method to ~1e-3,
  log=False and log=True). Prototype first to measure the tolerances.
- h2py (`nwcol/h2py`, `h2py/uncerts.py`): `compute_lrt_adjustment` uses `inv(HH)` with
  HH the Hessian of +LL, so the factor is negative (fix: `inv(-HH)`); global model cache
  keyed only by parameters. The user is telling the author.
