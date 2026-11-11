# Godambe / LRT calibration experiment: notes (updated 2026-11-11)

Where we stopped: the first full run on UW–Madison CHTC (193 of 200 datasets) is
summarized in "First full run: preliminary results". It is being redone from scratch
with infinite-sites mutations and chained optimizers (see "Second run"); the write-up
waits for that. The proposed fix for the LRT is described in "The LRT adjustment: method
and intuition".

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
| `submit_replicates.sub` | HTCondor submit file: one job per row of the queue file (switch the `queue` line between `queue_test.csv` and `queue_inputs.csv`), 4 CPUs, returns `data/{S}_{N}.pkl` and `results/{ld,sfs}_{S}_{N}.pkl`. |
| `run_replicate.sh` | Job script: sets one thread, runs `run --scenario $1 --seed $2 --regions $3 --boot-refits $4 --procs $5` (`$5` = `request_cpus`). |
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

One core per job would be 31–39 h, so jobs now use `request_cpus = 4`, passed to
`--procs`: ~10 h per mig/weak job, ~8 h per null job. Total ~7,000 CPU-hours; with 200
refits instead of 500, ~6,500.

Memory: 4 GB was too little for 4 processes. Full-run jobs went over the cgroup limit,
because each worker process simulating and parsing a region needs about as much memory
as the whole 1-CPU test job (~1 GB), so the peak grows with the number of processes.
Request more than 1 GB per process plus the parent.

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

## The LRT adjustment: method and intuition

This is the proposed fix for `moments.Godambe.LRT_adjust` and
`moments.LD.Godambe.LRT_adjust`. In the script it is `adjustment_inverse_blocks` with
`centre=True`.

### Notation

- θ = (ψ, τ): ψ are the q parameters fixed by the nested model (here ψ = M, q = 1, fixed
  at 0), and τ are the other parameters (nu0, nu1, T, Ne), which are refit in both
  models.
- cℓ(θ): the composite log-likelihood. The LD and SFS likelihoods both treat linked
  data as independent (LD: all bins and regions; SFS: all sites).
- U = ∇cℓ: the score.
- H = −∇²cℓ: the sensitivity matrix, the curvature of the likelihood.
- J = Var(U): the variability matrix, the variance of the score across datasets.
  Estimated from the scores of the bootstrap datasets, evaluated at the parameters of
  the nested fit.
- D = 2 [cℓ(θ̂) − cℓ(0, τ̂₀)]: the likelihood ratio statistic. θ̂ is the full fit and τ̂₀
  is the nested fit.

All of H, J and the mean score cU are computed **over all parameters**, at the nested
fit.

### Formula

    adj = q / tr[ (H^ψψ)⁻¹ G^ψψ ],    with  H^ψψ = (H⁻¹)_ψψ,  G^ψψ = (H⁻¹ J H⁻¹)_ψψ,
                                       and  J = mean(U Uᵀ) − cU cUᵀ  (centred)

Compare adj · D with the null distribution: ½χ²₀ + ½χ²₁ here, since M = 0 is on the
boundary (Self & Liang 1987); χ²_q in general. With q = 1 this is just
adj = H^MM / G^MM.

Reference: Varin, Reid & Firth (2011) §2.3. D tends to Σ λⱼ Zⱼ², where the λⱼ are the
eigenvalues of (H^ψψ)⁻¹ G^ψψ. Dividing by the mean eigenvalue is first-order moment
matching (Rotnitzky & Jewell 1990; Geys et al. 1999). For q = 1 it is exact
asymptotically.

### Intuition

1. **Why an adjustment is needed.** For a true likelihood, H = J (the information
   identity), and D ~ χ²_q (Wilks). A composite likelihood of linked data counts each
   correlated piece of data as independent evidence. The likelihood surface is then too
   sharply peaked for how much the estimates actually vary between datasets: J is
   "larger" than H. D inherits that overconfidence and is inflated by a factor λ. The
   adjustment divides it back out.

2. **What λ is.** For q = 1, D ≈ (ψ̂ − ψ₀)² / H^ψψ. The numerator varies across datasets
   with the true variance of the estimate, G^ψψ: the GIM (sandwich) variance. The
   denominator is the variance the likelihood's curvature claims, H^ψψ: the FIM-type
   variance. So

       λ = G^ψψ / H^ψψ = (SD_GIM(M) / SD_FIM(M))²,    adj = 1 / λ.

   (H and J here are at the same parameters, the nested fit; the CIs in the results
   use the full fit, so the match there is only approximate.)
   **The LRT adjustment is the same correction as using GIM instead of FIM for the
   uncertainty of M, applied to the test.** Median adj under the null is ~0.21 for LD
   and ~0.11 for the SFS, so the unadjusted D is 5–10 times too large. The FIM interval
   undercoverage in the results is the same effect.

3. **Why inverse blocks rather than submatrices.** The current `LRT_adjust` uses
   q / tr(J_ψψ H_ψψ⁻¹), the ψψ submatrices of H and J. That measures how inflated the
   score for M is **with the other parameters held fixed**. But D is a profile
   statistic: when M moves, T and the sizes are refit and absorb part of the signal (T
   and M have correlation −0.89 in H). The relevant score is the part of the M score
   the other parameters can't absorb, the efficient score
   U* = U_ψ − H_ψτ H_ττ⁻¹ U_τ. Its curvature is the Schur complement
   (H^ψψ)⁻¹ = H_ψψ − H_ψτ H_ττ⁻¹ H_τψ. Then λ = Var(U*) / (H^ψψ)⁻¹, which is exactly the
   inverse-block formula. (The ψ rows of H⁻¹ are H^ψψ [I, −H_ψτ H_ττ⁻¹], so
   G^ψψ = H^ψψ Var(U*) H^ψψ.) The submatrix version is right only when H_ψτ = 0, that
   is, when M is uncorrelated with the other parameters.

4. **Why J must be centred.** J should be the variance of the score: its noise. At the
   nested fit the model is not stationary in M, and the M score has a nonzero mean.
   That mean is the evidence for migration that the test is trying to detect. The
   uncentred mean(U Uᵀ) = Var(U) + cU cUᵀ counts that signal as noise. adj shrinks as the
   evidence grows, so the test loses power and is conservative even under the null,
   because at the boundary about half the null datasets have a positive M slope. It is
   also why `score_stat` with an uncentred J is always below 1 (Sherman–Morrison:
   s / (1 + s)).

5. **Neither change works alone.** In both the semi-parametric study and this run,
   centring J with submatrices, or inverse blocks with an uncentred J, stays
   conservative. Both are needed.

### Numerical details that matter

- Finite differences at M = 0: use second-order forward differences in M, never stepping
  to negative M. They match central differences to ~1e-6 for LD.
- SFS: evaluate "M = 0" at M = 1e-12 (`M_FLOOR`), so every evaluation uses the same
  integrator. See "SFS likelihood jumps at M = 0" above.
- Start the full fit from the nested optimum (plus a few M starts), or D comes out too
  small.

### Known weakness: H can be indefinite at the nested fit

In LD `mig` (M = 2), the inverse-block adj came out **negative in 30 of 46 datasets**.
The full-parameter H at the nested fit has one negative eigenvalue. Its eigenvector
trades nu1 against M (about (0, 0.70, 0, −0.71, 0)). When the nested model is badly
wrong, its fit can sit at a saddle of the full likelihood, and H⁻¹ then has the wrong
sign. The submatrix formula only needs H_MM > 0, so it doesn't see this.

- It happened only under strong alternatives, where D ≈ 3,000. It never happened under
  the null or with weak migration (0 of 147 datasets).
- A negative adj · D gives p = 1, so power in LD `mig` at α = 0.05 is 0.30 instead of
  1.00.
- Possible guards (not tested):
  - check that H is positive definite and warn or raise if not;
  - use H (and J) at the full fit, which is asymptotically equivalent under the null
    and local alternatives;
  - treat adj ≤ 0 as "reject" when D is large.

## First full run: preliminary results (2026-11-11)

Copied from CHTC to `results/` (test seeds 9001–9002 excluded). There are 46 `mig`,
47 `weak` and 100 `null` datasets for each method; 7 jobs were still running. Shown by
`python godambe_calibration.py summarize`. Rerun everything once the last jobs are in.

**LRT under the null.** Fraction of p-values below α = .01 / .05 / .10 / .25 (±~0.02
at 0.05):

| Adjustment | LD | SFS |
| --- | --- | --- |
| None | .14 / .21 / .23 / .38 | .25 / .32 / .36 / .40 |
| `LRT_adjust` (library) | .00 / .00 / .03 / .17 | .00 / .01 / .05 / .15 |
| Submatrix, centred J | .00 / .00 / .04 / .17 | .01 / .03 / .06 / .15 |
| Inverse blocks, centred J | .00 / .04 / .09 / .22 | .01 / .06 / .07 / .27 |

About 57% of null datasets have D = 0.

**Power:**

- Weak migration (M = 0.2), α = 0.01:
  - LD: inverse blocks 0.96, `LRT_adjust` 0.51.
  - SFS: `LRT_adjust` 0.00, inverse blocks 1.00.
- Strong migration (M = 2), SFS: `LRT_adjust` shrinks D ≈ 25,000 to adj · D ≈ 5, so
  p ≈ 0.011. It rejects at 0.05 but never at 0.01.
- Strong migration (M = 2), LD: inverse blocks only 0.30, because of the negative
  adjustments described above.

**95% CI coverage** (n ≈ 46; ±~0.03):

- FIM undercovers: 0.70–0.85 LD `mig`, 0.40–0.72 LD `weak`, 0.38–0.66 SFS.
- GIM covers 0.85–1.00 and is slightly conservative: its mean SD is up to ~35% above
  the spread of the estimates (SFS T: 0.0156 against 0.0115 of the truth).
- GIM (log) is similar; for the SFS it is 0.98–1.00.

**Issues to look into (hypotheses, not checked):**

- **SFS fits are slightly biased:** −0.2% to −0.5% relative. That is about 1 SD for Ne
  and nu0, which lowers SFS coverage for those (GIM: Ne 0.85 `mig`, 0.91 `weak`). Guess:
  msprime's finite-sites mutations against moments' infinite-sites model; recurrent
  hits would lower diversity and pull Ne down. Check: count multi-hit sites in a
  simulated region.
- **LD `weak` bootstrap is too narrow:** the bootstrap SD of M is 0.155 of the truth,
  against an empirical SD of 0.27 and a GIM SD of 0.24. The other parameters are ~18%
  too narrow. Bootstrap-SE coverage of M is 0.62. LD `mig` and SFS bootstraps are fine.
  Guess: the L-BFGS-B refits start from `p_full` and stop early where the surface is
  flat near small M. Check: refit some bootstrap sets from perturbed starts.
- **Negative inverse-block adjustment under strong migration** (above).

## Second run (prepared 2026-11-11)

The first run's results above are superseded: all datasets are being rerun from
scratch with these changes.

- **Infinite sites.** `sim_mutations(..., discrete_genome=False)`, as in the moments
  models, so there are no recurrent mutations. This may explain the SFS bias in the
  first run. VCFs are written with `position_transform="legacy"`, because ~9 sites per
  2 Mb region would otherwise share a rounded position.
- **Chained optimizers.** Every fit (nested, each start of the full fit, and each
  bootstrap refit) runs `FIT_CHAIN = ["lbfgsb", "fmin", "powell"]` on the log scale, each
  starting from the previous result. The first run used a single L-BFGS-B for the full
  fits and refits, and two rounds of fmin for the nested fit. The log-likelihood after
  each optimizer is saved (`fit_traces`, `boot_traces`) and printed in the job log.
- **200 bootstrap refits** instead of 500 (`queue_inputs.csv`). All 500 bootstrap sets
  are still used for the Godambe methods.
- **`request_memory = 8GB`** for 4 CPUs.

Smoke test (laptop, one `weak` dataset, 16 regions, 8 refits, 8 processes):

- **L-BFGS-B alone stops well short, especially for LD.** The fmin step after it
  improved the LD full fits by 16–336 log-likelihood units, and the LD refits by a
  median of 7 (max 72). For the SFS: the full fits by 0.05–26, the refits by about 0.1.
  The nested fits gained 18 (LD) and 0.1 (SFS). Powell after fmin changed nothing
  (< 1e-5). The first run's LD fits were probably under-converged too, but its data had
  500 regions, so the size of the effect there is unknown.
- **Time:** simulation was as before (206 s wall for 16 regions). Analysis took about
  twice as long as the single-optimizer pilot: LD 97 s, SFS 139 s, against 53 s and 77 s.
  Estimate on CHTC: ~40 s per LD refit and ~75 s per SFS refit, so 200 refits is
  ~6.5 CPU-h per mig/weak dataset, a bit less than the first run's 500 single-optimizer
  refits (~8.3 CPU-h).

Before submitting on CHTC: move the first run's `data/`, `results/` and `log/` aside,
for example to `run1_data/` etc. `simulate` and `analyze` skip outputs that already
exist, so old files would otherwise be kept and nothing would be rerun.

## Next steps

1. Pull the branch on CHTC, move the first run's outputs aside, and submit the second
   run (`queue_inputs.csv`). Then copy `results/` and run `summarize`.
2. Optionally check the two hypotheses above (SFS bias; LD bootstrap refits).
3. Library fix: inverse blocks with centred J in both `LRT_adjust` functions, with a
   guard for an indefinite H, and M = 0 handled without crossing the SFS integrator
   jump. Open the issues first (texts drafted).
4. Done 2026-11-10: seed per job; the SFS `M_FLOOR`; the CHTC test run. The SFS
   projection stayed at 40 haploids per deme.

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
