"""
Calibration of the moments uncertainty and likelihood ratio test methods, for LD
statistics and for the SFS, using msprime simulations.

For many independent simulated datasets from a two-population isolation-with-migration
model with known parameters, we run the full workflows: simulate regions with msprime
and write them to VCFs, then

- LD: compute LD statistics with moments.LD.Parsing.compute_ld_statistics, and fit
  moments.LD.Demographics2D.split_mig to the means and varcovs over regions;
- SFS: compute the joint SFS of each region with moments.Spectrum.from_vcf, and fit
  moments.Demographics2D.split_mig to their sum (with theta fit by multinomial scaling).

Both analyses use the same simulated data and the same bootstrap replicates (the same
regions resampled with replacement). We then check:

1. Confidence intervals. For datasets with migration, we compute 95% CIs for every
   parameter (nu0, nu1, T, M = 2 Ne m, and Ne) using
     - the Fisher information matrix (FIM_uncert),
     - the Godambe information matrix (GIM_uncert), on the linear and log scales, and
     - refits of the model to bootstrap replicates (percentile intervals and bootstrap
       standard errors),
   and report how often each interval contains the true value (should be about 95%),
   and the average interval widths.

2. Likelihood ratio tests for migration. For each dataset, we fit the nested model
   without migration and the full model, and compute D = 2 * (ll_full - ll_nested).
   Without migration ("null" scenario), the adjusted p-values should be approximately
   uniform (with a point mass at one, as the migration rate is on the boundary), so a
   fraction alpha of them are below alpha. We compare
     - no adjustment,
     - LRT_adjust (moments.LD.Godambe or moments.Godambe): first-order moment matching
       using the nested parameter's submatrices of H and an uncentred J, following
       Coffman et al. (2016),
     - the same with J centred,
     - first-order moment matching using the blocks of H^-1 and H^-1 J H^-1 for the
       nested parameter, with H and J over all parameters and J centred (the composite
       likelihood result; Varin, Reid and Firth 2011, Section 2.3).
   For the last two, H and J are computed here by finite differences: central
   differences, except for the migration rate at the nested fit (zero), which uses
   second-order forward differences, so that the model is never evaluated at a negative
   migration rate. The "weak" and "mig" scenarios measure the power of the tests.

Scenarios, and the number of independent datasets we suggest for each:
    mig   M = 2 (m = 1e-4)    50 datasets: CI coverage, and LRT power
    weak  M = 0.2 (m = 1e-5)  50 datasets: CI coverage, and LRT power
    null  M = 0               100 datasets: LRT false positive rate

Each dataset is identified by its scenario and a seed, from which all of its random
seeds are derived. Datasets are saved in ./data/{scenario}_{seed}.pkl, and the analysis
of each in ./results/{method}_{scenario}_{seed}.pkl. Stages skip outputs that already
exist, so they can be stopped and rerun.

    # one dataset: simulate and parse, then analyze with both methods
    python godambe_calibration.py run --scenario mig --seed 1

    # or the same in separate stages
    python godambe_calibration.py simulate --scenario mig --seed 1
    python godambe_calibration.py analyze --scenario mig --seed 1 --method ld
    python godambe_calibration.py analyze --scenario mig --seed 1 --method sfs

    # a quick pilot, to check the run time of each stage
    python godambe_calibration.py run --scenario mig --seed 1 --regions 16 \\
        --boot-refits 8 --procs 8

    # after all datasets are done, print CI coverage and LRT calibration
    python godambe_calibration.py summarize

On a cluster, run one job per dataset (seeds 1-50 for mig and weak, 1-100 for null).
submit_replicates.sub, run_replicate.sh and queue_inputs.csv do this on HTCondor (UW-
Madison CHTC); see NOTES.md.

--procs sets the number of processes (default 1). --regions (default REGIONS) needs at
least 16 regions, as bootstrap_data needs at least as many regions as statistics per
bin.

Requires msprime, demes, scikit-allel, h5py, and pandas, in addition to moments. The
log-scale LD GIM intervals need the fix for log=True step sizes in moments.LD.Godambe.
"""

import argparse
import contextlib
import glob
import io
import os
import sys
import pickle
import tempfile
import time
import warnings
from multiprocessing import Pool

import demes
import msprime
import numpy as np

import moments
import moments.LD

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------

# The demographic model: an ancestral population of size NE splits T_GEN generations
# ago into deme0 (size N0) and deme1 (size N1), which exchange migrants symmetrically
# at a per-generation rate m that depends on the scenario.
NE = 10000
N0 = 2000
N1 = 20000
T_GEN = 1500
SCENARIOS = {"mig": 1e-4, "weak": 1e-5, "null": 0.0}
POPS = ["deme0", "deme1"]

# Data: each dataset has REGIONS independent regions of length L, with NUM_SAMPLES
# diploid individuals sampled from each deme.
REGIONS = 500
L = 2e6
NUM_SAMPLES = 40
U = 1.25e-8
R = 1.00e-8

# Bootstrap replicates (regions resampled with replacement), shared by the LD and SFS
# analyses. All are used for the Godambe methods; the model is refit to the first
# --boot-refits of them.
NUM_BOOTSTRAPS = 500

# LD statistics: recombination bin edges, in Morgans.
R_BINS = np.concatenate(([0], np.logspace(-6, -3, 13)))

# SFS: the haploid sample size per deme that the SFS is projected to for inference.
# Integrating the 2D SFS at the full size (2 * NUM_SAMPLES = 80) is much slower.
SFS_SAMPLE_SIZE = 40

# Parameters (nu0, nu1, T, M, Ne) in genetic units, with M = 2 Ne m. In the LD model,
# Ne is the last parameter. In the SFS model, Ne is computed from
# theta = 4 Ne U L_total, where L_total is the total length of the regions.
PARAM_NAMES = ["nu0", "nu1", "T", "M", "Ne"]
MIG_IDX = 3
LOWER = [1e-3, 1e-3, 1e-3, 0, 100]
UPPER = [100, 100, 10, 100, 1e6]
DATA_DIR = "data"
RESULTS_DIR = "results"


def true_params(scenario):
    m = SCENARIOS[scenario]
    return np.array([N0 / NE, N1 / NE, T_GEN / 2 / NE, 2 * NE * m, NE])


def demographic_model(m):
    b = demes.Builder(time_units="generations")
    b.add_deme("anc", epochs=[dict(start_size=NE, end_time=T_GEN)])
    b.add_deme("deme0", ancestors=["anc"], epochs=[dict(start_size=N0)])
    b.add_deme("deme1", ancestors=["anc"], epochs=[dict(start_size=N1)])
    if m > 0:
        b.add_migration(demes=["deme0", "deme1"], rate=m)
    return b.resolve()


def sample_ids():
    # msprime creates the sampled individuals in the order of the samples dictionary,
    # and we name them ind_0, ind_1, ... in the VCF.
    return {
        pop: [f"ind_{j * NUM_SAMPLES + i}" for i in range(NUM_SAMPLES)]
        for j, pop in enumerate(POPS)
    }


# -----------------------------------------------------------------------------
# Simulating and parsing data
# -----------------------------------------------------------------------------


def _seed_sequence(scenario, seed, *keys):
    # Independent random streams for each scenario, dataset seed, and use within it.
    return np.random.SeedSequence([list(SCENARIOS).index(scenario), seed, *keys])


def _region_seeds(scenario, seed, region):
    # msprime ancestry and mutation seeds for one region
    ss = _seed_sequence(scenario, seed, 0, region)
    return [int(s) % (2**31 - 2) + 1 for s in ss.generate_state(2)]


def _map(func, args, procs):
    # map with a process pool, or in this process with procs = 1
    if procs == 1:
        return list(map(func, args))
    with Pool(procs) as pool:
        return pool.map(func, args)


def data_path(scenario, seed):
    return os.path.join(DATA_DIR, f"{scenario}_{seed}.pkl")


def results_path(method, scenario, seed):
    return os.path.join(RESULTS_DIR, f"{method}_{scenario}_{seed}.pkl")


def simulate_region(args):
    """
    Simulates one region and writes it to a VCF. Returns its LD statistics, from
    moments.LD.Parsing.compute_ld_statistics, and its joint SFS (at the full sample
    size), from moments.Spectrum.from_vcf. msprime writes the ancestral allele as REF,
    so the SFS is polarized.
    """
    scenario, seed, region, tmpdir = args
    anc_seed, mut_seed = _region_seeds(scenario, seed, region)
    demog = msprime.Demography.from_demes(demographic_model(SCENARIOS[scenario]))
    ts = msprime.sim_ancestry(
        {pop: NUM_SAMPLES for pop in POPS},
        demography=demog,
        sequence_length=L,
        recombination_rate=R,
        random_seed=anc_seed,
    )
    ts = msprime.sim_mutations(ts, rate=U, random_seed=mut_seed)

    vcf_file = os.path.join(tmpdir, f"region_{region}.vcf")
    names = [f"ind_{i}" for i in range(ts.num_individuals)]
    with open(vcf_file, "w") as fout:
        ts.write_vcf(fout, individual_names=names, allow_position_zero=True)
    quiet = io.StringIO()
    with (
        warnings.catch_warnings(),
        contextlib.redirect_stdout(quiet),
        contextlib.redirect_stderr(quiet),
    ):
        warnings.simplefilter("ignore")
        ld_stats = moments.LD.Parsing.compute_ld_statistics(
            vcf_file,
            rec_map_file=os.path.join(tmpdir, "rec_map.txt"),
            pop_file=os.path.join(tmpdir, "pops.txt"),
            pops=POPS,
            r_bins=R_BINS,
            use_h5=False,
            report=False,
        )
        sfs = moments.Spectrum.from_vcf(vcf_file, pop_mapping=sample_ids())
    os.remove(vcf_file)
    return ld_stats, np.asarray(sfs.data)


def write_pop_file_and_map(tmpdir):
    with open(os.path.join(tmpdir, "pops.txt"), "w") as fout:
        fout.write("sample\tpop\n")
        for pop, ids in sample_ids().items():
            for sample in ids:
                fout.write(f"{sample}\t{pop}\n")
    # A uniform recombination map, in cM.
    with open(os.path.join(tmpdir, "rec_map.txt"), "w") as fout:
        fout.write("pos\tMap(cM)\n")
        fout.write("0\t0\n")
        fout.write(f"{int(L)}\t{R * L * 100}\n")


def simulate(scenario, seed, regions, procs):
    """Simulates and parses one dataset, and saves it. Returns its path."""
    out = data_path(scenario, seed)
    if os.path.exists(out):
        print(f"{out}: already exists", flush=True)
        return out
    os.makedirs(DATA_DIR, exist_ok=True)
    t0 = time.time()
    with tempfile.TemporaryDirectory() as tmpdir:
        write_pop_file_and_map(tmpdir)
        args = [(scenario, seed, i, tmpdir) for i in range(regions)]
        outputs = _map(simulate_region, args, procs)
    region_ld = {i: o[0] for i, o in enumerate(outputs)}
    region_sfs = np.array([o[1] for o in outputs])

    # The same bootstrap replicates for LD and SFS: each resamples the regions with
    # replacement, as moments.LD.Parsing.get_bootstrap_sets does.
    rng = np.random.default_rng(_seed_sequence(scenario, seed, 1))
    boot_idx = rng.integers(regions, size=(NUM_BOOTSTRAPS, regions))
    stats = region_ld[0]["stats"]
    ld_boots = [
        moments.LD.Parsing.means_from_region_data(
            {i: region_ld[j] for i, j in enumerate(idx)}, stats
        )
        for idx in boot_idx
    ]
    data = {
        "scenario": scenario,
        "seed": seed,
        "regions": regions,
        "ld_region_stats": region_ld,
        "ld_mv": moments.LD.Parsing.bootstrap_data(region_ld),
        "ld_bootstraps": ld_boots,
        "sfs": region_sfs.sum(axis=0),
        "sfs_bootstraps": np.array([region_sfs[idx].sum(axis=0) for idx in boot_idx]),
        "boot_idx": boot_idx,
    }
    # write to a temporary name first, so an interrupted job leaves no partial file
    with open(out + ".tmp", "wb") as fout:
        pickle.dump(data, fout)
    os.replace(out + ".tmp", out)
    print(f"{out}: {regions} regions in {time.time() - t0:.0f} s", flush=True)
    return out


# -----------------------------------------------------------------------------
# Models and fits
# -----------------------------------------------------------------------------


class LDAnalysis:
    """Fits and Godambe methods for LD statistics. Parameters: (nu0, nu1, T, M, Ne)."""

    name = "ld"

    def __init__(self, data):
        self.means = data["ld_mv"]["means"]
        self.varcovs = data["ld_mv"]["varcovs"]
        self.bootstraps = data["ld_bootstraps"]
        # Statistics with the normalizing statistics removed, for the log-likelihood
        stats, self.means_r, self.varcovs_r, self.boots_r = (
            moments.LD.Godambe._remove_normalizing_statistics(
                moments.LD.Demographics2D.split_mig,
                true_params("mig"),
                False,
                self.means,
                self.varcovs,
                self.bootstraps,
                0,
                None,
            )
        )
        self.stats = stats
        self.cache = {}

    @staticmethod
    def model(params, rho=None, theta=0.001):
        return moments.LD.Demographics2D.split_mig(params, rho=rho, theta=theta)

    @staticmethod
    def nomig_model(params, rho=None, theta=0.001):
        return moments.LD.Demographics2D.split_mig([*params, 0], rho=rho, theta=theta)

    def fit(self, p_start, nested=False, data=None):
        means = self.means if data is None else data
        data = [means, self.varcovs]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if nested:
                p = np.delete(p_start, MIG_IDX)
                for _ in range(2):
                    p, neg_ll = moments.LD.Inference.optimize_log_fmin(
                        p, data, [self.nomig_model], rs=R_BINS
                    )
                return np.insert(p, MIG_IDX, 0), -neg_ll
            p, neg_ll = moments.LD.Inference.optimize_log_lbfgsb(
                p_start,
                data,
                [self.model],
                rs=R_BINS,
                lower_bound=LOWER,
                upper_bound=UPPER,
            )
            return np.array(p), -neg_ll

    def loglik(self, p, x):
        """Composite log-likelihood at p of data x (normalizing statistics removed)."""
        key = tuple(p)
        if key not in self.cache:
            y = moments.LD.Inference.bin_stats(
                self.model, p[:-1], rho=4 * p[-1] * R_BINS
            )
            y = moments.LD.Inference.sigmaD2(y)
            self.cache[key] = moments.LD.Inference.remove_nonpresent_statistics(
                y, self.stats
            )
        return moments.LD.Inference.ll_over_bins(x, self.cache[key], self.varcovs_r)

    def uncerts(self, p):
        G = moments.LD.Godambe
        return {
            "FIM": G.FIM_uncert(
                self.model, p, self.means, self.varcovs, r_edges=R_BINS
            ),
            "GIM": G.GIM_uncert(
                self.model, self.bootstraps, p, self.means, self.varcovs, r_edges=R_BINS
            ),
            "GIM (log)": G.GIM_uncert(
                self.model,
                self.bootstraps,
                p,
                self.means,
                self.varcovs,
                r_edges=R_BINS,
                log=True,
            ),
        }

    def lrt_adjust(self, p_nested):
        return moments.LD.Godambe.LRT_adjust(
            self.model,
            self.bootstraps,
            p_nested,
            [MIG_IDX],
            self.means,
            self.varcovs,
            r_edges=R_BINS,
        )

    def data_and_bootstraps(self):
        # data for loglik, with the normalizing statistics removed
        return self.means_r, self.boots_r

    def bootstrap_data(self, b):
        # data for fit, with all statistics, which the optimizer takes
        return self.bootstraps[b]


class SFSAnalysis:
    """
    Fits and Godambe methods for the joint SFS. The model functions take
    (nu0, nu1, T, M), and theta is fit by multinomial scaling. Parameter vectors here
    are (nu0, nu1, T, M, Ne), with Ne = theta / (4 U L_total).
    """

    name = "sfs"

    def __init__(self, data):
        ns = [SFS_SAMPLE_SIZE] * 2
        self.ns = ns
        self.L_total = data["regions"] * L
        self.fs = moments.Spectrum(data["sfs"], pop_ids=POPS).project(ns)
        self.bootstraps = [
            moments.Spectrum(b, pop_ids=POPS).project(ns)
            for b in data["sfs_bootstraps"]
        ]
        self.cache = {}

    @staticmethod
    def model(params, ns):
        return moments.Demographics2D.split_mig(params, ns)

    @staticmethod
    def nomig_model(params, ns):
        return moments.Demographics2D.split_mig([*params, 0], ns)

    def _with_Ne(self, p4, fs):
        theta = moments.Inference.optimal_sfs_scaling(self.model(p4, self.ns), fs)
        return np.append(p4, theta / (4 * U * self.L_total))

    def fit(self, p_start, nested=False, data=None):
        fs = self.fs if data is None else data
        p4 = np.array(p_start[:4], dtype=float)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if nested:
                p3 = np.delete(p4, MIG_IDX)
                for _ in range(2):
                    p3 = moments.Inference.optimize_log_fmin(p3, fs, self.nomig_model)
                p4 = np.insert(p3, MIG_IDX, 0)
            else:
                p4 = moments.Inference.optimize_log_lbfgsb(
                    p4, fs, self.model, lower_bound=LOWER[:4], upper_bound=UPPER[:4]
                )
            ll = moments.Inference.ll_multinom(self.model(p4, self.ns), fs)
        return self._with_Ne(np.array(p4), fs), ll

    def loglik(self, p, x):
        """Poisson log-likelihood at p = (nu0, nu1, T, M, Ne) of the spectrum x."""
        key = tuple(p[:4])
        if key not in self.cache:
            self.cache[key] = self.model(p[:4], self.ns)
        theta = 4 * p[4] * U * self.L_total
        return moments.Inference.ll(theta * self.cache[key], x)

    def uncerts(self, p):
        # The Godambe functions return uncertainties for (nu0, nu1, T, M, theta); the
        # relative uncertainty of theta is that of Ne.
        G = moments.Godambe
        theta = 4 * p[4] * U * self.L_total
        scale = np.array([1, 1, 1, 1, p[4] / theta])
        return {
            "FIM": G.FIM_uncert(self.model, p[:4], self.fs) * scale,
            "GIM": G.GIM_uncert(self.model, self.bootstraps, p[:4], self.fs) * scale,
            "GIM (log)": G.GIM_uncert(
                self.model, self.bootstraps, p[:4], self.fs, log=True
            ),
        }

    def lrt_adjust(self, p_nested):
        return moments.Godambe.LRT_adjust(
            self.model, self.bootstraps, p_nested[:4], self.fs, [MIG_IDX]
        )

    def data_and_bootstraps(self):
        return self.fs, self.bootstraps

    def bootstrap_data(self, b):
        return self.bootstraps[b]


ANALYSES = {"ld": LDAnalysis, "sfs": SFSAnalysis}


def fit_both(analysis, p_true):
    """
    Fits the nested model without migration, then the full model from the nested fit
    with a few values of the migration rate, and from the true parameters, keeping the
    best fit, so that its log-likelihood is at least that of the nested model.
    """
    p_nested, ll_nested = analysis.fit(p_true, nested=True)
    starts = [
        np.where(np.arange(5) == MIG_IDX, m0, p_nested) for m0 in [0.01, 0.1, 0.5]
    ]
    starts.append(np.maximum(p_true, [0, 0, 0, 0.01, 0]))
    best = (None, -np.inf)
    for p_start in starts:
        p, ll = analysis.fit(p_start)
        if ll > best[1]:
            best = (p, ll)
    return p_nested, ll_nested, best[0], best[1]


# -----------------------------------------------------------------------------
# Matrices for the LRT adjustment, over all parameters
# -----------------------------------------------------------------------------


def godambe_matrices(loglik, p0, x0, bootstraps, rel_step=1e-3, zero_step=1e-3):
    """
    Returns the hessian H of the negative log-likelihood at p0 for data x0, and the
    uncentred J and the mean score cU over the bootstrap replicates, with respect to all
    parameters. Uses central differences, except for parameters equal to zero (the
    migration rate at the nested fit), which use second-order forward differences, so
    the model is only evaluated at non-negative values.
    """
    p0 = np.array(p0, dtype=float)
    n = len(p0)
    h = np.where(p0 == 0, zero_step, rel_step * np.abs(p0))
    E = np.diag(h)
    forward = p0 == 0

    def d1(q, i, x):
        # first derivative in direction i at q
        if forward[i]:
            return (
                -3 * loglik(q, x) + 4 * loglik(q + E[i], x) - loglik(q + 2 * E[i], x)
            ) / (2 * h[i])
        return (loglik(q + E[i], x) - loglik(q - E[i], x)) / (2 * h[i])

    H = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            if i == j:
                if forward[i]:
                    d2 = (
                        2 * loglik(p0, x0)
                        - 5 * loglik(p0 + E[i], x0)
                        + 4 * loglik(p0 + 2 * E[i], x0)
                        - loglik(p0 + 3 * E[i], x0)
                    ) / h[i] ** 2
                else:
                    d2 = (
                        loglik(p0 + E[i], x0)
                        - 2 * loglik(p0, x0)
                        + loglik(p0 - E[i], x0)
                    ) / h[i] ** 2
            else:
                # derivative in direction a of the first derivative in direction b,
                # with a the forward direction, if there is one
                a, b = (i, j) if forward[i] else (j, i)
                if forward[a]:
                    g = [d1(p0 + s * E[a], b, x0) for s in range(3)]
                    d2 = (-3 * g[0] + 4 * g[1] - g[2]) / (2 * h[a])
                else:
                    d2 = (d1(p0 + E[a], b, x0) - d1(p0 - E[a], b, x0)) / (2 * h[a])
            H[i, j] = H[j, i] = -d2
    scores = np.array([[d1(p0, i, x) for i in range(n)] for x in bootstraps])
    J = scores.T @ scores / len(scores)
    return H, J, scores.mean(axis=0)


def adjustment_inverse_blocks(H, J, cU, nested, centre=True):
    """
    First-order moment matching with the blocks of H^-1 and H^-1 J H^-1 for the nested
    parameters: adj = q / tr[((H^-1)_nn)^-1 (H^-1 J H^-1)_nn].
    """
    if centre:
        J = J - np.outer(cU, cU)
    Hi = np.linalg.inv(H)
    A = Hi[np.ix_(nested, nested)]
    B = (Hi @ J @ Hi)[np.ix_(nested, nested)]
    return len(nested) / np.trace(np.linalg.inv(A) @ B)


def adjustment_submatrix(H, J, cU, nested, centre=False):
    """
    First-order moment matching with the submatrices of H and J for the nested
    parameters, as in LRT_adjust: adj = q / tr(J_nn H_nn^-1).
    """
    if centre:
        J = J - np.outer(cU, cU)
    sub = np.ix_(nested, nested)
    return len(nested) / np.trace(J[sub] @ np.linalg.inv(H[sub]))


# -----------------------------------------------------------------------------
# Analysis of one dataset
# -----------------------------------------------------------------------------


def _refit(args):
    method, path, b, p_full = args
    analysis = _load_analysis(method, path)
    return analysis.fit(p_full, data=analysis.bootstrap_data(b))[0]


_analysis_cache = {}


def _load_analysis(method, path):
    # each worker process loads the dataset once
    key = (method, path)
    if key not in _analysis_cache:
        with open(path, "rb") as fin:
            _analysis_cache[key] = ANALYSES[method](pickle.load(fin))
    return _analysis_cache[key]


def analyze_dataset(method, path, procs, boot_refits):
    with open(path, "rb") as fin:
        data = pickle.load(fin)
    scenario = data["scenario"]
    p_true = true_params(scenario)
    analysis = ANALYSES[method](data)
    res = {
        "method": method,
        "scenario": scenario,
        "seed": data["seed"],
        "regions": data["regions"],
        "p_true": p_true,
    }

    p_nested, ll_nested, p_full, ll_full = fit_both(analysis, p_true)
    res.update(p_nested=p_nested, ll_nested=ll_nested, p_full=p_full, ll_full=ll_full)
    res["D"] = max(2 * (ll_full - ll_nested), 0)

    # LRT adjustments at the nested fit
    res["adj_LRT_adjust"] = analysis.lrt_adjust(p_nested)
    x0, boots = analysis.data_and_bootstraps()
    res["H"], res["J"], res["cU"] = godambe_matrices(
        analysis.loglik, p_nested, x0, boots
    )

    # Uncertainties at the full fit, for datasets with migration (without migration,
    # the best fit is often on the boundary, where these intervals don't apply)
    if scenario != "null":
        res["sd"] = analysis.uncerts(p_full)
        args = [(method, path, b, p_full) for b in range(boot_refits)]
        res["p_boot"] = np.array(_map(_refit, args, procs))
    return res


def analyze(scenario, seed, method, procs, boot_refits):
    """Analyzes one saved dataset with one method, and saves the results."""
    path = data_path(scenario, seed)
    out = results_path(method, scenario, seed)
    if os.path.exists(out):
        print(f"{out}: already exists", flush=True)
        return
    os.makedirs(RESULTS_DIR, exist_ok=True)
    t0 = time.time()
    res = analyze_dataset(method, path, procs, boot_refits)
    with open(out + ".tmp", "wb") as fout:
        pickle.dump(res, fout)
    os.replace(out + ".tmp", out)
    print(f"{out}: analyzed in {time.time() - t0:.0f} s", flush=True)
    print_dataset_summary(res)


def print_dataset_summary(res):
    # a short summary of one analysis, for the job's log
    adj = {
        "LRT_adjust": res["adj_LRT_adjust"],
        "inverse blocks, centred J": adjustment_inverse_blocks(
            res["H"], res["J"], res["cU"], [MIG_IDX], True
        ),
    }
    np.set_printoptions(precision=4)
    print(f"  truth:        {res['p_true']}")
    print(f"  full fit:     {res['p_full']}")
    print(f"  nested fit:   {res['p_nested']}")
    print(f"  D = {res['D']:.4f}")
    for name, a in adj.items():
        p = moments.Godambe.sum_chi2_ppf(a * res["D"], weights=(0.5, 0.5))
        print(f"  adj ({name}) = {a:.4g}, p = {p:.4g}")
    for name, sd in res.get("sd", {}).items():
        print(f"  sd {name}: {sd}")
    if "p_boot" in res:
        print(f"  sd bootstrap: {res['p_boot'].std(axis=0, ddof=1)}")
    sys.stdout.flush()


def run(scenario, seed, regions, procs, boot_refits, methods):
    """Simulates one dataset and analyzes it with each method."""
    simulate(scenario, seed, regions, procs)
    for method in methods:
        analyze(scenario, seed, method, procs, boot_refits)


# -----------------------------------------------------------------------------
# Summaries
# -----------------------------------------------------------------------------


def load_results(method, scenario):
    results = []
    pattern = os.path.join(RESULTS_DIR, f"{method}_{scenario}_*.pkl")
    for path in sorted(glob.glob(pattern)):
        with open(path, "rb") as fin:
            results.append(pickle.load(fin))
    return results


def summarize_intervals(results, label):
    p_true = results[0]["p_true"]
    intervals = {}
    for r in results:
        p = r["p_full"]
        for name in ["FIM", "GIM"]:
            sd = r["sd"][name]
            intervals.setdefault(name, []).append((p - 1.96 * sd, p + 1.96 * sd))
        sd = r["sd"]["GIM (log)"]
        intervals.setdefault("GIM (log)", []).append(
            (p * np.exp(-1.96 * sd), p * np.exp(1.96 * sd))
        )
        sd = r["p_boot"].std(axis=0, ddof=1)
        intervals.setdefault("Bootstrap SE", []).append((p - 1.96 * sd, p + 1.96 * sd))
        intervals.setdefault("Bootstrap percentile", []).append(
            tuple(np.percentile(r["p_boot"], [2.5, 97.5], axis=0))
        )

    print(f"\n{label}: 95% confidence intervals, {len(results)} datasets")
    print("Coverage (should be about 0.95), with the mean width relative to the truth")
    print(f"{'':22s}" + "".join(f"{n:>14s}" for n in PARAM_NAMES))
    for name, ivs in intervals.items():
        lo = np.array([iv[0] for iv in ivs])
        hi = np.array([iv[1] for iv in ivs])
        cover = np.mean((lo <= p_true) & (p_true <= hi), axis=0)
        width = np.mean(hi - lo, axis=0) / p_true
        print(
            f"{name:22s}"
            + "".join(f"{c:7.2f} ({w:4.2f})" for c, w in zip(cover, width))
        )
    est = np.array([r["p_full"] for r in results])
    bias = (est.mean(0) - p_true) / est.std(0, ddof=1)
    print(f"{'Bias (sd units)':22s}" + "".join(f"{b:14.2f}" for b in bias))


def summarize_lrt(results, label):
    D = np.array([r["D"] for r in results])
    n = [MIG_IDX]
    adjustments = {
        "None": np.ones(len(results)),
        "LRT_adjust": np.array([r["adj_LRT_adjust"] for r in results]),
        "Submatrix, centred J": np.array(
            [adjustment_submatrix(r["H"], r["J"], r["cU"], n, True) for r in results]
        ),
        "Inverse blocks, centred J": np.array(
            [
                adjustment_inverse_blocks(r["H"], r["J"], r["cU"], n, True)
                for r in results
            ]
        ),
    }
    alphas = [0.01, 0.05, 0.1, 0.25]
    print(
        f"\n{label}: LRT for migration, {len(results)} datasets "
        f"(fraction with D = 0: {np.mean(D < 1e-8):.2f})"
    )
    print("Fraction of p-values below alpha (equals alpha under the null)")
    print(
        f"{'':28s}{'median adj':>12s}" + "".join(f"{'a=' + str(a):>9s}" for a in alphas)
    )
    for name, adj in adjustments.items():
        p = moments.Godambe.sum_chi2_ppf(adj * D, weights=(0.5, 0.5))
        print(
            f"{name:28s}{np.median(adj):12.4f}"
            + "".join(f"{np.mean(p < a):9.2f}" for a in alphas)
        )


def summarize():
    for method in ANALYSES:
        for scenario in SCENARIOS:
            results = load_results(method, scenario)
            if len(results) == 0:
                continue
            label = (
                f"{method.upper()}, {scenario} (M = {true_params(scenario)[MIG_IDX]})"
            )
            if scenario != "null":
                summarize_intervals(results, label)
            summarize_lrt(results, label)


# -----------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="stage", required=True)
    p_run = sub.add_parser("run", help="simulate one dataset and analyze it")
    p_sim = sub.add_parser("simulate", help="simulate and parse one dataset")
    p_ana = sub.add_parser("analyze", help="fit models, uncertainties, and LRTs")
    sub.add_parser("summarize", help="print CI coverage and LRT calibration")
    for p in [p_run, p_sim, p_ana]:
        p.add_argument("--scenario", choices=SCENARIOS, required=True)
        p.add_argument("--seed", type=int, required=True)
        p.add_argument("--procs", type=int, default=1)
    for p in [p_run, p_sim]:
        p.add_argument("--regions", type=int, default=REGIONS)
    for p in [p_run, p_ana]:
        p.add_argument("--boot-refits", type=int, default=NUM_BOOTSTRAPS)
    p_run.add_argument("--methods", nargs="+", choices=ANALYSES, default=list(ANALYSES))
    p_ana.add_argument("--method", choices=ANALYSES, required=True)
    args = parser.parse_args()

    if args.stage == "run":
        run(
            args.scenario,
            args.seed,
            args.regions,
            args.procs,
            args.boot_refits,
            args.methods,
        )
    elif args.stage == "simulate":
        simulate(args.scenario, args.seed, args.regions, args.procs)
    elif args.stage == "analyze":
        analyze(args.scenario, args.seed, args.method, args.procs, args.boot_refits)
    else:
        summarize()


if __name__ == "__main__":
    main()
