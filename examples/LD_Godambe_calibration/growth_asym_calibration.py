"""
A second calibration model for the Godambe LRT adjustments and the Hessian warning:
a split in which deme1 grows exponentially and migration is asymmetric.

The model is built with demes in this script (no library model functions):

    ancestor (size Ne) splits 2 Ne T generations ago into
      deme0, constant size nu0 Ne, and
      deme1, growing exponentially from nu1_0 Ne to nu1_F Ne,
    with migration from deme0 to deme1 at M01 = 2 Ne m01 and from deme1 to deme0 at
    M10 = 2 Ne m10 (demes source -> dest, forward in time).

Parameters (genetic units): (nu0, nu1_0, nu1_F, T, M01, M10, Ne). For LD, Ne is fit
(rho = 4 Ne r); for the SFS, theta is fit by multinomial scaling and reported as
Ne = theta / (4 u L_total).

Four likelihood ratio tests are run on every dataset, each against the full model:

    asym    nested M01 = M10              interior null, q = 1 (chi^2_1)
    mig     nested M01 = M10 = 0          two boundary parameters, q = 2 (chi-bar^2)
    oneway  nested M10 = 0                one boundary parameter, q = 1, with M01 free
    growth  nested nu1_F = nu1_0          interior null, q = 1 (chi^2_1)

For the interior tests, the Godambe matrices are computed in coordinates where the
nested parameter is a coordinate: (M, a) with M01 = M e^a and M10 = M e^-a for asym,
and g with nu1_F = nu1_0 e^g for growth (both 0 under the nested model, central
differences). For the boundary tests, they are computed in (M01, M10), with
second-order forward differences at zero.

Per test and dataset we save D, the nested fit, H, J and the mean score cU at the nested
fit (in the test's coordinates), the library's LRT_adjust, the conditioning measures
used by the warning and the fallback, the chi-bar-squared weights for the q = 2 test,
and a Monte Carlo p-value from the asymptotic null distribution of D with the
estimated H and J. At the full fit, FIM and GIM standard errors and bootstrap refits
give CIs and Wald tests.

Scenarios (true values; Ne = 10,000, nu0 = 0.2, T = 0.075):

    New datasets, simulated here (deme1 grows from 0.2 to 2):
      g_nomig    M01 = M10 = 0     null for mig and oneway
      g_sym      M01 = M10 = 1     null for asym; power for mig
      g_oneway   M01 = 1, M10 = 0  null for oneway; power for asym and mig
      g_asym     M01 = 1, M10 = 0.2  power for asym, mig and oneway; CIs
    Datasets from the second run of godambe_calibration.py (constant sizes, symmetric
    migration), reanalyzed with this model without simulating again:
      mig (M = 2), weak (M = 0.2): null for asym and growth
      null (M = 0): null for mig, oneway and growth

All four tests are run in every scenario, so each also gives power where its null is
false.

    python growth_asym_calibration.py simulate --scenario g_sym --seed 1 --procs 4
    python growth_asym_calibration.py analyze --scenario g_sym --seed 1 --method ld
    python growth_asym_calibration.py run --scenario g_sym --seed 1     # both stages
    python growth_asym_calibration.py analyze --scenario mig --seed 1 --method sfs
    python growth_asym_calibration.py summarize

Datasets are read from and written to ./data (the same format and names as
godambe_calibration.py, so its datasets can be reanalyzed), results go to
./results_growth/{method}_{scenario}_{seed}.pkl. Stages skip outputs that exist.
Imports godambe_calibration.py (in the same directory) for the data handling, fit
chains and adjustment formulas; that script is not changed.
"""

import argparse
import contextlib
import glob
import io
import os
import pickle
import sys
import tempfile
import time
import warnings
from itertools import combinations

import demes
import msprime
import numpy as np

import moments
import moments.LD

import godambe_calibration as base

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------

NE = base.NE
NREF = 10000  # reference size for building graphs from genetic-unit parameters
POPS = base.POPS
PARAM_NAMES = ["nu0", "nu1_0", "nu1_F", "T", "M01", "M10", "Ne"]
RATES = [4, 5]  # indices of the migration rates (boundary at zero)
LOWER = np.array([1e-3, 1e-3, 1e-3, 1e-3, 0, 0, 100])
UPPER = np.array([100, 100, 100, 10, 100, 100, 1e6])

# True values, in genetic units. "source" says whether the dataset is simulated here
# ("new") or comes from the second run of godambe_calibration.py ("run2").
_NU0, _T = base.N0 / NE, base.T_GEN / 2 / NE
SCENARIOS = {
    "g_nomig": dict(p=[_NU0, 0.2, 2, _T, 0, 0, NE], source="new"),
    "g_sym": dict(p=[_NU0, 0.2, 2, _T, 1, 1, NE], source="new"),
    "g_oneway": dict(p=[_NU0, 0.2, 2, _T, 1, 0, NE], source="new"),
    "g_asym": dict(p=[_NU0, 0.2, 2, _T, 1, 0.2, NE], source="new"),
}
for _name, _m in base.SCENARIOS.items():
    _M = 2 * NE * _m
    SCENARIOS[_name] = dict(p=[_NU0, 2, 2, _T, _M, _M, NE], source="run2")

# The number of bootstrap refits of the full model per dataset (fewer than in the
# first model: the growth model is slower to fit, and the CIs are a side result here)
BOOT_REFITS = 100

# Monte Carlo draws for the asymptotic null distribution of D
MC_DRAWS = 20000

DATA_DIR = base.DATA_DIR
RESULTS_DIR = "results_growth"


def true_params(scenario):
    return np.array(SCENARIOS[scenario]["p"], dtype=float)


def build_graph(d, nref=NREF, m_floor=0.0):
    """
    The demes graph for demographic parameters d = (nu0, nu1_0, nu1_F, T, M01, M10),
    with reference size nref. Migration rates below m_floor are set to m_floor.
    """
    nu0, nu1_0, nu1_F, T, M01, M10 = d
    b = demes.Builder(time_units="generations")
    b.add_deme("anc", epochs=[dict(start_size=nref, end_time=2 * nref * T)])
    b.add_deme("deme0", ancestors=["anc"], epochs=[dict(start_size=nu0 * nref)])
    b.add_deme(
        "deme1",
        ancestors=["anc"],
        epochs=[dict(start_size=nu1_0 * nref, end_size=nu1_F * nref)],
    )
    for source, dest, M in [("deme0", "deme1", M01), ("deme1", "deme0", M10)]:
        M = max(M, m_floor)
        if M > 0:
            b.add_migration(source=source, dest=dest, rate=M / (2 * nref))
    return b.resolve()


# -----------------------------------------------------------------------------
# The four tests: nested models and coordinates
# -----------------------------------------------------------------------------
#
# Each test gives
#   free_to_d(f), d_to_free(d): the nested model's free demographic parameters (no Ne)
#     to and from the full model's demographic parameters d (6 values);
#   to_test(p), from_test(q): full parameter vectors (with Ne) to and from the
#     coordinates in which the nested parameters psi are coordinates, fixed at 0;
#   psi: the indices of the nested parameters in the test coordinates;
#   boundary: whether psi = 0 is on the boundary (psi >= 0);
#   rates: indices (test coordinates) of parameters that must stay >= 0, which use
#     forward differences when they are zero.


def _asym_to_test(p):
    p = np.array(p, dtype=float)
    M01, M10 = p[4], p[5]
    q = p.copy()
    q[4] = np.sqrt(M01 * M10)
    q[5] = 0.5 * np.log(M01 / M10) if M01 > 0 and M10 > 0 else 0.0
    return q


def _asym_from_test(q):
    p = np.array(q, dtype=float)
    p[4], p[5] = q[4] * np.exp(q[5]), q[4] * np.exp(-q[5])
    return p


def _growth_to_test(p):
    q = np.array(p, dtype=float)
    q[2] = np.log(p[2] / p[1])
    return q


def _growth_from_test(q):
    p = np.array(q, dtype=float)
    p[2] = q[1] * np.exp(q[2])
    return p


TESTS = {
    "asym": dict(
        free_to_d=lambda f: np.array([f[0], f[1], f[2], f[3], f[4], f[4]]),
        d_to_free=lambda d: np.array([d[0], d[1], d[2], d[3], np.sqrt(d[4] * d[5])]),
        free_bounds=[0, 1, 2, 3, 4],
        to_test=_asym_to_test,
        from_test=_asym_from_test,
        psi=[5],
        boundary=False,
        rates=[4],
    ),
    "mig": dict(
        free_to_d=lambda f: np.array([f[0], f[1], f[2], f[3], 0, 0]),
        d_to_free=lambda d: np.array(d[:4]),
        free_bounds=[0, 1, 2, 3],
        to_test=lambda p: np.array(p, dtype=float),
        from_test=lambda q: np.array(q, dtype=float),
        psi=[4, 5],
        boundary=True,
        rates=[4, 5],
    ),
    "oneway": dict(
        free_to_d=lambda f: np.array([f[0], f[1], f[2], f[3], f[4], 0]),
        d_to_free=lambda d: np.array(d[:5]),
        free_bounds=[0, 1, 2, 3, 4],
        to_test=lambda p: np.array(p, dtype=float),
        from_test=lambda q: np.array(q, dtype=float),
        psi=[5],
        boundary=True,
        rates=[4, 5],
    ),
    "growth": dict(
        free_to_d=lambda f: np.array([f[0], f[1], f[1], f[2], f[3], f[4]]),
        d_to_free=lambda d: np.array([d[0], np.sqrt(d[1] * d[2]), d[3], d[4], d[5]]),
        free_bounds=[0, 1, 3, 4, 5],
        to_test=_growth_to_test,
        from_test=_growth_from_test,
        psi=[2],
        boundary=False,
        rates=[4, 5],
    ),
}


def null_is_true(test, p_true):
    nu1_0, nu1_F, M01, M10 = p_true[1], p_true[2], p_true[4], p_true[5]
    return {
        "asym": M01 == M10 and M01 > 0,
        "mig": M01 == 0 and M10 == 0,
        "oneway": M10 == 0,
        "growth": nu1_0 == nu1_F,
    }[test]


# -----------------------------------------------------------------------------
# Simulating and parsing data (new scenarios only)
# -----------------------------------------------------------------------------


def _seed_sequence(scenario, seed, *keys):
    # independent of the streams in godambe_calibration.py (which start at 0-2)
    return np.random.SeedSequence([100 + list(SCENARIOS).index(scenario), seed, *keys])


def _region_seeds(scenario, seed, region):
    ss = _seed_sequence(scenario, seed, 0, region)
    return [int(s) % (2**31 - 2) + 1 for s in ss.generate_state(2)]


def simulate_region(args):
    """As base.simulate_region, for the growth model."""
    scenario, seed, region, tmpdir = args
    anc_seed, mut_seed = _region_seeds(scenario, seed, region)
    p = true_params(scenario)
    demog = msprime.Demography.from_demes(build_graph(p[:6], nref=p[6]))
    ts = msprime.sim_ancestry(
        {pop: base.NUM_SAMPLES for pop in POPS},
        demography=demog,
        sequence_length=base.L,
        recombination_rate=base.R,
        random_seed=anc_seed,
    )
    ts = msprime.sim_mutations(
        ts, rate=base.U, random_seed=mut_seed, discrete_genome=False
    )
    vcf_file = os.path.join(tmpdir, f"region_{region}.vcf")
    names = [f"ind_{i}" for i in range(ts.num_individuals)]
    with open(vcf_file, "w") as fout:
        ts.write_vcf(
            fout,
            individual_names=names,
            allow_position_zero=True,
            position_transform="legacy",
        )
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
            r_bins=base.R_BINS,
            use_h5=False,
            report=False,
        )
        sfs = moments.Spectrum.from_vcf(vcf_file, pop_mapping=base.sample_ids())
    os.remove(vcf_file)
    return ld_stats, np.asarray(sfs.data)


def simulate(scenario, seed, regions, procs):
    """Simulates and parses one dataset of a new scenario, as base.simulate does."""
    if SCENARIOS[scenario]["source"] != "new":
        raise ValueError(f"{scenario} datasets come from godambe_calibration.py")
    out = base.data_path(scenario, seed)
    if os.path.exists(out):
        print(f"{out}: already exists", flush=True)
        return out
    os.makedirs(DATA_DIR, exist_ok=True)
    t0 = time.time()
    with tempfile.TemporaryDirectory() as tmpdir:
        base.write_pop_file_and_map(tmpdir)
        args = [(scenario, seed, i, tmpdir) for i in range(regions)]
        outputs = base._map(simulate_region, args, procs)
    region_ld = {i: o[0] for i, o in enumerate(outputs)}
    region_sfs = np.array([o[1] for o in outputs])
    rng = np.random.default_rng(_seed_sequence(scenario, seed, 1))
    boot_idx = rng.integers(regions, size=(base.NUM_BOOTSTRAPS, regions))
    # Compact storage (see load_data): each region's SFS projected to the inference
    # sample size, in float32, instead of 500 full-size bootstrap spectra (~26 MB); the
    # bootstraps, LD and SFS, are rebuilt from these and boot_idx when loading.
    ns = [base.SFS_SAMPLE_SIZE] * 2
    data = {
        "format": "compact",
        "scenario": scenario,
        "seed": seed,
        "regions": regions,
        "ld_region_stats": region_ld,
        "ld_mv": moments.LD.Parsing.bootstrap_data(region_ld),
        "sfs": region_sfs.sum(axis=0),
        "sfs_regions_projected": np.array(
            [moments.Spectrum(x).project(ns).data for x in region_sfs],
            dtype=np.float32,
        ),
        "boot_idx": boot_idx.astype(np.int16),
    }
    with open(out + ".tmp", "wb") as fout:
        pickle.dump(data, fout)
    os.replace(out + ".tmp", out)
    print(f"{out}: {regions} regions in {time.time() - t0:.0f} s", flush=True)
    return out


def load_data(path):
    """
    Loads a dataset. Compact datasets (from simulate here) store the projected SFS of
    each region; their LD and SFS bootstrap replicates are rebuilt here, in the format
    of godambe_calibration.py datasets. Since projection is linear, the summed
    projected spectra equal the projections of the summed spectra (up to float32
    rounding).
    """
    with open(path, "rb") as fin:
        data = pickle.load(fin)
    if data.get("format") != "compact":
        return data
    boot_idx = data["boot_idx"].astype(int)
    region_ld = data["ld_region_stats"]
    stats = region_ld[0]["stats"]
    data["ld_bootstraps"] = [
        moments.LD.Parsing.means_from_region_data(
            {i: region_ld[j] for i, j in enumerate(idx)}, stats
        )
        for idx in boot_idx
    ]
    proj = data["sfs_regions_projected"].astype(float)
    data["sfs_bootstraps"] = np.array([proj[idx].sum(axis=0) for idx in boot_idx])
    return data


# -----------------------------------------------------------------------------
# Models and fits
# -----------------------------------------------------------------------------


class LDModel:
    """LD statistics for the growth model. Parameter vectors are (d, Ne)."""

    name = "ld"

    def __init__(self, data):
        # data handling (normalizing statistics, bootstraps) as in the first model
        self.base = base.LDAnalysis(data)
        self.means, self.varcovs = self.base.means, self.base.varcovs
        self.bootstraps = self.base.bootstraps
        self.cache = {}
        self.scenario, self.seed = data["scenario"], data["seed"]

    @staticmethod
    def expectations(d, rho=None, theta=0.001):
        return moments.LD.LDstats.from_demes(build_graph(d), POPS, rho=rho, theta=theta)

    def fit(self, free_to_d, f_start, bounds, data=None, key=0, chain=base.FIT_CHAIN):
        """
        Fits the free parameters f (with Ne last) of a model given by free_to_d.
        bounds are the indices in d whose bounds apply to the free parameters.
        Returns the full parameter vector (d, Ne), the log-likelihood, and the trace.
        """
        means = self.means if data is None else data
        lower = np.append(LOWER[bounds], LOWER[6])
        upper = np.append(UPPER[bounds], UPPER[6])

        def model(f, rho=None, theta=0.001):
            return self.expectations(free_to_d(f), rho=rho, theta=theta)

        def optimize(name, f):
            optimizer = getattr(moments.LD.Inference, f"optimize_log_{name}")
            f, neg_ll = optimizer(
                f,
                [means, self.varcovs],
                [model],
                rs=base.R_BINS,
                lower_bound=lower,
                upper_bound=upper,
            )
            return f, -neg_ll

        f0 = np.clip(f_start, lower, upper)
        f, ll, lls = base.run_fit_chain(
            optimize, f0, lower, upper, _fit_rng(self, key), chain
        )
        return np.append(free_to_d(f[:-1]), f[-1]), ll, lls

    def loglik(self, p, x):
        """Composite log-likelihood at p = (d, Ne) of data x (normalizing removed)."""
        key = tuple(p)
        if key not in self.cache:
            y = moments.LD.Inference.bin_stats(
                self.expectations, p[:6], rho=4 * p[6] * base.R_BINS
            )
            y = moments.LD.Inference.sigmaD2(y)
            self.cache[key] = moments.LD.Inference.remove_nonpresent_statistics(
                y, self.base.stats
            )
        return moments.LD.Inference.ll_over_bins(
            x, self.cache[key], self.base.varcovs_r
        )

    def lrt_adjust(self, test, q_nested):
        # the library function, with the model in the test's coordinates
        from_test = TESTS[test]["from_test"]

        def model(qd, rho=None, theta=0.001):
            return self.expectations(
                from_test(np.append(qd, 1))[:6], rho=rho, theta=theta
            )

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return moments.LD.Godambe.LRT_adjust(
                model,
                self.bootstraps,
                q_nested,
                TESTS[test]["psi"],
                self.means,
                self.varcovs,
                r_edges=base.R_BINS,
            )

    def data_and_bootstraps(self):
        return self.base.means_r, self.base.boots_r

    def bootstrap_data(self, b):
        return self.bootstraps[b]


class SFSModel:
    """
    The joint SFS for the growth model. The optimizer fits d, with theta by
    multinomial scaling; parameter vectors are (d, Ne), Ne = theta / (4 u L_total).
    Migration rates are floored at base.M_FLOOR (see base.SFSAnalysis.model).
    """

    name = "sfs"

    def __init__(self, data):
        self.base = base.SFSAnalysis(data)
        self.ns, self.fs, self.bootstraps = (
            self.base.ns,
            self.base.fs,
            self.base.bootstraps,
        )
        self.L_total = self.base.L_total
        self.cache = {}
        self.scenario, self.seed = data["scenario"], data["seed"]

    @staticmethod
    def expectations(d, ns):
        g = build_graph(d, m_floor=base.M_FLOOR)
        return moments.Spectrum.from_demes(g, sampled_demes=POPS, sample_sizes=ns)

    def _with_Ne(self, d, fs):
        theta = moments.Inference.optimal_sfs_scaling(self.expectations(d, self.ns), fs)
        return np.append(d, theta / (4 * base.U * self.L_total))

    def fit(self, free_to_d, f_start, bounds, data=None, key=0, chain=base.FIT_CHAIN):
        """As LDModel.fit; f_start may include Ne last, which is ignored."""
        fs = self.fs if data is None else data
        lower, upper = LOWER[bounds], UPPER[bounds]

        def model(f, ns):
            return self.expectations(free_to_d(f), ns)

        def optimize(name, f):
            optimizer = getattr(moments.Inference, f"optimize_log_{name}")
            f = optimizer(f, fs, model, lower_bound=lower, upper_bound=upper)
            return f, moments.Inference.ll_multinom(model(f, self.ns), fs)

        f0 = np.clip(f_start[: len(bounds)], lower, upper)
        f, ll, lls = base.run_fit_chain(
            optimize, f0, lower, upper, _fit_rng(self, key), chain
        )
        return self._with_Ne(free_to_d(f), fs), ll, lls

    def loglik(self, p, x):
        key = tuple(p[:6])
        if key not in self.cache:
            self.cache[key] = self.expectations(p[:6], self.ns)
        theta = 4 * p[6] * base.U * self.L_total
        return moments.Inference.ll(theta * self.cache[key], x)

    def lrt_adjust(self, test, q_nested):
        from_test = TESTS[test]["from_test"]

        def model(qd, ns):
            return self.expectations(from_test(np.append(qd, 1))[:6], ns)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return moments.Godambe.LRT_adjust(
                model, self.bootstraps, q_nested[:6], self.fs, TESTS[test]["psi"]
            )

    def data_and_bootstraps(self):
        return self.fs, self.bootstraps

    def bootstrap_data(self, b):
        return self.bootstraps[b]


MODELS = {"ld": LDModel, "sfs": SFSModel}


def _fit_rng(model, key):
    method_idx = list(MODELS).index(model.name)
    ss = _seed_sequence(model.scenario, model.seed, 2, method_idx, key)
    return np.random.default_rng(ss)


def _start(f, floor=0.01):
    # log-scale optimizers can't start at zero
    return np.where(np.asarray(f) <= 0, floor, f)


def fit_nested(model, test, p_true, key):
    t = TESTS[test]
    f0 = _start(np.append(t["d_to_free"](p_true[:6]), p_true[6]))
    return model.fit(t["free_to_d"], f0, t["free_bounds"], key=key)


def fit_full(model, starts, key0):
    """Fits the full model from each start, and keeps the best fit."""
    identity = list(range(6))
    best, traces = (None, -np.inf), []
    for k, p in enumerate(starts):
        p_fit, ll, trace = model.fit(
            lambda f: np.asarray(f, dtype=float), _start(p), identity, key=key0 + k
        )
        traces.append(trace)
        if ll > best[1]:
            best = (p_fit, ll)
    return best[0], best[1], traces


# -----------------------------------------------------------------------------
# Godambe matrices, adjustments, and the asymptotic null distribution of D
# -----------------------------------------------------------------------------


def godambe_matrices(
    loglik, q0, x0, bootstraps, forward, rel_step=1e-3, zero_step=1e-3
):
    """
    As base.godambe_matrices, in any coordinates: H of the negative log-likelihood
    and the uncentred J and mean score cU over the bootstrap replicates, at q0.
    forward[i] marks the parameters at zero on a boundary (second-order forward
    differences, so they are never evaluated below zero); others use central
    differences, with step rel_step * |q| or zero_step when q is zero.
    """
    q0 = np.array(q0, dtype=float)
    n = len(q0)
    h = np.where(q0 == 0, zero_step, rel_step * np.abs(q0))
    E = np.diag(h)

    def d1(q, i, x):
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
                        2 * loglik(q0, x0)
                        - 5 * loglik(q0 + E[i], x0)
                        + 4 * loglik(q0 + 2 * E[i], x0)
                        - loglik(q0 + 3 * E[i], x0)
                    ) / h[i] ** 2
                else:
                    d2 = (
                        loglik(q0 + E[i], x0)
                        - 2 * loglik(q0, x0)
                        + loglik(q0 - E[i], x0)
                    ) / h[i] ** 2
            else:
                a, b = (i, j) if forward[i] else (j, i)
                if forward[a]:
                    g = [d1(q0 + s * E[a], b, x0) for s in range(3)]
                    d2 = (-3 * g[0] + 4 * g[1] - g[2]) / (2 * h[a])
                else:
                    d2 = (d1(q0 + E[a], b, x0) - d1(q0 - E[a], b, x0)) / (2 * h[a])
            H[i, j] = H[j, i] = -d2
    scores = np.array([[d1(q0, i, x) for i in range(n)] for x in bootstraps])
    return H, scores.T @ scores / len(scores), scores.mean(axis=0)


def conditioning(H, psi):
    """
    The smallest eigenvalue of the correlation-scaled H (the library warning's
    measure), and the smallest eigenvalue of H_pp^(-1/2) S H_pp^(-1/2), with S the
    Schur complement for psi (the share of the psi curvature that survives refitting
    the other parameters; S / H_psipsi for q = 1).
    """
    d = np.diag(H)
    if np.any(d <= 0) or not np.all(np.isfinite(H)):
        return np.nan, np.nan
    s = np.sqrt(d)
    min_eig = np.linalg.eigvalsh(H / np.outer(s, s))[0]
    tau = [i for i in range(len(H)) if i not in psi]
    Hpp = H[np.ix_(psi, psi)]
    S = Hpp - H[np.ix_(psi, tau)] @ np.linalg.solve(
        H[np.ix_(tau, tau)], H[np.ix_(tau, psi)]
    )
    try:
        Li = np.linalg.inv(np.linalg.cholesky(Hpp))
    except np.linalg.LinAlgError:
        return min_eig, np.nan
    return min_eig, np.linalg.eigvalsh(Li @ S @ Li.T)[0]


def psi_covariances(H, J, cU, psi):
    """The psi blocks of H^-1 and of H^-1 J H^-1 (J centred)."""
    Jc = J - np.outer(cU, cU)
    Hi = np.linalg.inv(H)
    sub = np.ix_(psi, psi)
    return Hi[sub], (Hi @ Jc @ Hi)[sub]


def chibar_weights(cov):
    """Chi-bar-squared weights (0, 1, 2 df) for two parameters with covariance cov."""
    # w2 = P(both estimates positive) = 1/4 + arcsin(rho) / (2 pi), w1 = 1/2
    rho = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
    w2 = 0.25 + np.arcsin(np.clip(rho, -1, 1)) / (2 * np.pi)
    return np.array([0.5 - w2, 0.5, w2])


def mc_null(H, J, cU, psi, boundary, rng, draws=MC_DRAWS):
    """
    Draws from the asymptotic null distribution of D: with Z ~ N(0, G^pp) and
    A = (H^pp)^-1, D = Z'AZ in the interior, and Z'AZ - min_{x >= 0} (Z-x)'A(Z-x) on
    the boundary (psi >= 0). Returns the draws, or None if A is not positive definite.
    """
    Hpp_inv, Gpp = psi_covariances(H, J, cU, psi)
    A = np.linalg.inv(Hpp_inv)
    if np.linalg.eigvalsh(A)[0] <= 0 or np.linalg.eigvalsh(Gpp)[0] <= 0:
        return None
    Z = rng.multivariate_normal(np.zeros(len(psi)), Gpp, size=draws)
    full = np.einsum("ni,ij,nj->n", Z, A, Z)
    if not boundary:
        return full
    # min over the orthant of (Z-x)'A(Z-x): the best feasible stationary point over the
    # faces (each subset of free coordinates, the others at zero)
    q = len(psi)
    best = full.copy()  # x = 0
    for k in range(1, q + 1):
        for free in combinations(range(q), k):
            free = list(free)
            fixed = [i for i in range(q) if i not in free]
            # minimize over x_free with x_fixed = 0: x_free = Z_free + A_ff^-1 A_fc Z_c
            Aff = A[np.ix_(free, free)]
            x = Z[:, free] + (
                np.linalg.solve(Aff, A[np.ix_(free, fixed)] @ Z[:, fixed].T).T
                if fixed
                else 0
            )
            X = np.zeros_like(Z)
            X[:, free] = x
            R = Z - X
            val = np.einsum("ni,ij,nj->n", R, A, R)
            ok = np.all(x >= 0, axis=1)
            best = np.where(ok & (val < best), val, best)
    return full - best


def test_quantities(model, test, p_nested, D, seed_key):
    """H, J, cU, adjustments, conditioning and Monte Carlo p-value for one test."""
    t = TESTS[test]
    psi = t["psi"]
    q0 = t["to_test"](p_nested)
    forward = [(i in t["rates"]) and q0[i] == 0 for i in range(len(q0))]
    x0, boots = model.data_and_bootstraps()

    def loglik(q, x):
        return model.loglik(t["from_test"](q), x)

    H, J, cU = godambe_matrices(loglik, q0, x0, boots, forward)
    out = dict(q_nested=q0, H=H, J=J, cU=cU)
    out["adj_LRT_adjust"] = model.lrt_adjust(test, q0)
    out["adj_submatrix"] = base.adjustment_submatrix(H, J, cU, psi, True)
    out["adj_inverse_blocks"] = base.adjustment_inverse_blocks(H, J, cU, psi, True)
    out["min_eig"], out["schur_ratio"] = conditioning(H, psi)
    Hpp_inv, Gpp = psi_covariances(H, J, cU, psi)
    if len(psi) == 2:
        out["chibar_w_G"] = chibar_weights(Gpp)
        out["chibar_w_H"] = chibar_weights(Hpp_inv)
    rng = np.random.default_rng(seed_key)
    null = mc_null(H, J, cU, psi, t["boundary"], rng)
    out["p_mc"] = np.nan if null is None else np.mean(null >= D)
    return out


# -----------------------------------------------------------------------------
# Analysis of one dataset
# -----------------------------------------------------------------------------

_model_cache = {}


def _load_model(method, path):
    key = (method, path)
    if key not in _model_cache:
        _model_cache[key] = MODELS[method](load_data(path))
    return _model_cache[key]


def _refit(args):
    method, path, b, p_full = args
    model = _load_model(method, path)
    return model.fit(
        lambda f: np.asarray(f, dtype=float),
        _start(p_full),
        list(range(6)),
        data=model.bootstrap_data(b),
        key=1000 + b,
        chain=base.REFIT_CHAIN,
    )


def analyze_dataset(method, path, procs, boot_refits):
    data = load_data(path)
    scenario = data["scenario"]
    p_true = true_params(scenario)
    model = MODELS[method](data)
    res = dict(
        method=method,
        scenario=scenario,
        seed=data["seed"],
        regions=data["regions"],
        p_true=p_true,
        fit_chain=base.FIT_CHAIN,
        refit_chain=base.REFIT_CHAIN,
    )

    # nested fits, then the full fit from each nested optimum and from the truth
    nested = {}
    for k, test in enumerate(TESTS):
        nested[test] = fit_nested(model, test, p_true, key=k)
    starts = [p for p, _, _ in nested.values()] + [p_true]
    p_full, ll_full, traces_full = fit_full(model, starts, key0=10)
    res.update(p_full=p_full, ll_full=ll_full, traces_full=traces_full)

    res["tests"] = {}
    for k, (test, (p_nested, ll_nested, trace)) in enumerate(nested.items()):
        D = max(2 * (ll_full - ll_nested), 0)
        out = dict(
            p_nested=p_nested,
            ll_nested=ll_nested,
            trace=trace,
            D=D,
            null_true=null_is_true(test, p_true),
        )
        seed_key = _seed_sequence(scenario, data["seed"], 3, k)
        out.update(test_quantities(model, test, p_nested, D, seed_key))
        res["tests"][test] = out

    # FIM and GIM at the full fit, in the model's coordinates, with forward
    # differences for migration rates at zero
    forward = [i in RATES and p_full[i] == 0 for i in range(7)]
    x0, boots = model.data_and_bootstraps()
    H, J, cU = godambe_matrices(model.loglik, p_full, x0, boots, forward)
    res.update(H_full=H, J_full=J, cU_full=cU)
    try:
        Hi = np.linalg.inv(H)
        res["sd"] = {
            "FIM": np.sqrt(np.diag(Hi)),
            "GIM": np.sqrt(np.diag(Hi @ (J - np.outer(cU, cU)) @ Hi)),
        }
    except np.linalg.LinAlgError:
        res["sd"] = None

    if boot_refits > 0:
        args = [(method, path, b, p_full) for b in range(boot_refits)]
        refits = base._map(_refit, args, procs)
        res["p_boot"] = np.array([r[0] for r in refits])
        res["boot_traces"] = np.array([r[2] for r in refits])
    return res


def results_path(method, scenario, seed):
    return os.path.join(RESULTS_DIR, f"{method}_{scenario}_{seed}.pkl")


def analyze(scenario, seed, method, procs, boot_refits):
    path = base.data_path(scenario, seed)
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


def _p_values(t, D, test):
    """p-values for one test under each adjustment."""
    if test == "mig":
        weights = {"1/4,1/2,1/4": (0.25, 0.5, 0.25), "G": t["chibar_w_G"]}
    elif TESTS[test]["boundary"]:
        weights = {"": (0.5, 0.5)}
    else:
        weights = {"": (0, 1)}
    out = {}
    for name in ["LRT_adjust", "submatrix", "inverse_blocks"]:
        for wname, w in weights.items():
            label = f"{name} ({wname})" if wname else name
            out[label] = moments.Godambe.sum_chi2_ppf(t[f"adj_{name}"] * D, w)
    out["Monte Carlo"] = t["p_mc"]
    return out


def print_dataset_summary(res):
    np.set_printoptions(precision=4, suppress=True)
    print(f"  truth:    {res['p_true']}")
    print(f"  full fit: {res['p_full']}  ll {res['ll_full']:.3f}")
    for trace in res["traces_full"]:
        print(f"    full fit start, ll per step: {np.array(trace)}")
    for test, t in res["tests"].items():
        print(
            f"  {test:7s} D = {t['D']:.3f} (null {'true' if t['null_true'] else 'false'}),"
            f" min eig {t['min_eig']:.4f}, Schur ratio {t['schur_ratio']:.4f}"
        )
        print(f"    nested fit {t['p_nested']}, ll per step {np.array(t['trace'])}")
        for name, p in _p_values(t, t["D"], test).items():
            print(f"    p {name}: {p:.4g}")
    if res["sd"] is not None:
        for name, sd in res["sd"].items():
            print(f"  sd {name}: {sd}")
    if "p_boot" in res:
        print(f"  sd bootstrap: {res['p_boot'].std(axis=0, ddof=1)}")
    sys.stdout.flush()


def run(scenario, seed, regions, procs, boot_refits, methods):
    if SCENARIOS[scenario]["source"] == "new":
        simulate(scenario, seed, regions, procs)
    for method in methods:
        analyze(scenario, seed, method, procs, boot_refits)


# -----------------------------------------------------------------------------
# Summaries
# -----------------------------------------------------------------------------

ALPHAS = [0.01, 0.05, 0.1]


def load_results(method, scenario):
    out = []
    for path in sorted(
        glob.glob(os.path.join(RESULTS_DIR, f"{method}_{scenario}_*.pkl"))
    ):
        with open(path, "rb") as fin:
            out.append(pickle.load(fin))
    return out


def summarize():
    for method in MODELS:
        for scenario in SCENARIOS:
            R = load_results(method, scenario)
            if not R:
                continue
            print(f"\n{method.upper()}, {scenario}, {len(R)} datasets")
            for test in TESTS:
                T = [r["tests"][test] for r in R]
                kind = "size" if T[0]["null_true"] else "power"
                ps = {}
                for t in T:
                    for k, v in _p_values(t, t["D"], test).items():
                        ps.setdefault(k, []).append(v)
                me = np.array([t["min_eig"] for t in T])
                sr = np.array([t["schur_ratio"] for t in T])
                print(
                    f"  {test} ({kind}): median D {np.median([t['D'] for t in T]):.3g};"
                    f" min eig min {np.nanmin(me):.4f}; Schur ratio min"
                    f" {np.nanmin(sr):.4f}, median {np.nanmedian(sr):.4f}"
                )
                for k, v in ps.items():
                    v = np.array(v, dtype=float)
                    rates = " ".join(f"{np.nanmean(v < a):.2f}" for a in ALPHAS)
                    print(f"    {k:32s} p < {ALPHAS}: {rates}")
            with_sd = [r for r in R if r["sd"] is not None]
            if with_sd:
                p_true = R[0]["p_true"]
                print("  95% CI coverage: " + " ".join(f"{n:>7s}" for n in PARAM_NAMES))
                for name in ["FIM", "GIM"]:
                    cover = np.mean(
                        [
                            np.abs(r["p_full"] - p_true) <= 1.96 * r["sd"][name]
                            for r in with_sd
                        ],
                        axis=0,
                    )
                    print(f"    {name:14s}" + " ".join(f"{c:7.2f}" for c in cover))
                if "p_boot" in R[0]:
                    lo_hi = [np.percentile(r["p_boot"], [2.5, 97.5], axis=0) for r in R]
                    cover = np.mean(
                        [(lo <= p_true) & (p_true <= hi) for lo, hi in lo_hi], axis=0
                    )
                    print("    Bootstrap pct " + " ".join(f"{c:7.2f}" for c in cover))


# -----------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="stage", required=True)
    p_run = sub.add_parser("run", help="simulate (new scenarios) and analyze")
    p_sim = sub.add_parser("simulate", help="simulate and parse one new dataset")
    p_ana = sub.add_parser("analyze", help="fit the models and run the tests")
    sub.add_parser("summarize", help="print test calibration and CI coverage")
    for p in [p_run, p_sim, p_ana]:
        p.add_argument("--scenario", choices=SCENARIOS, required=True)
        p.add_argument("--seed", type=int, required=True)
        p.add_argument("--procs", type=int, default=1)
    for p in [p_run, p_sim]:
        p.add_argument("--regions", type=int, default=base.REGIONS)
    for p in [p_run, p_ana]:
        p.add_argument("--boot-refits", type=int, default=BOOT_REFITS)
    p_run.add_argument("--methods", nargs="+", choices=MODELS, default=list(MODELS))
    p_ana.add_argument("--method", choices=MODELS, required=True)
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
