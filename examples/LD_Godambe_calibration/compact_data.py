"""
Shrinks the dataset pickles from godambe_calibration.py (data/{scenario}_{seed}.pkl),
which are ~30 MB each, almost all in the 500 full-size (81 x 81, float64) SFS bootstrap
spectra.

The analyses only use these spectra projected to SFS_SAMPLE_SIZE haploids per deme, so
each bootstrap spectrum is replaced by its projection (41 x 41) in float32, and boot_idx
is stored as int16. Everything else is kept, so both godambe_calibration.py and
growth_asym_calibration.py read the converted files unchanged (SFSAnalysis projects the
bootstraps again, which leaves them as they are). The full-size bootstraps are lost:
analyses at a larger projection would need the data to be simulated again.

Each file is rewritten in place through a temporary file, after checking that the
converted bootstraps match the projections of the originals. Files already converted are
skipped.

    python compact_data.py data/mig_1.pkl            # one file
    python compact_data.py data/*.pkl                # all of them
    python compact_data.py --dry-run data/*.pkl      # report sizes only
"""

import argparse
import os
import pickle

import numpy as np

import moments

import godambe_calibration as base


def compact(data):
    ns = [base.SFS_SAMPLE_SIZE] * 2
    full = data["sfs_bootstraps"]
    proj = np.array(
        [moments.Spectrum(b).project(ns).data for b in full], dtype=np.float32
    )
    # check: the analysis projects these again, and must get the same spectra
    for b, p in zip(full[:5], proj[:5]):
        a = moments.Spectrum(b).project(ns)
        c = moments.Spectrum(p.astype(float)).project(ns)
        if not np.allclose(a.data[~a.mask], c.data[~c.mask], rtol=1e-6, atol=1e-6):
            raise ValueError("projected bootstraps don't match")
    data = dict(data)
    data["sfs_bootstraps"] = proj
    data["sfs_bootstraps_projected_to"] = ns
    data["boot_idx"] = np.asarray(data["boot_idx"]).astype(np.int16)
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("paths", nargs="+")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    before = after = 0
    for path in args.paths:
        with open(path, "rb") as fin:
            data = pickle.load(fin)
        size = os.path.getsize(path)
        if "sfs_bootstraps_projected_to" in data or data.get("format") == "compact":
            print(f"{path}: already compact ({size / 1e6:.1f} MB)")
            before += size
            after += size
            continue
        new = compact(data)
        blob = pickle.dumps(new)
        before += size
        after += len(blob)
        print(f"{path}: {size / 1e6:.1f} MB -> {len(blob) / 1e6:.1f} MB")
        if not args.dry_run:
            with open(path + ".tmp", "wb") as fout:
                fout.write(blob)
            os.replace(path + ".tmp", path)
    print(f"total: {before / 1e9:.2f} GB -> {after / 1e9:.2f} GB")


if __name__ == "__main__":
    main()
