#!/usr/bin/env python3
"""Monte Carlo DP-ML synthetic data generation (parallel).

Trains differentially private ML generators from scratch for each
(scenario, method, seed) combination, then samples one synthetic
dataset per training run. Companion to 26-mc-ml-generate.py, which
covers the non-private generators. Supports the DP machine
learning comparison sentence in the main text's DP simulation
subsection (route tables in supplement S6): DP-SGD CTGAN and
PATE-CTGAN coefficient deviations at eps = 10 on the same six
scenarios as the MC all-method run.

Methods (smartnoise-synth):
    dpctgan   -- CTGAN trained with DP-SGD (opacus accountant)
    patectgan -- CTGAN trained under the PATE framework
                 [jordon-pate-gan-2019 applies PATE to GANs]

Privacy: total epsilon = 10.0 per synthesizer, of which
preprocessor_eps = 1.0 is spent learning continuous-column bounds.
delta is the smartnoise default 1/(n * sqrt(n)), about 1.1e-5 at
n = 2000 (recorded in the log).

Input:  sim_scenario_{01,03,04,06,08,09}_orig.csv
        (from 12-export-scenarios.r)
Output: results/mc-dpml/sc{01,03,04,06,08,09}_seed{seed}_{method}.csv

Environment (Python 3.13):
    uv venv .venv-dpml --python 3.13
    UV_HTTP_TIMEOUT=300 uv pip install smartnoise-synth pandas \
        packaging setuptools
    Unpinned resolve gives smartnoise-synth 1.0.8,
    opacus 0.14.0, torch 2.13.0+cu130, numpy 2.5.1. On a flaky
    connection fetch the large wheels with wget -c into a directory
    and add --find-links, verifying each wheel with
    python -m zipfile -t (a silent partial poisons resolution).
    GPU: verified with a full DPCTGAN and PATECTGAN fit on an
    RTX 6000 Ada (native cu13 kernels).
    batch_size 500 requires n >= 500 (opacus sample_rate = batch/n).

Usage on HPC:
    export PATH="$HOME/.local/bin:$PATH"
    MC_WORKERS=12 nohup .venv-dpml/bin/python 56-mc-dpml-generate.py &
"""

import multiprocessing
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

# ── Configuration ──────────────────────────────────────────────

B = 25
MASTER_SEED = 6489
SCENARIOS = [1, 3, 4, 6, 8, 9]
METHODS = ["dpctgan", "patectgan"]
EPSILON = 10.0
PREPROC_EPS = 1.0
MAX_WORKERS = int(os.environ.get("MC_WORKERS", "8"))
OUT_DIR = Path("results/mc-dpml")
LOG_FILE = Path("results/mc-dpml-generate.log")

CATEGORICAL = ["female", "race"]
CONTINUOUS = ["response", "age", "income"]


# ── Logging ────────────────────────────────────────────────────

class TeeWriter:
    """Write to both a file and the original stream."""
    def __init__(self, stream, path):
        self.stream = stream
        self.file = open(path, "a", buffering=1)

    def write(self, msg):
        try:
            self.stream.write(msg)
        except Exception:
            pass
        try:
            self.file.write(msg)
        except Exception:
            pass

    def flush(self):
        try:
            self.stream.flush()
        except Exception:
            pass
        try:
            self.file.flush()
        except Exception:
            pass


# ── Worker function ────────────────────────────────────────────

def run_one(task):
    """Train one DP synthesizer and sample one synthetic dataset.

    Each call is a fresh process (ProcessPoolExecutor), so no shared
    state or GIL concerns. OMP/MKL threads are pinned to 1 to avoid
    oversubscription; torch is seeded per task for reproducibility.
    """
    sc, method, seed = task

    for var in ("OMP_NUM_THREADS", "MKL_NUM_THREADS",
                "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[var] = "1"

    import warnings
    warnings.filterwarnings("ignore")

    import random
    import numpy as np
    import torch
    from snsynth import Synthesizer
    import pandas as pd

    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    torch.set_num_threads(1)

    f_in = f"results/sim_scenario_{sc:02d}_orig.csv"
    f_out = OUT_DIR / f"sc{sc:02d}_seed{seed}_{method}.csv"

    try:
        d = pd.read_csv(f_in)
        synth = Synthesizer.create(method, epsilon=EPSILON, verbose=False)
        synth.fit(d,
                  categorical_columns=CATEGORICAL,
                  continuous_columns=CONTINUOUS,
                  preprocessor_eps=PREPROC_EPS)
        d_syn = synth.sample(len(d))
        d_syn = d_syn[d.columns]
        d_syn.to_csv(f_out, index=False)
        return (sc, method, seed, "ok")
    except Exception as e:
        return (sc, method, seed, f"FAIL: {type(e).__name__}: {e}")


# ── Main ───────────────────────────────────────────────────────

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    sys.stdout = TeeWriter(sys.__stdout__, LOG_FILE)
    sys.stderr = TeeWriter(sys.__stderr__, LOG_FILE)

    # Per-rep seeds: deterministic from master seed, same construction
    # as 26-mc-ml-generate.py so reps pair across the two pipelines.
    import numpy as np
    rng = np.random.RandomState(MASTER_SEED)
    all_seeds = rng.randint(1, 2**31 - 1, size=B)

    n = 2000
    delta = 1.0 / (n * n**0.5)
    import torch
    print(f"=== MC DP-ML GENERATION: eps={EPSILON} "
          f"(preproc {PREPROC_EPS}), delta~{delta:.2e} (smartnoise "
          f"default), B={B} ===", flush=True)
    print(f"torch {torch.__version__}, cuda available: "
          f"{torch.cuda.is_available()}"
          + (f" ({torch.cuda.get_device_name(0)})"
             if torch.cuda.is_available() else ""), flush=True)

    tasks = [
        (sc, method, int(all_seeds[b]))
        for sc in SCENARIOS
        for method in METHODS
        for b in range(B)
    ]

    remaining = [
        t for t in tasks
        if not (OUT_DIR / f"sc{t[0]:02d}_seed{t[2]}_{t[1]}.csv").exists()
    ]

    total = len(tasks)
    done = total - len(remaining)
    print(f"{total} tasks, {done} done, {len(remaining)} remaining, "
          f"{MAX_WORKERS} workers\n", flush=True)

    if not remaining:
        print("All tasks complete.", flush=True)
        return

    t0 = time.time()
    n_ok = done
    n_fail = 0

    # CUDA cannot re-initialize in forked children.  Spawn is mandatory
    # for GPU workers (the parent touches CUDA for the banner above).
    with ProcessPoolExecutor(
            max_workers=MAX_WORKERS,
            mp_context=multiprocessing.get_context("spawn")) as pool:
        futures = {pool.submit(run_one, t): t for t in remaining}

        for future in as_completed(futures):
            sc, method, seed, status = future.result()
            if status == "ok":
                n_ok += 1
            else:
                n_fail += 1
                print(f"  FAIL sc{sc} {method} seed{seed}: {status}",
                      file=sys.stderr, flush=True)

            if (n_ok + n_fail) % 20 == 0 or (n_ok + n_fail) == total:
                elapsed = time.time() - t0
                rate = (n_ok + n_fail - done) / max(elapsed, 1)
                eta = (len(remaining) - (n_ok + n_fail - done)) / max(rate, 0.01)
                print(f"  [{n_ok + n_fail}/{total}] "
                      f"{elapsed:.0f}s elapsed, ~{eta:.0f}s remaining",
                      flush=True)

    elapsed = time.time() - t0
    print(f"\nDone: {n_ok} ok, {n_fail} failed, {elapsed:.0f}s total",
          flush=True)


if __name__ == "__main__":
    main()
