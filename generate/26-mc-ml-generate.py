#!/usr/bin/env python3
"""Monte Carlo ML synthetic data generation (parallel).

Trains SDV generators from scratch for each (scenario, method, seed)
combination, then samples one synthetic dataset per training run.
This captures the full training variability (SGD initialisation,
minibatch ordering) that a practitioner faces.

Input:  sim_scenario_{01,03,04,06,08,09}_orig.csv
        (from 12-export-scenarios.r)
Output: results/mc/sc{01,03,04,06,08,09}_seed{seed}_{method}.csv

Usage on HPC:
    export PATH="$HOME/.local/bin:$PATH"
    .venv/bin/python 26-mc-ml-generate.py
"""

import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

# ── Configuration ──────────────────────────────────────────────

B = 100
MASTER_SEED = 6489
SCENARIOS = [1, 3, 4, 6, 8, 9]
METHODS = ["ctgan", "tvae", "gauscop", "copgan"]
MAX_WORKERS = int(os.environ.get("MC_WORKERS", "8"))
OUT_DIR = Path("results/mc")
LOG_FILE = Path("results/mc-ml-generate.log")


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
    """Train one SDV model and sample one synthetic dataset.

    Each call is a fresh process (ProcessPoolExecutor), so no
    shared state or GIL concerns.  OMP/MKL threads are pinned to 1
    to avoid oversubscription on multi-core machines.
    """
    sc, method, seed = task

    # Pin threads per worker.  LOKY_MAX_CPU_COUNT matters as much as the BLAS
    # vars: SDV pulls in joblib, and without it each of the N workers spawns its
    # own nested loky pool sized to the whole machine.
    for var in ("OMP_NUM_THREADS", "MKL_NUM_THREADS",
                "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS",
                "LOKY_MAX_CPU_COUNT"):
        os.environ[var] = "1"

    import warnings
    warnings.filterwarnings("ignore")

    import numpy as np
    import pandas as pd
    import torch
    from sdv.metadata import Metadata
    from sdv.single_table import (
        CTGANSynthesizer,
        TVAESynthesizer,
        GaussianCopulaSynthesizer,
        CopulaGANSynthesizer,
    )

    # One thread also makes the torch reductions deterministic, so a seed
    # reproduces bit for bit rather than only in distribution.
    torch.set_num_threads(1)

    # Seeding, all three lines required.  SDV runs every generator through
    # ctgan/synthesizers/base.py, which falls back to an ambient RNG whenever the
    # model's own random_states is None (its default).  np.random.seed alone
    # therefore reaches neither the torch RNG that ctgan/tvae/copgan train under
    # Nor the fixed internal RNG that copulas gives GaussianCopula.  Omitting
    # torch.manual_seed makes the neural runs unreproducible.  Omitting
    # _set_random_state below makes GaussianCopula return the same sample for
    # every seed.
    np.random.seed(seed)
    torch.manual_seed(seed)

    cls_map = {
        "ctgan":   (CTGANSynthesizer,
                    dict(epochs=300, batch_size=500, verbose=False)),
        "tvae":    (TVAESynthesizer,
                    dict(epochs=300, batch_size=500)),
        "gauscop": (GaussianCopulaSynthesizer,
                    dict()),
        "copgan":  (CopulaGANSynthesizer,
                    dict(epochs=300, batch_size=500, verbose=False)),
    }

    f_in = f"results/sim_scenario_{sc:02d}_orig.csv"
    f_out = OUT_DIR / f"sc{sc:02d}_seed{seed}_{method}.csv"

    try:
        d = pd.read_csv(f_in)
        md = Metadata.detect_from_dataframe(data=d)
        cls, kwargs = cls_map[method]
        synth = cls(md, **kwargs)
        synth.fit(d)
        synth._set_random_state(seed)   # after fit: _model exists only then
        d_syn = synth.sample(num_rows=len(d))
        d_syn.to_csv(f_out, index=False)
        return (sc, method, seed, "ok")
    except Exception as e:
        return (sc, method, seed, f"FAIL: {e}")


# ── Main ───────────────────────────────────────────────────────

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    sys.stdout = TeeWriter(sys.__stdout__, LOG_FILE)
    sys.stderr = TeeWriter(sys.__stderr__, LOG_FILE)

    # Generate per-task seeds: deterministic from master seed
    import numpy as np
    rng = np.random.RandomState(MASTER_SEED)
    all_seeds = rng.randint(1, 2**31 - 1, size=B)

    # Build task list
    tasks = [
        (sc, method, int(all_seeds[b]))
        for sc in SCENARIOS
        for method in METHODS
        for b in range(B)
    ]

    # Skip already-completed tasks (for restartability)
    remaining = [
        t for t in tasks
        if not (OUT_DIR / f"sc{t[0]:02d}_seed{t[2]}_{t[1]}.csv").exists()
    ]

    total = len(tasks)
    done = total - len(remaining)
    print(f"=== MC ML GENERATION: {total} tasks, {done} done, "
          f"{len(remaining)} remaining, {MAX_WORKERS} workers ===\n",
          flush=True)

    if not remaining:
        print("All tasks complete.", flush=True)
        return

    t0 = time.time()
    n_ok = done
    n_fail = 0

    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as pool:
        futures = {pool.submit(run_one, t): t for t in remaining}

        for future in as_completed(futures):
            sc, method, seed, status = future.result()
            if status == "ok":
                n_ok += 1
            else:
                n_fail += 1
                print(f"  FAIL sc{sc} {method} seed{seed}: {status}",
                      file=sys.stderr, flush=True)

            if (n_ok + n_fail) % 50 == 0 or (n_ok + n_fail) == total:
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
