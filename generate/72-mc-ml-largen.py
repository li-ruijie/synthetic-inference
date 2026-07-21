#!/usr/bin/env python3
"""Phase 2 of the larger-n ML robustness check (see 72-mc-ml-largen-report.r).

Trains the four SDV generators on the scenario-1 draws written by
72-mc-ml-largen-gen.r, at n in {2000, 20000} over 3 reps (seeds 6489+k), and
samples n rows each.  Seeded exactly as 13-gen-ml-default.py.  Work dir is
$LARGEN_DIR (default largen/).

Usage: uv run --python .venv 72-mc-ml-largen.py
"""
import os, time
import pandas as pd, numpy as np, torch
from sdv.metadata import Metadata
from sdv.single_table import (CTGANSynthesizer, TVAESynthesizer,
                              GaussianCopulaSynthesizer, CopulaGANSynthesizer)

WD = os.environ.get("LARGEN_DIR", "results/largen")
NS = [2000, 20000]; REPS = [1, 2, 3]

def build(md):
    return {"gauscop": GaussianCopulaSynthesizer(md),
            "tvae":    TVAESynthesizer(md, epochs=300, batch_size=500),
            "ctgan":   CTGANSynthesizer(md, epochs=300, batch_size=500, verbose=False),
            "copgan":  CopulaGANSynthesizer(md, epochs=300, batch_size=500, verbose=False)}

for n in NS:
    for k in REPS:
        seed = 6489 + k
        d = pd.read_csv(os.path.join(WD, f"largen_sc01_n{n}_rep{k}_orig.csv"))
        md = Metadata.detect_from_dataframe(data=d)
        for name, synth in build(md).items():
            t = time.time()
            np.random.seed(seed); torch.manual_seed(seed)
            try:
                synth.fit(d)
                synth._set_random_state(seed)
                synth.sample(num_rows=len(d)).to_csv(
                    os.path.join(WD, f"largen_sc01_n{n}_rep{k}_{name}.csv"), index=False)
                print(f"n={n:6d} rep{k} {name:8s} {time.time()-t:7.1f}s", flush=True)
            except Exception as e:
                print(f"n={n:6d} rep{k} {name:8s} FAILED: {e}", flush=True)
print("ml done.", flush=True)
