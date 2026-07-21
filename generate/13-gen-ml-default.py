#!/usr/bin/env python3
"""Generate synthetic data using default-hyperparameter ML methods.

Methods: CTGAN, TVAE, GaussianCopula, CopulaGAN (all from SDV).
Scenarios: 1-10 simulation scenarios.
Hyperparameters: SDV defaults (CTGAN/TVAE: 300 epochs, batch 500).

Input:  sim_scenario_{01..10}_orig.csv  (from 12-export-scenarios.r)
Output: sim_scenario_{01..10}_{ctgan,tvae,gauscop,copgan}.csv

Usage:
    uv run --python .venv 13-gen-ml-default.py
"""
import sys
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

SEED = 6489
SCENARIOS = range(1, 11)   # ML baselines cover the ten continuous/count
                           # scenarios.  Scenario 11 (binomial) is compared
                           # against the projection and copula only.

METHODS = {
    "ctgan": lambda md: CTGANSynthesizer(md, epochs=300, batch_size=500,
                                          verbose=False),
    "tvae": lambda md: TVAESynthesizer(md, epochs=300, batch_size=500),
    "gauscop": lambda md: GaussianCopulaSynthesizer(md),
    "copgan": lambda md: CopulaGANSynthesizer(md, epochs=300, batch_size=500,
                                               verbose=False),
}


def generate_one(sc, method_name, make_synth):
    path_in = f"results/sim_scenario_{sc:02d}_orig.csv"
    path_out = f"results/sim_scenario_{sc:02d}_{method_name}.csv"

    try:
        d = pd.read_csv(path_in)
    except FileNotFoundError:
        print(f"  {path_in} not found, skipping")
        return

    # All three lines are required. np.random.seed reaches neither the torch RNG
    # that ctgan/tvae/copgan train under nor the fixed internal RNG that copulas
    # gives GaussianCopula, which returns one sample for every seed unless
    # _set_random_state is called after fit.
    np.random.seed(SEED)
    torch.manual_seed(SEED)

    metadata = Metadata.detect_from_dataframe(data=d)
    synth = make_synth(metadata)

    try:
        synth.fit(d)
        synth._set_random_state(SEED)   # after fit: _model exists only then
        d_syn = synth.sample(num_rows=len(d))
        d_syn.to_csv(path_out, index=False)
        print(f"  {method_name} -> {path_out} ({len(d_syn)} rows)")
    except Exception as e:
        print(f"  {method_name} FAILED: {e}", file=sys.stderr)


def main():
    for sc in SCENARIOS:
        print(f"\n=== Scenario {sc} ===")
        for name, factory in METHODS.items():
            generate_one(sc, name, factory)

    print("\nDone.")


if __name__ == "__main__":
    main()
