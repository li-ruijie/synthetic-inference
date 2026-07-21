#!/usr/bin/env python3
"""Generate synthetic data using tuned CTGAN hyperparameters.

Tests whether hyperparameter tuning resolves the structural failures
predicted by Proposition (no-guarantee). Tuning: more epochs (1000),
smaller batch (100), deeper generator and discriminator (three hidden
layers of 256 against the default two), lower learning rate.

Input:  sim_scenario_{01,03,04,06}_orig.csv
Output: sim_scenario_{01,03,04,06}_ctgan_tuned.csv

Usage:
    uv run --python .venv 23-gen-ml-tuned.py
"""
import sys
import numpy as np
import pandas as pd
import torch
from sdv.metadata import Metadata
from sdv.single_table import CTGANSynthesizer

SEED = 6489
SCENARIOS = [1, 3, 4, 6]


def main():
    for sc in SCENARIOS:
        path_in = f"results/sim_scenario_{sc:02d}_orig.csv"
        path_out = f"results/sim_scenario_{sc:02d}_ctgan_tuned.csv"

        try:
            d = pd.read_csv(path_in)
        except FileNotFoundError:
            print(f"Scenario {sc}: {path_in} not found, skipping")
            continue

        print(f"\n=== Scenario {sc} (n={len(d)}) ===")

        # CTGAN trains
        # under the torch RNG, which np.random.seed does not reach.
        np.random.seed(SEED)
        torch.manual_seed(SEED)

        metadata = Metadata.detect_from_dataframe(data=d)

        synth = CTGANSynthesizer(
            metadata,
            epochs=1000,
            batch_size=100,
            generator_dim=(256, 256, 256),
            discriminator_dim=(256, 256, 256),
            generator_lr=1e-4,
            discriminator_lr=1e-4,
            verbose=True,
        )

        try:
            synth.fit(d)
            synth._set_random_state(SEED)   # after fit: _model exists only then
            d_syn = synth.sample(num_rows=len(d))
            d_syn.to_csv(path_out, index=False)
            print(f"  Saved to {path_out} ({len(d_syn)} rows)")
        except Exception as e:
            print(f"  FAILED: {e}", file=sys.stderr)

    print("\nDone.")


if __name__ == "__main__":
    main()
