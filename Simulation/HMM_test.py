import msprime
import numpy as np
import tskit
import os

# Simulation parameters
n_samples = 25  # Sample size
sequence_length = 1e5  # 100,000 base pairs
recomb_rate = 1e-8
mutation_rate = 1e-8
Ne = 10000
sweep_position = sequence_length / 2

# --- Neutral Simulation ---
ts_neutral = msprime.sim_ancestry(
    samples=n_samples,
    sequence_length=sequence_length,
    recombination_rate=recomb_rate,
    population_size=Ne,
    random_seed=42
)
ts_neutral = msprime.sim_mutations(ts_neutral, rate=mutation_rate, random_seed=43)

# --- Selective Sweep Simulation ---
ts_sweep = msprime.sim_ancestry(
    samples=n_samples,
    sequence_length=sequence_length,
    recombination_rate=recomb_rate,
    population_size=Ne,
    model=msprime.SweepGenicSelection(
        position=sweep_position,
        start_frequency=0.01,
        end_frequency=0.99,
        s=0.1,        # selection coefficient
        dt=1e-4       # time step for selection modeling
    ),
    random_seed=44
)
ts_sweep = msprime.sim_mutations(ts_sweep, rate=mutation_rate, random_seed=45)

# --- Save trees with genomic positions ---
def save_newick_trees_with_positions(ts, filename_prefix):
    os.makedirs("Results", exist_ok=True)
    filepath = f"Results/{filename_prefix}_trees_with_pos.txt"
    with open(filepath, "w") as f:
        for tree in ts.trees():
            left = int(tree.interval.left)
            right = int(tree.interval.right)
            for i, root in enumerate(tree.roots):
                newick = tree.as_newick(root=root)
                f.write(f"{left},{right},root{i}\t{newick}\n")


save_newick_trees_with_positions(ts_neutral, "HMM_test_neutral")
save_newick_trees_with_positions(ts_sweep, "HMM_test_sweep")

print("Tree simulation complete. Files saved as:")
print("- neutral_trees_with_pos.txt")
print("- sweep_trees_with_pos.txt")