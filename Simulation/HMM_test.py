import msprime
import os
import pandas as pd
from tqdm import tqdm

# Parameters
Ne = 10_000
L = int(1e5)
num_samples = 25
num_simulations = 300
selection_s = 0.01
recomb_rate = 1.25e-8
mut_loc = int(L / 4)

output_dir = "Results/SelectedTreesWithPositions"
os.makedirs(output_dir, exist_ok=True)

def select_scenario(s, Ne, pos):
    return msprime.SweepGenicSelection(
        position=pos,
        start_frequency=1.0 / (2 * Ne),
        end_frequency=1.0 - (1.0 / (2 * Ne)),
        s=s,
        dt=1e-6
    )

for i in tqdm(range(num_simulations), desc="Simulating + saving with positions"):
    seed = 1000 + i
    model = [select_scenario(selection_s, Ne, mut_loc), msprime.StandardCoalescent()]
    
    ts = msprime.sim_ancestry(
        samples=num_samples,
        model=model,
        population_size=Ne,
        sequence_length=L,
        recombination_rate=recomb_rate,
        random_seed=seed,
    )

    records = []
    for j, tree in enumerate(ts.trees()):
        newick = tree.as_newick()
        left = tree.interval.left
        right = tree.interval.right
        records.append({
            "tree_index": j + 1,
            "left": left,
            "right": right,
            "newick": newick
        })

    df = pd.DataFrame(records)
    df.to_csv(os.path.join(output_dir, f"chrom_{i+1}.csv"), index=False)
