import msprime
import json
import os
import random


def get_file_path(selection_scenario, Ne, L, num_samples, folder="Results"):
    """Generate a dynamic file name for the output JSON file."""
    file_name = f"trees_newick_Ne_{Ne}_L_{L}_samples_{num_samples}_s_{selection_scenario}.json"
    os.makedirs(folder, exist_ok=True)
    return os.path.join(folder, file_name)


# Function to simulate a chromosome and save trees
def run_simple_simulation(selection_scenario, Ne, L, num_samples, recombination_rate=1.0 * 1e-7, seed = None):
    """Run msprime simulation and save trees in Newick format.""" 

    if selection_scenario==0:
        ts = msprime.sim_ancestry(
        samples=num_samples,
        model=msprime.StandardCoalescent(),
        population_size=Ne,
        recombination_rate=recombination_rate,
        sequence_length=L,
        random_seed=seed, 
        )
    else:
        ts = msprime.sim_ancestry(
            samples=num_samples,
            model=[msprime.SweepGenicSelection(position=L//4, start_frequency=1/(2*Ne), end_frequency=1-1/(2*Ne), s=selection_scenario, dt=1e-6),
                msprime.StandardCoalescent()],
            population_size=Ne,
            recombination_rate=recombination_rate,
            sequence_length=L,
            random_seed=seed, 
        )

    # Extract trees and their positions
    trees_data = [{"pos": tree.interval.left, "newick": tree.as_newick()} for tree in ts.trees()]

    output_file = get_file_path(selection_scenario, Ne, L, num_samples)

    # Save to JSON
    with open(output_file, "w") as f:
        json.dump(trees_data, f, indent=2)

    print(f"Saved tree sequence data to {output_file}")


if __name__ == "__main__":
    # Define parameters
    Ne = int(1e4)  # Effective population size
    L = int(1e5)   # Chromosome length
    num_samples = 10  # Number of samples
    selection_scenario = 0  # Selection coefficient
    seed = 42

    # Run and save simulation
    run_simple_simulation(selection_scenario, Ne, L, num_samples, seed=seed)
