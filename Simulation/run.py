import msprime
import json
import os
import random
import daiquiri
import logging
from tqdm import tqdm

def select_scenario(sc, Ne, pos):
    """Define msprime selective sweep model."""
    return msprime.SweepGenicSelection(
        position=pos,
        start_frequency=1.0 / (2 * Ne),
        end_frequency=1.0 - (1.0 / (2 * Ne)),
        s=sc,
        dt=1e-6
    )

def extract_trees(sim, indices):
    """Extract trees at specific indices and return Newick format strings."""
    return [sim.at_index(idx).as_newick() for idx in indices]

def get_file_path(Ne, L, num_samples, folder="Results"):
    """Generate a dynamic file name for the output JSON file."""
    file_name = f"results_Ne_{Ne}_L_{L}_samples_{num_samples}.json"
    os.makedirs(folder, exist_ok=True)
    return os.path.join(folder, file_name)

def run_simple_simulations(selection_scenario, Ne, L, num_samples, num_simulations, recombination_rate=1e-7, folder="results"):
    """
    Run msprime simulations for multiple selection scenarios and save the results.
    
    Parameters:
        selection_scenario (list): List of selection coefficients (e.g., [0.01, 0.5, 1.0]).
        Ne (int): Effective population size.
        L (int): Sequence length.
        num_samples (int): Number of samples in each simulation.
        num_simulations (int): Number of simulations per selection coefficient.
        recombination_rate (float): Recombination rate per base pair.
        folder (str): Folder to save the results.

    Returns:
        str: Path to the output JSON file.
    """
    # Create the file path
    file_path = get_file_path(Ne, L, num_samples, folder=folder)

    # Load existing data if file exists
    if os.path.exists(file_path):
        with open(file_path, "r") as f:
            results = json.load(f)
    else:
        results = {}

    # Run simulations
    for sc in selection_scenario:
        scenario_key = f"selection_scenario_{sc}"
        if scenario_key not in results:
            results[scenario_key] = []

        for _ in tqdm(range(num_simulations), colour="green", desc=f"Scenario {sc}"):
            mut_loc = int(L / 4)
            sim = msprime.sim_ancestry(
                samples=num_samples,
                model=[select_scenario(sc, Ne, mut_loc), msprime.StandardCoalescent()],
                population_size=Ne,
                recombination_rate=recombination_rate,
                sequence_length=L,
            )
            m = sim.at(mut_loc).index
            num_trees = sim.num_trees
            indices = [m, m+1, int((m + num_trees) / 2), sim.num_trees - 1]
            trees = extract_trees(sim, indices)
            
            # Append to scenario data
            results[scenario_key].append(trees)

    # Save updated results as JSON
    with open(file_path, "w") as f:
        json.dump(results, f, indent=2)

    return file_path

def simulate_multiple_sweeps(selection_scenario, Ne, L, num_samples, num_simulations, num_sweeps, recombination_rate=1e-7, folder="results", log_simulations=None):
    """
    Run msprime simulations with multiple sweeps for various selection scenarios and save the results.
    
    Parameters:
        selection_scenario (list): List of selection coefficients (e.g., [0.01, 0.5, 1.0]).
        Ne (int): Effective population size.
        L (int): Sequence length.
        num_samples (int): Number of samples in each simulation.
        num_simulations (int): Number of simulations per selection coefficient.
        num_sweeps (int): Number of selective sweeps to simulate.
        recombination_rate (float): Recombination rate per base pair.
        folder (str): Folder to save the results.
        log_simulations (list of tuples): Specific scenarios and simulation indices to log, e.g., [(1, 0), (2, 5)].

    Returns:
        str: Path to the output JSON file.
    """
    import logging
    
    # random.seed(123)  # Seed for reproducibility
    # Create the file path
    file_path = get_file_path_multiple_sweeps(Ne, L, num_samples, folder=folder)

    # Load existing data if file exists
    if os.path.exists(file_path):
        with open(file_path, "r") as f:
            results = json.load(f)
    else:
        results = {}

    # Configure logging
    log_simulations = log_simulations or []
    logger = daiquiri.getLogger("msprime_simulation")

    # Run simulations
    for scenario_idx, sc in enumerate(selection_scenario):
        scenario_key = f"selection_scenario_{sc}_sweeps_{num_sweeps}"
        if scenario_key not in results:
            results[scenario_key] = []

        for sim_idx in tqdm(range(num_simulations), colour="green", desc=f"Scenario {sc} with {num_sweeps} sweeps"):
            # Enable logging for specific simulations
            if (scenario_idx, sim_idx) in log_simulations:
                daiquiri.setup(level="INFO")
            else:
                daiquiri.setup(level="WARNING")  # Suppress detailed logs

            # Generate models with multiple sweeps
            mut_loc = int(L / 4)
            models = []
            for _ in range(num_sweeps):
                models.append(msprime.StandardCoalescent(duration=random.uniform(0, 100)))  # Random coalescent
                models.append(
                    msprime.SweepGenicSelection(
                        position=mut_loc,
                        start_frequency=1 / (2 * Ne),
                        end_frequency=1.0 - 1 / (2 * Ne),
                        s=sc,
                        dt=1e-6,
                    )
                )
            # Add a final coalescent to ensure coalescence
            models.append(msprime.StandardCoalescent())

            # Simulate ancestry
            sim = msprime.sim_ancestry(
                samples=num_samples,
                model=models,
                population_size=Ne,
                recombination_rate=recombination_rate,
                sequence_length=L,
            )
            
            # Extract trees
            m = sim.at(mut_loc).index
            num_trees = sim.num_trees
            if m >= num_trees - 1:
                pass
            else:
                indices = [m, m + 1, int((m + num_trees) / 2), num_trees - 1]
                trees = extract_trees(sim, indices)
            
                # Append to scenario data
                results[scenario_key].append(trees)

    # Save updated results as JSON
    with open(file_path, "w") as f:
        json.dump(results, f, indent=2)

    return file_path


# Utility Functions
def get_file_path_multiple_sweeps(Ne, L, num_samples, folder="results"):
    """Generate a dynamic file name for the output JSON file."""
    file_name = f"results_Ne_{Ne}_L_{L}_samples_{num_samples}_sweeps.json"
    os.makedirs(folder, exist_ok=True)
    return os.path.join(folder, file_name)

# def extract_trees(sim, indices):
#     """Extract trees at specific indices and return Newick format strings."""
#     return [sim.at_index(idx).as_newick() for idx in indices]

# Example usage
if __name__ == "__main__":
    selection_scenario = [0.001, 0.01, 0.1]
    Ne = int(1e4)
    L = int(1e5)
    num_samples = 25
    num_simulations = 300
    num_sweeps = 1
    recombination_rate = 1.25 * 1e-8
    log_setting = [(1, 0), (1, 1), (1, 2)]
    
    # output_file = simulate_multiple_sweeps(
    #     selection_scenario=selection_scenario,
    #     Ne=Ne,
    #     L=L,
    #     num_samples=num_samples,
    #     num_simulations=num_simulations,
    #     num_sweeps=num_sweeps,
    #     recombination_rate=recombination_rate,
    #     # log_simulations=log_setting
    # )
    # print(f"Results saved to {output_file}")
    
    output_file = run_simple_simulations(
        selection_scenario=selection_scenario,
        Ne=Ne,
        L=L,
        num_samples=num_samples,
        num_simulations=num_simulations,
        recombination_rate=recombination_rate
    )
    print(f"Results saved to {output_file}")
