import msprime
import json
import os
import random
import daiquiri
import logging
from tqdm import tqdm
import tskit
from IPython.display import SVG, display

def select_scenario(sc, Ne, pos):
    """Define msprime selective sweep model."""
    if sc == 0:
        return msprime.StandardCoalescent()  # Neutral evolution
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
    random.seed(123)

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
            num_trees = sim.num_trees - 1
            if m == num_trees:
                indices = [m, m, m, m]
            else:
                indices = [m, m+1, int((m + num_trees) / 2), num_trees]
            trees = extract_trees(sim, indices)
            
            # Append to scenario data
            results[scenario_key].append(trees)

    # Save updated results as JSON
    with open(file_path, "w") as f:
        json.dump(results, f, indent=2)

    return file_path


# Example usage
if __name__ == "__main__":

    # ############################
    # # Run Multiple Simulations #
    # ############################
    selection_scenario = [0.01]
    Ne = int(1e4)
    L = int(1e5)
    num_samples = 25
    num_simulations = 300
    num_sweeps = 1
    recombination_rate = 1.25 * 1e-8
    log_setting = [(1, 0), (1, 1), (1, 2)]
    
    # Simulate multiple sweeps
    # # output_file = simulate_multiple_sweeps(
    # #     selection_scenario=selection_scenario,
    # #     Ne=Ne,
    # #     L=L,
    # #     num_samples=num_samples,
    # #     num_simulations=num_simulations,
    # #     num_sweeps=num_sweeps,
    # #     recombination_rate=recombination_rate,
    # #     # log_simulations=log_setting
    # # )
    # # print(f"Results saved to {output_file}")
    
    # Simulate single sweep
    output_file = run_simple_simulations(
        selection_scenario=selection_scenario,
        Ne=Ne,
        L=L,
        num_samples=num_samples,
        num_simulations=num_simulations,
        recombination_rate=recombination_rate,
        folder="Results/Two_sample_results"
    )
    print(f"Results saved to {output_file}")


    #########################################
    # Simulate Single Instance and Visualize#
    #########################################

    # Define parameters
    # Ne = int(1e4)  # Effective population size
    # L = int(1e5)  # Chromosome length
    # num_samples = 10  # Number of samples
    # selection_scenario = 0.1  # Example selection coefficient

    # # Run the simulation
    # ts = run_simple_simulation(selection_scenario, Ne, L, num_samples)

    # # Save the visualization to an SVG file
    # with open("tree_sequence.svg", "w") as f:
    #     f.write(ts.draw_svg(y_axis=True))

    # print("SVG saved as tree_sequence.svg. Open the file to view the visualization.")
