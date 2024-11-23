import msprime
from IPython.display import SVG

# Set up selective sweep model with parameters
def select_scenario(sc, Ne, pos):
    sweep_model = msprime.SweepGenicSelection(
        position=pos,
        start_frequency=1.0 / (2 * Ne),
        end_frequency=1.0 - (1.0 / (2 * Ne)),
        s=sc,
        dt=1e-6
    )
    return sweep_model

# Parameters for a small sample size simulation
selection_scenario = 0.5   # Selection coefficient for the example
Ne = 1e3                    # Effective population size
L = 1e6                     # Chromosome length
sample_size = 4            # Small sample size for visualization
mut_loc = int(L / 4)        # Mutation location at 1/4 of chromosome length

# Run a single simulation with selective sweep and coalescent model
sim = msprime.sim_ancestry(
    samples=sample_size,
    model=[select_scenario(selection_scenario, Ne, mut_loc), msprime.StandardCoalescent()],
    population_size=Ne,
    recombination_rate=1e-7,
    sequence_length=L,
)

# Visualize a single tree in the simulation
tree_index = mut_loc
tree = sim.at(tree_index)

# Save SVG to file if not in Jupyter
with open("genealogy_tree.svg", "w") as f:
    f.write(tree.draw_svg(y_axis=True))