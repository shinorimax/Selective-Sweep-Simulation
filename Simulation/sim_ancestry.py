import msprime
import os

# === Parameters ===
Ne = 10_000
L = int(1e5)
num_individuals = 25            # 25 diploid individuals = 50 haploid genomes
selection_s = 0.1
recomb_rate = 1.25e-8
mut_rate = 1e-8
mut_loc = int(L / 4)            # Selected site at 25,000

# === Define selective sweep model ===
sweep = msprime.SweepGenicSelection(
    position=mut_loc,
    start_frequency=1 / (2 * Ne),
    end_frequency=1.0 - 1 / (2 * Ne),
    s=selection_s,
    dt=1e-6,
)

# === Step 1: Simulate ancestry with diploid individuals ===
ts = msprime.sim_ancestry(
    samples=num_individuals,
    sequence_length=L,
    recombination_rate=recomb_rate,
    population_size=Ne,
    model=[sweep, msprime.StandardCoalescent()],
    random_seed=42,
)

# === Step 2: Simulate mutations ===
mutated_ts = msprime.sim_mutations(ts, rate=mut_rate, random_seed=43)

# === Step 3: Write VCF ===
output_dir = "Results/Two_Sample_Test_ARG_0.1"
os.makedirs(output_dir, exist_ok=True)

with open(os.path.join(output_dir, "simulated_data.vcf"), "w") as vcf_file:
    mutated_ts.write_vcf(vcf_file)
