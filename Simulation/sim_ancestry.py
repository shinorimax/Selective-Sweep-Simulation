import msprime
import os

# Parameters
Ne = 10_000
L = int(1e5)
num_samples = 25
selection_s = 0.01
recomb_rate = 1.25e-8
mut_rate = 1e-8
mut_loc = int(L / 4)  # 25_000

# Setup selection at one site
sweep = msprime.SweepGenicSelection(
    position=mut_loc,
    start_frequency=1 / (2 * Ne),
    end_frequency=1.0 - 1 / (2 * Ne),
    s=selection_s,
    dt=1e-6,
)

# Step 1: Simulate ancestry with selection
ts = msprime.sim_ancestry(
    samples=num_samples,
    sequence_length=L,
    recombination_rate=recomb_rate,
    population_size=Ne,
    model=[sweep],
    random_seed=42,
)

# Step 2: Simulate mutations with explicit mutation rate
mutated_ts = msprime.sim_mutations(ts, rate=mut_rate, random_seed=43)

# Export to VCF format for ARGweaver
os.makedirs("Results/Two_Sample_Test_ARG", exist_ok=True)
with open("Results/Two_Sample_Test_ARG/simulated_data.vcf", "w") as vcf_file:
    mutated_ts.write_vcf(vcf_file)
    # Export to .sites format for ARGweaver
with open("Results/Two_Sample_Test_ARG/simulated_data.sites", "w") as f:
    print("NAMES", " ".join(mutated_ts.individual(i).metadata.get("name", f"tsk_{i}")
                           for i in range(mutated_ts.num_samples)), file=f)
    print("REGION 0", int(mutated_ts.sequence_length), file=f)
    print("SITES", file=f)

    for var in mutated_ts.variants():
        alleles = var.alleles
        assert len(alleles) == 2  # Ensure biallelic
        print(var.site.position, *var.genotypes, sep="\t", file=f)
