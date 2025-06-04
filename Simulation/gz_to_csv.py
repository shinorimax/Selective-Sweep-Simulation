import gzip
import csv
import os

def parse_smc_file(file_path, output_csv_path):
    trees = []
    tree_index = 1

    with gzip.open(file_path, "rt") as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("TREE"):
            parts = line.split(maxsplit=3)
            if len(parts) < 4:
                i += 1
                continue  # skip malformed
            left = int(parts[1])
            right = int(parts[2])
            newick = parts[3]

            # Continue to next lines if the Newick is multiline
            while not newick.strip().endswith(";") and i + 1 < len(lines):
                i += 1
                newick += lines[i].strip()

            trees.append([tree_index, left, right, newick, "unknown"])
            tree_index += 1
        i += 1

    # Write CSV
    with open(output_csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["tree_index", "left", "right", "newick", "state"])
        writer.writerows(trees)

# Run the script for all .smc.gz files
for i in range(0, 1001, 10):
    smc_path = f"/Users/yagishinnosuke/Documents/2024-2025 Stanford/Research/Selective-Sweep-Simulation/Results/Two_Sample_Test_ARG_0.1/run2.{i}.smc.gz"
    csv_path = f"/Users/yagishinnosuke/Documents/2024-2025 Stanford/Research/Selective-Sweep-Simulation/Results/Two_Sample_Test_ARG_0.1/CSVs_run2/trees_run2.{i}.csv"
    if os.path.exists(smc_path):
        print(f"Processing {smc_path} → {csv_path}")
        parse_smc_file(smc_path, csv_path)
    else:
        print(f"File not found: {smc_path}")
