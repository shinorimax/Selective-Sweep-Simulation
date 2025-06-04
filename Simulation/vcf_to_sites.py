import os

# === INPUT AND OUTPUT ===
input_vcf = "/Users/yagishinnosuke/Documents/2024-2025 Stanford/Research/Selective-Sweep-Simulation/Results/Two_Sample_Test_ARG_0.1/simulated_data.vcf"
output_sites = "/Users/yagishinnosuke/Documents/2024-2025 Stanford/Research/Selective-Sweep-Simulation/Results/Two_Sample_Test_ARG_0.1/simulated_data.sites"


def parse_diploid_vcf(vcf_path, output_path, chrom="chr", region_start=1, region_end=100000):
    sample_names = []
    haplotype_names = []
    site_data = []

    with open(vcf_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("##"):
                continue
            elif line.startswith("#CHROM"):
                parts = line.split('\t')
                sample_names = parts[9:]
                # Create names for haplotypes: tsk_0_a, tsk_0_b, ...
                haplotype_names = [f"{s}_a" for s in sample_names] + [f"{s}_b" for s in sample_names]
            else:
                parts = line.split('\t')
                pos = parts[1]
                ref = parts[3]
                alt = parts[4]
                gts = parts[9:]

                # Skip non-biallelic SNPs
                if len(ref) != 1 or len(alt) != 1 or ',' in alt:
                    continue

                alleles = [ref, alt]
                hap_row = []

                for gt in gts:
                    gt_val = gt.split(":")[0]  # Only take genotype field, e.g., "0|1"
                    if '|' in gt_val:
                        hap1, hap2 = gt_val.split('|')
                    elif '/' in gt_val:
                        hap1, hap2 = gt_val.split('/')
                    else:
                        hap1 = hap2 = '.'

                    # First haplotype
                    if hap1 == ".":
                        hap_row.append("N")
                    else:
                        try:
                            hap_row.append(alleles[int(hap1)])
                        except:
                            hap_row.append("N")

                    # Second haplotype
                    if hap2 == ".":
                        hap_row.append("N")
                    else:
                        try:
                            hap_row.append(alleles[int(hap2)])
                        except:
                            hap_row.append("N")

                site_data.append((pos, ''.join(hap_row)))

    with open(output_path, "w") as out:
        out.write("NAMES\t" + "\t".join(haplotype_names) + "\n")
        out.write(f"REGION\t{chrom}\t{region_start}\t{region_end}\n")
        for pos, seq in site_data:
            out.write(f"{pos}\t{seq}\n")

# Run the function
parse_diploid_vcf(input_vcf, output_sites)
