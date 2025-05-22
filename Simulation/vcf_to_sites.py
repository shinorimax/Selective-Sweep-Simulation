from cyvcf2 import VCF
import sys

# === INPUT AND OUTPUT ===
input_vcf = "/Users/yagishinnosuke/Documents/2024-2025 Stanford/Research/Selective-Sweep-Simulation/Results/Two_Sample_Test_ARG/simulated_data.vcf"         # change this path if needed
output_sites = "/Users/yagishinnosuke/Documents/2024-2025 Stanford/Research/Selective-Sweep-Simulation/Results/Two_Sample_Test_ARG/simulated_data.sites"    # desired output path


def parse_vcf(vcf_path, output_path, chrom="chr", region_start=1, region_end=100000):
    sample_names = []
    site_data = []

    with open(vcf_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("##"):
                continue
            elif line.startswith("#CHROM"):
                parts = line.split('\t')
                sample_names = parts[9:]
            else:
                parts = line.split('\t')
                pos = parts[1]
                ref = parts[3]
                alt = parts[4]
                gts = parts[9:]

                # Skip multi-allelic or indel sites
                if len(ref) != 1 or len(alt) != 1 or ',' in alt:
                    continue

                alleles = [ref, alt]
                row = []

                for gt in gts:
                    gt_val = gt.split(":")[0]  # get only genotype (e.g., "0|1")
                    hap = gt_val.split('|')[0]  # take first haplotype
                    if hap == ".":
                        row.append("N")
                    else:
                        try:
                            row.append(alleles[int(hap)])
                        except:
                            row.append("N")
                site_data.append((pos, ''.join(row)))

    with open(output_path, "w") as out:
        out.write("NAMES\t" + "\t".join(sample_names) + "\n")
        out.write(f"REGION\t{chrom}\t{region_start}\t{region_end}\n")
        for pos, seq in site_data:
            out.write(f"{pos}\t{seq}\n")

parse_vcf(input_vcf, output_sites)
