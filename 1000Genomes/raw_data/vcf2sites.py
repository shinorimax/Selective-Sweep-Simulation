#!/usr/bin/env python3
# vcf2sites.py — convert a phased bcftools-subsetted VCF to ARGweaver SITES format
import sys, gzip

vcf_file = sys.argv[1]   # input .vcf or .vcf.gz
chrom    = sys.argv[2]   # e.g. "2"
start    = int(sys.argv[3])
end      = int(sys.argv[4])

opener = gzip.open if vcf_file.endswith(".gz") else open

with opener(vcf_file, "rt") as fh:
    samples = []
    hap_names = []
    rows = []

    for line in fh:
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            samples = line.strip().split("\t")[9:]
            hap_names = [f"{s}_1" for s in samples] + [f"{s}_2" for s in samples]
            # we'll interleave haplotypes below
            hap_names = []
            for s in samples:
                hap_names += [f"{s}_1", f"{s}_2"]
            continue

        fields = line.strip().split("\t")
        pos = int(fields[1])
        ref, alt = fields[3], fields[4]

        # biallelic SNPs only
        if len(ref) != 1 or len(alt) != 1 or "," in alt:
            continue
        if pos < start or pos > end:
            continue

        gts = []
        for g in fields[9:]:
            gt = g.split(":")[0]
            if "|" not in gt:
                continue  # skip unphased
            a1, a2 = gt.split("|")
            gts.append(ref if a1 == "0" else alt)
            gts.append(ref if a2 == "0" else alt)

        if len(gts) != len(hap_names):
            continue
        if len(set(gts)) == 1:
            continue  # invariant after subsetting

        rows.append((pos, "".join(gts)))

# write SITES format
print("NAMES\t" + "\t".join(hap_names))
print(f"REGION\t{chrom}\t{start}\t{end}")
for pos, alleles in rows:
    print(f"{pos}\t{alleles}")
    