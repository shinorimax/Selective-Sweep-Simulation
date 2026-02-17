import argweaver.smc as smc

outdir = "argweaver_chr2_130_140Mb"

for i in range(0, 201, 10):
    infile  = f"{outdir}/ceu.{i}.smc.gz"
    outfile = f"{outdir}/ceu.{i}.trees"
    try:
        ts = smc.smc2tskit(infile)
        ts.dump(outfile)
        print(f"wrote {outfile}  ({ts.num_trees} trees)")
    except Exception as e:
        print(f"failed {infile}: {e}")