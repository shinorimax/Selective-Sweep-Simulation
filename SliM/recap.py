# import sys, os
# import pyslim
# import numpy as np
# import msprime as msp
# import tskit

# # ----- helpers (API-compat) -----
# def recapitate_compat(ts, Ne, r, seed):
#     # Prefer msprime.recapitate if ts.recapitate doesn't exist
#     if hasattr(ts, "recapitate"):
#         return ts.recapitate(recombination_rate=r, Ne=Ne, random_seed=seed % 2**32)
#     else:
#         return msp.recapitate(ts, Ne=Ne, recombination_rate=r, random_seed=seed % 2**32)

# def mutate_compat(ts, mu, seed):
#     # Try msprime.mutate, else sim_mutations
#     if hasattr(msp, "mutate"):
#         return msp.mutate(ts, rate=mu, random_seed=seed % 2**32)
#     else:
#         return msp.sim_mutations(ts, rate=mu, random_seed=seed % 2**32, model="infinite_sites")

# def subsample(ts, N, seed):
#     rng = np.random.default_rng(seed % 2**32)
#     samples = ts.samples()
#     if len(samples) < N:
#         raise ValueError(f"Requested N={N} but only {len(samples)} samples in ts.")
#     keep = rng.choice(samples, N, replace=False)
#     return ts.simplify(keep)

# if __name__ == "__main__":
#     # Args: simID N Ne r mu
#     simID = int(sys.argv[1])
#     N     = int(sys.argv[2])
#     Ne    = int(sys.argv[3])
#     r     = float(sys.argv[4])
#     mu    = float(sys.argv[5])

#     in_path  = f"trees/{simID}.trees"       # must match your SLiM script output
#     out_path = f"trees/r{simID}.trees"      # recap’d+mutated+subsampled output

#     os.makedirs("trees", exist_ok=True)

#     print("recap is starting")
#     # ts = pyslim.load(in_path).simplify()
#     # NEW
#     ts = tskit.load(in_path).simplify()

#     ts = recapitate_compat(ts, Ne=Ne, r=r, seed=simID)
#     print("recap is done")

#     ts = mutate_compat(ts, mu=mu, seed=simID)
#     print("mutate is done")

#     ts = subsample(ts, N=N, seed=simID)
#     print("subsample is done")

#     ts.dump(out_path)
#     print("it's done")

#     # OPTIONAL: move original to scratch; comment out if you don't have that path
#     # scratch_prefix = os.environ.get("SCRATCH_PREFIX", None)  # e.g., "/scratch/my_lab/user"
#     # if scratch_prefix:
#     #     dst_dir = os.path.join(scratch_prefix, "trees")
#     #     os.makedirs(dst_dir, exist_ok=True)
#     #     dst_path = os.path.join(dst_dir, f"{simID}.trees")
#     #     import shutil
#     #     shutil.move(in_path, dst_path)
#     #     print(f"Moved original to {dst_path}")

import sys, os
import numpy as np
import tskit
import msprime as msp
import pyslim  # keep imported; we’ll use pyslim.recapitate & helpers

# ----- helpers (API-compat) -----
def recapitate_compat(ts, Ne, r, seed):
    # Prefer pyslim.recapitate (recommended)
    if hasattr(pyslim, "recapitate"):
        # In pyslim, use ancestral_Ne (diploid) and recombination_rate
        return pyslim.recapitate(
            ts,
            recombination_rate=r,
            ancestral_Ne=Ne,
            random_seed=seed % (2**32),
        )
    # Try TreeSequence method (older stacks)
    if hasattr(ts, "recapitate"):
        return ts.recapitate(recombination_rate=r, Ne=Ne, random_seed=seed % (2**32))
    # Try top-level msprime function (some versions)
    if hasattr(msp, "recapitate"):
        return msp.recapitate(ts, Ne=Ne, recombination_rate=r, random_seed=seed % (2**32))
    # Fallback
    raise RuntimeError(
        "No recapitate API found. Install a compatible stack, e.g.: "
        "pip install 'msprime==1.2.*' 'tskit==0.5.*' 'pyslim>=1.0'"
    )

def mutate_compat(ts, mu, seed):
    # Use msprime.sim_mutations; if reloading into SLiM is possible, keep SLiM model + IDs
    next_id = None
    try:
        next_id = pyslim.next_slim_mutation_id(ts)
    except Exception:
        pass
    try:
        model = msp.SLiMMutationModel(type=0, next_id=next_id) if next_id is not None else "infinite_sites"
    except Exception:
        model = "infinite_sites"

    return msp.sim_mutations(
        ts,
        rate=mu,
        model=model,
        keep=True,  # keep any existing (SLiM) muts, if present
        random_seed=seed % (2**32),
    )

def subsample_genomes(ts, N, seed):
    # Keep N genomes (nodes) at time 0
    rng = np.random.default_rng(seed % (2**32))
    samples = ts.samples()
    if len(samples) < N:
        raise ValueError(f"Requested N={N} but only {len(samples)} sample nodes in ts.")
    keep = rng.choice(samples, N, replace=False)
    # After recap, regular simplify is fine
    return ts.simplify(keep)

def subsample_diploids(ts, N_diploids, seed):
    # Keep N_diploids individuals (=> 2*N_diploids genomes)
    rng = np.random.default_rng(seed % (2**32))
    alive_inds = pyslim.individuals_alive_at(ts, 0)
    if len(alive_inds) < N_diploids:
        raise ValueError(f"Requested {N_diploids} diploids but only {len(alive_inds)} alive.")
    keep_inds = rng.choice(alive_inds, N_diploids, replace=False)
    keep_nodes = []
    for i in keep_inds:
        keep_nodes.extend(ts.individual(i).nodes)
    # Keep input roots by default here is not necessary post-recap, but harmless
    return ts.simplify(keep_nodes, keep_input_roots=True)

if __name__ == "__main__":
    # Args: simID N Ne r mu
    simID = int(sys.argv[1])
    N     = int(sys.argv[2])      # interpret as *genomes*; see note below to switch to diploids
    Ne    = int(sys.argv[3])
    r     = float(sys.argv[4])
    mu    = float(sys.argv[5])

    in_path  = f"trees/{simID}.trees"       # SLiM raw output
    out_path = f"trees/r{simID}.trees"      # recap’d + mutated + subsampled

    os.makedirs("trees", exist_ok=True)
    print("recap is starting")

    # IMPORTANT: do NOT simplify before recapitation
    ts = tskit.load(in_path)

    # Recapitate (backwards coalescent history)
    ts = recapitate_compat(ts, Ne=Ne, r=r, seed=simID)
    print("recap is done")

    # Add neutral mutations (discrete sites by default in msprime>=1.0)
    ts = mutate_compat(ts, mu=mu, seed=simID)
    print("mutate is done")

    # Subsample:
    #   If N means genomes: use subsample_genomes(ts, N, ...)
    #   If you want N diploid individuals: use subsample_diploids(ts, N, ...)
    ts = subsample_genomes(ts, N, seed=simID)
    # ts = subsample_diploids(ts, N, seed=simID)  # <-- use this instead if N = #diploids
    print("subsample is done")

    ts.dump(out_path)
    print("it's done")
