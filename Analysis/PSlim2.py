import os, re, sys, json, subprocess, time
import numpy as np
import pandas as pd
import tskit
import subprocess
import shutil

# ---- optional deps used for ROC plotting ----
try:
    from sklearn.metrics import roc_curve, roc_auc_score
    _HAVE_SKLEARN = True
except Exception:
    _HAVE_SKLEARN = False

import sys
# add the parent of the *inner* bim/ to sys.path
sys.path.insert(0, "/opt/anaconda3/envs/sweep312/bin/bim")
from bim.utils import InferEta

# ---- optional: demography inference helper (only if you actually train eta) ----
# If you won't use train_eta(), you can ignore this import and leave it commented.
# sys.path.append('/home/enes/bim/')
# from utils import InferEta


# ----------------------------- utilities -----------------------------

def _ensure_dirs(*paths):
    for p in paths:
        os.makedirs(p, exist_ok=True)

def get_simIDs(cID, setID, nrep):
    """
    Deterministic simIDs & filename regex, matching the original.
    cID: class ID (0..9), setID: index of model within class (0..99), nrep: replicates
    """
    if any([cID > 9, setID > 99, nrep > 10000]):
        raise ValueError('Wrong IDs')
    start = cID * 1_000_000 + setID * 10_000
    ids = np.arange(start, start + nrep)

    setID_str = str(setID)
    cID_str = str(cID)
    if len(setID_str) == 1:
        setID_str = '0' + setID_str
    regex = 'r' + cID_str + setID_str + '[0-9]{4}.trees'
    return ids, regex

def treepaths(regex, folder='trees'):
    r = re.compile(regex)
    files = os.listdir(folder) if os.path.isdir(folder) else []
    tpaths = []
    for f in files:
        m = r.match(f)
        if m is not None:
            tpaths.append(os.path.join(folder, m.group()))
    return tpaths

def get_eps(Ne, gr, t, constant_gen=0):
    dipNe = 2 * Ne
    rep = gr + 1
    epshist = []
    for i in range(t):
        dipNe1 = dipNe * rep if i > constant_gen else dipNe
        epshist.append(dipNe)
        dipNe = dipNe1
    return epshist

def ROC(ax, y0, y1, score_ascending=False, label=None):
    """
    y0: scores for class 0, y1: scores for class 1
    """
    if not _HAVE_SKLEARN:
        raise ImportError("scikit-learn is required for ROC plotting.")
    if not score_ascending:
        y0 = -y0
        y1 = -y1
    leny1 = len(y1)
    leny0 = len(y0)
    y_true = np.r_[np.zeros(leny0), np.ones(leny1)]
    y_score = np.r_[y0, y1]
    fpr, tpr, _ = roc_curve(y_true, y_score)
    auc = roc_auc_score(y_true, y_score)
    lab = label if label is not None else ""
    ax.plot(fpr, tpr, label=lab + ' (' + str(round(auc, 3)) + ')', linewidth=3)
    return auc

def ROCAUC(y0, y1, score_ascending=False):
    if not _HAVE_SKLEARN:
        raise ImportError("scikit-learn is required for ROC AUC.")
    if not score_ascending:
        y0 = -y0
        y1 = -y1
    y_true = np.r_[np.zeros(len(y0)), np.ones(len(y1))]
    y_score = np.r_[y0, y1]
    return roc_auc_score(y_true, y_score)

def hist(y, color='blue', label=None):
    import matplotlib.pyplot as plt
    plt.hist(y, bins=50, density=True, alpha=0.5, color=color, label=label)


# ----------------------------- core class -----------------------------

class PSlim:
    """
    Thin wrapper to:
      1) render Slim.txt -> .slim (placeholder filling)
      2) run SLiM with a seed & simID
      3) run recap.py to recapitate/mutate/subsample
      4) (optionally) run infertrees.py to compute extra per-tree stats
    """

    def __init__(self, Args):
        """
        Args: dict of parameters used by Slim.txt placeholders AND by recap.py
          required keys used by Slim.txt:
            'slimTxt','L','h','s','r','Ne','extmut','Until','start','rep','Freq','reset_lost'
          required keys used by recap.py:
            'N','Ne','r','mu'
          optional HPC keys:
            'srun' (None or a Slurm job runner)
        """
        self.Args = Args
        _ensure_dirs('slimfiles', 'trees', 'checkpoints', 'outs')
        self.gen_slimFile()

    def gen_slimFile(self):
        """
        Render Slim.txt into a runnable .slim script.
        Only replace [Identifier]-style placeholders to avoid eating Eidos bracket syntax.
        """
        Args = self.Args
        with open(Args['slimTxt'], 'r') as f:
            slimFile = f.read()

        # Split into literal + placeholder chunks; only identifiers inside brackets count
        L = re.split(r'(\[[A-Za-z_][A-Za-z0-9_]*\])', slimFile)

        # Build a short, stable output filename: <base>_<name or tag>.slim
        base = os.path.splitext(os.path.basename(Args['slimTxt']))[0]
        tag = Args.get('name', f"{Args.get('s','s')}_{Args.get('Freq','F')}")
        safe_tag = re.sub(r'[^A-Za-z0-9._-]+', '_', str(tag))
        out_path = os.path.join('slimfiles', f"{base}_{safe_tag}.slim")

        # Fill placeholders
        for i in range(1, len(L), 2):
            varname = L[i][1:-1]
            if varname not in Args:
                raise KeyError(f"Placeholder [{varname}] has no value in Args.")
            L[i] = str(Args[varname])

        with open(out_path, 'w') as f:
            f.write(''.join(L))

        Args['slimFile'] = out_path
        self.Args = Args

    def _run(self, cmd, check=True):
        print(f"running cmd with _run: {cmd}")
        p = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        print(f"cmd with _run done: {p.stdout}")
        if check and p.returncode != 0:
            raise RuntimeError(f"Command failed:\nCMD: {cmd}\nSTDOUT:\n{p.stdout}\nSTDERR:\n{p.stderr}")
        return p.stdout

    def _run(self, args, check=True, log_path=None, cwd=None):
        """
        Run a subprocess with live streaming, no shell.
        args: list[str] (e.g., ["slim","-s","3000000",...])
        log_path: optional file to tee output into
        cwd: optional working directory
        """
        if not isinstance(args, (list, tuple)):
            raise TypeError(f"_run expects list/tuple, got: {args!r}")

        # ensure log dir exists if logging
        lf = None
        if log_path:
            os.makedirs(os.path.dirname(log_path), exist_ok=True)
            lf = open(log_path, "w")

        print(f"running cmd with _run: {args}")
        t0 = time.time()
        p = subprocess.Popen(
            args,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            cwd=cwd,
        )
        try:
            for line in p.stdout:
                sys.stdout.write(line)   # show live
                if lf: lf.write(line)    # also save
            p.wait()
        finally:
            if lf: lf.close()

        dt = time.time() - t0
        print(f"cmd with _run done in {dt:.1f}s (returncode={p.returncode})")
        if check and p.returncode != 0:
            raise RuntimeError(f"Command failed ({args})")
        return p.returncode


    def sim(self, simID, itrees=True):
        """
        Run one replicate:
        SLiM -> trees/<simID>.trees
        recap.py -> trees/r<simID>.trees
        (optional) infertrees.py -> additional stats
        """
        Args = self.Args.copy()
        slimFile = Args['slimFile']
        simID_str = str(simID)

        # Build commands as argument lists (no shell splitting)
        cmd1 = ["slim", "-s", simID_str, "-d", f"simID={simID_str}", slimFile]
        # cmd2 = [sys.executable, os.path.abspath("recap.py"),
        #         simID_str, str(Args['N']), str(Args['Ne']),
        #         str(Args['r']), str(Args['mu'])]

        cmd2 = ["python", "recap.py", simID_str, str(self.Args['N']), str(self.Args['Ne']), str(self.Args['r']), str(self.Args['mu'])]

        print("commands set up")

        srun = Args.get('srun', None)

        if srun is None:
            # Run locally
            print("running cmd1")
            self._run(cmd1)   # SLiM
            print("running cmd2")
            self._run(cmd2)   # recap
            if itrees:
                rtrees = f"trees/r{simID_str}.trees"
                cmd3 = [sys.executable, os.path.abspath("infertrees.py"), rtrees]
                try:
                    self._run(cmd3)
                except Exception as e:
                    print(f"[warn] infertrees.py not run or failed for {rtrees}: {e}")
        else:
            # HPC path: submit combined jobs as strings for Slurm
            jobs = [
                " ".join(cmd1),
                " ".join(cmd2),
            ]
            if itrees:
                rtrees = f"trees/r{simID_str}.trees"
                jobs.append(f"{sys.executable} infertrees.py {rtrees}")
            srun.run(jobs)

        return "ok"


# ----------------------------- experiment orchestration -----------------------------

class Experiment:
    """
    Orchestrates:
      - generating simIDs per set
      - running simulations per set
      - computing SFS (calc_sfs)
      - (optional) training eta (train_eta)
      - estimating stats over many trees with BIM script (est)
      - merging outputs (merge_outs)
    """

    def __init__(self, cID, nrep, Args):
        """
        cID: class id (int)
        nrep: replicates per condition
        Args: dict mapping setid -> arg dict used by PSlim + metadata:
              - must include 'N' for SFS shape
              - optionally: 'etapath' for train_eta output
        """
        _ensure_dirs('trees', 'outs')
        self.cID = cID
        self.nrep = nrep
        self.Args = Args
        self.setids = setids = list(Args.keys())
        self.neutrals = [sid for sid in setids
                         if (Args[sid].get('s', None) == 0) and (Args[sid].get('h', None) == 0.5)]

        self.outs = {}
        self.simIDs = {}
        self.setReg = {}
        self.df = {}

        for i, setid in enumerate(setids):
            self.simIDs[setid], self.setReg[setid] = get_simIDs(cID, i, nrep)
            self.outs[setid] = []
            self.df[setid] = os.path.join('outs', f"{cID}{setid}.csv")

        self.AFS = {sid: np.zeros(Args[sid]['N'] - 1, dtype='int') for sid in setids}

    def sim(self):
        setids = self.setids
        simIDs = self.simIDs
        simJobs = {}
        for setid in setids:
            print(f"running setid: {setid}")
            pslim = PSlim(self.Args[setid])
            print(f"pslim set up for {setid}")
            simJobs[setid] = []
            for simID in simIDs[setid]:
                simJobs[setid].append(pslim.sim(simID))
                print(f"simID: {simID} done")
        self.simJobs = simJobs
        print('If you are using HPC (srun is not None) check the jobs!')

    def calc_sfs(self):
        setReg = self.setReg
        setids = self.setids
        Args = self.Args
        AFS = {sid: np.zeros(Args[sid]['N'] - 1, dtype='int') for sid in setids}

        for setid in setids:
            regex = setReg[setid]
            tpaths = treepaths(regex)
            for tpath in tpaths:
                ts = tskit.load(tpath)
                # [1:-1] excludes singletons of 0 and fixed-of-1
                AFS[setid] += ts.allele_frequency_spectrum(
                    span_normalise=False, polarised=True
                )[1:-1].astype('int')

        self.AFS = AFS

    def train_eta(self):
        """
        Fit eta (demography) from neutral SFS and save to JSON at Args[setid]['etapath'].
        Requires utils.InferEta; skip if you don't use it.
        """
        if 'InferEta' not in globals():
            raise ImportError("InferEta not available. Uncomment its import if you need train_eta().")

        Args = self.Args
        AFS = self.AFS
        neutrals = self.neutrals

        # breakpoints (generations)
        t = np.logspace(np.log10(1), np.log10(20000), 100)
        t = np.concatenate((np.array([0]), t))

        a1 = 0.  # L1 penalty
        a2 = 1.  # L2 penalty

        self.ebl = {}
        for setid in neutrals:
            sfs = AFS[setid]
            inferEta = InferEta(Args[setid]['N'], t, a1=a1, a2=a2)
            etaout = inferEta.predict(sfs, maxiter=1000)
            eps = etaout.x
            ebl = inferEta.get_esfs(eps)
            self.ebl[setid] = ebl

            ti = t
            ai = np.float16(etaout.x)
            diff = np.diff(ai).nonzero()[0]
            ti = np.r_[ti[0], ti[diff + 1]]
            ai = np.r_[ai[0], ai[diff + 1]]

            outname = Args[setid]['etapath']
            with open(outname, 'w') as fp:
                json.dump({0: {'t': ti.tolist(), 'a': (1 / ai).tolist()}}, fp)
        print('Done!')

    def est(self, BIM, setid, now=2, srun=None, arg="", onlysimloc=False):
        import os, shutil, sys, numpy as np, subprocess

        Arg = self.Args[setid]
        eta = Arg['etapath']

        # Build list from disk: include r*.trees and ir*.trees if present
        simIDs = self.simIDs[setid]
        want = [f"../Simulation/Slim(primary)/trees/r{i}.trees" for i in simIDs] + [f"../Simulation/Slim(primary)/trees/ir{i}.trees" for i in simIDs]
        tree_paths = [p for p in want if os.path.exists(p)]
        if not tree_paths:
            raise RuntimeError(f"No trees found for setid={setid}. Looked for r/ir files.")

        n = len(tree_paths)
        chunk_sizes = [n // now] * now
        for k in range(n - sum(chunk_sizes)):
            chunk_sizes[k] += 1

        # normalize flags
        if isinstance(arg, str):
            arg = arg.replace("--stats=", "--stat=")   # old -> new
            arg = arg.replace("--treew=", "--weights=")
            extra_args = arg.split()
        else:
            extra_args = [a.replace("--stats=", "--stat=").replace("--treew=", "--weights=") for a in arg]

        # how to invoke
        use_cli = (os.path.sep in BIM and os.access(BIM, os.X_OK)) or shutil.which(BIM)
        prefix = [BIM] if use_cli else [sys.executable, BIM]

        msgs = []
        self.outs[setid] = []
        it = iter(tree_paths)
        for sz in chunk_sizes:
            if sz == 0:
                continue
            if onlysimloc:
                out = os.path.join('outs_onlysimloc', f"{np.random.randint(1e8, 1_000_000_000-1):09d}.csv")
            else:
                out = os.path.join('outs', f"{np.random.randint(1e8, 1_000_000_000-1):09d}.csv")
            self.outs[setid].append(out)
            chunk = [next(it) for _ in range(sz)]
            argv = prefix + chunk + [
                str(Arg['N']),
                "--stat=all",
                f"--eta={eta}",
                f"--out={out}",
                *extra_args,
            ]
            if srun is None:
                p = subprocess.run(argv, capture_output=True, text=True)
                if p.returncode != 0:
                    raise RuntimeError(f"[BIM failed]\nARGV: {argv}\nSTDOUT:\n{p.stdout}\nSTDERR:\n{p.stderr}")
                msgs.append(p.stdout)
            else:
                msgs.append(srun.run(" ".join(argv)))
        print("If you are using HPC (srun is not None) check the jobs!")
        return msgs

    # def merge_outs(self, setid, name=None):
    #     outs = self.outs.get(setid, [])
    #     if not outs:
    #         raise ValueError(f"No outs recorded for setid={setid}. Did you run est()?")

    #     if name is None:
    #         name = self.df[setid]

    #     df = pd.concat([pd.read_csv(out, comment='#') for out in outs]).reset_index(drop=True)

    #     dfr = df[[p[0] == 'r' for p in df['path']]].copy()
    #     dfi = df[[p[0] == 'i' for p in df['path']]].copy()

    #     dfr['path'] = [int(re.findall('[0-9]+', p)[0]) for p in dfr['path']]
    #     dfi['path'] = [int(re.findall('[0-9]+', p)[0]) for p in dfi['path']]
    #     dfi = dfi.rename(columns={'Colless': 'iColless', 'btree': 'ibtree'})[['path', 'iColless', 'ibtree']]
    #     dfm = dfr.merge(dfi, on='path')
    #     dfm.sort_values('path').to_csv(name, index=False)

    #     # cleanup shards
    #     for out in outs:
    #         try:
    #             os.remove(out)
    #         except FileNotFoundError:
    #             pass

    def merge_outs_nonstrict(self, setid, name=None, onlysimloc=False):
        outs = self.outs.get(setid, [])
        if not outs:
            raise ValueError(f"No outs recorded for setid={setid}. Did you run est()?")

        if name is None:
            if onlysimloc:
                # Extract the filename from self.df[setid] and use outs_onlysimloc directory
                base_name = os.path.basename(self.df[setid])
                name = os.path.join('outs_onlysimloc', base_name)
                _ensure_dirs('outs_onlysimloc')
            else:
                name = self.df[setid]

        df = pd.concat([pd.read_csv(out, comment='#') for out in outs]).reset_index(drop=True)

        # rows from recap (r*.trees)
        dfr = df[df['path'].str.startswith('r')].copy()
        dfr['path'] = dfr['path'].str.extract(r'(\d+)').astype(int)

        # rows from inferred (i*.trees) — may be absent
        has_inferred = (df['path'].str.startswith('i')).any()
        if has_inferred:
            dfi = df[df['path'].str.startswith('i')].copy()
            dfi['path'] = dfi['path'].str.extract(r'(\d+)').astype(int)
            dfi = dfi.rename(columns={'Colless':'iColless','btree':'ibtree'})[['path','iColless','ibtree']]
            dfm = dfr.merge(dfi, on='path', how='left')
        else:
            dfm = dfr  # no merge; just keep recap stats

        dfm.sort_values('path').to_csv(name, index=False)

        # cleanup shard files
        for out in outs:
            try: os.remove(out)
            except FileNotFoundError: pass

    # monkey-patch
    # Experiment.merge_outs = merge_outs_nonstrict