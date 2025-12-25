What it simulates
	•	A single panmictic population (p1) of size Ne evolving forward in time with uniform recombination across a chromosome of length L.
	•	A de novo beneficial mutation is forcibly introduced at the chromosome midpoint (pmut = L/2) at tick 1. Its selection coefficient is s and dominance h.
	•	Optionally, a few additional linked selected mutations can be added symmetrically around the focal site (extmut flag), to mimic a cluster of selected/linked variants.
	•	No neutral mutations are generated during the forward run (initializeMutationRate(0)), so the output is a clean tree sequence with selection + recombination only; you typically add neutral mutations later (e.g., during recapitation/decoration in msprime).

Key pieces

initialize() block
	•	initializeTreeSeq();
Turns on tree-sequence recording so you get a compact .trees file capturing the full local genealogy along the genome.
	•	initializeMutationRate(0); and initializeMutationType("m1", …)
Defines a neutral type m1, but sets the forward neutral mut rate to 0 (neutral mutations added later off-line).
	•	initializeMutationType("m3", [h], "f", [s]);
Defines the selected mutation type with dominance h and selection s.
	•	initializeGenomicElementType("g1", m1, 1.0); initializeGenomicElement(g1, 0, L);
Covers the whole chromosome with a neutral element type (again, mutation rate 0 right now).
	•	initializeRecombinationRate([r]);
Uniform recombination at rate r across length L.

1 early()
	•	sim.addSubpop("p1", Ne);
Creates the population at tick 1 with size Ne.

1 late()
	•	Seeds the sweep: randomly picks one haplosome (haploid chromosome) and adds a new drawn mutation of type m3 at position pmut (the midpoint).
	•	Optional linked selected mutations (extmut): adds a few more m3 mutations at positions pmut ± 5*i (alternating sides), giving a small cluster around the target.
	•	Checkpoint: immediately writes checkpoints/slim_<simID>.trees.
	•	This snapshot is used to restart the attempt if the selected mutation fixes too fast or is lost (see “reset/restart” below).

1:[Until] late()

Runs every tick from 1 up to Until and controls the run/stop logic:
	•	Demography (optional growth/decline)
If community.tick > [start], resets p1 size to round(rep^(tick - start) * Ne).
→ This is exponential change after start, controlled by rep (e.g., growth if rep>1, contraction if <1).
	•	Sweep progress check
Computes the total frequency of all m3 mutations:
freqs = sum(sim.mutationFrequencies(NULL, m3muts));
If freqs >= [Freq], it writes trees/<simID>.trees and terminates.
→ This gives you a tree sequence stopped exactly when the beneficial allele reaches a target frequency (e.g., 0.7, 0.9), which is what we want for benchmarking selection-detection methods.
	•	Reset/restart logic ([reset_lost])
If there are no m3 mutations present:
	•	If there is exactly one m3 substitution, the allele fixed (too far); else, it was lost.
	•	In either case, it reloads the checkpoint (checkpoints/slim_<simID>.trees) and reseeds the RNG, then keeps trying until it hits the frequency target or Until.
→ This ensures you only keep replicates where the sweep reaches your desired frequency window.
	•	Progress logging every 1000 ticks.

[Until] late()
	•	Hard stop at the maximum time Until: prints the current frequency, writes trees/<simID>.trees, and exits.
→ Guarantees an output even if the target frequency wasn’t reached.

What files you get (and why they’re useful)
	•	checkpoints/slim_<simID>.trees: tree sequence right after introducing the selected mutation. Used internally for fast restart upon loss/fixation.
	•	trees/<simID>.trees: the final tree sequence for the replicate. It contains:
	•	Local genealogies along the genome under your recombination rate r.
	•	The selected mutation(s) on branches (type m3) at the specified position(s).
	•	No forward neutral mutations, keeping the file compact and clean for post hoc recapitation (extend ancestral history back in time) and neutral mutation overlay in msprime/pyslim.

What this provides for the selection-detection project
	•	Ground-truth sweep replicates at a controlled allele frequency stage (pre-fixation, near-fixation, etc.), under configurable demography and recombination.
	•	Compact tree-sequence outputs ideal for:
	1.	Recapitation (to add deep ancestral history consistent with coalescent)
	2.	Neutral mutation overlay (to create sequence data / SFS)
	3.	ARG or local-tree inference benchmarking (if you infer trees from the sequences)
	4.	Window scans comparing neutral vs. selective regions using your tree-distance or β-based statistics.

Parameter placeholders (the bracketed ones)
	•	[L] chromosome length; [r] recombination rate; [Ne] census size at tick 1.
	•	[s], [h] selection coefficient and dominance of the selected allele.
	•	[extmut] number of extra linked selected mutations to add around the focal site (0 = none).
	•	[Until] maximum number of ticks to simulate.
	•	[start], [rep] control the onset and rate of exponential size change.
	•	[Freq] the target total frequency of selected mutations that triggers finishing.
	•	[reset_lost] whether to auto-restart from the checkpoint if the allele fixes or is lost.
	•	simID is just a run label for reproducible file naming; reseeding uses it during restarts.