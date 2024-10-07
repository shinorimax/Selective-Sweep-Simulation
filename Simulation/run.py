import tskit
import msprime
from tqdm import tqdm 
import pandas as pd
import numpy as np
import io
import random

def select_scenario(sc, Ne, pos):
    #msprime model 
    sweep_model = msprime.SweepGenicSelection(
        position= pos,  
        start_frequency=1.0 / (2 * Ne),
        end_frequency=1.0 - (1.0 / (2 * Ne)),
        s = sc,
        dt=1e-6 ##double check what this does
    )
    return sweep_model

def get_random_tree ():
    num_trees = sim.num_trees
    start = int (np.round (num_trees / 2))
    rand_index = sim.at_index(random.randrange(num_trees-start-2))
    return rand_index

# main block, simulations here
selection_scenario = [0.01, 0.5, 1.0]   #change these accordingly, based on what selection coefficient we want to explore
file_name = "all_scenarios.txt"  #name of output file

Ne = 1e3 #initial population size
L = 1e6 #chromosome length
num_simulations = 3 #can change, warning: more simulations --> analysis takes significantly more time

o = open (file_name, "a")
for i in range (0, len (selection_scenario)):
    o.write ("NEW SELECTION SCENARIO #" + str (i) + "COEFFICIENT: " + str (selection_scenario [i]) + "\n"*3)
    for j in tqdm (range (0, num_simulations), colour = "green"):
        mut_loc = int (L/4)
        sim = msprime.sim_ancestry (
            samples = 25, #samples are by default diploid.
            model = [select_scenario (selection_scenario[i], Ne, mut_loc), msprime.StandardCoalescent()], 
            population_size = Ne,
            recombination_rate = 1e-7,
            sequence_length = L,
        )
        m = sim.at(mut_loc).index
        extract_trees = [sim.at_index (m), sim.at_index(m + 1), get_random_tree(), sim.last()]
        for tr in extract_trees:
            tr_str = str (tr.as_newick() + "\n\n\n\n\n")
            o.write(tr_str)
    o.write ("\n"*4)


o.close()