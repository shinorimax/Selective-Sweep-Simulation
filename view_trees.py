from ete3 import Tree
import re

# Read the file content
file_name = "all_scenarios.txt"

with open(file_name, "r") as file:
    content = file.read()

# Extract individual Newick trees (assuming they are separated by multiple new lines)
newick_trees = re.findall(r"\(.*?;\)", content, re.DOTALL)

# Iterate over the trees and display them
for newick_tree in newick_trees:
    t = Tree(newick_tree)
    t.show()  # This will open an interactive tree viewer