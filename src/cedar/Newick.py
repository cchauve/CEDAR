"""
CEDAR: manipulating phylogenetic rooted trees representations as vectors
Conversion to/from Newick files
"""

__author__ = "Cedric Chauve"
__credits__ = ["Cedric Chauve", "Louxin Zhang"]
__license__ = "GPL"
__version__ = "1.2.0"
__maintainer__ = "Cedric Chauve"
__email__ = "cedric.chauve@sfu.ca"
__status__ = "Release"

from cedar.utils import (
    __read_file,
    __write_file
)
from cedar.LeavesOrder import (
    read_leaves_order_file
)
from cedar.TreeVec import (
    TreeVec,
    read_TreeVec_file,
    write_TreeVec_file
)

def read_Newick_file(in_Newick_file):
    _Newick_trees = __read_file(in_Newick_file)
    Newick_trees = []
    for tree in _Newick_trees:
        if tree[-1] != ";":
            Newick_trees.append(f"{tree};")
        else:
            Newick_trees.append(tree)
    return Newick_trees

def convert_Newick2TreeVec(in_Newick_file, out_TreeVec_file, leaves_order_file=None):
    """
    Converts the trees in in_Newick_file into TreeVec strings written in out_TreeVec_file
    If leaves_order_file is not None, it is used to define the leaves order
    """
    # Reading the Newick trees
    in_Newick_trees = read_Newick_file(in_Newick_file)
    # Determining the leaves order if provided
    if leaves_order_file is not None:
        leaf2idx,idx2leaf = read_leaves_order_file(leaves_order_file)
    else:
        leaf2idx,idx2leaf = None,None
    # Converting trees
    out_TreeVec_trees = []
    for Newick_tree in in_Newick_trees:
        TreeVec_tree = TreeVec(newick_str=Newick_tree, leaf2idx=leaf2idx, idx2leaf=idx2leaf)
        if leaf2idx is None or idx2leaf is None:
            leaf2idx,idx2leaf = TreeVec_tree.extract_leaves_order()
        out_TreeVec_trees.append(TreeVec_tree)
    # Writing TreeVec trees
    write_TreeVec_file(out_TreeVec_trees, idx2leaf, out_TreeVec_file)

def convert_TreeVec2Newick(in_TreeVec_file, out_Newick_file, Newick_format=0):
    """
    Converts the trees in in_TreeVec_file into Newick strings written in out_Newick_file
    """
    TreeVec_trees,leaf2idx,idx2leaf = read_TreeVec_file(in_TreeVec_file)
    Newick_trees = []
    for TreeVec_tree in TreeVec_trees:
        Newick_trees.append(
            TreeVec_tree.treevec2newick(newick_format=Newick_format)
        )
    __write_file(Newick_trees, out_Newick_file)

