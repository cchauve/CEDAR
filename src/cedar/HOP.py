"""
CEDAR: manipulating phylogenetic rooted trees representations as vectors
HOP functions
"""

__author__ = "Cedric Chauve"
__credits__ = ["Cedric Chauve", "Louxin Zhang"]
__license__ = "GPL"
__version__ = "1.2.0"
__maintainer__ = "Cedric Chauve"
__email__ = "cedric.chauve@sfu.ca"
__status__ = "Release"


from cedar.utils import (
    __write_file
)
from cedar.LeavesOrder import (
    order2str
)
from cedar.TreeVec import (
    read_TreeVec_file,
    write_TreeVec_file,
)

def __hop_similarity_list(in_TreeVec_trees, mode="sequence"):
    """
    Given a list of TreeVec trees, compute the hop similarity between them
    - in_TreeVec_trees: list(TreeVec) list of TreeVec objects
    - if mode == "sequence", computes the similarity between successive trees
    - if mode == "pairwise", computes the similarity between all pairs of different trees
    - if mode == "first", computes the similarity betwen the first tree and all other trees
    Output:
    - list(int,int,int): (tree1,tree2,similarity)
    """
    similarity_list = []
    nb_trees = len(in_TreeVec_trees)
    for i in range(0,nb_trees-1):
        tree1 = in_TreeVec_trees[i]
        if mode == "sequence":
            range_j = [i+1]
        elif mode == "pairwise":
            range_j = [j for j in range(i+1,nb_trees)]
        elif mode == "first":
            range_j = [j for j in range(1,nb_trees) if i == 0]
        for j in range_j:
            tree2 = in_TreeVec_trees[j]
            similarity = tree1.hop_similarity(tree2, compute_seq=False)
            similarity_list.append([i+1,j+1,similarity])
    return similarity_list

def hop_similarity(in_TreeVec_file, out_dist_file, mode="sequence"):
    """
    Reads a file of TreeVec trees and computes the hop similarity between them
    If mode=="sequence" computes the similarity between successive trees
    If mode=="pairwise" computes the similarity between all pairs of trees
    If mode=="first", computes the similarity betwen the first tree and all other trees
    Input:
    - in_treevec_file (str): path to an existing treevec file
    - out_dist_file (str): path to the CSV distances file
      format: tree1,tree2,similarity
    - mode (str): "sequence" or "pairwise"
    Output:
    - None (writes file out_dist_file)
    """
    TreeVec_trees,leaf2idx,idx2leaf = read_TreeVec_file(in_TreeVec_file)
    # Computing the hop similarity
    similarity = __hop_similarity_list(TreeVec_trees, mode=mode)
    # Writing the output
    out_str = [
        f"#order {order2str(idx2leaf)}", f"#tree1,tree2,similarity"
    ] + [
        f"{i},{j},{s}" for [i,j,s] in similarity
    ]
    __write_file(out_str, out_dist_file)

def hop_neighbourhood_size(in_TreeVec_file, out_size_file):
    """
    Computes the size of the hop neighbourhood for all trees in a TreeVec file
    Output format: tree_id,size
    """
    TreeVec_trees,leaf2idx,idx2leaf = read_TreeVec_file(in_TreeVec_file)
    ngb_size = []
    for TreeVec_tree in TreeVec_trees:
        ngb_size.append(
            TreeVec_tree.hop_neighbourhood_size()
        )
    # Writing the output
    out_str = [
        f"#order {order2str(idx2leaf)}", f"#tree,neighbourhood_size"
    ] + [
        f"{i+1},{ngb_size[i]}" for i in range(0,len(ngb_size))
    ]
    __write_file(out_str, out_size_file)

def hop_neighbourhood(in_TreeVec_file, out_ngb_file_prefix="CEDAR_HOP_neighbourhood"):
    """
    Computes the hop neighbourhood for all trees in a TreeVec file
    """
    TreeVec_trees,leaf2idx,idx2leaf = read_TreeVec_file(in_TreeVec_file)
    i = 1
    for TreeVec_tree in TreeVec_trees:
        ngb = TreeVec_tree.hop_neighbourhood()
        out_TreeVec_file = f"{out_ngb_file_prefix}.{i}.vec"
        write_TreeVec_file(ngb, idx2leaf, out_TreeVec_file)
        i += 1

def hop_path(in_TreeVec_file, out_TreeVec_file):
    """
    Given a file of TreeVec trees, inserts between any pair of successive trees a list of
    trees forming a hop path; exports in TreeVec format
    """
    in_TreeVec_trees,leaf2idx,idx2leaf = read_TreeVec_file(in_TreeVec_file)
    nb_trees = len(in_TreeVec_trees)
    out_TreeVec_trees = []
    for i in range(0,nb_trees-1):
        tree1,tree2 = in_TreeVec_trees[i],in_TreeVec_trees[i+1]
        if i==0:
            out_TreeVec_trees.append(tree1)
        j = 0
        tree3 = tree1.hop_next(tree2)
        while tree3 is not None:
            out_TreeVec_trees.append(tree3)
            tree1 = tree3
            tree3 = tree1.hop_next(tree2)
            j += 1
    write_TreeVec_file(out_TreeVec_trees, idx2leaf, out_TreeVec_file)

def hop_random(in_TreeVec_file, out_TreeVec_file, in_seed):
    """
    Given a file of TreeVec trees, performs a random hop to each tree
    """
    in_TreeVec_trees,leaf2idx,idx2leaf = read_TreeVec_file(in_TreeVec_file)
    out_TreeVec_trees = []
    for TreeVec_tree in in_TreeVec_trees:
        new_TreeVec_tree = TreeVec_tree.random_hop(seed=in_seed, inplace=False)
        out_TreeVec_trees.append(new_TreeVec_tree)
    write_TreeVec_file(out_TreeVec_trees, idx2leaf, out_TreeVec_file)
