"""
CEDAR: manipulating phylogenetic rooted trees representations as vectors
Convergence of an MCMC using Gelman Rubin statistics
"""

__author__ = "Cedric Chauve"
__credits__ = ["Cedric Chauve", "Louxin Zhang"]
__license__ = "GPL"
__version__ = "1.2.0"
__maintainer__ = "Cedric Chauve"
__email__ = "cedric.chauve@sfu.ca"
__status__ = "Release"

import numpy as np

from cedar.Newick import (
    read_Newick_file
)
from cedar.LeavesOrder import (
    order2str
)
from cedar.TreeVec import (
    TreeVec,
    get_nb_taxa
)

def __hop_min_distance(in_TreeVec_tree_1, in_TreeVec_tree_2, nb_taxa):
    """
    Computing a distance between two trees by taking the min HOP distance
    for a list of encodings of each tree with different leaves orders.
    Input:
    - in_TreeVec_tree_1 (list(TreeVec)): list of TreeVec objects
    - in_TreeVec_tree_2 (list(TreeVec)): list of TreeVec objects
    - nb_taxa (int): number of taxa in trees
    Assumptions:
    - lists in_TreeVec_tree_1 and in_TreeVec_tree_2 are of same length
    - trees in_TreeVec_tree_1[i] and in_TreeVec_tree_2[i] are encoded using the
      same leaves order
    Output:
    - int: HOP distance
    """
    nb_orders = len(in_TreeVec_tree_1)
    distances = [
        nb_taxa - in_TreeVec_tree_1[i].hop_similarity(
            in_TreeVec_tree_2[i], compute_seq=False
        )
        for i in range(nb_orders)
    ]
    return min(distances)

def __hop_distance_within(in_TreeVec_trees_1):
    """
    Given a list of trees, compute the hop distance between all pairs of trees
    Input:
    - in_TreeVec_trees_1: dict(int: list(TreeVec)) dict of list of TreeVec objects
    Assumptions:
    - in_TreeVec_trees_1.keys() is a range from 0 to number of trees
    - all in_TreeVec_trees_1[i] have the same length
    - in_TreeVec_trees_1[i] contains encodingof the same tree with different random orders
    - in_TreeVec_trees_1[i1][j] and in_TreeVec_trees_1[i2][j] are encoded with the same leaves order
    Output:
    - np.array where entry [i][j] is the distance between trees i from list 1 and j from list 1 (0-base index)
    """
    nb_trees = len(in_TreeVec_trees_1.keys())
    nb_taxa = get_nb_taxa(in_TreeVec_trees_1[0][0])
    distances = np.zeros([nb_trees, nb_trees])
    for i in range(0,nb_trees-1):
        tree_1 = in_TreeVec_trees_1[i]
        range_j = [j for j in range(i+1,nb_trees)]
        for j in range_j:
            tree_2 = in_TreeVec_trees_1[j]
            distances[i][j] = __hop_min_distance(tree_1, tree_2, nb_taxa)
            distances[j][i] = distances[i][j]
    return distances

def __hop_distance_between(in_TreeVec_trees_1, in_TreeVec_trees_2):
    """
    Given two list of TreeVec trees, compute the hop distance between trees of list 1 an trees of list 2
    Input:
    - in_TreeVec_trees_1: dict(int: list(TreeVec)) dict of list of TreeVec objects
    - in_TreeVec_trees_2: dict(int: list(TreeVec)) dict of list of TreeVec objects
    Assumptions:
    - in_TreeVec_trees_1.keys() is a range from 0 to number of trees
    - in_TreeVec_trees_2.keys() == in_TreeVec_trees_1.keys()
    - all in_TreeVec_trees_?[i] have the same length
    - in_TreeVec_trees_?[i] contains encodingof the same tree with different random orders
    - in_TreeVec_trees_?[i1][j] and in_TreeVec_trees_?[i2][j] are encoded with the same leaves order
    Output:
    - np.array where entry [i][j] is the distance between trees i from list 1 and j from list 2 (0-base index)
    """
    nb_trees = len(in_TreeVec_trees_1.keys())
    nb_taxa = get_nb_taxa(in_TreeVec_trees_1[0][0])
    distances = np.zeros([nb_trees, nb_trees])
    for i in range(0,nb_trees-1):
        tree_1 = in_TreeVec_trees_1[i]
        for j in range(0,nb_trees):
            tree_2 = in_TreeVec_trees_2[j]
            distances[i][j] = __hop_min_distance(tree_1, tree_2, nb_taxa)
    return distances

def _GelmanRubin(in_TreeVec_trees_1, in_TreeVec_trees_2):
    """
    Given two lists of trees (chains), compute the Gelman Rubin diagnostic value for each tree of each chain
    Input:
    - in_TreeVec_trees_1: dict(int: list(TreeVec)) dict of list of TreeVec objects
    - in_TreeVec_trees_2: dict(int: list(TreeVec)) dict of list of TreeVec objects
    Assumptions:
    - in_TreeVec_trees_1.keys() is a range from 0 to number of trees
    - in_TreeVec_trees_2.keys() == in_TreeVec_trees_1.keys()
    - all in_TreeVec_trees_?[i] have the same length
    - in_TreeVec_trees_?[i] contains encoding of the same tree with different random leaves orders
    - in_TreeVec_trees_?[i1][j] and in_TreeVec_trees_?[i2][j] are encoded with the same leaves order
    Output:
    - dict(c: dict(i: float) c in [1,2]): entry[c][i] = GelmanRubin value for tree i in chain c
    """
    nb_trees = len(in_TreeVec_trees_1.keys())
    in_trees = {1: in_TreeVec_trees_1, 2: in_TreeVec_trees_2}
    # GelmanRubin value for tree of rank i in each chain
    GR =  {c: np.zeros(nb_trees) for c in [1,2]}
    for c1 in [1,2]:
        # -- Computing the needed distances
        # squared distances within chain c1
        squared_distances_within = np.square(__hop_distance_within(in_trees[c1]))
        # squared distances between chain c1 (index 1) and chain c2 (index 2)
        if c1 == 1:
            squared_distances_between = np.square(__hop_distance_between(in_trees[c1], in_trees[(c1%2)+1]))
        else:
            squared_distances_between = np.transpose(squared_distances_between)
        # -- Processing trees of chain c1
        # sum from 0 to i of row s of squared_distances_within and squared_distances_between
        sum_s_within,sum_s_between = {0: np.float64(0.0)},{0: np.float64(0.0)}
        for i in range(1,nb_trees):
            # Processing tree i
            for s in range(i+1):
                if s == i:
                    # Computing the sums as index s is seen for the first time
                    sum_s_within[s] = np.sum(squared_distances_within[s][:i+1])
                    sum_s_between[s] = np.sum(squared_distances_between[s][:i+1])
                else:
                    # Updating sum_s_* by adding the squared distance for (s,i)
                    sum_s_within[s] += squared_distances_within[s][i]
                    sum_s_between[s] += squared_distances_between[s][i]
                PSRF = np.sqrt(sum_s_between[s]/sum_s_within[s])
                GR[c1][i] += PSRF
            GR[c1][i] /= np.float64(i+1)
    return GR

def GelmanRubin(
        in_Newick_trees_file_1,
        in_Newick_trees_file_2,
        nb_trees,
        nb_orders,
        random_seed,
        out_gr_file,
        out_orders_file
):
    """
    Computes and writes in a file the Gelman Rubin values for the tail of two Newick trees files
    Input:
    - in_Newick_trees_file_1 (str): path to file for first lit of Newick trees (chain 1)
    - in_Newick_trees_file_1 (str): path to file for second list of Newick trees (chain 2)
    - nb_trees (int): number of trees (trailing trees) to conside in each file
    - nb_orders (int): number of random leaves orders
    - random_seed (int): random number generator seed
    - out_gr_file (str): path to output GelmanRubin statistics file
    - out_orders_file (str): path to output leaves orders file
    Format of out_gr_file:
    - line 1: index<TAB>chain1<TAB>chain2
    - line 2+: <index i><TAB><GR value for tree i in chain 1><TAB><GR value for tree i in chain 2>
    Format of out_orders_file:
    - line 1: seed<TAB>index<TAB>order
    - line 2+: <random generator seed><TAB><order index><TAB><leaves order>
    """
    # Extracting one list of Newick strings from each file; adding a ";" at the end
    in_Newick_trees = {
        1: read_Newick_file(in_Newick_trees_file_1),
        2: read_Newick_file(in_Newick_trees_file_2)
    }
    assert len(in_Newick_trees[1])>=nb_trees, "Chain 1 too short"
    assert len(in_Newick_trees[2])>=nb_trees, "Chain 2 too short"

    # Generating from each list of Newick strings a structure
    # dict(int: list(TreeVec)) indexed by number of trees in same order than in Newick files
    # and where each list is a list of TreeVec objects obtained with a different order on leaves
    in_TreeVec_trees = {c: {i: [] for i in range(nb_trees)} for c in [1,2]}
    leaves_orders = []
    rng = np.random.default_rng(random_seed)
    _tree = TreeVec(newick_str=in_Newick_trees[1][0])
    for _ in range(nb_orders):
        # Generating a random leaves order
        _reordered_tree = _tree.reorder_leaves(rng)
        leaf2idx,idx2leaf = _reordered_tree.extract_leaves_order()
        leaves_orders.append(order2str(idx2leaf))
        # Populating the dictionary of lis of trees
        for c in [1,2]:
            i = 0
            for in_Newick_tree in in_Newick_trees[c][-nb_trees:]:
                _tree = TreeVec(newick_str=in_Newick_tree, leaf2idx=leaf2idx)
                in_TreeVec_trees[c][i].append(_tree)
                i += 1
    # Computing Gelman Rubin statistics
    GR = _GelmanRubin(in_TreeVec_trees[1], in_TreeVec_trees[2])
    # Writing output files
    with open(out_orders_file, "w") as _out_file:
        _out_file.write("seed\tindex\torder\n")
        i = 0
        for order_str in leaves_orders:
            _out_file.write(f"{random_seed}\t{i}\t{order_str}\n")
            i += 1
    with open(out_gr_file, "w") as _out_file:
        _out_file.write("index\tchain1\tchain2\n")
        for i in range(nb_trees):
            _out_file.write(f"{i}\t{GR[1][i]}\t{GR[2][i]}\n")
