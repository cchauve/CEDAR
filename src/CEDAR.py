"""
CEDAR: manipulating phylogenetic rooted trees representations as vectors
"""

__author__ = "Cedric Chauve"
__credits__ = ["Cedric Chauve", "Louxin Zhang"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Cedric Chauve"
__email__ = "cedric.chauve@sfu.ca"
__status__ = "Release"

import os
import sys
import argparse
from numpy import random
import numpy as np

from TreeVec import TreeVec

# Auxiliary functions for manipulating leaves orders

# Separator between leaves names in a leaves order
SEP_ORDER = ","

def str2order(s, sep=SEP_ORDER):
    """
    Reads a leaves order from a string
    Input:
    - s (str)
    Output:
    - leaf2idx (dict str -> int): leaf name to index in a total order on leaves
      (1-base)  
    - idx2leaf (dict int -> str): reverse dictionary
    """
    leaves = s.rstrip().split(sep)
    n = len(leaves)
    leaf2idx = {
        leaves[i]: i+1 for i in range(0,n)
    }
    idx2leaf = {
        i: l for l,i in leaf2idx.items() 
    }
    return leaf2idx,idx2leaf

def order2str(idx2leaf, sep=SEP_ORDER):
    """
    Write a leaves order into a string
    Input:
    - leaves_order (dict int -> str)
    Output:
    - (str): string with leaf names in increasing order separated by SEP_ORDER
    """
    n = len(idx2leaf.keys())
    return sep.join([str(idx2leaf[i]) for i in range(1,n+1)])

def read_leaves_order_file(in_leaves_order_file, sep=SEP_ORDER):
    """
    Reads a leaves order from a file
    Input:
    - in_leaves_order_file (str): path to a file whose first line contains
      a leaves order encoded as the leaves names in increasing order separated
      by SEP_ORDER
    Output:
    - leaf2idx (dict str -> int): leaf name to index in a total order on leaves
      (1-base)  
    - idx2leaf (dict int -> str): reverse dictionary
    """
    with open(in_leaves_order_file) as in_file:        
        line = in_file.readlines()[0]
    leaf2idx,idx2leaf = str2order(line, sep=sep)
    return leaf2idx,idx2leaf

def random_leaves_order(in_TreeVec_file, nb_orders=1, out_file_prefix="leaves_order"):
    """
    Generates nb_orders random leaves orders and write them in files 
    {out_prefix_file}_{nb order}.txt
    """
    _,_,idx2leaf = __read_TreeVec_file(in_TreeVec_file)
    nb_leaves = len(idx2leaf.keys())
    idx = np.array(list(idx2leaf.values()))
    for i in range(1,nb_orders+1):
        random.shuffle(idx)
        _idx2leaf = {j: idx[j-1] for j in range(1,nb_leaves+1)}
        with open(f"{out_file_prefix}_{i}.txt", "w") as out_file:
            out_file.write(f"{order2str(_idx2leaf)}\n")
        
# Converting between Newick and TreeVec

def __read_file(input_file):
    """
    Read a file and returns a list of strings, one per line
    """
    lines = []
    with open(input_file) as in_file:
        for line in in_file.readlines():
            lines.append(line.rstrip())
    return lines

def __read_TreeVec_file(in_TreeVec_file):
    # Reading the TreeVec trees
    in_TreeVec_trees = __read_file(in_TreeVec_file)
    # Determining the leaves order
    leaf2idx,idx2leaf = str2order(in_TreeVec_trees[0].split()[1])
    # Creating TreeVec objects
    TreeVec_trees = []
    for TreeVec_tree in in_TreeVec_trees[1:]:
        TreeVec_trees.append(
            TreeVec(
                treevec_str=TreeVec_tree, idx2leaf=idx2leaf, format=1, compact=True
            )
        )
    return TreeVec_trees,leaf2idx,idx2leaf

def __write_file(in_lines, output_file):
    """
    Write a list of strings into a file
    """
    with open(output_file, "w") as out_file:
        for line in in_lines:
            out_file.write(f"{line}\n")

def __write_TreeVec_file(in_TreeVec_trees, idx2leaf, out_TreeVec_file):
    out_str = [f"#order {order2str(idx2leaf)}"]
    for TreeVec_tree in in_TreeVec_trees:
        out_str.append(
            TreeVec_tree.treevec2str(format=1, compact=True)
        )
    __write_file(out_str, out_TreeVec_file)
            
def convert_Newick2TreeVec(in_Newick_file, out_TreeVec_file, leaves_order_file=None):
    """
    Converts the trees in in_Newick_file into TreeVec strings written in out_TreeVec_file
    If leaves_order_file is not None, it is used to define the leaves order
    """
    # Reading the Newick trees
    in_Newick_trees = __read_file(in_Newick_file)
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
    __write_TreeVec_file(out_TreeVec_trees, idx2leaf, out_TreeVec_file)

def convert_TreeVec2Newick(in_TreeVec_file, out_Newick_file, Newick_format=0):
    """
    Converts the trees in in_TreeVec_file into Newick strings written in out_Newick_file
    """
    TreeVec_trees,leaf2idx,idx2leaf = __read_TreeVec_file(in_TreeVec_file)
    Newick_trees = []
    for TreeVec_tree in TreeVec_trees:
        Newick_trees.append(
            TreeVec_tree.treevec2newick(newick_format=Newick_format)
        )
    __write_file(Newick_trees, out_Newick_file)

# Hop similarity

def __hop_similarity(in_TreeVec_trees, mode="sequence"):
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
    TreeVec_trees,leaf2idx,idx2leaf = __read_TreeVec_file(in_TreeVec_file)
    # Computing the hop similarity
    similarity = __hop_similarity(TreeVec_trees, mode=mode)
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
    TreeVec_trees,leaf2idx,idx2leaf = __read_TreeVec_file(in_TreeVec_file)
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
    TreeVec_trees,leaf2idx,idx2leaf = __read_TreeVec_file(in_TreeVec_file)
    i = 1
    for TreeVec_tree in TreeVec_trees:
        ngb = TreeVec_tree.hop_neighbourhood()
        out_TreeVec_file = f"{out_ngb_file_prefix}.{i}.vec"
        __write_TreeVec_file(ngb, idx2leaf, out_TreeVec_file)
        i += 1

def hop_path(in_TreeVec_file, out_TreeVec_file):
    """
    Given a file of TreeVec trees, inserts between any pair of successive trees a list of
    trees forming a hop path; exports in TreeVec format
    """
    in_TreeVec_trees,leaf2idx,idx2leaf = __read_TreeVec_file(in_TreeVec_file)
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
    __write_TreeVec_file(out_TreeVec_trees, idx2leaf, out_TreeVec_file)

def _parse_arguments():
    description = "CEDAR: manipulating phylogenetic rooted trees representations as vectors"

    argparser = argparse.ArgumentParser(prog="CEDAR", description=description)
    
    subparsers = argparser.add_subparsers(title="commands", help="command help")

    # Converting Newick trees file to CEDAR trees file
    convert1 = subparsers.add_parser("fromNewick", help="Convert Newick file to CEDAR file")
    convert1.set_defaults(cmd="fromNewick")
    convert1.add_argument("--input_file", type=str, help="Input Newick file")
    convert1.add_argument("--output_file", type=str, help="Output CEDAR file")
    convert1.add_argument("--order_file", type=str, default=None, help="[OPTIONAL] leaves order file")

    # Converting CEDAR trees file to Newick trees file
    convert2 = subparsers.add_parser("toNewick", help="Convert CEDAR file to Newick file")
    convert2.set_defaults(cmd="toNewick")
    convert2.add_argument("--input_file", type=str, help="Input CEDAR file")
    convert2.add_argument("--output_file", type=str, help="Output Newick file")
    convert2.add_argument("--format", type=int, default=1, help="[OPTIONAL] Newick format (default 1)")

    # Creating random taxa orders
    orders = subparsers.add_parser("orders", help="Create random leaves oder files")
    orders.set_defaults(cmd="orders")
    orders.add_argument("--input_file", type=str, help="Input CEDAR file")
    orders.add_argument("--output_dir", type=str, help="Output directory")
    orders.add_argument("--nb_orders", type=int, default=1, help="[OPTIONAL] Number of random orders to generate (default=1)")
    orders.add_argument("--output_prefix", type=str, default="CEDAR_random_order", help="[OPTIONAL] Prefix of random order files")    

    # Computing the HOP similarity between trees
    hop_sim = subparsers.add_parser("HOP_sim", help="HOP similarity between trees")
    hop_sim.set_defaults(cmd="HOP_sim")
    hop_sim.add_argument("--input_file", type=str, help="Input CEDAR file")
    hop_sim.add_argument("--output_file", type=str, help="Output similarity file")
    hop_sim.add_argument("--mode", type=str, default="sequence", help="[OPTIONAL] Mode of similarity (sequence/pairwise/first, default=sequence)")

    # Computing the HOP neighbourhood sizes of a set of trees
    hop_ngb1 = subparsers.add_parser("HOP_ngb1", help="HOP neighbourhood size")
    hop_ngb1.set_defaults(cmd="HOP_ngb1")
    hop_ngb1.add_argument("--input_file", type=str, help="Input CEDAR file")
    hop_ngb1.add_argument("--output_file", type=str, help="Output neighbourhood sizes file")

    # Computing the HOP neighbourhoods of a set of trees
    hop_ngb2 = subparsers.add_parser("HOP_ngb2", help="HOP neighbourhood sizes")
    hop_ngb2.set_defaults(cmd="HOP_ngb2")
    hop_ngb2.add_argument("--input_file", type=str, help="Input CEDAR file")
    hop_ngb2.add_argument("--output_dir", type=str, help="Output directory")
    hop_ngb2.add_argument("--output_prefix", type=str, default="CEDAR_HOP_neighbourhood", help="[OPTIONAL] Prefix of output files")    

    # Computing a HOP path beween successive trees
    hop_path = subparsers.add_parser("HOP_path", help="Create HOP path between successive trees")
    hop_path.set_defaults(cmd="HOP_path")
    hop_path.add_argument("--input_file", type=str, help="Input CEDAR file")
    hop_path.add_argument("--output_file", type=str, help="Output CEDAR file")

    return argparser.parse_args()

def fromNewick(args):
    convert_Newick2TreeVec(
        args.input_file,
        args.output_file,
        leaves_order_file=args.order_file
    )

def toNewick(args):
    convert_TreeVec2Newick(
        args.input_file,
        args.output_file,
        Newick_format=args.format
    )

def orders(args):
    random_leaves_order(
        args.input_file,
        nb_orders=args.nb_orders,
        out_file_prefix=os.path.join(args.output_dir,args.output_prefix)
    )

def HOP_sim(args):
    if args.mode not in ["sequence", "pairwise", "first"]:
        raise Exception("The HOP similar mode must be eiher sequence, or pairwise or first")
    hop_similarity(
        args.input_file,
        args.output_file,
        mode=args.mode
    )

def HOP_ngb1(args):
    hop_neighbourhood_size(
        args.input_file,
        args.output_file
    )

def HOP_ngb2(args):
    hop_neighbourhood(
        args.input_file,
        out_ngb_file_prefix=os.path.join(args.output_dir,args.output_prefix)
    )

def HOP_path(args):
    hop_path(
        args.input_file,
        args.output_file
    )
    
def main(args):
    
    if args.cmd == "fromNewick":
        fromNewick(args)
        
    elif args.cmd == "toNewick":
        toNewick(args)

    elif args.cmd == "orders":
        orders(args)

    elif args.cmd == "HOP_sim":
        HOP_sim(args)

    elif args.cmd == "HOP_ngb1":
        HOP_ngb1(args)

    elif args.cmd == "HOP_ngb2":
        HOP_ngb2(args)

    elif args.cmd == "HOP_path":
        HOP_path(args)
        
    else:
        raise Exception("ERROR: Unknown command")

if __name__ == "__main__":
    
    args = _parse_arguments()

    main(args)
