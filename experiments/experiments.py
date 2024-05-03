"""
Experiments for the paper "A Vecor Representation for Phylogenetic Trees
"""

import os
import sys
import gzip
import shutil

from os.path import dirname
sys.path.append("../src")

from ete3 import Tree
from SPR import random_spr

from TreeVec import (
    random_leaves_order,
    convert_Newick2TreeVec,
    convert_TreeVec2Newick,
    hop_similarity,
    hop_neighbourhood_size,
    hop_path
)

def gzip_file(in_file, out_file):
    with open(in_file, "rb") as f_in:
        with gzip.open(out_file, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

def experiment1():
    """
    Using 1000 random trees on 100 taxa and comparing the disk space for storing them in Newick format
    and in TreeVec format, both uncompressed and gzipped.
    """
    
    in_Newick_file = "random.trees.1000.nwk"
    out_TreeVec_file = "random.trees.1000.vec"

    print(f"EXP1\tConverting {in_Newick_file} into {out_TreeVec_file}")
    convert_Newick2TreeVec(in_Newick_file, out_TreeVec_file, leaves_order_file=None)
    size_Newick_file = os.stat(in_Newick_file).st_size
    size_TreeVec_file = os.stat(out_TreeVec_file).st_size
    print(f"EXP1\tUncompressed Newick file {in_Newick_file}:    \t{size_Newick_file} Bytes")
    print(f"EXP1\tUncompressed TreeVec file {out_TreeVec_file}: \t{size_TreeVec_file} Bytes")
    print(f"EXP1\tSpace usage ratio:\t{round(size_TreeVec_file/size_Newick_file,3)}")

    in_Newick_file_gz = f"exp1_{in_Newick_file}.gz"
    out_TreeVec_file_gz = f"exp1_{out_TreeVec_file}.gz"

    print(f"EXP1\tCompressing {in_Newick_file} into {in_Newick_file_gz}")
    gzip_file(in_Newick_file, in_Newick_file_gz)
    print(f"EXP1\tCompressing {out_TreeVec_file} into {out_TreeVec_file_gz}")
    gzip_file(out_TreeVec_file, out_TreeVec_file_gz)

    size_Newick_file_gz = os.stat(in_Newick_file_gz).st_size
    size_TreeVec_file_gz = os.stat(out_TreeVec_file_gz).st_size
    print(f"EXP1\tCompressed Newick file {in_Newick_file_gz}:     \t{size_Newick_file_gz} Bytes")
    print(f"EXP1\tCompressed TreeVec file {out_TreeVec_file_gz}:  \t{size_TreeVec_file_gz} Bytes")
    print(f"EXP1\tSpace usage ratio:\t{round(size_TreeVec_file_gz/size_Newick_file_gz,3)}")

def experiment2():
    """
    Using the first 50 random trees on 100 taxa of experiment 1, adding hop trees beween them, and
    comparing the disk space for storing them in Newick format and in vector format, both uncompressed
    and gzipped.
    """

    in_TreeVec_file = "random.trees.50.vec"
    out_TreeVec_file = "exp2_random.trees.50.path.vec"
    out_Newick_file = "exp2_random.trees.50.path.nwk"
    out_sim_file = "exp2_random.trees.50.sim.txt"

    print(f"EXP2\tComputing hop similarity in file {out_sim_file}")
    hop_similarity(in_TreeVec_file, out_sim_file, mode="sequence")
    
    print(f"EXP2\tComputing hop paths in file {out_TreeVec_file}")
    hop_path(in_TreeVec_file, out_TreeVec_file)
    print(f"EXP2\tConverting {out_TreeVec_file} into {out_Newick_file}")
    convert_TreeVec2Newick(out_TreeVec_file, out_Newick_file, Newick_format=0)
    
    size_TreeVec_file = os.stat(out_TreeVec_file).st_size
    size_Newick_file = os.stat(out_Newick_file).st_size
    print(f"EXP2\tUncompressed Newick file {out_Newick_file}:   \t{size_Newick_file} Bytes")
    print(f"EXP2\tUncompressed TreeVec file {out_TreeVec_file}: \t{size_TreeVec_file} Bytes")
    print(f"EXP2\tSpace usage ratio:\t{round(size_TreeVec_file/size_Newick_file,3)}")

    out_Newick_file_gz = f"{out_Newick_file}.gz"
    out_TreeVec_file_gz = f"{out_TreeVec_file}.gz"

    print(f"EXP2\tCompressing {out_Newick_file} into {out_Newick_file_gz}")
    gzip_file(out_Newick_file, out_Newick_file_gz)
    print(f"EXP2\tCompressing {out_TreeVec_file} into {out_TreeVec_file_gz}")
    gzip_file(out_TreeVec_file, out_TreeVec_file_gz)

    size_Newick_file_gz = os.stat(out_Newick_file_gz).st_size
    size_TreeVec_file_gz = os.stat(out_TreeVec_file_gz).st_size
    print(f"EXP2\tCompressed Newick file {out_Newick_file_gz}:     \t{size_Newick_file_gz} Bytes")
    print(f"EXP2\tCompressed TreeVec file {out_TreeVec_file_gz}:  \t{size_TreeVec_file_gz} Bytes")
    print(f"EXP2\tSpace usage ratio:\t{round(size_TreeVec_file_gz/size_Newick_file_gz,3)}")

def experiment3():
    """
    Computing the hop neighbourhood sizes
    """

    in_TreeVec_file = "random.trees.1000.vec"
    out_size_file = "exp3_results.txt"

    print(f"EXP3\tComputing hop neighbourhood size in file {out_size_file}")
    hop_neighbourhood_size(in_TreeVec_file, out_size_file)


def __exp4_compute_SPR(in_Newick_file, out_prefix, spr_min, spr_max, spr_step):
    out_dict = {}
    with open(in_Newick_file) as in_file:
        i = 1
        for Newick_str in in_file.readlines():
            out_Newick_file = os.path.join("exp4", f"{out_prefix}.{i}.nwk")
            out_dict[i] = out_Newick_file
            out_Newick_trees = [Newick_str.rstrip()]
            tree = Tree(Newick_str.rstrip(), format=1)
            p = 1
            for node in tree.traverse("postorder"):
                if not node.is_leaf():
                    node.name = str(p)
                    p += 1
            for j in range(spr_min,spr_max+1,spr_step):
                for k in range(0,spr_step):
                    nb_tries = random_spr(tree)
                out_Newick_trees.append(tree.write(format=1))
            with open(out_Newick_file, "w") as out_file:
                for tree_str in out_Newick_trees:
                    out_file.write(f"{tree_str}\n")
            i += 1
    return out_dict

def __exp4_compute_RF(in_Newick_file):
    out_RF = {}
    Newick_trees = []
    with open(in_Newick_file) as in_file:
        for Newick_str in in_file.readlines():
            Newick_trees.append(Tree(Newick_str.rstrip(), format=1))
    tree1 = Newick_trees[0]
    j = 1
    for tree2 in Newick_trees[1:]:
        [rf,max_rf] = tree1.robinson_foulds(tree2)[0:2]
        out_RF[j] = [rf,max_rf]
        j += 1
    return out_RF

def __exp4_compute_LCS(out_Newick_files):
    LCS_files = {}
    for i,out_Newick_file in out_Newick_files.items():
        LCS_files[i] = {}
        for j in range(1,11):
            order_file = os.path.join("exp4", f"order_{j}.txt")
            out_TreeVec_file = out_Newick_file.replace(".nwk",f".{j}.vec")
            convert_Newick2TreeVec(out_Newick_file, out_TreeVec_file, leaves_order_file=order_file)
            out_dist_file = f"{out_TreeVec_file}.dist"
            LCS_files[i][j] = out_dist_file
            hop_similarity(out_TreeVec_file, out_dist_file, mode="first")
    # i: starting ree, j: order
    return LCS_files

def __exp4_aggregate_results(RF, LCS_files, out_dist_file, spr_step, nb_taxa):
    with open(out_dist_file, "w") as out_file:
        out_file.write("tree,SPR_dist,order,dist_type,dist_value\n")
        for i,LCS_files_tree in LCS_files.items():
            # i: starting tree
            for k,RF_dist in RF[i].items():
                out_file.write(f"{i},{k*spr_step},1,RF,{RF_dist[0]}\n")
            for j,LCS_file_tree_order in LCS_files_tree.items():
                # j: order
                with open(LCS_file_tree_order) as in_file:
                    lines = in_file.readlines()[2:]
                    for line in lines:
                        tree1,tree2,sim = line.rstrip().split(",")
                        out_file.write(f"{i},{spr_step*int(tree2)},{j},HOP,{nb_taxa-int(sim)}\n")
                

    
def experiment4():
    """
    Comparing RF and LCS distances on trees with known SPR distance from a starting tree.
    We consider 50 starting trees, each of 100 taxa, and for each we build a sequence of trees by
    doing SPR, from 5 SPR to 100 SPR by steps of 5.
    We compute the RF distance and the LCS distance, for 10 random leaves orders.
    """
    in_Newick_file = "random.trees.50.nwk"
    in_TreeVec_file = "random.trees.50.vec"

    print(f"EXP4\tComputing SPR trees")
    out_Newick_files = __exp4_compute_SPR(in_Newick_file, "random.trees.50.spr", 5, 100, 5)
    print(f"EXP4\tComputing RF distance")
    RF = {}
    for i,out_Newick_file in out_Newick_files.items():
        RF[i] = __exp4_compute_RF(out_Newick_file)
    print(f"EXP4\tComputing random orders")
    random_leaves_order(
        in_TreeVec_file, nb_orders=10, out_file_prefix=os.path.join("exp4", "order")
    )
    print(f"EXP4\tComputing LCS distance")
    LCS_files = __exp4_compute_LCS(out_Newick_files)
    print(f"EXP4\tAggregating results")
    out_dist_file = "exp4_results.txt"
    __exp4_aggregate_results(RF, LCS_files, out_dist_file, 5, 100)


if __name__ == "__main__":
    experiment = sys.argv[1]

    if experiment == "exp1":
        experiment1()
    elif experiment == "exp2":
        experiment2()
    elif experiment == "exp3":
        experiment3()
    elif experiment == "exp4":
        experiment4()
