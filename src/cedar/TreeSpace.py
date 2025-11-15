"""
CEDAR: exploration of the tree space in a maximum likelihood framework
"""

__author__ = "Cedric Chauve"
__credits__ = ["Cedric Chauve", "Louxin Zhang"]
__license__ = "GPL"
__version__ = "1.2.0"
__maintainer__ = "Cedric Chauve"
__email__ = "cedric.chauve@sfu.ca"
__status__ = "Release"

import os
import numpy as np
from Bio import SeqIO
from ete3 import Tree

from cedar._raxml import (
    raxml_loss
)
from cedar.TreeVec import (
    TreeVec
)

def _read_fasta(data_path):
    """
    Read a fasta file

    Input:
    - data_path (str): Path to fasta file
    Output:
    - Generator
        An iterator for all sequences
        Each element has "id" (taxon name) and "seq" (sequence) attributes
    """
    return SeqIO.parse(data_path, "fasta")

def _create_random_tree(leaves):
    """
    Create a random tree as a TreeVec object

    Input:
    - leaves (list(str)): List of leaves labels
    Output:
    - TreeVec: TeeVec object
    """
    tree_aux = Tree()
    nb_leaves = len(leaves)
    tree_aux.populate(nb_leaves, names_library=leaves)
    tree = TreeVec(tree=tree_aux)
    return tree

def _clean_tree(tree_folder_path, tree_file, clean_tree):
    _tree_file = os.path.join(tree_folder_path, tree_file)
    if clean_tree and os.path.exists(_tree_file):
        os.remove(_tree_file)

def _neighborhood_ML_scores(
        current_tree,
        fasta_path, DNA_model,
        tree_folder_path, iteration_counter, clean_tree
):
    """
    Compute the likelihood of all trees in the neighbourhood of a tree using Raxml

    Input:
    - current_tree (TreeVec): current tree
    - fasta_path (str): pah to FASTA file
    - DNA_model (str): DNA evolution model to use in Raxml
    - tree_folder_path (str): pah to folder where Newick trees are written
    - iteration_counter (int): iteration counter in hill-climbing process
    - clean_tree (bool): if True the created tree file is deleted after Raxml runs
    Output:
    - ngb_scores (list((int,int),float)):
      list of the hop (i,j) and likelihood score of all trees in the neighbourhood
    """
    # Set of hops defining the neigbourhood of current tree
    ngb_hops = current_tree.hop_neighbourhood(export_list=True)
    # Output list
    ngb_scores = []
    # Loop on possible hops
    for (i,j) in ngb_hops:
        # New tree obtained by prforming hop (i,j)
        ngb_tree = current_tree.hop(i,j,inplace=False)
        # Liklihood score of new_tree
        _tree_file = f"tmp_{iteration_counter}_{i}_{j}.tree"
        ngb_tree_score = raxml_loss(
            fasta_path, ngb_tree, DNA_model,
            tree_folder_path, outfile=_tree_file
        )
        _clean_tree(tree_folder_path, _tree_file, clean_tree)
        ngb_scores.append(((i,j),ngb_tree_score))
    return ngb_scores

def best_ML_neighbour(
        current_tree, current_score,
        fasta_path, DNA_model,
        tree_folder_path, iteration_counter, clean_tree
):
    """
    Chose the next tree as the best neighbouring tree in terms of likelihood score

    Input:
    - current_tree (TreeVec): current tree
    - current_score (float): score of current tree
    - fasta_path (str): pah to FASTA file
    - DNA_model (str): DNA evolution model to use in Raxml
    - tree_folder_path (str): pah to folder where Newick trees are written
    - iteration_counter (int): iteration counter in hill-climbing process
    - clean_tree (bool): if True the created tree file is deleted after Raxml runs
    Output:
    - best_ngb_tree (TreeVec): best neighbour tree
    - best_ngb_score (float): score of best neighbour tree
    """
    # Computing the likelihood of all trees in the neighbourhood
    ngb_scores_and_hops = _neighborhood_ML_scores(
        current_tree, fasta_path, DNA_model,
        tree_folder_path, iteration_counter, clean_tree
    )
    # Isolating liklihood scores in the same order than in ngb_scores_and_hops
    ngb_scores = [ngb_result[1] for ngb_result in ngb_scores_and_hops]
    # Index of the best likelihood score
    best_score_idx = np.argmin(ngb_scores)
    # HOP defining the best tree
    (i,j) = ngb_scores_and_hops[best_score_idx][0]
    # Creating the best tree by perforing the corresponding hop
    best_ngb_tree = current_tree.hop(i,j,inplace=False)
    best_ngb_score = ngb_scores[best_score_idx]
    return best_ngb_tree,best_ngb_score

def _write_tree(tree, score, prefix, out_file, sep):
    """
    Write "<prefix><sep><score><sep><tree newick><sep><tree treevec>" in out_file
    """
    tree_newick = tree.treevec2newick()
    tree_treevec = tree.treevec2str(format_str=1,compact=False)
    out_file.write(f"{prefix}{sep}{score}{sep}{tree_newick}{sep}{tree_treevec}\n")

def _hill_climbing(
        fasta_path, DNA_model, tree_folder_path,
        first_tree, rng,
        tol, max_patience, max_nb_iterations,
        out_file_path,
        sep,
        clean_tree
):
    """
    Using a Hill-Climbing heuristic to find a maximum likelihood tree from a set of
    FASTA sequences.

    Algorithm:
    current tree = first_tree
    iterate
    - reorder randomly the leaves of the current tree
    - compute the likelihood of all trees in the hop-neighbourhood
    - select the highest likelihood tree
    - if its likelihood is within tol of the best tree so far:
      - decrease the patience counter
    - otherwise
      - the best tree becomes the current tree
      - the patience counter is set to max_patience
    until the max number of iterations is reached or the patience counter is 0

    Input:
    - fasta_path (str): pah to FASTA file
    - DNA_model (str): DNA evolution model to use in Raxml
    - tree_folder_path (str): pah to folder where Newick trees are written
    - first_tree (TreeVec): starting tree
    - rng (Generator): numpy random number generator
    - tol (float): tolerance used to define a stopping criterion
    - max_patience (int): maximum number of iterations of rordering leaves
    - max_nb_iterations (int): maximum number of iterations
    - out_file_path (str): path of the file where all explored trees are recorded
      together with their likelihood score
    - sep (str): separator between fields in the output file
    - clean_tree (bool): if True the created tree file is deleted after Raxml runs
    Output: Explored trees written in out_file_path
    """
    out_file = open(out_file_path, "w")
    current_tree = first_tree
    _tree_file = "tmp_0_.tree"
    current_score = raxml_loss(
        fasta_path, current_tree, DNA_model,
        tree_folder_path, outfile=_tree_file
    )
    _clean_tree(tree_folder_path, _tree_file, clean_tree)
    _write_tree(
        current_tree, current_score, "START", out_file, sep=sep
    )
    stop = False
    patience_counter = max_patience
    iteration_counter = 1
    while not stop:
        current_tree = current_tree.reorder_leaves(rng)
        best_ngb_tree,best_ngb_score = best_ML_neighbour(
            current_tree, current_score,
            fasta_path, DNA_model, tree_folder_path,
            iteration_counter, clean_tree
        )
        if current_score - best_ngb_score > tol:
            current_tree = best_ngb_tree
            current_score = best_ngb_score
            patience_counter = max_patience
            prefix = "NEIGHBOUR"
        else:
            patience_counter -= 1
            prefix = "REORDER"
        stop = (
            (patience_counter==0)
            or
            (
                (max_nb_iterations is not None)
                and
                (iteration_counter>max_nb_iterations)
            )
        )
        iteration_counter += 1
        _write_tree(
            current_tree, current_score, prefix, out_file, sep=sep
        )

def _random_walk(
        fasta_path, DNA_model, tree_folder_path,
        first_tree, rng,
        max_nb_iterations,
        out_file_path,
        sep,
        clean_tree
):
    """
    Performing a random walk in the tree space for a given number of iterations

    Input:
    - fasta_path (str): pah to FASTA file
    - DNA_model (str): DNA evolution model to use in Raxml
    - tree_folder_path (str): pah to folder where Newick trees are written
    - first_tree (TreeVec): starting tree
    - rng (Generator): numpy random number generator
    - max_nb_iterations (int): maximum number of iterations (steps of the walk)
    - out_file_path (str): path of the file where all explored trees are recorded
      together with their likelihood score
    - sep (str): separator between fields in the output file
    - clean_tree (bool): if True the created tree file is deleted after Raxml runs
    Output: Explored trees written in out_file_path
    """
    out_file = open(out_file_path, "w")
    current_tree = first_tree
    current_score = None
    _tree_file = "tmp_0_.tree"
    current_score = raxml_loss(
        fasta_path, current_tree, DNA_model,
        tree_folder_path, outfile=_tree_file
    )
    _clean_tree(tree_folder_path, _tree_file, clean_tree)
    _write_tree(
        current_tree, current_score, "START", out_file, sep=sep
    )
    for i in range(max_nb_iterations):
        new_tree = current_tree.random_hop(rng, inplace=False)
        current_tree = new_tree.reorder_leaves(rng)
        _tree_file = f"tmp_0_{i}.tree"
        current_score = raxml_loss(
            fasta_path, current_tree, DNA_model,
            tree_folder_path, outfile=_tree_file
        )
        _clean_tree(tree_folder_path, _tree_file, clean_tree)
        _write_tree(
            current_tree, current_score, "RANDOM", out_file, sep=sep
        )

def hill_climbing(
        fasta_path, DNA_model, tree_folder_path,
        tol, max_patience, max_nb_iterations,
        random_seed,
        out_file_path
):
    # Read a FASTA file to record sequences (taxa) names
    leaves = [record.id for record in _read_fasta(fasta_path)]
    # Instantiate a random numbers generators
    rng = np.random.default_rng(random_seed)
    # Create a random tree to start from
    first_tree = _create_random_tree(leaves)
    _hill_climbing(
        fasta_path, DNA_model, tree_folder_path,
        first_tree, rng,
        tol, max_patience, max_nb_iterations,
        out_file_path,
        sep="\t",
        clean_tree=True
    )
