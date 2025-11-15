"""
CEDAR: manipulating phylogenetic rooted trees representations as vectors
"""

__author__ = "Cedric Chauve"
__credits__ = ["Cedric Chauve", "Louxin Zhang"]
__license__ = "GPL"
__version__ = "1.2.0"
__maintainer__ = "Cedric Chauve"
__email__ = "cedric.chauve@sfu.ca"
__status__ = "Release"

import os
import argparse

from cedar.TreeVec import (
    random_leaves_order
)
from cedar.Newick import (
    convert_Newick2TreeVec,
    convert_TreeVec2Newick
)
from cedar.HOP import (
    hop_similarity,
    hop_neighbourhood_size,
    hop_neighbourhood,
    hop_path,
    hop_random
)
from cedar.TreeSpace import (
    hill_climbing
)
from cedar.GelmanRubin import (
    GelmanRubin
)

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
    orders.add_argument("--seed", type=int, default=0, help="[OPTIONAL] Random generator seed")

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

    # Computing a random HOP for each tree of a file
    hop_random = subparsers.add_parser("HOP_random", help="Create new trees by applying a rando HOP to each input trees")
    hop_random.set_defaults(cmd="HOP_random")
    hop_random.add_argument("--input_file", type=str, help="Input CEDAR file")
    hop_random.add_argument("--output_file", type=str, help="Output CEDAR file")
    hop_random.add_argument("--seed", type=int, default=0, help="[OPTIONAL] Random number generator seed")

    # Hill-Climbing
    hop_hc = subparsers.add_parser("HOP_hc", help="Hill-Climbing heuristic")
    hop_hc.set_defaults(cmd="HOP_hc")
    hop_hc.add_argument("--fasta_path", type=str, help="Input FASTA file")
    hop_hc.add_argument("--DNA_model", type=str, help="Input DNA model for Raxml")
    hop_hc.add_argument("--tree_folder_path", type=str, help="Path to folder where trees are written for Raxml (not created)")
    hop_hc.add_argument("--tol", type=float, default=0.001, help="[OPTIONAL] Tolerance to dcide if a best tree has been found")
    hop_hc.add_argument("--max_patience", type=int, default=5, help="[OPTIONAL] Max number of leaves reordering steps if no best neighbour is found")
    hop_hc.add_argument("--max_nb_iterations", type=int, default=None, help="[OPTIONAL] Maximum number of iterations")
    hop_hc.add_argument("--seed", type=int, default=0, help="[OPTIONAL] Random number generator seed")
    hop_hc.add_argument("--out_file_path", type=str, help="Path to output file")

    # Gelman Rubin convegence statistics
    GR = subparsers.add_parser("GR", help="Gelman Rubin diagnostic values")
    GR.set_defaults(cmd="GR")
    GR.add_argument("--Newick_file_1", type=str, help="Input Newick file 1")
    GR.add_argument("--Newick_file_2", type=str, help="Input Newick file 2")
    GR.add_argument("--nb_trees", type=int, help="Number of trees (tail of Newick files to consider)")
    GR.add_argument("--nb_orders", type=int, default=5, help="[OPTIONAL] Number of random leaves orders")
    GR.add_argument("--seed", type=int, default=0, help="[OPTIONAL] Random number generator seed")
    GR.add_argument("--output_gr_file", type=str, help="Output file containing GR statistics")
    GR.add_argument("--output_orders_file", type=str, help="Output file containing leaves orders")

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
        in_seed=args.seed,
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

def HOP_random(args):
    hop_random(
        args.input_file,
        args.output_file,
        in_seed=args.seed
    )

def HOP_hc(args):
    hill_climbing(
        args.fasta_path, args.DNA_model, args.tree_folder_path,
        args.tol, args.max_patience, args.max_nb_iterations, args.seed,
        args.out_file_path
    )

def GR(args):
    GelmanRubin(
        args.Newick_file_1,
        args.Newick_file_2,
        args.nb_trees,
        args.nb_orders,
        args.seed,
        args.output_gr_file,
        args.output_orders_file
    )

def main():
    args = _parse_arguments()

    match args.cmd:
        case "fromNewick":
            fromNewick(args)
        case "toNewick":
            toNewick(args)
        case "orders":
            orders(args)
        case "HOP_sim":
            HOP_sim(args)
        case "HOP_ngb1":
            HOP_ngb1(args)
        case "HOP_ngb2":
            HOP_ngb2(args)
        case "HOP_path":
            HOP_path(args)
        case "HOP_random":
            HOP_random(args)
        case "HOP_hc":
            HOP_hc(args)
        case "GR":
            GR(args)
        case _:
            raise Exception("ERROR: Unknown command")


if __name__ == "__main__":
    main()
