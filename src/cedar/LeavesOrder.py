"""
CEDAR: manipulating phylogenetic rooted trees representations as vectors
Auxiliary functions for manipulating leaves orders
"""

__author__ = "Cedric Chauve"
__credits__ = ["Cedric Chauve", "Louxin Zhang"]
__license__ = "GPL"
__version__ = "1.2.0"
__maintainer__ = "Cedric Chauve"
__email__ = "cedric.chauve@sfu.ca"
__status__ = "Release"

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


