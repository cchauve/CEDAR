"""
Implementing SPR tree rearrangement in ETE3
"""
import collections

from ete3 import Tree
import random

def spr_inplace(t, a, b):
    """
    Modifies a tree t by pruning the subtree rooted at node named b
    above node names a
    - t (Tree)
    - a,b (str,str)
    """
    # A, B
    B = t.search_nodes(name=b)[0]
    A = t.search_nodes(name=a)[0]
    aa = A.dist
    # B parent
    Y = B.up
    y = Y.name
    # B sister
    X = [x for x in Y.get_children() if x != B][0]
    # A parent
    U = A.up
    # Prune subtree rooted at node B
    _ = B.detach()#__detach(t,B)
    if not Y.is_root():
        X.dist += Y.dist
        Y.delete()#__delete(t,Y)
    else:
        t.name = X.name
        t.children = X.children
        for v in t.children: v.up = t
    # Regraft above node A
    U.add_child(name=y, dist=aa/2.0)
    A.dist = aa/2.0
    AA = A.detach()
    Y = t.search_nodes(name=y)[0]
    Y.add_child(B)
    Y.add_child(AA)

def check_spr(t, a, b):
    """
    Check if the SPR pruning the subtree at node named b and
    regrafting above node b is correct and modifies the topology
    of t
    - t (Tree)
    - a,b (str,str)
    """
    if a == b:
        return False
    B = t.search_nodes(name=b)[0]
    A = t.search_nodes(name=a)[0]
    if A.is_root():
        # f"{a} = root"
        return False
    if B.is_root():
        # f"{b} = root"
        return False
    ancestor = t.get_common_ancestor(a,b)
    if ancestor.name == b:
        # f"{b} ancestor of {a}"
        return False
    Y = B.up
    if Y.name == a:        
        # f"{b} child of {a}"
        return False
    if a in [V.name for V in Y.get_children()]:
        # f"{b} sister of {a}"
        return False
    return True
    
def random_spr(t):
    """
    Generates random SPR coordinates until a valid one is found
    and do the SPR
    - t (Tree)
    """
    __nodes = [node.name for node in t.traverse("postorder")]
    nb_nodes = len(__nodes)
    nodes = {i: __nodes[i] for i in range(0,nb_nodes)}
    nb_tries = 1
    while True:
        i = random.randint(0,nb_nodes-1)
        j = random.randint(0,nb_nodes-1)
        a,b = nodes[i],nodes[j]
        if check_spr(t, a, b):
            spr_inplace(t, a, b)
            return nb_tries
        else:
            nb_tries += 1
