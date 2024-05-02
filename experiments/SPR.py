"""
Implementing SPR tree rearrangement in ETE3
"""

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
    _ = B.detach()
    if not Y.is_root():
        X.dist += Y.dist
        Y.delete()
    else:
        t.name = X.name
        t.children = X.children
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
    __nodes = [node.name for node in t.traverse("postorder")]
    if a == b:
        return False
    B = t.search_nodes(name=b)[0]
    A = t.search_nodes(name=a)[0]
    if A.is_root():
        # f"{a} = root"
        return False
    ancestor = t.get_common_ancestor(a,b)
    if ancestor.name == b:
        # f"{b} ancestor of {a}"
        return False
    Y = B.up
    if Y.name == a:        
        # f"{b} child of {a}"
        return False
    if A in Y.get_children():
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


# #t = Tree('((((H:1.0,K:1.0)D:1.0,(F:1.0,I:1.0)G:1.0)B:1.0,E:1.0)A:1.0,((L:1.0,(N:1.0,Q:1.0)O:1.0)J:1.0,(P:1.0,S:1.0)M:1.0)C:1.0)X;', format=1)

# def print_tree(t):
#     print(t.get_ascii(show_internal=True, attributes=["dist","name"]))

# x_random = True
            
# if not x_random:
#     print_tree(t)
#     spr_inplace(t, "N", "S")
#     print_tree(t)
#     spr_inplace(t, "M", "F")
#     print_tree(t)
#     spr_inplace(t, "A", "S")
#     print_tree(t)
#     spr_inplace(t, "J", "M")
#     print_tree(t)
#     a,b = "Q","X"
#     if check_spr(t, a, b):
#         spr_inplace(t, a, b)
#     else:
#         print(f"SPR {b} pruned and regrafted above {a} invalid")
# else:
#     for i in range(0,5):
#         nb_tries = random_spr(t)
#         print(f"{nb_tries} tries")
#         print_tree(t)
