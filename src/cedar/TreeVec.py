"""
CEDAR: manipulating phylogenetic rooted trees representations as vectors
TreeVec class for vector representation of a tree
"""

__author__ = "Cedric Chauve"
__credits__ = ["Cedric Chauve", "Louxin Zhang"]
__license__ = "GPL"
__version__ = "1.2.0"
__maintainer__ = "Cedric Chauve"
__email__ = "cedric.chauve@sfu.ca"
__status__ = "Release"

import numpy as np
from ete3 import Tree
from cedar.LIS import LIS_len, LIS_seq
from cedar.utils import (
    __read_file,
    __write_file
)
from cedar.LeavesOrder import (
    str2order,
    order2str
)

# Separator between a label and a node name in a tree representation
SEP_NODE = ":"
# Separators between tree representation elements
SEP_VEC = ","

# Main Class
class TreeVec:
    """
    Vector representation of a tree.

    The topology of a tree with n leaves, ordered from 1 to n is encoded by a
    list of 2n integer labels
    - start with 1,
    - ends with n,
    - contains 2 occurrences of every integer in {1,...,n}
    - the first occurrence of i>1 appears before the second copy of i-1
    - the second occurrence of i>1 appears after the second occurrence of i-1
    - the first occurrence of i encodes an internal node
    - the second occurrence of i encodes a leaf
    The tree is augmented by a root labeld 1 and with a single child called the
    dummy root.

    Data structure: list([int,str,float,bool])
    - field 0 (int): label
    - field 1 (str): name of the node in the tree
    - field 2 (float): length of the branch to the parent;
      the root and the dummy root have a branch length equal to 0.0
    - field 3 (bool): True if second occurrence (leaf)
                      False if first occurrence (internal node)

    String encoding
    A tree vector representation can be written in format 1 or 2 and in compact or
    non-compact writing:
    - nodes are separated by SEP_VEC
    - format 1.non-compact: each node is written as label:name:dist
    - format 2.non-compact: each node is written as label:name
    - format 1.compact:
      - each internal node is written as label:name:dist
      - each leaf is written as dist
        to be decoded this requires a mapping idx2leaf (dict int -> str) that
        defines a total order on leaves and allows to recover the leaf name and label
        associated to positions in the vector encoding leaves
    - format 2.compact:
      - each internal node is written as label:name
      - each leaf is written as an empty string whose name and labels can be recovered
        from the mapping idx2leaf as described above
    """

    def __init__(
            self,
            treevec_vec=None,
            tree=None,
            newick_str=None,
            treevec_str=None,
            leaf2idx=None,
            idx2leaf=None,
            format_str=None,
            compact=None
    ):
        """
        Instantiate a vector representation for a tree on n leaves
        - If treevec_vec is not None, the vector is created using it as vector
        - If tree is not None, tree is a Tree object and the vector is created from it
          using leaf2idx
        - If newick_str is not None it is created from newick_str using idx2leaf and
          expected in Newick format=1
        - If treevec_str is not None it is created from treevec_str using idx2leaf and
          expected in format defined by format_str and compact
        - Otherwise an empty vector is created
        - leaf2idx (dict str -> int): leaf name to index in a total order on leaves
          (1-base)
        - idx2leaf (dict int -> str): reverse dictionary
        - format_str (int in [1,2])
        - compact (bool)
        """
        self.vector = []
        if treevec_vec is not None:
            self.vector = treevec_vec
        elif tree is not None:
            self.vector = self.tree2treevec(tree, leaf2idx=leaf2idx)
        elif newick_str is not None:
            self.vector = self.newick2treevec(newick_str, leaf2idx=leaf2idx)
        elif treevec_str is not None:
            self.vector = self.str2treevec(
                treevec_str, idx2leaf,
                format_str=format_str, compact=compact
            )

    def copy(self):
        """
        Creates a copy of self
        """
        def __copy_vector(v):
            """
            Creates a copy of vector v
            """
            return [node.copy() for node in v]
        return TreeVec(treevec_vec=__copy_vector(self.vector))

    def check_vector(self):
        """
        Check that a tree vector is a proper tree representation
        Ouput:
        - (bool): True if it is, false otherwise
        """
        v = self.vector
        n = int(len(v) / 2)
        # Positions in v of all labels
        positions = {i: [] for i in range(1,2*n+1)}
        for i in range(0,2*n):
            positions[v[i][0]].append(i)
        for x in range(2,n+1):
            if len(positions[x]) != 2: return False
            elif (x<n) and positions[x][0] > positions[x-1][1]: return False
            elif (x<n) and positions[x][1] < positions[x-1][1]: return False
        return True

    def extract_leaves_order(self):
        """
        Given a vector representation of length 2n, computes
        - leaf2idx (dict str -> int): leaf name to index in a total order on leaves
          indexed in 1-base
        - idx2leaf (dict int -> str): reverse dictionary
        numbers and the implicit order on leaves
        """
        leaf2idx,idx2leaf,j = {},{},1
        for node in self.vector:
            leaf = node[3]
            if leaf:
                name = node[1]
                leaf2idx[name],idx2leaf[j] = j,name
                j += 1
        return leaf2idx,idx2leaf

    def treevec2tree(self):
        """
        Given a tree vector representation, compute a Tree object
        Ouput:
        - (Tree) Tree object
        """
        v = self.vector
        n = int(len(v) / 2)
        # Computing the label and name of each node
        __label, __name, __dist = {}, {}, {}
        for i in range(0,2*n):
            [label,name,dist,leaf] = v[i]
            __label[i] = label
            if leaf:
                __name[label] = name
                __dist[label] = dist
            else:
                __name[n+label] = name
                __dist[n+label] = dist
        # Decoding into edges
        edges = {i: [] for i in range(1,2*n+1)}
        o = {j: (2 if v[j][3] else 1) for j in range(0,2*n)}
        for j in range(0,2*n-1):
            if o[j] == 1 and o[j+1] == 1:
                edges[n+__label[j]].append(n+__label[j+1])
            elif o[j] == 1 and o[j+1] == 2:
                edges[n+__label[j]].append(__label[j+1])
            elif o[j] == 2 and o[j+1] == 1:
                k = j+2
                while o[k] == 1: k += 1
                edges[n+__label[k]].append(n+__label[j+1])
            elif o[j] == 2 and o[j+1] == 2:
                edges[n+__label[j+1]].append(__label[j+1])
        # Creating Tree structure
        nodes = {i: Tree(name=__name[i],dist=__dist[i]) for i in range(1,2*n+1)}
        for node_from_idx,nodes_to_idx in edges.items():
            node_from = nodes[node_from_idx]
            for node_to_idx in nodes_to_idx:
                node_from.add_child(nodes[node_to_idx])
        root = nodes[n+1].children[0]
        return root

    def tree2treevec(self, tree, leaf2idx=None):
        """
        Compute the vector representation of a tree with n leaves from a Tree objet
        Input:
        - t: Tree object with features "name" and "dist" (branch length)
        - leaf2idx: dict(str -> int) leaf name to leaf label
        if None: leaf labels added during a postorder traversal in order of visit.
        """
        # Adding a root labeled 1 and named ""
        T = Tree(name="")
        T.add_feature("label", 1)
        T.add_child(tree)
        # Labeling nodes
        label = 1
        for node in tree.traverse("postorder"):
            if node.is_leaf() and leaf2idx is not None:
                node.add_feature("label", leaf2idx[node.name])
                node.add_feature("min_label", node.label)
            elif node.is_leaf():
                node.add_feature("label", label)
                node.add_feature("min_label", node.label)
                label += 1
            else:
                children_min_label = [child.min_label for child in node.children]
                node.add_feature("min_label", min(children_min_label))
                node.add_feature("label", max(children_min_label))
        # Computing a dictionary from label to leaf
        label2leaf = {}
        for node in T.traverse("postorder"):
            if node.is_leaf():
                label2leaf[node.label] = node
                n = len(label2leaf.keys())
                # Computing paths from leaves to the internal node of same label
        paths = {}
        for i in range(1,n+1):
            paths[i] = []
            node = label2leaf[i].up
            while node.label != i:
                paths[i].append([node.label, node.name, node.dist, False])
                node = node.up
        # Concatenating reversed paths
        v = [[1, T.name, 0.0, False]]
        for i in range(1,n+1):
            leaf = label2leaf[i]
            v += paths[i][::-1] + [[i,leaf.name,leaf.dist,True]]
        return v

    def treevec2str(self, format_str=1, compact=True):
        """
        Transform a tree vector representation into a string in format
        defined by format_str (1 or 2) and compact (True or False)
        """
        sep_vec,sep_node = SEP_VEC,SEP_NODE,
        _,idx2leaf = self.extract_leaves_order()
        out_str = []
        for node in self.vector:
            [label,name,dist,leaf] = node
            if format_str == 1 and (not compact):
                out_str.append(f"{label}{sep_node}{name}{sep_node}{dist}")
            elif format_str == 2 and (not compact):
                out_str.append(f"{label}{sep_node}{name}")
            elif format_str == 1 and compact and (not leaf):
                out_str.append(f"{label}{sep_node}{name}{sep_node}{dist}")
            elif format_str == 1 and compact and leaf:
                out_str.append(f"{dist}")
            elif format_str == 2 and compact and (not leaf):
                out_str.append(f"{label}{sep_node}{name}")
            elif format_str == 2 and compact and leaf:
                out_str.append("")
        return f"{sep_vec.join(out_str)};"

    def str2treevec(self, s, idx2leaf, format_str=1, compact=True):
        """
        Reads a tree vector representation from a string s
        Input:
        - s (str): string encoding of a tree vector representation
        - format_str (int in [1,2]): expected encoding format
        - compact (bool): if True, compact writing
        """
        sep_vec,sep_node = SEP_VEC,SEP_NODE
        s1 = s.rstrip()[:-1].split(sep_vec)
        n = int(len(s1)/2)
        occurrences = {i: 0 for i in range(1,n+1)}
        v,j = [],1
        for node_str in s1:
            node = node_str.split(sep_node)
            if format_str == 1 and (not compact):
                label,name,dist = int(node[0]),node[1],float(node[2])
            elif format_str == 2 and (not compact):
                label,name,dist = int(node[0]),node[1],1.0
            elif format_str == 1 and compact and len(node) == 3:
                label,name,dist = int(node[0]),node[1],float(node[2])
            elif format_str == 1 and compact and len(node) == 1:
                label,name,dist = j,idx2leaf[j],float(node[0])
                j += 1
            elif format_str == 2 and compact and len(node) == 2:
                label,name,dist = int(node[0]),node[1],1.0
            elif format_str == 2 and compact and len(node) == 1:
                label,name,dist = j,idx2leaf[j],1.0
                j += 1
            leaf = (False if occurrences[label]==0 else True)
            v.append([label,name,dist,leaf])
            occurrences[label] += 1
        return v

    def treevec2newick(self, newick_format=0):
        """
        Transform a vector representation into a Newick string in newick_format
        """
        tree = self.treevec2tree()
        tree1 = tree.get_tree_root().children[0]
        return tree1.write(format=newick_format)

    def newick2treevec(self, newick_str, leaf2idx=None):
        """
        Instantiate a vector representation from a Newick string
        - leaf2idx: dict(str -> int) leaf name to leaf label
        if None: leaf labels added during a postorder traversal in order of visit.
        """
        tree = Tree(newick_str, format=1)
        return self.tree2treevec(tree, leaf2idx=leaf2idx)

    def reorder_leaves(self, rng):
        """
        Reorder the leaves of the current tree randomly
        Input:
        - rng (Generator): numpy random numbers generator
        Output:
        - TreeVec
        """
        leaf2idx,_ = self.extract_leaves_order()
        leaves_idx = list(leaf2idx.values())
        leaves_labels = list(leaf2idx.keys())
        rng.shuffle(leaves_idx)
        leaf2idx = {
            leaves_labels[i]: leaves_idx[i]
            for i in range(len(leaves_idx))
        }
        new_tree = TreeVec(tree=self.treevec2tree(), leaf2idx=leaf2idx)
        return new_tree

    """ Hop-related functions
    A hop is an SPR, so a subtree is pruned then regrafted.
    Not every SPR is a hop as the following restriction applies to hop:
    if a is the smallest (in the leaves order associated to a tree) leaf in the pruned subtree,
    it can only be regrafted on an edge on one of the LTSv paths for leaves 1 to a-1.

    In terms of vector encoding, a hop is as follows:
    - pruning: removing the first occurrence of some label a > 1: this element corresponds to the
      unique internal node labelled by a and the pruned subtree is rooted at its child on the path
      to the leaf a; 
    - regrafting: inserting it anywhere between the first occurrence of 1 and the second occurrence 
      of a-1; if it is inserted before an element labelled by b, then the pruned subtree is regrafted
      on the branch from the corresponding node labelled b to its parent.
    Branch lengths are modified as follows: the branch created by the pruning has length equal to 
    the sum of the two branches defining it, while the branch where the prune subtree is regrafted 
    leads to wo branches with length half the length of the initial branch.
    
    In terms of implementation, pruning is implemented using the pop function and regrafting using 
    the insert function, both functions applying on python lists.

    Internally, a hop is encoded with a pair (i,j) where:
    - i>0 is the position in self.vector of the popped element (labeled a);
    - j is the position in self.vector where the popped element will be inserted;
    both positions are in 0-index.

    By definition of the vector encoding, if the second ocurrence of a-1 is in position k, then
    i<k and j must satisfy j<k. So j must be in [1,k]. Pairs that do not define a new tree and
    must be excluded are
    - j=i so (i,i)
    - j=i+1 so (i,i+1)
    - j=i+2 if i+1 is not a leaf as (i,i+2) is equivalent to (i+1,i)
 
    """

    def __hop_update_dist(self, v, i, j, x, y, z):
        """
        Update the branch lengths of a vector representation after having done a hop
        i: position in vector of node whose branch to parent was subdivided by the pruning
           (pruned subree root sibling)
        j: position in vector of node whose branch to parent contains the regrafting
           (regrafting node)
           implies that j-1 is the position of its parent, the new node created by regrafing
        x: length of branch from pruned subtree root sibling to parent
        y: length of branch from pruned subtree root parent to parent
        z: length of branch from regrafting node to parent
        """
        v[i][2] = x+y
        v[j-1][2] = z/2
        v[j][2] = z/2

    def __hop_init(self, i, j):
        """
        Create a new TreeVec object by a hop (i,j)
        Input:
        - i,j (int): hop coordinates
        WARNING: does not check i and j define a valid hop
        """
        v = self.copy().vector
        new_v = []
        x,y,z = v[i+1][2],v[i][2],v[j][2]
        if j <= i-1:
            new_v = v[0:j]+[v[i]]+v[j:i]+v[i+1:]
            self.__hop_update_dist(new_v,i+1,j+1,x,y,z)
        elif j > i+1:
            new_v = v[0:i]+v[i+1:j]+[v[i]]+v[j:]
            self.__hop_update_dist(new_v,i,j,x,y,z)

        return TreeVec(treevec_vec=new_v)

    def __hop_inplace(self, i, j):
        """
        Modify self.vector by a hop (i,j)
        Input:
        - i,j (int): hop coordinates
        WARNING: does not check i and j define a hop leading to a valid tree
        representation
        """
        v = self.vector
        pruned_subtree = v.pop(i)
        x,y,z = v[i+1][2],v[i][2],v[j][2]
        if j <= i-1:
            v.insert(j,pruned_subtree)
            self.__hop_update_dist(v,i+1,j+1,x,y,z)
        elif j > i+1:
            v.insert(j-1,pruned_subtree)
            self.__hop_update_dist(v,i,j,x,y,z)

    def hop(self, i, j, inplace=False):
        """
        Implement the hop rearrangement (i,j)
        Input:
        - i,j (int): hop coordinates
        - inplace (bool): if True modifies the current object oherwise returns a new object
        WARNING: does not check i and j define a hop leading to a valid tree
        representation
        """
        v = self.vector
        if not (i>0 and j>0 and (j<=i-1 or j>i+1) and (not v[i][3])):
            raise Exception("ERROR: improper HOP coordinates")
        if inplace:
            self.__hop_inplace(i, j)
            return None
        else:
            return self.__hop_init(i, j)

    def _compute_leaves_positions(self):
        """
        Computes a dictionary leafpos that records he position in the
        vector of the second occurrence of each label
        """
        v = self.vector
        leafpos = {}
        for k in range(1,len(v)):
            node = v[k]
            if node[3]:
                leafpos[node[0]] = k
        return leafpos

    def _compute_hops_range(self, i, k, size=False):
        """
        Compute the number/set of positions elment in position i can be moved to by a hop
        if size=True: count number of positions, otherwise cmpute the set of positions
        """
        range_j = list(range(1,k+1)) # 1 to exclude first position 0
        # excluding (i,i)
        range_j.remove(i)
        # excluding (i,i+1) if i+1 in the range of possible hops
        if i+1 <= k:
            range_j.remove(i+1)
        # excluding (i,i+2) if i+2 in the range of possible hops and v[i+1] not a leaf
        # as it is equivalent to (i+1,i) that will not be excluded
        if i+2<=k and not self.vector[i+1][3]:
            range_j.remove(i+2)
        return len(list(range_j)) if size else range_j

    def hop_neighbourhood_size(self):
        """
        Compute the size of the hop neighbourhood of a tree
        - Output: (int), neighbourhood size
        """
        leafpos = self._compute_leaves_positions()
        ngb_size = 0
        for i in range(1,len(self.vector)):
            moved_node = self.vector[i]
            if not moved_node[3]:
                # k: position of second occurrence of a-1 if moved_node has label a
                k = leafpos[moved_node[0]-1]
                ngb_size += self._compute_hops_range(i, k, size=True)
        return ngb_size

    def hop_neighbourhood(self, export_list=False):
        """
        Compute the hop neighbourhood of a tree
        - Output:
        if export_list is False: list(TreeVec)
        else: list((int,int)), list of HOP defining the neighbourhood
        """
        leafpos = self._compute_leaves_positions()
        ngb = []
        for i in range(1,len(self.vector)):
            moved_node = self.vector[i]
            if not moved_node[3]:
                k = leafpos[moved_node[0]-1]
                range_j = list(self._compute_hops_range(i, k, size=False))
                ngb += [(i,j) for j in range_j] if export_list else [
                    self.hop(i,j, inplace=False) for j in range_j
                ]
        return ngb

    def random_hop(self, rng, inplace=False):
        """
        Computes a new TreeVec representing a tree differing from self by a single random hop
        under the uniform distribution.
        Input:
        - rng: random numbers generator
        - inplace(bool): If True modifies the curren object, otherwise returns a new object
        Output:
        - TreeVec
        """
        # Number of possible random hops
        hop_ngb_size = self.hop_neighbourhood_size()
        # Index of the hop to apply
        hop_rank = rng.integers(1,high=hop_ngb_size+1)
        leafpos = self._compute_leaves_positions()
        # Number of possible hops already considered
        nb_hops = 0
        for i in range(1,len(self.vector)):
            candidate_moved_node = self.vector[i]
            if not candidate_moved_node[3]:
                k = leafpos[candidate_moved_node[0]-1]
                range_j = list(self._compute_hops_range(i, k, size=False))
                nb_possible_hops = len(range_j)
                if  nb_possible_hops >= hop_rank:
                    j = range_j[hop_rank-1]
                    return self.hop(i, j, inplace=inplace)
                else:
                    nb_hops += nb_possible_hops
                    hop_rank -= nb_possible_hops

    # Hop-similarity related functions

    def hop_similarity(self, t2, compute_seq=False):
        """
        Compute the hop smilarity to another tree representations
        Input:
        - t2 (TreeVec)
        assumption: both are on the same leaves order (asserted)
        - compute_seq (bool): if True, returns an actual LCS, if False, returns
          the similarity value
        Output:
        - compute_seq=False: (int) in [0,n]
        - compute_seq=True: list((int,bool)) list of (integers,True if leaf)
          encoding the LCS between v1 and v2
        """

        def __relabel_segment(segment1, segment2):
            """
            Relabel the labels of segment1 increasingly from 0
            and the labels of segment2 according to the relabeling of
            segment1, excluding labels not present in segment1
            """
            # map1[x] = position of label x in segment1
            map1 = {segment1[i1]: i1 for i1 in range(0,len(segment1))}
            # Relabeling segment2 according to __map1,
            # excluding labels not in segment1
            relabeled_segment2 = []
            for i2 in range(0,len(segment2)):
                if segment2[i2] in map1.keys():
                    relabeled_segment2.append(map1[segment2[i2]])
            return relabeled_segment2

        v1,v2 = self.vector,t2.vector
        n = int(len(v1)/2)
        second_occ_order = [
            v1[i][0]
            for i in range(0,2*n) if v1[i][3]
        ]
        # Compute a list of pairs of subsequences to compare pairwise
        # boundaries1[i] = [j,k]: boundaries of the segment of internal nodes
        # in v1 before leaf i+1 similar for boundaries2 and v2
        # if j>k: empty segment
        boundaries = {1: [], 2: []}
        i1,i2 = 0,0
        for i in range(0,2*n):
            leaf1,leaf2 = v1[i][3],v2[i][3]
            if leaf1:
                boundaries[1].append([i1,i-1])
                i1 = i+1
            if leaf2:
                boundaries[2].append([i2,i-1])
                i2 = i+1
        # Computes an LCS for each pair of segments using an LIS algorithm
        lcs_len,lcs_seq = 0,[]
        for j in range(0,n):
            b1_start,b1_end = boundaries[1][j][0], boundaries[1][j][1]
            b2_start,b2_end = boundaries[2][j][0], boundaries[2][j][1]
            # Checking that both segments are non-empty (otherwise, no LCS)
            if (b1_end>=b1_start) and (b2_end>=b2_start):
                # Segments of v1 and v2 to consider
                __segment1 = [v1[k][0] for k in range(b1_start, b1_end+1)]
                __segment2 = [v2[k][0] for k in range(b2_start, b2_end+1)]
                # Relabeling __segment2 according to __map1,
                # excluding labels not in __segment1
                segment2 = __relabel_segment(__segment1, __segment2)
                # Computing an LIS in segment2
                if compute_seq: lcs_seq += [
                        (__segment1[i2],False)
                        for i2 in LIS_seq(segment2)
                ]
                else: lcs_len += LIS_len(segment2)
            lcs_seq += [(second_occ_order[j],True)]
        return (lcs_seq if compute_seq else lcs_len)

    def hop_next(self, t2, inplace=False):
        """
        Computes a new TreeVec representing a tree differing
        from self=t1 by a single hop on the hop-path from t1 to t2
        Input:
        - t2 (TreeVec)
        Output:
        - TreeVec
        """
        v1,v2 = self.vector,t2.vector
        alignment = self.hop_similarity(t2, compute_seq=True)
        # If v1==v2: return None
        if len(alignment) == len(v1):
            return None
        def __equal(a,s):
            return (a[0]==s[0] and a[1]==s[3])
        # Find first misaligned element in v2
        k = 0
        l = len(v2)
        while __equal(alignment[k],v2[k]): k += 1
        x,y = (v2[k][0],v2[k][3]),alignment[k-1]
        # Find x,y in v1
        i = 0
        while not __equal(x,v1[i]): i+=1
        j = 0
        while not __equal(y,v1[j]): j+=1
        # Move x after y
        return self.hop(i,j+1, inplace=inplace)

def read_TreeVec_file(in_TreeVec_file):
    # Reading the TreeVec trees
    in_TreeVec_trees = __read_file(in_TreeVec_file)
    # Determining the leaves order
    leaf2idx,idx2leaf = str2order(in_TreeVec_trees[0].split()[1])
    # Creating TreeVec objects
    TreeVec_trees = []
    for TreeVec_tree in in_TreeVec_trees[1:]:
        TreeVec_trees.append(
            TreeVec(
                treevec_str=TreeVec_tree, idx2leaf=idx2leaf, format_str=1, compact=True
            )
        )
    return TreeVec_trees,leaf2idx,idx2leaf

def write_TreeVec_file(in_TreeVec_trees, idx2leaf, out_TreeVec_file):
    out_str = [f"#order {order2str(idx2leaf)}"]
    for TreeVec_tree in in_TreeVec_trees:
        out_str.append(
            TreeVec_tree.treevec2str(format_str=1, compact=True)
        )
    __write_file(out_str, out_TreeVec_file)

def get_nb_taxa(in_TreeVec_tree):
    leaf2idx,_ = in_TreeVec_tree.extract_leaves_order()
    nb_taxa = len(list(leaf2idx.keys()))
    return nb_taxa

def random_leaves_order(in_TreeVec_file, nb_orders=1, in_seed=0, out_file_prefix="leaves_order"):
    """
    Generates nb_orders random leaves orders and write them in files
    {out_prefix_file}_{nb order}.txt
    """
    _,_,idx2leaf = read_TreeVec_file(in_TreeVec_file)
    nb_leaves = len(idx2leaf.keys())
    idx = np.array(list(idx2leaf.values()))
    rng = np.random.default_rng(in_seed)
    for i in range(1,nb_orders+1):
        rng.shuffle(idx)
        _idx2leaf = {j: idx[j-1] for j in range(1,nb_leaves+1)}
        with open(f"{out_file_prefix}_{i}.txt", "w") as out_file:
            out_file.write(f"{order2str(_idx2leaf)}\n")
