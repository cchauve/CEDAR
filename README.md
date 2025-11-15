# CEDAR: Encoding phylogenetic trees as vectors

CEDAR is a program aimed at manipulating rooted phylogenetic trees encoded as vecors, as described in the paper <a href="https://doi.org/10.1098/rstb.2024.0226">A Vector Representation for Phylogenetic Trees</a>.

**WARNING.** CEDAR is still in development and comes with no warranty.

CEDAR can be used in command line (described below) or wihin python progams using the class `TreeVec` implemented in the file [TreeVec.py](src/cedar/TreeVec.py).

*Dependencies:* <a href="https://numpy.org/">numpy</a> and <a href="http://etetoolkit.org/docs/latest/index.html">ete3</a>; <a href="https://github.com/amkozlov/raxml-ng">RAxML-NG</a> for the hill-climbing tree space exploration heuristic.

The folder [example](example/) contains an example of using CEDAR through command-line.

The folder [experiments](experiments) contains the code to reproduce the experiments described in the  paper <a href="https://doi.org/10.1098/rstb.2024.0226">A Vector Representation for Phylogenetic Trees</a>.

The folder [experiments_hc](experiments_hc) contains the code to reproduce the experiments using the hill-climbing tree space exploration heuristic implemented in CEDAR.  

The folder [experiments_gr](experiments_gr) contains the code to reproduce the experiments using the Gelman-Rubin MCMC convergence statstics implemented in CEDAR.

## Installation

This package can be installed via `pip`. (Missing a released .whl for this 
atm, then pip install https://github.com/cchauve/CEDAR/releases/download/...
/*.whl)

To build from source, you have to first clone the repository.
For an installation its enough to run `pip install .` after cloning.

You can build a `.whl` via `python -m build` (requires `build` package).
The wheel can then be installed via `pip install dist/*.whl`.

## Vector encoding of phylogenetic trees

CEDAR implements a vector representation of rooted phylogenetic trees. 
Consider the set of all rooted binary phylogenetic trees with totally ordered leaves, labelled $1,\dots,n$. 
We augment all trees with a dummy root connected to the root of $T$ by a branch of length $0.0$.

A vector representation of a tree is a vector $\mathbf{v}=(v_1,\dots,v_{2n})$ of $2n$ entries that satisfies the following constraints:
- $v_1=1$;
- $\forall i \in \{1,\dots,n\}$, $i$ appears exactly twice in $v$;
- $\forall i \in \{2,\dots,n\}$, the first occurrence of $i$ in $v$ appears before the second occurrence of $i-1$;
- $\forall i \in \{2,\dots,n\}$, the second occurrence of $i$ in $v$ appears after the second occurrence of $i-1$.

There is a one-to-one correspondence between the set of rooted binary phylogenetic trees with leaves labelled $1,\dots,n$ and the set of vector representations on $\{1,\dots,n\}$ defined above.
In a vector representation, the second occurrence of $i$ encodes the leaf labelled $i$ and the first occurrence of $i$ encodes an internal node on the path from the leaf $i$ to the root of the tree.

Given a rooted binary phylogenetic tree $T$ on $n$ leaves, the vector representation of $T$ is written as follows in a CEDAR tree:
- the entries of the vector are separated by commas (`,`);
- if the entry $i$ corresponds to an internal node, it is written as `i:N:b` where `N` is the name of the node and `b` the length of the branch to its parent;
- if the entry $i$ corresponds to a leaf it is written as `b` where here too `b` is the length of the branch to its parent.  
Note that the label $i$ of leaves does not need to be written as it is implicit that the $i^{th}$ leaf is labelled by $i$.
To recover the names of leaves, a CEDAR file always starts with a first line that contains the order on the leaves used to define the vector representation.

For example, a CEDAR file encoding the two trees 
```
((b:1.0,d:2.0)e:2.0,(a:1.0,c:1.0)g:1.0)f;
(((b:1.0,d:2.0)x:2.0,a:1.0)y:1.0,c:1.0)z;
```
where the first tree has internal nodes named `e,f,g` and the second tree has internal nodes named `x,y,z`, and where leaves are ordered as
`1=a,2=b,3=c,4=d`, is
```
#a,b,c,d
1::0.0,2:f:0.0,3:g:1.0,1.0,4:e:2.0,1.0,1.0,2.0;
1::0.0,3:z:0.0,2:y:1.0,1.0,4:x:2.0,1.0,1.0,2.0;
```
corresponding to vectors $(1,2,3,1,4,2,3,4)$ and $(1,3,2,1,4,2,3,4)$.

## Command-line usage

The command-line script is [CEDAR.py](src/cedar/CEDAR.py), and allows to perform the following tasks:
-  Converting trees written in Newick format into the CEDAR format:
   ```
   python src/CEDAR.py fromNewick --input_file Newick_file --output_file CEDAR_file [--order_file order]
   ```
   Encodes the files from `Newick_file` in the CEDAR format in file `CEDAR_file`.
   If the file `order` is provided, it is used to determine the taxa order used for the encoding;
   otherwise, the taxa order is defined by a postorder traversal of the first tree in `Newick_file`.
   The format of a taxa order file is a single line, starting conaining the list of all
   leaves, in the desired order, separated by commas.

   **Assumption**: all trees in `Newick_file` are rooted phylogenetic trees on the same set of taxa.

- Creating random taxa orders:
  ```
  python src/CEDAR.py orders --input_file CEDAR_file --output_dir out_dir [--nb_orders N] [--output_prefix out_pref] [--seed random_number_generator_seed]
  ```
  Creates `N` random taxa orders in files `out_dir/out_pref_I.txt` for `I` from `1` to `N`.
  Parameter `random_number_generator_seed` is the seed of the random generator and has default value `0`.
  Default values: `N=1`, `out_pref=CEDAR_random_order`, `seed=0`.

  **Assumption**: in all commands, all trees in a `CEDAR_file` are rooted phylogenetic trees on the same set of taxa
  and having been encoded into vector using the same order on taxa, written in the first line of the file.

- Converting trees written in CEDAR format into the Newick format:
  ```
  python src/CEDAR.py toNewick --input_file CEDAR_file --output_file Newick_file [--format Newick_format]
  ```
  The optional `format` argument defines the Newick format of the created Newick file (default: `1`).
  
- Computing the HOP similarity between trees:
  ```
  python src/CEDAR.py HOP_sim --input_file CEDAR_file --output_file similarity_file [--mode distance_mode]
  ```
  Creates a file that records the HOP similarity between the trees in CEDAR_file.
  - If `mode=pairwise` the HOP similarity is computed between all pairs of trees.
  - If `mode=sequence` the HOP similarity is computed between successive pairs of trees.
  - If `mode=first` the HOP similarity is computed between the first tree and all other trees.
  - The default mode is `sequence`.
    
- Computing the HOP neighbourhood size:
  ```
  python src/CEDAR.py HOP_ngb1 --input_file CEDAR_file --output_file size_file
  ```
  Computes the size of the HOP neighbourhood for every tree in `CEDAR_file`

- Computing the HOP neighbourhoods:
  ```
  python src/CEDAR.py HOP_ngb2 --input_file CEDAR_file --output_dir out_dir [--output_prefix out_pref]
  ```
  Computes the trees in the HOP neighbourhood for every tree in `CEDAR_file`. The HOP neighbourhood for
  tree `i` is in the file `out_dir/out_pref.i.vec`. The default value for `out_pref` is `CEDAR_HOP_neighbourhood`.

- Computing HOP paths between successive trees:
  ```
  python src/CEDAR.py HOP_path --input_file CEDAR_file --output_file CEDAR_path_file --mode sequence
  ```
  Creates a new CEDAR file where a sequence of trees, each differing from the previous one by a single
  HOP, is inserted between all pairs of successive trees in the input file `CEDAR_file`.

- Peforming a random HOP to each tree in a file
  ```
  python src/CEDAR.py HOP_random --input_file CEDAR_file --output_file CEDAR_path_file [--seed random_number_generator_seed]
  ```
  Creates a new CEDAR file containing on each line a tree differing by a single random HOP from the tree
  in the same line in the input file.
  Parameter `random_number_generator_seed` is the seed of the random generator and has default value `0`.

- Hill-climbing heuristic exploration of the tree space using HOPs:  
  starting from a random tree, the heuristic iterates the following steps
  - reorder randomly the leaves of the current tree
  - compute the likelihood of all trees in the HOP neighbourhood of the current tree using <a hef="https://github.com/amkozlov/raxml-ng">RAxML-NG</a> 
  - select the highest likelihood tree  
  - if its likelihood is within a given tolerance of the best tree so far:  
    - decrease a patience counter [patience step]  
  - otherwise:  
    - the best tree becomes the current tree  
    - the patience counter is set to max_patience
    - 
  until the maximum number of iterations is reached or the patience counter is 0.  
 
  The command `raxml-ng` is assumed to be available in the default path.

  ```
  python src/CEDAR.py HOP_hc --fasta_path FASTA_file --DNA_model DNA_model --tree_folder_path tree_folder_path \
  --out_file_path out_file \
  [--tol tolerance] [--max_patience] max_nb_patience_steps [--max_nb_iterations max_iter] \
  [--seed random_number_generator_seed]			     
  ```	
  All visited trees are recorded in the TSV file `out_file` in format `<S[TART,NEIGHBOUR,REORDER]>TAB<likelihood score>TAB<Newick tree>TAB<TreeVec tree>`  
  where `START` indicates the starting tree, `NEIGHBOUR` a step where a better neighbour was found, `REORDER` a patince step where
  the leaves order was randomly shuffled.
  Parameter `DNA_model` is a DNA model recognized by `raxml-ng`.  
  Parameter `random_number_generator_seed` is the random seed used to reorder randomly leaves in patience steps and has default value `0`.  
  Parameter `max_iter` limts the number of iterations and has value default `None` (no limit).  
  Parameter `max_nb_patience_steps` limits the number of times a patience steps is repeated consecutively and has default value `5`.  
  Parameter `tolerance` is used to determine if a better neighbour was found (if the likelihood of a neighbour tree is at least
  `tolerance` larger than the likelihood of the current tree) and has default value `0.001`.  
  The results of `raxml-ng` are stored in folder `tree_folder_path`.

- Gelman-Rubin MCMC convergence statistics:  
  Computes the Gelman-Rubin convegence test described in <a href="https://doi.org/10.1109/TCBB.2024.3457875">An Automated Convergence Diagnostic for Phylogenetic MCMC Analyses</a>
  using a distance between pais of trees defined as the minimum HOP distance over a specified number of random leaves orders.
  ```
  python src/CEDAR.py GR --Newick_file_1 input_file_1 --Newick_file_2 [--nb_trees nb_trees] \
  --nb_orders nb_orders [--seed random_number_generator_seed] \
  --output_gr_file output_file_1 --output_orders_file output_file_2
  ```	
  Parameters `input_file_1,input_file_2`: Newick files containing the trees of wo parallel MCMC runs.   
  Parameter `nb_trees`: number of trees to consider in the input file; the last trees of bth files are considered, previous trees being discarded as "burn-in" trees.    
  Parameter `nb_orders`: number of random leaves orders used to compute the distance between pairs of trees and has defaut value `5`; the same orders are used for all distance computations.   
  Parameter `random_number_generator_seed` is the random seed used to reorder randomly leaves in patience steps and has default value `0`.  
  Parameter `output_file_1`: main output file recording the Gelman-Rubin statistics for the considered trees in TSV format:
  `<index i>TAB<GR value tree i in chain 1>TAB<GR value tree i in chain 2>`.    
  Parameter `output_file_2`: file recording the random leaves orders.

## Class TreeVec

The python class `TreeVec` is implemented in the file  [TreeVec.py](src/cedar/TreeVec.py). 
A `TreeVec` object, that encodes a rooted phylogenetic tree as a vector, can be instantiated from
- a Newick string, expected to be in format 1, or
- an <a href="http://etetoolkit.org/docs/latest/index.html">ete3</a> `Tree` object, or
- a string encoding a `TreeVec` vector itself.

The main methods of the class, aside of the constuctor, are:
- `copy`: creates a copy of a `TreeVec` object;
- `check_vector`: checks that a TreVec object is a valid vector representation of a phylogenetic tree;
- `extract_leaves_order`: extracts from a `TreeVec` object the taxa order used to create the vector representation of the tree;
- `treevec2tree` and `tree2treevec`: conversion between a `TreeVec` object and an <a href="http://etetoolkit.org/docs/latest/index.html">ete3</a> `Tree` object;
- `treevec2str` and `str2treevec`: conversion between a `TreeVec` object and a string encoding the vector;
- `treevec2newick` and `newick2treevec`: conversion between a `TreeVec` object and a Newick string;
- `hop`: performs a HOP tree rearrangement on the tree encoded by the object;
- `hop_neighbourhood_size`: computes the size of the HOP neighbourhood of the tree encoded by the object;
- `hop_neighbourhood`: computes the `TreeVec` representations of the trees in the HOP neighbourhood of the tree encoded by the object;
- `hop_similarity`: computes the HOP similarity or a HOP LCS (see paper for an explanation) with another tree;
- `hop_next`: creates a tree one HOP closer to a target tree;
- `random_hop`: peform a uniform random HOP on the current tree;
- `reorder_leaves`: create a new `TreeVec` object fom the current one differing by a random reordering of the leaves.

