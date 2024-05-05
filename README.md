# CEDAR: Encoding phylogenetic trees as vectors

CEDAR is a program aimed at manipulating rooted phylogenetic trees encoded as vecors, as described in
the paper *"A Vector Representation for Phylogenetic Trees"*.

**WARNING.** CEDAR is still in development and comes with no warranty.

CEDAR can be used in command line (described below) or wihin python progams using the class `TreeVec` implemented in the file 
[TreeVec.py](src/TreeVec.py).

*Dependencies:* <a href="https://numpy.org/">numpy</a> and <a href="http://etetoolkit.org/docs/latest/index.html">ete3</a>.

The directory [example](example/) contains an example of using CEDAR through command-line.

The directory [experiments](experiments) contains the code to reproduce the experiments described in the 
paper *"A Vector Representation for Phylogenetic Trees"*.

## Command-line usage

The command-line script is [CEDAR.py](src/CEDAR.py), and allows to perform the following tasks:
-  Converting trees written in Newick format into the CEDAR format:
   ```
   python src/CEDAR.py fromNewick --input_file Newick_file --output_file CEDAR_file [--order_file orde]
   ```
   Encodes the files from `Newick_file` in the CEDAR format in file `CEDAR_file`.
   If the file `order` is provided, it is used to determine the taxa order used for the encoding;
   otherwise, the taxa order is defined by a postorder travrsal of the first tree in `Newick_file`.
   The format of a taxa order file is a single line, starting by `#order ` followed by the list of all
   leaves separated by commas.

   **Assumption**: all trees in `Newick_file` are rooted phylogenetic trees on the same set of taxa.

- Creating random taxa orders:
  ```
  python src/CEDAR.py orders --input_file CEDAR_file --output_dir out_dir [--nb_orders N] [--output_prefix out_pref]
  ```
  Creates `N` random taxa orders in files `out_dir/out_pref_I.txt` for `I` from `1` to `N`.
  Default values: `N=1`, `out_pref=CEDAR_random_order`.

  **Assumption**: in all commands, all trees in a `CEDAR_file` are rooted phylogenetic trees on the same set of taxa
  and having been encoded into vector using the same order on taxa, written in the first line of the file.

- Converting trees written in CEDAR format into the Newick format:
  ```
  python src/CEDAR.py toNewick --input_file CEDAR_file --output_file Newick_file [--format Newick_format]
  ```
  The optional `format` agument definesthe Newick format of the created Newick file (default: `1`).
  
- Computing the HOP similarity between trees:
  ```
  python src/CEDAR.py HOP_sim --input_file CEDAR_file --output_file similarity_file [--mode distance_mode]
  ```
  Creates a file that records the HOP similarity between the trees in CEDAR_file.
  - If `mode=pairwise` the HOP similariy is computed between all pairs of trees.
  - If `mode=sequence` the HOP similariy is computed between successive pairs of trees.
  - If `mode=first` the HOP similariy is computed between the first tree and all other trees.
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
  python src/CEDAR.py HOP_path --input_file CEDAR_file --output_file CEDAR_path_file
  ```
  Creates a new CEDAR file where a sequence of trees, each differing from the previous one by a single
  HOP, is inserted between all pairs of successive trees in the input file `CEDAR_file`.

## Class TreeVec

The python class `TreeVec` implemented in the file  [TreeVec.py](src/TreeVec.py). A `TreeVec` object, that
encodes a rooted phylogenetic tree as a vector can be instantiated from
- a Newick string, expcted to be in format 1,
- an <a href="http://etetoolkit.org/docs/latest/index.html">ete3</a> `Tree` object,
- a string encoding a TeeVec vector.

The main mehods of the class, aside of the constuctor, are:
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
- `random_hop`: peform a andom HOP on the current tree.

