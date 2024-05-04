# CEDAR: Encoding phylogenetic trees as vectors

CEDAR is a program aimed at manipulating rooted phylogenetic trees encoded as vecors, as described in
the paper "A Vector Representation for Phylogenetic Trees".

**WARNING.** CEDAR is still in development.

CEDAR allows to perform the following tasks:
-  Converting trees written in Newick format into the CEDAR format:
   ```
   python ../src/CEDAR.py fromNewick --input_file Newick_file --output_file CEDAR_file [--order_file orde]
   ```
   Encodes the files from `Newick_file` in the CEDAR format in file `CEDAR_file`.
   If the file `order` is provided, it is used to determine the taxa order used for the encoding;
   otherwise, the taxa order is defined by a postorder travrsal of the first tree in `Newick_file`.
   The format of a taxa order file is a single line, starting by `#order ` followed by the list of all
   leaves separated by commas.

   **Assumption**: all trees in `Newick_file` are rooted phylogenetic trees on the same set of taxa.

- Creating random taxa orders:
  ```
  python ../src/CEDAR.py orders --input_file CEDAR_file --output_dir out_dir [--nb_orders N] [--output_prefix out_pref]
  ```
  Creates `N` random taxa orders in files `out_dir/out_pref_I.txt` for `I` from `1` to `N`.
  Default values: `N=1`, `out_pref=CEDAR_random_order`.

  **Assumption**: all trees in `CEDAR_file` are rooted phylogenetic trees on the same set of taxa.

- Converting trees written in CEDAR format into the Newick format:
  ```
  python ../src/CEDAR.py toNewick --input_file CEDAR_file --output_file Newick_file [--format Newick_format]
  ```
  The optional `format` agument definesthe Newick format of the created Newick file (default: `1`).
  
  **Assumption**: all trees in `CEDAR_file` are rooted phylogenetic trees on the same set of taxa.

- Computing the HOP similarity between trees:
  ```
  python ../src/CEDAR.py HOP_sim --input_file CEDAR_file --output_file similarity_file [--mode distance_mode]
  ```
  Creates a file that records the HOP similarity between the trees in CEDAR_file.
  - If `mode=pairwise` the HOP similariy is computed between all pairs of trees.
  - If `mode=sequence` the HOP similariy is computed between successive pairs of trees.
  - If `mode=first` the HOP similariy is computed between the first tree and all other trees.
  - The default mode is `sequence`.

  **Assumption**: all trees in `CEDAR_file` are rooted phylogenetic trees on the same set of taxa.

-  Computing the HOP neighbourhood size:
   ```
   python ../src/CEDAR.py HOP_ngb1 --input_file CEDAR_file --output_file size_file
   ```
   Computes the size of the HOP neighbourhood for every tree in `CEDAR_file`

  **Assumption**: all trees in `CEDAR_file` are rooted phylogenetic trees on the same set of taxa.

- Computing the HOP neighbourhoods:
  ```
  python ../src/CEDAR.py HOP_ngb2 --input_file CEDAR_file --output_dir out_dir [--output_prefix out_pref]
  ```
  Computes the trees in the HOP neighbourhood for every tree in `CEDAR_file`. The HOP neighbourhood for
  tree `i` is in the file `out_dir/out_pref.i.vec`. The default value for `out_pref` is `CEDAR_HOP_neighbourhood`.

  **Assumption**: all trees in `CEDAR_file` are rooted phylogenetic trees on the same set of taxa.

- Computing HOP paths between successive trees:
  ```
  python ../src/CEDAR.py HOP_path --input_file CEDAR_file --output_file CEDAR_path_file
  ```
  Creates a new CEDAR file where a sequence of trees, each differing from the previous one by a single
  HOP, is inserted between all pairs of successive trees in the input file `CEDAR_file`.

  **Assumption**: all trees in `CEDAR_file` are rooted phylogenetic trees on the same set of taxa.

