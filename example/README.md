# CEDAR: example

The file `random.trees.5.nwk` contains 5 trees on 100 taxa in Newick format.

Converting the trees in CEDAR format:
```
python ../src/CEDAR.py fromNewick --input_file random.trees.5.nwk --output_file random.trees.5.vec
```

Creating a random leaves order file `random_order_1.txt`:
```
python ../src/CEDAR.py orders --input_file random.trees.5.vec --output_dir . --nb_orders=1 --output_prefix random_order
```

Converting the trees in CEDAR format using the created random leaves order:
```
python ../src/CEDAR.py fromNewick --input_file random.trees.5.nwk --output_file random.trees.5.2.vec --order_file random_order_1.txt
```

Converting back `random.trerandom.trees.5.2.veces.5.2.vec` into Newick format:
```
python ../src/CEDAR.py toNewick --input_file random.trees.5.2.vec --output_file random.trees.5.2.nwk
```

Computing the HOP similarity:
```
python ../src/CEDAR.py HOP_sim --input_file random.trees.5.vec --output_file random.trees.5.vec.dist.sequence.txt --mode sequence
python ../src/CEDAR.py HOP_sim --input_file random.trees.5.vec --output_file random.trees.5.vec.dist.pairwise.txt --mode pairwise
python ../src/CEDAR.py HOP_sim --input_file random.trees.5.vec --output_file random.trees.5.vec.dist.first.txt --mode first
```

Computing the HOP neighbourhood size:
```
python ../src/CEDAR.py HOP_ngb1 --input_file random.trees.5.vec --output_file random.trees.5.ngb.txt
```

Computing the HOP neighbourhoods:
```
python ../src/CEDAR.py HOP_ngb2 --input_file random.trees.5.vec --output_dir . --output_prefix random.trees.5.ngb
```

Computing HOP paths between successive trees:
```
python ../src/CEDAR.py HOP_path --input_file random.trees.5.vec --output_file random.trees.5.paths.vec
```