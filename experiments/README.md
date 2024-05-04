# Experiments for the paper "A Vector Representation for Phylogenetic Trees"

## Experiment 1

From 1000 random trees on 100 taxa, comparing the disk space for storing them in Newick format and in vector format, both uncompressed and gzipped.

```
python experiments.py exp1 > exp1_results.txt
```

## Experiment 2

From 50 random trees on 100 taxa, adding HOP trees beween them, and comparing the disk space for storing them in Newick format and in vector format, both uncompressed and gzipped.

```
python experiments.py exp2 > exp2_results.txt
```

## Experiment 3

Computing the size of the HOP neighbourhood for 1000 random trees on 100 taxa.

```
python experiments.py exp3 > exp3_results.txt
```

## Experiment 4

Comparing RF and HOP distances on trees with known SPR distance from a starting tree.
We consider 50 starting trees, each of 100 taxa, and for each we build a sequence of trees by doing random SPRs, from 5 SPRs to 100 SPRs by steps of 5.
We compute the RF distance and the HOP distance, for 10 random orders.

```
python experiments.py exp4
python figures.py exp4_results.txt figures 100 50 10
```
