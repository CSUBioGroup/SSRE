# SSRE 
SSRE (Sparse Subspace Representation and similarity Enhancement) is a single-cell clustering framework based on single-cell RNA sequencing (scRNA-seq) data.

## Requirements of Python (version 3.5) packages 
* numpy 
* sklearn 
* pandas 

## Run
To use SSRE on given dataset Engel, please run
```
  cd path_to_SSRE-py/src/
  python main.py -d ../data/Test_Engel.txt -l ../data/Test_Engel_label.txt -p 10
```

- **-d the path of data**\
        data is a gene expression matrix in which each column denotes a gene and each row denotes a cell
- **-l the path of label**\
        label is the cells' true labels for calculating NMI, ARI
- **-p**\
the penalty coefficient
