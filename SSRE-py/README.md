SSRE
SSRE (Sparse Subspace Representation and similarity Enhancement) is a single-cell clustering framework based on single-cell RNA sequencing (scRNA-seq) data.

Requirements of Python (version 3.5) packages
numpy
sklearn
pandas
time
matplotlib

Example
To use SSRE on given dataset Engel, please run
cd path_to_SSRE-py/scr/; python main.py  -d ../data/Test_Engel.txt -l ../data/Test_Engel_label.txt -p 10

# Test_Engel.txt is the gene expression matrix, each column denotes a gene and each row denotes a cell
# Test_Engel_label.txt is the cells'true labels for calculating NMI, ARI
# The parameter p is the penalty coefficient
