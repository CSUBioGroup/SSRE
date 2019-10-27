# SSRE-master
Cell type detection based on sparse-subspace representation and similarity enhancement
## Description
% > Matlab R2016b <br />
% SSRE(data,label) <br />
% data is the gene expression matrix, each column denotes a gene and each row denotes a cell <br />
% label is the cells'true labels for calculating NMI, ARI <br />
% omiga is the penalty coefficient <br />

[NMI,ARI,cluster] = SSRE(data,label,omiga);

## Example
run demo.m

## Acknowlege

SubKit (https://github.com/sjtrny/SubKit))
