# SSRE-matlab
Cell type detection based on sparse-subspace representation and similarity enhancement
## Description
% > Matlab R2016b <br />
% SSRE(data,label) <br />
% data is the gene expression matrix, each column denotes a gene and each row denotes a cell <br />
% label is the cells'true labels for calculating NMI, ARI <br />
% omiga is the penalty coefficient <br />
% you can choose different combination between SSR and three correlation methods (pearson, spearman, consine) by running demo_combine.m <br 

[NMI,ARI,cluster] = SSRE(data,label,omiga);

## Example
run demo.m

## Acknowlege

SubKit (https://github.com/sjtrny/SubKit))
