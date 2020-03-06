clear
clc
addpath('Data')
load(['Test_Engel']);
%the values in NMI and ARI are obtained from in order SSR+pearson, SSSR+spearman, SSR+consine, SR+pearson+spearman, SSR+pearson+consine, SSR+spearman+consine
[NMI,ARI,cluster]=SSRE_combine(in_X,true_labs,10);
