clear
clc
addpath('Data')
load(['Test_Engel']);
[NMI,ARI,cluster]=SSRE(in_X,true_labs,10);
