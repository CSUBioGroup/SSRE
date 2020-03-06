function [NMI,ARI,cluster]=SSRE_combine(data,label,omiga)
n = length(unique(label));
data = FilterGenesZero(data);%each row is a cell
%calculate pairwise similarity
pairwise_data0=data;
[pear_sim0,spear_sim0,cos_sim0]=ssc_cal_sim(pairwise_data0');%each column of input x should be a cell
%gene score
[pear_score] = LaplacianScore(pairwise_data0, pear_sim0);
[spear_score] = LaplacianScore(pairwise_data0, spear_sim0);
[cos_score] = LaplacianScore(pairwise_data0, cos_sim0);

%calculate sparse similarity
ssc_data0=matrixNormalize(data');%each column is a cell
CMat3 = admmLasso_mat_func(ssc_data0,false,omiga);
C0 = abs(CMat3)+abs(CMat3')+eps;
%gene score
[ssc_score] = LaplacianScore(ssc_data0', C0);

%gene selection with different combination
[gene_select,gene_slect_combine]=l_gene_select_combine(ssc_score,pear_score,spear_score,cos_score);

%test 6 different combination of 4 similarities, ssc+pearson, ssc+spearman,
%ssc+cosine, ssc+pearsaon+spearman, ssc+pearsaon+cosine,
%ssc+spearman+cosine
NMI=zeros(6,1)';
ARI=zeros(6,1)';
for i=1:6
    gene_select1=sort(gene_slect_combine{i});
    data_select=data(:,gene_select1);

    %calculate sparse similarity again
    ssc_data1=matrixNormalize(data_select');%each column is a cell
    CMat3 = admmLasso_mat_func(ssc_data1,false,omiga);
    C1 = abs(CMat3)+abs(CMat3');

    %calculate pairwise similarity again
    pairwise_data=data_select';
    [pear_sim1,spear_sim1,cos_sim1]=ssc_cal_sim(pairwise_data);%each column of input x should be a cell

    [selected_edge_combine]=l_choose_edge_combine(pear_sim1,spear_sim1,cos_sim1);

    elected_edges=selected_edge_combine{i};
    [test_data]=l_enhance(C1,elected_edges);

    cluster = SpectralClustering(test_data,n);
    NMI(i)=Cal_NMI(label, cluster);
    ARI(i)=RandIndex(label, cluster);
end
end