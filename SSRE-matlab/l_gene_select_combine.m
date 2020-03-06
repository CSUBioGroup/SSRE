function [gene_select,gene_slect_combine,gene_inter_num]=l_gene_select_combine(ssc_score,pear_score,spear_score,cos_score)
score_set={ssc_score,pear_score,spear_score,cos_score};%gene_num*1
for i=1:4
    [score,index]=sort(score_set{i},'descend');
    gene_num=length(score_set{i});
    thred1=round(0.1*gene_num);
    thred2=round(0.5*gene_num);
    for  j=thred1:1:thred2
        score1=score(1:j);
        score2=score(j+1:gene_num);
        var1=var(score1);
        var2=var(score2);
        gene_var(j)=var1+var2;
    end
    gene_var(1:thred1-1)=inf;
    [min_gene_var,select_index]=min(gene_var);
    
    gene_inter{i}=index(1:select_index);
    gene_inter_num{i}=length(gene_inter{i});
    
end

gene_slect_ssc_pear=intersect(gene_inter{1},gene_inter{2});
gene_slect_ssc_spear=intersect(gene_inter{1},gene_inter{3});
gene_slect_ssc_cos=intersect(gene_inter{1},gene_inter{4});
%combination of ssc and each pairwise similarity
gene_slect_combine{1}=gene_slect_ssc_pear;
gene_slect_combine{2}=gene_slect_ssc_spear;
gene_slect_combine{3}=gene_slect_ssc_cos;
%combination of ssc and any two pairwise similarities
gene_slect_combine{4}=intersect(gene_slect_ssc_pear,gene_inter{3});
gene_slect_combine{5}=intersect(gene_slect_ssc_pear,gene_inter{4});
gene_slect_combine{6}=intersect(gene_slect_ssc_spear,gene_inter{4});
%combination of ssc and three pairwise similarities
gene_select=intersect(gene_slect_combine{4},gene_inter{4});


