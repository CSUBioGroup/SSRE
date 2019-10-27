function [gene_select,gene_inter,gene_inter_num]=l_gene_select(ssc_score,pear_score,spear_score,cos_score)
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


gene_select0=intersect(gene_inter{1},gene_inter{3});
gene_select1=intersect(gene_inter{2},gene_inter{4});
gene_select=intersect(gene_select0,gene_select1);


