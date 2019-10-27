function [pear_sim,spear_sim,cos_sim,eucli_sim]=ssc_cal_sim(x)
dist_cos=pdist(x','cosine');
dist_cos=squareform(dist_cos);
cos_sim0=1-dist_cos;
dist_eucli=pdist(x','euclidean');
dist_eucli=squareform(dist_eucli);
eucli_sim0=1-dist_eucli;
dist_pearson=corr(x,'type','pearson');
dist_spearman=corr(x,'type','spearman');
if min(min(dist_pearson))<0
    display('minimize of pearson less than 0');
    dist_pearson(dist_pearson<0)=0;
    pear_sim=dist_pearson;
end
if min(min(dist_spearman))<0
    display('minimize of spearman less than 0');
    dist_spearman(dist_spearman<0)=0;
    spear_sim=dist_spearman;
end
if max(max(cos_sim0))>1
    display('minimize of cosine less than 0');
    cos_sim0(cos_sim0>1)=0;
    cos_sim=cos_sim0;
end
if max(max(eucli_sim0))>1
    display('minimize of cosine less than 0');
    eucli_sim0(eucli_sim0>1)=0;
    eucli_sim=eucli_sim0;
end
pear_sim=dist_pearson;
spear_sim=dist_spearman;
cos_sim0(cos_sim0<0)=0;
eucli_sim0(eucli_sim0<0)=0;
cos_sim=cos_sim0;
eucli_sim=eucli_sim0;

end