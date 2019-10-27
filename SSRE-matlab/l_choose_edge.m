function [elected_edge]=l_choose_edge(pear_sim,spear_sim,cos_sim)
[m,n]=size(pear_sim);
if m>5000
    edge_num=100;
else
    edge_num=round(m*0.1);
end
[data1,index1]=sort(pear_sim,2,'descend');
pear_index=index1(:,1:edge_num);
[data2,index2]=sort(spear_sim,2,'descend');
spear_index=index2(:,1:edge_num);
[data1,index3]=sort(cos_sim,2,'descend');
cos_index=index3(:,1:edge_num);
edges={};
for i=1:m
    edges{i}=union(union(pear_index(i,:),spear_index(i,:)),cos_index(i,:));
    elected_edge{i}=sort(edges{i},2);
end



