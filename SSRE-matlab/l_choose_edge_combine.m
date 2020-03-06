function [selected_edge_combine]=l_choose_edge_combine(pear_sim,spear_sim,cos_sim)
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
%selected edge with different combination of similarities
selected_edge_combine{1}=num2cell(pear_index);
selected_edge_combine{2}=num2cell(spear_index);
selected_edge_combine{3}=num2cell(cos_index);

edges={};
elected_edge={};
for i=1:m
    edges{i}=union(pear_index(i,:),spear_index(i,:));
    elected_edge{i}=sort(edges{i},2);
end
selected_edge_combine{4}=elected_edge;

edges={};
elected_edge={};
for i=1:m
    edges{i}=union(pear_index(i,:),cos_index(i,:));
    elected_edge{i}=sort(edges{i},2);
end
selected_edge_combine{5}=elected_edge;

edges={};
elected_edge={};
for i=1:m
    edges{i}=union(spear_index(i,:),cos_index(i,:));
    elected_edge{i}=sort(edges{i},2);
end
selected_edge_combine{6}=elected_edge;
