function [test_data]=l_enhance(data,selected_edge)
data=abs(data);
%index=[data~=0];
[m,n]=size(data);
test_data=data;
%complete some links nonexisting in the graph
data=test_data;
RA_score=zeros(m,n);
WRA_score=zeros(m,n);
for i=1:m
    edge_len=length(selected_edge{i});
    for j=1:edge_len
        if data(i,selected_edge{i}(j))==0
            RA=0;
            WRA=0;
            for z=1:m
                if data(i,z)~=0&&data(selected_edge{i}(j),z)~=0
                    neighbor_num=sum(sum(data(z,:)~=0));
                    if neighbor_num==0
                        disp('a error happening!')
                    else
                        RA=RA+1/neighbor_num;
                        WRA=WRA+(data(i,z)+data(selected_edge{i}(j),z))/neighbor_num;
                        %WRA=WRA+min(data(i,z)/neighbor_num,data(selected_edge{i}(j),z)/neighbor_num);
                    end
                end
            end
            RA_score(i,selected_edge{i}(j))=RA;
            WRA_score(i,selected_edge{i}(j))=WRA;
        end
    end
end
%number of exsiting links
nonzeros=sum(sum(data~=0));
links_num=nonzeros/2;

WRA_value=WRA_score+WRA_score';
test_data0=data+WRA_value;
test_data=test_data0-diag(diag(test_data0));
end
