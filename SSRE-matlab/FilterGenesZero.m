function [out_X] = FilterGenesZero(in_X)
[m,n]=size(in_X);
flag=(in_X~=0);%m*n 0-1 matrix
flag_count=sum(flag);
out_X=in_X(:,flag_count>0);
end