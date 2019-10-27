function [K1] = Estimate_Number_of_Clusters_Gap(W,Num)

    warning off;
    N = size(W,1);

    % Normalized spectral clustering according to Ng & Jordan & Weiss
    % using Normalized Symmetric Laplacian L = I - D^{-1/2} W D^{-1/2}

    DN = diag(1./sqrt(sum(W)+eps) );
    LapN = speye(N) - DN * W * DN;
   
    [~,U,vN] = svd(LapN);
    eigvalue1 = diag(U);
    %eig1_select = eigvalue1((N-Num+1:N));
	eig1_select = eigvalue1((N-Num+1:N));
    eigvalue2 = diag(U);
    %eig2_select = eigvalue2(N-Num:N-1);
	eig2_select = eigvalue2(N-Num:N-1);
  
    tmp = abs(eig2_select-eig1_select);
    [sort_id,index]= sort(tmp,1,'descend');
    
    K1=Num-index(1)+1
end