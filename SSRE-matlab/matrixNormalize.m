
function [Yn,n] = matrixNormalize(Y)

for i = 1:size(Y,2)
    n(i) = norm(Y(:,i));
    Yn(:,i) = Y(:,i) ./ n(i);
end