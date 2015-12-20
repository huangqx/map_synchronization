function [corres] = func_2_point(sBasis, tBasis, X)
%
nP = size(sBasis.vecs,1);
corres = [1:nP;1:nP];

for i = 1:nP
    t = X*sBasis.vecs(i,:)';
    tfunc = tBasis.vecs*t;
    [s,id] = max(tfunc);
    corres(2, i) = id;
end