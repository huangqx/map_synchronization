function [X] = fm_fmap_fitting2(S, T, weights, sEigs, tEigs, lambda)
% 
ms = size(S, 1);
mt = size(T, 1);
% w = tEigs;
% w(1) = tEigs(2);
% w = sqrt(1./w);
% w = w/w(2);
sEigs = sEigs/max(sEigs);
tEigs = tEigs/max(tEigs);
W = sparse(1:length(weights), 1:length(weights), weights);
H1 = S*W*S';
G = T*W*S';

X = zeros(mt, ms);
for i = 1:mt
    d = (tEigs(i) - sEigs);
    H = H1 + lambda*diag(d.*d) ;
    g = G(i,:);
    X(i,:) = g*inv(H); 
end
