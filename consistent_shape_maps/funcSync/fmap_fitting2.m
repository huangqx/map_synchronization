function [X] = fmap_fitting2(S, T, sEigs, tEigs, lambda)
% 
ms = size(S, 1);
mt = size(T, 1);
w = tEigs;
w(1) = tEigs(2);
w = sqrt(1./w);
w = w/w(2);

H1 = S*S';
G = T*S';

X = zeros(mt, ms);
for i = 1:mt
    d = w(i)*(tEigs(i) - sEigs);
    H = H1 + lambda*diag(d.*d) ;
    g = G(i,:);
    X(i,:) = g*inv(H); 
end
