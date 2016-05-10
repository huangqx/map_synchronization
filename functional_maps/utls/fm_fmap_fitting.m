function [X, eij] = fm_fmap_fitting(Fs, Ft, corresW, Ys, Yt,...
    eigsS, eigsT, lambda)
% Compute the functional map in order to minimize the following objective
% function
% min_{X, D} \|diag(corresW)(X*Fs - DFt)\|_{F}^2 + \|X*Ys - Yt\|_{F}^2 + lambda\|XLs -
% Lt*X\|_{F}^2
% Here D is a diagonal matrix
W = sparse(1:length(corresW), 1:length(corresW), corresW);
H1 = Fs*W*Fs' + Ys*Ys';
G = Ft*W*Fs' + Yt*Ys';

ms = size(Fs, 1);
mt = size(Ft, 1);

X = zeros(mt, ms);

eigsS = eigsS/max(eigsS);
eigsT = eigsT/max(eigsT);

for i = 1:mt
    d = eigsT(i) - eigsS;
    H = H1 + lambda*diag(d.*d) ;
    g = G(i,:);
    X(i,:) = g*inv(H); 
end

W1 = sparse(1:length(corresW), 1:length(corresW), sqrt(corresW));
tp1 = (X*Fs - Ft)*W1;
tp2 = X*Ys - Yt;
tp3 = X*diag(eigsS) - diag(eigsT)*X;
eij(1) = sum(sum(tp1.*tp1));
eij(2) = sum(sum(tp2.*tp2));
eij(3) = lambda*sum(sum(tp3.*tp3));
