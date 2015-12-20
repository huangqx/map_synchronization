function [X, eij] = fmap_fitting(Fs, Ft, corresW, Ys, Yt,...
    eigsS, eigsT, lambda)
% Compute the functional map in order to minimize the following objective
% function
% min_{X, D} \|diag(corresW)(X*Fs - DFt)\|_{F}^2 + \|X*Ys - Yt\|_{F}^2 + lambda\|XLs -
% Lt*X\|_{F}^2
% Here D is a diagonal matrix
if 0
    inner = sum(Fs.*Ft);
    D = Fs - Ft.*kron(ones(size(Fs,1),1), inner);
    D = D.*kron(ones(size(D,1),1), corresW);
    H1 = D*D' + Ys*Ys';
    G = Yt*Ys';
else
    Fs = Fs.*kron(ones(size(Fs,1),1), corresW);
    Ft = Ft.*kron(ones(size(Ft,1),1), corresW);
    H1 = Fs*Fs' + Ys*Ys';
    G = Ft*Fs' + Yt*Ys';
end

ms = size(Fs, 1);
mt = size(Ft, 1);

X = zeros(mt, ms);

w = eigsT;
w(1) = eigsT(2);
w = sqrt(1./w);
w = w/w(2);

for i = 1:mt
    d = w(i)*(eigsT(i) - eigsS);
    H = H1 + lambda*diag(d.*d) ;
    g = G(i,:);
    X(i,:) = g*inv(H); 
end

W = sparse(1:length(w),1:length(w),w);

if 0
    tp1 = X*D;
else
    tp1 = W*(X*Fs - Ft);
end
tp2 = W*(X*Ys - Yt);
tp3 = W*(X*diag(eigsS) - diag(eigsT)*X);
eij = sum(sum(tp1.*tp1)) + sum(sum(tp2.*tp2)) + lambda*sum(sum(tp3.*tp3));
