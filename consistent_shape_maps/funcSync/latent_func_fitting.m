function [Y] = latent_func_fitting(fmaps, n, m, nB)
%
dim = n*m;
X = zeros(dim, dim);
for id = 1:length(fmaps)
    map = fmaps{id};
    sIds = ((map.sId-1)*m+1):(map.sId*m);
    tIds = ((map.tId-1)*m+1):(map.tId*m);
    X(sIds, sIds) = X(sIds, sIds) + map.X'*map.X;
    X(tIds, sIds) = X(tIds, sIds) - map.X;
    X(sIds, tIds) = X(sIds, tIds) - map.X';
    X(tIds, tIds) = X(tIds, tIds) + eye(m);
end

[U,V] = eig(X);
v = diag(V)';
[s,ids] = sort(v);
Y = U(:, ids(1:nB));