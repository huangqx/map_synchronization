function [consistentVertexIds] = csm_sdp_rounding(Data, X, Para)
%
[U,V] = eig(X);
m0 = Para.m0;
m = Para.m;
n = length(Data.shapes);
v = diag(V)';
[s,ids] = sort(v);
ids = ids((length(ids)-m0+1):length(ids));
X_lr = U(:, ids)*V(ids,ids)*U(:,ids)';

ids = (1+(Para.rootId-1)*m):((Para.rootId-1)*m+m0);
Y_lr = X_lr(ids,:);
consistentSampleIds = kron(ones(1,n), (1:m0)');
for i = 1:(Para.rootId-1)
    colIds = (1+(i-1)*m):(i*m);
    Y = Y_lr(:, colIds);
    Y_exp = [1-Y; 10*ones(m-m0,m)];
    sol = lapjv(Y_exp, 0.0);
    sol = sol(1:m0);
    consistentSampleIds(:,i) = sol';
end

for i = (Para.rootId+1):n
    colIds = (1+(i-2)*m+m0):((i-1)*m+m0);
    Y = Y_lr(:, colIds);
    Y_exp = [1-Y; 10*ones(m-m0,m)];
    sol = lapjv(Y_exp, 0.0);
    sol = sol(1:m0);
    consistentSampleIds(:,i) = sol';
end

consistentVertexIds = consistentSampleIds;
for i = 1:n
    consistentVertexIds(:,i) = Data.SAMPLE{i}.sampleIds(consistentSampleIds(:,i)); 
end
