function [R_out, lambda] = lap_rot_sync(I, RR, edgeWeights)

numObjects = max(max(I));
numEdges = size(I, 2);
dim = 3*numObjects;


% Generate the lower part of the data matrix
% Row indices of the matrix entries
rowIdsV = kron(3*(I(2,:)-1), ones(3,1)) + kron(ones(1, numEdges), (1:3)');
rowIdsV = kron(rowIdsV, ones(1,3));
rowIdsV = reshape(rowIdsV, [9, numEdges]);

% Column indices of the matrix entries
colIdsV = kron(3*(I(1,:)-1), ones(1,3)) + kron(ones(1, numEdges), (1:3));
colIdsV = kron(colIdsV, ones(3,1));
colIdsV = reshape(colIdsV, [9, numEdges]);

% Values of the matrix entries
RR = reshape(RR, [9, numEdges]);
if length(edgeWeights) > 0
    RR = RR.*kron(ones(9,1), edgeWeights);
end
% Generate the data matrix
R = sparse(rowIdsV, colIdsV, RR, dim, dim);

% From the data matrix to the connection Laplacian
% Generate the adjacency matrix
Adj = sparse(I(1,:), I(2,:), edgeWeights, numObjects, numObjects);
Adj = Adj + Adj';
% Compute the vertex degrees
d = full(sum(Adj));

% Generate the connection Laplacian
Laplacian_connection = sparse(1:dim, 1:dim, kron(d, ones(1,3))) - R - R';

% Perform the eigen-decomposition, and compute the top eigenvectors
[B, Lambda] = eigs(Laplacian_connection, 3, -1e-7);
R_out = zeros(3, 3, numObjects);
lambda = max(diag(Lambda));

for i = 1:numObjects
    [U,V,W] = svd(B((3*i-2):(3*i),:));
    if det(U*W') > 0
        R_out(:,:,i) = U*W';
    else
        R_out(:,:,i) = -U*W';
    end
end

