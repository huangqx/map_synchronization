function [] = consistent_basis(Data, Para, Camera)
% Input arguments:
%       Data.shapes:   input shapes
%       Data.basis:    the basis per shape
%       Data.fmaps:    the pre-computed pair-wise fmaps
% Output argument:
%       Segs:          the face segment ids
%
n = length(Data.shapes);
m = length(Data.basis{1}.eigVals);

dim = n*m;
X = zeros(dim, dim);
for id = 1:length(Data.consistent_fmaps)
    map = Data.consistent_fmaps{id};
    sIds = ((map.sId-1)*m+1):(map.sId*m);
    tIds = ((map.tId-1)*m+1):(map.tId*m);
    X(sIds, sIds) = X(sIds, sIds) + map.X'*map.X;
    X(tIds, sIds) = X(tIds, sIds) - map.X;
    X(sIds, tIds) = X(sIds, tIds) - map.X';
    X(tIds, tIds) = X(tIds, tIds) + eye(m);
end
%
X = X*Para.lambda_consistency;
for i = 1:n
    ids = ((i-1)*m+1):(i*m);
    eigVals = Data.basis{i}.eigVals/max(Data.basis{i}.eigVals);
    X(ids, ids) = X(ids, ids) + diag(eigVals);
end
%
[U,V] = eigs(X, 8, 1e-10);
for basisId = 1:8
    figure(basisId);
    for id = 1 : n
        ids = ((id-1)*m+1):(id*m);
        coord = U(ids, 9-basisId);
        func = Data.basis{id}.eigVecs*coord;
        images{id} = fm_render_func(Data.shapes_ori{id}, func', Para.Camera);
    end
    for shapeId = 1: length(Data.shapes)
        subplot(4,5, shapeId);
        imshow(images{shapeId});
    end
end
