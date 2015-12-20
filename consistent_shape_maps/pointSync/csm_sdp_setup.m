function [SDP] = csm_sdp_setup(Data, Para)
% Generate the sdp optimization problem given the input data
% Input parameters:
%  Para.numMatchedPoints: the number of consistent points across the
%                         dataset. The default number is 32
%  Para.numSamples: the number of samples per shape used for the
%                         computation. The default number is 128. This is
%                         used to ensure that we can find a consistent
%                         subset. Increase this number for shapes that have
%                         large distortion
% Output :
%  SDP.C: objective function
%  SDP.A, SDP.b: the linear constraint: SDP.A(X) = SDP.b
%  SDP.IdsZ: indices of elements that are non negative.

m = Para.m;
m0 = Para.m0;
n = length(Data.shapes);
dim = 1+m*n;

% Compute the maximum diameters
for i = 1:length(Data.shapes)
    diat(i) = max(max(Data.SAMPLE{i}.distMat));
end
mean_diat = mean(diat);

C = zeros(dim, dim);
for mapId = 1:length(Data.initial_maps)
    map = Data.initial_maps{mapId};
%    [C_st] = pairwise_data_matrix(distMats{map.tId},...
%        map.corres, m0, m);
    [C_st] = pairwise_data_matrix(Data.SAMPLE{map.sId},...
        Data.SAMPLE{map.tId},...
        map.corres, m0, m, mean_diat);

    rowIds = ((map.sId-1)*m + 2):(map.sId*m+1);
    colIds = ((map.tId-1)*m + 2):(map.tId*m+1);
    C(rowIds, colIds) = C_st;
end

rootId = Para.rootId;
IDX = [1:(1+(Para.rootId-1)*m+m0), (2+Para.rootId*m):(1+m*n)];
SDP.C = C(IDX, IDX);
SDP.C = (SDP.C+SDP.C')/2;

[SDP.A, SDP.b, SDP.ids_geq_0, SDP.ids_leq_1] = gen_constraints(...
    n, m0, m, rootId);
% 
% function [sampleMaps, distMats] = vertexMap_2_sampleMaps(Data, m)
% % Convert maps between vertices into maps between sample points (m per
% % shape)
% 
% for shapeId = 1:length(Data.SAMPLE)
%     Sample = Data.SAMPLE{shapeId};
%     % Extract distance matrix and closestSampleIds
%     distMats{shapeId} = Sample.distMat(1:m, Sample.sampleIds(1:m));
%     [s, closestSampleIds{shapeId}] = min(Sample.distMat(1:m,:));
% end

% for mapId = 1:length(Data.initial_maps)
%     map = Data.initial_maps{mapId};
%     sampleIds_s = 1:m;
%     vertexIds_s = Data.SAMPLE{map.sId}.sampleIds(sampleIds_s);
%     vertexIds_t = map.corres(2, vertexIds_s);
%     sampleIds_t = closestSampleIds{map.tId}(vertexIds_t);
%     map.corres = [sampleIds_s;sampleIds_t];
%     sampleMaps{mapId} = map;
% end

function [C_st] = pairwise_data_matrix(Sample_s, Sample_t, corres_st,...
    m0, m, mean_diat)
%
%C_st = min(1, distMat_t(corres_st(2,:),:));
sVIds = Sample_s.sampleIds(1:m);
tVIds = corres_st(2, sVIds);
C_st = min(mean_diat/3, Sample_t.distMat(1:m, tVIds)');

% Favor points that are sampled early
t = (floor(m/m0)+1)^2;
v = 1:((t-1)/(m-1)):t;
v = sqrt(v);
d1 = power(v, 0.5);
d2 = power(v, 0.5);
C_st = diag(d1)*C_st*diag(d2);

% sVIds = Sample_s.sampleIds;
% tVIds = corres_st(2, sVIds);
% tSIds = Sample_t.closestSampleIds(tVIds);
% sDis = Sample_s.distMat(:, sVIds);
% tDis = Sample_t.distMat(tSIds, Sample_t.sampleIds(tSIds));
% scale = mean(mean(tDis))/mean(mean(sDis));
% sDis = sDis*scale;
% sDis = sDis + 0.1;
% tDis = tDis + 0.1;
% dif = abs(sDis - tDis);
% ave = (sDis + tDis)/2;
% ratio = dif./ave;
% ratio = sqrt(mean(mean(ratio.*ratio)));
% w = exp(-ratio^2/0.08);
% %C_st = C_st*w;

function [A,b, ids_geq_0, ids_leq_1] = gen_constraints(n, m0, m, rootId)
% Generate the constraints
dim = 1 + (n-1)*m + m0;
F = ones(dim, dim);
for i = 1:(rootId-1)
    ids = (2+(i-1)*m):(1+i*m);
    F(ids,ids) = 2;
end
ids = (2+(rootId-1)*m):(1+(rootId-1)*m+m0);
F(ids,ids) = 2;
for i = (rootId+1):n
    ids  = (2+(i-2)*m+m0):(1+(i-1)*m+m0);
    F(ids,ids) = 2;
end
for i = 1:dim
    F(i, 1:i) = 0;
end
F(1,:) = 0;
F = sparse(F);

% ids_geq_zero
[rows, cols, vals] = find(F == 1);
ids_geq_0 = (cols-1)*dim + rows;

% ids_leg_one
ids_x = [2:(1+(rootId-1)*m),(2+(rootId-1)*m+m0):(1+(n-1)*m+m0)];
ids_leq_1 = (ids_x'-1)*dim + ids_x';

% Ax ==b
%X(1,1) == m0
rowsA = [1];
colsA = [1];
valsA = [1];
b = [m0];

% diag(X_{rr}) == 1
ids_1 = (2+(rootId-1)*m):(1+(rootId-1)*m+m0);
rowsA = [rowsA, 2:(m0+1)];
colsA = [colsA, (ids_1-1)*dim + ids_1];
valsA = [valsA, ones(1, m0)];
b = [b, ones(1,m0)];

% x_{r} == 1
rowsA = [rowsA, (m0+2):(2*m0+1)];
colsA = [colsA, (ids_1-1)*dim + 1];
valsA = [valsA, ones(1, m0)];
b = [b, ones(1,m0)];

% diag(X_ii) = x_i, 1\leq i \leq n
tp = (2*m0+2):(2*m0+1+length(ids_x));
rowsA = [rowsA, kron(tp, ones(1,2))];
tp = [(ids_x-1)*dim+1;(ids_x-1)*dim+ids_x];
colsA = [colsA, reshape(tp, [1, size(tp,1)*size(tp,2)])];
valsA = [valsA, kron(ones(1,size(tp,2)), [1,-1])];
b = [b, zeros(1, size(tp,2))];
rowId = max(rowsA);
for i = 1:(rootId-1)
    ids = (2+(i-1)*m):(1+i*m);
    rowId = rowId + 1;
    rowsA = [rowsA, rowId*ones(1,m)];
    colsA = [colsA, (ids-1)*dim+ids];
    valsA = [valsA, ones(1,m)];
    b = [b, m0];
end
for i = (rootId+1):n
    ids = (2+(i-2)*m+m0):(1+(i-1)*m+m0);
    rowId = rowId + 1;
    rowsA = [rowsA, rowId*ones(1,m)];
    colsA = [colsA, (ids-1)*dim+ids];
    valsA = [valsA, ones(1,m)];
    b = [b, m0];
end

[rows, cols, vals] = find(F==2);
rowsA = [rowsA, (max(rowsA)+1):(max(rowsA) + length(rows))];
colsA = [colsA, (cols'-1)*dim+rows'];
valsA = [valsA, ones(1, length(rows))];
b = [b, zeros(1, length(rows))];


b = b';
A = sparse(rowsA, colsA, valsA);