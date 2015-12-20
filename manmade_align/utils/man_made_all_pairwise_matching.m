function [unaryTerm, binaryTerm, edges, numOfStates, rootId] =...
    man_made_all_pairwise_matching(...
    SAMPLE,...% The input shapes
    G_knn,... % The shape graph
    Para_align) % The parameters
%
% This function performs all pair-wise affine matching along pairs of shapes
% specified by a shape graph
%
d = full(sum(G_knn));
[s, rootId] = max(d);
%
dimB = 4*Para_align.numRotSamples;
numShapes = length(SAMPLE);
numOfStates = ones(1, numShapes-1)*dimB;
%
IDX = [1:(rootId-1), (rootId+1):numShapes];
%
top1 = G_knn(IDX, rootId);
% Fill in the unary term
unaryTerm = zeros(1, (numShapes-1)*dimB);
ids = find(top1);
fprintf('Computing unary...\n');
for i = 1:length(ids)
    sId = rootId;
    tId = IDX(ids(i));
    pairScores = man_made_pairwise_affine_matching(...
        SAMPLE{sId},...
        SAMPLE{tId},...
        Para_align.numRotSamples);
    subIds = ((ids(i)-1)*dimB+1):(ids(i)*dimB);
    unaryTerm(subIds) = pairScores(:,1)';
end

%
top2 = G_knn(IDX, IDX);
[rows, cols, vals] = find(top2);
ids = find(rows < cols);
sVIds = rows(ids)';
tVIds = cols(ids)';
edges = [sVIds; tVIds];
numEdges = length(sVIds);

binaryTerm = zeros(1, dimB*dimB*numEdges);
fprintf('Computing pairwise...\n');
for eId = 1: numEdges;
    sId = edges(1, eId);
    tId = edges(2, eId);
    pairScores = man_made_pairwise_affine_matching(...
        SAMPLE{sId},...
        SAMPLE{tId},...
        Para_align.numRotSamples);
    Dst = reshape(pairScores, [1, dimB*dimB]);
    ids = ((eId-1)*dimB*dimB+1):(eId*dimB*dimB);
    binaryTerm(ids) = Dst;
    if mod(eId, 100) == 0
        fprintf('  Finished %f.\n', eId*100/numEdges);
    end
end