function [vertexMaps] = batch_feature_2_vertex(...
    Data,...
    consVertexIds,...
    lambda)

numShapes = length(Data.shapes);

mapId = 0;
for sId = 1:numShapes
    for tId = 1:numShapes
        if sId == tId
            continue;
        end
        sVIds = consVertexIds(:, sId);
        tVIds = consVertexIds(:, tId);
        S = Data.basis{sId}.vecs(sVIds,:)';
        S = S./kron(ones(size(S,1),1), sqrt(sum(S.*S)));
        T = Data.basis{tId}.vecs(tVIds,:)';
        T = T./kron(ones(size(T,1),1), sqrt(sum(T.*T)));
        vMap.sId = sId;
        vMap.tId = tId;
        X = fmap_fitting2(S,T,...
            Data.basis{sId}.vals,...
            Data.basis{tId}.vals,...
            lambda);
        vMap.corres =...
            func_2_point(Data.basis{sId}, Data.basis{tId}, X);
        
        mapId = mapId + 1;
        vertexMaps{mapId} = vMap;
        fprintf('Finished converting %d to %d\n', sId, tId); 
    end
end
