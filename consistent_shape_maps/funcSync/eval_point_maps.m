function [curve] = eval_point_maps(Data, maps, numBins)
%
for id = 1:length(Data.shapes)
    model = Data.shapes{id};
    lc = min(model.vertexPoss')';
    uc = max(model.vertexPoss')';
    diat = norm(uc-lc);
    dist(id) = diat;
end

curve = zeros(1, numBins);
maxDiat = mean(dist);
for id = 1:length(maps)
    sId = maps{id}.sId;
    tId = maps{id}.tId;
    corres = maps{id}.corres;
    sModel = Data.shapes{sId};
    tModel = Data.shapes{tId};
    for j = 1:length(sModel.featurePointIds)
        tId_gt = tModel.featurePointIds(j);
        tId = corres(2, sModel.featurePointIds(j));
        d = tModel.vertexPoss(:, tId_gt) - tModel.vertexPoss(:, tId);
        binId = min(numBins, max(1, floor(sqrt(d'*d)*numBins/maxDiat)));
        curve(binId) = curve(binId) + 1;
    end
end
for i = 2:numBins
    curve(i) = curve(i) + curve(i-1);
end
curve = curve/max(curve);