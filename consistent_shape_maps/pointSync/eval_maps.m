function [curve] = eval_maps(maps, feaDisMat)
%
numBins = 256;
curve = zeros(1, 256);

for mapId = 1:length(maps)
    map = maps{mapId};
    [featurePointIds, s] = find(feaDisMat{map.sId}' == 0);
    tVIds = map.corres(2, featurePointIds);
    for j = 1:length(tVIds)
        d = feaDisMat{map.tId}(j, tVIds(j));
        binId = min(numBins, max(1, floor(d*numBins/3.0)));
        curve(binId) = curve(binId) + 1;
    end
end

for i = 2:numBins
    curve(i) = curve(i) + curve(i-1);
end
curve = curve/max(curve);

