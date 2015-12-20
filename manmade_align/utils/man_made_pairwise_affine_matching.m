function [pairScores] = man_made_pairwise_affine_matching(...
        sourceShape,...
        targetShape,...
        numRotSamples)
% This function calculates the distance between rotated copies of two
% shapes
dimX = numRotSamples;
dimY = numRotSamples*4;
%
aff_trans = zeros(dimX, dimY);
%
for i = 1 : dimX
    thetaX = (i-1)*pi/dimX/2;
    Rx = [cos(thetaX), 0, -sin(thetaX);
        0, 1, 0;
        sin(thetaX), 0, cos(thetaX)];
    sP = Rx*sourceShape;
    for j = 1 : dimY
        thetaY = (j-1)*pi*2/dimY;
        Ry = [cos(thetaY), 0, -sin(thetaY);
              0, 1, 0;
              sin(thetaY), 0, cos(thetaY)];
        tP = Ry*targetShape;
        aff_trans(i,j) = calculate_error(sP, tP);
    end
end

pairScores = zeros(dimY, dimY);
pairScores(:, 1:dimX) = aff_trans';
ids = [(dimX+1):(4*dimX), 1:dimX];
pairScores(:, (dimX+1):(2*dimX)) = aff_trans(:, ids)';
ids = [(2*dimX+1):(4*dimX), 1:(2*dimX)];
pairScores(:, (2*dimX+1):(3*dimX)) = aff_trans(:, ids)';
ids = [(3*dimX+1):(4*dimX), 1:(3*dimX)];
pairScores(:, (3*dimX+1):(4*dimX)) = aff_trans(:, ids)';


function [error] = calculate_error(sP, tP)
%
sBox = max(sP')' - min(sP')';
tBox = max(tP')' - min(tP')';

meanBox = (sBox + tBox)/2;
sScale = meanBox./sBox;
tScale = meanBox./tBox;
%
sP(1,:) = sP(1,:)*sScale(1);
sP(2,:) = sP(2,:)*sScale(2);
sP(3,:) = sP(3,:)*sScale(3);
%
tP(1,:) = tP(1,:)*tScale(1);
tP(2,:) = tP(2,:)*tScale(2);
tP(3,:) = tP(3,:)*tScale(3);
%
[IDX, dist_st] = knnsearch(tP', sP');
[IDX, dist_ts] = knnsearch(sP', tP');
error = sqrt((sum(dist_st.*dist_st) + sum(dist_ts.*dist_ts))/2/size(sP, 2));