function [angleDistCurve] = gt_eval_relative(I, Rgt, RR)
% This function evaluates optimized rotations with the ground-truth
% rotations

angleDistCurve = zeros(1, 180);
for edgeId = 1 : size(I,2)
    R1 = Rgt(:,:, I(2, edgeId))*Rgt(:,:,I(1,edgeId))';
    R2 = RR(:,:, edgeId);
    angleErr =...
            acos(max(min((R1(1,:)*R2(1,:)'+R1(2,:)*R2(2,:)'+R1(3,:)*R2(3,:)'-1)/2,1),-1));
    binId = min(180, floor(angleErr*360/2/pi)+1);
    angleDistCurve(binId:180) = angleDistCurve(binId:180) + 1;
end
angleDistCurve = angleDistCurve/size(I,2);
