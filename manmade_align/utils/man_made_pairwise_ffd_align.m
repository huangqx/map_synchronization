function [ffd_opt, medianDis] = man_made_pairwise_ffd_align(...
    sourceFFD, sourceSample, targetSample, Para)
% Optimize a free-form deformation that aligns the source shape to the
% target shape

% The FFD coefficient of the source points
C_source_ori = sp_ffd_basis_coeff(sourceFFD, sourceSample);

% Applies non-rigid 
for iter = 1:Para.numIterations_pairwise
    sourcePoints = FFD.ctrlPos_cur*C_source_ori';
    % Compute correspondences
    [Corres, medianDis] = sp_closest_point_corres(sourcePoints,...
        targetSample);
    
    % Deform the shape accordingly
    nC = size(Corres, 2);
    W = sparse(1:nC, 1:nC, Corres(3,:));
    
    Ds = C_source_ori(Corres(1,:),:)';
    Pt = targetSample(:, Corres(2,:));
    
    dimX = size(Ds, 1);
    A = Ds*W*Ds' + Para.lambda_first*eye(dimX) + Para.lambda_smooth*FFD.H_smooth;
    b = Ds*W*Pt' + Para.lambda_first*FFD.ctrlPos_ori';
    sourceFFD.ctrlPos_cur = (A\b)';
end
%
ffd_opt = sourceFFD;

sourceShape_def = sourceShape;
C_source_ori = sp_ffd_basis_coeff(FFD, sourceShape.vertexPoss);
sourceShape_def.vertexPoss = FFD.ctrlPos_cur*C_source_ori';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute closest point pairs between three point clouds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Corres, medianDis] = sp_closest_point_corres(sourcePoints,...
    targetPoints)
%
% Compute cloest point pairs in both directions
[footIds_s_in_t, d_s_in_t] = knnsearch(targetPoints', sourcePoints');
[footIds_t_in_s, d_t_in_s] = knnsearch(sourcePoints', targetPoints');

dis = [d_s_in_t', d_t_in_s'];
Corres_s = [1:length(footIds_s_in_t), footIds_t_in_s'];
Corres_t = [footIds_s_in_t', 1:length(footIds_t_in_s)];

Corres = [Corres_s; Corres_t];
d = sort(dis);
sigma = d(floor(length(d)*0.85));
sigma = max(5e-2, sigma);
weights = sigma*sigma./(sigma*sigma + dis.*dis);

Corres = [Corres; weights];
medianDis = sqrt(mean(dis.*dis));