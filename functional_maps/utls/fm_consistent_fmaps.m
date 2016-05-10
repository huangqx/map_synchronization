function [fmaps] = fm_consistent_fmaps(Data, Para)
% Input arguments:
%       Data.shapes:       the input shapes      
%       Data.basis:        the precomputed basis
%       Para.lambda_regu:  the regularization prior for pair-wise
%                          functional maps
%       Para.lambda_cons:  the consistency we enforce on the shape
%                          network
%
% Convert point-based maps into functional maps
% Initialize the pair-wise functional maps
fmaps = cell(1, length(Data.input_pair_maps));
for id = 1:length(Data.input_pair_maps)
    sId = Data.input_pair_maps{id}.sId;
    tId = Data.input_pair_maps{id}.tId;
    corres = Data.input_pair_maps{id}.corres;
    Fs = normalize_prob_funcs(Data.basis{sId}.eigVecs(corres(1,:),:)');
    Ft = normalize_prob_funcs(Data.basis{tId}.eigVecs(corres(2,:),:)');
    corresW = corres(3,:)*size(Fs,1)/length(corres(3,:));
    fmaps{id}.sId = sId;
    fmaps{id}.tId = tId;
    fmaps{id}.X = fm_fmap_fitting2(Fs, Ft, corresW, ...
        Data.basis{sId}.eigVals,...
        Data.basis{tId}.eigVals, Para.lambda_regu);
    n = size(corres, 2);
end

n = length(Data.shapes);
m = length(Data.basis{1}.eigVals);

% Perform alternating optimization
for iIter = 1:Para.numIterations_alternate
    % Compute the latent function
    [Y] = fm_latent_func_fitting(fmaps, n, m, Para.dimBasis);
    % Fit the pair-wise funcmaps again
    objVal = zeros(1,3);
    for id = 1:length(Data.input_pair_maps)
        sId = Data.input_pair_maps{id}.sId;
        tId = Data.input_pair_maps{id}.tId;
        sIds = ((sId-1)*m+1):(sId*m);
        tIds = ((tId-1)*m+1):(tId*m);
        corres = Data.input_pair_maps{id}.corres;
        Fs = normalize_prob_funcs(Data.basis{sId}.eigVecs(corres(1,:),:)');
        Ft = normalize_prob_funcs(Data.basis{tId}.eigVecs(corres(2,:),:)');
        corresW = corres(3,:)*size(Fs,1)/length(corres(3,:));
        % Weight the functions
        [fmaps{id}.X, eij] = fm_fmap_fitting(Fs, Ft, corresW,...
            Y(sIds,:)*sqrt(Para.lambda_consistency),...
            Y(tIds,:)*sqrt(Para.lambda_consistency),...
            Data.basis{sId}.eigVals,...
            Data.basis{tId}.eigVals, Para.lambda_regu);
        objVal = objVal + eij;
    end
    fprintf(' iIter = %d, objVal = %.2f, [%.2f, %.2f, %.2f].\n', iIter, sum(objVal), objVal(1), objVal(2), objVal(3));
end

function [S] = normalize_prob_funcs(S_in)
% Normalize so that each function is a normalized vector
S = S_in./kron(ones(size(S_in,1),1), sqrt(sum(S_in.*S_in)));
