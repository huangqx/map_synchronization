function [fmaps] = joint_fmap_opt_lb(Data, Para)
% Initialize the pair-wise functional maps
for id = 1:length(Data.initial_maps)
    sId = Data.initial_maps{id}.sId;
    tId = Data.initial_maps{id}.tId;
    corres = Data.initial_maps{id}.corres;
    S = Data.basis{sId}.vecs(corres(1,:),:)';
    S = S./kron(ones(size(S,1),1), sqrt(sum(S.*S)));
    T = Data.basis{tId}.vecs(corres(2,:),:)';
    T = T./kron(ones(size(T,1),1), sqrt(sum(T.*T)));
    fmaps{id}.sId = sId;
    fmaps{id}.tId = tId;
    fmaps{id}.X = fmap_fitting2(S,T,...
        Data.basis{sId}.vals,...
        Data.basis{tId}.vals, Para.lambda);
    n = size(corres, 2);
    Data.initial_maps{id}.corres = [Data.initial_maps{id}.corres; ones(1,n)];
end

n = length(Data.shapes);
m = length(Data.basis{1}.vals);

% Perform alternating optimization
for oIter = 1:Para.nIters_outer
    fprintf('oIter = %d.....\n', oIter);
    for iIter = 1:Para.nIters_inner
        % Compute the latent function
        [Y] = latent_func_fitting(fmaps, n, m, Para.nB);
        % Fit the pair-wise funcmaps again
        objVal = 0;
        for id = 1:length(Data.initial_maps)
            sId = Data.initial_maps{id}.sId;
            tId = Data.initial_maps{id}.tId;
            sIds = ((sId-1)*m+1):(sId*m);
            tIds = ((tId-1)*m+1):(tId*m);
            corres = Data.initial_maps{id}.corres;
            Fs = normalize_prob_funcs(Data.basis{sId}.vecs(corres(1,:),:)');
            Ft = normalize_prob_funcs(Data.basis{tId}.vecs(corres(2,:),:)');
            corresW = sqrt(corres(3,:));
            corresW = corresW*size(Fs,1)/length(corresW);
            % Weight the functions
            [fmaps{id}.X, eij] = fmap_fitting(Fs, Ft, corresW,...
                Y(sIds,:)*sqrt(Para.mu),...
                Y(tIds,:)*sqrt(Para.mu),...
                Data.basis{sId}.vals,...
                Data.basis{tId}.vals, Para.lambda);
            objVal = objVal + eij;
        end
        fprintf(' iIter = %d, objVal = %.2f.\n', iIter, objVal);
    end
    % Reweight the correspondences
    numCorres = 0;
    for id = 1:length(Data.initial_maps)
        numCorres = numCorres + size(Data.initial_maps{id}.corres,2);
    end
    dist2 = zeros(1, numCorres);
    numCorres = 0;
    for id = 1:length(Data.initial_maps)
        sId = Data.initial_maps{id}.sId;
        tId = Data.initial_maps{id}.tId;
        corres = Data.initial_maps{id}.corres;
        Fs = normalize_prob_funcs(Data.basis{sId}.vecs(corres(1,:),:)');
        Ft = normalize_prob_funcs(Data.basis{tId}.vecs(corres(2,:),:)');
        sqrRes2 = squared_corres_residual(fmaps{id}.X, Fs, Ft, Data.basis{tId}.vals);
        Data.initial_maps{id}.corres(3,:) = sqrRes2;
        dist2((numCorres+1):(numCorres+length(sqrRes2))) = sqrRes2;
        numCorres = numCorres + length(sqrRes2);
    end
    sigma2 = median(dist2);
    for id = 1:length(Data.initial_maps)
        Data.initial_maps{id}.corres(3,:) =...
            sigma2./(sigma2 + Data.initial_maps{id}.corres(3,:));
    end
    fprintf('sigma2 = %f.\n', sigma2);
end

function [S] = normalize_prob_funcs(S_in)
% Normalize so that each function is a normalized vector
S = S_in./kron(ones(size(S_in,1),1), sqrt(sum(S_in.*S_in)));

function [sqrDist] = squared_corres_residual(X, S, T, tEigs)

w = tEigs;
w(1) = tEigs(2);
w = sqrt(1./w);
w = w/w(2);
W = sparse(1:length(w), 1:length(w), w);

T1 = W*X*S;
T = W*T;

if 0
    inner = sum(T1.*T)./sum(T.*T);
    DMat = T1 - kron(ones(size(T1,1),1), inner).*T;
else
    DMat = T1 - T;
end
sqrDist = sum(DMat.*DMat);

