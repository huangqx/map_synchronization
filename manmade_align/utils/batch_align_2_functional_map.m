function [FMAP] = batch_align_2_functional_map(...
    Shapes,...
    FFDs_opt,...
    Para_align)
% 
% Compute functional basis
for shapeId = 1 : length(Shapes)
    BASIS{shapeId} = manmade_basis(Shapes{shapeId}, Para_align.dimBasis, 1e-2);
    C_source_ori = sp_ffd_basis_coeff(FFDs_opt{shapeId},...
        Shapes{shapeId}.vertexPoss); % FFD coefficient
    Shapes{shapeId}.vertexPoss = FFDs_opt{shapeId}.ctrlPos_cur*C_source_ori';
    % Current point positions

    fprintf('Finished computing basis for shape %d\n', shapeId);
end

for sId = 1 : length(Shapes)
    for tId = 1 : length(Shapes)
        if sId == tId
            continue;
        end
        FMAP{sId, tId} = nearest_neighbor_2_fmap(...
            Shapes{sId}, Shapes{tId}, BASIS{sId}, BASIS{tId});
    end
end

function [fmap] = nearest_neighbor_2_fmap(Shape_s, Shape_t, Basis_s, Basis_t)
%
[IDX_st, d_st] = knnsearch(Shape_s.vertexPoss', Shape_t.vertexPoss');
[IDX_ts, d_ts] = knnsearch(Shape_t.vertexPoss', Shape_s.vertexPoss');

corres = [[IDX_st',1:length(IDX_ts)];[1:length(IDX_st), IDX_ts']];
d_all = [d_st', d_ts'];
sigma = median(d_all);
ids = find(d_all < 3*sigma);
corres = corres(:, ids);
Bs = Basis_s.eigVecs(corres(1,:),:)';
Bt = Basis_t.eigVecs(corres(2,:),:)';

fmap.X = zeros(size(Bt, 1), size(Bs,1));
H1 = Bs*Bs';
G = Bt*Bs';
for i = 1 : size(Bs, 1)
    d = (Basis_t.eigVals(i) - Basis_s.eigVals);
    H = H1 + 1e-1*diag(d.*d) ;
    g = G(i,:);
    fmap.X(i,:) = g*inv(H); 
end