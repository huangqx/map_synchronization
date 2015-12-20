function [pMaps] = batch_func_2_point(Data, fMaps)
%
for id = 1:length(fMaps)
    sId = fMaps{id}.sId;
    tId = fMaps{id}.tId;
    pMaps{id}.sId = sId;
    pMaps{id}.tId = tId;
    pMaps{id}.corres =...
        func_2_point(Data.basis{sId}, Data.basis{tId}, fMaps{id}.X);
    fprintf('Finished %d/%d\n', id, length(fMaps));
end