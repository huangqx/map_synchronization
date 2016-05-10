function [Image] = fm_visualize_basis(Shape, Basis, Camera)
%
for id = 1:size(Basis.eigVecs, 2)
    images{id} = fm_render_func(Shape, Basis.eigVecs(:, id), Camera);
end
%
[dimX, dimY, k] = size(images{1});
dimX2 = dimX*5;
dimY2 = dimY*6;
Image = 255 - uint8(zeros(dimX2, dimY2,3));
for rowId = 1 : 5
    rowIds = ((rowId-1)*dimX + 1):(rowId*dimX);
    for colId = 1 : 6
        colIds = ((colId-1)*dimY + 1):(colId*dimY);
        basisId = (rowId-1)*6 + colId;
        Image(rowIds, colIds, :) = images{basisId};
    end
end