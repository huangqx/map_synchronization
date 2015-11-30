function [Shape_out] = clean_mesh(Shape_in)
%
[subIds, idt] = extract_sub_indices(Shape_in.vertexPoss, 64);
Shape_out = Shape_in;
Shape_out.vertexPoss = Shape_in.vertexPoss(:, idt);
for i = 1:length(Shape_out.meshes)
    Shape_out.meshes{i}.vertexIds = subIds(Shape_out.meshes{i}.vertexIds);
end
Shape_out.faceVIds = subIds(Shape_out.faceVIds);

tp = sort(Shape_out.faceVIds);
[subIds2, idt2] = extract_sub_indices(double(tp), 16);
Shape_out.faceVIds = Shape_out.faceVIds(:, idt2);


function [subIds, idt] = extract_sub_indices(vertexPoss, k)
%
[IDX,D] = knnsearch(vertexPoss', vertexPoss', 'K', k);
subIds = zeros(1, size(vertexPoss,2));
idt = zeros(1, size(vertexPoss, 2));
off = 0;
for i = 1:size(IDX, 1)
    if subIds(i) > 0
        continue;
    end
    tp = find(D(i,:) < 1e-15);
    tp = IDX(i, tp);
    off = off + 1;    
    subIds(i) = off;
    idt(off) = i;
    if length(tp) > 0
        subIds(tp) = off;
    end
end
idt = idt(1:off);