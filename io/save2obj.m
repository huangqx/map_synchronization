function [] = save2obj(model, filename)

box_max = max(model.vertexPoss')';
box_min = min(model.vertexPoss')';
scale = 1/(box_max(2) - box_min(2));
box_center = (box_min+box_max)/2;
box_center(2) = box_min(2);
%model.vertexPoss = model.vertexPoss - box_center*single(ones(1, size(model.vertexPoss,2)));
%model.vertexPoss = scale*model.vertexPoss;

f_id = fopen(filename, 'w');

for vId = 1:size(model.vertexPoss, 2)
    pos = model.vertexPoss(:, vId);
    fprintf(f_id, 'v %f %f %f\n', pos(1), pos(2), pos(3)); 
end

for fId = 1:size(model.faceVIds, 2)
    ids = model.faceVIds(:, fId);
    fprintf(f_id, 'f %d %d %d\n', ids(1), ids(2), ids(3));
end

fclose(f_id);

