function [model] = loadfobj(filename)

nv = 0;
nf = 0;
vertex_pos = zeros(3, 200000);
face_vids = zeros(3, 400000);
f_id = fopen(filename, 'r');
while 1
    tline = fgetl(f_id);
    if length(tline) == 0
        continue;
    end
    if tline(1) == -1
        break;
    end
    if tline(1) == 'v' && tline(2) == ' '
        nv = nv + 1;
        vertex_pos(:, nv) = sscanf(tline, 'v %lf %lf %lf\n');
    end
    if tline(1) == 'f'
        nf = nf + 1;
        tp = sscanf(tline, 'f %d//%d %d//%d %d//%d\n');
        face_vids(:, nf) = tp(1:2:6)';
    end
end
fclose(f_id);

model.vertexPoss = vertex_pos(:, 1:nv);
model.faceVIds = face_vids(:, 1:nf);
% model.vertexPoss = model.vertexPoss - mean(model.vertexPoss')'*ones(1,nv);
% box_min = min(model.vertexPoss')';
% box_max = max(model.vertexPoss')';
% scale = max(box_max - box_min);
% model.vertexPoss = model.vertexPoss/scale;

model.vertexPoss = single(model.vertexPoss);
model.faceVIds = uint32(model.faceVIds);

% Remove redundant vertices
nv = size(model.vertexPoss, 2);
flags = zeros(1, nv);
flags(model.faceVIds) = 1;
ids = find(flags);
flags(ids) = 1:length(ids);
model.vertexPoss = model.vertexPoss(:, ids);
model.faceVIds = flags(model.faceVIds);

% center = (max(model.vertexPoss')' + min(model.vertexPoss')')/2;
% model.vertexPoss = model.vertexPoss - center*ones(1, size(model.vertexPoss,2));
% box = max(model.vertexPoss')' - min(model.vertexPoss')';
% s = norm(box);
% model.vertexPoss = model.vertexPoss/s;

% [vertexNors, faceNors] = compute_normals(model.vertexPoss, model.faceVIds);
% model.faceNors = faceNors;
% model.vertexNors = vertexNors;
% model.facePoss = (model.vertexPoss(:, model.faceVIds(1,:))...
%     + model.vertexPoss(:, model.faceVIds(2,:))...
%     + model.vertexPoss(:, model.faceVIds(3,:)))/3;

%scale = 1.0/(max(model.vertexPoss(2,:)) - min(model.vertexPoss(2,:)));
%model.vertexPoss = model.vertexPoss/scale;
%model.vertexPoss(2,:) = model.vertexPoss(2,:) + 0.5;

function [vertexNors, faceNors] = compute_normals(vertexPoss, faceVIds)

num_vertices = size(vertexPoss, 2);
num_faces = size(faceVIds, 2);

vertexNors = zeros(3, num_vertices);
faceNors = zeros(3, num_faces);

for fId = 1:num_faces
    v1Id = faceVIds(1, fId);
    v2Id = faceVIds(2, fId);
    v3Id = faceVIds(3, fId);
    p1 = vertexPoss(:, v1Id);
    p2 = vertexPoss(:, v2Id);
    p3 = vertexPoss(:, v3Id);
    faceNors(:, fId) = face_nor(p1, p2, p3);
    vertexNors(:, v1Id) = vertexNors(:, v1Id) + faceNors(:, fId);
    vertexNors(:, v2Id) = vertexNors(:, v2Id) + faceNors(:, fId);
    vertexNors(:, v3Id) = vertexNors(:, v3Id) + faceNors(:, fId);
end


for vId = 1:num_vertices
    l = norm(vertexNors(:, vId));
    if l > 1e-20
        vertexNors(:, vId) = vertexNors(:, vId)/l;
    end
end

function [n] = face_nor(p1, p2, p3)

e1 = p2 - p1;
e2 = p3 - p1;

n = -[e1(2)*e2(3) - e1(3)*e2(2); e1(3)*e2(1) - e1(1)*e2(3); e1(1)*e2(2) - e1(2)*e2(1)];
if norm(n) > 1e-20
    n = n/norm(n);
end
