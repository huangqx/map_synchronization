function [Data] = load_dataset(foldername)
% Load the dataset stored in folder 'foldername'.
% The dataset contains shapes in subfolder 'shapes'
% and pairwise maps in subfolder 'initial_maps'.
% The data is stored in Data.shapes and Data.initial_maps

temp = dir([foldername, '\shapes\']);
id = 1;
for i = 1:(length(temp)-2)
    name = temp(i+2).name;
    sExt = name((length(name)-3):length(name));
    name = name(1:(length(name)-4));
    if strcmp(sExt, '.obj') == 1
        objNames{id} = name;
        id = id + 1;
    end
end

% load shapes
fprintf('Loading shapes...\n');
for id = 1:length(objNames)
    shapeFilePathName = [foldername, '\shapes\', objNames{id}, '.obj'];
    Data.shapes{id} = load_shape(shapeFilePathName);
    featureFilePathName = [foldername, '\shapes\', objNames{id}, '.txt'];
    Data.shapes{id}.name = objNames{id};
    Data.shapes{id}.featurePointIds = load_features(featureFilePathName);
    fprintf('Finished loading shape %d\n', id);
end

% load maps
fprintf('Loading maps...\n');
fprintf('Finished loading maps %d\n', id);
mapId = 0;
for sId = 1:length(objNames)
    for tId = 1:length(objNames)
        if sId == tId
            continue;
        end
        mapId = mapId + 1;
        Data.initial_maps{mapId}.sId = sId;
        Data.initial_maps{mapId}.tId = tId;
        mapFilePathName = [foldername, '\initial_maps\',...
            objNames{sId}, '_', objNames{tId}, '.txt'];
        Data.initial_maps{mapId}.corres = load_map(mapFilePathName);
        fprintf('Finished loading map %d to %d\n', sId, tId);
    end
end

function [Shape] = load_shape(filename)
% Load a 3D model from a file in obj wavefront

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
    if tline(1) == 'v'
        nv = nv + 1;
        vertex_pos(:, nv) = sscanf(tline, 'v %lf %lf %lf\n');
    end
    if tline(1) == 'f'
        nf = nf + 1;
        face_vids(:, nf) = sscanf(tline, 'f %d %d %d\n');
    end
end
fclose(f_id);

Shape.vertexPoss = vertex_pos(:, 1:nv);
Shape.faceVIds = face_vids(:, 1:nf);

function [vIds] = load_features(filename)
vIds = [];

f_id = fopen(filename, 'r');
if f_id == -1
    return;
end
fclose(f_id);

vIds = load(filename);
vIds = vIds(2:length(vIds))'+1;

function [corres] = load_map(filename)
nums = load(filename);
corres = nums' + 1;