function [] = save_corres(Data, consistentSampleIds, foldername)

for i = 1:length(Data.shapes)
    name = Data.shapes{i}.name;
    vIds = Data.SAMPLE{i}.sampleIds(consistentSampleIds(:,i));
    f_id = fopen([foldername, name, '.txt'], 'w');
    fprintf(f_id, '%d\n', length(vIds));
    for j = 1:length(vIds)
        fprintf(f_id, '%d\n', vIds(j)-1);
    end
    fclose(f_id);
end
