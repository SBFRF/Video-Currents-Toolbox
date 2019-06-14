% script to remove stack from saved data
globstring = sprintf('data/processed/*.mat');
flist = glob(globstring);
for i=1:length(flist);
    file = flist{i};
    try
        data=load(file);
        fprintf('loaded file %s\n', file)
        try
            dataSave = rmfield(data.dataSave, 'stack');
            save(file, 'dataSave')
        catch
        end
    catch
    end
end