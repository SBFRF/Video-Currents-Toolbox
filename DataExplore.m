w = warning('off');
%% find available data to process
gauge = 200;
months = {'Oct'}; 
year = '2015';
% focus on storms in Oct/sept
days = string([13:31]); 
stackInterpResolution = 0.02; % 2 cm is peak alongshore resolution of the image 
targetYFRF = 860;   %alongshore location of postion to process
%% lists to loop through
%%%%%% loop through below lists
stackResolutionRange = stackInterpResolution;    % interpolated stack resolution in y
dyStepRange = fliplr([5, 10]);                        % output resolution from OCM wrapper
dyWindowRange =  fliplr([10, 20, 30, 40, 50]);   % alongshore distance over which to sample stack
%: the time length of the FFT window (in points)
%       For 2 Hz data, Twin = 128 will yield a 64s average current
tWinRange =  fliplr([32, 64, 128, 256]);  % from paper, but expanded
% time length to step the window (in points) Smoothing in time
%       For 2 Hz data, Tstep = 64 will yield a current estimate every 32 s
tStepRange = [256, 64];
% radial filter step
radialStep = [10, 20, 30];   % radial filter used for radon transfer
%% loop through time
for mm=1:length(months)
    for dd=1:length(days)
        day = days{dd};
        month = months{mm};
        globstring = sprintf('/mnt/gaia/peeler/argus/argus02b/%s/cx/*/*%s.%s*vbar%d*', year, month, day, gauge);
        flist = sort(glob(globstring));  % generate files to process through
%         flist=sort(flist);         % sort so files are processed in temporal order
        if isempty(flist)
            fprintf('No Files found for %s %s for gauge %dm\n', month, day, gauge);
            continue
        end 
        for TWinIdx = 1:length(tWinRange)
           for dyWinIdx = 1:length(dyWindowRange)
                for TStepIdx = 1:length(tStepRange)
                    for dyIdx = 1:length(dyStepRange)
                            fprintf('\nworking on %s %s gauge %d\nTstep= %d, Twin=%d dyOut=%d, dyWindow=%d\n', ...
                                     month, day, gauge, tStepRange(TStepIdx), tWinRange(TWinIdx), ...
                                     dyStepRange(dyIdx), dyWindowRange(dyWinIdx))
                            dataIn = cell(length(flist),1);
                            param = cell(length(flist),1);
                            dataOut = cell(length(flist),1);
                            parfor i=1:length(flist)
                                try
                                    data = load(flist{i});
                                catch; continue; end 
                                dataIn{i}.x= data.XYZ(:,1);     % parse out x Locations
                                dataIn{i}.y = data.XYZ(:,2);    % parse out y Locations
                                dataIn{i}.stack = data.RAW;     % parse out Timestack data
                                dataIn{i}.time = data.T;        % load time
                             
                                %% set parameters
                                param{i}.dyStep = dyStepRange(dyIdx);              % output resolution from OCM wrapper
                                param{i}.dyWindow = dyWindowRange(dyWinIdx);          % alongshore distance over which to sample stack
                                param{i}.dyInterpOut = stackInterpResolution;     % interpolated stack resolution in y
                                param{i}.Tstep =  tStepRange(TStepIdx);
                                param{i}.Twin =  tWinRange(TWinIdx);
                                param{i}.radialFilterThresh = radialStep;
                                param{i}.plotFlag = 0;
                                param{i}.plotFnameBase = 'figures/base';
                                param{i}.vB = [];   % use defaults
                                param{i}.fkB = [];  % use defaults                            
                                % make output location list to FRF pipes
%                                 slop = param{i}.dyStep(dyIdx);  % generates 3 outputs at each gauge location (one north, one co-located, and one south)
                                gaugeOutputLocs = [targetYFRF]; % (769-slop):param.dyOut:(769+slop) (945-slop):param.dyOut:(945+slop)];
                                param{i}.outputYlocations = gaugeOutputLocs;
                            
                                %% run OCM on alongshore stack
                                dataOut{i} = OCMwrapper(param{i}, dataIn{i});
                                                       
                            end % i (globbed flist)
           
                           fname = sprintf('data/processed/OCM_%dm_%s_%s_tWin%d_dyWin%d_tStep%d_dyStep%d.mat', gauge, month, day, ...
                                                       param{1}.Twin(1), param{1}.dyWindow(1), param{1}.Tstep(1),  param{1}.dyStep(1)); 
                            % create a save structure for the day's data 
                        [dataSave] = makeSaveableData(dataOut, param, dataIn, flist);
                        % save  daily file after looping through day's data     
                        save(fname, 'dataSave');                        
                        
                    end % dyOutRange
                end % dyWindowRange
            end % tWinRange
        end % tStepRange
    end % months
end % days

function dataSave = makeSaveableData(dataOut, param, dataIn, flist)
%
% This function makes saveable data from the pre-processor
%%
ny = length(param{1,1}.outputYlocations);
nWin = floor((size(dataIn{1}.time,2)-(param{1}.Twin-param{1}.Tstep))/param{1}.Tstep);
nRadon = max(size(dataOut{1}.radonData(1).time'));
nt = squeeze(length(flist)' * nWin);
nRadonSmooths = size(dataOut{1}.radonData,1);
% initalize variables
dataSave.param = param{1};                      % don't bother looping stays static
dataSave.y =  param{1}.outputYlocations;        % don't bother looping stays static
dataSave.t =   NaN(nt, ny);
dataSave.meanV =  NaN(nt, ny);
dataSave.meanI =  NaN(nt, ny);
dataSave.QCspan =  NaN(nt, ny);
dataSave.stdV =  NaN(nt, ny);
dataSave.prob =  NaN(nt, ny);
dataSave.ci =  NaN(nt, 2, ny);
dataSave.cispan =  NaN(nt, ny);
dataSave.SNR =  NaN(nt, ny);
dataSave.stack =  cell(nt, ny);  
dataSave.xy =  cell(nt, ny);  
dataSave.rd.time = NaN(nRadon, ny, nRadonSmooths);
dataSave.rd.v = NaN(nRadon, ny, nRadonSmooths);
dataSave.rd.stdI = NaN(nRadon, ny, nRadonSmooths);
dataSave.rd.QCspan = NaN(nRadon, ny, nRadonSmooths);
dataSave.rd.meanI = NaN(nRadon, ny, nRadonSmooths);
dataSave.rd.AngPixlIntensDensity = NaN(nRadon, ny, nRadonSmooths, 180);
dataSave.rd.radialFilterThresh = NaN(nRadon, ny, nRadonSmooths);
       
for ii = 0:length(dataOut)-1
    idxMinSave = ii*nWin+1;
    idxMaxSave = idxMinSave + nWin-1;

    dataSave.t(idxMinSave:idxMaxSave, :) = dataOut{ii+1}.t;    
    dataSave.meanV(idxMinSave:idxMaxSave, :) = dataOut{ii+1}.meanV;    
    dataSave.meanI(idxMinSave:idxMaxSave, :) = dataOut{ii+1}.meanI;    
    dataSave.QCspan(idxMinSave:idxMaxSave, :) = dataOut{ii+1}.QCspan;    
    dataSave.stdV(idxMinSave:idxMaxSave, :) = dataOut{ii+1}.stdV;    
    dataSave.prob(idxMinSave:idxMaxSave, :) = dataOut{ii+1}.prob;    
    dataSave.ci(idxMinSave:idxMaxSave, :, :) = dataOut{ii+1}.ci;    
    dataSave.cispan(idxMinSave:idxMaxSave, :) = dataOut{ii+1}.cispan;    
    dataSave.SNR(idxMinSave:idxMaxSave, :) = dataOut{ii+1}.SNR;     
    %[dataSave.stack{ii+1, :}] = deal(dataOut{ii+1}.Raw.stack);
    %[dataSave.xy{ii+1, :}] = deal(dataOut{ii+1}.Raw.xy);     
    
    for kk = 1:size(dataOut{ii+1}.radonData,1)
       idxMaxSaveRD = idxMinSave +  length(dataOut{ii+1}.radonData(kk).time')-1;
       dataSave.rd.time(idxMinSave: idxMaxSaveRD, :, kk) = dataOut{ii+1}.radonData(kk).time';
       dataSave.rd.v(idxMinSave: idxMaxSaveRD, :, kk) = dataOut{ii+1}.radonData(kk).v;
       dataSave.rd.stdI(idxMinSave: idxMaxSaveRD, :, kk) = dataOut{ii+1}.radonData(kk).stdI;
       dataSave.rd.QCspan(idxMinSave: idxMaxSaveRD, :, kk) = dataOut{ii+1}.radonData(kk).QCspan; 
       dataSave.rd.meanI(idxMinSave: idxMaxSaveRD, :, kk) = dataOut{ii+1}.radonData(kk).meanI;
       dataSave.rd.AngPixlIntensDensity(idxMinSave: idxMaxSaveRD, :, kk, :) = dataOut{ii+1}.radonData(kk).AngPixlIntensDensity;
       dataSave.rd.radialFilterThresh(idxMinSave: idxMaxSaveRD, :, kk) = dataOut{ii+1}.radonData(kk).radialFilterThresh;
    end

end


end
