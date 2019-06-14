w = warning('off');
%% find available data to process
gauge = 150;
months = {'Oct'}; 
% focus on storms in Oct/sept
days = {'13' '14' '15' '16' '17' '18' '19'}; % '26' '27' '28' '29'}; '20', '21', '24''25''12'
%%%%%% loop through below lists
stackResolutionRange = fliplr([0.25, 0.5, 1, 2]); % interpolated stack resolution in y
dyOutRange = fliplr([2.5, 5, 7.5, 10]);           % output resolution from OCM wrapper
dyWindowRange =  fliplr([10, 20, 30]);            % alongshore distance over which to sample stack
%: the time length of the FFT window (in points)
%       For 2 Hz data, Twin = 128 will yield a 64s average current
tWinRange =  fliplr([32, 64, 128, 256]);  % from paper, but expanded
% time length to step the window (in points) Smoothing in time
%       For 2 Hz data, Tstep = 64 will yield a current estimate every 32 s
tStepRange = [256, 128, 64,32];
%% loop through time
for mm=1:length(months)
    for dd=1:length(days)
        day = days{dd};
        month = months{mm};
        globstring = sprintf('/mnt/gaia/peeler/argus/argus02b/2017/cx/*/*%s.%s*vbar%d*', month, day, gauge);
        flist = glob(globstring);  % generate files to process through
        flist=sort(flist);         % sort so files are processed in temporal order
        if isempty(flist)
            fprintf('No Files found for %s %s for gauge %dm\n', month, day, gauge);
            continue
        end 
        for stackRidx=1:length(stackResolutionRange)
            for dyIdx = 1:length(dyOutRange)
                for dyWinIdx = 1:length(dyWindowRange)
                    for TWinIdx = 1:length(tWinRange)
                        for TStepIdx = 1:length(tStepRange)
                            fprintf('\nworking on %s %s gauge %d\nTstep= %d, Twin=%d dyOut=%d, dyWindow=%d, stackRes=%d m\n', ...
                                     month, day, gauge, tStepRange(TStepIdx), tWinRange(TWinIdx), ...
                                     dyOutRange(dyIdx), dyWindowRange(dyWinIdx), stackResolutionRange(stackRidx))
                            dataIn = cell(length(flist),1);
                            param = cell(length(flist),1);
                            dataOut = cell(length(flist),1);
                            parfor i=1:length(flist)
                                try
                                    data = load(flist{i});
                                catch
                                    continue
                                end
                                dataIn{i}.x= data.XYZ(:,1);     % parse out x Locations
                                dataIn{i}.y = data.XYZ(:,2);    % parse out y Locations
                                dataIn{i}.stack = data.RAW;     % parse out Timestack data
                                dataIn{i}.time = data.T;        % load time
                             
                                %% set parameters
                                param{i}.dyOut = dyOutRange(dyIdx);              % output resolution from OCM wrapper
                                param{i}.dyWindow = dyWindowRange(dyWinIdx);          % alongshore distance over which to sample stack
                                param{i}.dyInterpOut = stackResolutionRange(stackRidx);     % interpolated stack resolution in y
                                param{i}.vB = [];   % use defaults
                                param{i}.fkB = [];  % use defaults
                                param{i}.Tstep =  tStepRange(TStepIdx);
                                param{i}.Twin =  tWinRange(TWinIdx);
                                param{i}.plotFlag = 0;
                                param{i}.plotFnameBase = 'figures/base';
                            
                                % make output location list to FRF pipes
                                slop = param{i}.dyOut; % generates 3 outputs at each gauge location (one north, one co-located, and one south)
                                gaugeOutputLocs = [(860-slop):param{i}.dyOut:(860+slop)];% (769-slop):param.dyOut:(769+slop) (945-slop):param.dyOut:(945+slop)];
                                param{i}.outputYlocations = gaugeOutputLocs;
                            
                                %% run OCM on alongshore stack
                                dataOut{i} = OCMwrapper(param{i}, dataIn{i});
                                                       
                            end % i (globbed flist)
                           fname = sprintf('data/processed/OCM_%dm_%s_%s_Tstep%d_Twin%d_dy%.1f_dyWin%d_StackRes%.2f.mat', gauge, month, day, ...
    param{1}.Tstep(1), param{1}.Twin(1), param{1}.dyOut(1), param{1}.dyWindow(1), param{1}.dyInterpOut(1)); 
                            % create a save structure for the day's data 
                        [dataSave] = makeSaveableData(dataOut, param, dataIn, flist);
                        % save  daily file after looping through day's data     
                        save(fname, 'dataSave');                        
                        end % stackRidx
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
nt = squeeze(length(flist)' * nWin);
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


for ii = 0:length(dataOut)-1
    idxMinSave = ii*nWin+1;
    idxMaxSave = idxMinSave + nWin-1;
%     fprintf('idx min: %d idxMax: %d\n', idxMinSave, idxMaxSave);
    dataSave.t(idxMinSave:idxMaxSave, :) = dataOut{ii+1}.t;    
    dataSave.meanV(idxMinSave:idxMaxSave, :) = dataOut{ii+1}.meanV;    
    dataSave.meanI(idxMinSave:idxMaxSave, :) = dataOut{ii+1}.meanI;    
    dataSave.QCspan(idxMinSave:idxMaxSave, :) = dataOut{ii+1}.QCspan;    
    dataSave.stdV(idxMinSave:idxMaxSave, :) = dataOut{ii+1}.stdV;    
    dataSave.prob(idxMinSave:idxMaxSave, :) = dataOut{ii+1}.prob;    
    dataSave.ci(idxMinSave:idxMaxSave, :, :) = dataOut{ii+1}.ci;    
    dataSave.cispan(idxMinSave:idxMaxSave, :) = dataOut{ii+1}.cispan;    
    dataSave.SNR(idxMinSave:idxMaxSave, :) = dataOut{ii+1}.SNR;     
    [dataSave.stack{ii+1, :}] = deal(dataOut{ii+1}.Raw.stack);
    [dataSave.xy{ii+1, :}] = deal(dataOut{ii+1}.Raw.xy);     
end


end
