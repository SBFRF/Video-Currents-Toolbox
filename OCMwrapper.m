function totalDataOut = OCMwrapper(parametersIn, data)
% UNTITLED2 This fuction is a wrapper for video Current Gen to window in
%   the alongshore direction, will sort, interpolate, and section raw
%   argus stack
%
%    INPUTS:
%       parametersIn (struct): must have fields of
%           'dyWindow' -- alongshore window for stack, the width over which
%               to average to get output point [m]
%           'dyOut' -- alongshore spacing of output points [m]  --- REMOVED dependency
%           'dyInterpOut' -- interpolate the stack to this resolution [m]
%           'vB'  -- velocity bounds see videoCurrentGen for description
%           'fkB' -- frequency wave number bounds see videoCurrentGen for description
%           'Twin' -- time window see videoCurrentGen for description
%           'Tstep'-- time step see videoCurrentGen for description
%           'plotFlag' -- plot flage for turning plots on/off see videoCurrentGen for description
%       data (struct): must have fields of:
%           'x' -- cross shore location data (will take median and assume
%              cross-shore uniform
%           'y' -- alongshore location of data (will interpolate to uniform
%              spacing in y at a resolution of dyInterpOut)
%           'stack' -- time stack of imagery data
%           'time' -- epoch time in seconds
%
%    OUTPUTS:
%        a mat structure with the below fields 
%        totalDataOut.t -- UTC time step of data (as input)
%        totalDataOut.y -- the specific alonghosre location of processed
%               data
%        totalDataOut.meanI -- mean intensity of time stack
%        totalDataOut.QCspan -- the 95th percentile minus the 50th 
%              percentile of the timestack histogram, used to measure the 
%              amount of video "texture"
%        totalDataOut.meanV -- mean surface current velocity 
%        totalDataOut.stdV -- the width (std. dev.) of the energy in velocity spectrum
%        totalDataOut.prob --  the probability of the model fit
%        totalDataOut.ci -- the 95% conf. interval around meanV
%        totalDataOut.cispan --  the width of ci
%        totalDataOut.SNR -- model guess at signal to noise ratio
%        totalDataOut.Raw.stack -- interpolated stack actually used to
%               process
%        totalDataOut.Raw.timeIn -- time input for processing 
%        totalDataOut.Raw.xy -- xy vector used for processing 
%        totalDataOut.Raw.Twin -- time window used for processing 
%        totalDataOut.Raw.Tstep -- time step used for processing
%
%  This code requires:
%     Optimization Toolbox (lsqcurvefit.m)
%     Statistics and Machine Learning Toolbox (nlparci.m)
%     Signal Processing Toolbox (bartlett.m)
%
%% Some preprocesssing stuff, 
% unpack input structures 
% dyOut = parametersIn.dyOut;
dyWindow = parametersIn.dyWindow;
dyInterpOut = parametersIn.dyInterpOut;
vB = parametersIn.vB;
fkB = parametersIn.fkB;
Tstep = parametersIn.Tstep; 
Twin = parametersIn.Twin; 
plotFlag = parametersIn.plotFlag;
fnameOutBase = parametersIn.plotFnameBase;
outputYlocations = parametersIn.outputYlocations; 
% fprintf('using restricted ylocations of %d to %d\n', max(outputYlocations), min(outputYlocations));
pierWindow = [495, 520];  % bounds for the pier, will not process data in this window
%% initalize output 
nToutput = floor((size(data.time)-(Twin-Tstep))/Tstep); % number of timesteps for output
outDataDummy = NaN(nToutput(2), length(outputYlocations));
% ouptut structure initalized
totalDataOut.meanI = outDataDummy;
totalDataOut.QCspan = outDataDummy;
totalDataOut.meanV = outDataDummy;
totalDataOut.stdV = outDataDummy;
totalDataOut.prob = outDataDummy;
totalDataOut.ci = outDataDummy;
totalDataOut.cispan = outDataDummy;
totalDataOut.SNR = outDataDummy;
totalDataOut.meanI = outDataDummy;
totalDataOut.QCspan = outDataDummy;
totalDataOut.meanV = outDataDummy;
totalDataOut.stdV = outDataDummy;
totalDataOut.prob = outDataDummy;
totalDataOut.ci =NaN(nToutput(2), 2, length(outputYlocations));
totalDataOut.cispan = outDataDummy;
totalDataOut.SNR = outDataDummy;
%% pre process for loop in alongshore
% sort data (raw argus files are not alongshore sorted and unique)
[data.y, sortIdxY, ~] = unique(data.y, 'sorted');
data.stack = data.stack(:,sortIdxY);
% data.y = data.y(sortIdxY);
% define output window and locations
% nWindowY = floor((max(data.y) - min(data.y))/dyWindow);  % how many output points in the alongshore
%% loop each alongshore window
for yy=1:length(outputYlocations)
    %% interpolate input to get constant alongshore spacing
    slush = 10; % extra data to ensure that interp has more than enough data 
    yminWindow = outputYlocations(yy) - dyWindow/2;  % output min window
    ymaxWindow = outputYlocations(yy) + dyWindow/2;  % output max window
    %fprintf('processing Alongshore location of yFRF = %d m\n', outputYlocations(yy));
    % find the indicies of data of interest for stack and y coordinte
    idxYdata = find(data.y >= yminWindow - slush & data.y <=ymaxWindow + slush);

    % interpolate to constant alongshore spacing
    [yInInterp, tIn] = meshgrid(data.y(idxYdata), data.time);      % don't interpolate in time
    % opticalStack coordinates that I want to interpolate to
    yOutCoord = yminWindow:dyInterpOut:ymaxWindow;
    if any(yOutCoord < pierWindow(2)) && any(yOutCoord > pierWindow(1))
        %ismember(pierWindow(1), floor(yOutCoord))
        fprintf('your output location %d is using data within the pier window [%d %d]\n', outputYlocations(yy), pierWindow);
        continue
    end

    [yOutInterp, tOut] = meshgrid(yOutCoord, data.time);
    % now interpolate and smooth
    interpType = 'nearest'; 
    stackNew = interp2(yInInterp, tIn, data.stack(:,idxYdata), yOutInterp, tOut, interpType);
    stackNew = smoothdata(stackNew, 2, 'sgolay', 'degree', 10);
    
     if plotFlag  % plot data to see what interpolation did to data
        plotTime = datenum(datetime(data.time, 'convertfrom', 'posixtime'));
        figure('renderer', 'painters', 'position', [10, 10, 900, 1000]); 
        sgtitle(sprintf('Output location yFRF %.2f m -- interp: %s',  outputYlocations(yy), interpType))
        ax1 = subplot(6,6, [1:12]);
        imagesc(data.y, plotTime, data.stack); 
        shading flat; colormap gray; hold on
        plot([data.y(idxYdata(1)), data.y(idxYdata(1))], [plotTime(1), plotTime(end)], 'k-', 'linewidth', 3)
        plot([data.y(idxYdata(end)), data.y(idxYdata(end))], [plotTime(1), plotTime(end)], 'k-', 'linewidth', 3)
        datetick('y', 13, 'keeplimits')
        ylabel('time')
        
        ax2 = subplot(6,6, [13,14,15,19,20,21, 25,26,27,31,32,33]);
        imagesc(data.y(idxYdata), plotTime, data.stack(:, idxYdata));
        shading flat; colormap gray; title('Pre-interpolation'); hold on;
        xlim([yOutCoord(1), yOutCoord(end)])
        datetick('y', 13, 'keeplimits')
        xlabel('yFRF [m]')
        ylabel('time')

        ax3 = subplot(6,6,[16,17,18,22,23,24,28,29,30,34,35,36] );
        imagesc(yOutCoord, data.time, stackNew);
        shading flat; colormap gray; title('Post-interpolation');
        xlabel('yFRF [m]')
        set(gca, 'ytick', [])
        set(gca, 'yticklabel',[])
        %linkaxes([ax1, ax2, ax3], 'x')
        
        if isempty(fnameOutBase)
            close()
        else
            fnameEnd = sprintf('_temporalWindow_y%.2gm.png',  outputYlocations(yy));
            saveas(gcf, strcat(fnameOutBase, fnameEnd)); close();
        end
    end
    xy = [yOutInterp(1,:)', yOutInterp(1,:)'];
    xy(:,1) = median(data.x);           % assume cross-shore position is always constant
    timeStart = data.time(1);           % logging to shift output back to real time
    timeIn = data.time - timeStart;     % change starting point in time to zero
    
    %% Run OCM codes
    radonData = myRadonCurrent(stackNew, timeIn, xy(:,2), parametersIn); % Twin, Tstep, plotFlag);  %median(diff(timeIn)),dyInterpOut,dyWindow,stackNew);
    windowedDataOut = videoCurrentGen(stackNew, timeIn, xy, vB, ...
        fkB, Twin, Tstep, plotFlag, fnameOutBase);
    
    %% save alongshore window out to larger output Structure
    % save data
    if plotFlag && ~all(isnan(windowedDataOut.t))
        plotOCM(windowedDataOut, fnameOutBase, outputYlocations(yy)) % plot summary Data out for particular alongshore window
    end
    idxGoodData = ~isnan(windowedDataOut.t); % index of non-NaN'd data 
    
    if any(idxGoodData)
        totalDataOut.t(:, yy) = windowedDataOut.t + timeStart;  % reasign UTC time step, whole time stack will be same, let it overwrite
        totalDataOut.y(yy) = outputYlocations(yy);                    % save the specific point i just processed 
        totalDataOut.meanI(:, yy) = windowedDataOut.meanI;
        totalDataOut.QCspan(:, yy) = windowedDataOut.QCspan;
        totalDataOut.meanV(:, yy) = windowedDataOut.meanV;
        totalDataOut.stdV(:, yy) = windowedDataOut.stdV;
        totalDataOut.prob(:, yy) = windowedDataOut.prob;
        totalDataOut.ci(:,:, yy) = windowedDataOut.ci;
        totalDataOut.cispan(:, yy) = windowedDataOut.cispan;
        totalDataOut.SNR(: , yy) = windowedDataOut.SNR;
        % save input data for comparisons later 
        totalDataOut.Raw.stack{yy} = stackNew;
        totalDataOut.Raw.timeIn{yy} = timeIn;
        totalDataOut.Raw.xy{yy} = xy;
        totalDataOut.Raw.Twin{yy} = Twin;
        totalDataOut.Raw.Tstep{yy} = Tstep;
        totalDataOut.radonV(:,yy) = radonData.v;
        totalDataOut.radonT(:,yy) = radonData.time+timeStart;
        
    else
        totalDataOut.t(:,yy) = windowedDataOut.t + timeStart;   % reasign UTC time step -- will come out NaNs
        totalDataOut.y(yy) = outputYlocations(yy);                    % save the specific point i'm looking for 
        totalDataOut.meanI(:, yy) =  ones('like',windowedDataOut.t) * NaN;
        totalDataOut.QCspan(:, yy) = ones('like',windowedDataOut.t) * NaN;
        totalDataOut.meanV(:, yy) = ones('like',windowedDataOut.t) * NaN;
        totalDataOut.stdV(:, yy) = ones('like',windowedDataOut.t) * NaN;
        totalDataOut.prob(:, yy) = ones('like',windowedDataOut.t) * NaN;
        totalDataOut.ci(:, :, yy) =ones('like',windowedDataOut.t) * NaN;
        totalDataOut.cispan(:, yy) = ones('like',windowedDataOut.t) * NaN;
        totalDataOut.SNR(:, yy) = ones('like',windowedDataOut.t) * NaN;
        % save input data for comparisons later 
        totalDataOut.Raw.stack{yy} = stackNew;
        totalDataOut.Raw.timeIn{yy} = timeIn;
        totalDataOut.Raw.xy{yy} = xy;
        totalDataOut.Raw.Twin{yy} = Twin;
        totalDataOut.Raw.Tstep{yy} = Tstep;
    end
    
end
% disp('remember to remove pier window data, if you change the outputYlocations')

end

%%%%%%%%%%%%%%%%%%%%%%%%%
% function plotImageToAnalyze(param, dataIn)
% 
%     figure()
%     suplot(121)
%     imshow(dataIn.stack/256)
%     
%     
% end
function plotOCM(dataOut, fnameOutBase, outputYlocation)
% Plots OCM output into a five panel plot assumes output of:
%   INPUT:
%       Dataoutput from videocurrentGen function
%
%
%       fnameOut (str): file name for saving plot if interested if empty
%           will not save
%
figure();
subplot(511)
errorbar(dataOut.t, dataOut.meanV, dataOut.stdV, '-')
ylabel('v [m/s]')
subplot(514)
plot(dataOut.t, dataOut.meanI, '.')
ylabel('mean Intensity')
subplot(513)
plot(dataOut.t, dataOut.SNR, '.')
ylabel('SNR')
subplot(512)
plot(dataOut.t, dataOut.QCspan, '.')
ylabel('QC span')
subplot(515)
plot(dataOut.t, dataOut.cispan, '.')
ylabel('Confidence Interval Span')
xlabel('time [s]')

% save if data are available
if isempty(fnameOutBase)
    close()
else
    fnameEnd = sprintf('_vBarTemporalOutput_y%.2gm.png',  outputYlocation);
    saveas(gcf, strcat(fnameOutBase, fnameEnd)); close();
end
end