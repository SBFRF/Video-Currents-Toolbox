function totalDataOut  = OCMwrapper(parametersIn, data)
% UNTITLED2 This fuction is a wrapper for video Current Gen to window in
%   the alongshore direction
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
% sort data (raw argus files are not alongshore sorted)
[~, sortIdxY] = sort(data.y);
data.y = data.y(sortIdxY);
data.stack = data.stack(:,sortIdxY);
% define output window and locations
nWindowY = floor((max(data.y) - min(data.y))/dyWindow);  % how many output points in the alongshore
outputYlocations = parametersIn.outputYlocations; % 750:dyOut:950; %nWindowY*dyWindow;
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
%% loop each alongshore window
for yy=1:length(outputYlocations)
    %% interpolate input to get constant alongshore spacing
    slush = 10; % extra data to ensure that interp has more than enough data 
    yminWindow = outputYlocations(yy) - dyWindow/2;  % output min window
    ymaxWindow = outputYlocations(yy) + dyWindow/2;  % output max window
    %fprintf('processing Alongshore location of yFRF = %d m\n', outputYlocations(yy));
    % find the indicies of data of interest for stack and y coordinte
    idxYdataPre = find(data.y >= yminWindow - slush & data.y <=ymaxWindow + slush );
    % remove duplicates from y
    [yInInterp, idxYdata] = unique(data.y(idxYdataPre));
    idxYdata = idxYdataPre(idxYdata);         % index of an index 
    
    % interpolate to constant alongshore spacing
    [yInInterp, tIn] = meshgrid(yInInterp, data.time);      % don't interpolate in time
    % opticalStack coordinates that I want to interpolate to
    yOutCoord = yminWindow:dyInterpOut:ymaxWindow;
    if any(yOutCoord < pierWindow(2)) && any(yOutCoord > pierWindow(1))
        %ismember(pierWindow(1), floor(yOutCoord))
        fprintf('your output location %d is using data within the pier window [%d %d]\n', outputYlocations(yy), pierWindow);
        continue
    end

    [yOutInterp, tOut] = meshgrid(yOutCoord, data.time);
    % now interpolate 
    stackNew = interp2(yInInterp, tIn, data.stack(:,idxYdata), yOutInterp, tOut);

    if plotFlag  % plot data to see what interpolation did to data
        figure();
        sgtitle(sprintf('Output location yFRF %d m',  outputYlocations(yy)))
        ax1 = subplot(311);
        pcolor(data.y(idxYdataPre), data.time, data.stack(:, idxYdataPre));
        shading flat; colormap gray; title('preUnique');
        ax2 = subplot(312);
        pcolor(data.y(idxYdata), data.time(:), data.stack(:, idxYdata));
        shading flat; colormap gray; title('postUnique');
        ax3 = subplot(313);
        pcolor(yOutInterp, tOut, stackNew);
        shading flat; colormap gray; title('postInterp');
        xlabel('yFRF [m]')
        linkaxes([ax1, ax2, ax3], 'x')
        if isempty(fnameOutBase)
            close()
        else
            fnameEnd = sprintf('_temporalWindow_y%dm.png',  outputYlocations(yy));
            saveas(gcf, strcat(fnameOutBase, fnameEnd)); close();
        end
    end
    xy = [yOutInterp(1,:)', yOutInterp(1,:)'];
    xy(:,1) = median(data.x);           % assume cross-shore position is always constant
    timeStart = data.time(1);           % logging to shift output back to real time
    timeIn = data.time - timeStart;     % change starting point in time to zero
    
    %% Run OCM codes
    radonData = myRadonCurrent(stackNew, timeIn, xy(:,2), parametersIn) % Twin, Tstep, plotFlag);  %median(diff(timeIn)),dyInterpOut,dyWindow,stackNew);
    windowedDataOut = videoCurrentGen(stackNew, timeIn, xy, vB, ...
        fkB, Twin, Tstep, plotFlag);
    
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
    fnameEnd = sprintf('_vBarTemporalOutput_y%dm.png',  outputYlocation);
    saveas(gcf, strcat(fnameOutBase, fnameEnd)); close();
end
end