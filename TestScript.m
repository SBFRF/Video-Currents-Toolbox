% load Example file 
load 1458475140.Sun.Mar.20_11_59_00.GMT.2016.argus02b.cx.vbar125.mat
data.x= XYZ(:,1);   % parse out x Locations 
data.y = XYZ(:,2);  % parse out y Locations 
data.stack = RAW;   % parse out Timestack data 
data.time = T;      % load time 
%% set input parameters for videoCurrentGen function 
%
% set plotFlag:  0 -- no plotting; 1 -- turn plotting on; 2 -- pause
%       plotting with each step
plotFlag = 1;      % turn plotting on with no pause
%
% Set Twin: the time length of the FFT window (in points)
%       For 2 Hz data, Twin = 128 will yield a 64s average current
Twin = 128;        % recommended ____________
%
% Set Tstep: time length to step the window (in points)
%       For 2 Hz data, Tstep = 64 will yield a current estimate every 32 s
Tstep = 32;        % recommended _____________
%
% Set vBounds: [minV maxV], units m/s or vector of desired velocity steps,
%       Set this to empty, [], to use defaults [-3 3]
vB = [];
%
% Set fkBounds = [fmin fmax kmin kmax], vector of frequency and wavenumber
%       bounds energy out of side of these bounds will be set to 0.  Useful
%       to eliminate some of the wave contamination that leaks in.  Set this
%       to empty, [], to use defaults.
fkB = [];          % using defaults 
%% interpolate input to get constant alongshore spacing
ymin = 500;   % take output from starting a local y of 500 
ymax = 1400;  % take output ending at local 1400 m 
dyOut = 0.25; % interpolate to 25 cm spacing

yOutCoord = ymin:dyOut:ymax;  % interpolated output coordinates
% remove duplicates from y
[yIn, idx] = unique(data.y);

% interpolate to constant alongshore spacing 
[yIn, tIn] = meshgrid(yIn, data.time);      % don't interpolate in time 
[yOut, tOut] = meshgrid(yOutCoord, data.time);
stackNew = interp2(yIn, tIn, data.stack(:,idx), yOut, tOut);

if plotFlag  % plot data to see what interpolation did to data 
    figure(); 
    ax1 = subplot(311); 
    pcolor(data.y, data.time, data.stack); shading flat; colormap gray; title('preUnique');
    ax2 = subplot(312); 
    pcolor(data.y(idx), data.time, data.stack(:, idx)); shading flat; colormap gray; title('postUnique');
    ax3 = subplot(313); 
    pcolor(yOut, tOut, stackNew); shading flat; colormap gray; title('postInterp'); 
    linkaxes([ax1, ax2, ax3], 'x')
end

xy = [yOut(1,:)', yOut(1,:)']; xy(:,1) = median(data.x);  % assume cross-shore position is always constant
data.time = data.time - data.time(1); % change starting point in time to zero

%% Run OCM code 
dataOut = videoCurrentGen(stackNew, data.time, xy, vB, ...
        fkB, Twin, Tstep, plotFlag);
%% plot data out
figure();
subplot(511)
errorbar(dataOut.t, dataOut.meanV, dataOut.stdV, '-.')
ylabel('v [m/s]')
xlabel('time [s]')
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

