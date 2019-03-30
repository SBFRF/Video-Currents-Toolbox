% load FRF Example file 
load data/1458476940.Sun.Mar.20_12_29_00.GMT.2016.argus02b.cx.vbar200.mat

dataIn.x= XYZ(:,1);     % parse out x Locations 
dataIn.y = XYZ(:,2);    % parse out y Locations 
dataIn.stack = RAW;     % parse out Timestack data 
dataIn.time = T;        % load time 

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
Tstep = 64;        % recommended _____________
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
%% set input parameters for what will be wrapper  
% test below sensitivities 
dyOut = 5;              % output resolution from OCM wrapper 
dyWindow = 20;          % alongshore distance over which to sample stack
dyInterpOut = 0.25;     % interpolated stack resolution in y

%% process one time period 
param.dyOut = dyOut;
param.dyWindow = dyWindow;
param.dyInterpOut = dyInterpOut;
param.vB = vB;
param.fkB = fkB;
param.Tstep = Tstep; 
param.Twin = Twin; 
param.plotFlag = plotFlag; 
param.plotFnameBase = 'figures/base';

myData = OCMwrapper(param, dataIn);
%plot data
time = datetime(myData.t(1,:), 'ConvertFrom', 'posixtime');
figure(); 
pcolor(myData.y, datenum(time),  myData.meanV');
contourf(myData.y, datenum(time),  myData.meanV');
datetick('y', 'HH:MM:SS', 'keepticks')
shading flat;
colorbar();