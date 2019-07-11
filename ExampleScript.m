% load FRF Example file 
load data/1446152340.Thu.Oct.29_20_59_00.GMT.2015.argus02b.cx.vbar200.mat
% ADV=load('ExampleScriptObs.mat');
%%% measured velocity (at y=945 x=300)  is -0.3191844 m/s  from aquadopp

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
%% for Radon Current method
radialFilterThresh = 5; % radial distance in pixels to filter (figure 2 in almar et al 2016)

%% set input parameters for what will be wrapper  
% test below sensitivities 
dyOut = 10;              % output resolution from OCM wrapper 
dyWindow = 50;          % alongshore distance over which to sample stack
dyInterpOut = 0.025;     % interpolated stack resolution in y
yFRF = 800;             % alongshore location of interest (just do one for speed's sake) 
filePrefix = join(['figures/ExampleScript/Single_' string(radialFilterThresh)],'');
%% process one time period
% setup input structure
param.dyOut = dyOut;
param.dyWindow = dyWindow;
param.dyInterpOut = dyInterpOut;
param.vB = vB;
param.fkB = fkB;
param.Tstep = Tstep; 
param.Twin = Twin; 
param.plotFlag = plotFlag;
param.radialFilterThresh = radialFilterThresh;
%param.dyOut = 0.5; % m resolution for interpolation
slop = param.dyOut; % generates 3 outputs at each gauge location (one north, one co-located, and one south)
gaugeOutputLocs = [yFRF]; %(yFRF-slop):param.dyOut:(yFRF+slop)];% (769-slop):param.dyOut:(769+slop) (945-slop):param.dyOut:(945+slop)];
param.outputYlocations = gaugeOutputLocs;
param.plotFnameBase = filePrefix;

%% run the wrapper code, runs both methods 
vb = OCMwrapper(param, dataIn);

%% plot data
%%%%%%%%%%%%%%% load obs 
url="https://chldata.erdc.dren.mil/thredds/dodsC/frf/projects/bathyduck/data/BathyDuck-ocean_currents_p23_201510.nc";
ADV.yLoc = ncread(url, 'alongVelYloc');
ADV.xLoc = ncread(url, 'alongVelXloc');
ADV.time = datenum(datetime(ncread(url, 'time'), 'convertfrom', 'posixtime'));
ADV.v = ncread(url, 'alongVel');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yLoc = 1;
timeVB = datenum(datetime(vb.t(:, yLoc), 'ConvertFrom', 'posixtime'));
radonTime = datenum(datetime(vb.radonT(:, yLoc), 'ConvertFrom', 'posixtime'));
[idxObsR, idxRadon] = timeMatchOCM(ADV.time, radonTime);
[idxObsV, idxVBar] = timeMatchOCM(ADV.time, timeVB);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(); 
subplot(121)
plot(timeVB,  vb.meanV(:, yLoc), '.', 'displayname', 'vBar Method'); hold on;
plot(radonTime,  vb.radonV(:, yLoc), 'g.','displayname', 'Radon Method');
plot(ADV.time, ADV.v, 'rx--', 'displayname', 'Observations')
xlim([min(radonTime)-0.01, max(radonTime)+0.01])
legend()
datetick('x', 'HH:MM:SS', 'keepticks')
subplot(122)
plot(ADV.v(idxObsV), vb.meanV(idxVBar, yLoc),'.', 'displayname', 'Vbar'); hold on;
plot(ADV.v(idxObsV), vb.radonV(idxRadon, yLoc),'.', 'displayname', 'Radon')
plot([-3,3], [-3,3], 'k-') %, 'dislayname', 'unity')
legend()
ylabel('camera derived velocities [m/s]')
xlabel('observed velocities');

function [idxObs, idxOptical] = timeMatchOCM(obsTime, opticalTime)
%    """
%    function looks though optical time for the closest value that is below the time step of the optical time
%    the obs sample period.  It does this so multiple OCM measurements can be compared to the observations
%    Args:
%        obsTime:  epochtime (or some other numeric)
%        opticalTime: epoch time (or some other numeric) -- this is master time to be matched to
%
%
%% 
rc=1;
obsSamplePeriod = median(diff(obsTime));
for idxOT = 1:length(opticalTime)
    oTime = opticalTime(idxOT);
    [~, idxMaybe] = min(abs(oTime - obsTime));
    if abs(oTime - obsTime(idxMaybe)) < obsSamplePeriod
        idxObs(rc) = idxMaybe;
        idxOptical(rc) = idxOT;
    end
    rc=rc+1;
end
end


