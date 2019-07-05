% RadonCurrent_20141129.m
function [CRadmoy]=myRadonCurrent(In, timeIn, y, parametersIn)
%Inputs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   In: input stack (nt, nx)
%   timeIn: time input (seconds since start of stack) .. ie (0..2048)
%   xy: x and y pixal locations of data of stack (size=(nx,2))
%   tWin: temporal window to do analysis over 
%   tStep: temporal step to include (not currently used)
%   varargin: plotting 1 = true, 0 = false 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright 2017 Rafael Almar (IRD, France)- rafael.almar@ird.fr
%
%   Citation:
%   Almar, Rafael, et al. "On the use of the Radon transform to estimate
%   longshore currents from video imagery." Coastal Engineering 114
%   (2016): 301-308.
%
%    Requirements:
%       Curve Fitting Toolbox (smooth function)
%            used move mean when not available
%Outputs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% an: Radon outputs better if (angle in degrees) is between 30 and 70 degrees
% Hx: elevation of peaks detected for the calculation of celerity
% Hm: elevation of the hollows detected for the calculation of celerity
% Tim: peaks detected for calculating the clerity
% C: peaks of detected crests
%Outputs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
plotFlag = parametersIn.plotFlag;
fnameOutBase = parametersIn.plotFnameBase;
tStep = parametersIn.Tstep;
tWin = parametersIn.Twin;

%pdx: Spatial resolution degradation (in point)  -- sub sampling 
pdx=1;
dx = median(diff(y));  % calculate the pixel size in size over ground [m]
%freq_x : spatial frequency (1/dx)
freq_x = 1./(dx.*pdx);  %(here 2 pixel/m)
%Wx: Window size in space for Radon computation (sensity test : better if WX > 10)
Wx = abs(y(1) - y(end));  % alongshore window size in meters
% 2nd: time -- time of input
dt = median(diff(timeIn));
%pdt: Temporal resolution degradation (in point)
pdt=1;
%freq_t : temporal frequency (1/dt) (taking into account for pdt)
freq_t = 1./dt;  %
%Wt: window size in time for Radon computation (in sec) %Wt=30;
iang=1:1:180;  % define angular resolution
disp('need to incorporate tstep, right now operates as 0% overlap')
% this subsamples in x and t (as defined above by pdx, pdt)
M=In(1:pdt:length(In(:,1)),1:pdx:length(In(1,:)));
clear Mat VER XHFR an Tim Hx Hm C CC2 Tan

for tt = 1:tWin:size(squeeze(timeIn),2)-(tWin+1)   % loop on time window
    for ix=Wx+1:Wx:size(M,2)-(Wx+1)       % loop on x positions
       
        MR=M(tt:tt+tWin, ix-Wx:ix+Wx);                      % window selection
        %     try
        %      
        %     catch
        %     [fr1 gtf1]=lmax(movmean(MR(:,round(size(MR,2)/2)),10));  % search for wave crests
        %    end 
        %     [fr2 gtf2]=lmax((MR(:,round(size(MR,2)/2))));            % search for wave crests
        %
        %     clear gtf fr
        %     for ip=1:length(gtf1)
        %         [fre hgy]=min(abs(gtf1(ip)-gtf2));
        %         gtf(ip)=gtf2(hgy);
        %         fr(ip)=fr2(hgy);
        %     end
        %     [kp tp]=sort(fr);
        %     fr=kp(round(length(fr)/3:length(fr)));
        %     gtf=gtf(tp(round(length(kp)/3:length(kp))));
        %     [gtf tp]=sort(gtf);
        %     fr=fr(tp);
        %       try
        %      [frfr gtffr]=lmin(smooth(MR(:,round(size(MR,2)/2)),5)); % search for
        %     hollows (troughs?)
        %     catch 
        %     [frfr gtffr]=lmin(movmean(MR(:,round(size(MR,2)/2)),5)); % search for
        %     hollows (troughs?)
        %      end
        %      
        %
        nt=size(MR,1);
        [R, Xp] =radon(detrend(double(MR'))', iang); % radon transform
        if plotFlag == 1  % make plot for QA/QC to see what's happening with transform 
            figure('Renderer', 'painters', 'Position', [10 10 1500 500]);  clf;
            subplot(131)
            imshow(MR/256);ylabel('time'); xlabel('pixels alongshore')
            subplot(132)
            pcolor(iang, Xp, R); shading flat; ylabel('radial distance');
            xlabel('angle'); colorbar()
            subplot(133)
            plot(iang, sum(R,1)); xlabel('angle'); ylabel('integrated radial concentration')

            if isempty(fnameOutBase)
                pause
                close()
            else
                fnameEnd = sprintf('_RadonTransformQAQC_%ds_%dm.png', tt, ix*Wx);
                saveas(gcf, strcat(fnameOutBase, fnameEnd)); close();
            end
        end
        tr=abs(cosd(90:270));
        nr=size(R,1);
        amp=nr/nt;
        
        k=nt;   % -round(nt./2);
        % trk=round(nr/2+tr*((k-(nr+1)/2))*amp);
        trk=floor((size(MR,1))/2)-floor((0*cosd(iang) + ((size(MR,1))/2-k*amp)*sind(iang)));
        trk=trk-min(trk);
        res=(nt*dt)./(trk.*2);
        
        R2=R;
        for i=iang
            try
                R2(:,i)=R(:,i)-smooth(R(:,i),round(1+20./(res(i))));
            catch
                R2(:,i)=R(:,i)-movmean(R(:,i),round(1+20./(res(i))));
            end

        end
        
        [frd, a2] = max(std(R2(round(size(R2,1)/4:3*size(R2,1)./4),:)));  %mean celerity calculation
        
        if length(freq_x)==1
            C2=(1/mean(freq_x))/(tand(90-a2)*(1/mean(freq_t)));   % mean celerity
        else
            C2=(1/mean(freq_x(ix-Wx:ix+Wx)))/(tand(90-a2)*(1/mean(freq_t)));% mean celerity
        end
        %calculate celerity wave to wave
        % figure(16);clf;pcolor(R);shading flat; hold on;
        %  for k=1:length(gtf)
        %      try
        %
        % % hold on;plot(1:length(trk),trk,'k')
        % vag=[];
        % for i=-5:5 % we take the results on the points around for less noise
        %      try
        %     vag=[vag smooth(smooth(diag(R(trk(ang)+i,ang)),3),10)];
        %      catch
        %     vag=[vag movmean(movemean(diag(R(trk(ang)+i,ang)),3),10)];
        %      end
        % end
        % res=0.1;
        % vec=max(vag');
        % try
        % vec=smooth(smooth(interp1(1:length(vec),vec,1:res:length(vec)),15),32);
        % catch
        % vec=movmean(movmean(interp1(1:length(vec),vec,1:res:length(vec)),15),32);
        % [frd g11]=max(vec(a2/res-10/res:min([a2/res+10/res length(vec)])));g11=g11+(a2-10)/res;
        % g11=g11*res+1;
        %
        % % figure(39);clf;pcolor(vag');shading flat;axis equal%caxis([-160 160]);
        % % hold on;plot(g11,5,'ko')
        % % pause
        % %
        % % if abs((1/freq_x)/(tand(90-g11)*(1/freq_t))-C2)<1
        %
        % % figure(37);plot(vec);
        % % hold on;plot(g11,vec(g11),'rsq');
        % % pause
        %
        %         co=co+1;
        % qual(co,ix)=frd;%Crit�re qualit�
        %
        % if length(freq_x)==1
        % CRad(co,ix)=(1/freq_x)/(tand(90-g11)*(1/freq_t));%celerites (angle in degrees)
        % else
        % CRad(co,ix)=(1/mean(freq_x(ix-Wx:ix+Wx)))/(tand(90-g11)*(1/freq_t));%celerites (angle in degree)
        % end
        %
        % TCRad(co,ix)=nt-k;%celerites (angle en degree)
        % % end
        % %  end
        %  end
        %  end
        
        %  mean(VELHF(:,ix)) % celerity data
        
        %  [Tok indok]=sort(XHF(:,ix));
        % figure(378);clf
        % plot(Tan(ix,:),an(ix,:),'k');hold on;plot(Tok,VELHF(indok,ix),'r')
        % pause
        
        CRadmoy(ix, tt)=C2; % mean celerity
    end
end
% interpolate output before returning
CRadmoy=interp1(pdx.*(Wx+1:Wx:size(M,2)-(Wx+1)), CRadmoy(Wx+1:Wx:size(M,2)-(Wx+1)), 1:size(In,2));
