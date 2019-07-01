% RadonCurrent_20141129.m
function [CRadmoy]=RadonCurrent_20141129(dt,dx,Wx,In)
%Inputs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%In: Input spatio-temporal matrix, format: (nt,nx)
% In=ZHF;
%pdx: Spatial resolution degradation (in point)
pdx=5;
%pdt: Temporal resolution degradation (in point)
pdt=1;
%Wx: Window size in space for Radon computation (sensity test : better if WX > 10)
% Wx= 30; %
%Wt: window size in time for Radon computation (in sec)
%Wt=30;
%freq_t : temporal frequency (1/dt) (taking into account for pdt)
freq_t = 1./dt;%
%freq_x : spatial frequency (1/dx)
freq_x = 1./(dx.*pdx);%(here 2 pixel/m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright 2017 Rafael Almar (IRD, France)- rafael.almar@ird.fr

%Outputs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% an: sorties brutes radon (angle en °) mieux si entre 30 et 70
% Hx: élévation des crêtes détectées pour le calcul de la célérité
% Hm: élévation des creux détectées pour le calcul de la célérité
% Tim: crêtes détectées pour le calcul de la célérité
% C: célérités des crêtes détectées
%Outputs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Wx=round(Wx./pdx);

M=In(1:pdt:length(In(:,1)),1:pdx:length(In(1,:)));

clear Mat VER XHFR an Tim Hx Hm C CC2 Tan
length(Wx+1:Wx:size(M,2)-(Wx+1))
for ix=Wx+1:Wx:size(M,2)-(Wx+1)%Boucle sur les positions X


MR=M(1:length(M(:,1)),ix-Wx:ix+Wx);%selection de la fenetre

% [fr1 gtf1]=lmax(smooth(MR(:,round(size(MR,2)/2)),10));%Recherche des crêtes
% [fr2 gtf2]=lmax((MR(:,round(size(MR,2)/2))));%Recherche des crêtes
% 
% clear gtf fr
% for ip=1:length(gtf1)
%  [fre hgy]=min(abs(gtf1(ip)-gtf2));
%  gtf(ip)=gtf2(hgy);
%  fr(ip)=fr2(hgy);
% end
% [kp tp]=sort(fr);
% fr=kp(round(length(fr)/3:length(fr)));
% gtf=gtf(tp(round(length(kp)/3:length(kp))));
% [gtf tp]=sort(gtf);
% fr=fr(tp);

% [frfr gtffr]=lmin(smooth(MR(:,round(size(MR,2)/2)),5));%Recherche des creux

nt=size(MR,1);
R=radon(detrend(double(MR'))',1:1:180);%transformée de Radon
% figure(128);clf;pcolor(R);shading flat
% pause
tr=abs(cosd(90:270));
nr=size(R,1);
amp=nr/nt;
iang=1:1:180;

k=nt;%-round(nt./2);
% trk=round(nr/2+tr*((k-(nr+1)/2))*amp);
trk=floor((size(MR,1))/2)-floor((  0*cosd(iang)+  ((size(MR,1))/2-k*amp)*sind(iang)));trk=trk-min(trk);
res=(nt*dt)./(trk.*2);

R2=R;
for i=iang
   R2(:,i)=R(:,i)-smooth(R(:,i),round(1+20./(res(i)))); 
end

[frd a2]=max(std(R2(round(size(R2,1)/4:3*size(R2,1)./4),:)));%calcul celerite moyenne

if length(freq_x)==1
C2=(1/mean(freq_x))/(tand(90-a2)*(1/mean(freq_t)));% celerite moyenne
else
C2=(1/mean(freq_x(ix-Wx:ix+Wx)))/(tand(90-a2)*(1/mean(freq_t)));% celerite moyenne

end
%calcul celerite vague a vague
% figure(16);clf;pcolor(R);shading flat; hold on;
%  for k=1:length(gtf)
%      try
% 
% % hold on;plot(1:length(trk),trk,'k')
% vag=[];
% for i=-5:5% on prend les résultats sur les points autours pour moins de bruit
%     vag=[vag smooth(smooth(diag(R(trk(ang)+i,ang)),3),10)];
% end
% res=0.1;
% vec=max(vag');vec=smooth(smooth(interp1(1:length(vec),vec,1:res:length(vec)),15),32);
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
% qual(co,ix)=frd;%Critére qualité
% 
% if length(freq_x)==1
% CRad(co,ix)=(1/freq_x)/(tand(90-g11)*(1/freq_t));%celerites (angle en degree) 
% else
% CRad(co,ix)=(1/mean(freq_x(ix-Wx:ix+Wx)))/(tand(90-g11)*(1/freq_t));%celerites (angle en degree)    
% end
% 
% TCRad(co,ix)=nt-k;%celerites (angle en degree)  
% % end
% %  end
%  end
%  end
 
%  mean(VELHF(:,ix))% celerite données


%  [Tok indok]=sort(XHF(:,ix));
% figure(378);clf
% plot(Tan(ix,:),an(ix,:),'k');hold on;plot(Tok,VELHF(indok,ix),'r')
% pause
 
CRadmoy(ix)=C2; %celerite moyenne
end
CRadmoy=interp1(pdx.*(Wx+1:Wx:size(M,2)-(Wx+1)),CRadmoy(Wx+1:Wx:size(M,2)-(Wx+1)),1:size(In,2));
