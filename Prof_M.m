% function Prof_M(stp) 
% Plot tracer profiles from the time series with step 'stp'
clear

stp = 40;
if exist('stp')==0;
disp('Prof should be called with a step variable, e.g. Prof(2)')
return
end

%-----
load('OutThilda_M')

% Data ordered accordingly:
% sLL(tracer,level,time) contains low latitude ocean data
% sHL(tracer,level,time) contains high latitude ocean data
% sAT(tracer,zone,time) contains atmospheric data
% st(time) contains model time
%-----
% Ocean tracers are: 
%                     1:   - temperature [oC]
%                     2:   - salinity []
%                     3:   - phosphate [mol/m3
%                     4:   - DIC [mol/m3] 
%                     5:   - DI^13C [mol/m3]
%                     6:   - DI^14C [mol/m3]
%                     7:   - alkalinity [eq/m3] 
%                     8:   - oxygen  [mol/m3] 
%                     9:   - Isotope ration d18O [-]
%                     10:  - CH4 [mol/m3]
%                     11:  - CH4^13C  [mol/m3]
%                     12:  - nitrate [mol/m3] 
%                     131:  - NH_4 [mol/m3] 
%                     14:  - H_2S [mol/m3] 
% Zones are 
% 1:low latitude (LL, red lines), 2:high latitude (HL, blue lines) 


%-----
% Get parameter values
ParVal_M    % Activate global parameters
global nto dm d n R13pdb 
%-----

%-----
zcent= [dm/2 dm+(d:d:(n-1)*d)-d/2]/1e3;   % Vertical center of boxes
%-----

% Calculate delta13 and Delta 14 values in permil for ocean and atmosphere
%-----

d13LLo = squeeze( (sLL(5,:,:)./sLL(4,:,:)/R13pdb-1)*1e3 );
d13HLo = squeeze( (sHL(5,:,:)./sHL(4,:,:)/R13pdb-1)*1e3 );
d13MLLo = squeeze( (sLL(11,:,:)./sLL(10,:,:)/R13pdb-1)*1e3 );
d13MHLo = squeeze( (sHL(11,:,:)./sHL(10,:,:)/R13pdb-1)*1e3 );

O2LL =squeeze(sLL(8,:,:));
O2HL =squeeze(sHL(8,:,:));
for i = 1:n
    for j=1:length(st)
      if O2LL(i,j)<0
         O2LL(i,j)=0;
      end
     if O2HL(i,j)<0
        O2HL(i,j)=0;
     end
    end
end

%-----
% Loop through the time series with step 'stp'

% disp(strcat([int2str(round(length(st)/stp)) '  images will be shown with short intervals.']))
% disp(strcat(['The model time between images is ' int2str(round(st(2)-st(1))*stp) ' years.' ]))
% if length(st)/stp>10 disp('!! Increase input step variable to Prof(step) for fewer images. !!'); end

% for i=1:stp:length(st)
for i=1

% Iterate Carbonate system for the vertical profile of 
% the CO3 concentration and saturation state 
%-----
[LLCO3,LLCO3s,LLCO2,LLHCO3,LLHp]=CarSysPres_M(sLL(:,:,i));
[HLCO3,HLCO3s,HLCO2,HLHCO3,HLHp]=CarSysPres_M(sHL(:,:,i));
%-----
%-----
% Plotting
%-----

figure("Position",[0,0,750,700])
ha = tight_subplot(2,5,[0.06,0.06],[0.06,0.02],[0.06,0.02]);

axes(ha(1))
plot([sLL(1,1,i) sLL(1,1,i) sLL(1,2:n,i)],[0 dm/1e3 zcent(2:end)],'r-');axis ij; hold on
plot([sHL(1,1,i) sHL(1,1,i) sHL(1,2:n,i)],[0 dm/1e3 zcent(2:end)],'b-');
ylabel('z (km)')
xlabel('T (^oC)')
axis([-2 50 0 n*d*1e-3])
ax(1)=gca;

axes(ha(2))
plot([sLL(2,1,i) sLL(2,1,i) sLL(2,2:n,i)],[0 dm/1e3 zcent(2:end)],'r-');axis ij; hold on
plot([sHL(2,1,i) sHL(2,1,i) sHL(2,2:n,i)],[0 dm/1e3 zcent(2:end)],'b-');
ylabel('z (km)')
xlabel('S')
axis([31 40 0 n*d*1e-3])
ax(2)=gca;

axes(ha(3))
plot(10^3*[sLL(10,1,i) sLL(10,1,i) sLL(10,2:n,i)],[0 dm/1e3 zcent(2:end)],'r-');axis ij; hold on
plot(10^3*[sHL(10,1,i) sHL(10,1,i) sHL(10,2:n,i)],[0 dm/1e3 zcent(2:end)],'b-');
ylabel('z (km)')
xlabel('Methane (mmol m^{-3})')
axis([0 max(max(squeeze(10^3*sLL(10,:,:)))) 0 n*d*1e-3])
ax(3)=gca;

axes(ha(4))
plot(10^3*[sLL(12,1,i) sLL(12,1,i) sLL(12,2:n,i)],[0 dm/1e3 zcent(2:end)],'r-');axis ij; hold on
plot(10^3*[sHL(12,1,i) sHL(12,1,i) sHL(12,2:n,i)],[0 dm/1e3 zcent(2:end)],'b-');
ylabel('z (km)')
xlabel('Nitrate (mmol m^{-3})')
axis([0 max(max(squeeze(10^3*sLL(12,:,:)))) 0 n*d*1e-3])
ax(4)=gca;

axes(ha(5))
plot([O2LL(1,i) O2LL(1,i) O2LL(2:n,i)'],[0 dm/1e3 zcent(2:end)],'r-');axis ij; hold on
plot([O2HL(1,i) O2HL(1,i) O2HL(2:n,i)'],[0 dm/1e3 zcent(2:end)],'b-');
ylabel('z (km)')
xlabel('O_2 (mol m^{-3})')
axis([-.05 .4 0 n*d*1e-3])
ax(5)=gca;

axes(ha(6))
plot([sLL(4,1,i) sLL(4,1,i) sLL(4,2:n,i)],[0 dm/1e3 zcent(2:end)],'r-');axis ij; hold on
plot([sHL(4,1,i) sHL(4,1,i) sHL(4,2:n,i)],[0 dm/1e3 zcent(2:end)],'b-');
plot([sLL(7,1,i) sLL(7,1,i) sLL(7,2:n,i)],[0 dm/1e3 zcent(2:end)],'r:');axis ij;
plot([sHL(7,1,i) sHL(7,1,i) sHL(7,2:n,i)],[0 dm/1e3 zcent(2:end)],'b:');
ylabel('z (km)')
xlabel('DIC (mol m^{-3}) & Alk (eq m^{-3})')
axis([1.7 2.4 0 n*d*1e-3])
ax(6)=gca;

axes(ha(7))
plot([d13LLo(1,i) d13LLo(1,i) d13LLo(2:end,i)'],[0 dm/1e3 zcent(2:end)],'r-');axis ij; hold on
plot([d13HLo(1,i) d13HLo(1,i) d13HLo(2:end,i)'],[0 dm/1e3 zcent(2:end)],'b-');
xlim([-5,3]);
xticks(-5:2:3);
ylabel('z (km)');
xlabel('\delta^{13}C (o/oo)');
axis([-8 5 0 n*d*1e-3]);
ax(7)=gca;

axes(ha(8))
plot([LLCO3(1)  LLCO3(1)  LLCO3(2:n) ],[0 dm/1e3 zcent(2:end)],'r-');axis ij; hold on
plot([HLCO3(1)  HLCO3(1)  HLCO3(2:n) ],[0 dm/1e3 zcent(2:end)],'b-');
plot([LLCO3s(1) LLCO3s(1) LLCO3s(2:n)],[0 dm/1e3 zcent(2:end)],'r--');  %LL saturation
plot([HLCO3s(1) HLCO3s(1) HLCO3s(2:n)],[0 dm/1e3 zcent(2:end)],'b--');  %HL saturation
plot([LLCO2(1) LLCO2(1) LLCO2(2:n)],[0 dm/1e3 zcent(2:end)],'r:');      %LL CO2  
plot([HLCO2(1) HLCO2(1) HLCO2(2:n)],[0 dm/1e3 zcent(2:end)],'b:');      %HL CO2
ylabel('z (km)')
xlabel('CO_3 & CO_2 (mol m^{-3})')
axis([0 .3 0 n*d*1e-3])
ax(8)=gca;

if nto>=9
  axes(ha(9))
  plot([sLL(9,1,i) sLL(9,1,i) sLL(9,2:n,i)]*1000,[0 dm/1e3 zcent(2:end)],'r-');axis ij;hold on
  plot([sHL(9,1,i) sHL(9,1,i) sHL(9,2:n,i)]*1000,[0 dm/1e3 zcent(2:end)],'b-');
  ylabel('z (km)')
  xlabel('\delta^{18}O (o/oo)')
  axis([-0.8 2 0 n*d*1e-3])
  ax(9)=gca;
end

axes(ha(10))
plot([dwcLL(1,i) dwcLL(1,i) dwcLL(2:end,i)'],[0 dm/1e3 zcent(2:end)],'r-');axis ij; hold on
plot([dwcHL(1,i) dwcHL(1,i) dwcHL(2:end,i)'],[0 dm/1e3 zcent(2:end)],'b-');
ylabel('z (km)')
xlabel('CaCO3 dwf, low (r) and high (b)latitude')
axis([0 1 0 n*d*1e-3])
ax(10)=gca;

drawnow
pause(.1)

set(ax,'nextplot','replace')

end % Loop

% return
