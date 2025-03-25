% function Thilda_M

% Integrates the ordinary differential equations in ODE.m using 
% a fourth order Runge Kutta scheme with fixed timestep (h).
tic
clear all

% Read in initial conditions
%--------------------------------------------
[LL,HL,AT,LB,LLCO3,HLCO3,LLCO3s,HLCO3s,pCalLL,pCalHL,pOrgLL,pOrgHL,dwpCalLL,dwpCalHL,...
 fimLL,fimHL,dwpCalssLL,dwpCalssHL,dwpOrgLL,dwpOrgssLL,dwpOrgHL,dwpOrgssHL,wsedLL,wsedHL ]=InitCnd_M;

global sy n lmdcorg lmdorg dm d rcp lmdcar nto 
global aLL aHL LfLL LfHL ULL UHL GAw GAc NCFLL NCFHL CAFLL CAFHL TSL TSH NPPIL NPPIH 
%% set initial parameter and inject carbon

tend = 100000;   % End-time of integration [yr]
input.tend = tend;
dtout  = 100;     % Output interval [yr]
h      = 1/25;    % Timestep [yr], st val 1/25
ko     = 1/h;     % Multiple of timestep at which CO3 and time dependent sediment model are calculated
dt    = h*sy;    % Timestep [sec]
t     = 0;

dataf = readtable("data\Forcing.xlsx","Sheet","Tune1");
input.time = dataf.time;
input.meth.rate  =dataf.meth/2;
input.mett.rate = dataf.mett/2;
input.org.rate = dataf.org/2;
input.lip.rate = dataf.lip/2;

input.meth.total = input_total(input.time,input.meth.rate);
input.mett.total = input_total(input.time,input.mett.rate);
input.org.total = input_total(input.time,input.org.rate);
input.lip.total = input_total(input.time,input.lip.rate);

input.meth.mfrac  = 0;
input.mett.mfrac  = 1;

t_serial = 0:0.1:100; % unit kyr
for i = 1: length(t_serial)
input.meth.flux(i) = enforcing(t_serial(i), input.time, input.meth.rate);
input.mett.flux(i) = enforcing(t_serial(i), input.time, input.mett.rate);
input.org.flux(i) = enforcing(t_serial(i), input.time, input.org.rate);
input.lip.flux(i) = enforcing(t_serial(i), input.time, input.lip.rate);
end
%% 
figure("Position",[0,0,500,450])
[ha,~] = tight_subplot(2,1,[0.15,0],[0.15,0.02],[0.1,0.02]);
axes(ha(1))
plot(t_serial,input.meth.flux*2,"LineWidth",1); hold on
plot(t_serial,input.mett.flux*2,"LineWidth",1);
plot(t_serial,input.org.flux*2,"LineWidth",1);
plot(t_serial,input.lip.flux*2,"LineWidth",1);
ylim([0,1]);
ylimits = ylim;
line([22.3,22.3],[ylimits(1),ylimits(2)],"LineStyle","--","Color","k","LineWidth",1);
line([18.5,18.5],[ylimits(1),ylimits(2)],"LineStyle","--","Color","k","LineWidth",1);
xlim([0,max(t_serial)]);
xticks(0:5:max(t_serial));
ticklabels("x",xticks,2,"F");
xtickangle(0);
text_norm(0.8,0.9,"Meth:"+ input.meth.total*2,10);
text_norm(0.8,0.80,"Mett:"+ input.mett.total*2,10);
xlabel('Time (kyr)');
ylabel("Carbon input (Gtyr^{-1})");
box off
BorderLine(gca,1);
set(gca,"linewidth", 1,"FontSize",10,"FontName", "Times","TickLength",[0.01,0.015],"Layer","top");

axes(ha(2))
data1 = readtable("C:\Users\19328\Nutstore\1\我的坚果云\PETM\Data\Xiong-2023-PETM_1.xlsx","Sheet","C and O isotope");
scatter(data1.age_westerhold(1:end-3)+22.3,data1.dC(1:end-3)-mean(data1.dC(1:60)),7,"Marker","o","MarkerEdgeColor","#404040","MarkerFaceColor","#404040","LineWidth",0.01);
xlim([0,100]);
xticks(0:5:100);
ticklabels("x",xticks,2,"F");
xtickangle(0);
xlabel('Time (kyr)');
ylim([-7.5,3]);
ylimits = ylim;
yticks(-7:1:3);
line([22.3,22.3],[ylimits(1),ylimits(2)],"LineStyle","--","Color","k","LineWidth",1);
line([18.5,18.5],[ylimits(1),ylimits(2)],"LineStyle","--","Color","k","LineWidth",1);
ticklabels("y",yticks,2,"F");
ylabel('\Delta\delta^{13}C (‰)');
limY = get(gca, 'Ylim');
box off
BorderLine(gca,1);
set(gca,"linewidth", 1,"FontSize",10,"FontName", "Times","TickLength",[0.01,0.015],"Layer","top");
print(gcf,"Figure\Input","-dpng","-r600")
%% 
%--------------------------------------------
% Integration
%------------------------------


for c =1:tend/dtout
  for cc=0:dt:(dtout*sy-dt)
%%---------------
      
    if mod(t/sy,ko*h)==0
      %---------------
      % Sediment model
      %---------------
    
      % Iterate the carbonate system for the vertical profile of 
      % the CO3 concentration and saturation state 
      
      % org_vel = [0,0];

      [LLCO3,LLCO3s,LLCO2,LLHCO3,LLHp]=CarSysPres_M(LL);
      [HLCO3,HLCO3s,HLCO2,HLHCO3,HLHp]=CarSysPres_M(HL);
      
      [QEBLL,QEBHL,QLL,QHL,aHLNI,Fw,Td] = AtmEnerBal_M(LL(1,1),HL(1,1),AT);
      [RcarLL,RcarHL,RorgLL,RorgHL,Wcarb,Wsil,Worg]=ExtForce_M(AT);
      % Get surface CO2
      [asLL,P14C,CO2LL,CO3LL,HCO3LL,GAMLL] = GasExc_M(LL,AT,ULL,aLL);
      [asHL,P14C,CO2HL,CO3HL,HCO3HL,GAMHL] = GasExc_M(HL,AT,UHL,aHLNI);
      
      % Get New Production
      
      [srcLL,rpl,fnpl,fcpl,fCO2l,fCO3l,fHCO3l] = OrgFlx_M(TSL,NPPIL,LL,AT,aLL  ,LfLL,1,CO2LL,CO3LL,HCO3LL,GAMLL,GAw,pCalLL,pOrgLL,RcarLL,RorgLL,1);
      [srcHL,rph,fnph,fcph,fCO2h,fCO3h,fHCO3h] = OrgFlx_M(TSH,NPPIH,HL,AT,aHLNI,LfHL,0,CO2HL,CO3HL,HCO3HL,GAMHL,GAc,pCalHL,pOrgHL,RcarHL,RorgHL,0);
               
      NPL = (srcLL(3,1)-RorgLL)/aLL;
      NPHM = (srcHL(3,1)-RorgHL)/aHL;
            
      %Time dependent sediment model,
      
      [pCalLL,pOrgLL,dwpCalLL,dwpOrgLL,fimLL,wsedLL]=...
         SMtdCorgNew2_M(TSL,LL,rpl,NCFLL,CAFLL,LLCO3,LLCO3s,NPL,dwpCalLL,dwpCalssLL,dwpOrgLL,dwpOrgssLL,fimLL,wsedLL,ko,h);
      [pCalHL,pOrgHL,dwpCalHL,dwpOrgHL,fimHL,wsedHL]=...
         SMtdCorgNew2_M(TSH,HL,rph,NCFHL,CAFHL,HLCO3,HLCO3s,NPHM,dwpCalHL,dwpCalssHL,dwpOrgHL,dwpOrgssHL,fimHL,wsedHL,ko,h);
     
    
      RMorg  =  exp(-lmdorg* ( ( dm:d:n ...
          *d )-dm)); 
      RMcorg =  exp(-lmdcorg*( ( dm:d:n*d )-dm));
      RMcar  =  exp(-lmdcar* ( ( dm:d:n*d )-dm));
     
                 
      RMorgLL  = aLL*NPL*       ...
               [-RMorg(1) (RMorg(1:n-2)-RMorg(2:n-1)).*GAw(2:n-1)+...
               RMorg(1:n-2).*(GAw(1:n-2)-GAw(2:n-1)).*pOrgLL(1:n-2) ...
               (RMorg(n-1)-RMorg(n)).*GAw(n)+...
               RMorg(n-1)*(GAw(n-1)-GAw(n))*pOrgLL(n-1)+...
               RMorg(n)*GAw(n)*pOrgLL(n)];
           
      RMorgHL  = aHL*NPHM*       ...
               [-RMorg(1) (RMorg(1:n-2)-RMorg(2:n-1)).*GAc(2:n-1)+...
               RMorg(1:n-2).*(GAc(1:n-2)-GAc(2:n-1)).*pOrgHL(1:n-2) ...
               (RMorg(n-1)-RMorg(n)).*GAc(n)+...
               RMorg(n-1)*(GAc(n-1)-GAc(n))*pOrgHL(n-1)+...
               RMorg(n)*GAc(n)*pOrgHL(n)];
           
      RMcorgLL  = aLL*NPL*rcp*       ...
               [-RMcorg(1) (RMcorg(1:n-2)-RMcorg(2:n-1)).*GAw(2:n-1)+...
               RMcorg(1:n-2).*(GAw(1:n-2)-GAw(2:n-1)).*pOrgLL(1:n-2) ...
               (RMcorg(n-1)-RMcorg(n)).*GAw(n)+...
               RMcorg(n-1)*(GAw(n-1)-GAw(n))*pOrgLL(n-1)+...
               RMcorg(n)*GAw(n)*pOrgLL(n)];
           
      RMcorgHL  = aHL*NPHM*rcp*       ...
               [-RMcorg(1) (RMcorg(1:n-2)-RMcorg(2:n-1)).*GAc(2:n-1)+...
               RMcorg(1:n-2).*(GAc(1:n-2)-GAc(2:n-1)).*pOrgHL(1:n-2) ...
               (RMcorg(n-1)-RMcorg(n)).*GAc(n)+...
               RMcorg(n-1)*(GAc(n-1)-GAc(n))*pOrgHL(n-1)+...
               RMcorg(n)*GAc(n)*pOrgHL(n)];
           
      RMcarLL  = aLL*NPL*rcp*rpl*...
               [-RMcar(1) (RMcar(1:n-2)-RMcar(2:n-1)).*GAw(2:n-1)+...
               RMcar(1:n-2).*(GAw(1:n-2)-GAw(2:n-1)).*pCalLL(1:n-2) ...
               (RMcar(n-1)-RMcar(n)).*GAw(n)+...
               RMcar(n-1)*(GAw(n-1)-GAw(n))*pCalLL(n-1)+...
               RMcar(n)*GAw(n)*pCalLL(n)];
           
      RMcarHL  = aHL*NPHM*rcp*rph*...
               [-RMcar(1) (RMcar(1:n-2)-RMcar(2:n-1)).*GAc(2:n-1)+...
               RMcar(1:n-2).*(GAc(1:n-2)-GAc(2:n-1)).*pCalHL(1:n-2) ...
               (RMcar(n-1)-RMcar(n)).*GAc(n)+...
               RMcar(n-1)*(GAc(n-1)-GAc(n))*pCalHL(n-1)+...
               RMcar(n)*GAc(n)*pCalHL(n)];
                      
      BurorgLL= sum(RMorgLL);
      BurorgHL= sum(RMorgHL);        
      BurcorgLL= sum(RMcorgLL);
      BurcorgHL= sum(RMcorgHL);
      BurcarLL= sum(RMcarLL);
      BurcarHL= sum(RMcarHL); 
           
    end  
   
    t=t+dt;

    [k1LL,k1HL,k1AT,k1LB,srLL1,srHL1,car1] =...
        ODE_M(t,LL,HL,input,AT,LB,pCalLL,pCalHL,pOrgLL,pOrgHL);
      nLL = LL+k1LL*dt/2; 
      nHL = HL+k1HL*dt/2;
      nAT = AT+k1AT*dt/2;
      nLB = LB+k1LB*dt/2;
    [k2LL,k2HL,k2AT,k2LB,srLL2,srHL2,car2] =...
        ODE_M(t+dt/2,nLL,nHL,input,nAT,nLB,pCalLL,pCalHL,pOrgLL,pOrgHL);
      nLL = LL+k2LL*dt/2; 
      nHL = HL+k2HL*dt/2; 
      nAT = AT+k2AT*dt/2;
      nLB = LB+k2LB*dt/2;
    [k3LL,k3HL,k3AT,k3LB,srLL3,srHL3,car3] =...
        ODE_M(t+dt/2,nLL,nHL,input,nAT,nLB,pCalLL,pCalHL,pOrgLL,pOrgHL);
      nLL = LL+k3LL*dt; 
      nHL = HL+k3HL*dt; 
      nAT = AT+k3AT*dt; 
      nLB = LB+k3LB*dt;
    [k4LL,k4HL,k4AT,k4LB,srLL4,srHL4,car4] =...
        ODE_M(t+dt,nLL,nHL,input,nAT,nLB,pCalLL,pCalHL,pOrgLL,pOrgHL);
 

    LL = LL+1/6*dt*( k1LL+2*k2LL+2*k3LL+k4LL );
    HL = HL+1/6*dt*( k1HL+2*k2HL+2*k3HL+k4HL );
    AT = AT+1/6*dt*( k1AT+2*k2AT+2*k3AT+k4AT );
    LB = LB+1/6*dt*( k1LB+2*k2LB+2*k3LB+k4LB );
    
   end % counter cc
  car = 1/6*(car1+2*car2+2*car3+car4);
  meth_vel_12C = car(1)
  mett_vel_12C = car(2)

  % sorgLL_input = RorgLL_input;
  scar(:,c) = car; 
  st(  c    )   = t/sy;
  sAT( :,:,c)   = AT;  
  sLL( :,:,c)   = LL;
  sHL( :,:,c)   = HL;
  sMT( 1,:,c)   = .84*sLL(1,:,c)+ .16*sHL(1,:,c);

  sLB( :,:,c)   = LB; 
  srLL(:,:,c)   = 1/6*(srLL1+2*srLL2+2*srLL3+srLL4);
  srHL(:,:,c)   = 1/6*(srHL1+2*srHL2+2*srHL3+srHL4);
  
  dwcLL(:,c)    = dwpCalLL;
  dwcHL(:,c)    = dwpCalHL;
  dworgCLL(:,c) = dwpOrgLL;
  dworgCHL(:,c) = dwpOrgHL;
  fiLL(:,c)     = fimLL;
  fiHL(:,c)     = fimHL;
  pCalCLL(:,c)  = pCalLL;
  pCalCHL(:,c)  = pCalHL;
  pOrgCLL(:,c)  = pOrgLL;
  pOrgCHL(:,c)  = pOrgHL;
  wsLL(:,c)     = wsedLL;
  wsHL(:,c)     = wsedHL;
  
  QeLL(:,c)     = QEBLL;
  QeHL(:,c)     = QEBHL;
  ASLL(:,c)     = asLL;
  ASHL(:,c)     = asHL;
  
  RoLL(1,1,c)   = RorgLL;
  RcLL(1,1,c)   = RcarLL;
  RoHL(1,1,c)   = RorgHL;
  RcHL(1,1,c)   = RcarHL;
  
  BurcLL(:,c)   = BurcorgLL;
  BurcHL(:,c)   = BurcorgHL;
  BurcaLL(:,c)  = BurcarLL;
  BurcaHL(:,c)  = BurcarHL;
  BurpLL(:,c)   = BurorgLL;
  BurpHL(:,c)   = BurorgHL;
 

  disp(strcat(['Integration at year ' num2str(c*dtout) ...      
		   ' (tend = ',num2str(tend) ')']))
 
 
end % counter c
%------------------------------
%
toc
%% 

save OutThilda_M st sLL sHL sMT sAT sLB srLL srHL dwcLL dwcHL dworgCLL dworgCHL...
    fiLL fiHL pCalCLL pCalCHL pOrgCLL pOrgCHL wsLL wsHL BurcLL BurcHL BurcaLL ...
    BurcaHL BurpLL BurpHL RoLL RoHL RcLL RcHL QeLL QeHL ASLL ASHL input scar
%% 

save("dworkspace.mat")

msgbox('Program has finished running!', 'Notification');
beep;

%------------------------------
%
% Copyright � 2017 Danish Center for Earth System Science
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software 
% and associated documentation files (the "Software"), to deal in the Software without restriction, 
% including without limitation the rights to use, copy, and modify copies of the Software, and 
% to permit persons to whom the Software is furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all copies or 
% substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING 
% BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

