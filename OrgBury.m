function [sr,rp,fnp,fcp,fCO2,fCO3,fHCO3] = OrgFlx_PETMIsoM_DN(LH,Ao,Lf,CO2,GAM,GA,pOrg)
% Calculates tracer sources for one sector due to 
% POM and CaCO_3 production in the euphotic zone, 
% remineralization and dissolution at depth, 
% nitrogen fixation, nitrate reduction and sulphate 
% reduction. Sinks are defined negative.   
% 
% Input : LH  -  Low or High latitude vertical tracer distribution
%                see ODE.m for the data structure.
%         AT -Atmospheric tracer concentrations
%         Ao  - Ocean surface area [m2]
%         Lf  - Limitation factor for new production [-]
%               (eg. light/iron limitation, unity means no limitation)
%         NF  - Nitrogen fixation on/off = 1/0
%         CO2 - surface ocean CO2 concentration [mol/m3] 
%         GA  - Ocean area at depth relative to surface
%         pCal- Percentage of Calcite dissolution in the sediment
%        
% Output: sr -  Vertical distribution of total tracer sources (object? ocean?)
%               sr(3,:)   - phosphate [mol P/s]   
%               sr(4,:)   - DIC [mol DIC/s]   
%               sr(7,:)   - alkalinity [eq/s]   
%               sr(8,:)   - oxygen [mol O_2/s]   
%               sr(10,:)  - nitrogen [mol NO_3/s]   
%               sr(11,:)  - ammonium [mol NH_4/s]
%               sr(12,:)  - hydrogen sulphide [mol H_2S/s]

% Activate global paramter values
%ParVal_PETMIsoM;
global  d dm n sy rnp  lmdcar lmdorg lmdcorg LfHL rpm

%-------

% Nondimensional (No unit) flux profiles of carbonate and organic POM (Production)
% at the lower level of model layers

% 1. New production form with CaCo3 formation (Coral, Coccolithophore);
% 2. Remineralization of new production cause the decreased oganic matter but increased phosphate and dissolved organic carbon. 
% Notely, The profile of phosphate and dissolved are natural exponential distribution, but increased rate is difference, which can expressed by lmdorg and lmdcorg.
% 3. lmdorg and lmdcorg is the inverse of e-fold depth for phosphate and dissolved organic carbon
% 4. The profile of CaCO3 and new production is contrast with that of DIC,
% phoshate and dissoleved organic matter.


%-------
RMorg  =  exp(-lmdorg* ( ( dm:d:n*d )-dm)); %  The flux profiles of new production based on phosphate verticle distribution;   
RMcorg =  exp(-lmdcorg*( ( dm:d:n*d )-dm)); %  The flux profiles of new production based on dissolved organic carbon verticle distribution;

%-------

% New production of organic matter (euphotic zone phosphate uptake)
%-------
bpe   = 1/sy;                                       % Bio-Production Efficiency [s-1]
Phlf  = 1.6e-5;                                     % Half-saturation constant [molP/m3] 
NPP   = bpe*Lf*dm*LH(3,1)^2/(Phlf+LH(3,1));         % [mol P/m2 s] 
Nhlf  = 16*Phlf;                                    % Half-saturation constant [molP/m3], st. val. 1e-6 
NPN   = 1/rnp*bpe*Lf*dm*LH(12,1)^2/(Nhlf+LH(12,1)); % [mol P/m2 s]
NP    = min(NPP,NPN);                               % [mol P/m2 s] % New production (organic matter) 
                        
% Production ratio 'rp' ("Rain ratio")
%-------

Tr  = 10;                                           % [oC] Reference temperatur, st. val 10
p1  = 1;                                            % [-] O.Marchal 98,Maier-Reimer 93
p2  = 0.1565;                                        % [oC^-1] Maier-Reimer 93
                                                    %recalibrated   22/8/2014
if GAM<1
    rp=0;
else
rp  = rpm* (p1*exp(p2*(LH(1,1)-Tr))/(1+p1*exp(p2*(LH(1,1)-Tr))))*((GAM-1)/(1+(GAM-1)));
end

%-------
% Carbon ion fractionation for organic and carbonate new production
%-------

fcp=1;

if LH(3,1)>0.5e-3 %0.2e-3
fnp    = ((-(25-(116.96*1000*LfHL*LH(3,1)+81.42)/(1000*CO2)))/1000+1);   
else
fnp    = ((-(25-(116.96*1000*LH(3,1)+81.42)/(1000*CO2)))/1000+1);
end

%-------

% RMorg: The remineralization flux of new production in different layers.
% Notably, the first layer (0-100 m) is formation of new production (negative), not considering remineralization;
% The other layers is remineralization of new production (Positive), not consider formation
% RMorg are compoposed of remineralizaion overlying the water and seafloor.
% The partial organic matter input into seafloor in each layer (such as slope)
RMorg  = Ao*NP*...
     [-RMorg(1) (RMorg(1:n-2)-RMorg(2:n-1)).*GA(2:n-1)+...
               RMorg(1:n-2).*(GA(1:n-2)-GA(2:n-1)).*pOrg(1:n-2) ...
               (RMorg(n-1)-RMorg(n)).*GA(n)+...
               RMorg(n-1)*(GA(n-1)-GA(n))*pOrg(n-1)+...
               RMorg(n)*GA(n)*pOrg(n)];   % [mol P/s] 
end