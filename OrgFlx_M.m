
function [sr,rp,fnp,fcp,fCO2,fCO3,fHCO3] = OrgFlx_PETMIsoM_DN(TS,NPPI,LH,AT,Ao,Lf,NF,CO2,CO3,HCO3,GAM,GA,pCal,pOrg,Rcar,Rorg,l)
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
global nto d dm n dv sy rcp rcop rno rnp rcd rda rca lmdcar lmdorg lmdcorg 
global minO minN rDNno3 rDNnh4 rRBonh4 rRBnnh4 rSRh2s rSRnh4 rRBoh2s R13pdb 
global RCarbC13 swf LfHL rpm aLL aHL

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
RMcar  =  exp(-lmdcar* ( ( dm:d:n*d )-dm)); %  The flux profiles of CaCO3 based on the dissolved DIC;  
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
           
RMcorg  = Ao*NP*       ...
  [-RMcorg(1) (RMcorg(1:n-2)-RMcorg(2:n-1)).*GA(2:n-1)+...
               RMcorg(1:n-2).*(GA(1:n-2)-GA(2:n-1)).*pOrg(1:n-2) ...
               (RMcorg(n-1)-RMcorg(n)).*GA(n)+...
               RMcorg(n-1)*(GA(n-1)-GA(n))*pOrg(n-1)+...
               RMcorg(n)*GA(n)*pOrg(n)]; % [mol P/s]


% [mol Carbonate/s]

RMcar  = Ao*NP*rcp*rp*...
    [-RMcar(1) (RMcar(1:n-2)-RMcar(2:n-1)).*GA(2:n-1)+...
               RMcar(1:n-2).*(GA(1:n-2)-GA(2:n-1)).*pCal(1:n-2) ...
               (RMcar(n-1)-RMcar(n)).*GA(n)+...
               RMcar(n-1)*(GA(n-1)-GA(n))*pCal(n-1)+...
               RMcar(n)*GA(n)*pCal(n)];
          

f13dg = 1+((0.014*CO3/LH(4,1)-0.107)*LH(1,1)+10.53)/1000;
fCO2=(0.99869+4.9e-6*LH(1,1))/f13dg;
fCO3=(1.00722-52e-6*LH(1,1))/f13dg;
fHCO3=(1.01078-114e-6*LH(1,1))/f13dg;

CarbC13=LH(5,1)/LH(4,1)*fcp*fHCO3;
OrgC13=LH(5,1)/LH(4,1)*fnp*fCO2;
RcarbC13=(CarbC13/R13pdb-1)*1e3;
RorgC13=(OrgC13/R13pdb-1)*1e3;

%-------

% Oxidation of NH_4 and H_2S produces by denitrification and sulphate reduction is 
% expressed as a decay of the two compounds at oxygen rich levels. Decay scale set to 200 days 
%-------
  dc    = (1/(200*24*60^2));                                          % [s-1]
  Rnh4  = -LH(13,:).*dc*Ao.*GA.*dv.*(LH(8,:)>minO);                       % [mol NH4/s] 
  Rh2s  = -LH(14,:).*dc*Ao.*GA.*dv.*(LH(8,:)>minO);                       % [mol H2S/s]
% --------------------------------------------------------------------------

%-------
% Nitrogen fixation expressed as a nitrate source term in the LL mixed layer
%-------
nfo   = 1e6;            % [mol NO_3/s]                                          % [mol NO_3/s], st.val. 1e6

NFF=1; % Fijacion de Nitrogeno para ambas zonas
  nf    = nfo*NFF*(exp(NPP/NPN-1)-1)*(NPP>NPN);                      % [mol NO_3/s]


% Calculate/combine interior source terms
%-------

sr(3,:) = RMorg;                                                    % [mol P/s]
sr(3,1) = sr(3,1) + Rorg;                          % [mol P/s]

sr(4,:) = (RMcorg*rcp + RMcar*rcd);                                 % [mol D/s]
sr(4,1) = sr(4,1) + 2*Rcar;                          %  [mol C/s]

sr(5,:) = LH(5,1)/LH(4,1)* ...
        (RMcorg*rcp*fnp*fCO2 + RMcar*rcd*fcp);
sr(5,1) = sr(5,1) + (Rcar/(1+swf)*((R13pdb*(RCarbC13+1))...
          +(1+2*swf)*AT(5,1)/AT(4,1)));              % [molC13/s]
 sr(6,:) = LH(6,1)/LH(4,1)* ...
         (RMcorg*(1-2*(1-fnp))*(1-2*(1- fCO2))*rcp + ...
          RMcar*rcd*(1-2*(1-fcp)));                  % [molD14/s]
 sr(7,:) = RMcar*rca + RMorg*rda +...                               % [eq/s]
    RMcorg*rcp.*(LH(8,:)<=minO).*(LH(12,:)<=minN)+...
    RMorg*rDNnh4.*(LH(8,:)<=minO).*(LH(12,:)>minN)+... 
    RMorg*(rSRnh4-12).*(LH(8,:)<=minO).*(LH(12,:)<=minN) + Rnh4 + Rh2s;                              % [eq/s]
sr(7,1) = sr(7,1) +(2*Rcar + Rorg*rda);                                       % [eq/s]
sr(8,:) = ((RMorg*rno+RMcorg*rcop).*(LH(8,:)>minO)-rRBonh4*Rnh4   -rRBoh2s*Rh2s).*(LH(8,:)>minO);    % [mol O/s]

% The profile of N and P concentration is similar, thus RMorg based on P profile is equal to RMNorg based on P profile (shaffer et al., 1996) 
% The profile of DIC associate with remineralization is eauqal to that of DOC associate with remineralization 
% RMorg*rno+RMcorg*rcop: considering oxygen consume of natrate and DIC associated with remineralization, respectively

sr(8,1) = RMorg(1)*(rno+rcop)                    -rRBonh4*Rnh4(1)-rRBoh2s*Rh2s(1); % [mol O/s]       % The created O2 during the organic matter production in the surface water

sr(10,:) = 0;
sr(10,1) = 0;

% if nto>=12
  
  %sr(12,:)= RMorg*rnp + ...
  %          RMorg*rDNno3.*(LH(8,:)<=minO).*(LH(12,:)>minN)-rRBnnh4*Rnh4;   % [mol N0_3/s]
  
  %% trc 12,13,14,15 revised 20/07/2016
  sr(12,1)= RMorg(1)*rnp+nf+Rorg*rnp-rRBnnh4*Rnh4(1);
  
  sr(12,2:n)= RMorg(2:n)*rnp.*(LH(8,2:n)>minO) + ...
              RMorg(2:n)*rDNno3.*(LH(8,2:n)<=minO).*(LH(12,2:n)>minN)-...
              rRBnnh4*Rnh4(2:n);   % [mol N0_3/s]      
          
%   sr(12,1)= sr(12,1) + nf +Rorg*rnp;                                                 % [mol N0_3/s]

  sr(13,:)= RMorg*rDNnh4.*(LH(8,:)<=minO).*(LH(12,:)>minN) + ...
            RMorg*rSRnh4.*(LH(8,:)<=minO).*(LH(12,:)<=minN) + Rnh4;        % [mol NH_4/s]

  sr(14,1)=Rh2s(1);      
  sr(14,2:n)= RMorg(2:n)*rSRh2s.*(LH(8,2:n)<=minO).*(LH(12,2:n)<=minN) + Rh2s(2:n);        % [mol H_2S/s]
              % transformacion de sulfato a H2S


%   sr(15,1)   =-rRBnnh4*Rnh4(1).*facN2O(1)./(0.5*facN2O(1)+1); 
%   sr(15,2:n) = RMorg(2:n)*rnp.*facN2O(2:n)./(0.5*facN2O(2:n)+1)...
%               -rRBnnh4*Rnh4(2:n).*facN2O(2:n)./(0.5*facN2O(2:n)+1)...
%               -DN(2:n);
%           sr(15,1): % no hay consumpcion de N2O en capa superficial




% end
%-------
% rcp     =  106;                       % Carbon to phosphate, organic matter [-]
% rno     =  -32;                       % Oxygen to nitrogen, organic matter [-],    

%*****************************
% rcop    = -118;                       % which should be 106?
%*****************************

% rDNno3  = -84.8;                      % NO_3 to PO_4 in denitrification [] 
% rDNnh4  =  16;                        % NH_4 to PO_4 in denitrification [] 
% rRBonh4 = -2;                         % O_2 to NH_4 at redox boundaries [] 
% rRBnnh4 =  1;                         % NO_3 to NH_4 at redox boundaries [] 
% 
% rSRh2s  =  53;                        % H_2S to PO_4 in sulphate reduction [] 
% rSRnh4  =  16;                        % NH_4 to PO_4 in sulphate reduction [] 
% rRBoh2s = -2;       
return
