function [al]=LandExcnew_M(AT,LB)

% Calculates Land Biomass changes and land-sea gas exchanges of CO_2 and CH_4
% Input : 
%         AT    -  Atmospheric tracers, see ODE.m for the data structure.
%         LB    -  Land Biomasses 
%
% Output: 
%         al(1) - increment change in leafy biomass 12C (GtC/s)
%         al(2) - increment change in woody biomass 12C (GtC/s)
%         al(3) - increment change in litter biomass 12C (GtC/s)
%         al(4) - increment change in soil biomass 12C (GtC/s)
%         al(5) - increment change in leafy biomass 13C(GtC/s)
%         al(6) - increment change in woody biomass 13C (GtC/s)
%         al(7) - increment change in litter biomass 13C (GtC/s)
%         al(8) - increment change in soil biomass 13C (GtC/s)
%         al(9) - land air exchange of DIC [mol/s]  
%         al(10) - land air exchange of DI^13C [mol/s]
%         al(11) - land air exchange of CH4-C [mol/s]  
%         al(12) - land air exchange of CH4-^13C [mol/s]
%         al(13) - N20 production [mol/s]

% Activate global parameters
global sy fdiv Q10 CO2fer mgt rVa Q10met mdts
global Htem pCH4int pN2Oint pCO2int dc14 FacVeg


% Carbon ion fractionation factors for land air exchange
%-------
% C13

f13al   = 0.9819;
f13la   = 1;

% Carbon ion fractionation factor during land biosphere methane production

f13lm = 0.9700;
f14lm = 1-2*(1-f13lm);

% C14
f14al   = 1-2*(1-f13al);
f14la   = 1;
%-------

%-------------
% four box Land biosphere model derived from Siegenthaler and Oeschger (1987), includes CO2 fertilisation
% effect, and bacterial respiration in litter and soil as a functions of global mean temperature. Biomasses refer to one hemisphere
% revised according to the results of Gerber et al 2004

% revised values, 1/12/2006 75,375,200,700

Gro = FacVeg*50;                                                           %Pre-Anthropogenic (PA) leafy biomass, GtC
Woo = FacVeg*250;                                                          %PA woody biomass, GtC
Lio = FacVeg*60;                                                           %PA litter biomass, GtC
Slo = FacVeg*750;                                                          %PA soil biomass, GtC
LBtoto = Gro+Woo+Lio+Slo;

NPPLo = 30;                                                         %PA primary production on land,  GtC/yr
LBMPo = pCH4int/(rVa*mgt*mdts);                                     %PA land biosphere methane production, GtC/yr minus 1/2"PA",
                                                                    %anthropogenic input
N2OPo = pN2Oint/(rVa*150*sy);                                       %PA N2O production, mol/s 
 
ATem = AT(1,1)*sin(fdiv)+AT(1,2)*(1-sin(fdiv));                     %Mean atmospheric temperature
       
AGr = (35/60)*NPPLo/Gro;                                            %decay rate for leafy biomass, GtC/yr
AWo = (25/60)*NPPLo/Woo;                                            %decay rate for woody biomass, GtC/yr 
ALi = ((55/60)*NPPLo/Lio)*Q10^((ATem-Htem)/10);                     %decay rate for litter biomass with Q10 T dependence, GtC/yr  
ASl = ((15/60)*NPPLo-(LBMPo))/Slo*Q10^((ATem-Htem)/10);             %decay rate for soil biomass with Q10 T dependence, GtC/yr  

NPPL = NPPLo*...
          (1+CO2fer*log(AT(4,1)/pCO2int));                      
LBMP = LBMPo*LB(4,1)/Slo*Q10met^((ATem-Htem)/10);                   % land biosphere methane production with Q10 T dependence
                                                                    % land biosphere methane production with Q10 T dependence
N2OP = N2OPo*LB(4,1)/Slo*Q10^((ATem-Htem)/10);                      % N2O production 

al(1) = ((35/60)*NPPL - AGr*LB(1,1))/sy;                            % 12C rate of change for leafy biomass, GtC/s
al(2) = ((25/60)*NPPL - AWo*LB(2,1))/sy;                            % 12C rate of change for wood, Gt/s
al(3) = (AGr*LB(1,1) + (20/25)*AWo*LB(2,1) - ALi*LB(3,1))/sy;       % 12C rate of change for litter, GtC/s
al(4) = ((10/55)*ALi*LB(3,1) + (5/25)*AWo*LB(2,1)...
         - ASl*LB(4,1)-LBMP)/sy;                                    % 12C rate of change for soil, GtC/s

al(5) = ((35/60)*NPPL*AT(5,1)/AT(4,1)*f13al - AGr*LB(5,1))/sy;      % 13C rate of change for leafy biomass, GtC/s   
al(6) = ((25/60)*NPPL*AT(5,1)/AT(4,1)*f13al - AWo*LB(6,1))/sy;      % 13C rate of change for wood, GtC/s
al(7) = (AGr*LB(5,1) + (20/25)*AWo*LB(6,1) - ALi*LB(7,1))/sy;       % 13C rate of change for litter, GtC/s
al(8) = ((10/55)*ALi*LB(7,1) + (5/25)*AWo*LB(6,1)...
          - ASl*LB(8,1)-LBMP*LB(8,1)/LB(4,1)*f13lm)/sy;             % 13C rate of change for soil biomass, GtC/s 
  

al(9) = (-NPPL + (45/55)*ALi*LB(3,1) + ASl*LB(4,1))*mgt/sy;         %pCO2-12C sink/source to atmosphere from changes in LB, mol/s

al(10) = (-NPPL*AT(5,1)/AT(4,1)*f13al + (45/55)*ALi*LB(7,1) ...     %pCO2-13C sink/source to atmosphere from changes in LB, mol/s 
         + ASl*LB(8,1))*mgt/sy;
al(11) = LBMP*mgt/sy;                                               %CH4-12C source to atmosphere
  
al(12) = LBMP*LB(8,1)/LB(4,1)*f13lm *mgt/sy;                        %CH4-13C source to atmosphere

al(13) = N2OP;                                                      %N2O source to atmosphere 
al(14) = ((35/60)*NPPL*AT(6,1)/AT(4,1)*f14al - (AGr+dc14*sy)*LB(9,1))/sy;        % 14C rate of change for leafy biomass, GtC/s   
al(15) = ((25/60)*NPPL*AT(6,1)/AT(4,1)*f14al - (AWo+dc14*sy)*LB(10,1))/sy;       % 14C rate of change for wood, GtC/s
al(16) = (AGr*LB(9,1) + (20/25)*AWo*LB(10,1) - (ALi+dc14*sy)*LB(11,1))/sy;       % 14C rate of change for litter, GtC/s
al(17) = ((10/55)*ALi*LB(11,1) + (5/25)*AWo*LB(10,1)...
          - (ASl+dc14*sy)*LB(12,1)-LBMP*LB(12,1)/LB(4,1)*f14lm)/sy;             % 14C rate of change for soil biomass, GtC/s 
al(18) = (-NPPL*AT(6,1)/AT(4,1)*f14al + (45/55)*ALi*LB(11,1) ...    %pCO2-14C sink/source to atmosphere from changes in LB, mol/s 
         + ASl*LB(12,1))*mgt/sy;
al(19) = LBMP*LB(8,1)/LB(4,1)*f14lm *mgt/sy;                        %CH4-14 source to atmosphere

return

