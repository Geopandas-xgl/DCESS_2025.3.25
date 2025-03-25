function [CO3,CO3s,CO2,HCO3,Hp]=CarSysPres_M(LH)

% Calculates solubilities and iterate the pressure dependent carbonate system 
% Input : LH  - Low or High latitude vertical tracer distribution
%               see ODE.m for the data structure.
%               The sub-set used in CarSysPres are    
%                     1:   - temperature [oC]
%                     2:   - salinity []
%                     4:   - DIC [mol/m3] 
%                     7:   - alkalinity [eq/m3]
% Output variables
%      K0         [mol/m3/atm]
%      CO2        [mol/m3]
%      CO3        [mol/m3]
%      HCO3       [mol/m3]
%
%

warning off

% Activate global Parameters

global d dm n Cao So

% Convert input DIC and ALK to [mol/kg] from [mol/m3] using a fixed density
Ro   = 1028;    % [kg/m3]  
DIC  = LH(4,:)/Ro;  % [mol/kg]
ALK  = LH(7,:)/Ro;  % [mol/kg]

% Convert temperature to kelvin scale
t    = LH(1,:)  ; % [oC]
T    = t +273.15; % [K]


% Define salinity for convenience
S    = LH(2,:);

% Solubility of CO2 [mol/kg/atm] given by Weiss 1974, see Millero 1995 eq. 26
K0   =  exp(- 60.2409 + 93.4517*(100./T) + 23.3585*log(T./100) + ...
            S.*( 0.023517  - 0.023656*(T./100) + 0.0047036*(T./100).^2));

% Thermodynamic constants for the dissociation of carbonic acid, 
% K1 and K2 [mol/kg] (Mehrbach et al. 1973), see Millero 95 eq 35/36
% The Mehrbach constants are recomended by Lee and Millero 1997 
% K1  = 10.^(-(3670.7./T-62.008+9.7944*log(T)-0.0118*S+0.000116*S.^2));
K1mod   = 10.^(-(3670.7./T-62.008+9.7944*log(T)-0.0118*S+0.000116*S.^2));
K1     =   (1+0.155*((29.9-52.82)/52.82))*K1mod;        %revised from Panchuk 2008

% K2   = 10.^(-(1394.7./T+ 4.777              -0.0184*S+0.000118*S.^2));
K2mod   = 10.^(-(1394.7./T+ 4.777              -0.0184*S+0.000118*S.^2));
K2     =   (1+0.442*((29.9-52.82)/52.82))*K2mod;


% Dissociation constant for sea-water [mol/kg], Millero 95 (eq 63) 
Kw   = exp(148.9802-13847.26./T-23.6521*log(T) + ...
           (-5.977+118.67./T+1.0495*log(T)).*S.^.5 - ...
           0.01615*S);

% Dissociation constant for boric acid [mol/kg], Dickson 1990, recommended by Millero 95 (eq 52) 
Kb   = exp((-8966.90-2890.51*S.^.5-77.942*S+1.726*S.^1.5-0.0993*S.^2)./T + ...
           148.0248+137.194*S.^.5+1.62247*S + ...
           (-24.4344-25.085*S.^.5-0.2474*S).*log(T) + ...
           0.053105*S.^.5.*T);   

% Dissociation constant for calcite [mol/kg], Mucci 1983
Kc   = 10.^(-171.9065-0.077993*T+2839.319./T+71.595*log10(T)+...
       (-0.77712+0.0028426*T+178.34./T).*S.^0.5-0.07711*S+0.0041249*S.^1.5);

% Total boric acid in sea-water, [mol/kg] Millero 1982, see Millero 1995 below eq 53
BR   = 0.000416*S/So;

%---------------------------------------------
% Pressure dependence of system thermodynamics
R = 83.131;                         %[bar cm3 mol-1 K-1] Gas constant
P = 1+([dm/2 dm+(1:n-1)*d-d/2])/10;   %[bar] Pressure

% Molar volume and compressibility changes of dissociation reactions
%Kw: H2O
MolarVol  =  (-2.00E01) + ( 1.120E-01)*t + (-1.410E-03)*t.^2; 
DissCompr =  (-5.13E00) + ( 7.940E-02)*t                     ; 
ratio     = -(MolarVol./(R*T)).*P+(.5E-03*DissCompr./(R*T)).*P.^2;
Kw_p      =  Kw.*exp(ratio);

%K1: H2CO3
MolarVol  =  (-2.550E+01) + ( 1.271E-01)*t + (    0E-03)*t.^2; 
DissCompr =  (-3.080E+00) + ( 8.770E-02)*t                   ; 
ratio     = -(MolarVol./(R*T)).*P+(.5E-03*DissCompr./(R*T)).*P.^2;
K1_p      =  K1.*exp(ratio);

%K2: HCO3(-)
MolarVol  =  (-1.582E+01) + (-2.190E-02)*t + (    0E-03)*t.^2 ; 
DissCompr =  ( 1.130E+00) + (-1.475E-01)*t                    ; 
ratio     = -(MolarVol./(R*T)).*P+(.5E-03*DissCompr./(R*T)).*P.^2;
K2_p      =  K2.*exp(ratio);

%Kb: B(OH)3
MolarVol  =  (-2.948E+01) + ( 1.622E-01)*t + ( 2.608E-03)*t.^2; 
DissCompr =  (-2.840E+00) + (     0E-02)*t                    ; 
ratio     = -(MolarVol./(R*T)).*P+(.5E-03*DissCompr./(R*T)).*P.^2;
Kb_p      =  Kb.*exp(ratio);

%Kc: CaCO3(cal)
MolarVol  =  (-4.876E+01) + ( 5.304E-01)*t + (     0E-03)*t.^2; 
DissCompr =  (-1.176E+01) + ( 3.692E-01)*t                    ; 
ratio     = -(MolarVol./(R*T)).*P+(.5E-03*DissCompr./(R*T)).*P.^2;
%Kc_p      =  Kc.*exp(ratio);
Kc_pmod      =  Kc.*exp(ratio);
Kc_p = Kc_pmod - 3.655e-8*(52.82/10.25-29.9/18.2);  %revised from Panchuk 2008

% Use a salinity corrected value for the concentration of Ca(2+) 
% [Millero 1982] to  estimate the saturation concentration of CO3(2-) 
CO3s = Kc_p./(Cao*S/So);

%---------------------------------------------

% Iterate the carbon dioxide in solution CO_2 using the recursive formulation 
% of Antonie and Morel 1995. The values for the hydrogen ion activity, H+, 
% and the carbonate alkalinity are also computed, CA = HCO3- + 2 CO3--
Kr   = K1_p./K2_p;
A    = ALK;
imx  = 500;
for i=1:imx
Z    = ( (DIC.*Kr).^2+A.*Kr.*(2*DIC-A).*(4-Kr) ).^.5;
CO2  = DIC-A+(A.*Kr-DIC.*Kr-4*A+Z)./(2*(Kr-4));
Hp   = CO2.*K1_p./(2*A) + (( ((CO2.*K1_p).^2+8*A.*CO2.*K1_p.*K2_p) ).^.5)./(2*A);
CA   = ALK-BR.*Kb_p./(Kb_p+Hp)-Kw_p./Hp+Hp;
if max(abs(CA-A))<1e-8 break; end
A    = CA;
end
if i==imx; disp('!! Carbonate system not converging, CarSysPres !!'); end

CO3  = K1_p.*K2_p.*(  DIC./( (1+K1_p./Hp+K1_p.*K2_p./Hp.^2) )  )./Hp.^2;  % [mol/kg]

% Convert output values to model units

CO3  = CO3 *Ro;                       % [mol/m3] 
CO3s = CO3s*Ro;                       % [mol/m3] 
CO2  = CO2 *Ro;                       % [mol/m3] 
HCO3  = DIC*Ro-CO3-CO2;
return

