function [as,P14C,CO2,CO3,HCO3,GAM]=GasExc_M(LH,AT,U,Ao)

% Calculates the air-sea gas exchange of CO_2 and O_2 for a single zone
% Input : LH    -  Low or High latitude vertical tracer distribution
%                  see ODE.m for the data structure.
%         AT    -  Atmospheric tracers, see ODE.m for the data structure.
%         U     - long term averaged wind speed [m2/s]
%         Ao    - ocean surface area [m2]
%
% Output: as(4) - air sea exchange of DIC [mol/s]  
%         as(5) - air sea exchange of DI^13C [mol/s]
%         as(6) - air sea exchange of DI^14C [mol/s]
%         as(8) - air sea exchange of O2 [mol/s] 
%         P14C  - cosmogenic production of ^14C [atoms/(sm^2)]
%         CO2   - surface ocean CO2 concentration  [mol/m3] 

% Activate global Parameters
global R13pdb R14oas kwo sy

%-------
% Iterate Carbonate system for pCO2 - partial pressure of surface CO2, 
% and K0 - solubility of CO2 

[pCO2w,K0,CO2,CO3,HCO3,GAM]=CarSys_M(LH(1,1),LH(2,1),LH(4,1),LH(7,1));
%-------
%LH(1,1)
%LH(2,1)
%LH(4,1)
%LH(7,1)
% Carbon ion fractionation factors for air-sea exchange
%-------
% C13
%LH(1,1)=25;
w13w    = 0.99912;
f13a    =  0.99869+4.9e-6*LH(1,1);
a13hco3 = (0.99869+4.9e-6*LH(1,1))/(1.01078-114e-6*LH(1,1));
a13co3  = (0.99869+4.9e-6*LH(1,1))/(1.00722-52e-6*LH(1,1));
f13T    = (CO2 + a13hco3*HCO3 + a13co3*CO3)/LH(4,1);
f13as   = w13w*f13a;
f13sa   = w13w*f13T;

w13CH4=0.9992;

% C14
w14w    = 1-2*(1-w13w);
f14a    = 1-2*(1-f13a);
a14hco3 = 1-2*(1-a13hco3);
a14co3  = 1-2*(1-a13co3);
f14T    = (CO2 + a14hco3*HCO3 + a14co3*CO3)/LH(4,1);
f14as   = f14a*w14w;
f14sa   = f14T*w14w;

%From Zhang et al, 1994
f13dg = 1+((0.014*CO3/LH(4,1)-0.107)*LH(1,1)+10.53)/1000;
f13sa = f13as/f13dg;
f14as = 1-2*(1-f13as);
f14sa = 1-2*(1-f13sa);

%-------

%-------------
% Gas exchange coefficient kw (m/s), formulated in the long-term averaged 
% winds U (m/s) (Wanninkhof 1992) and temperature dependent Schmidt number 
% for CO2 (Groger 2011).

Sc      =  1955.9*exp(-0.0663147*LH(1,1))+142.653;    %Groger
Kw      = .39/(100*60^2)*U^2*(Sc/660)^(-.5);

% Air sea C02 exchange
as(4) = Ao*Kw*K0*(                       AT(4,1) -                     pCO2w ); % [mol/s]
as(5) = Ao*Kw*K0*( f13as*AT(5,1)/AT(4,1)*AT(4,1) - f13sa*LH(5,1)/LH(4,1)*pCO2w ); % [mol/s]
as(6) = Ao*Kw*K0*( f14as*AT(6,1)/AT(4,1)*AT(4,1) - f14sa*LH(6,1)/LH(4,1)*pCO2w ); % [mol/s]
%-------------

%-------------
% Gas exchange coefficient kw (m/s), formulated in the long-term averaged 
% winds U (m/s) (Wanninkhof 1992) and temperature dependent Schmidt number 
% for O2 (Keeling et al 1998).
Sc   = 1638-81.83*LH(1,1)+1.483*LH(1,1)^2-0.008004*LH(1,1)^3; % T in oC
Kw      = .39/(100*60^2)*U^2*(Sc/660)^(-.5);

% Bunsen solubility coefficient for oxygen (Weiss 1970)
Tk   = LH(1,1)+273.15; % Kelvin temperature
beta = exp( -58.3877+85.8079*(100/Tk)+23.8439*log(Tk/100)+ ...
	    LH(2,1)*(-0.034892+0.015568*(Tk/100)-0.0019387*(Tk/100)^2) );
% Convertsion of Bunsen coefficient [atm-1] to solubility K0 [mol/m3/atm]. 
% The conversion ignores the the correction for non-ideality, ok for oxygen.
MVi  = 22.41361;         % Molar volume for an ideal gas [liter/mol]
K0   = beta/MVi*1e3;     % [mol/m3/atm]
 
% Air sea 02 exchange
as(8)= Ao*Kw*(K0*AT(8,1)-LH(8,1));  % [mol/s]
%-------------

% Methane
Scm      =  1771*exp(-0.06501*LH(1,1))+128.8;    %Groger
Kwm      = .39/(100*60^2)*U^2*(Scm/660)^(-.5);
% Bunsen solubility coefficient for methane (Weiss 1970)
Tk   = LH(1,1)+273.15; % Kelvin temperature
betam = exp(-68.8862+101.4956*(100/Tk)+28.7314*log(Tk/100)+ ...
	    LH(2,1)*(-0.076146+0.043970*(Tk/100)-0.0068672*(Tk/100)^2) );
% Convertsion of Bunsen coefficient [atm-1] to solubility K0 [mol/m3/atm]. 
% The conversion ignores the the correction for non-ideality, ok for oxygen.
MVi  = 22.41361;         % Molar volume for an ideal gas [liter/mol]
K0m   = betam/MVi*1e3;     % [mol/m3/atm]
% as(2)= Ao*Kwm*(K0m*AT(2,1)-LH(10,1));  % [mol/s] 
as(2)= Ao*Kwm*(K0m*AT(2,1)-LH(10,1));  % [mol/s] 
as(7)=Ao*Kwm*w13CH4*(K0m*AT(7,1)-LH(11,1));

%-------------
% Calculate delta^13C and Delta^14C values in permil for the atmosphere 
% and set ^14C production to maintain Delta^14C approx = 0;
d13a     = (AT(5,1)./AT(4,1)/R13pdb-1)*1e3;
D14a     = (AT(6,1)./AT(4,1)./R14oas.*(R13pdb*.975./(AT(5,1)./AT(4,1))).^2-1)*1e3;
P14C     = max(0,1e5*(0-D14a));   % Production [atoms/(sm^2)]
%-------------

return
