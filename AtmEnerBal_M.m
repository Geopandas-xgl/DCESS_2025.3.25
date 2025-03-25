function [QEBLL,QEBHL,QLL,QHL,aHLNI,Fw,PTad] = AtmEnerBal_M(ToLL,ToHL,AT)

% Redefine
Ta=AT(1,:); 
pCO2=AT(4,1);
pCH4=AT(2,1);
pN2O=AT(3,1);

% Activate global parameters
global olfHL aH aLL fdiv fland pCO2int pCH4int pN2Oint
global Tsnow SCred rc dm albch

% Local parameter values
%----------------
% Albedo, function of latitude (Hartmann 1994)
al0   =  0.7*(1-albch);
al2   = -0.175;
ali   =  .62; % Albedo of ice covered area,  st. val.  0.62

% Solar radiation, function of latitude, Nakamura, J Clim 1994, st. val. 1365
Qo    =  1365*(1-SCred);                 %[W/m2]
Q2    = -0.482;   

% Zero degree backgrod radiation for pCO2 = 280,
% pCH4 = 700, pN2O = 270; 

Ao=241;                                  %St val 241.0 for B=1  

ACO2=5.32*log(pCO2/pCO2int)+ 0.39*log(pCO2/pCO2int)^2; 

% Constants for new CH4 forcing fit at high CH4
a = -4.054; b = -6482; c = 9.542; d = -60.99;   

pCH4int2=pCH4int-0.005e-6;
M1=2.5*1e-6;
if pCH4<2.5e-6
  ACH4 = 1173*(pCH4^(.5)-pCH4int2^(.5))-71636*(pCH4^(.5)-pCH4int2^(.5))^2;
elseif pCH4>100e-6
  ACH4= a*exp(b*pCH4) + c*exp(d*pCH4)-0.867;
else
 ACH4=0.824+0.8*log(pCH4/M1)+0.2*log(pCH4/M1)^2;
end

 AN2O=3899*(pN2O^(.5)- pN2Oint^(.5))+38256*(pN2O^(.5)- pN2Oint^(.5))^2;
 
 AOverlap= -16.16*exp(-0.036*(log(pCO2-pCO2int)-0.0024)^2 ...
    -0.05*(log(pN2O-pN2Oint)+6.5)^2)-24*exp(-0.02*(log(pCH4-pCH4int2)-0.01)^2 ...
    -0.044*(log(pN2O-pN2Oint)+7.73)^2);

% Zero degree background radiation for variable pCO2, pCH4 and pN2O 
Atot  =  Ao - ACO2 - ACH4 - AN2O - AOverlap;        % [W/m2]


% Sensitivity of LW radiation to temp.

B=1;                      % st val 1.0 for Ao =241
  
% Atmospheric transports
Lv    =  2.25e6*1e3;  % Latent heat of vaporization [J/1000kg] = [J/m^3]   
Cfw   =  1.1*10e9;    % Atmospheric exchange coefficient for water, st. val 10e9    
Cs    =  1.237*2.5e11;% Atmospheric exchange coefficient for sensible heat, st.val2.5e11     
N     =  2.5;         % Gradient exponent

% Haney coefficient, st. val 30
lambda= 30;           % [W/m2/oC] 

% Direct solar heating of ocean surface, LL only 
Dheat = 30;           % [W/m2] 
%----------------

% Establish the atmospheric temperature profile, 2.order legendre pol. in 
% sine of lat, Ta(lat) = PTa(1) + .5*PTa(2) * ( 3*sin(lat)^2 -1 )
% Coefficients are calculated so that the area weighted mean of the profile
% matches the surface mean temperatures in each sector.
%----------------
CTa(1,:) = [ 1 .5*(sin(fdiv)^2-1)                   ];
CTa(2,:) = [ 1 .5*(sin(fdiv)-sin(fdiv)^3)/(1-sin(fdiv)) ];
RTa      = [ Ta(1) Ta(2)]';
PTa      = CTa\RTa;
%----------------

% Sea-Ice
%----------------
  
% Calculate the latitude fi where Ta(fseaice) =Tseaice & Ta(fsnow)   =Tsnow 
  fseaiceo = pi/2-0.0001;
  fseaice15C = 0.6922*pi/2;
  fseaice = fland; %pi/2-0.0001;    %Calculate with no ice or snow
  fsnow=pi/2-0.0001;
  fmax     = max(fsnow , fseaice); % Used for albedo calculations

%----------------

% Calculate mean temperature of the high latitude atmosphere over ice free ocean
%----------------
THNI  =  1/(sin(fseaice)-sin(fdiv))*...
         [(PTa(1)-PTa(2)*.5)*(sin(fseaice)-sin(fdiv))+...
         0.5*PTa(2)*(sin(fseaice)^3-sin(fdiv)^3 )];

% Calculate ice free ocean area 
%----------------
aHLNI = olfHL*aH*(sin(fseaice)-sin(fdiv)); 

% Apply extra high latitude, longwave cloud forcing
aHLLWCF=aH*(sin(fseaiceo)-sin(fseaice15C));  
HLLWCF=30; %W/m2
QHLLWCF=HLLWCF*aHLLWCF;

if fseaice < fseaice15C
QHLLWCF=0;
end

% Air sea heat exchange, positive upwards
%----------------
QLL   = (lambda*(ToLL-Ta(1))-Dheat); % [W/m2]
QHL   =  lambda*(ToHL-THNI);         % [W/m2]

% Meridional temperature gradient at latitude = fdiv
%----------------
PTay  = 3*PTa(2)*sin(fdiv)*cos(fdiv);

% Temperature at latitude = fdiv
%----------------
PTad  = PTa(1) + PTa(2)*.5*(3*sin(fdiv).^2-1);

% Integrated mean incomming shortwave radiation for the low latitude sector 
%----------------
Aw    = 1/sin(fdiv)*Qo/4 * [ ...
         (al0-.5*al2)*(1-.5*Q2)*sin(fdiv) + ...
         .5*(Q2*(al0-.5*al2)+al2*(1-.5*Q2))*sin(fdiv)^3 + ...
         9/20*al2*Q2*sin(fdiv)^5 ];      %[W/m2]
% Subtract zero degree outgoing longwave radiation 
Aw    = Aw-Atot;                         %[W/m2]   

% Integrated mean incomming radiation for the high latitude sector, north of fd, 
% including icecover north of fseaice of fixed albedo ali. 
%----------------
Ac1 = (al0-.5*al2)*(1-.5*Q2) - (1-ali)*(1-.5*Q2);
Ac2 = .5*(Q2*(al0-.5*al2)+al2*(1-.5*Q2)) - (1-ali)*Q2*.5;
Ac3 = 9/20*al2*Q2;
Ac4 = -(al0-.5*al2)*(1-.5*Q2)*sin(fdiv) ...
      -.5*(Q2*(al0-.5*al2)+al2*(1-.5*Q2))*sin(fdiv)^3 ...
      - 9/20*al2*Q2*sin(fdiv)^5 ...
      +(1-ali)*(1-.5*Q2) + (1-ali)*Q2*.5;
% 'Ocean zone' (to pole): albedo changes at fmax
Aco = 1/(1-sin(fdiv))*Qo/4*...
      ( sin(fmax )*Ac1 + sin(fmax )^3*Ac2 + sin(fmax )^5*Ac3 + Ac4 );
% Land zone (to pole): albedo changes at fsnow
Acl = 1/(1-sin(fdiv))*Qo/4*...
      ( sin(fsnow)*Ac1 + sin(fsnow)^3*Ac2 + sin(fsnow)^5*Ac3 + Ac4 );
% Combine the two estimates

Ac  = Aco*olfHL+Acl*(1-olfHL) - Atot;
% Atmospheric freshwater flux into the cold sector
%----------------
Fw    = Cfw * exp(-5420/(PTad+273))*abs(PTay)^N;    %[m3/s]

% Atmospheric energy transport, sensible and latent
%----------------
Hash  = Cs * abs(PTay)^N;         %[W]
Halh  = Lv*Fw;                    %[W]
Ha    = Hash + Halh;              %[W]

%----------------
% Heat convergence/divergence
%----------------
QEBLL = -Ha + QLL*aLL   + (Aw - B*Ta(1))*aH*sin(fdiv)     ; %[W]
QEBHL =  Ha + QHL*aHLNI + (Ac - B*Ta(2))*aH*(1-sin(fdiv)) + QHLLWCF; %[W]
%----------------



return
