function ParVal_M

% Model settings
global nto 
%-----
nto     = 14;      % Number of ocean tracers in the list (1 to nto)
                  
disp(strcat(['Number of ocean tracers: ' int2str(nto)]))

% Useful constants
global  sy rc mgt
%-----
sy      =  31536000;	              % Seconds per year [s/yr] 
rc      =  4.186e6;	                  % Scaled heat capacity [J/m3/K]  
mgt     =  8.326e13;                  % Moles C per GtC
%-----

% Geometry
global  dm d n dv olf olfHL r aH fdiv fland aHL aLL Lxf Lyf
%-----
dm      =  100;                           % Surface  layer thickness [m]
d       =  100;                           % Interior layer thickness [m]
n       =  55;                            % Total number of layers [-]
dv      =  [dm d*ones(1,n-1)];            % Vector of layer thichnesses [m]
olf     =  1.0173*270/360;                % Ocean to land fraction[-],, for%extension of old topo
r       =  6.371e6;	                      % Earth Radius [m]
aH      =  2*pi*r^2;                      % Hemispheric area [m2]
fdiv    =  52*pi/180;                     % High - low latitude division [rad], st. val 52
fland   =  71*pi/180;                     % High latitude poleward ocean limit , for extension of old topo 
olfHL   = olf*(sin(fland)-sin(fdiv))/(sin(fland)-sin(fdiv));%revised 11/06/13
aHL     =  olfHL*aH*(sin(fland)-sin(fdiv)); % High latitude ocean area [m2), 
aLL     =  olf*aH*sin(fdiv);	          % Low latutude ocean area [m2], revised 11/06/13 
                                          % Note: olf, fdiv and HL and LL ocean surface areas
                                          % adjusted slightly to give internal agreement. 
Lxf     =  2*pi*r*olf*cos(fdiv);          % Zonal ocean length scale at fdiv [m]
Lyf     =  0.5*r;                         % Meridional ocean length scale at fdiv [m]
%-----



global dc14 gLL gHL C14a rVa Avg ULL UHL R13pdb R14oas kwo Tgas P14Co rpm
%-----

rVa     =  1/(1.7676e20/2);           % Inverse atmospheric volume (one hemisphere) meassure [atm/mol], update Jan 2012  
Avg     =  6.024e23;                  % Avogadro number [atoms/mol]
kwo     =  1.5*4.72e-5;               % piston velocity for gas exchange, st. val 4.72 or 17.1 cm/hr
ULL     =  8;                         % Wind speed, low lat. [m/s], st. val 5.6, Jan 29, 2006
UHL     =  8;                         % Wind speed, high lat. [m/s], st. val 5.6, Jan 29 2006 
dc14    =  3.84e-12;                  % C14 decay rate [s-1]
gLL     =  2.32e-7;                   % Air-sea gas exchange velocety [m2/s]
gHL     =  2.32e-7;                   % Air-sea gas exchange velocety [m2/s]
Tgas    =  0;                         % Temperature offset for calculating HL gas exchange
C14a    =  1;                         % Fixed (arbitrary) atmospheric 14^C concentration
R13pdb  = 0.0112372;                  % 'Pee-Dee Belemnite' 13C/12C standard []
R14oas  = 1.176e-12;                  % 'Oxalic Acid Standard' for the 14C/12C ratio []
P14Co   = 1.8752e+004;                % Steady state, Holocene C14 production value [atoms/m^2]
%-----


% Redfield atomic ratios, remineralization scales and new production 
global rcp rno rcop rnp rcd rda rca minO minN rDNno3 rDNnh4 rRBonh4...
    rRBnnh4 rSRh2s rSRnh4 rRBoh2s LfLL LfHL lmdcar lmdorg lmdcorg rf methc13 orgc13 lipc13 mettc13
%-----
rcp     =  106;                       % Carbon to phosphate, organic matter [-]
rno     =  -32;                       % Oxygen to nitrogen, organic matter [-],    
rcop    = -118;                       % C+H to phosphate, organic matter [-], 
rnp     =  16;                        % Nitrate to phosphate, organic matter [-]
rcd     =  1;                         % Carbonate to DIC [-]
rca     =  2;                         % Carbonate to ALK [-]
rda     = -16;                        % DIC to ALK [-]
rDNno3  = -94.4;                      % NO_3 to PO_4 in denitrification []; revised from -84.8 
rDNnh4  =  16;                        % NH_4 to PO_4 in denitrification [] 
rRBonh4 = -2;                         % O_2 to NH_4 at redox boundaries [] 
rRBnnh4 =  1;                         % NO_3 to NH_4 at redox boundaries [] 

rSRh2s  =  59;                        % H_2S to PO_4 in sulphate reduction []; revised from 53 
rSRnh4  =  16;                        % NH_4 to PO_4 in sulphate reduction [] 
rRBoh2s = -2;                         % O_2 to H_2S at redox boundaries [] 


%************************
methc13 = -60;                        %del13C of methane input deault: -40   methane hydrate: -60
mettc13 = -40;                                 % methermogenic methane: -40;
orgc13 = -25;                         %del13C of organic matter
lipc13 = -15;                         %del13C associated with lip 
%************************

lmdcar  = 1/3000;                     %water column dissolution scale for CaCO3 , st. val. 2000 [m^-1]
lmdorg  = 1/750;                      %water column remineralization scale for P, st. val 780 [m-1]
lmdcorg = 1/1050;                     %water column remineralization scale for Org C, st. val 1000 [m-1]

minO    =  3e-3;                   % Min. O_2 conc. for oxidation of POM, implies onset of DN [mol/m3] default: 3e-3 mol/m3
minN    =  3e-5;                      % Min. NO_3 conc. for denitrification of POM, implies onset of SR [mol/m3]
LfLL    =  1;                         % Limitation factor for new production, low lat. [-], st. val 1
LfHL    =  0.36;                      % Limitation factor for new production, high lat. [-]
rf      =  0.9;                       % remineralization fraction of organic matter falling on ocean bottom, 
                                      % st. val 0.9, equal to 1 minus the burial fraction
rpm     = 0.4517;                     % maximum rain ratio  


% Atmosphere   
global SCred Tseaice Tsnow weLL weHL pCO2int pCH4int pN2Oint Htem mdts mdts2 albch

%-----
SCred   = 0.005;                      % Fractional reduction of solar constant, 0.005 at 50 Myr BP
albch = -0.05;                        % Change in global-mean background albedo, 0.05 at 50 Myr BP; changed to 0.02 Aug. 5, 2013
Tseaice = -5;                         % Fixed ice line temperature,  st. val. -5
Tsnow   = 0;                          % Fixed snow line temperature, st. val.  0 
weLL    =  5;                         % Atmospheric heat capacity (water equivalent) [m]
weHL    =  20;                        % Atmospheric heat capacity (water equivalent) [m]
pCO2int = 278e-6;                     % Inital mean atmposhere pCO2
pCH4int = 0.72e-6;                    % Initial mean atmospheric CH4
pN2Oint = 0.27e-6;                    % Initial mean atmospheric N2O
Htem = 15.0;                          % Initial atmospher temperature
mdts  = 6.9;                          % Initial methane decay time scale (yr) 
mdts2 = 9.5;                          % Revised initial methane decay scale

% Ocean
global So Cao Pint DICint ALKint DelCO3LL DelCO3HL TSL TSH NPPIL NPPIH
%-----
So      =  33.8;                      % Mean reference salinity for PETM/Eocene[]
Cao     =  0.0182;                    % Ocean mean calcium concentration, modern value 0.01028, P/E value 0.0182 
Pint    = 2.11e-3;                    % Initial ocean mean phosphate concentration
DICint  = 2.33;                       % Initial ocean mean DIC concentration
ALKint  = 2.443;                      % Initial ocean mean alkalinity concentration
DelCO3LL = 0.2070;
DelCO3HL = 0.0831;
TSL=21.254;                           % low lat ocean surface temp , PI value 21.234, used only to calculate rp 
TSH=5.0;                              % high lat ocean surface temp, PI value 0.108; value of 5 gives double rph(PI)
NPPIL= 2.1488e-10;
NPPIH=1.6916e-09;

% Land Biosphere
global Q10 CO2fer Q10met aNSint aHLNIint FacVeg

Q10      = 2;                       % Q10 for soil = litter respiration rate, st. val 2.0
Q10met   = 2;                       % Q10 for terrestrial methane release, st. val 2.0
CO2fer   = 0.37;                    % New value based on C4MIP results, Friedlingstein et al 2006, Changed January 2012)
fsnowint = 55.829*pi/180;           % PI snowline latitude
fseaiceint = 63.551*pi/180;         % PI snowline latitude
aNSint   = (1-olf)*aH*(sin(fsnowint));  %PI snowfree land area
aHLNIint = olf*aH*(sin(fseaiceint)-sin(fdiv)); 
FacVeg = 1.3;                       % ratio of vegetation area to PI vegetaion area, 1.3 for Eocene  

%Sediment model

global kcalcite NCFLL NCFHL CAFLL CAFHL Mca rhom rhoc Moc compc k_Poxidw

Mca   = 100;                         % Molar weight of CaCO3 (40+12+16x3) [g/mol]
rhom  = 2.7;                         % Density of mineral fraction in sediments [g/cm3]
rhoc  = 1.1;                         % density if organic fraction in sediments [g/cm3]
Moc   = 12;                          % Molar weight of C
compc = 2.7;                         % total organic to org C ration in organic matter [g/g]
kcalcite   = 0.0015;                 % CaCO3 dissolution rate coefficient, st. val 0.001 [1/day];
NCFLL  = 0.3;                        % Open ocean, non-Calcite flux to sediment,Low Lat., st. val 0.3 [g/(cm2*kyr)]
NCFHL  = 0.3;                        % Open ocean, non-Calcite flux to sediment,High Lat., st. val 0.6 [g/(cm2*kyr)] 
CAFLL  = 20;                         % Amplification factor for NCF at the "coast", Low Lat., st. val. 20
CAFHL  = 20;                         % Amplification factor for NCF at the "coast", High Lat., st. val. 10 

% External inputs
global BCarbPA BCorgPA RCarbC13 RCorgC13 RVolC13 RCfosC13 RPPA RAtmC13 swf cwf Volo wfac wfacP Volin wfacC
swf      =     0.85;                 % Silicate weathering fraction, st. val 0.15
RPPA     =  0.83168e3;               % Pre-Anthropogenic river input of P, based on -sum(RMorg), see Orgfix
BCarbPA  =   2.6522e5; %4.4028       % PA burial of carbonate C, based on -sum(RMcar),see OrgFlx, mol/s
BCorgPA  =  0.95988e5;               % PA burial of organic C, based on -sum(RMorg)
RCarbC13  =  0.0011522;              % assumed delta C13 of carbonate C, calculated from burial

RCorgC13 = -0.023168;                % "    "     "     "   of organic C, calculated from burial 
RVolC13  =  -0.005;                  % "    "     "     "   volcanic CO2 entering the atmosphere, st.val.-5 per mil 
RCfosC13 = -0.028;                   % "    "     "     "   of fossil fuel, -28 per mil

RAtmC13   = -0.0064; 

Volo      =  swf*BCarbPA/(1+swf)*(RCarbC13-RCorgC13)/(RVolC13-RCorgC13);
tfac      =1.3;                      % see line below,
Volin    =  tfac*Volo;               % assumed volcanic input of CO2
wfac      = 0.65;                    % multiple of PI weathering
wfacP     = 0.8468; %0.65;           % multiple of PI weathering
wfacC     = 1.0086;

%***************************************
k_Poxidw = 0;                     % Oxidative fraction of P weathering
%***************************************

return
