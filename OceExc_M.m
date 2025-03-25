function [q,wLL,wHL,kvLL,kvHL,kh] = OceExc_PETMIsoM

% Output
%          q         - Horizontal exchange between layers, posetive poleward [m3/s]
%          wLL/wHL   - Vertical exchange between layers, posetive upwards [m3/s] 
%          kvLL/kvHL - Vertical diffusion [m2/s]
%          kh        - Horizontal diffusion [m2/s]

% Activate global parameters
global d dm n aLL


mpl=1;            %Multiplier of physical exchange, Nov. 2012
% Ocean overturning, poleward in the surface layers 
% and equatorward at the ocean bottom - zero in between.
qmax  =  mpl*aLL*1.8e-8;                % [m3/s], st. val. 2.7.5e-8
q(1)  =  qmax;      
q(n)  = -qmax;

% Upwelling/downwelling given by continuity, zero at the ocean bottom
wLL(n:-1:1) =  [0 -cumsum(q(n:-1:2))]; % [m3/s]
wHL(n:-1:1) =  [0  cumsum(q(n:-1:2))]; % [m3/s]

% Kv given at layer interfaces 
kvLL(1:n-1) = mpl*2e-5.*(1+5.5*(1-exp(-100*(1:n-1)/4000)));                        % [m^2/s], st. val 3.8e-5
kvHL(1:n-1) = mpl*230e-5;              % [m^2/s], st. val 240e-5

% Horizontal diffusion
khd   =  mpl*1.65e3;                   % Constant (deep) diff. [m2/s], st. val 2.6e3    
khs   =  0*khd;                        % Surface intensified diff. [m2/s]
zkhs  =  200;                          % Depth scale for surface intens. diff. [m] 
z     = [-dm/2 -dm/2-((1:n-1)*d)];     % Center of layers [m]
kh    =  (khd+(khs*exp(z/zkhs)));      % Combined horizontal diffusion [m2/s]

return

