function [mpr,mdr,MRR] = AtmMet_M(t,input ...
    ,AT)

% Input : t - time
%         At - atmospheric tracers
% Output: mpr(1) - MH release to atmosphere, 12C. mol/sec 
%         mpr(2) - MH release to atmosphere, 13C. mol/sec 
%         mpr(3) - MH release to one ocean layer, 12C. mol/sec 
%         mpr(4) - MH release to one ocean layer, 13C. mol/sec
%         mdr(1) - methane oxydation in the atmosphere , 12C, mol/sec
%         mdr(2) - methane oxydation in the atmosphere , 13C, mol/sec


global sy rVa mgt R13pdb mdts2 methc13 mettc13

pCH4o =  0.72e-6;                           %Pre-industrial methane concentration
pCH4  =  AT(2,1);
fatm=1.0039;                                % C fractionation in atmospheric oxidation
M     = (pCH4-pCH4o)/pCH4o;

                                             
RCH4o   =  1/(rVa*mdts2*sy);                %PA decay rate for methane, residence time 8.4 yrs
RCH4    =  RCH4o*(1-0.78*M/(M+11));         % Decay rate is the inverse of residence time 

%***********************************
nlay  = 30-3+1;                                 % number of ocean layers receiving methane 
%**********************************

% if (t/sy>=MH(2))
%     MRR  =  (mgt*MH(3)*(t/sy-MH(2))^4*exp(-MH(4)*(t/sy-MH(2))))/sy; 
% 
% else
%     MRR =0;
% end         
% input_end = MH(5)+MH(2);

%------------------------------------------------

%%%%%%% Hydrate C release rate derive from oxidation mol/sy
% input_end = input.meth.sta+ input.meth.dur;
% if (t/sy>=input.meth.sta && t/sy<= input_end)
%     [~,minIndex] = min(abs(input.meth.rate(:,1)-t/sy));
%     MRRh  = input.meth.rate(minIndex,2)*mgt/sy;
% else
%     MRRh = 0;
% end
%------------------------------------------------

MRRh = enforcing(t/sy/1000,input.time, input.meth.rate)*mgt/sy;

MRRt = enforcing(t/sy/1000,input.time, input.mett.rate)*mgt/sy;


%%%%%%% Thermogenic C release rate derive from oxidation mol/sy
% input_end = input.mett.sta+ input.mett.dur;
% if (t/sy>=input.mett.sta && t/sy<= input_end)
%     [~,minIndex] = min(abs(input.mett.rate(:,1)-t/sy));
%     MRRt  = input.mett.rate(minIndex,2)*mgt/sy;
% else
%     MRRt = 0;
% end

MRR = [MRRh, MRRt];

mpr(1) =  input.meth.mfrac*MRRh + input.mett.mfrac*MRRt;
mpr(2) =  (input.meth.mfrac*MRRh*(methc13*1e-3+1))*R13pdb + (input.mett.mfrac*MRRt*(mettc13*1e-3+1))*R13pdb;   % methc13 is del 13c of input methane
mpr(3) =  (1-input.meth.mfrac)/nlay*MRRh + (1-input.mett.mfrac)/nlay*MRRt;      
mpr(4) =  (1-input.meth.mfrac)/nlay*MRRh*R13pdb*(methc13*1e-3+1)+ (1-input.mett.mfrac)/nlay*MRRt*R13pdb*(mettc13*1e-3+1);



mdr(1) =  pCH4*RCH4;  % atmospheric methane oxidation to CO2; pCH4 is the concentration of CH4 in atmosphere
mdr(2) =  (pCH4*RCH4*AT(7,1)/pCH4)*fatm;



return


