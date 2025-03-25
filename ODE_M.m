function [kLL,kHL,kAT,kLB,srLL,srHL,car]=...
    ODE_M(t,LL,HL,input,AT,LB,pCalLL,pCalHL,pOrgLL,pOrgHL)

% ----------------------------------------------------------------------------------
% ODE_PETMIsoM_New: Adición condición OrgFlx O2>minO y oxidación anaeróbica de CH4
% ----------------------------------------------------------------------------------

% Ordinary Differential Equations for Thilda
%
% Input :  t        - time
%          LL/HL    - ocean tracers
%          AT       - atmospheric tracers
% Data ordered accordingly:
% LL(tracer,level) contains low latiude ocean data
% HL(tracer,level) contains high latitude ocean data
% AT(tracer,zone) contains atmospheric data
% LB(tracer
% Tracers are: 
%                     1:   - temperature [oC]
%                     2:   - salinity []
%                     3:   - phosphate [mol/m3
%                     4:   - DIC [mol/m3] 
%                     5:   - DI^13C [mol/m3]
%                     6:   - DI^14C [mol/m3]
%                     7:   - alkalinity [eq/m3] 
%                     8:   - oxygen  [mol/m3] 
%                     9:   - Isotope ration d18O [-]
%                     10:  - CH_4 [mol/m3]
%                     11:  - 13CH_4 [mol/m3]
%                     12:  - nitrate [mol/m3] 
%                     13:  - NH_4 [mol/m3] 
%                     14:  - H_2S [mol/m3] 
%         AT        - Atmospheric tracers
%                     1:   - temperature [oC]
%                     2:   - pCH_4 [atm]
%                     3:   - pN_2O [atm]
%                     4:   - pCO_2 [atm]
%                     5:   - p^13CO_2 [atm]
%                     6:   - p^14CO_2 [atm]
%                     7:   - p^13CH_4
%                     8:   - pO_2 [atm]  
% Zones are (only temperature is not well mixed):    
%                     1:   - low latitude (LL)
%                     2:   - high latitude (HL)
%        LB          - Land Biosphere tracers
%                     1-4: ^12C
%                     5-8: ^13C
%                     9_12: ^14C 
%
%
% Output: kLL/kHL   - Rate of change (per sec) of low (kLL) and 
%                     high (kHL) oceanic tracers. Used for time-stepping.
%                     Dimensions are as input ocean data.
%         kAT       - Rate of change (per sec) of atmospheric tracers
%                     Dimensions are as input atmospheric data.  
%         kLB       - Rate of change (per sec) of Land Biosphere tracers
%                     dimensions are as input land biosphere data   
%         srLL/srHL - ocean source and sink terms. Used for diagnostics.


% Activate global Parameters
global nto sim GAw GAc RVolC13 RCorgC13 RCarbC13 Volin 
global d dm n dv aH aLL aHL Lxf Lyf rc fdiv weLL weHL sy rno rcp rcop
global dc14 R13pdb rVa Avg ULL UHL LfLL LfHL TSL TSH NPPIL NPPIH 
global minO minN k_Poxidw RorgLL_input mgt orgc13 lipc13
% Ocean and atmospheric exchanges (see functions) 
%-------------
[q,wLL,wHL,kvLL,kvHL,kh]          = OceExc_M;
%-------------

[QEBLL,QEBHL,QLL,QHL,aHLNI,Fw,Td] = AtmEnerBal_M(LL(1,1),HL(1,1),AT);
%-------------
% External inputs

% Organic matter oxidation (wildfire and thawing permafrost)
%*************************************************
% org_vel = AtmOrg_M(t,input);
% lip_vel = AtmLip_M(t,input);

org_rate = enforcing(t/sy/1000,input.time,input.org.rate)*mgt/sy;
org_vel(1) = org_rate;  %
org_vel(2) = org_rate*(orgc13*1e-3+1)*R13pdb;

lip_rate = enforcing(t/sy/1000,input.time,input.lip.rate)*mgt/sy;
lip_vel(1) = lip_rate;  %
lip_vel(2) = lip_rate*(lipc13*1e-3+1)*R13pdb;

[RcarLL,RcarHL,RorgLL,RorgHL,Wcarb,Wsil,Worg]=...
    ExtForce_M(AT);
%************************************************

% d18O ratio in atmospheric water vapor
%-------------

  [A18O] = WaVa18O(LL(1,1),Td,LL(9,1));
%-------------

% Air-sea exchange of DIC/pCO_2 (incl. isotopes) and O_2 
% Cosmogenic production of ^14C is also estimated (see function) 
%-------------

[asLL,P14C,CO2LL,CO3LL,HCO3LL,GAMLL] = GasExc_M(LL,AT,ULL,aLL);
[asHL,P14C,CO2HL,CO3HL,HCO3HL,GAMHL] = GasExc_M(HL,AT,UHL,aHLNI);
%-------------

% Land biosphere component and land-air exchangespOrgLL
%----------------
[al]=LandExcnew_M(AT,LB);



% Methane hydrate and CO2 sources and atmospheric methane sinks
%---------------------------------------------------------------
[mpr,mdr,MRR] = AtmMet_M(t,input,AT);
%------------------------------------------------------------



% Tracer sources and sinks (see function) 
%-------------

[srLL,rpl,fnpl,fcpl,fCO2l,fCO3l,fHCO3l] = OrgFlx_M(TSL,NPPIL,LL,AT,aLL  ,LfLL,1,CO2LL,CO3LL,HCO3LL,GAMLL,GAw,pCalLL,pOrgLL,RcarLL,RorgLL,1);
[srHL,rph,fnph,fcph,fCO2h,fCO3h,fHCO3h] = OrgFlx_M(TSH,NPPIH,HL,AT,aHLNI,LfHL,1,CO2HL,CO3HL,HCO3HL,GAMHL,GAc,pCalHL,pOrgHL,RcarHL,RorgHL,0);
%-------------

% Establish source terms for all tracers
%-------------
% Temperature and salinity only have surface sources/sinks
srLL(1,:)    = zeros(1,n); srLL(1,1) = -QLL*aLL/rc;
srHL(1,:)    = zeros(1,n); srHL(1,1) = -QHL*aHLNI/rc;
srLL(2,:)    = zeros(1,n); srLL(2,1) = 0;     %   +Fw*So;
srHL(2,:)    = zeros(1,n); srHL(2,1) = 0;     % -Fw*So;

% Combine surface fluxes and interior sources for DIC

dCH4=0.7e-9;
dCH4a=0.1*dCH4;                      %revised, 
fOH=1.01;

srLL(4,1) = srLL(4,1)+asLL(4);
srLL(4,:) = srLL(4,:)+dCH4*LL(10,:)*aLL.*GAw.*dv.*(LL(8,:)>minO)...
                     +dCH4*LL(10,:)*aLL.*GAw.*dv.*(LL(8,:)<=minO).*(LL(12,:)>minN)...
                     +dCH4a*LL(10,:)*aLL.*GAw.*dv.*(LL(8,:)<=minO).*(LL(12,:)<=minN);
srHL(4,1) = srHL(4,1)+asHL(4);
srHL(4,:) = srHL(4,:)+dCH4*HL(10,:)*aHL.*GAc.*dv.*(HL(8,:)>minO)...
                     +dCH4*HL(10,:)*aHL.*GAc.*dv.*(HL(8,:)<=minO).*(HL(12,:)>minN)...
                     +dCH4a*HL(10,:)*aHL.*GAc.*dv.*(HL(8,:)<=minO).*(HL(12,:)<=minN);

% Combine surface fluxes and interior sources for DI^13C
srLL(5,1) = srLL(5,1)+asLL(5);
srLL(5,:)= srLL(5,:)+dCH4*fOH*LL(11,:)*aLL.*GAw.*dv.*(LL(8,:)>minO)...
                    +dCH4*fOH*LL(11,:)*aLL.*GAw.*dv.*(LL(8,:)<=minO).*(LL(12,:)>minN)...
                    +dCH4a*fOH*LL(11,:)*aLL.*GAw.*dv.*(LL(8,:)<=minO).*(LL(12,:)<=minN);
srHL(5,1) = srHL(5,1)+asHL(5);
srHL(5,:) = srHL(5,:)+dCH4*fOH*HL(11,:)*aHL.*GAc.*dv.*(HL(8,:)>minO)...
                     +dCH4*fOH*HL(11,:)*aHL.*GAc.*dv.*(HL(8,:)<=minO).*(HL(12,:)>minN)...
                     +dCH4a*fOH*HL(11,:)*aHL.*GAc.*dv.*(HL(8,:)<=minO).*(HL(12,:)<=minN);

% Combine surface fluxes and interior sources for DI^14C
% Radiocarbon decay is given via the decay constant (dc14).

srLL(6,:)  = srLL(6,:)-dc14*LL(6,:)*aLL.*GAw.*dv;
srLL(6,1)  = srLL(6,1)+asLL(6);
srHL(6,:)  = srHL(6,:)-dc14*HL(6,:)*aHL.*GAc.*dv;
srHL(6,1)  = srHL(6,1)+asHL(6);

%Alkalinity,   added!!!

srLL(7,:) = srLL(7,:)+2*dCH4a*LL(10,:)*aLL.*GAw.*dv.*(LL(8,:)<=minO).*(LL(12,:)<=minN);
srHL(7,:) = srHL(7,:)+2*dCH4a*HL(10,:)*aHL.*GAc.*dv.*(HL(8,:)<=minO).*(HL(12,:)<=minN);

% Combine surface fluxes and interior sources for oxygen
srLL(8,1) = srLL(8,1)+asLL(8);
srLL(8,:) = srLL(8,:)-2*dCH4*LL(10,:)*aLL.*GAw.*dv.*(LL(8,:)>minO);
srHL(8,1) = srHL(8,1)+asHL(8);
srHL(8,:) = srHL(8,:)-2*dCH4*HL(10,:)*aHL.*GAc.*dv.*(HL(8,:)>minO);

% The d18O isotopic ratio is forced only by surface freshwater fluxes

srLL(9,:)   = zeros(1,n); srLL(9,1) = -Fw*A18O;
srHL(9,:)   = zeros(1,n); srHL(9,1) = +Fw*A18O;

% Anerobic oxidation of methane with nitrate
% 5CH4+8NO3 = 5CO2 +4N2 +14H2O
srLL(12,:) = srLL(12,:)-1.6*dCH4*LL(10,:)*aLL.*GAw.*dv.*(LL(8,:)<=minO).*(LL(12,:)>minN);
srHL(12,:) = srHL(12,:)-1.6*dCH4*HL(10,:)*aHL.*GAc.*dv.*(HL(8,:)<=minO).*(HL(12,:)>minN);

%Anerobic oxidation of methane with sulphate
%CH4+SO4 = CO2 +H2S +2H20
srLL(14,:) = srLL(14,:) + dCH4a*LL(10,:)*aLL.*GAw.*dv.*(LL(8,:)<=minO).*(LL(12,:)<minN);
srHL(14,:) = srHL(14,:) + dCH4a*HL(10,:)*aHL.*GAc.*dv.*(HL(8,:)<=minO).*(HL(12,:)<minN);


% methane C12
srLL(10,:)  =-dCH4*LL(10,:)*aLL.*GAw.*dv.*(LL(8,:)>minO)...
    -dCH4*LL(10,:)*aLL.*GAw.*dv.*(LL(8,:)<=minO).*(LL(12,:)>minN)...
    -dCH4a*LL(10,:)*aLL.*GAw.*dv.*(LL(8,:)<=minO).*(LL(12,:)<=minN);
srLL(10,1) = srLL(10,1)+asLL(2);
srHL(10,:) = -dCH4*HL(10,:)*aHL.*GAc.*dv.*(HL(8,:)>minO)...
    -dCH4*HL(10,:)*aHL.*GAc.*dv.*(HL(8,:)<=minO).*(HL(12,:)>minN)...
     -dCH4a*HL(10,:)*aHL.*GAc.*dv.*(HL(8,:)<=minO).*(HL(12,:)<=minN);
srHL(10,1) = srHL(10,1)+asHL(2);
%-------------

%methane C13
srLL(11,:)  =-dCH4*fOH*LL(11,:)*aLL.*GAw.*dv.*(LL(8,:)>minO)...
             -dCH4*fOH*LL(11,:)*aLL.*GAw.*dv.*(LL(8,:)<=minO).*(LL(12,:)>minN)...
             -dCH4a*fOH*LL(11,:)*aLL.*GAw.*dv.*(LL(8,:)<=minO).*(LL(12,:)<=minN);
srLL(11,1) = srLL(11,1)+asLL(7);

srHL(11,:) = -dCH4*fOH*HL(11,:)*aHL.*GAc.*dv.*(HL(8,:)>minO)...
             -dCH4*fOH*HL(11,:)*aHL.*GAc.*dv.*(HL(8,:)<=minO).*(HL(12,:)>minN)...
             -dCH4a*fOH*HL(11,:)*aHL.*GAc.*dv.*(HL(8,:)<=minO).*(HL(12,:)<=minN);
             
srHL(11,1) = srHL(11,1)+asHL(7);


for m=3:30   
    
    srLL(10,m)=srLL(10,m)+0.84*mpr(3);
    srHL(10,m)=srHL(10,m)+0.16*mpr(3);
    
    srLL(11,m)=srLL(11,m)+0.84*mpr(4);
    srHL(11,m)=srHL(11,m)+0.16*mpr(4);

end

% Calculate time derivative of ocean tracers
%-------------

kvLL(1:n-1)=kvLL(1:n-1).*GAw(1:n-1);
kvHL(1:n-1)=kvHL(1:n-1).*GAc(1:n-1);
kh(1:n)=kh.*GAw(1:n);

%

    kLL     = 1./repmat((aLL*[dm d*GAw(2:n)]),nto,1).*( ...
             + [-kvLL(1)*aLL/(3*d)*(8*LL(:,1)-9*LL(:,2)+LL(:,3)) ...
                 repmat(kvLL(:,2:n-1),nto,1)*aLL/d.*(LL(:,3:n)-LL(:,2:n-1)) zeros(nto,1)] ...
	     - [zeros(nto,1) -kvLL(1)*aLL/(3*d)*(8*LL(:,1)-9*LL(:,2)+LL(:,3)) ...
                 repmat(kvLL(:,2:n-1),nto,1)*aLL/d.*(LL(:,3:n)-LL(:,2:n-1))]...
             - repmat(kh.*Lxf/Lyf.*dv,nto,1).*(LL-HL) ...
	     + [wLL(1)*LL(:,1) repmat(wLL(2:n-1),nto,1).*( LL(:,2:n-1)+LL(:,3:n) )*.5 zeros(nto,1)] ...
             - [zeros(nto,1) wLL(1)*LL(:,1) repmat(wLL(2:n-1),nto,1).*( LL(:,2:n-1) + LL(:,3:n) )*.5] ...
             - repmat((q>0).*(q-Fw),nto,1).*LL ...
             - repmat((q<0).*q,nto,1).*HL ...
             + srLL );

    kHL     = 1./repmat((aHL*[dm d*GAc(2:n)]),nto,1).*( ...         
             + [-kvHL(1)*aHL/(3*d)*(8*HL(:,1)-9*HL(:,2)+HL(:,3)) ...
                 repmat(kvHL(:,2:n-1),nto,1)*aHL/d.*(HL(:,3:n)-HL(:,2:n-1)) zeros(nto,1)] ...
	     - [zeros(nto,1) -kvHL(1)*aHL/(3*d)*(8*HL(:,1)-9*HL(:,2)+HL(:,3)) ...
                 repmat(kvHL(:,2:n-1),nto,1)*aHL/d.*(HL(:,3:n)-HL(:,2:n-1))]...
             + repmat(kh.*Lxf/Lyf.*dv,nto,1).*(LL-HL) ...
	     + [wHL(1)*HL(:,1) repmat(wHL(2:n-1),nto,1).*( HL(:,2:n-1)+HL(:,3:n) )*.5 zeros(nto,1)] ...
             - [zeros(nto,1) wHL(1)*HL(:,1) repmat(wHL(2:n-1),nto,1).*( HL(:,2:n-1) + HL(:,3:n) )*.5] ...
             + repmat((q>0).*(q-Fw),nto,1).*LL ...
             + repmat((q<0).*q,nto,1).*HL ...
             + srHL );
    
% Calculate time derivative of atmospheric tracers

nta = 12; 

for i=1:nta
% Atmospheric temperature, gasses and tracers
  if i==1
      kAT(i,:)  = [QEBLL QEBHL]./(rc*aH*[weLL*sin(fdiv) weHL*(1-sin(fdiv)) ]); 
  elseif i==2
      kAT(i,:)  = -rVa*(asLL(2) + asHL(2)-al(11) + mdr(1) - mpr(1));% 
  elseif i==3
      kAT(i,:)  = -rVa*(asLL(3) + asHL(3)-al(13) + AT(3,1)*1/(rVa*150*sy));                    % N2O lifetime 150 yrs;  
  elseif i==4
      kAT(i,:)  =-rVa*(asLL(4) + asHL(4) - al(9) - mdr(1)+...%-(asLL(2) + asHL(2)
           Wcarb + 2*Wsil -Worg - Volin- org_vel(1)-lip_vel(1)); 
      %last term is volcanic input     
  elseif i==5 
      kAT(i,:)  =-rVa*( + asLL(5) + asHL(5) - al(10) - mdr(2) +...%
         (Wcarb+2*Wsil)*AT(5,1)/AT(4,1) - ...
          Worg*(R13pdb*(RCorgC13+1))-Volin*(R13pdb*(RVolC13+1))-org_vel(2) -lip_vel(2));   
  elseif i==6 
      kAT(i,:)  =-rVa*(asLL(6) + asHL(6)-al(18)+(Wcarb+2*Wsil)*AT(6,1)/AT(4,1)) - dc14*AT(6,1) + aH*P14C*rVa/Avg;
   elseif i==7
      kAT(i,:)  = -rVa*(asLL(7)+asHL(7)-al(12) + mdr(2) - mpr(2)); 
  elseif i==8 
       kAT(i,:)  = -rVa*(asLL(8) + asHL(8) + 1.1*(al(9)+al(11)) + 2*(mdr(1)-al(11)) - ...     %revised 19/10/2012 to contain CH4 oxidation
     (Worg+Volin*(RCarbC13-RVolC13)/(RCarbC13-RCorgC13))*rcop/rcp-rno*(RorgLL+RorgHL));
  end 
end

%---------
%***************************************************************
car = nan([5,1]);                                   
car(1) = MRR(1);                                   % C mol/s associated with methane hydrate
car(2) = MRR(2);                                   % C mol/s associated with thermogenic hydrate
car(3) = org_vel(1);                               % C mol/s organic matter oxidation (wildfire and thawing permofrost)
car(4) = lip_vel(1);                               % C mol/s associated with LIP
car(5) = 2*Wsil;                                   % C mol/s silicate weathering bury


%***************************************************************
                                                              
% Calculate time derivative of C12 and C13 of land biomasses
% kLB=zeros(1:12,:);
for i=1:12
  if i==1
     kLB(i,:)  = al(1); 
 elseif i==2
     kLB(i,:)  = al(2);    
 elseif i==3
     kLB(i,:)  = al(3);   
 elseif i==4
     kLB(i,:)  = al(4);
 elseif i==5
     kLB(i,:)  = al(5); 
 elseif i==6
     kLB(i,:)  = al(6);    
 elseif i==7
     kLB(i,:)  = al(7);   
 elseif i==8
     kLB(i,:)  = al(8); 
 elseif i==9
     kLB(i,:)  = al(14); 
 elseif i==10
     kLB(i,:)  = al(15);    
 elseif i==11
     kLB(i,:)  = al(16);   
 elseif i==12
     kLB(i,:)  = al(17);      
 end 
end 

return

