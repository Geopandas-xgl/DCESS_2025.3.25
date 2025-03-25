function  [pCal,pOrg,dwpCal,dwpOrg,fim,wsed]=SMtdCorgNew2_PETMIsoM(TS,LH,rp,NCF,CAF,LHCO3,LHCO3s,NP,dwpCal,dwpCalss,dwpOrg,dwpOrgss,fim,wsed,ko,h)

% "Complete" time dependent sediment model in a 10 cm thick, bioturbated
% sediment layer, driven by water column concentrations from the model
% and carbonate flxes,simpliest version with dissolution not depending
% on organic carbon remineralization, approximate carbonate solution
% (ignoring [CO2]),all 5 sediment layers, Carbonate and Non-carbonate,
% fractions only.

% Modified August 15, 2005 by Gary

warning off all

global d dm sy n rcp rca rda kcalcite NCFLL NCFHL CAFLL CAFHL lmdcar lmdorg lmdcorg Mca rhom rpm

TLH   = LH(1,:)+273.16;
CO3LH = LHCO3*1e-6;
CO3LHs= LHCO3s*1e-6;
kc    = kcalcite/(3600*24); % Calcite dissolution rate constant,
                            % [s-1], st. val, 0.0003
NCFo  = NCF/1000/sy;        % Non-carbonate rain [g/cm^2/s]

%carbon burial parameters
Moc =  12;
rhoc = 1.1;              % organic matter density, g per cm3
compc = 2.7;             % g organic matter per g C    
Dbo  = 0.0232*18.774/sy; % Archer et al 2002, surface value
                         %18.774 facter is 31.5^^0.85
jo2o=1e-9;                                
R = 1.4;                 %O to C ratio

% Layers (5)
%-------------
zb   = [0 0.2 0.5 1 1.8 3.2 6 10]; %Chosen layer boundaries in sediment layer, cm 
difzb =[zb(2)-zb(1) zb(3)-zb(2) zb(4)-zb(3) zb(5)-zb(4)...
    zb(6)-zb(5) zb(7)-zb(6) zb(8)-zb(7)];
dtc   = ko*h*sy;         %time step, must be consistent with time step in Thilda!
fzo  = 1;                %porosity at the sediment surface

% Surface flux of CaCO3 and organic carbon 
RMcar = -NP*rcp*rp/(1e4);  % [mol/cm2/s]
RMorgc=-NP*rcp/(1e4);

%-------------
% Loop layers
%-------------
for i = 1:n
%-------------

  % Non carbonate flux, rapidly decreasing with depth to
  % mimic high levels in coastal areas ?
  NC   = CAF*NCFo*exp((1-i)/2)+ NCFo;    % [g/cm2/s]
 % Flux of calcite, exponential decay
  Fcal  = RMcar*exp(-lmdcar*(i*d-dm));      % [mol/cm2/s]
  
  %organic carbon flux, exponential decay
  FC = RMorgc*exp(-lmdcorg*(i*d-dm));    %Org C flux in mole C cm-2 s-1  
%   if LH(8,i)>1e-8
%      IniO2 = LH(8,i)*1e-6;
%   else
%       IniO2 = 1e-8;
%   end
if LH(8,i)>5e-3
     IniO2 = LH(8,i)*1e-6;
  else
      IniO2 = 5e-9;
end

   Db=Dbo/18.774*((FC*1e6*sy)^0.85)*IniO2/(IniO2+20*1e-9);
  %Db=Dbo;
  jo2=jo2o*Db/Dbo;
  beta=0.1*(FC*1e6*sy/31.5)^-0.3;
  %beta=0.025;

  % Temperature dependent diffusion of O2 and CO3
  My    = 2.31./(1 + 0.036 *(TLH(i)-273.16) + 1.85e-4*(TLH(i)-273.16).^2);
  Do2   = 1.25e-5*TLH(i)./278.16 * 1.95./My; 
  DCO3  = 7.9E-6 * (TLH(i)/291.16) * (1.3525/My);  

      PerCaCO3o=dwpCalss(i);
      alf = 0.25*PerCaCO3o+3*(1-PerCaCO3o);
      fimax = 1-(0.483+0.45*PerCaCO3o)/2.5;
      fimeano = fimax-alf*(1-fimax)*(exp(-zb(8)/alf)-1)/zb(8);  %mean porosity in 10 cm bioturbated layer
      fiboto=fimax-alf*(fzo-fimax)*(exp(-zb(8)/alf)-...         %mean porosity of bottom layer
                  exp(-zb(7)/alf))/(zb(8)-zb(7)); 
          
      PerOrgo = dwpOrgss(i);
      rhmeano =(1-PerOrgo)*rhom+PerOrgo*rhoc;   
    
     PerCaCO3 = dwpCal(i);
     PerOrg = dwpOrg(i);
     rhmean =(1-PerOrg)*rhom+PerOrg*rhoc;   
     fimeanb=fim(i);
     w=wsed(i);

      
    if CO3LH(i) >= CO3LHs(i)
      % Case: super saturated, no calcite dissolution:  
      %-------------
     alf = 0.25*PerCaCO3+3*(1-PerCaCO3);
     fimax = 1-(0.483+0.45*PerCaCO3)./2.5; 
     fimean = fimax-alf*(fzo-fimax)*(exp(-zb(8)/alf)-1)/(zb(8));
     fibot = fimax+(fzo-fimax)*exp(-zb(8)/alf);
     Tcaldis=0;Tcaldis1=0;Tcaldis2=0;Tcaldis3=0;Tcaldis4=0;Tcaldis5=0;
     Tcaldis6=0;Tcaldis7=0;s3=0;
      
     % Calculate mean sediment formation factor F (fi to the minus three)
    
     for k =1:7
    
     % calculate mean sediment form factor F for each sediment model layer   
     F(k)=(zb(k+1)-zb(k))/(fimax.^3*(zb(k+1)-zb(k))-3*fimax.^2*(fzo-fimax).*alf.*...
     (exp(-zb(k+1)./alf) -exp(-zb(k)./alf))-3/2*fimax.*(fzo-fimax).^2*alf.*...
     (exp(-2*zb(k+1)./alf) -exp(-2*zb(k)./alf))-1/3*(fzo-fimax).^3*alf.*...
     (exp(-3*zb(k+1)./alf) -exp(-3*zb(k)./alf)));
     fi(k)=fimax-alf*(fzo-fimax)*(exp(-zb(k+1)/alf)-...
                  exp(-zb(k)/alf))/(zb(k+1)-zb(k)); 
     end
     
      fimean = sum(fi.*difzb)/zb(8);
      
    else
      % Case: under saturated, calcite dissolution:  
      %-------------
           alf = 0.25*PerCaCO3+3*(1-PerCaCO3);
           fimax = 1-(0.483+0.45*PerCaCO3)/2.5;
           fimean = fimax-alf*(fzo-fimax)*(exp(-zb(8)/alf)-1)/(zb(8));
           fibot = fimax+(fzo-fimax)*exp(-zb(8)/alf);
           
            for k =1:7
    
              F(k)=(zb(k+1)-zb(k))/(fimax^3*(zb(k+1)-zb(k))-...
              3*fimax^2*(fzo-fimax)*alf*(exp(-zb(k+1)/alf)...
              -exp(-zb(k)/alf))-3/2*fimax*(fzo-fimax)^2*alf*...
             (exp(-2*zb(k+1)/alf)-exp(-2*zb(k)/alf))-...
             1/3*(fzo-fimax)^3*alf*(exp(-3*zb(k+1)/alf)...
             -exp(-3*zb(k)/alf)));
            fi(k)=fimax-alf*(fzo-fimax)*(exp(-zb(k+1)/alf)-...
                  exp(-zb(k)/alf))/(zb(k+1)-zb(k)); 
            end
            
             fimean = sum(fi.*difzb)/zb(8);
             
          %Apply analytical solutions for CO3 concentration,
      
           
      s3 = ((kc*PerCaCO3*rhmean*(1-fi(1))*F(1))/(CO3LHs(i)*Mca*DCO3))^0.5;
 
      
          alf2=((1-fi(2))*F(2)/((1-fi(1))*F(1)))^0.5;
          alf3=((1-fi(3))*F(3)/((1-fi(1))*F(1)))^0.5;
          alf4=((1-fi(4))*F(4)/((1-fi(1))*F(1)))^0.5;
          alf5=((1-fi(5))*F(5)/((1-fi(1))*F(1)))^0.5;
          alf6=((1-fi(6))*F(6)/((1-fi(1))*F(1)))^0.5;
          alf7=((1-fi(7))*F(7)/((1-fi(1))*F(1)))^0.5;
          
          
          co=(CO3LH(i)- CO3LHs(i));
          
          c1=exp(s3*zb(2));
          c2=exp(-s3*zb(2));
          c3=exp(alf2*s3*zb(2));
          c4=exp(-alf2*s3*zb(2));
          c5=F(1)/F(2)*alf2*c3; 
          c6=F(1)/F(2)*alf2*c4;
          c7=exp(alf2*s3*zb(3));
          c8=exp(-alf2*s3*zb(3));
          c9=exp(alf3*s3*zb(3));
          c10=exp(-alf3*s3*zb(3));
          c11=F(2)/F(3)*alf3*c9;
          c12=F(2)/F(3)*alf3*c10;
          c13=exp(alf3*s3*zb(4));
          c14=exp(-alf3*s3*zb(4));
          c15=exp(alf4*s3*zb(4)); 
          c16=exp(-alf4*s3*zb(4)); 
          c17=F(3)/F(4)*alf4*c15;
          c18=F(3)/F(4)*alf4*c16;
          c19=exp(alf4*s3*zb(5));
          c20=exp(-alf4*s3*zb(5));
          c21=exp(alf5*s3*zb(5));
          c22=exp(-alf5*s3*zb(5));
          c23=F(4)/F(5)*alf5*c21; 
          c24=F(4)/F(5)*alf5*c22;
          c25=exp(alf5*s3*zb(6));
          c26=exp(-alf5*s3*zb(6));
          c27=exp(alf6*s3*zb(6)); 
          c28=exp(-alf6*s3*zb(6));
          c29=F(5)/F(6)*alf6*c27; 
          c30=F(5)/F(6)*alf6*c28;
          c31=exp(alf6*s3*zb(7)); 
          c32=exp(-alf6*s3*zb(7));
          c33=exp(alf7*s3*zb(7));
          c34=exp(-alf7*s3*zb(7));
          c35=F(6)/F(7)*alf7*c33;
          c36=F(6)/F(7)*alf7*c34;
          c37=exp(alf7*s3*zb(8));
          c38=exp(-alf7*s3*zb(8));
                
            
 A =  [1    1    0    0    0    0    0    0    0    0    0    0    0    0;...
      c1    c2  -c3  -c4   0    0    0    0    0    0    0    0    0    0;...
      c1   -c2  -c5   c6   0    0    0    0    0    0    0    0    0    0;...
       0    0    c7   c8  -c9  -c10  0    0    0    0    0    0    0    0;...
       0    0 alf2*c7 -alf2*c8 -c11 c12 0 0    0    0    0    0    0    0 ;...
       0    0    0    0   c13   c14 -c15 -c16  0    0    0    0    0    0;...
       0    0    0    0 alf3*c13 -alf3*c14 -c17 c18 0 0  0    0    0    0;...
       0    0    0    0    0    0   c19  c20 -c21 -c22   0    0    0    0;...
       0    0    0    0    0    0 alf4*c19 -alf4*c20 -c23 c24 0 0  0    0;...
       0    0    0    0    0    0    0    0   c25  c26 -c27 -c28   0    0;...
       0    0    0    0    0    0    0    0 alf5*c25 -alf5*c26 -c29 c30 0 0;...
       0    0    0    0    0    0    0    0    0    0   c31 c32 -c33 -c34;...
       0    0    0    0    0    0    0    0    0    0 alf6*c31 -alf6*c32 -c35 c36;...
       0    0    0    0    0    0    0    0    0    0    0    0   c37  c38];  %BBC: CO3s-C03=0 
                             
              b = [co; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
              x = linsolve(A,b);
              
          
          kcal= kc*PerCaCO3*rhmean/(Mca*CO3LHs(i));     
           
         Tcaldis1=kcal*...
        ((1-fi(1))*(-x(1)/s3*(exp(s3*zb(2))-1)-x(2)/(-s3)*(exp(-s3*zb(2))-1)));
        Tcaldis2 =  kcal*...
        ((1-fi(2))*(-x(3)/(alf2*s3)*(exp(alf2*s3*zb(3))-exp(alf2*s3*zb(2)))-...
         x(4)/(-alf2*s3)*(exp(-alf2*s3*zb(3))-exp(-alf2*s3*zb(2)))));
        Tcaldis3 =  kcal*...
        ((1-fi(3))*(-x(5)/(alf3*s3)*(exp(alf3*s3*zb(4))-exp(alf3*s3*zb(3)))-...
        x(6)/(-alf3*s3)*(exp(-alf3*s3*zb(4))-exp(-alf3*s3*zb(3)))));
        Tcaldis4 =  kcal*...
        ((1-fi(4))*(-x(7)/(alf4*s3)*(exp(alf4*s3*zb(5))-exp(alf4*s3*zb(4)))-...
         x(8)/(-alf4*s3)*(exp(-alf4*s3*zb(5))-exp(-alf4*s3*zb(4))))); 
        Tcaldis5 =  kcal*...
        ((1-fi(5))*(-x(9)/(alf5*s3)*(exp(alf5*s3*zb(6))-exp(alf5*s3*zb(5)))-...
         x(10)/(-alf5*s3)*(exp(-alf5*s3*zb(6))-exp(-alf5*s3*zb(5)))));
        Tcaldis6 =  kcal*...
        ((1-fi(6))*(-x(11)/(alf6*s3)*(exp(alf6*s3*zb(7))-exp(alf6*s3*zb(6)))-...
         x(12)/(-alf6*s3)*(exp(-alf6*s3*zb(7))-exp(-alf6*s3*zb(6)))));
        Tcaldis7 =  kcal*...
        ((1-fi(7))*(-x(13)/(alf7*s3)*(exp(alf7*s3*zb(8))-exp(alf7*s3*zb(7)))-...
         x(14)/(-alf7*s3)*(exp(-alf7*s3*zb(8))-exp(-alf7*s3*zb(7)))));
     
        Tcaldis = Tcaldis1+Tcaldis2+Tcaldis3+Tcaldis4+Tcaldis5+Tcaldis6+Tcaldis7;     
                  
         % Total calcite dissolution rate in the model sediment.
        
    end   
           
 del=R*jo2*(1-fi).*F./Do2;
 w=[w w w w w w w];

[s,x,Rm,C1,D1,Res1]=corg1a7w_eq2(w,Db,jo2,beta,FC,fi,F,del,IniO2,zb,zb(2));
[s,x,Rm,C1,D2,D1,C2,Res2]=corg2a7w_eq2(w,Db,jo2,beta,FC,fi,F,del,IniO2,zb,zb(3));
[s,x,Rm,C1,D3,D2,D1,C2,C3,Res3]=corg3a7w_eq2(w,Db,jo2,beta,FC,fi,F,del,IniO2,zb,zb(4));
[s,x,Rm,C1,D4,D3,D2,D1,C2,C3,C4,Res4]=...
    corg4a7w_eq2(w,Db,jo2,beta,FC,fi,F,del,IniO2,zb,zb(5));
[s,x,Rm,C1,D5,D4,D3,D2,D1,C2,C3,C4,C5,Res5]=...
    corg5a7w_eq2(w,Db,jo2,beta,FC,fi,F,del,IniO2,zb,zb(6));
[s,x,Rm,C1,D6,D5,D4,D3,D2,D1,C2,C3,C4,C5,C6,Res6]=...
      corg6a7w_eq2(w,Db,jo2,beta,FC,fi,F,del,IniO2,zb,zb(7));
[s,x,Rm,C1,D7,D6,D5,D4,D3,D2,D1,C2,C3,C4,C5,C6,C7,Res7]=...
      corg7a7w_eq2(w,Db,jo2,beta,FC,fi,F,del,IniO2,zb,zb(8)); 
if Res1>0
  corg_h1=corg1a7w_eq(w,Db,jo2,beta,FC,fi,F,del,IniO2,zb);
    zo=fzero(corg_h1,[0.001 zb(2)]);
   [s,x,Rm,C1,D1,Res1]=corg1a7w_eq2(w,Db,jo2,beta,FC,fi,F,del,IniO2,zb,zo);
  pt=1;
elseif Res2>0
  corg_h2=corg2a7w_eq(w,Db,jo2,beta,FC,fi,F,del,IniO2,zb);
  zo=fzero(corg_h2,[zb(2) zb(3)]);   
  [s,x,Rm,C1,D2,D1,C2,Res2]=corg2a7w_eq2(w,Db,jo2,beta,FC,fi,F,del,IniO2,zb,zo); 
  pt=2;
elseif Res3>0
  corg_h3=corg3a7w_eq(w,Db,jo2,beta,FC,fi,F,del,IniO2,zb);
  zo=fzero(corg_h3,[zb(3) zb(4)]);
  [s,x,Rm,C1,D3,D2,D1,C2,C3,Res3]=corg3a7w_eq2(w,Db,jo2,beta,FC,fi,F,del,IniO2,zb,zo);
  pt=3;
elseif Res4>0
  corg_h4=corg4a7w_eq(w,Db,jo2,beta,FC,fi,F,del,IniO2,zb);
  zo=fzero(corg_h4,[zb(4) zb(5)]);
  [s,x,Rm,C1,D4,D3,D2,D1,C2,C3,C4,Res4]=...
      corg4a7w_eq2(w,Db,jo2,beta,FC,fi,F,del,IniO2,zb,zo);
  pt=4;
elseif Res5>0
  corg_h5=corg5a7w_eq(w,Db,jo2,beta,FC,fi,F,del,IniO2,zb);
  zo=fzero(corg_h5,[zb(5) zb(6)]);
  [s,x,Rm,C1,D5,D4,D3,D2,D1,C2,C3,C4,C5,Res5]=...
      corg5a7w_eq2(w,Db,jo2,beta,FC,fi,F,del,IniO2,zb,zo); 
  pt=5;
elseif Res6>0
  corg_h6=corg6a7w_eq(w,Db,jo2,beta,FC,fi,F,del,IniO2,zb);
  zo=fzero(corg_h6,[zb(6) zb(7)]);
  [s,x,Rm,C1,D6,D5,D4,D3,D2,D1,C2,C3,C4,C5,C6,Res6]=...
      corg6a7w_eq2(w,Db,jo2,beta,FC,fi,F,del,IniO2,zb,zo); 
  pt=6;
  elseif Res7>0
  corg_h7=corg7a7w_eq(w,Db,jo2,beta,FC,fi,F,del,IniO2,zb);
  zo=fzero(corg_h7,[zb(7) zb(8)]);
  [s,x,Rm,C1,D7,D6,D5,D4,D3,D2,D1,C2,C3,C4,C5,C6,C7,Res7]=...
      corg7a7w_eq2(w,Db,jo2,beta,FC,fi,F,del,IniO2,zb,zo); 
  pt=7;
else  
    zo=zb(8);
    pt=8;
end

Tremin=Rm(1)+Rm(2)+Rm(3)+Rm(4)+Rm(5)+Rm(6)+Rm(7);  


    %Total organic C remineralization in the sediment. This is equal
    %to total source rate to the bottom of the water column of CO2 ion,
    %i.e. one times the source for DIC 
    %
  
         
    if (Fcal-Tcaldis)*Mca +(FC-Tremin)*Moc*compc+NC-...
        (fimeanb-fimean)*zb(8)*rhmean/dtc>=0
    
    w=((Fcal-Tcaldis)*Mca +(FC-Tremin)*Moc*compc+NC-...
        (fimeanb-fimean)*zb(8)*rhmean/dtc)/((1-fi(7))*rhmean);           
    
   
    dP=((Fcal-Tcaldis)*Mca*dtc-(w*(1-fi(7))+zb(8)*(fimeanb-fimean)/dtc)*...
        rhmean*PerCaCO3*dtc)/((1-fimean)*zb(8)*rhmean);
    
       
    dOrg=((FC-Tremin)*Moc*compc*dtc-(w*(1-fi(7))+zb(8)*(fimeanb-fimean)/dtc)*...
        rhmean*PerOrg*dtc)/((1-fimean)*zb(8)*rhmean);    
    
     
    else
        
     w=((Fcal-Tcaldis)*Mca +(FC-Tremin)*Moc*compc+NC-...
        (fimeanb-fimean)*zb(8)*rhmean/dtc)/((1-fiboto)*rhmeano);        
        
     dP=((Fcal-Tcaldis)*Mca*dtc-(w*(1-fiboto)*PerCaCO3o*rhmeano +...
     zb(8)*(fimeanb-fimean)/dtc*PerCaCO3*rhmean)*dtc)/((1-fimean)*zb(8)*rhmean);
 
        
     dOrg=((FC-Tremin)*Moc*compc*dtc-(w*(1-fiboto)*PerOrgo*rhmeano +...
     zb(8)*(fimeanb-fimean)/dtc*PerOrg*rhmean)*dtc)/((1-fimean)*zb(8)*rhmean);
     
    end
 
    PerCaCO3=PerCaCO3+dP;
    PerOrg=PerOrg+dOrg;
    fimeanb = fimean;


  % Output parameters
  %-------------
  % - percentage of calcite dissolution
  pCal(i)   = Tcaldis/Fcal;
  % - percentage organic carbon dissolution
  pOrg(i) =Tremin/FC;
  % - dry weight fraction of calcite in sediment
  dwpCal(i) = PerCaCO3;
  dwpOrg(i) = PerOrg;
  wsed(i) = w;
  fim(i) = fimeanb;
end
  
   
end

           
