function [RcarLL,RcarHL,RorgLL,RorgHL,Wcarb,Wsil,Worg]=...
    ExtForce_M(AT)

% Calculates external inputs to atmosphere and ocean
% Input : 
%         AT    -  Atmospheric tracers, see ODE.m for the data structure.
%         
% Output: 
%         RCarb - River input of carbonate carbon (mol/s)
%         RP - River input of inorganic phosphorus (mol/s)


global fdiv Q10 BCarbPA BCorgPA RPPA Volo swf Htem wfac wfacP wfacC % Get parameter values

ATem = AT(1,1)*sin(fdiv)+AT(1,2)*(1-sin(fdiv));          %Mean atmospheric temperature
       
RcarLL = wfac*0.80*(BCarbPA)*exp(0.1*log(Q10)*(ATem-Htem));          %LL river input of C, Q1o T dependency             
RcarHL = wfac*0.20*(BCarbPA)*exp(0.1*log(Q10)*(ATem-Htem));          %HL river inout of C, Q10 T dependency 

%**********************************************************
% RorgLL_input = 0.8*org_vel(1)*k_Poxidw;
% RorgLL = 0.80*(wfacP*(RPPA)*exp(0.1*log(Q10)*(ATem-Htem)) + org_vel(1)*k_Poxidw);           %Adding org_input
% RorgHL = 0.20*(wfacP*(RPPA)*exp(0.1*log(Q10)*(ATem-Htem)) + org_vel(1)*k_Poxidw);

RorgLL = 0.80*(wfacP*(RPPA)*exp(0.1*log(Q10)*(ATem-Htem)));           %Adding org_input
RorgHL = 0.20*(wfacP*(RPPA)*exp(0.1*log(Q10)*(ATem-Htem)));

%*********************************************************

Wcarb = wfac*(BCarbPA)/(1+swf)*exp(0.1*log(Q10)*(ATem-Htem)); 
Wsil  = swf*Wcarb; 
Worg  = wfacC*((BCorgPA)+swf*(BCarbPA)/(1+swf)-Volo)*exp(0.1*log(Q10)*(ATem-Htem));%
return
