% Subroutine for calculating zo, the level in the sediment where dissolved
%oxygen concentration goes to zero, anoxix remineralization included, 
%including the oxidation at zo of H2S and NH4 producted in the anoxic layer
%below zo.


function [s,x,Rm,C1,D1,Res1]=corg1a7w_eq2(w,Db,jo2,beta,FC,fi,F,del,IniO2,zb,zo)

 
  s(1) = 0.5*(w(1)/Db + (w(1)^2/(Db^2)+4*jo2/Db)^0.5);
  s(2) = 0.5*(w(1)/Db - (w(1)^2/(Db^2)+4*jo2/Db)^0.5);
  s(3) = 0.5*(w(2)/Db + (w(2)^2/(Db^2)+4*jo2/Db)^0.5);
  s(4) = 0.5*(w(2)/Db - (w(2)^2/(Db^2)+4*jo2/Db)^0.5);
  s(5) = 0.5*(w(3)/Db + (w(3)^2/(Db^2)+4*jo2/Db)^0.5);
  s(6) = 0.5*(w(3)/Db - (w(3)^2/(Db^2)+4*jo2/Db)^0.5);
  s(7) = 0.5*(w(4)/Db + (w(4)^2/(Db^2)+4*jo2/Db)^0.5);
  s(8) = 0.5*(w(4)/Db - (w(4)^2/(Db^2)+4*jo2/Db)^0.5);
  s(9) = 0.5*(w(5)/Db + (w(5)^2/(Db^2)+4*jo2/Db)^0.5);
  s(10) = 0.5*(w(5)/Db - (w(5)^2/(Db^2)+4*jo2/Db)^0.5);

  s(11) = 0.5*(w(1)/Db + (w(1)^2/(Db^2)+4*beta*jo2/Db)^0.5);
  s(12) = 0.5*(w(1)/Db - (w(1)^2/(Db^2)+4*beta*jo2/Db)^0.5);
  s(13) = 0.5*(w(2)/Db + (w(2)^2/(Db^2)+4*beta*jo2/Db)^0.5);
  s(14) = 0.5*(w(2)/Db - (w(2)^2/(Db^2)+4*beta*jo2/Db)^0.5);
  s(15) = 0.5*(w(3)/Db + (w(3)^2/(Db^2)+4*beta*jo2/Db)^0.5);
  s(16) = 0.5*(w(3)/Db - (w(3)^2/(Db^2)+4*beta*jo2/Db)^0.5);
  s(17) = 0.5*(w(4)/Db + (w(4)^2/(Db^2)+4*beta*jo2/Db)^0.5);
  s(18) = 0.5*(w(4)/Db - (w(4)^2/(Db^2)+4*beta*jo2/Db)^0.5);
  s(19) = 0.5*(w(5)/Db + (w(5)^2/(Db^2)+4*beta*jo2/Db)^0.5);
  s(20) = 0.5*(w(5)/Db - (w(5)^2/(Db^2)+4*beta*jo2/Db)^0.5);
  s(21) = 0.5*(w(6)/Db + (w(6)^2/(Db^2)+4*beta*jo2/Db)^0.5);
  s(22) = 0.5*(w(6)/Db - (w(6)^2/(Db^2)+4*beta*jo2/Db)^0.5);
  s(23) = 0.5*(w(7)/Db + (w(7)^2/(Db^2)+4*beta*jo2/Db)^0.5);
  s(24) = 0.5*(w(7)/Db - (w(7)^2/(Db^2)+4*beta*jo2/Db)^0.5);
  
 
  c1=(w(1)-Db*s(1))*(1-fi(1));
  c2=(w(1)-Db*s(2))*(1-fi(1));
  c3=exp(s(1)*zo);
  c4=exp(s(2)*zo);
  c5=-exp(s(11)*zo);
  c6=-exp(s(12)*zo);
  c7=-(w(1)-Db*s(1))*exp(s(1)*zo);
  c8=-(w(1)-Db*s(2))*exp(s(2)*zo);
  c9=(w(1)-Db*s(11))*exp(s(11)*zo);
  c10=(w(1)-Db*s(12))*exp(s(12)*zo);
  c11=exp(s(11)*zb(2));
  c12=exp(s(12)*zb(2));
  c13=-exp(s(13)*zb(2));
  c14=-exp(s(14)*zb(2));
  %c15=-(Db*s(11))*exp(s(11)*zb(2))*(1-fi(1));
  %c16=-(Db*s(12))*exp(s(12)*zb(2))*(1-fi(1));
  %c17=(Db*s(13))*exp(s(13)*zb(2))*(1-fi(2));
  %c18=(Db*s(14))*exp(s(14)*zb(2))*(1-fi(2));
  c15=-(w(1)-Db*s(11))*exp(s(11)*zb(2))*(1-fi(1));
  c16=-(w(1)-Db*s(12))*exp(s(12)*zb(2))*(1-fi(1));
  c17=(w(2)-Db*s(13))*exp(s(13)*zb(2))*(1-fi(2));
  c18=(w(2)-Db*s(14))*exp(s(14)*zb(2))*(1-fi(2));
  c19=exp(s(13)*zb(3));
  c20=exp(s(14)*zb(3));
  c21=-exp(s(15)*zb(3));
  c22=-exp(s(16)*zb(3));
  c23=-(w(2)-Db*s(13))*exp(s(13)*zb(3))*(1-fi(2));
  c24=-(w(2)-Db*s(14))*exp(s(14)*zb(3))*(1-fi(2));
  c25=(w(3)-Db*s(15))*exp(s(15)*zb(3))*(1-fi(3));
  c26=(w(3)-Db*s(16))*exp(s(16)*zb(3))*(1-fi(3));
  c27=exp(s(15)*zb(4));
  c28=exp(s(16)*zb(4));
  c29=-exp(s(17)*zb(4));
  c30=-exp(s(18)*zb(4));
  c31=-(w(3)-Db*s(15))*exp(s(15)*zb(4))*(1-fi(3));
  c32=-(w(3)-Db*s(16))*exp(s(16)*zb(4))*(1-fi(3));
  c33=(w(4)-Db*s(17))*exp(s(17)*zb(4))*(1-fi(4));
  c34=(w(4)-Db*s(18))*exp(s(18)*zb(4))*(1-fi(4));
  c35=exp(s(17)*zb(5));
  c36=exp(s(18)*zb(5));
  c37=-exp(s(19)*zb(5));
  c38=-exp(s(20)*zb(5));
  c39=-(w(4)-Db*s(17))*exp(s(17)*zb(5))*(1-fi(4));
  c40=-(w(4)-Db*s(18))*exp(s(18)*zb(5))*(1-fi(4));
  c41=(w(5)-Db*s(19))*exp(s(19)*zb(5))*(1-fi(5));
  c42=(w(5)-Db*s(20))*exp(s(20)*zb(5))*(1-fi(5));
  c43=exp(s(19)*zb(6));
  c44=exp(s(20)*zb(6));
  c45=-exp(s(21)*zb(6));
  c46=-exp(s(22)*zb(6));
  c47=-(w(5)-Db*s(19))*exp(s(19)*zb(6))*(1-fi(5));
  c48=-(w(5)-Db*s(20))*exp(s(20)*zb(6))*(1-fi(5));
  c49=(w(6)-Db*s(21))*exp(s(21)*zb(6))*(1-fi(6));
  c50=(w(6)-Db*s(20))*exp(s(20)*zb(6))*(1-fi(6));
  c51=exp(s(21)*zb(7));
  c52=exp(s(22)*zb(7));
  c53=-exp(s(23)*zb(7));
  c54=-exp(s(24)*zb(7));
  c55=-(w(6)-Db*s(21))*exp(s(21)*zb(7))*(1-fi(6));
  c56=-(w(6)-Db*s(22))*exp(s(22)*zb(7))*(1-fi(6));
  c57=(w(7)-Db*s(23))*exp(s(23)*zb(7))*(1-fi(7));
  c58=(w(7)-Db*s(24))*exp(s(24)*zb(7))*(1-fi(7));
  c59=(Db*s(23))*exp(s(23)*zb(8));
  c60=(Db*s(24))*exp(s(24)*zb(8));
  
 
  
  A =[  c1  c2  0   0   0   0   0   0   0   0   0   0   0   0   0   0;...
        c3  c4  c5  c6  0   0   0   0   0   0   0   0   0   0   0   0;...
        c7  c8  c9  c10 0   0   0   0   0   0   0   0   0   0   0   0;...
        0   0   c11 c12 c13 c14 0   0   0   0   0   0   0   0   0   0;...
        0   0   c15 c16 c17 c18 0   0   0   0   0   0   0   0   0   0;...
        0   0   0   0   c19 c20 c21 c22 0   0   0   0   0   0   0   0;...
        0   0   0   0   c23 c24 c25 c26 0   0   0   0   0   0   0   0;...
        0   0   0   0   0   0   c27 c28 c29 c30 0   0   0   0   0   0;...
        0   0   0   0   0   0   c31 c32 c33 c34 0   0   0   0   0   0;...
        0   0   0   0   0   0   0   0   c35 c36 c37 c38 0   0   0   0;...
        0   0   0   0   0   0   0   0   c39 c40 c41 c42 0   0   0   0;...
        0   0   0   0   0   0   0   0   0   0   c43 c44 c45 c46 0   0;...
        0   0   0   0   0   0   0   0   0   0   c47 c48 c49 c50 0   0;...
        0   0   0   0   0   0   0   0   0   0   0   0   c51 c52 c53 c54;...
        0   0   0   0   0   0   0   0   0   0   0   0   c55 c56 c57 c58;...
        0   0   0   0   0   0   0   0   0   0   0   0   0   0   c59 c60];...
        
  
 B = [FC; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];         
           
 x = linsolve(A,B);
 
 Rm1=(1-fi(1))*jo2*(x(1)*(exp(s(1)*zo)-1)/s(1)+...
                    x(2)*(exp(s(2)*zo)-1)/s(2)+...
              beta*(x(3)*(exp(s(11)*zb(2))-exp(s(11)*zo))/s(11)+...
                    x(4)*(exp(s(12)*zb(2))-exp(s(12)*zo))/s(12)));
 Rm2=(1-fi(2))*beta*jo2*(x(5)*(exp(s(13)*zb(3))-exp(s(13)*zb(2)))/s(13)...
                        +x(6)*(exp(s(14)*zb(3))-exp(s(14)*zb(2)))/s(14));                   
 Rm3=(1-fi(3))*beta*jo2*(x(7)*(exp(s(15)*zb(4))-exp(s(15)*zb(3)))/s(15)...
                        +x(8)*(exp(s(16)*zb(4))-exp(s(16)*zb(3)))/s(16));
 Rm4=(1-fi(4))*beta*jo2*(x(9)*(exp(s(17)*zb(5))-exp(s(17)*zb(4)))/s(17)...
                        +x(10)*(exp(s(18)*zb(5))-exp(s(18)*zb(4)))/s(18));
 Rm5=(1-fi(5))*beta*jo2*(x(11)*(exp(s(19)*zb(6))-exp(s(19)*zb(5)))/s(19)+...
                         x(12)*(exp(s(20)*zb(6))-exp(s(20)*zb(5)))/s(20));
 Rm6=(1-fi(6))*beta*jo2*(x(13)*(exp(s(21)*zb(7))-exp(s(21)*zb(6)))/s(21)+...
                         x(14)*(exp(s(22)*zb(7))-exp(s(22)*zb(6)))/s(22));
 Rm7=(1-fi(7))*beta*jo2*(x(15)*(exp(s(23)*zb(8))-exp(s(23)*zb(7)))/s(23)+...
                         x(16)*(exp(s(24)*zb(8))-exp(s(24)*zb(7)))/s(24));
 Rm8= (1-fi(1))*jo2*(beta*(x(3)*(exp(s(11)*zb(2))-exp(s(11)*zo))/s(11)+...
                    x(4)*(exp(s(12)*zb(2))-exp(s(12)*zo))/s(12)))+Rm2+Rm3+Rm4+Rm5+Rm6+Rm7;

 Rm = [Rm1 Rm2 Rm3 Rm4 Rm5 Rm6 Rm7 Rm8];
 
 
  C1 = IniO2 - del(1)*(x(1)/(s(1)^2)+x(2)/(s(2)^2));
  
  D1 = -del(1)*((x(1)/s(1)*exp(s(1)*zo)+x(2)/s(2)*exp(s(2)*zo)...
  +beta*(x(3)*(exp(s(11)*zb(2))-exp(s(11)*zo))/s(11)+x(4)*(exp(s(12)*zb(2))-exp(s(12)*zo))/s(12)...
    +del(2)/del(1)*(x(5)*(exp(s(13)*zb(3))-exp(s(13)*zb(2)))/s(13)+x(6)*(exp(s(14)*zb(3))-exp(s(14)*zb(2)))/s(14))...
    +del(3)/del(1)*(x(7)*(exp(s(15)*zb(4))-exp(s(15)*zb(3)))/s(15)+x(8)*(exp(s(16)*zb(4))-exp(s(16)*zb(3)))/s(16))...
    +del(4)/del(1)*(x(9)*(exp(s(17)*zb(5))-exp(s(17)*zb(4)))/s(17)+x(10)*(exp(s(18)*zb(5))-exp(s(18)*zb(4)))/s(18))...
    +del(5)/del(1)*(x(11)*(exp(s(19)*zb(6))-exp(s(19)*zb(5)))/s(19)+x(12)*(exp(s(20)*zb(6))-exp(s(20)*zb(5)))/s(20))...
    +del(6)/del(1)*(x(13)*(exp(s(21)*zb(7))-exp(s(21)*zb(6)))/s(21)+x(14)*(exp(s(22)*zb(7))-exp(s(22)*zb(6)))/s(22))...
    +del(7)/del(1)*(x(15)*(exp(s(23)*zb(8))-exp(s(23)*zb(7)))/s(23)+x(16)*(exp(s(24)*zb(8))-exp(s(24)*zb(7)))/s(24)))));
    
    %D1 = - del(1)*(x(1)/s(1)*(exp(s(1)*zo))+x(2)/s(2)*(exp(s(2)*zo)));
  
  Res1 = -(C1 + D1*zo + del(1)*(x(1)*exp(s(1)*zo)./(s(1)^2)+x(2)*exp(s(2)*zo)/(s(2)^2)));

end
