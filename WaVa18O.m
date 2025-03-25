function [AO18c]=WaVa18O(Tw,Td,Ow)
% Water Vapour 18O

h      = .75;
k18    = .006;
alp18w = (.9884 + 1.025e-4*(Tw)   - 3.57e-7*(Tw).^2)^-1;
alp18d = (.9884 + 1.025e-4*(Td)   - 3.57e-7*(Td).^2)^-1;
alp18c = .5*(alp18w+alp18d);
eww    = 10^(9.4-2.35e3/(Tw+273));
ewd    = 10^(9.4-2.35e3/(Td+273));
Fv     = ewd/eww;
AO18w  = ( 1/alp18w * (Ow+1)*(1-k18)/(1-h*k18)-1 );
AO18c  = Fv^(alp18c-1)*(1+AO18w)-1;

return