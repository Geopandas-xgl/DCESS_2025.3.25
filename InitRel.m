function [A,B] = InitRel(tot,frac,tfrac)

% Set release parameters
% MH or CO2 release has the form A*t^n*exp(-B*t), n = 4
% The parameter B is found so that a set fraction of the total release of
% MH (frac)occurs before a given time (tfrac) after the onset of release (tCH4).

% Initial guess for B 
B=4./tfrac;

fA=frac*24;
pol=24-((tfrac*B)^4+4*(tfrac*B)^3+12*(tfrac*B)^2+ ...
    24*tfrac*B+24)*exp(-B*tfrac);
delta=0.0000010;

while(abs(fA-pol) > 0.1)

if (fA<pol)
    B=B-delta;
else
    B=B+delta;
end
pol=24-((tfrac*B)^4+4*(tfrac*B)^3+12*(tfrac*B)^2+...
    24*tfrac*B+24)*exp(-B*tfrac);
end

delta=delta/10.;
while(abs(fA-pol) > 0.01)

if (fA<pol)
    B=B-delta;
else
    B=B+delta;
end
pol=24-((tfrac*B)^4+4*(tfrac*B)^3+12*(tfrac*B)^2+...
    24*tfrac*B+24)*exp(-B*tfrac);
end

A=tot/24*B^5;

return

