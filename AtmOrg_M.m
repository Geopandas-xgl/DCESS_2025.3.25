function org_vel = AtmOrg_M(t,input)

global sy mgt R13pdb orgc13

input_end = input.org.sta + input.org.dur;

if (t/sy>=input.org.sta && t/sy<= input_end)
    [~,minIndex] = min(abs(input.org.rate(:,1)-t/sy));
    rate  = input.org.rate(minIndex,2)*mgt/sy;
else
    rate = 0;
end

org_vel(1) = rate;  %
org_vel(2) = rate*(orgc13*1e-3+1)*R13pdb;
end