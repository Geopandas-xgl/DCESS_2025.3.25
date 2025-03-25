function lip_vel = AtmLip_M(t,input)

global sy mgt R13pdb lipc13

input_end = input.lip.sta + input.lip.dur;

if (t/sy>=input.lip.sta && t/sy<= input_end)
    [~,minIndex] = min(abs(input.lip.rate(:,1)-t/sy));
    rate  = input.lip.rate(minIndex,2)*mgt/sy;
else
    rate = 0;
end

lip_vel(1) = rate;  %
lip_vel(2) = rate*(lipc13*1e-3+1)*R13pdb;
end