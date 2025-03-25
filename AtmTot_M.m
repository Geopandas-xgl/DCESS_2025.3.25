function tot_vel = AtmTot_M(input)

input.meth.end = input.meth.sta + input.meth.dur;
input.mett.end = input.mett.sta + input.mett.dur;
input.org.end = input.org.sta + input.org.dur;
input.lip.end = input.lip.sta + input.lip.dur;

t = 0:10:max(input.lip.end);
num = length(t);
meth_rate = nan(num,1);
mett_rate = nan(num,1);
org_rate = nan(num,1);
lip_rate = nan(num,1);

for i = 1:num
    if (t(i)>=input.meth.sta && t(i)<= input.meth.end)
        [~,minIndex] = min(abs(input.meth.rate(:,1)-t(i)));
        meth_rate(i)  = input.meth.rate(minIndex,2);
    else
        meth_rate(i) = 0;
    end

     if (t(i)>=input.mett.sta && t(i)<= input.mett.end)
        [~,minIndex] = min(abs(input.mett.rate(:,1)-t(i)));
        mett_rate(i)  = input.mett.rate(minIndex,2);
    else
        mett_rate(i) = 0;
    end
    
    if (t(i)>=input.org.sta && t(i)<= input.org.end)
        [~,minIndex] = min(abs(input.org.rate(:,1)-t(i)));
        org_rate(i)  = input.org.rate(minIndex,2);
    else
        org_rate(i) = 0;
    end
    
    if (t(i)>=input.lip.sta && t(i)<= input.lip.end)
        [~,minIndex] = min(abs(input.lip.rate(:,1)-t(i)));
        lip_rate(i)  = input.lip.rate(minIndex,2);
    else
        lip_rate(i) = 0;
    end
end

t = t';
rate = meth_rate + mett_rate + org_rate + lip_rate;

tot_vel = [t, rate];

end


