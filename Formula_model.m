%% exponential function
figure
x = 0:0.1:2;

for k = 0:0.1:1
    plot(x,exp(k.*(x-1))); hold on
end
