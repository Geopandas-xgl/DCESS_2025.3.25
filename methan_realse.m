clear all
% Parameters
a = 19;      % Transition point 17.65
b = 0.28;      % Power function exponent
k = 25;      % Gamma shape parameter
theta = 0.77;  % Gamma scale parameter
t = 13;

% Define the power function for x < a

power_func = @(x) 0.25*(1 - 1 ./ (1 + exp(-b * (x - 25))));
% Compute the matching value at x = a
power_value = power_func(a);

% Define the Gamma distribution function for x >= a, scaled to match continuity
gamma_func = @(x) 0.6*(gampdf(x, k, theta)+0.01);

gamma_value = gamma_func(a);

gamma_value_t = gamma_func(t);

line_func = @(x) gamma_value_t;

% Define the domain
x = 0:0.01:23;

% Define the piecewise function
f = @(x) (t>x).*line_func(x)+ (t <=x & x < a) .* (gamma_func(x)) + (x >= a).* (power_func(x)-power_value + gamma_value);

f1 = @(x) 6660*(f(x)- f(max(x)));


% Evaluate the function
y = f1(x) ;

S = integral(f1, min(x), max(x));

y1 = power_func(x)-power_value + gamma_value;
% Plot the function
figure
% plot(x, y1, 'LineWidth', 2); hold on
plot(x, y, 'LineWidth', 2);
xlim([0,23]);
% ylim([0,0.12]);
xlabel('x');
% xticks(3:5:23);
% xticklabels(-20:5:0);
ylabel('f(x)');
grid on;
