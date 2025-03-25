function c_flux = enforcing(t,x,y) 
% Preallocate arrays for slope and intercept
x = x(~isnan(y));
y = y(~isnan(y));
num_segments = length(x) - 1;
slopes = zeros(1, num_segments);
intercepts = zeros(1, num_segments);

% Calculate slopes and intercepts for each segment
for i = 1:num_segments
    slopes(i) = (y(i+1) - y(i)) / (x(i+1) - x(i));
    intercepts(i) = y(i) - slopes(i) * x(i);
end
%% 
% Function to evaluate piecewise linear equation
evaluate_piecewise = @(x_val) piecewise_eval(x_val, x, slopes, intercepts);

% Example: Evaluate for a set of x values
c_flux = arrayfun(evaluate_piecewise, t); %unit: kyr


% Piecewise evaluation function
function y = piecewise_eval(x_val, x_points, slopes, intercepts)
    % Find the segment corresponding to x_val
    for i = 1:length(slopes)
        if x_val >= x_points(i) && x_val <= x_points(i+1)
            % Apply the linear equation for the current segment
            y = slopes(i) * x_val + intercepts(i);
            return;
        end
    end
    % Handle out-of-bound x values (optional)
    error('Input x value %.2f is outside the range of defined segments.', x_val);
end
end
