clc; clear; close all;

% Load the image
img = imread('Figure\Methane release rate_1.png'); % Replace with your image file
grayImg = rgb2gray(img);  % Convert to grayscale
bwImg = imbinarize(grayImg, 'adaptive'); % Convert to binary (black & white)

% Extract edge points
[y, x] = find(bwImg); % Get (x, y) coordinates of the curve

% Convert pixel coordinates to meaningful values
x = -max(x) + x; % Flip to match the image's x-axis
y = max(y) - y;  % Flip the y-axis if needed

% Sort x-values for proper plotting
[x, idx] = sort(x);
y = y(idx);

% Plot extracted points
figure;
imshow(img);
hold on;
scatter(x, y, 10, 'r', 'filled'); % Show extracted points
title('Extracted Curve from Image');
xlabel('X');
ylabel('Y');
grid on;

figure
plot(x, y)


% % Fit a Gamma distribution
% a0 = 5; % Initial shape parameter
% b0 = 1; % Initial scale parameter
% params = mle(y, 'distribution', 'gamma'); % Maximum likelihood estimation
% 
% % Generate fitted curve
% xx = linspace(min(x), max(x), 1000);
% yy_fit = gampdf(xx, params(1), params(2));
% 
% % Plot the fitted curve
% figure;
% plot(x, y, 'ro', 'MarkerSize', 5, 'DisplayName', 'Extracted Data');
% hold on;
% plot(xx, yy_fit, 'b-', 'LineWidth', 2, 'DisplayName', 'Gamma Fit');
% legend;
% xlabel('X');
% ylabel('Y');
% title('Gamma Distribution Fit');
% grid on;
