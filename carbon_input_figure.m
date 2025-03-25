function C_input_rate = carbon_input_figure(C_total, start_time, duration, pathway)

% Read and process the image
img = imread(pathway);  % Load image

% Convert to grayscale and binarize
grayImg = rgb2gray(img); 
bwImg = imbinarize(grayImg, 'adaptive'); % Adaptive thresholding

% Extract edge points
[y, x] = find(bwImg); % Get (x, y) coordinates

% Convert pixel coordinates to meaningful values
x = -max(x) + x; % Flip x-axis
y = max(y) - y;  % Flip y-axis

% Sort x-values for proper plotting
[x, idx] = sort(x);
y = y(idx);

% Normalize x-axis and y-aixs
xNorm = (x - min(x)) / (max(x) - min(x)); 
yNorm = (y-min(y))/ (max(y)-min(y)); 

% Scale x-axis to match physical units
TC = start_time + xNorm * duration; 

% Apply LOESS smoothing
span = max(0.05, min(0.1, 10/length(TC))); % Auto-adjust span
ySmooth = smooth(TC, yNorm, span, 'loess');

% Calculate area under the curve using trapezoidal rule
area_trapz = trapz(TC, ySmooth);

% Normalize ySmooth to match a reference area
RC = (C_total / area_trapz) * ySmooth;

C_input_rate = [TC, RC];

% % Plot results
% figure("Position",[0,0,500,200])
% plot(TC, RC, 'LineWidth', 1); hold on;
% xlabel('Time (yr)');
% ylabel("Carbon input (Gtyr^{-1})");
% set(gca,"linewidth", 0.5,"FontSize",10,"FontName", "Times New Roman","TickLength",[0.01,0.012],"Layer","top");
% print(gcf,"Figure\Methane volume","-dpng","-r600");

end

