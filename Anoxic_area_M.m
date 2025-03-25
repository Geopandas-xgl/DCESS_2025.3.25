clear 

global dm d n aHL aLL

load('OutThilda_M')
load ResThilda_M.mat
GAw = load('GAw_Eo.txt');
GAc = load('GAc_Eo.txt');

% Get parameter values
ParVal_M    % Activate global parameters
%%
fHL = aHL/(aHL+aLL);
fLL = aLL/(aHL+aLL);
% Define the range for x and y
zcent= [dm/2 dm+(d:d:(n-1)*d)-d/2]/1e3;   % Vertical center of boxes
x = st/1e3;
y = zcent;
% Create the grid
[X, Y] = meshgrid(x, y);
% Compute the function values
ZLL = permute(sLL(8,:,:),[2,3,1])*1000;
ZHL = permute(sHL(8,:,:),[2,3,1])*1000;

RS_LL = nan(55,1);
RS_HL =  nan(55,1);

RS_LL(1) = 0;
RS_HL(1) = 0;

for i = 1:54
   RS_LL(i+1) = (GAw(i)-GAw(i+1))/GAw(1); % obtain seafloor area of each layer
   RS_HL(i+1) = (GAc(i)-GAc(i+1))/GAc(1);
end

RS_LL_inv = RS_LL';
RS_HL_inv = RS_HL';

LLO_init = LL(8,:)*1000; %mmol/m3
HLO_init = HL(8,:)*1000; %mmol/m3LHO_init

fanox_init = 0.002;
crit = 8.93; % mmol/m3 equal to anoxic criteria of 0.2 ml/mol

syms k

eqn = fLL*sum(exp(-k.*(LLO_init-crit)).*RS_LL_inv) + fHL*sum(exp(-k.*(HLO_init-crit)).*RS_HL_inv) == fanox_init;

k = double(vpasolve(eqn, k));

RS_LL_matr = repmat(RS_LL,1,length(x));
RS_HL_matr = repmat(RS_HL,1,length(x));

ZLL(ZLL<=crit) = crit;
deg_LL = exp(-k.*(ZLL-crit));

ZHL(ZHL<=crit) = crit;
deg_HL = exp(-k.*(ZHL-crit));

fanox_LL = sum(deg_LL.*RS_LL_matr)*100;
fanox_HL = sum(deg_HL.*RS_HL_matr)*100;
fanox_global = fLL*fanox_LL + fHL*fanox_HL;

fanox_LL_init = sum(exp(-k.*(LLO_init-crit)).*RS_LL_inv);
fanox_HL_init = sum(exp(-k.*(HLO_init-crit)).*RS_HL_inv);
%%
figure("Position",[0,0,800,300])

[ha,~] = tight_subplot(1,2,[0,0.06],[0.16,0.08],[0.06,0.01]);

axes(ha(1))
[c, h] = contour(X, Y, ZLL);
clabel(c, h, 'FontSize', 6)
colorbar;
xlabel('Time (kyr)');
ylabel('Depth (km)');
title("Low latitude");
set(gca,"YDir","reverse","linewidth", 1,"FontSize",10,"FontName", "Times","TickLength",[0.02,0.03],"Layer","top");


axes(ha(2))
[c, h] = contour(X, Y, ZHL);
clabel(c, h, 'FontSize', 6)
colorbar;
xlabel('Time (kyr)');
ylabel('Depth (km)');
title("High Latitude");
set(gca,"YDir","reverse","linewidth", 1,"FontSize",10,"FontName", "Times","TickLength",[0.02,0.03],"Layer","top");
print(gcf,"Figure\Oxygen concentration_Contour");
%% 

figure("Position",[0,0,400,180])
[ha,~] = tight_subplot(1,1,[0,0],[0.22,0.05],[0.1,0.01]);
axes(ha(1))
plot(st/1e3,fanox_LL, "LineWidth",1.4); hold on
plot(st/1e3,fanox_HL, "LineWidth",1.4);
plot(st/1e3,fanox_global, "LineWidth",1.4); 
xlabel("Time (kyr)");
ylabel("Area fraction (%)");
legend(["Low latitude", "High latitude","Global"],"Location","northwest","FontSize",8)
print(gcf,"Figure\anoxic fraction","-dpng","-r600");
set(gca,"linewidth", 0.8,"FontSize",10,"FontName", "Times","TickLength",[0.01,0.015],"Layer","top");


% axes(ha(4))
% plot(st/1e3,fanox_LL/fanox_LL_init, "LineWidth",1.2); hold on
% plot(st/1e3,fanox_HL/fanox_HL_init, "LineWidth",1.2);
% plot(st/1e3,fanox_global/fanox_init, "LineWidth",1.2); 
% xlabel("Time (kyr)");
% ylabel("Anoxic expansion times");
%%
% 
% figure("Position",[0,0,1000,300])
% [h,~]= tight_subplot(1,3,[0,0.08],[0.15,0.1],[0.08,0.01]);
% axes(h(1))
% for i = 20:20:120
% area_layer_matrix_1 = area_layer_matrix;
% area_layer_matrix_1(ZLL>=i) = 0; % selected data under the threshold of O2 concentration; % 0.2 ml/mol: 8.93 mmol/m3; 0.5 ml/mol: 22.3 mmol/m3
% area_anox_1 = sum(area_layer_matrix_1)*100; % calclated total seafloor anoxic area in specifield time;
% plot(st/1e3,area_anox_1, "LineWidth",1.2); hold on
% xlabel("Time (kyr)");
% ylabel("Area fraction (%)");
% end
% lgd = legend(h(1),string(20:20:120));
% lgd.Title.String = "O_{2}";
% title("Figure 2a");
% set(gca,"linewidth", 1,"FontSize",12,"FontName", "Times New Roman","TickLength",[0.02,0.025],"Layer","top");
% 
% 
% axes(h(2))
% for i = 20
% area_layer_matrix_1 = area_layer_matrix;
% area_layer_matrix_1(ZLL>=i) = 0; % 0.2 ml/mol: 8.93 mmol/m3; 0.5 ml/mol: 22.3 mmol/m3
% area_anox_1 = sum(area_layer_matrix_1)*100;
% plot(st/1e3,area_anox_1, "LineWidth",1.2); hold on
% xlabel("Time (kyr)");
% ylabel("Area fraction (%)");
% end
% lgd = legend(h(2),string(20));
% lgd.Title.String = "O_{2}";
% title("Figure 2b");
% set(gca,"linewidth", 1,"FontSize",12,"FontName", "Times New Roman","TickLength",[0.02,0.025],"Layer","top");
% 
% axes(h(3))
% oxy = 60:20:80;
% clr = ["#F2C965","#7E2F8E"];
% for i = 1: length(oxy)
% oxy_i = oxy(i); 
% area_layer_matrix_1 = area_layer_matrix;
% area_layer_matrix_1(ZLL>=oxy_i) = 0; % 0.2 ml/mol: 8.93 mmol/m3; 0.5 ml/mol: 22.3 mmol/m3
% area_anox_1 = sum(area_layer_matrix_1)*100;
% fold_anox_1 = area_anox_1/area_anox_1(1);
% plot(st/1e3,fold_anox_1,"Color",clr(i),"LineWidth",1.2); hold on
% xlabel("Time (kyr)");
% ylabel("expanding folds");
% end
% legend(h(3),string(60:20:80));
% title("Figure 2c");
% set(gca,"linewidth", 1,"FontSize",12,"FontName", "Times New Roman","TickLength",[0.02,0.025],"Layer","top");
% 
% 

% oxy_layer_matrix = Z.*area_layer_matrix;
% oxy_sum = sum(oxy_layer_matrix(3:30,:));
% oxy_norm = oxy_sum/oxy_sum(1); % In the modern seawater, the fraction of  seafloor anoxia (<0.2 ml/mol) is 0.2% 
% anox_area = 1./oxy_norm;
%% 


% %%
% figure
% plot(con_o,danox)
% xlabel("O_2 concentration (mmol/m^{3})")
% ylabel("anoxic degree (%)")
% box off
% BorderLine(gca,1);
% set(gca,"linewidth", 1,"FontSize",10,"FontName", "Times","TickLength",[0.02,0.03],"Layer","top");


