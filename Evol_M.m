clear
% % function Evol_M
global sy rcp n fdiv R13pdb rpm mgt methc13 mettc13 orgc13 lipc13 dm d n aHL aLL pOrgLL pOrgHL

load('OutThilda_M')
load ResThilda_M.mat

warning off

ParVal_M

GAw = load('GAw_Eo.txt');
GAc = load('GAc_Eo.txt');

% data1 = readtable("C:\Users\19328\Nutstore\1\我的坚果云\PETM\Data\Xiong-2023-PETM_1.xlsx","Sheet","C and O isotope");
% data2 = readtable("C:\Users\19328\Nutstore\1\我的坚果云\PETM\Data\Xiong-2023-PETM_1.xlsx","Sheet","Leach and disolution");
% data3 = readtable("C:\Users\19328\Nutstore\1\我的坚果云\PETM\Data\recovered_trends.csv");
% color = readtable("C:\Users\19328\Documents\Corel\Corel Content\Palettes\Isson-2018-Science.xml");

data1 = readtable("data\Xiong-2025-PETM.xlsx","Sheet","C and O isotope");
data2 = readtable("data\Xiong-2025-PETM.xlsx","Sheet","U and Ba isotope");
data3 = readtable("data\Xiong-2025-PETM.xlsx","Sheet","Anoxic area");
color = readtable("data\Xiong-2025-PETM.xlsx","Sheet","Color");

data.dC1 = data1.x_13C(1:end-3);
data.age1 = data1.age_westerhold(1:end-3);

data.Ba = data2(1:end-3, ["Age_westerhold","x_138Ba","x2SD_1"]);
data.Ba = data.Ba(~isnan(data.Ba.x_138Ba),:);
data.Ba.Age_westerhold = data.Ba.Age_westerhold + 22.3;

data.age3 = data3.time/1000;
data.dU3 = data3(:,3:102);
data.fanox3 = data3(:,203:end);
numericCols3 = varfun(@isnumeric,data.fanox3,"OutputFormat","uniform");
data.fanox3{:,numericCols3} = 100*10.^data.fanox3{:,numericCols3};

data.dU_age3 = [[data.age3,data.dU3.d238U_0_16];[flip(data.age3),flip(data.dU3.d238U_0_84)]];
data.fanox_age3 = [[data.age3,data.fanox3.fanox_0_16];[flip(data.age3),flip(data.fanox3.fanox_0_84)]];
data.num_idx3 = 1:size(data.dU_age3,1);


% data.dC1 = data1.dC(1:end-3);
% data.age1 = data1.age_westerhold(1:end-3)+22.3;
% data.Ba = data2(5:end-3, ["Age_westerhold","d138Ba_cor","dBa2SD"]);
% data.Ba = data.Ba(~isnan(data.Ba.d138Ba_cor),:);
% data.Ba.Age_westerhold = data.Ba.Age_westerhold + 22.3;
% 
% data.age3 = data3.time/1000;
% data.dU3 = data3(:,3:102);
% data.fanox3 = data3(:,203:end);
% numericCols3 = varfun(@isnumeric,data.fanox3,"OutputFormat","uniform");
% data.fanox3{:,numericCols3} = 100*10.^data.fanox3{:,numericCols3};
% 
% data.dU_age3 = [[data.age3,data.dU3.d238U_0_16];[flip(data.age3),flip(data.dU3.d238U_0_84)]];
% data.fanox_age3 = [[data.age3,data.fanox3.fanox_0_16];[flip(data.age3),flip(data.fanox3.fanox_0_84)]];
% data.num_idx3 = 1:size(data.dU_age3,1);

for i = 1 : length(color.rgb)
    data.clr(i,:) = sscanf(string(color.rgb(i)),"%f,%f,%f")';
end

label = ["A","B","C","D","E","F","G","H","I","J","K","L"];
row = 1:4;
col = 1:3;
label = reshape(label,length(row),length(col));
%% 
% Data ordered accordingly:
% sLL(tracer,level,time) contains low latitude ocean data
% sHL(tracer,level,time) contains high latitude ocean data
% sAT(tracer,zone,time) contains atmospheric data
% st(time) contains model time
% Zones are 
% 1:low latitude (LL), 2:high latitude (HL)
% Ocean tracers are 
% 1:T, 2:S, 3:PO4, 4:DIC, 5:DI13C, 6:DI14C, 7:Alk, 8:O, 9:18O, 
%10:CH4, 11:13CH4, 12: NO3, 13:NH_4, 14:H2S
%Atmosphere tracers are
%1:T, 2:pCH4, 3:pN20, 4:pCO2, 5:pC(13)O2, 6: pC(14)O2, 7:pC(13)H4, 8:pO2
%Land biomasses are
%1:Leaves, 2: Wood, 3:Litter, 4. Soil


for ii=1:length(sAT(4,1,:))
    st(ii)=100*(ii);
           
    [pCO2,K0,CO2,CO3,HCO3,GAM]=CarSys_M(sLL(1,1,ii),sLL(2,1,ii),sLL(4,1,ii),sLL(7,1,ii));
    GAMLL(1,1,ii)=GAM;
    [pCO2,K0,CO2,CO3,HCO3,GAM]=CarSys_M(sHL(1,1,ii),sHL(2,1,ii),sHL(4,1,ii),sHL(7,1,ii));
    GAMHL(1,1,ii)=GAM;
     
       for i=2:n
         if dwcLL(i,ii)<0.1
             CCDLL10(ii)=(i-1)*100+(dwcLL(i-1,ii)-0.1)./(dwcLL(i-1,ii)-dwcLL(i,ii))*100; break
         end
        end
       for i=2:n
         if dwcHL(i,ii)<0.1
             CCDHL10(ii)=(i-1)*100+(dwcHL(i-1,ii)-0.1)./(dwcHL(i-1,ii)-dwcHL(i,ii))*100; break
         else
             CCDHL10(ii)=5500;
         end
        end
end     

% Calculate delta13 and Delta 14 values in permil for ocean and atmosphere
%-----
t_kyr = st/1e3;
 
[~,idx] = min(abs(t_kyr-19.5));

d13a   = squeeze( (sAT(5,1,:)./sAT(4,1,:)/R13pdb-1)*1e3 );
d13ao   = squeeze( (sAT(5,1,idx)./sAT(4,1,idx)/R13pdb-1)*1e3 );      
d13ad = d13a - d13ao;                                            %atmosphere

d13LL100mo = squeeze( (sLL(5,1,idx)./sLL(4,1,idx)/R13pdb-1)*1e3 );
d13LL1000mo = squeeze( (sLL(5,10,idx)./sLL(4,10,idx)/R13pdb-1)*1e3 );
d13LL3000mo = squeeze( (sLL(5,30,idx)./sLL(4,30,idx)/R13pdb-1)*1e3 );
d13LL5000mo = squeeze( (sLL(5,50,idx)./sLL(4,50,idx)/R13pdb-1)*1e3 );

d13LL100m = squeeze( (sLL(5,1,:)./sLL(4,1,:)/R13pdb-1)*1e3 );
d13LL1000m = squeeze( (sLL(5,10,:)./sLL(4,10,:)/R13pdb-1)*1e3 );
d13LL3000m = squeeze( (sLL(5,30,:)./sLL(4,30,:)/R13pdb-1)*1e3 );
d13LL5000m = squeeze( (sLL(5,50,:)./sLL(4,50,:)/R13pdb-1)*1e3 );

d13LL100md = d13LL100m - d13LL100mo;                       %ocean, 100m
d13LL1000md = d13LL1000m - d13LL1000mo;                    %ocean, 1000m
d13LL3000md = d13LL3000m - d13LL3000mo;                    %ocean, 3000m 
d13LL5000md = d13LL5000m - d13LL5000mo;

d18cLL100mo = squeeze( sLL(9,1,1).*1000+(16.5-sLL(1,1,1))./4.8 );
d18cLL1000mo = squeeze( sLL(9,10,1).*1000+(16.5-sLL(1,10,1))./4.8 );
d18cLL3000mo = squeeze( sLL(9,30,1).*1000+(16.5-sLL(1,30,1))./4.8 );

d18cLL100m = squeeze( sLL(9,1,:).*1000+(16.5-sLL(1,1,:))./4.8 );
d18cLL1000m = squeeze( sLL(9,10,:).*1000+(16.5-sLL(1,10,:))./4.8 );
d18cLL3000m = squeeze( sLL(9,30,:).*1000+(16.5-sLL(1,30,:))./4.8 );

d18cLL100md = d18cLL100m-d18cLL100mo;                      %ocean, 100m
d18cLL1000md = d18cLL1000m-d18cLL1000mo;                   %ocean, 1000m
d18cLL3000md = d18cLL3000m-d18cLL3000mo;                   %ocean, 3000m
%-----
%%
% Establish the atmospheric temperature profile, 2.order legendre pol. in 
% sine of lat, Ta(lat) = PTa(1) + .5*PTa(2) * ( 3*sin(lat)^2 -1 )
% Coefficients are calculated so that the area weighted mean of the profile
% matches the surface mean temperatures in each sector.
for ii=1:length(sAT(1,1,:))
  CTa(1,:) = [ 1 .5*(sin(fdiv)^2-1)                   ];
  CTa(2,:) = [ 1 .5*(sin(fdiv)-sin(fdiv)^3)/(1-sin(fdiv)) ];
  RTa      = [ sAT(1,1,ii) sAT(1,2,ii)]';
  PTa(1:2,ii)= CTa\RTa;
end
%% Calculate time-dependent changes in the anoxic seafloor area based on dissolved oxygen levels
%**********************************************************************************************%
% Setting initial parameters 
OLL = (LL(8,:)*1000)'; % 18O in LL mmol/m3
OHL = (HL(8,:)*1000)'; % 18O in HL mmol/m3
fanox_init = 0.0007; % The proportion of anoxic seafloor area in the total seafloor area 
crit = 3; % the threshold of O2min for anoxic condition, is equal to that of denitrification (mmol/m3)  
fHL = aHL/(aHL+aLL); % the area proportion of high latitude in ocean area
fLL = aLL/(aHL+aLL); % % the area proportion of low latitude in ocean area

% Calculate proportion of each layer in total seafloor area
RS_LL = nan(55,1);
RS_HL =  nan(55,1);
RS_LL(1) = 0;
RS_HL(1) = 0;
for i = 1:54
   RS_LL(i+1) = (GAw(i)-GAw(i+1))/GAw(1); 
   RS_HL(i+1) = (GAc(i)-GAc(i+1))/GAc(1);
end

% Calculate k1 values
syms k1 % variability rate（k < 0）
eqn = fLL*sum(exp(k1.*(OLL/crit-1)).*RS_LL) + fHL*sum(exp(k1.*(OHL/crit-1)).*RS_HL) == fanox_init; % a negative exponential relationship between dissolved oxygen concentration and seafloor anoxia 
k1 = double(vpasolve(eqn, k1));  

% Extract the dissolved oxygen concentration 
ZLL = permute(sLL(8,:,:),[2,3,1])*1000; % dissolved oxygen concentration in LL mmol/m3
ZHL = permute(sHL(8,:,:),[2,3,1])*1000; % dissolved oxygen concentration in HL mmol/m3

%
RS_LL_matr = repmat(RS_LL,1,length(t_kyr));
RS_HL_matr = repmat(RS_HL,1,length(t_kyr));

% Transfer dissolved oxygen concentration into anoxic seafloor area 
ZLL_1 = ZLL;
ZLL_1(ZLL_1<=crit) = crit; % The O2min are utilized to repalce the dissolved oxygenation concentration if the dissolved oxygenation concentration are less than O2min     
deg_LL = exp(k1.*(ZLL_1/crit-1)); % Transfer dissolved oxygen concentration into seaflooor anoxic degree for each layer  
ZHL_1 = ZHL;
ZHL(ZHL<=crit) = crit;
deg_HL = exp(k1.*(ZHL/crit-1));
fanox_LL = sum(deg_LL.*RS_LL_matr)*100; % Calculate seafloor aoxic area in LL （unit: %）
fanox_HL = sum(deg_HL.*RS_HL_matr)*100; % Calculate seafloor aoxic area in HL unit: %）
fanox_global = fLL*fanox_LL + fHL*fanox_HL; % % Calculate global seafloor aoxic area in HL unit: %）
%**************************************************************************************************%
%% 
xmin = 0;
xmax = 100;
xtick_inter = 10*xmax/100; 
onset = 22.3; % The time are correspond to the onset of PETM
bef = 18.5; % The time are correspond to the first peak of anoxic seafloor area 

figure("Position",[0,0,1100,650])
ha = tight_subplot(4,3,[0.04,0.06],[0.08,0.04],[0.06,0.02]);

axes(ha(3*0+1))
plt1 = plot(t_kyr,2*scar(1,:)/mgt*sy,"LineWidth",1.4); hold on
plt2 = plot(t_kyr,2*scar(2,:)/mgt*sy,"LineWidth",1.4);
xlim([xmin,xmax]);
xticks(xmin:xtick_inter:xmax);
ylim([0,1.1]);
limY = get(gca, 'Ylim');
ylabel('Carbon input (GtCyr-1)')
line([onset,onset,NaN,bef,bef],[limY(1),limY(2), NaN, limY(1),limY(2)],"LineStyle",'--',"Color","#B4B5B5","LineWidth",1.2);
text_norm(0.03,0.92,label(1,1),10);
box off
BorderLine(gca,1);
lgd = legend([plt1,plt2],["Hydrate methane ("+ input.meth.total*2 + "Gt )","Thermogenic methane ("+ input.mett.total*2 + "Gt )"],"FontSize",8,"EdgeColor","none");
lgd.ItemTokenSize = [18,1];
set(gca,"xticklabel",[],"linewidth", 1,"FontSize",10,"FontName", "Times","TickLength",[0.02,0.03],"Layer","top");

axes(ha(3*1+1))
lin1 = line([xmin,xmax],[methc13,methc13],"LineWidth",1.4); hold on
lin2 = line([xmin,xmax],[mettc13,mettc13],"LineWidth",1.4);
xlim([xmin,xmax]);
xticks(xmin:xtick_inter:xmax);
ylim([-70,0]);
yticks(-70:10:0);
ticklabels("y",yticks,2,"F");
ylabel("\delta^{13}C");
limY = get(gca, 'Ylim');
line([onset,onset,NaN,bef,bef],[limY(1),limY(2), NaN, limY(1),limY(2)],"LineStyle",'--',"Color","#B4B5B5","LineWidth",1.2);
text_norm(0.03,0.92,label(2,1),10);
box off
BorderLine(gca,1);
lgd = legend([lin1,lin2],["Hydrate methane","Thermogenic methane"],"FontSize",8,"EdgeColor","none");
lgd.ItemTokenSize = [18,1];
set(gca,"xticklabel",[],"linewidth", 1,"FontSize",10,"FontName", "Times","TickLength",[0.02,0.03],"Layer","top");

axes(ha(3*2+1))
plt1 = plot(t_kyr,squeeze(sAT(2,1,:))*1e6,"LineWidth",1.4); hold on
plt2 = plot(t_kyr,squeeze(sAT(4,1,:))*1e6/100,"LineWidth",1.4);
xlim([xmin,xmax]);
xticks(xmin:xtick_inter:xmax);
ylim([0,37])
ylabel({'pCO_2 (× 10^{-2} ppm)','pCH4 (ppm)'})
limY = get(gca, 'Ylim');
line([onset,onset,NaN,bef,bef],[limY(1),limY(2), NaN, limY(1),limY(2)],"LineStyle",'--',"Color","#B4B5B5","LineWidth",1.2);
text_norm(0.03,0.92,label(3,1),10);
box off
BorderLine(gca,1);
lgd = legend([plt1,plt2],["CO2","CH4"],"FontSize",8,"EdgeColor","none","Location","northeast");
lgd.ItemTokenSize = [18,1];
set(gca,"xticklabel",[],"linewidth", 1,"FontSize",10,"FontName", "Times","TickLength",[0.02,0.03],"Layer","top");

axes(ha(3*3+1))
plt1 = plot(t_kyr,squeeze(sLL(1,1,:)),"LineWidth",1.4); hold on
plt2 = plot(t_kyr,PTa(1,:),"LineWidth",1.4); 
xlim([xmin,xmax]);
xticks(xmin:xtick_inter:xmax);
xlabel("t (kyr)");
ylabel('T (^oC)')
limY = get(gca, 'Ylim');
line([onset,onset,NaN,bef,bef],[limY(1),limY(2), NaN, limY(1),limY(2)],"LineStyle",'--',"Color","#B4B5B5","LineWidth",1.2);
text_norm(0.03,0.92,label(4,1),10);
box off
BorderLine(gca,1);
lgd = legend([plt1,plt2],["Ocean (100 m)","Average atmosphere"],"EdgeColor","none","FontSize",8,"Location","southeast");
lgd.ItemTokenSize = [18,1]; 
set(gca,"linewidth", 1,"FontSize",10,"FontName", "Times","TickLength",[0.02,0.03],"Layer","top");

axes(ha(3*0+2))
plt1 = plot(t_kyr,d13LL100md,"LineWidth",1.4); hold on
plt2 = plot(t_kyr,d13LL5000md,"LineWidth",1.4);
scatter(data.age1,data.dC1-mean(data.dC1(1:60)),7,"Marker","o","MarkerEdgeColor","#404040","MarkerFaceColor","#404040","LineWidth",0.01);
xlim([xmin,xmax]);
xticks(xmin:xtick_inter:xmax);
ylim([-7.5,3.5]);
yticks(-7:1:3);
ticklabels("y",yticks,2,"F");
ylabel('\Delta\delta^{13}C (‰)');
limY = get(gca, 'Ylim');
line([onset,onset,NaN,bef,bef],[limY(1),limY(2), NaN, limY(1),limY(2)],"LineStyle",'--',"Color","#B4B5B5","LineWidth",1.2);
text_norm(0.03,0.92,label(1,2),10);
box off
BorderLine(gca,1);
lgd = legend([plt1,plt2],["100 m","5000 m"],"EdgeColor","none","FontSize",8);
lgd.ItemTokenSize = [18,1];
set(gca,"xticklabel",[],"linewidth", 1,"FontSize",10,"FontName", "Times","TickLength",[0.02,0.03],"Layer","top");

axes(ha(3*1+2))
plt1 = plot(t_kyr,squeeze(sLL(8,1,:))*1000,"LineWidth",1.4);  hold on
plt2 = plot(t_kyr,squeeze(sLL(8,10,:))*1000,"LineWidth",1.4);
plt3 = plot(t_kyr,squeeze(sLL(8,30,:))*1000,"LineWidth",1.4);
xlim([xmin,xmax]);
xticks(xmin:xtick_inter:xmax);
ylim([0,300]);
ylabel('O_2 (mmol m^{-3})');
limY = get(gca, 'Ylim');
lin1 = line([xmin,xmax],[3,3],"LineStyle","--","Color","#000000","LineWidth",1.4);
line([onset,onset,NaN,bef,bef],[limY(1),limY(2), NaN, limY(1),limY(2)],"LineStyle",'--',"Color","#B4B5B5","LineWidth",1.2);
text_norm(0.03,0.92,label(2,2),10);
box off
BorderLine(gca,1);
lgd = legend([plt1,plt2,plt3,lin1],["100 m","1000 m","3000 m","O_{min}"],"EdgeColor","none","FontSize",8);
lgd.ItemTokenSize = [18,1];
set(gca,"xticklabel",[],"linewidth", 1,"FontSize",10,"FontName", "Times","TickLength",[0.02,0.03],"Layer","top");

axes(ha(3*2+2))
plt1 = plot(t_kyr,squeeze(sLL(12,10,:))*1000,"LineWidth",1.4); hold on
plt2 = plot(t_kyr,squeeze(sLL(12,30,:))*1000,"LineWidth",1.4);
xlim([xmin,xmax]);
xticks(xmin:xtick_inter:xmax);
ylim([0,50]);
limY = get(gca, 'Ylim');
ylabel('NO_{3}^{-2} (mmol/m3)');
lin1 = line([xmin,xmax],[0.03,0.03],"LineStyle","--","Color","#000000","LineWidth",1.4);
line([onset,onset,NaN,bef,bef],[limY(1),limY(2), NaN, limY(1),limY(2)],"LineStyle",'--',"Color","#B4B5B5","LineWidth",1.2);
text_norm(0.03,0.92,label(3,2),10);
box off
BorderLine(gca,1);
lgd = legend([plt1,plt2,lin1],["1000 m","3000 m","N_{min}"],"EdgeColor","none","FontSize",8);
lgd.ItemTokenSize = [18,1];
set(gca,"xticklabel",[],"linewidth", 1,"FontSize",10,"FontName", "Times","TickLength",[0.02,0.03],"Layer","top");

axes(ha(3*3+2))
plt3 = plot(t_kyr,fanox_global, "LineWidth",1.4);
patch("Faces",data.num_idx3,"Vertices",data.fanox_age3,"FaceColor",data.clr(2, :),"EdgeColor","None","FaceAlpha",0.8,"LineWidth",1.4); hold on
plt4 = plot(data.age3, data.fanox3.fanox_0_5,"-","Color",data.clr(1, :),"LineWidth",1.4); 
xlim([xmin,xmax]);
xticks(xmin:xtick_inter:xmax);
xlabel("t (kyr)");
limY = get(gca, 'Ylim');
ylabel("Seafloor anoxic (%)");
line([onset,onset,NaN,bef,bef],[limY(1),limY(2), NaN, limY(1),limY(2)],"LineStyle",'--',"Color","#B4B5B5","LineWidth",1.2);
text_norm(0.03,0.92,label(4,2),10);
box off
BorderLine(gca,1);
lgd = legend([plt3,plt4],["Model", "Observed"],"Location","southeast","FontSize",8,"EdgeColor","none");
lgd.ItemTokenSize = [18,1];
set(gca,"linewidth", 1,"FontSize",10,"FontName", "Times","TickLength",[0.01,0.015],"Layer","top");
axes(ha(3*0+3))
plot(data.Ba.Age_westerhold,data.Ba.x_138Ba,"LineWidth",1,"Color","#000000"); hold on
errorbar(data.Ba.Age_westerhold,data.Ba.x_138Ba,data.Ba.x2SD_1, "o", "vertical","MarkerSize",4.7,"MarkerEdgeColor","#000000","MarkerFaceColor","#93BBDB","LineWidth",0.3,"CapSize",0,"Color","#757475"); 
xlim([xmin,xmax]);
xticks(xmin:xtick_inter:xmax);
ylim([0.3,0.9]);
limY = get(gca, 'Ylim');
line([onset,onset,NaN,bef,bef],[limY(1),limY(2), NaN, limY(1),limY(2)],"LineStyle",'--',"Color","#B4B5B5","LineWidth",1.2);
yticks(0.3:0.15:0.9);
ticklabels("y",0.3:0.15:0.9,2,"F");
ylabel("\delta^{138}Ba (‰)");
text_norm(0.03,0.92,label(1,3),10);
box off
BorderLine(gca);
set(gca,"xticklabel",[],"linewidth", 1,"FontSize",10,"FontName", "Times","TickLength",[0.02,0.03],"Layer","top");

axes(ha(3*1+3))
Tr  = 10;         % [oC] Reference temperatur, st. val 10
p1  = 1;          % [-] O.Marchal 98,Maier-Reimer 93
p2  = 0.1565;     % [oC^-1] Maier-Reimer 93; revised 22/08/2014
rpLL(1,1,:)  = rpm* p1*exp(p2*(sLL(1,1,:)-Tr)) ...
    ./(1+p1*exp(p2*(sLL(1,1,:)-Tr))).*(GAMLL(1,1,:)-1)./(1+(GAMLL(1,1,:)-1));  
rpHL(1,1,:)  = rpm* p1*exp(p2*(sHL(1,1,:)-Tr)) ...
    ./(1+p1*exp(p2*(sHL(1,1,:)-Tr))).*(GAMHL(1,1,:)-1)./(1+(GAMHL(1,1,:)-1));
NPLL(1,1,:)=-2*(srLL(3,1,:)-RoLL(1,1,:))*rcp*12.011*sy/1e15;
NPHL(1,1,:)=-2*(srHL(3,1,:)-RoHL(1,1,:))*rcp*12.011*sy/1e15;
plt1 = plot(t_kyr,squeeze(NPLL(1,1,:)),"LineWidth",1.4); hold on   %LL OrgC 
plt2 = plot(t_kyr,squeeze(NPHL(1,1,:)),"LineWidth",1.4);           %HL OrgC 

xlim([xmin,xmax]);
xticks(xmin:xtick_inter:xmax);
ylim([2.5,4.5]);
ylabel(['New production',"(Gton C/yr)"]);
limY = get(gca, 'Ylim');
line([onset,onset,NaN,bef,bef],[limY(1),limY(2), NaN, limY(1),limY(2)],"LineStyle",'--',"Color","#B4B5B5","LineWidth",1.2);
text_norm(0.03,0.92,label(2,3),10);
box off
BorderLine(gca,1);
lgd = legend([plt1,plt2],["Low latitude","High latitude"],"EdgeColor","none","FontSize",8,"Location","southeast");
lgd.ItemTokenSize = [18,1];
set(gca,"xticklabel",[],"linewidth", 1,"FontSize",10,"FontName", "Times","TickLength",[0.02,0.03],"Layer","top");

axes(ha(3*2+3))
plt1 = plot(t_kyr,2*(BurcLL+BurcHL)/mgt*sy/rcp,"LineWidth",1.4); hold on
plt2 = plot(t_kyr,2*(squeeze(RoLL(1,1,:))+squeeze(RoHL(1,1,:)))/1e15*30.974*sy,"LineWidth",1.4); 
xlim([xmin,xmax]);
xticks(xmin:xtick_inter:xmax);
ylim([0,6e-3])
limY = get(gca, 'Ylim');
ylabel({'P Flux,(GtP yr^{-1})'})
line([onset,onset,NaN,bef,bef],[limY(1),limY(2), NaN, limY(1),limY(2)],"LineStyle",'--',"Color","#B4B5B5","LineWidth",1.2);
text_norm(0.03,0.92,label(3,3),10);
box off
BorderLine(gca,1);
lgd = legend([plt1,plt2],["burial P","Input P"],"EdgeColor","none","FontSize",8,"Location","southeast");
lgd.ItemTokenSize = [18,1];
set(gca,"XTickLabel",[],"LineWidth", 1,"FontSize",10,"FontName", "Times","TickLength",[0.02,0.03],"Layer","top");

axes(ha(3*3+3))
plt2 = plot(t_kyr,squeeze(sLL(3,10,:))*1000,"LineWidth",1.4); hold on
plt3 = plot(t_kyr,squeeze(sLL(3,30,:))*1000,"LineWidth",1.4);
xlim([xmin,xmax]);
xticks(xmin:xtick_inter:xmax);
xlabel('t (kyr)');
ylabel('P mmol/m^3');
limY = get(gca, 'Ylim');
line([onset,onset,NaN,bef,bef],[limY(1),limY(2), NaN, limY(1),limY(2)],"LineStyle",'--',"Color","#B4B5B5","LineWidth",1.2);
text_norm(0.03,0.92,label(4,3),10);
box off
BorderLine(gca,1);
lgd = legend([plt2, plt3],["1000 m", "3000 m"],"EdgeColor","none","FontSize",8, "Location","southeast");
lgd.ItemTokenSize = [18,1];
set(gca,"linewidth", 1,"FontSize",10,"FontName", "Times","TickLength",[0.02,0.03],"Layer","top");

print(gcf,"Figure\DCESS_Evol","-dpng","-r600")
%% Carlibate organic carbon dissolution based on dissolved oxygenation concentration 
figure("Position",[0,0,600,400])
[ha, ~] = tight_subplot(1,3,[0,0.07],[0.11,0.02],[0.08,0.02]);

axes(ha(1))
k2 = 0; % varibility rate
for i = 1:50:1000
oxy_norm = ZLL(:,i)./(LL(8,:)'*1000); % Normalized dissolved oxygenation concentration  
pOrgCLL_alt = exp(k2*(oxy_norm-1)).*pOrgCLL(:,i); % positive exponential function between organic carbon dissolution and dissolved oxygenation concentration 
plot(pOrgCLL_alt,1:55); hold on
end
xlim([0,1]);
xticks(0:0.2:1)
ticklabels("x",xticks,2,"T");
ylim([0,55]);
ylabel("Depth (km)");
text_norm(0.05,0.05,"k = "+k2,12);
set(gca,"YDir","reverse","LineWidth", 1,"FontSize",10,"FontName", "Times","TickLength",[0.02,0.03],"Layer","top");

axes(ha(2))
k2 = 0.5;
for i = 1:50:1000
oxy_norm = ZLL(:,i)./(LL(8,:)'*1000);
pOrgCLL_alt = exp(k2*(oxy_norm-1)).*pOrgCLL(:,i);
plot(pOrgCLL_alt,1:55); hold on
end
xlim([0,1])
xticks(0:0.2:1)
ticklabels("x",xticks,2,"T");
xlabel("Orgnic carbon dissolution (%)")
ylim([0,55]);
text_norm(0.05,0.05,"k = "+k2,12);
set(gca,"YDir","reverse","LineWidth", 1,"FontSize",10,"FontName", "Times","TickLength",[0.02,0.03],"Layer","top");

axes(ha(3))
k2 = 1;
for i = 1:50:1000
oxy_norm = ZLL(:,i)./(LL(8,:)'*1000);
pOrgCLL_alt = exp(k2*(oxy_norm-1)).*pOrgCLL(:,i);
plot(pOrgCLL_alt,1:55); hold on
end
xlim([0,1])
xticks(0:0.2:1)
ticklabels("x",xticks,2,"T");
ylim([0,55]);
text_norm(0.05,0.05,"k = "+k2,12);
set(gca,"YDir","reverse","LineWidth", 1,"FontSize",10,"FontName", "Times","TickLength",[0.02,0.03],"Layer","top");

print(gcf,"Figure/Organic carbon dissolution","-dpng","-r600");
%%
% Define the range for x and y
zcent= [dm/2 dm+(d:d:(n-1)*d)-d/2]/1e3;   % Vertical center of boxes

% Create the grid
[X, Y] = meshgrid(t_kyr, zcent);

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

print(gcf,"Figure\Oxygen concentration_Contour","-dpng","-r600");