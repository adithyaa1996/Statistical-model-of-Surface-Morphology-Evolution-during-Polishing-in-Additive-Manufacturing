%% Model Parameters Fitting with Shilan Jin's data - Nov 11th 2021 - Mar 21st 2022

% Consider a 0.25mm x 0.25mm surface for the purpose of modelling.
% Surface packing density is assumed to be in the range of 0.7 for Ti-6Al-4V
% The minimun and maximum radii of asperities for the circle packing
% algorithm are specified accordingly.

% Stage 0 has already been realized.
% Load Stage0_Aug13.mat from the folder Aug_13_Stage0
% Then run the following blocks.

%% Initial - Network Formation Plot - DeLaunay Triangulation

figure;
DT = delaunayTriangulation(x_pos, y_pos);
triplot(DT, '--');
title('Surface Asperities - Network Formation');
xlabel('X (10^{-6})m')
ylabel('Y (10^{-6})m')
pub_fig;

%% Stored Variables from Stage0_Aug13.mat file

% grid_length = 250;  % Units in micrometers
% [Locations,radii_m0] = random_circle_packing_rectangle([grid_length grid_length],3.14,25,true);
% x_pos = Locations(:,1);
% y_pos = Locations(:,2);
% 
% N = size(Locations,1);  % number of asperities
% pd = makedist('Weibull',39.3, 4.27);
% heights_m0 = random(pd, N, 1);       % initial height distribution - weibull 
% 
% asperity_area = pi * sum(radii_m0.^2);  % units um^2
% packing_density = asperity_area / (grid_length^2);
% Temp = 20*ones(N,1);  % baseline temperature
% params = cell(1,2); 
% E = 0.2;
% dd = 0.01;
% params{1} = E;
% params{2} = dd;

%% Other variables to be Initialized

Force = 53.38*ones(1,10);
kappa_1 = 15*ones(1,10);
kappa_2 = 0.674*ones(1,10);
speed = 2.36*ones(1,10);

alpha = [0.30456 0.108642 1.372173 2.013025 1.343399 1.650823 0.413329 2.34497973 0.04668737 1.420183099];
sigma = [0.001827 0.000103 0.003548 0.002703 0.004522 0.007872 0.000083 0.00000016 0.00000042 0.00000000025];
beta = [0.000242 0.000029 0.000233 0.000038 0.001900 0.001363 0.001472 0.00008597 0.00000021 0.0000631721];
Fiedler = zeros(1,10);

%% Stage 0 - Data and Model Comparison - Plotting

% Table0 = readtable('matrix_stage0.csv');
% data_stage0 = table2array(Table0);
% data_stage0 = reshape(data_stage0, [3264 1]);
% data0 = mean(data_stage0) - data_stage0;
% count = 0;
% 
% % Stage 0 Data Statistics
% max_stage0 = max(data0);
% mean_stage0 = mean(data0);
% stdev_stage0 = std(data0);
% min_stage0 = min(data0);
% Hurst_stage0 = Gen_hurst(data0);
% Sa_stage0 = 0;
% for i = 1: 3264
%     Sa_stage0 = Sa_stage0 + ((1/3264)*(abs(data0(i) - mean_stage0)));
% end
% P2P_stage0 = (sum(maxk(data_stage0,320,1),'all') - sum(mink(data_stage0,320,1), 'all'))/320;
% 
% % Model 0 Statistics
% max_model_0 = max(heights_m0);
% min_model_0 = min(heights_m0);
% active_nodes = 0;
% [Model_H0, mean_model_0, stdev_model_0, Sa_model_0, Hurst_model_0, KS_test_0, p_value_0] = surface_roughness(radii_m0, initial_asperity_radii, heights_m0, N, data_stage0, packing_density);
% P2P_model_0 = (sum(maxk(heights_m0,20,1),'all') - sum(mink(heights_m0,20,1), 'all'))/20;

figure
subplot(2,3,1)
histogram(data_stage0, 'Normalization', 'probability', 'BinWidth', 1);
pub_fig;
ylim([0 0.06]);
xlim([0 50]);
title('Histogram of Stage 0 Data');
hold on
plot([0:0.5:50], wblpdf([0:0.5:50], 27.4283, 2.6496), 'r');
xlabel('Heights (10^{-6} m)')
ylabel('Frequency')

subplot(2,3,2)
H = histogram(max_model_0 - Model_H0, 'Normalization', 'probability', 'BinWidth', 1);
title('Histogram of generated Initial surface');
hold on
plot([0:0.5:50], wblpdf([0:0.5:50], 27.4283, 2.6496), 'r');
pub_fig;
ylim([0 0.06]);
xlim([0 50]);
xlabel('Heights (10^{-6} m)')
ylabel('Frequency')

subplot(2,3,4)
[Values_s, edges] = histcounts(data_stage0, linspace(0, 50, 51), 'Normalization', 'cdf');
centers = (edges(1:end-1)+edges(2:end))/2;
plot(Values_s, centers, 'r-', 'LineWidth', 2.5)
pub_fig;
xticks([0 0.2 0.4 0.6 0.8 1]);
set(gca,'ydir','reverse');
hold on
%subplot(2,3,5)
[Values_m, edges] = histcounts(max_model_0 - Model_H0, linspace(0, 50, 51), 'Normalization', 'cdf');
centers = (edges(1:end-1)+ edges(2:end))/2;
plot(Values_m, centers, 'b-', 'LineWidth', 2.5);
title('Bearing Area Curves');
pub_fig;
xticks([0 0.2 0.4 0.6 0.8 1]);
ylabel('Heights (10^{-6} m)')
xlabel('Quantile')
set(gca,'ydir','reverse');
legend('Stage 0 data','Model surface')

for i = 2:1:length(Values_s)
    Values_s(i) = Values_s(i) - Values_s(i-1);
    Values_m(i) = Values_m(i) - Values_m(i-1);
end
KL_div_test_stage0 = KLDiv(Values_s, Values_m);


subplot(2,3,3)
pos = get(subplot(2,3,3), 'Position');
Statistics = {'Mean';'Std dev'; 'Sa value'; 'Max'; 'Min'; 'Mean P2P'; '(KS & KL)'};
Stage_0 = [mean_stage0; stdev_stage0; Sa_stage0; max_stage0; min_stage0; P2P_stage0; p_value_0];
Model_0 = [mean_model_0; stdev_model_0; Sa_model_0; max_model_0; min_model_0; P2P_model_0; KL_div_test_stage0];
T = table(Stage_0, Model_0, 'RowNames', Statistics);
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,...
    'RowName',T.Properties.RowNames,'Units', 'Normalized','ColumnWidth', {70}, 'Position', pos);
pub_fig;
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('Surface Characteristics');
xlabel('All units in µm');

subplot(2,3,5)
qqplot(max_model_0 - Model_H0, data_stage0)
pub_fig;
title("Q-Q plot")
ylim([0 50]);
xlim([0 50]);
xlabel("Model Quantiles")
ylabel("Stage 0 - Data Quantiles")

%% Stage 0 to Stage 1 Parameter Fitting 
%%%%%%% Asperity network %%%%%%%%%

MRR(1) = 0;

Time = 300; 
k = 2; dt = 0.2;
a_0 = 2.5;
Horizon_01 = 0:dt:Time;
Horizon = Horizon_01;
T = zeros(N,length(Horizon));  T(:,1) = Temp;
s = zeros(N,length(Horizon));  s(:,1) = sigma(1);
h = zeros(N,length(Horizon));  h(:,1) = heights_m0;
r = zeros(N,length(Horizon));  r(:,1) = radii_m0;
F = zeros(N,length(Horizon));                   % Force carried by asperities
Rho = zeros(N,1); 
diff_heights = zeros(N,1);
count_active_nodes = zeros(length(Horizon),1);

st_d_crit = zeros(length(Horizon),1);
st_active_nodes = zeros(length(Horizon),1);
vector = 0: params{2}: max(heights_m0);
vectorlength = length(vector);

for t = dt:dt:Time
    
    d_crit = solve_for_d(Force(1),r(:,k-1),h(:,k-1),params,vectorlength);
    active_nodes = find(h(:,k-1) > d_crit); % load bearing asperities
    F(active_nodes, k) = (2/3)*E* sqrt(abs(r(active_nodes,k-1))).*((abs(h(active_nodes,k-1) - d_crit.*ones(length(active_nodes),1))).^(3/2));   % force carried by each node. 
    count_active_nodes(k-1,1) = length(active_nodes);
    % Temp dynamics
    T(:,k) = T(:,k-1) + dt*(-alpha(1)*(T(:,k-1)- Temp) + kappa_1(1)*F(:, k-1) + kappa_2(1)*speed(1)); 
    
    % Sigma dynamics
    s(:,k) = s(:,k-1) + (beta(1).*(T(:,k) - T(:,k-1))./Temp); 

%     Indicator = (T(:,k-1) > (Temp));
%     if (k==2)
%        s(:,k) = s(:,k-1) + (dt.*s(:,k-1).*Indicator.*T(:,k-1)*(0.1)./Temp);
%     else
%        s(:,k) = s(:,k-1) + (dt.*s(:,k-1).*Indicator.*(T(:,k-1) - T(:,k-2))./Temp); 
%     end
    
    % Differential height dynamics
    diff_heights(:,1) = (dt.*s(:,k-1).*T(:,k-1).*F(:,k-1));
    h(:,k) = h(:,k-1) - diff_heights;    
    
    % Radii update and MRR computation
    Rho(:,1) = h(:,k-1)./(4*r(:,k-1));
    MR_coeff_eta = exp(-(h(:,k-1)/max(heights_m0).^2)).*exp(-(a_0./r(:,k-1)).^0.0833);
    r(:,k) = r(:,k-1) + (MR_coeff_eta.*Rho(:,1).*(diff_heights));
    MRR(1) = MRR(1) + sum(pi.*(diff_heights.^2).*h(:,k-1).*MR_coeff_eta.*(1 - MR_coeff_eta)/4);
    
    % store
    st_d_crit(k-1) = d_crit;
    st_active_nodes(k-1) = length(active_nodes);
    k = k+1;
end

heights_m1 = h(:,length(Horizon));
radii_m1 = r(:,length(Horizon));
MRR(1) = MRR(1) / (250*250*Time/60);

% height, sigma and Temp change plots
figure
subplot(1,3,1)
h_01 = h(active_nodes,:);
plot(Horizon,h(active_nodes,:),'LineWidth',2)
pub_fig;
hold on
% title('Height change dynamics of active asperities')
xlabel('Time (sec)')
ylabel('Heights of active asperities (μm)')

subplot(1,3,2)
plot(Horizon,r(active_nodes,:),'LineWidth',2)
pub_fig;
hold on
% title('Radii change of asperities')
xlabel('Time (sec)')
ylabel('Radii of active asperities (μm)')

subplot(1,3,3)
Temp_01 = T(active_nodes,:);
plot(Horizon,T(active_nodes,:),'LineWidth',2)
pub_fig;
hold on
% title('Temp variation: Stage 0 - 1')
xlabel('Time (sec)')
ylabel('Temperature of active asperities (degree C)')

%%

% figure
% plot(Horizon,count_active_nodes,'LineWidth',1.6)
% pub_fig;
% hold on
% xlabel('Time (sec)')
% ylabel('Count of active asperities')


%% Stage 1 - Data and Model Comparison - Plotting

Table1 = readtable('matrix_stage1.csv');
data_stage1 = table2array(Table1);
data_p2p_stage1 = reshape(data_stage1, [102 32]);
data_stage1 = reshape(data_stage1, [3264 1]);
data1 = mean(data_stage1) - data_stage1;

% Stage 1 Data Statistics
max_stage1 = max(data1);
mean_stage1 = mean(data1);
stdev_stage1 = std(data1);
min_stage1 = min(data1);
Hurst_stage1 = Gen_hurst(data1);
Sa_stage1 = 0;
for i = 1: 3264
    Sa_stage1 = Sa_stage1 + ((1/3264)*(abs(data_stage1(i) - mean_stage1)));
end
P2P_stage1 = (sum(maxk(data_stage1,320,1),'all') - sum(mink(data_stage1,320,1), 'all'))/320;

% Model 1 Statistics
max_model_1 = max(heights_m1);
min_model_1 = min(heights_m1);
[Model_H1, mean_model_1, stdev_model_1, Sa_model_1, Hurst_model_1, KS_test_1, p_value_1] = surface_roughness(radii_m1, initial_asperity_radii, heights_m1, N, data_stage1, packing_density);
P2P_model_1 = (sum(maxk(heights_m1,20,1),'all') - sum(mink(heights_m1,20,1), 'all'))/20;

vector = 0: params{2}: max(heights_m0);
vectorlength = length(vector);
[AdjMatrix_1, Fiedler(1), E1, bins_CC1] = graph_evolution(x_pos, y_pos, packing_density, Force(1), radii_m1, heights_m1, heights_m0, params, vectorlength, 1);

figure
subplot(2,3,1)
histogram(data_stage1, 'Normalization', 'probability', 'BinWidth', 0.5);
title('Histogram of Stage 1 Data');
pub_fig;
ylim([0 0.4]);
xlim([0 20]);
xlabel('Heights (10^{-6} m)')
ylabel('Frequency')

subplot(2,3,4)
[Values_s, edges] = histcounts(data_stage1, linspace(0, 25, 30), 'Normalization', 'cdf');
centers = (edges(1:end-1)+edges(2:end))/2;
plot(Values_s, centers, 'r-', 'LineWidth', 2.5)
pub_fig;
xticks([0 0.2 0.4 0.6 0.8 1]);
hold on
x = linspace(0, Values_s(1), 10);
y = [zeros(1,length(x)-1) centers(1)];
plot(x,y,'r-', 'LineWidth', 2.5)
set(gca,'ydir','reverse');
hold on
norm_heights_m1 = max_model_1 - Model_H1;
norm_heights_m1 = norm_heights_m1(norm_heights_m1 < 26);
[Values_m, edges] = histcounts(norm_heights_m1, linspace(0, 25, 30), 'Normalization', 'cdf');
centers = (edges(1:end-1)+ edges(2:end))/2;
plot(Values_m, centers, 'b-', 'LineWidth', 2.5);
title('Stage 1 Bearing Area Curves: λ = 0.8');
pub_fig;
xticks([0 0.2 0.4 0.6 0.8 1]);
x = linspace(0, Values_m(1), 10);
y = [zeros(1,length(x)-1) centers(1)];
plot(x,y,'b-', 'LineWidth', 2.5)
ylabel('Heights (10^{-6} m)')
xlabel('Quantile')
pub_fig;
set(gca,'ydir','reverse');
legend('Data', '', 'Model surface')

subplot(2,3,2)
histogram(max_model_1 - Model_H1, 'Normalization', 'probability', 'BinWidth', 0.5);
title('Histogram of Model 1 Sim');
pub_fig;
ylim([0 0.4]);
xlim([0 20]);
xlabel('Heights (10^{-6} m)')
ylabel('Frequency')

for i = 2:1:length(Values_s)
    Values_s(i) = Values_s(i) - Values_s(i-1);
    Values_m(i) = Values_m(i) - Values_m(i-1);
end
KL_div_test_stage1 = KLDiv(Values_s, Values_m);

subplot(2,3,3)
pos = get(subplot(2,3,3), 'Position');
Statistics = {'Mean';'Std dev'; 'Sa value'; 'Max'; 'Min';'Mean P2P'; '(KS & KL)'};
Stage_1 = [mean_stage1; stdev_stage1; Sa_stage1; max_stage1; min_stage1; P2P_stage1; p_value_1];
Model_1 = [mean_model_1; stdev_model_1; Sa_model_1; max_model_1; min_model_1; P2P_model_1; KL_div_test_stage1];
T = table(Stage_1, Model_1, 'RowNames', Statistics);
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,...
    'RowName',T.Properties.RowNames,'Units', 'Normalized','ColumnWidth', {70}, 'Position', pos);
pub_fig;
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('Surface Characteristics');
xlabel('All units in µm');

subplot(2,3,6)
pos = get(subplot(2,3,6), 'Position');
Model_Parameters = {'Force (N)'; 'alpha';'beta'; 'kappa 1'; 'kappa 2'; 'sigma0'; 'λ2'; 'MRR (um/min)'};
Values_1 = [Force(1); alpha(1); beta(1); kappa_1(1); kappa_2(1); sigma(1); Fiedler(1); MRR(1)];
T = table(Values_1, 'RowNames', Model_Parameters);
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,...
    'RowName',T.Properties.RowNames,'Units', 'Normalized','ColumnWidth', {70}, 'Position', pos);
pub_fig;
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('Model Parameters - Stage 1');


%% Stage 1 to Stage 2 Parameter Fitting 
%%%%%%% Asperity network %%%%%%%%%

MRR(2) = 0;

Time = 300; 
k = 2; dt = 0.2;
a_0 = 2.5;
Horizon_12 = 0:dt:Time;
Horizon = Horizon_12;
Horizon_12 = Horizon_12 + 300;

T = zeros(N,length(Horizon));  T(:,1) = Temp;
s = zeros(N,length(Horizon));  s(:,1) = sigma(2);
h = zeros(N,length(Horizon));  h(:,1) = heights_m1;
r = zeros(N,length(Horizon));  r(:,1) = radii_m1;
F = zeros(N,length(Horizon));   % Force carried by asperities
Rho = zeros(N,1); 
diff_heights = zeros(N,1);

st_d_crit = zeros(length(Horizon),1);
st_active_nodes = zeros(length(Horizon),1);
vector = 0: params{2}: max(heights_m0);
vectorlength = length(vector);

for t = dt:dt:Time
    
    d_crit = solve_for_d(Force(2),r(:,k-1),h(:,k-1),params,vectorlength);
    active_nodes = find(h(:,k-1) > d_crit); % load bearing asperities
    F(active_nodes, k) = (2/3)*E* sqrt(abs(r(active_nodes,k-1))).*((abs(h(active_nodes,k-1) - d_crit.*ones(length(active_nodes),1))).^(3/2));   % force carried by each node. 
    
    % Temp dynamics
    T(:,k) = T(:,k-1) + dt*(-alpha(2)*(T(:,k-1)- Temp) + kappa_1(2)*F(:, k-1) + kappa_2(2)*speed(2)); 
    
    % Sigma dynamics
    s(:,k) = s(:,k-1) + (beta(2).*(T(:,k) - T(:,k-1))./Temp); 
    
    % Differential height dynamics
    diff_heights(:,1) = (dt.*s(:,k-1).*T(:,k-1).*F(:,k-1));
    h(:,k) = h(:,k-1) - diff_heights;    
    
    % Radii update and MRR computation
    Rho(:,1) = h(:,k-1)./(4*r(:,k-1));
    MR_coeff_eta = exp(-(h(:,k-1)/max(heights_m0).^2)).*exp(-(a_0./r(:,k-1)).^0.0833);
    r(:,k) = r(:,k-1) + (MR_coeff_eta.*Rho(:,1).*(diff_heights));
    MRR(2) = MRR(2) + sum(pi.*(diff_heights.^2).*h(:,k-1).*MR_coeff_eta.*(1 - MR_coeff_eta)/4);
    
    % store
    st_d_crit(k-1) = d_crit;
    st_active_nodes(k-1) = length(active_nodes);
    k = k+1;
end

heights_m2 = h(:,length(Horizon));
radii_m2 = r(:,length(Horizon));
MRR(2) = MRR(2) / (250*250*Time/60);

% height, sigma and Temp change plots
figure
subplot(1,3,1)
h_12 = h(active_nodes,:);
plot(Horizon,h(active_nodes,:),'LineWidth',2)
pub_fig;
hold on
title('Height change dynamics of active asperities')
xlabel('Time (sec)')
ylabel('Height (um)')

subplot(1,3,2)
plot(Horizon,r(active_nodes,:),'LineWidth',2)
pub_fig;
hold on
title('Radii change of asperities')
xlabel('Time (sec)')
ylabel('Radii (um)')

subplot(1,3,3)
Temp_12 = T(active_nodes,:);
plot(Horizon,T(active_nodes,:),'LineWidth',2)
pub_fig;
hold on
title('Temp variation: Stage 1 - 2')
xlabel('Time (sec)')
ylabel('Temp (degree C)')


%% Stage 2 - Data and Model Comparison - Plotting

Table2 = readtable('matrix_stage2.csv');
data_stage2 = table2array(Table2);
data_p2p_stage2 = reshape(data_stage2, [102 32]);
data_stage2 = reshape(data_stage2, [3264 1]);
data2 = mean(data_stage2) - data_stage2;

% Stage 2 Data Statistics
max_stage2 = max(data2);
mean_stage2 = mean(data2);
stdev_stage2 = std(data2);
min_stage2 = min(data2);
Hurst_stage2 = Gen_hurst(data2);
Sa_stage2 = 0;
for i = 1: 3264
    Sa_stage2 = Sa_stage2 + ((1/3264)*(abs(data_stage2(i) - mean_stage2)));
end
P2P_stage2 = (sum(maxk(data_stage2,320,1),'all') - sum(mink(data_stage2,320,1), 'all'))/320;

% Model 2 Statistics
max_model_2 = max(heights_m2);
min_model_2 = min(heights_m2);
[Model_H2, mean_model_2, stdev_model_2, Sa_model_2, Hurst_model_2, KS_test_2, p_value_2] = surface_roughness(radii_m2, initial_asperity_radii, heights_m2, N, data_stage2, packing_density);
P2P_model_2 = (sum(maxk(heights_m2,20,1),'all') - sum(mink(heights_m2,20,1), 'all'))/20;

vector = 0: params{2}: max(heights_m0);
vectorlength = length(vector);
[AdjMatrix_2, Fiedler(2), E2, bins_CC2] = graph_evolution(x_pos, y_pos, packing_density, Force(2), radii_m2, heights_m2, heights_m0, params, vectorlength, 2);

figure
subplot(2,3,1)
histogram(data_stage2, 'Normalization', 'probability', 'BinWidth', 0.5);
title('Histogram of Stage 2 Data');
pub_fig;
ylim([0 0.4]);
xlim([0 20]);
xlabel('Heights (10^{-6} m)')
ylabel('Frequency')

subplot(2,3,4)
[Values_s, edges] = histcounts(data_stage2, linspace(0, 20, 30), 'Normalization', 'cdf');
centers = (edges(1:end-1)+edges(2:end))/2;
plot(Values_s, centers, 'r-', 'LineWidth', 2.5)
pub_fig;
xticks([0 0.2 0.4 0.6 0.8 1]);
hold on
x = linspace(0, Values_s(1), 10);
y = [zeros(1,length(x)-1) centers(1)];
plot(x,y,'r-', 'LineWidth', 2.5)
set(gca,'ydir','reverse');
hold on
norm_heights_m2 = max_model_2 - Model_H2;
norm_heights_m2 = norm_heights_m2(norm_heights_m2 < 21);
[Values_m, edges] = histcounts(norm_heights_m2, linspace(0, 20, 30), 'Normalization', 'cdf');
centers = (edges(1:end-1)+ edges(2:end))/2;
plot(Values_m, centers, 'b-', 'LineWidth', 2.5);
title('Bearing Area Curves - Stage 2');
pub_fig;
xticks([0 0.2 0.4 0.6 0.8 1]);
x = linspace(0, Values_m(1), 10);
y = [zeros(1,length(x)-1) centers(1)];
plot(x,y,'b-', 'LineWidth', 2.5)
ylabel('Heights (10^{-6} m)')
xlabel('Quantile')
pub_fig;
set(gca,'ydir','reverse');
legend('Data', '', 'Model surface')

subplot(2,3,2)
histogram(max_model_2 - Model_H2, 'Normalization', 'probability', 'BinWidth', 0.5);
title('Histogram of Model 2 Sim');
pub_fig;
ylim([0 0.4]);
xlim([0 20]);
xlabel('Heights (10^{-6} m)')
ylabel('Frequency')

for i = 2:1:length(Values_s)
    Values_s(i) = Values_s(i) - Values_s(i-1);
    Values_m(i) = Values_m(i) - Values_m(i-1);
end
KL_div_test_stage2 = KLDiv(Values_s, Values_m);

subplot(2,3,3)
pos = get(subplot(2,3,3), 'Position');
Statistics = {'Mean';'Std dev'; 'Sa value'; 'Max'; 'Min';'Mean P2P'; '(KS & KL)'};
Stage_2 = [mean_stage2; stdev_stage2; Sa_stage2; max_stage2; min_stage2; P2P_stage2; p_value_2];
Model_2 = [mean_model_2; stdev_model_2; Sa_model_2; max_model_2; min_model_2; P2P_model_2; KL_div_test_stage2];
T = table(Stage_2, Model_2, 'RowNames', Statistics);
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,...
    'RowName',T.Properties.RowNames,'Units', 'Normalized','ColumnWidth', {70}, 'Position', pos);
pub_fig;
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('Surface Characteristics');
xlabel('All units in µm');

subplot(2,3,6)
pos = get(subplot(2,3,6), 'Position');
Model_Parameters = {'Force (N)'; 'alpha';'beta'; 'kappa 1'; 'kappa 2'; 'sigma0'; 'λ2'; 'MRR (um/min)'};
Values_1 = [Force(2); alpha(2); beta(2); kappa_1(2); kappa_2(2); sigma(2); Fiedler(2); MRR(2)];
T = table(Values_1, 'RowNames', Model_Parameters);
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,...
    'RowName',T.Properties.RowNames,'Units', 'Normalized','ColumnWidth', {70}, 'Position', pos);
pub_fig;
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('Model Parameters - Stage 2');


%% Stage 2 to Stage 3 Parameter Fitting 
%%%%%%% Asperity network %%%%%%%%%

MRR(3) = 0;

Time = 100; 
k = 2; dt = 0.2;
a_0 = 2.5;
Horizon_23 = 0:dt:Time;
Horizon = Horizon_23;
Horizon_23 = Horizon_23 + 600;

T = zeros(N,length(Horizon));  T(:,1) = Temp;
s = zeros(N,length(Horizon));  s(:,1) = sigma(3);
h = zeros(N,length(Horizon));  h(:,1) = heights_m2;
r = zeros(N,length(Horizon));  r(:,1) = radii_m2;
F = zeros(N,length(Horizon));                   % Force carried by asperities
Rho = zeros(N,1); 
diff_heights = zeros(N,1);

st_d_crit = zeros(length(Horizon),1);
st_active_nodes = zeros(length(Horizon),1);
vector = 0: params{2}: max(heights_m0);
vectorlength = length(vector);

for t = dt:dt:Time
    
    d_crit = solve_for_d(Force(3),r(:,k-1),h(:,k-1),params,vectorlength);
    active_nodes = find(h(:,k-1) > d_crit); % load bearing asperities
    F(active_nodes, k) = (2/3)*E* sqrt(abs(r(active_nodes,k-1))).*((abs(h(active_nodes,k-1) - d_crit.*ones(length(active_nodes),1))).^(3/2));   % force carried by each node. 
    
    % Temp dynamics
    T(:,k) = T(:,k-1) + dt*(-alpha(3)*(T(:,k-1)- Temp) + kappa_1(3)*F(:, k-1) + kappa_2(3)*speed(3)); 
    
    % Sigma dynamics
    s(:,k) = s(:,k-1) + (beta(3).*(T(:,k) - T(:,k-1))./Temp); 
    
    % Differential height dynamics
    diff_heights(:,1) = (dt.*s(:,k-1).*T(:,k-1).*F(:,k-1));
    h(:,k) = h(:,k-1) - diff_heights;    
    
    % Radii update and MRR computation
    Rho(:,1) = h(:,k-1)./(4*r(:,k-1));
    MR_coeff_eta = exp(-(h(:,k-1)/max(heights_m0).^2)).*exp(-(a_0./r(:,k-1)).^0.0833);
    r(:,k) = r(:,k-1) + (MR_coeff_eta.*Rho(:,1).*(diff_heights));
    MRR(3) = MRR(3) + sum(pi.*(diff_heights.^2).*h(:,k-1).*MR_coeff_eta.*(1 - MR_coeff_eta)/4);
    
    % store
    st_d_crit(k-1) = d_crit;
    st_active_nodes(k-1) = length(active_nodes);
    k = k+1;
end

heights_m3 = h(:,length(Horizon));
radii_m3 = r(:,length(Horizon));
MRR(3) = MRR(3) / (250*250*Time/60);

% height, sigma and Temp change plots
figure
subplot(1,3,1)
h_23 = h(active_nodes,:);
plot(Horizon,h(active_nodes,:),'LineWidth',2)
pub_fig;
hold on
title('Height change dynamics of active asperities')
xlabel('Time (sec)')
ylabel('Height (um)')

subplot(1,3,2)
plot(Horizon,r(active_nodes,:),'LineWidth',2)
pub_fig;
hold on
title('Radii change of asperities')
xlabel('Time (sec)')
ylabel('Radii (um)')

subplot(1,3,3) 
Temp_23 = T(active_nodes,:);
plot(Horizon,T(active_nodes,:),'LineWidth',2)
pub_fig;
hold on
title('Temp variation: Stage 2 - 3')
xlabel('Time (sec)')
ylabel('Temp (degree C)')

%% Stage 3 - Data and Model Comparison - Plotting

Table3 = readtable('matrix_stage3.csv');
data_stage3 = table2array(Table3);
data_p2p_stage3 = reshape(data_stage3, [102 32]);
data_stage3 = reshape(data_stage3, [3264 1]);
data3 = mean(data_stage3) - data_stage3;

% Stage 3 Data Statistics
max_stage3 = max(data3);
mean_stage3 = mean(data3);
stdev_stage3 = std(data3);
min_stage3 = min(data3);
Hurst_stage3 = Gen_hurst(data3);
Sa_stage3 = 0;
for i = 1: 3264
    Sa_stage3 = Sa_stage3 + ((1/3264)*(abs(data_stage3(i) - mean_stage3)));
end
P2P_stage3 = (sum(maxk(data_stage3,320,1),'all') - sum(mink(data_stage3,320,1), 'all'))/320;

% Model 3 Statistics
max_model_3 = max(heights_m3);
min_model_3 = min(heights_m3);
[Model_H3, mean_model_3, stdev_model_3, Sa_model_3, Hurst_model_3, KS_test_3, p_value_3] = surface_roughness(radii_m3, initial_asperity_radii, heights_m3, N, data_stage3, packing_density);
P2P_model_3 = (sum(maxk(heights_m3,20,1),'all') - sum(mink(heights_m3,20,1), 'all'))/20;

vector = 0: params{2}: max(heights_m0);
vectorlength = length(vector);
[AdjMatrix_3, Fiedler(3), E3, bins_CC3] = graph_evolution(x_pos, y_pos, packing_density, Force(3), radii_m3, heights_m3, heights_m0, params, vectorlength, 3);

figure
subplot(2,3,1)
histogram(data_stage2, 'Normalization', 'probability', 'BinWidth', 0.5);
title('Histogram of Stage 3 Data');
pub_fig;
ylim([0 0.7]);
xlim([0 15]);
xlabel('Heights (10^{-6} m)')
ylabel('Frequency')

subplot(2,3,4)
[Values_s, edges] = histcounts(data_stage3, linspace(0, 15, 30), 'Normalization', 'cdf');
centers = (edges(1:end-1)+edges(2:end))/2;
plot(Values_s, centers, 'r-', 'LineWidth', 2.5)
pub_fig;
xticks([0 0.2 0.4 0.6 0.8 1]);
hold on
x = linspace(0, Values_s(1), 10);
y = [zeros(1,length(x)-1) centers(1)];
plot(x,y,'r-', 'LineWidth', 2.5)
set(gca,'ydir','reverse');
hold on
norm_heights_m3 = max_model_3 - Model_H3;
norm_heights_m3 = norm_heights_m3(norm_heights_m3 < 16);
[Values_m, edges] = histcounts(norm_heights_m3, linspace(0, 15, 30), 'Normalization', 'cdf');
centers = (edges(1:end-1)+ edges(2:end))/2;
plot(Values_m, centers, 'b-', 'LineWidth', 2.5);
title('Bearing Area Curves - Stage 3');
pub_fig;
xticks([0 0.2 0.4 0.6 0.8 1]);
x = linspace(0, Values_m(1), 10);
y = [zeros(1,length(x)-1) centers(1)];
plot(x,y,'b-', 'LineWidth', 2.5)
ylabel('Heights (10^{-6} m)')
xlabel('Quantile')
pub_fig;
set(gca,'ydir','reverse');
legend('Data', '', 'Model surface')

subplot(2,3,2)
histogram(max_model_3 - Model_H3, 'Normalization', 'probability', 'BinWidth', 0.5);
title('Histogram of Model 3 Sim');
pub_fig;
ylim([0 0.7]);
xlim([0 15]);
xlabel('Heights (10^{-6} m)')
ylabel('Frequency')

for i = 2:1:length(Values_s)
    Values_s(i) = Values_s(i) - Values_s(i-1);
    Values_m(i) = Values_m(i) - Values_m(i-1);
end
KL_div_test_stage3 = KLDiv(Values_s, Values_m);

subplot(2,3,3)
pos = get(subplot(2,3,3), 'Position');
Statistics = {'Mean';'Std dev'; 'Sa value'; 'Max'; 'Min';'Mean P2P'; '(KS & KL)'};
Stage_3 = [mean_stage3; stdev_stage3; Sa_stage3; max_stage3; min_stage3; P2P_stage3; p_value_3];
Model_3 = [mean_model_3; stdev_model_3; Sa_model_3; max_model_3; min_model_3; P2P_model_3; KL_div_test_stage3];
T = table(Stage_3, Model_3, 'RowNames', Statistics);
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,...
    'RowName',T.Properties.RowNames,'Units', 'Normalized','ColumnWidth', {70}, 'Position', pos);
pub_fig;
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('Surface Characteristics');
xlabel('All units in µm');

subplot(2,3,6)
pos = get(subplot(2,3,6), 'Position');
Model_Parameters = {'Force (N)'; 'alpha';'beta'; 'kappa 1'; 'kappa 2'; 'sigma0'; 'λ2'; 'MRR (um/min)'};
Values_1 = [Force(3); alpha(3); beta(3); kappa_1(3); kappa_2(3); sigma(3); Fiedler(3); MRR(3)];
T = table(Values_1, 'RowNames', Model_Parameters);
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,...
    'RowName',T.Properties.RowNames,'Units', 'Normalized','ColumnWidth', {70}, 'Position', pos);
pub_fig;
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('Model Parameters - Stage 3');



%% Stage 3 to Stage 4 Parameter Fitting 
%%%%%%% Asperity network %%%%%%%%%

MRR(4) = 0;

Time = 450; 
k = 2; dt = 0.2;
a_0 = 2.5;
Horizon_34 = 0:dt:Time;
Horizon = Horizon_34;
Horizon_34 = Horizon_34 + 1200;

T = zeros(N,length(Horizon));  T(:,1) = Temp;
s = zeros(N,length(Horizon));  s(:,1) = sigma(4);
h = zeros(N,length(Horizon));  h(:,1) = heights_m3;
r = zeros(N,length(Horizon));  r(:,1) = radii_m3;
F = zeros(N,length(Horizon));                   % Force carried by asperities
Rho = zeros(N,1); 
diff_heights = zeros(N,1);

st_d_crit = zeros(length(Horizon),1);
st_active_nodes = zeros(length(Horizon),1);
vector = 0: params{2}: max(heights_m0);
vectorlength = length(vector);

for t = dt:dt:Time
    
    d_crit = solve_for_d(Force(4),r(:,k-1),h(:,k-1),params,vectorlength);
    active_nodes = find(h(:,k-1) > d_crit); % load bearing asperities
    F(active_nodes, k) = (2/3)*E* sqrt(abs(r(active_nodes,k-1))).*((abs(h(active_nodes,k-1) - d_crit.*ones(length(active_nodes),1))).^(3/2));   % force carried by each node. 
    
    % Temp dynamics
    T(:,k) = T(:,k-1) + dt*(-alpha(4)*(T(:,k-1)- Temp) + kappa_1(4)*F(:, k-1) + kappa_2(4)*speed(4)); 
    
    % Sigma dynamics
    s(:,k) = s(:,k-1) + (beta(4).*(T(:,k) - T(:,k-1))./Temp); 
    
    % Differential height dynamics
    diff_heights(:,1) = (dt.*s(:,k-1).*T(:,k-1).*F(:,k-1));
    h(:,k) = h(:,k-1) - diff_heights;    
    
    % Radii update and MRR computation
    Rho(:,1) = h(:,k-1)./(4*r(:,k-1));
    MR_coeff_eta = exp(-(h(:,k-1)/max(heights_m0).^2)).*exp(-(a_0./r(:,k-1)).^0.0833);
    r(:,k) = r(:,k-1) + (MR_coeff_eta.*Rho(:,1).*(diff_heights));
    MRR(4) = MRR(4) + sum(pi.*(diff_heights.^2).*h(:,k-1).*MR_coeff_eta.*(1 - MR_coeff_eta)/4);
    
    % store
    st_d_crit(k-1) = d_crit;
    st_active_nodes(k-1) = length(active_nodes);
    k = k+1;
end

heights_m4 = h(:,length(Horizon));
radii_m4 = r(:,length(Horizon));
MRR(4) = MRR(4) / (250*250*Time/60);

% height, sigma and Temp change plots
figure
subplot(1,3,1)
h_34 = h(active_nodes,:);
plot(Horizon,h(active_nodes,:),'LineWidth',2)
pub_fig;
hold on
title('Height change dynamics of active asperities')
xlabel('Time (sec)')
ylabel('Height (um)')

subplot(1,3,2)
plot(Horizon,r(active_nodes,:),'LineWidth',2)
pub_fig;
hold on
title('Radii change of asperities')
xlabel('Time (sec)')
ylabel('Radii (um)')

subplot(1,3,3)
Temp_34 = T(active_nodes,:);
plot(Horizon,T(active_nodes,:),'LineWidth',2)
pub_fig;
hold on
title('Temp variation: Stage 3 - 4')
xlabel('Time (sec)')
ylabel('Temp (degree C)')

%% Stage 4 - Data and Model Comparison - Plotting

Table4 = readtable('matrix_stage4.csv');
data_stage4 = table2array(Table4);
data_p2p_stage4 = reshape(data_stage4, [102 32]);
data_stage4 = reshape(data_stage4, [3264 1]);
data4 = mean(data_stage4) - data_stage4;

% Stage 4 Data Statistics
max_stage4 = max(data4);
mean_stage4 = mean(data4);
stdev_stage4 = std(data4);
min_stage4 = min(data4);
Hurst_stage4 = Gen_hurst(data4);
Sa_stage4 = 0;
for i = 1: 3264
    Sa_stage4 = Sa_stage4 + ((1/3264)*(abs(data_stage4(i) - mean_stage4)));
end
P2P_stage4 = (sum(maxk(data_stage4,320,1),'all') - sum(mink(data_stage4,320,1), 'all'))/320;

% Model 4 Statistics
max_model_4 = max(heights_m4);
min_model_4 = min(heights_m4);
[Model_H4, mean_model_4, stdev_model_4, Sa_model_4, Hurst_model_4, KS_test_4, p_value_4] = surface_roughness(radii_m4, initial_asperity_radii, heights_m4, N, data_stage4, packing_density);
P2P_model_4 = (sum(maxk(heights_m4,20,1),'all') - sum(mink(heights_m4,20,1), 'all'))/20;

vector = 0: params{2}: max(heights_m0);
vectorlength = length(vector);
[AdjMatrix_4, Fiedler(4), E4, bins_CC4] = graph_evolution(x_pos, y_pos, packing_density, Force(4), radii_m4, heights_m4, heights_m0, params, vectorlength, 4);

figure
subplot(2,3,1)
histogram(data_stage4, 'Normalization', 'probability', 'BinWidth', 0.5);
title('Histogram of Stage 4 Data');
pub_fig;
ylim([0 0.8]);
xlim([0 15]);
xlabel('Heights (10^{-6} m)')
ylabel('Frequency')

subplot(2,3,4)
[Values_s, edges] = histcounts(data_stage4, linspace(0, 15, 30), 'Normalization', 'cdf');
centers = (edges(1:end-1)+edges(2:end))/2;
plot(Values_s, centers, 'r-', 'LineWidth', 2.5)
pub_fig;
xticks([0 0.2 0.4 0.6 0.8 1]);
hold on
x = linspace(0, Values_s(1), 10);
y = [zeros(1,length(x)-1) centers(1)];
plot(x,y,'r-', 'LineWidth', 2.5)
set(gca,'ydir','reverse');
hold on
norm_heights_m4 = max_model_4 - Model_H4;
norm_heights_m4 = norm_heights_m4(norm_heights_m4 < 15);
[Values_m, edges] = histcounts(norm_heights_m4, linspace(0, 15, 30), 'Normalization', 'cdf');
centers = (edges(1:end-1)+ edges(2:end))/2;
plot(Values_m, centers, 'b-', 'LineWidth', 2.5);
title('Bearing Area Curves - Stage 4');
pub_fig;
xticks([0 0.2 0.4 0.6 0.8 1]);
x = linspace(0, Values_m(1), 10);
y = [zeros(1,length(x)-1) centers(1)];
plot(x,y,'b-', 'LineWidth', 2.5)
ylabel('Heights (10^{-6} m)')
xlabel('Quantile')
pub_fig;
set(gca,'ydir','reverse');
legend('Data', '', 'Model surface')

subplot(2,3,2)
histogram(max_model_4 - Model_H4, 'Normalization', 'probability', 'BinWidth', 0.5);
title('Histogram of Model 4 Sim');
pub_fig;
ylim([0 0.8]);
xlim([0 15]);
xlabel('Heights (10^{-6} m)')
ylabel('Frequency')

for i = 2:1:length(Values_s)
    Values_s(i) = Values_s(i) - Values_s(i-1);
    Values_m(i) = Values_m(i) - Values_m(i-1);
end
KL_div_test_stage4 = KLDiv(Values_s, Values_m);

subplot(2,3,3)
pos = get(subplot(2,3,3), 'Position');
Statistics = {'Mean';'Std dev'; 'Sa value'; 'Max'; 'Min';'Mean P2P'; '(KS & KL)'};
Stage_4 = [mean_stage4; stdev_stage4; Sa_stage4; max_stage4; min_stage4; P2P_stage4; p_value_4];
Model_4 = [mean_model_4; stdev_model_4; Sa_model_4; max_model_4; min_model_4; P2P_model_4; KL_div_test_stage4];
T = table(Stage_4, Model_4, 'RowNames', Statistics);
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,...
    'RowName',T.Properties.RowNames,'Units', 'Normalized','ColumnWidth', {70}, 'Position', pos);
pub_fig;
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('Surface Characteristics');
xlabel('All units in µm');

subplot(2,3,6)
pos = get(subplot(2,3,6), 'Position');
Model_Parameters = {'Force (N)'; 'alpha';'beta'; 'kappa 1'; 'kappa 2'; 'sigma0'; 'λ2'; 'MRR (um/min)'};
Values_1 = [Force(4); alpha(4); beta(4); kappa_1(4); kappa_2(4); sigma(4); Fiedler(4); MRR(4)];
T = table(Values_1, 'RowNames', Model_Parameters);
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,...
    'RowName',T.Properties.RowNames,'Units', 'Normalized','ColumnWidth', {70}, 'Position', pos);
pub_fig;
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('Model Parameters - Stage 4');


%% Stage 4 to Stage 5 Parameter Fitting 
%%%%%%% Asperity network %%%%%%%%%

MRR(5) = 0;

Time = 600; 
k = 2; dt = 0.2;
a_0 = 2.5;
Horizon_45 = 0:dt:Time;
Horizon = Horizon_45;
Horizon_45 = Horizon_45 + 1800;

T = zeros(N,length(Horizon));  T(:,1) = Temp;
s = zeros(N,length(Horizon));  s(:,1) = sigma(5);
h = zeros(N,length(Horizon));  h(:,1) = heights_m4;
r = zeros(N,length(Horizon));  r(:,1) = radii_m4;
F = zeros(N,length(Horizon));                   % Force carried by asperities
Rho = zeros(N,1); 
diff_heights = zeros(N,1);

st_d_crit = zeros(length(Horizon),1);
st_active_nodes = zeros(length(Horizon),1);
vector = 0: params{2}: max(heights_m0);
vectorlength = length(vector);

for t = dt:dt:Time
    
    d_crit = solve_for_d(Force(5),r(:,k-1),h(:,k-1),params,vectorlength);
    active_nodes = find(h(:,k-1) > d_crit); % load bearing asperities
    F(active_nodes, k) = (2/3)*E* sqrt(abs(r(active_nodes,k-1))).*((abs(h(active_nodes,k-1) - d_crit.*ones(length(active_nodes),1))).^(3/2));   % force carried by each node. 
    
    % Temp dynamics
    T(:,k) = T(:,k-1) + dt*(-alpha(5)*(T(:,k-1)- Temp) + kappa_1(5)*F(:, k-1) + kappa_2(5)*speed(5)); 
    
    % Sigma dynamics
    s(:,k) = s(:,k-1) + (beta(5).*(T(:,k) - T(:,k-1))./Temp); 
    
    % Differential height dynamics
    diff_heights(:,1) = (dt.*s(:,k-1).*T(:,k-1).*F(:,k-1));
    h(:,k) = h(:,k-1) - diff_heights;    
    
    % Radii update and MRR computation
    Rho(:,1) = h(:,k-1)./(4*r(:,k-1));
    MR_coeff_eta = exp(-(h(:,k-1)/max(heights_m0).^2)).*exp(-(a_0./r(:,k-1)).^0.0833);
    r(:,k) = r(:,k-1) + (MR_coeff_eta.*Rho(:,1).*(diff_heights));
    MRR(5) = MRR(5) + sum(pi.*(diff_heights.^2).*h(:,k-1).*MR_coeff_eta.*(1 - MR_coeff_eta)/4);
    
    % store
    st_d_crit(k-1) = d_crit;
    st_active_nodes(k-1) = length(active_nodes);
    k = k+1;
end

heights_m5 = h(:,length(Horizon));
radii_m5 = r(:,length(Horizon));
MRR(5) = MRR(5) / (250*250*Time/60);

% height, sigma and Temp change plots
figure
subplot(1,3,1)
h_45 = h(active_nodes,:);
plot(Horizon,h(active_nodes,:),'LineWidth',2)
pub_fig;
hold on
title('Height change dynamics of active asperities')
xlabel('Time (sec)')
ylabel('Height (um)')

% subplot(1,3,2)
% plot(Horizon,r(active_nodes,:),'LineWidth',2)
% pub_fig;
% hold on
% title('Radii change of asperities')
% xlabel('Time (sec)')
% ylabel('Radii (um)')

subplot(1,3,2)
Temp_45 = T(active_nodes,:);
plot(Horizon,T(active_nodes,:),'LineWidth',2)
pub_fig;
hold on
title('Temp variation: Stage 4 - 5')
xlabel('Time (sec)')
ylabel('Temp (degree C)')


%% Stage 5 - Data and Model Comparison - Plotting

Table5 = readtable('matrix_stage5.csv');
data_stage5 = table2array(Table5);
data_p2p_stage5 = reshape(data_stage5, [102 32]);
data_stage5 = reshape(data_stage5, [3264 1]);
data5 = mean(data_stage5) - data_stage5;

% Stage 5 Data Statistics
max_stage5 = max(data5);
mean_stage5 = mean(data5);
stdev_stage5 = std(data5);
min_stage5 = min(data5);
Hurst_stage5 = Gen_hurst(data5);
Sa_stage5 = 0;
for i = 1: 3264
    Sa_stage5 = Sa_stage5 + ((1/3264)*(abs(data_stage5(i) - mean_stage5)));
end
P2P_stage5 = (sum(maxk(data_stage5,320,1),'all') - sum(mink(data_stage5,320,1), 'all'))/320;

% Model 5 Statistics
max_model_5 = max(heights_m5);
min_model_5 = min(heights_m5);
[Model_H5, mean_model_5, stdev_model_5, Sa_model_5, Hurst_model_5, KS_test_5, p_value_5] = surface_roughness(radii_m5, initial_asperity_radii, heights_m5, N, data_stage5, packing_density);
P2P_model_5 = (sum(maxk(heights_m5,20,1),'all') - sum(mink(heights_m5,20,1), 'all'))/20;

vector = 0: params{2}: max(heights_m0);
vectorlength = length(vector);
[AdjMatrix_5, Fiedler(5), E5, bins_CC5] = graph_evolution(x_pos, y_pos, packing_density, Force(5), radii_m5, heights_m5, heights_m0, params, vectorlength, 5);

figure
subplot(2,3,1)
histogram(data_stage5, 'Normalization', 'probability', 'BinWidth', 0.5);
title('Histogram of Stage 5 Data');
pub_fig;
ylim([0 1]);
xlim([0 5]);
xlabel('Heights (10^{-6} m)')
ylabel('Frequency')

subplot(2,3,4)
[Values_s, edges] = histcounts(data_stage5, linspace(0, 5, 30), 'Normalization', 'cdf');
centers = (edges(1:end-1)+edges(2:end))/2;
plot(Values_s, centers, 'r-', 'LineWidth', 2.5)
pub_fig;
xticks([0 0.2 0.4 0.6 0.8 1]);
hold on
x = linspace(0, Values_s(1), 10);
y = [zeros(1,length(x)-1) centers(1)];
plot(x,y,'r-', 'LineWidth', 2.5)
set(gca,'ydir','reverse');
hold on
norm_heights_m5 = max_model_5 - Model_H5;
norm_heights_m5 = norm_heights_m5(norm_heights_m5 < 5);
[Values_m, edges] = histcounts(norm_heights_m5, linspace(0, 5, 30), 'Normalization', 'cdf');
centers = (edges(1:end-1)+ edges(2:end))/2;
plot(Values_m, centers, 'b-', 'LineWidth', 2.5);
title('Bearing Area Curves - Stage 5');
pub_fig;
xticks([0 0.2 0.4 0.6 0.8 1]);
x = linspace(0, Values_m(1), 10);
y = [zeros(1,length(x)-1) centers(1)];
plot(x,y,'b-', 'LineWidth', 2.5)
ylabel('Heights (10^{-6} m)')
xlabel('Quantile')
pub_fig;
set(gca,'ydir','reverse');
legend('Data', '', 'Model surface')

subplot(2,3,2)
histogram(max_model_5 - Model_H5, 'Normalization', 'probability', 'BinWidth', 0.5);
title('Histogram of Model 5 Sim');
pub_fig;
ylim([0 1]);
xlim([0 5]);
xlabel('Heights (10^{-6} m)')
ylabel('Frequency')

for i = 2:1:length(Values_s)
    Values_s(i) = Values_s(i) - Values_s(i-1);
    Values_m(i) = Values_m(i) - Values_m(i-1);
end
KL_div_test_stage5 = KLDiv(Values_s, Values_m);

subplot(2,3,3)
pos = get(subplot(2,3,3), 'Position');
Statistics = {'Mean';'Std dev'; 'Sa value'; 'Max'; 'Min';'Mean P2P'; '(KS & KL)'};
Stage_5 = [mean_stage5; stdev_stage5; Sa_stage5; max_stage5; min_stage5; P2P_stage5; p_value_5];
Model_5 = [mean_model_5; stdev_model_5; Sa_model_5; max_model_5; min_model_5; P2P_model_5; KL_div_test_stage5];
T = table(Stage_5, Model_5, 'RowNames', Statistics);
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,...
    'RowName',T.Properties.RowNames,'Units', 'Normalized','ColumnWidth', {70}, 'Position', pos);
pub_fig;
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('Surface Characteristics');
xlabel('All units in µm');

subplot(2,3,6)
pos = get(subplot(2,3,6), 'Position');
Model_Parameters = {'Force (N)'; 'alpha';'beta'; 'kappa 1'; 'kappa 2'; 'sigma0'; 'λ2'; 'MRR (um/min)'};
Values_1 = [Force(5); alpha(5); beta(5); kappa_1(5); kappa_2(5); sigma(5); Fiedler(5); MRR(5)];
T = table(Values_1, 'RowNames', Model_Parameters);
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,...
    'RowName',T.Properties.RowNames,'Units', 'Normalized','ColumnWidth', {70}, 'Position', pos);
pub_fig;
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('Model Parameters - Stage 5');


%% Stage 5 to Stage 6 Parameter Fitting 
%%%%%%% Asperity network %%%%%%%%%

MRR(6) = 0;

Time = 900; 
k = 2; dt = 0.2;
a_0 = 2.5;
Horizon_56 = 0:dt:Time;
Horizon = Horizon_56;
Horizon_56 = Horizon_56 + 2400;

T = zeros(N,length(Horizon));  T(:,1) = Temp;
s = zeros(N,length(Horizon));  s(:,1) = sigma(6);
h = zeros(N,length(Horizon));  h(:,1) = heights_m5;
r = zeros(N,length(Horizon));  r(:,1) = radii_m5;
F = zeros(N,length(Horizon));                   % Force carried by asperities
Rho = zeros(N,1); 
diff_heights = zeros(N,1);

st_d_crit = zeros(length(Horizon),1);
st_active_nodes = zeros(length(Horizon),1);
vector = 0: params{2}: max(heights_m0);
vectorlength = length(vector);

for t = dt:dt:Time
    
    d_crit = solve_for_d(Force(6),r(:,k-1),h(:,k-1),params,vectorlength);
    active_nodes = find(h(:,k-1) > d_crit); % load bearing asperities
    F(active_nodes, k) = (2/3)*E* sqrt(abs(r(active_nodes,k-1))).*((abs(h(active_nodes,k-1) - d_crit.*ones(length(active_nodes),1))).^(3/2));   % force carried by each node. 
    
    % Temp dynamics
    T(:,k) = T(:,k-1) + dt*(-alpha(6)*(T(:,k-1)- Temp) + kappa_1(6)*F(:, k-1) + kappa_2(6)*speed(6)); 
    
    % Sigma dynamics
    s(:,k) = s(:,k-1) + (beta(6).*(T(:,k) - T(:,k-1))./Temp); 
    
    % Differential height dynamics
    diff_heights(:,1) = (dt.*s(:,k-1).*T(:,k-1).*F(:,k-1));
    h(:,k) = h(:,k-1) - diff_heights;    
    
    % Radii update and MRR computation
    Rho(:,1) = h(:,k-1)./(4*r(:,k-1));
    MR_coeff_eta = exp(-(h(:,k-1)/max(heights_m0).^2)).*exp(-(a_0./r(:,k-1)).^0.0833);
    r(:,k) = r(:,k-1) + (MR_coeff_eta.*Rho(:,1).*(diff_heights));
    MRR(6) = MRR(6) + sum(pi.*(diff_heights.^2).*h(:,k-1).*MR_coeff_eta.*(1 - MR_coeff_eta)/4);
    
    % store
    st_d_crit(k-1) = d_crit;
    st_active_nodes(k-1) = length(active_nodes);
    k = k+1;
end

heights_m6 = h(:,length(Horizon));
radii_m6 = r(:,length(Horizon));
MRR(6) = MRR(6) / (250*250*Time/60);

% height, sigma and Temp change plots
figure
subplot(1,3,1)
h_56 = h(active_nodes,:);
plot(Horizon,h(active_nodes,:),'LineWidth',2)
pub_fig;
hold on
title('Height change dynamics of active asperities')
xlabel('Time (sec)')
ylabel('Height (um)')

% subplot(1,3,2)
% plot(Horizon,r(active_nodes,:),'LineWidth',2)
% pub_fig;
% hold on
% title('Radii change of asperities')
% xlabel('Time (sec)')
% ylabel('Radii (um)')

subplot(1,3,2)
Temp_56 = T(active_nodes,:);
plot(Horizon,T(active_nodes,:),'LineWidth',2)
pub_fig;
hold on
title('Temp variation: Stage 5 - 6')
xlabel('Time (sec)')
ylabel('Temp (degree C)')



%% Stage 6 - Data and Model Comparison - Plotting

Table6 = readtable('matrix_stage6.csv');
data_stage6 = table2array(Table6);
data_p2p_stage6 = reshape(data_stage6, [102 32]);
data_stage6 = reshape(data_stage6, [3264 1]);
data6 = mean(data_stage6) - data_stage6;

% Stage 6 Data Statistics
max_stage6 = max(data6);
mean_stage6 = mean(data6);
stdev_stage6 = std(data6);
min_stage6 = min(data6);
Hurst_stage6 = Gen_hurst(data6);
Sa_stage6 = 0;
for i = 1: 3264
    Sa_stage6 = Sa_stage6 + ((1/3264)*(abs(data_stage6(i) - mean_stage6)));
end
P2P_stage6 = (sum(maxk(data_stage6,320,1),'all') - sum(mink(data_stage6,320,1), 'all'))/320;

% Model 6 Statistics
max_model_6 = max(heights_m6);
min_model_6 = min(heights_m6);
[Model_H6, mean_model_6, stdev_model_6, Sa_model_6, Hurst_model_6, KS_test_6, p_value_6] = surface_roughness(radii_m6, initial_asperity_radii, heights_m6, N, data_stage6, packing_density);
P2P_model_6 = (sum(maxk(heights_m6,20,1),'all') - sum(mink(heights_m6,20,1), 'all'))/20;

vector = 0: params{2}: max(heights_m0);
vectorlength = length(vector);
[AdjMatrix_6, Fiedler(6), E6, bins_CC6] = graph_evolution(x_pos, y_pos, packing_density, Force(6), radii_m6, heights_m6, heights_m0, params, vectorlength, 6);

figure
subplot(2,3,1)
histogram(data_stage6, 'Normalization', 'probability', 'BinWidth', 0.5);
title('Histogram of Stage 6 Data');
pub_fig;
ylim([0 1]);
xlim([0 5]);
xlabel('Heights (10^{-6} m)')
ylabel('Frequency')

subplot(2,3,4)
[Values_s, edges] = histcounts(data_stage6, linspace(0, 5, 30), 'Normalization', 'cdf');
centers = (edges(1:end-1)+edges(2:end))/2;
plot(Values_s, centers, 'r-', 'LineWidth', 2.5)
pub_fig;
xticks([0 0.2 0.4 0.6 0.8 1]);
hold on
x = linspace(0, Values_s(1), 10);
y = [zeros(1,length(x)-1) centers(1)];
plot(x,y,'r-', 'LineWidth', 2.5)
set(gca,'ydir','reverse');
hold on
norm_heights_m6 = max_model_6 - Model_H6;
norm_heights_m6 = norm_heights_m6(norm_heights_m6 < 5);
[Values_m, edges] = histcounts(norm_heights_m6, linspace(0, 5, 30), 'Normalization', 'cdf');
centers = (edges(1:end-1)+ edges(2:end))/2;
plot(Values_m, centers, 'b-', 'LineWidth', 2.5);
title('Bearing Area Curves - Stage 6');
pub_fig;
xticks([0 0.2 0.4 0.6 0.8 1]);
x = linspace(0, Values_m(1), 10);
y = [zeros(1,length(x)-1) centers(1)];
plot(x,y,'b-', 'LineWidth', 2.5)
ylabel('Heights (10^{-6} m)')
xlabel('Quantile')
pub_fig;
set(gca,'ydir','reverse');
legend('Data', '', 'Model surface')

subplot(2,3,2)
histogram(max_model_6 - Model_H6, 'Normalization', 'probability', 'BinWidth', 0.5);
title('Histogram of Model 6 Sim');
pub_fig;
ylim([0 1]);
xlim([0 5]);
xlabel('Heights (10^{-6} m)')
ylabel('Frequency')

for i = 2:1:length(Values_s)
    Values_s(i) = Values_s(i) - Values_s(i-1);
    Values_m(i) = Values_m(i) - Values_m(i-1);
end
KL_div_test_stage6 = KLDiv(Values_s, Values_m);

subplot(2,3,3)
pos = get(subplot(2,3,3), 'Position');
Statistics = {'Mean';'Std dev'; 'Sa value'; 'Max'; 'Min';'Mean P2P'; '(KS & KL)'};
Stage_6 = [mean_stage6; stdev_stage6; Sa_stage6; max_stage6; min_stage6; P2P_stage6; p_value_6];
Model_6 = [mean_model_6; stdev_model_6; Sa_model_6; max_model_6; min_model_6; P2P_model_6; KL_div_test_stage6];
T = table(Stage_6, Model_6, 'RowNames', Statistics);
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,...
    'RowName',T.Properties.RowNames,'Units', 'Normalized','ColumnWidth', {70}, 'Position', pos);
pub_fig;
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('Surface Characteristics');
xlabel('All units in µm');

subplot(2,3,6)
pos = get(subplot(2,3,6), 'Position');
Model_Parameters = {'Force (N)'; 'alpha';'beta'; 'kappa 1'; 'kappa 2'; 'sigma0'; 'λ2'; 'MRR (um/min)'};
Values_1 = [Force(6); alpha(6); beta(6); kappa_1(6); kappa_2(6); sigma(6); Fiedler(6); MRR(6)];
T = table(Values_1, 'RowNames', Model_Parameters);
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,...
    'RowName',T.Properties.RowNames,'Units', 'Normalized','ColumnWidth', {70}, 'Position', pos);
pub_fig;
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('Model Parameters - Stage 6');


%% Stage 6 to Stage 7 Parameter Fitting 
%%%%%%% Asperity network %%%%%%%%%

MRR(7) = 0;

Time = 900; 
k = 2; dt = 0.2;
a_0 = 2.5;
Horizon_67 = 0:dt:Time;
Horizon = Horizon_67;
Horizon_67 = Horizon_67 + 3300;

T = zeros(N,length(Horizon));  T(:,1) = Temp;
s = zeros(N,length(Horizon));  s(:,1) = sigma(7);
h = zeros(N,length(Horizon));  h(:,1) = heights_m6;
r = zeros(N,length(Horizon));  r(:,1) = radii_m6;
F = zeros(N,length(Horizon));                   % Force carried by asperities
Rho = zeros(N,1); 
diff_heights = zeros(N,1);

st_d_crit = zeros(length(Horizon),1);
st_active_nodes = zeros(length(Horizon),1);
vector = 0: params{2}: max(heights_m0);
vectorlength = length(vector);

for t = dt:dt:Time
    
    d_crit = solve_for_d(Force(7),r(:,k-1),h(:,k-1),params,vectorlength);
    active_nodes = find(h(:,k-1) > d_crit); % load bearing asperities
    F(active_nodes, k) = (2/3)*E* sqrt(abs(r(active_nodes,k-1))).*((abs(h(active_nodes,k-1) - d_crit.*ones(length(active_nodes),1))).^(3/2));   % force carried by each node. 
    
    % Temp dynamics
    T(:,k) = T(:,k-1) + dt*(-alpha(7)*(T(:,k-1)- Temp) + kappa_1(7)*F(:, k-1) + kappa_2(7)*speed(7)); 
    
    % Sigma dynamics
    s(:,k) = s(:,k-1) + (beta(7).*(T(:,k) - T(:,k-1))./Temp); 
    
    % Differential height dynamics
    diff_heights(:,1) = (dt.*s(:,k-1).*T(:,k-1).*F(:,k-1));
    h(:,k) = h(:,k-1) - diff_heights;    
    
    % Radii update and MRR computation
    Rho(:,1) = h(:,k-1)./(4*r(:,k-1));
    MR_coeff_eta = exp(-(h(:,k-1)/max(heights_m0).^2)).*exp(-(a_0./r(:,k-1)).^0.0833);
    r(:,k) = r(:,k-1) + (MR_coeff_eta.*Rho(:,1).*(diff_heights));
    MRR(7) = MRR(7) + sum(pi.*(diff_heights.^2).*h(:,k-1).*MR_coeff_eta.*(1 - MR_coeff_eta)/4);
    
    % store
    st_d_crit(k-1) = d_crit;
    st_active_nodes(k-1) = length(active_nodes);
    k = k+1;
end

heights_m7 = h(:,length(Horizon));
radii_m7 = r(:,length(Horizon));
MRR(7) = MRR(7) / (250*250*Time/60);

% height, sigma and Temp change plots
figure
subplot(1,3,1)
h_67 = h(active_nodes,:);
plot(Horizon,h(active_nodes,:),'LineWidth',2)
pub_fig;
hold on
title('Height change dynamics of active asperities')
xlabel('Time (sec)')
ylabel('Height (um)')

% subplot(1,3,2)
% plot(Horizon,r(active_nodes,:),'LineWidth',2)
% pub_fig;
% hold on
% title('Radii change of asperities')
% xlabel('Time (sec)')
% ylabel('Radii (um)')

subplot(1,3,2)
Temp_67 = T(active_nodes,:);
plot(Horizon,T(active_nodes,:),'LineWidth',2)
pub_fig;
hold on
title('Temp variation: Stage 6 - 7')
xlabel('Time (sec)')
ylabel('Temp (degree C)')


%% Stage 7 - Data and Model Comparison - Plotting

Table7 = readtable('matrix_stage7.csv');
data_stage7 = table2array(Table7);
data_p2p_stage7 = reshape(data_stage7, [102 32]);
data_stage7 = reshape(data_stage7, [3264 1]);
data7 = mean(data_stage7) - data_stage7;

% Stage 7 Data Statistics
max_stage7 = max(data7);
mean_stage7 = mean(data7);
stdev_stage7 = std(data7);
min_stage7 = min(data7);
Hurst_stage7 = Gen_hurst(data7);
Sa_stage7 = 0;
for i = 1: 3264
    Sa_stage7 = Sa_stage7 + ((1/3264)*(abs(data_stage7(i) - mean_stage7)));
end
P2P_stage7 = (sum(maxk(data_stage7,320,1),'all') - sum(mink(data_stage7,320,1), 'all'))/320;

% Model 7 Statistics
max_model_7 = max(heights_m7);
min_model_7 = min(heights_m7);
[Model_H7, mean_model_7, stdev_model_7, Sa_model_7, Hurst_model_7, KS_test_7, p_value_7] = surface_roughness(radii_m7, initial_asperity_radii, heights_m7, N, data_stage7, packing_density);
P2P_model_7 = (sum(maxk(heights_m7,20,1),'all') - sum(mink(heights_m7,20,1), 'all'))/20;

vector = 0: params{2}: max(heights_m0);
vectorlength = length(vector);
[AdjMatrix_7, Fiedler(7), E7, bins_CC7] = graph_evolution(x_pos, y_pos, packing_density, Force(7), radii_m7, heights_m7, heights_m0, params, vectorlength, 7);

figure
subplot(2,3,1)
histogram(data_stage7, 'Normalization', 'probability', 'BinWidth', 0.5);
title('Histogram of Stage 7 Data');
pub_fig;
ylim([0 1]);
xlim([0 5]);
xlabel('Heights (10^{-6} m)')
ylabel('Frequency')

subplot(2,3,4)
[Values_s, edges] = histcounts(data_stage7, linspace(0, 5, 30), 'Normalization', 'cdf');
centers = (edges(1:end-1)+edges(2:end))/2;
plot(Values_s, centers, 'r-', 'LineWidth', 2.5)
pub_fig;
xticks([0 0.2 0.4 0.6 0.8 1]);
hold on
x = linspace(0, Values_s(1), 10);
y = [zeros(1,length(x)-1) centers(1)];
plot(x,y,'r-', 'LineWidth', 2.5)
set(gca,'ydir','reverse');
hold on
norm_heights_m7 = max_model_7 - Model_H7;
norm_heights_m7 = norm_heights_m7(norm_heights_m7 < 5);
[Values_m, edges] = histcounts(norm_heights_m7, linspace(0, 5, 30), 'Normalization', 'cdf');
centers = (edges(1:end-1)+ edges(2:end))/2;
plot(Values_m, centers, 'b-', 'LineWidth', 2.5);
title('Bearing Area Curves - Stage 7');
pub_fig;
xticks([0 0.2 0.4 0.6 0.8 1]);
x = linspace(0, Values_m(1), 10);
y = [zeros(1,length(x)-1) centers(1)];
plot(x,y,'b-', 'LineWidth', 2.5)
ylabel('Heights (10^{-6} m)')
xlabel('Quantile')
pub_fig;
set(gca,'ydir','reverse');
legend('Data', '', 'Model surface')

subplot(2,3,2)
histogram(max_model_7 - Model_H7, 'Normalization', 'probability', 'BinWidth', 0.5);
title('Histogram of Model 7 Sim');
pub_fig;
ylim([0 1]);
xlim([0 5]);
xlabel('Heights (10^{-6} m)')
ylabel('Frequency')

for i = 2:1:length(Values_s)
    Values_s(i) = Values_s(i) - Values_s(i-1);
    Values_m(i) = Values_m(i) - Values_m(i-1);
end
KL_div_test_stage7 = KLDiv(Values_s, Values_m);

subplot(2,3,3)
pos = get(subplot(2,3,3), 'Position');
Statistics = {'Mean';'Std dev'; 'Sa value'; 'Max'; 'Min';'Mean P2P'; '(KS & KL)'};
Stage_7 = [mean_stage7; stdev_stage7; Sa_stage7; max_stage7; min_stage7; P2P_stage7; p_value_7];
Model_7 = [mean_model_7; stdev_model_7; Sa_model_7; max_model_7; min_model_7; P2P_model_7; KL_div_test_stage7];
T = table(Stage_7, Model_7, 'RowNames', Statistics);
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,...
    'RowName',T.Properties.RowNames,'Units', 'Normalized','ColumnWidth', {70}, 'Position', pos);
pub_fig;
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('Surface Characteristics');
xlabel('All units in µm');

subplot(2,3,6)
pos = get(subplot(2,3,6), 'Position');
Model_Parameters = {'Force (N)'; 'alpha';'beta'; 'kappa 1'; 'kappa 2'; 'sigma0'; 'λ2'; 'MRR (um/min)'};
Values_1 = [Force(7); alpha(7); beta(7); kappa_1(7); kappa_2(7); sigma(7); Fiedler(7); MRR(7)];
T = table(Values_1, 'RowNames', Model_Parameters);
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,...
    'RowName',T.Properties.RowNames,'Units', 'Normalized','ColumnWidth', {70}, 'Position', pos);
pub_fig;
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('Model Parameters - Stage 7');


%% Stage 7 to Stage 8 Parameter Fitting 
%%%%%%% Asperity network %%%%%%%%%

MRR(8) = 0;

Time = 900; 
k = 2; dt = 0.2;
a_0 = 2.5;
Horizon_78 = 0:dt:Time;
Horizon = Horizon_78;
Horizon_78 = Horizon_78 + 4200;

T = zeros(N,length(Horizon));  T(:,1) = Temp;
s = zeros(N,length(Horizon));  s(:,1) = sigma(8);
h = zeros(N,length(Horizon));  h(:,1) = heights_m7;
r = zeros(N,length(Horizon));  r(:,1) = radii_m7;
F = zeros(N,length(Horizon));                   % Force carried by asperities
Rho = zeros(N,1); 
diff_heights = zeros(N,1);

st_d_crit = zeros(length(Horizon),1);
st_active_nodes = zeros(length(Horizon),1);
vector = 0: params{2}: max(heights_m0);
vectorlength = length(vector);

for t = dt:dt:Time
    
    d_crit = solve_for_d(Force(8),r(:,k-1),h(:,k-1),params,vectorlength);
    active_nodes = find(h(:,k-1) > d_crit); % load bearing asperities
    F(active_nodes, k) = (2/3)*E* sqrt(abs(r(active_nodes,k-1))).*((abs(h(active_nodes,k-1) - d_crit.*ones(length(active_nodes),1))).^(3/2));   % force carried by each node. 
    
    % Temp dynamics
    T(:,k) = T(:,k-1) + dt*(-alpha(8)*(T(:,k-1)- Temp) + kappa_1(8)*F(:, k-1) + kappa_2(8)*speed(8)); 
    
    % Sigma dynamics
    s(:,k) = s(:,k-1) + (beta(8).*(T(:,k) - T(:,k-1))./Temp); 
    
    % Differential height dynamics
    diff_heights(:,1) = (dt.*s(:,k-1).*T(:,k-1).*F(:,k-1));
    h(:,k) = h(:,k-1) - diff_heights;    
    
    % Radii update and MRR computation
    Rho(:,1) = h(:,k-1)./(4*r(:,k-1));
    MR_coeff_eta = exp(-(h(:,k-1)/max(heights_m0).^2)).*exp(-(a_0./r(:,k-1)).^0.0833);
    r(:,k) = r(:,k-1) + (MR_coeff_eta.*Rho(:,1).*(diff_heights));
    MRR(8) = MRR(8) + sum(pi.*(diff_heights.^2).*h(:,k-1).*MR_coeff_eta.*(1 - MR_coeff_eta)/4);
    
    % store
    st_d_crit(k-1) = d_crit;
    st_active_nodes(k-1) = length(active_nodes);
    k = k+1;
end

heights_m8 = h(:,length(Horizon));
radii_m8 = r(:,length(Horizon));
MRR(8) = MRR(8) / (250*250*Time/60);

% height, sigma and Temp change plots
figure
subplot(1,3,1)
h_78 = h(active_nodes,:);
plot(Horizon,h(active_nodes,:),'LineWidth',2)
pub_fig;
hold on
title('Height change dynamics of active asperities')
xlabel('Time (sec)')
ylabel('Height (um)')

subplot(1,3,2)
plot(Horizon,r(active_nodes,:),'LineWidth',2)
pub_fig;
hold on
title('Radii change of asperities')
xlabel('Time (sec)')
ylabel('Radii (um)')

subplot(1,3,3)
Temp_78 = T(active_nodes,:);
plot(Horizon,T(active_nodes,:),'LineWidth',2)
pub_fig;
hold on
title('Temp variation: Stage 7 - 8')
xlabel('Time (sec)')
ylabel('Temp (degree C)')


%% Stage 8 - Data and Model Comparison - Plotting

Table8 = readtable('matrix_stage8.csv');
data_stage8 = table2array(Table8);
data_p2p_stage8 = reshape(data_stage8, [102 32]);
data_stage8 = reshape(data_stage8, [3264 1]);
data8 = mean(data_stage8) - data_stage8;

% Stage 8 Data Statistics
max_stage8 = max(data8);
mean_stage8 = mean(data8);
stdev_stage8 = std(data8);
min_stage8 = min(data8);
Hurst_stage8 = Gen_hurst(data8);
Sa_stage8 = 0;
for i = 1: 3264
    Sa_stage8 = Sa_stage8 + ((1/3264)*(abs(data_stage8(i) - mean_stage8)));
end
P2P_stage8 = (sum(maxk(data_stage8,320,1),'all') - sum(mink(data_stage8,320,1), 'all'))/320;

% Model 8 Statistics
max_model_8 = max(heights_m8);
min_model_8 = min(heights_m8);
[Model_H8, mean_model_8, stdev_model_8, Sa_model_8, Hurst_model_8, KS_test_8, p_value_8] = surface_roughness(radii_m8, initial_asperity_radii, heights_m8, N, data_stage8, packing_density);
P2P_model_8 = (sum(maxk(heights_m8,20,1),'all') - sum(mink(heights_m8,20,1), 'all'))/20;

vector = 0: params{2}: max(heights_m0);
vectorlength = length(vector);
[AdjMatrix_8, Fiedler(8), E8, bins_CC8] = graph_evolution(x_pos, y_pos, packing_density, Force(8), radii_m8, heights_m8, heights_m0, params, vectorlength, 8);

figure
subplot(2,3,1)
histogram(data_stage8, 'Normalization', 'probability', 'BinWidth', 0.5);
title('Histogram of Stage 8 Data');
pub_fig;
ylim([0 0.4]);
xlim([0 20]);
xlabel('Heights (10^{-6} m)')
ylabel('Frequency')

subplot(2,3,4)
[Values_s, edges] = histcounts(data_stage8, linspace(0, 20, 30), 'Normalization', 'cdf');
centers = (edges(1:end-1)+edges(2:end))/2;
plot(Values_s, centers, 'r-', 'LineWidth', 2.5)
pub_fig;
xticks([0 0.2 0.4 0.6 0.8 1]);
hold on
x = linspace(0, Values_s(1), 10);
y = [zeros(1,length(x)-1) centers(1)];
plot(x,y,'r-', 'LineWidth', 2.5)
set(gca,'ydir','reverse');
hold on
norm_heights_m8 = max_model_8 - Model_H8;
norm_heights_m8 = norm_heights_m8(norm_heights_m8 < 21);
[Values_m, edges] = histcounts(norm_heights_m8, linspace(0, 20, 30), 'Normalization', 'cdf');
centers = (edges(1:end-1)+ edges(2:end))/2;
plot(Values_m, centers, 'b-', 'LineWidth', 2.5);
title('Bearing Area Curves - Stage 8');
pub_fig;
xticks([0 0.2 0.4 0.6 0.8 1]);
x = linspace(0, Values_m(1), 10);
y = [zeros(1,length(x)-1) centers(1)];
plot(x,y,'b-', 'LineWidth', 2.5)
ylabel('Heights (10^{-6} m)')
xlabel('Quantile')
pub_fig;
set(gca,'ydir','reverse');
legend('Data', '', 'Model surface')

subplot(2,3,2)
histogram(max_model_8 - Model_H8, 'Normalization', 'probability', 'BinWidth', 0.5);
title('Histogram of Model 8 Sim');
pub_fig;
ylim([0 0.4]);
xlim([0 20]);
xlabel('Heights (10^{-6} m)')
ylabel('Frequency')

for i = 2:1:length(Values_s)
    Values_s(i) = Values_s(i) - Values_s(i-1);
    Values_m(i) = Values_m(i) - Values_m(i-1);
end
KL_div_test_stage8 = KLDiv(Values_s, Values_m);

subplot(2,3,3)
pos = get(subplot(2,3,3), 'Position');
Statistics = {'Mean';'Std dev'; 'Sa value'; 'Max'; 'Min';'Mean P2P'; '(KS & KL)'};
Stage_8 = [mean_stage8; stdev_stage8; Sa_stage8; max_stage8; min_stage8; P2P_stage8; p_value_8];
Model_8 = [mean_model_8; stdev_model_8; Sa_model_8; max_model_8; min_model_8; P2P_model_8; KL_div_test_stage8];
T = table(Stage_8, Model_8, 'RowNames', Statistics);
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,...
    'RowName',T.Properties.RowNames,'Units', 'Normalized','ColumnWidth', {70}, 'Position', pos);
pub_fig;
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('Surface Characteristics');
xlabel('All units in µm');

subplot(2,3,6)
pos = get(subplot(2,3,6), 'Position');
Model_Parameters = {'Force (N)'; 'alpha';'beta'; 'kappa 1'; 'kappa 2'; 'sigma0'; 'λ2'; 'MRR (um/min)'};
Values_1 = [Force(8); alpha(8); beta(8); kappa_1(8); kappa_2(8); sigma(8); Fiedler(8); MRR(8)];
T = table(Values_1, 'RowNames', Model_Parameters);
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,...
    'RowName',T.Properties.RowNames,'Units', 'Normalized','ColumnWidth', {70}, 'Position', pos);
pub_fig;
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('Model Parameters - Stage 8');


%% Stage 8 to Stage 9 Parameter Fitting 
%%%%%%% Asperity network %%%%%%%%%
MRR(9) = 0;

Time = 900; 
k = 2; dt = 0.2;
a_0 = 2.5;
Horizon_89 = 0:dt:Time;
Horizon = Horizon_89;
Horizon_89 = Horizon_89 + 5100;

T = zeros(N,length(Horizon));  T(:,1) = Temp;
s = zeros(N,length(Horizon));  s(:,1) = sigma(9);
h = zeros(N,length(Horizon));  h(:,1) = heights_m8;
r = zeros(N,length(Horizon));  r(:,1) = radii_m8;
F = zeros(N,length(Horizon));                   % Force carried by asperities
Rho = zeros(N,1); 
diff_heights = zeros(N,1);

st_d_crit = zeros(length(Horizon),1);
st_active_nodes = zeros(length(Horizon),1);
vector = 0: params{2}: max(heights_m0);
vectorlength = length(vector);

for t = dt:dt:Time
    
    d_crit = solve_for_d(Force(9),r(:,k-1),h(:,k-1),params,vectorlength);
    active_nodes = find(h(:,k-1) > d_crit); % load bearing asperities
    F(active_nodes, k) = (2/3)*E* sqrt(abs(r(active_nodes,k-1))).*((abs(h(active_nodes,k-1) - d_crit.*ones(length(active_nodes),1))).^(3/2));   % force carried by each node. 
    
    % Temp dynamics
    T(:,k) = T(:,k-1) + dt*(-alpha(9)*(T(:,k-1)- Temp) + kappa_1(9)*F(:, k-1) + kappa_2(9)*speed(9)); 
    
    % Sigma dynamics
    s(:,k) = s(:,k-1) + (beta(9).*(T(:,k) - T(:,k-1))./Temp); 
    
    % Differential height dynamics
    diff_heights(:,1) = (dt.*s(:,k-1).*T(:,k-1).*F(:,k-1));
    h(:,k) = h(:,k-1) - diff_heights;    
    
    % Radii update and MRR computation
    Rho(:,1) = h(:,k-1)./(4*r(:,k-1));
    MR_coeff_eta = exp(-(h(:,k-1)/max(heights_m0).^2)).*exp(-(a_0./r(:,k-1)).^0.0833);
    r(:,k) = r(:,k-1) + (MR_coeff_eta.*Rho(:,1).*(diff_heights));
    MRR(9) = MRR(9) + sum(pi.*(diff_heights.^2).*h(:,k-1).*MR_coeff_eta.*(1 - MR_coeff_eta)/4);
    
    % store
    st_d_crit(k-1) = d_crit;
    st_active_nodes(k-1) = length(active_nodes);
    k = k+1;
end

heights_m9 = h(:,length(Horizon));
radii_m9 = r(:,length(Horizon));
MRR(9) = MRR(9) / (250*250*Time/60);

% height, sigma and Temp change plots
figure
subplot(1,3,1)
h_89 = h(active_nodes,:);
plot(Horizon,h(active_nodes,:),'LineWidth',2)
pub_fig;
hold on
title('Height change dynamics of active asperities')
xlabel('Time (sec)')
ylabel('Height (um)')

subplot(1,3,2)
plot(Horizon,r(active_nodes,:),'LineWidth',2)
pub_fig;
hold on
title('Radii change of asperities')
xlabel('Time (sec)')
ylabel('Radii (um)')

subplot(1,3,3)
Temp_89 = T(active_nodes,:);
plot(Horizon,T(active_nodes,:),'LineWidth',2)
pub_fig;
hold on
title('Temp variation: Stage 8 - 9')
xlabel('Time (sec)')
ylabel('Temp (degree C)')



%% Stage 9 - Data and Model Comparison - Plotting

Table9 = readtable('matrix_stage9.csv');
data_stage9 = table2array(Table9);
data_p2p_stage9 = reshape(data_stage9, [102 32]);
data_stage9 = reshape(data_stage9, [3264 1]);
data9 = mean(data_stage9) - data_stage9;

% Stage 9 Data Statistics
max_stage9 = max(data9);
mean_stage9 = mean(data9);
stdev_stage9 = std(data9);
min_stage9 = min(data9);
Hurst_stage9 = Gen_hurst(data9);
Sa_stage2 = 0;
for i = 1: 3264
    Sa_stage9 = Sa_stage9 + ((1/3264)*(abs(data_stage9(i) - mean_stage9)));
end
P2P_stage9 = (sum(maxk(data_stage9,320,1),'all') - sum(mink(data_stage9,320,1), 'all'))/320;

% Model 9 Statistics
max_model_9 = max(heights_m9);
min_model_9 = min(heights_m9);
[Model_H9, mean_model_9, stdev_model_9, Sa_model_9, Hurst_model_9, KS_test_9, p_value_9] = surface_roughness(radii_m9, initial_asperity_radii, heights_m9, N, data_stage9, packing_density);
P2P_model_9 = (sum(maxk(heights_m9,20,1),'all') - sum(mink(heights_m9,20,1), 'all'))/20;

vector = 0: params{2}: max(heights_m0);
vectorlength = length(vector);
[AdjMatrix_9, Fiedler(9), E9, bins_CC9] = graph_evolution(x_pos, y_pos, packing_density, Force(9), radii_m9, heights_m9, heights_m0, params, vectorlength, 2);

figure
subplot(2,3,1)
histogram(data_stage9, 'Normalization', 'probability', 'BinWidth', 0.5);
title('Histogram of Stage 9 Data');
pub_fig;
ylim([0 0.4]);
xlim([0 20]);
xlabel('Heights (10^{-6} m)')
ylabel('Frequency')

subplot(2,3,4)
[Values_s, edges] = histcounts(data_stage9, linspace(0, 20, 30), 'Normalization', 'cdf');
centers = (edges(1:end-1)+edges(2:end))/2;
plot(Values_s, centers, 'r-', 'LineWidth', 2.5)
pub_fig;
xticks([0 0.2 0.4 0.6 0.8 1]);
hold on
x = linspace(0, Values_s(1), 10);
y = [zeros(1,length(x)-1) centers(1)];
plot(x,y,'r-', 'LineWidth', 2.5)
set(gca,'ydir','reverse');
hold on
norm_heights_m9 = max_model_9 - Model_H9;
norm_heights_m9 = norm_heights_m9(norm_heights_m9 < 21);
[Values_m, edges] = histcounts(norm_heights_m9, linspace(0, 20, 30), 'Normalization', 'cdf');
centers = (edges(1:end-1)+ edges(2:end))/2;
plot(Values_m, centers, 'b-', 'LineWidth', 2.5);
title('Bearing Area Curves - Stage 9');
pub_fig;
xticks([0 0.2 0.4 0.6 0.8 1]);
x = linspace(0, Values_m(1), 10);
y = [zeros(1,length(x)-1) centers(1)];
plot(x,y,'b-', 'LineWidth', 2.5)
ylabel('Heights (10^{-6} m)')
xlabel('Quantile')
pub_fig;
set(gca,'ydir','reverse');
legend('Data', '', 'Model surface')

subplot(2,3,2)
histogram(max_model_9 - Model_H9, 'Normalization', 'probability', 'BinWidth', 0.5);
title('Histogram of Model 9 Sim');
pub_fig;
ylim([0 0.4]);
xlim([0 20]);
xlabel('Heights (10^{-6} m)')
ylabel('Frequency')

for i = 2:1:length(Values_s)
    Values_s(i) = Values_s(i) - Values_s(i-1);
    Values_m(i) = Values_m(i) - Values_m(i-1);
end
KL_div_test_stage9 = KLDiv(Values_s, Values_m);

subplot(2,3,3)
pos = get(subplot(2,3,3), 'Position');
Statistics = {'Mean';'Std dev'; 'Sa value'; 'Max'; 'Min';'Mean P2P'; '(KS & KL)'};
Stage_9 = [mean_stage9; stdev_stage9; Sa_stage9; max_stage9; min_stage9; P2P_stage9; p_value_9];
Model_9 = [mean_model_9; stdev_model_9; Sa_model_9; max_model_9; min_model_9; P2P_model_9; KL_div_test_stage9];
T = table(Stage_9, Model_9, 'RowNames', Statistics);
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,...
    'RowName',T.Properties.RowNames,'Units', 'Normalized','ColumnWidth', {70}, 'Position', pos);
pub_fig;
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('Surface Characteristics');
xlabel('All units in µm');

subplot(2,3,6)
pos = get(subplot(2,3,6), 'Position');
Model_Parameters = {'Force (N)'; 'alpha';'beta'; 'kappa 1'; 'kappa 2'; 'sigma0'; 'λ2'; 'MRR (um/min)'};
Values_1 = [Force(9); alpha(9); beta(9); kappa_1(9); kappa_2(9); sigma(9); Fiedler(9); MRR(9)];
T = table(Values_1, 'RowNames', Model_Parameters);
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,...
    'RowName',T.Properties.RowNames,'Units', 'Normalized','ColumnWidth', {70}, 'Position', pos);
pub_fig;
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('Model Parameters - Stage 9');


%% Stage 9 to Stage 10 Parameter Fitting 
%%%%%%% Asperity network %%%%%%%%%

MRR(10) = 0;

Time = 900; 
k = 2; dt = 0.2;
a_0 = 2.5;
Horizon_910 = 0:dt:Time;
Horizon = Horizon_910;
Horizon_910 = Horizon_910 + 6000;

T = zeros(N,length(Horizon));  T(:,1) = Temp;
s = zeros(N,length(Horizon));  s(:,1) = sigma(10);
h = zeros(N,length(Horizon));  h(:,1) = heights_m9;
r = zeros(N,length(Horizon));  r(:,1) = radii_m9;
F = zeros(N,length(Horizon));                   % Force carried by asperities
Rho = zeros(N,1); 
diff_heights = zeros(N,1);

st_d_crit = zeros(length(Horizon),1);
st_active_nodes = zeros(length(Horizon),1);
vector = 0: params{2}: max(heights_m0);
vectorlength = length(vector);

for t = dt:dt:Time
    
    d_crit = solve_for_d(Force(10),r(:,k-1),h(:,k-1),params,vectorlength);
    active_nodes = find(h(:,k-1) > d_crit); % load bearing asperities
    F(active_nodes, k) = (2/3)*E* sqrt(abs(r(active_nodes,k-1))).*((abs(h(active_nodes,k-1) - d_crit.*ones(length(active_nodes),1))).^(3/2));   % force carried by each node. 
    
    % Temp dynamics
    T(:,k) = T(:,k-1) + dt*(-alpha(10)*(T(:,k-1)- Temp) + kappa_1(10)*F(:, k-1) + kappa_2(10)*speed(10)); 
    
    % Sigma dynamics
    s(:,k) = s(:,k-1) + (beta(10).*(T(:,k) - T(:,k-1))./Temp); 
    
    % Differential height dynamics
    diff_heights(:,1) = (dt.*s(:,k-1).*T(:,k-1).*F(:,k-1));
    h(:,k) = h(:,k-1) - diff_heights;    
    
    % Radii update and MRR computation
    Rho(:,1) = h(:,k-1)./(4*r(:,k-1));
    MR_coeff_eta = exp(-(h(:,k-1)/max(heights_m0).^2)).*exp(-(a_0./r(:,k-1)).^0.0833);
    r(:,k) = r(:,k-1) + (MR_coeff_eta.*Rho(:,1).*(diff_heights));
    MRR(10) = MRR(10) + sum(pi.*(diff_heights.^2).*h(:,k-1).*MR_coeff_eta.*(1 - MR_coeff_eta)/4);
    
    % store
    st_d_crit(k-1) = d_crit;
    st_active_nodes(k-1) = length(active_nodes);
    k = k+1;
end

heights_m10 = h(:,length(Horizon));
radii_m10 = r(:,length(Horizon));
MRR(10) = MRR(10) / (250*250*Time/60);

% height, sigma and Temp change plots
figure
subplot(1,3,1)
h_910 = h(active_nodes,:);
plot(Horizon,h(active_nodes,:),'LineWidth',2)
pub_fig;
hold on
title('Height change dynamics of active asperities')
xlabel('Time (sec)')
ylabel('Height (um)')

subplot(1,3,2)
plot(Horizon,r(active_nodes,:),'LineWidth',2)
pub_fig;
hold on
title('Radii change of asperities')
xlabel('Time (sec)')
ylabel('Radii (um)')

subplot(1,3,3)
Temp_910 = T(active_nodes,:);
plot(Horizon,T(active_nodes,:),'LineWidth',2)
pub_fig;
hold on
title('Temp variation: Stage 9 - 10')
xlabel('Time (sec)')
ylabel('Temp (degree C)')



%% Stage 10 - Data and Model Comparison - Plotting

Table10 = readtable('matrix_stage10.csv');
data_stage10 = table2array(Table10);
data_p2p_stage10 = reshape(data_stage10, [102 32]);
data_stage10 = reshape(data_stage10, [3264 1]);
data10 = mean(data_stage10) - data_stage10;

% Stage 10 Data Statistics
max_stage10 = max(data10);
mean_stage10 = mean(data10);
stdev_stage10 = std(data10);
min_stage10 = min(data10);
Hurst_stage10 = Gen_hurst(data10);
Sa_stage10 = 0;
for i = 1: 3264
    Sa_stage10 = Sa_stage10 + ((1/3264)*(abs(data_stage10(i) - mean_stage10)));
end
P2P_stage10 = (sum(maxk(data_stage10,320,1),'all') - sum(mink(data_stage10,320,1), 'all'))/320;

% Model 10 Statistics
max_model_10 = max(heights_m10);
min_model_10 = min(heights_m10);
[Model_H10, mean_model_10, stdev_model_10, Sa_model_10, Hurst_model_10, KS_test_10, p_value_10] = surface_roughness(radii_m10, initial_asperity_radii, heights_m10, N, data_stage10, packing_density);
P2P_model_10 = (sum(maxk(heights_m10,20,1),'all') - sum(mink(heights_m10,20,1), 'all'))/20;

vector = 0: params{2}: max(heights_m0);
vectorlength = length(vector);
[AdjMatrix_10, Fiedler(10), E10, bins_CC10] = graph_evolution(x_pos, y_pos, packing_density, Force(10), radii_m10, heights_m10, heights_m0, params, vectorlength, 10);

figure
subplot(2,3,1)
histogram(data_stage10, 'Normalization', 'probability', 'BinWidth', 0.5);
title('Histogram of Stage 10 Data');
pub_fig;
ylim([0 0.4]);
xlim([0 20]);
xlabel('Heights (10^{-6} m)')
ylabel('Frequency')

subplot(2,3,4)
[Values_s, edges] = histcounts(data_stage10, linspace(0, 20, 30), 'Normalization', 'cdf');
centers = (edges(1:end-1)+edges(2:end))/2;
plot(Values_s, centers, 'r-', 'LineWidth', 2.5)
pub_fig;
xticks([0 0.2 0.4 0.6 0.8 1]);
hold on
x = linspace(0, Values_s(1), 10);
y = [zeros(1,length(x)-1) centers(1)];
plot(x,y,'r-', 'LineWidth', 2.5)
set(gca,'ydir','reverse');
hold on
norm_heights_m10 = max_model_10 - Model_H10;
norm_heights_m10 = norm_heights_m10(norm_heights_m10 < 21);
[Values_m, edges] = histcounts(norm_heights_m10, linspace(0, 20, 30), 'Normalization', 'cdf');
centers = (edges(1:end-1)+ edges(2:end))/2;
plot(Values_m, centers, 'b-', 'LineWidth', 2.5);
title('Bearing Area Curves - Stage 10');
pub_fig;
xticks([0 0.2 0.4 0.6 0.8 1]);
x = linspace(0, Values_m(1), 10);
y = [zeros(1,length(x)-1) centers(1)];
plot(x,y,'b-', 'LineWidth', 2.5)
ylabel('Heights (10^{-6} m)')
xlabel('Quantile')
pub_fig;
set(gca,'ydir','reverse');
legend('Data', '', 'Model surface')

subplot(2,3,2)
histogram(max_model_10 - Model_H10, 'Normalization', 'probability', 'BinWidth', 0.5);
title('Histogram of Model 10 Sim');
pub_fig;
ylim([0 0.4]);
xlim([0 20]);
xlabel('Heights (10^{-6} m)')
ylabel('Frequency')

for i = 2:1:length(Values_s)
    Values_s(i) = Values_s(i) - Values_s(i-1);
    Values_m(i) = Values_m(i) - Values_m(i-1);
end
KL_div_test_stage10 = KLDiv(Values_s, Values_m);

subplot(2,3,3)
pos = get(subplot(2,3,3), 'Position');
Statistics = {'Mean';'Std dev'; 'Sa value'; 'Max'; 'Min';'Mean P2P'; '(KS & KL)'};
Stage_10 = [mean_stage10; stdev_stage10; Sa_stage10; max_stage10; min_stage10; P2P_stage10; p_value_10];
Model_10 = [mean_model_10; stdev_model_10; Sa_model_10; max_model_10; min_model_10; P2P_model_10; KL_div_test_stage10];
T = table(Stage_10, Model_10, 'RowNames', Statistics);
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,...
    'RowName',T.Properties.RowNames,'Units', 'Normalized','ColumnWidth', {70}, 'Position', pos);
pub_fig;
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('Surface Characteristics');
xlabel('All units in µm');

subplot(2,3,6)
pos = get(subplot(2,3,6), 'Position');
Model_Parameters = {'Force (N)'; 'alpha';'beta'; 'kappa 1'; 'kappa 2'; 'sigma0'; 'λ2'; 'MRR (um/min)'};
Values_1 = [Force(10); alpha(10); beta(10); kappa_1(10); kappa_2(10); sigma(10); Fiedler(10); MRR(10)];
T = table(Values_1, 'RowNames', Model_Parameters);
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,...
    'RowName',T.Properties.RowNames,'Units', 'Normalized','ColumnWidth', {70}, 'Position', pos);
pub_fig;
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('Model Parameters - Stage 10');
