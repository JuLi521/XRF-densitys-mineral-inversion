%% XRF-density joint mineral inversion: ICDP-DIVE 5071_1_B
% Reproducible 10-sample demonstration (150-250 m).
% This main script retains the mineral model, likelihood weighting,
% MCMC settings, and visualization style of the original research script.
%
% MCMC method background:
% Cheng, L., Jin, G., Michelena, R., & Tura, A. (2022). 
% Practical Bayesian Inversions for Rock Composition and Petrophysical Endpoints in Multimineral Analysis. 
% SPE Reservoir Evaluation & Engineering, 25(04), 849–865. https://doi.org/10.2118/210576-PA

% IMPORTANT: The bundled mcmces_xrf1.m is copied verbatim from the
% author-supplied code. Its exact software provenance has NOT been verified.
% Its acceptance rule differs from the standard Goodman-Weare stretch move.
%
% Requires: MATLAB + Statistics and Machine Learning Toolbox (ksdensity).
% Run this script in the directory containing the CSV and mcmces_xrf1.m.

clear; clc; close all;
rng(2026,'twister');

%% 1. Read previously prepared, depth-aligned example data
% No processing of the full original XRF, lithology, density or log files.
data = readtable('DIVE_5071_1_B_KINZ_10samples_150-250m.csv');
dep = data.Depth_m;
den_kin = data.Density_g_cm3;
oxide_names = {'SiO2','K2O','Al2O3','CaO','MgO','MnO','Na2O','FeO','Fe2O3'};
oxide_columns = strcat(oxide_names,'_wt_pct');
XRF_Kin = data{:,oxide_columns};
[num_depths,nn] = size(XRF_Kin);
assert(num_depths == 10 && all(dep >= 150 & dep <= 250), ...
    'The example dataset must contain 10 samples from 150-250 m.');
assert(all(isfinite(XRF_Kin(:))) && all(isfinite(den_kin)), ...
    'All XRF and density entries must be finite.');
[dep,sortidx] = sort(dep);
XRF_Kin = XRF_Kin(sortidx,:);
den_kin = den_kin(sortidx);

%% 2. Mineral reference compositions (unchanged from original script)
%         SiO2   K2O   Al2O3 CaO   MgO   MnO  Na2O  FeO  Fe2O3
G_kin = [43.20  11.29  12.22  0     28.97 0     0     0     0; ... % Phlogopite
         35.21   9.20   9.96  0      0    0     0    42.10  0; ... % Annite
        100.00   0      0     0      0    0     0     0     0; ... % Quartz
         68.74   0     19.44  0      0    0    11.82  0     0; ... % Albite
         43.20   0     36.65 20.15   0    0     0     0     0; ... % Anorthite
         38.03   0     20.98  4.08   3.85 3.06  0.02 29.93 0.53; ... % Garnet
         64.76  16.92  18.32  0      0    0     0     0     0];    % K-feldspar

den_min = [2.784 3.3 2.65 2.62 2.76 4.10 2.56]; % g/cm3
mineral_names = {'Phl','Ant','Qtz','Ab','An','Grt','Kfs'};
N = 7;

%% 3. Original MCMC settings and likelihood weighting
% NOTE: These likelihood weights follow the original example exactly.
% They are not independently calibrated analytical uncertainties.
sigma = std(XRF_Kin,0,1);
if any(~isfinite(sigma(1:8))) || any(sigma(1:8) <= 0)
    error('The first 8 XRF columns need nonzero sample standard deviations.');
end
Nwalkers = 100;
Nstep = 2000;
BurnIn = 0.6;
StepSize = 5;
ThinChain = 1;
crop = ceil(Nstep*BurnIn);

Weight_mean = zeros(num_depths,N);
Weight_max = zeros(num_depths,N); % Marginal posterior KDE mode (not joint MAP)
devi = zeros(num_depths,N);
W_samples_all = cell(num_depths,1);
acceptance_rate = zeros(num_depths,1);

%% 4. Invert each sample using the UNMODIFIED mcmces_xrf1.m
for k = 1:num_depths
    % Original program uses the first 8 oxides in the likelihood.
    % The ninth column (Fe2O3) is retained for forward reconstruction.
    Likelihood = @(m) ...
       -0.5*sum(((XRF_Kin(k,1:8) - m'*G_kin(:,1:8)).^2)./sigma(1:8)) ...
       -0.5*(den_kin(k) - 1/sum(m'./den_min)).^2/0.00010;
    Prior = @(m) all(m >= 0 & m <= 1);

    [Mchain,~,rate,~] = mcmces_xrf1( ...
        Prior,Likelihood,N,Nwalkers,Nstep,StepSize,ThinChain);
    acceptance_rate(k) = rate;
    ModelsBurn = Mchain(:,:,crop+1:end);
    samples = reshape(ModelsBurn,N,[]);
    W_samples_all{k} = samples;

    for j = 1:N
        Weight_mean(k,j) = mean(samples(j,:));
        [f,xi] = ksdensity(samples(j,:));
        [~,imax] = max(f);
        Weight_max(k,j) = xi(imax);
        devi(k,j) = std(samples(j,:));
    end
    fprintf('Sample %2d/10  depth = %7.2f m  acceptance = %.3f\n', ...
        k, dep(k), rate);
end

%% 5. Figure 1: posterior PDFs with depth (original surface-map style)
n_bins = 100;
x_bins = linspace(0,1,n_bins);
xlm = [zeros(N,1) ones(N,1)];
clm = [zeros(N,1) ones(N,1)*0.3];

cmap_turbo = jet(256);       % Original script uses jet, not turbo
cmap_turbo(1:64,:) = [];
cmap_turbo(end-31:end,:) = [];
fade_len = 64;
white = [1 1 1];
blue = cmap_turbo(1,:);
fade = zeros(fade_len,3);
for i = 1:3
    fade(:,i) = linspace(white(i),blue(i),fade_len)';
end
cmap_trimmed = [fade; cmap_turbo];

figure('Color','w','Name','Figure 1: Mineral posterior PDFs', ...
       'Position',[80 80 1550 620]);
for j = 1:N
    density_map = zeros(num_depths,n_bins);
    for d = 1:num_depths
        f = ksdensity(W_samples_all{d}(j,:),x_bins);
        density_map(d,:) = f/sum(f); % As in the original script
    end

    subplot(1,N,j);
    [X,Y] = meshgrid(x_bins,dep);
    surface(X,Y,zeros(size(density_map)),density_map, ...
        'EdgeColor','none','FaceColor','interp');
    view(2); colormap(cmap_trimmed);
    set(gca,'YDir','reverse'); hold on;
    caxis(clm(j,:));
    plot(Weight_mean(:,j),dep,'k*','LineWidth',0.2,'MarkerSize',5);
    grid minor;
    xlabel('Weight Fraction');
    if j == 1, ylabel('Depth (m)'); end
    xlim(xlm(j,:)); ylim([150 250]);
    title(mineral_names{j});
end
cb = colorbar;
ylabel(cb,'Normalized KDE intensity');

%% 6. Figure 2: measured versus reconstructed oxide contents
% Preserve original: black stars (inversion) and red circles (data).
XRF_Inv = Weight_mean*G_kin;
figure('Color','w','Name','Figure 2: XRF reconstruction', ...
       'Position',[60 120 1660 550]);
for j = 1:nn
    subplot(1,nn,j); hold on;
    h1 = plot(XRF_Inv(:,j),dep,'k*','MarkerSize',5);
    h2 = plot(XRF_Kin(:,j),dep,'ro','MarkerSize',5);
    ylim([150 250]); set(gca,'YDir','reverse');
    title(oxide_names{j}); xlabel('wt.%'); grid minor; box on;
    if j == 1
        ylabel('Depth (m)');
        legend([h1 h2],{'Inversion','Data'},'Location','best');
    end
end

%% 7. Figure 3: measured versus reconstructed bulk density
Den_Inv = 1./sum(Weight_mean./den_min,2);
figure('Color','w','Name','Figure 3: Density reconstruction', ...
       'Position',[180 120 480 620]);
hold on;
h1 = plot(Den_Inv,dep,'k*','MarkerSize',6);
h2 = plot(den_kin,dep,'ro','MarkerSize',6);
set(gca,'YDir','reverse'); ylim([150 250]);
xlabel('Bulk Density (g/cm^3)'); ylabel('Depth (m)');
title('Bulk density reconstruction'); grid minor; box on;
legend([h1 h2],{'Inversion','Data'},'Location','best');

%% 8. Convert weight fractions to volume fractions
% For mass fractions w_i, V_i = (w_i/rho_i)/sum_j(w_j/rho_j).
Volume_frac = (Weight_mean./den_min) ./ sum(Weight_mean./den_min,2);

% Preserve original petrological grouping:
% Biotite (Phl+Ant), Quartz, Plagioclase (Ab+An), Garnet, K-feldspar
Volume_group = [Volume_frac(:,1)+Volume_frac(:,2), ...
                Volume_frac(:,3), ...
                Volume_frac(:,4)+Volume_frac(:,5), ...
                Volume_frac(:,6), ...
                Volume_frac(:,7)];
volume_group_names = {'Bt','Qtz','Pl','Grt','Kfs'};

%% 9. Figure 4: grouped mineral volume fractions (original depth-track style)
figure('Color','w','Name','Figure 4: Mineral volume fractions', ...
       'Position',[80 80 1320 620]);
for j = 1:numel(volume_group_names)
    subplot(1,numel(volume_group_names),j); hold on;
    plot(Volume_group(:,j),dep,'Marker','o', ...
        'MarkerFaceColor','k','MarkerSize',5,'LineStyle','none', ...
        'Color','k');
    set(gca,'YDir','reverse');
    xlabel('Volume Fraction');
    if j == 1, ylabel('Depth (m)'); end
    xlim([0 1]); ylim([150 250]);
    title(volume_group_names{j});
    grid minor; box on;
end

%% 10. Export numerical results (not needed by the inversion itself)
result = table(dep,den_kin,Den_Inv,acceptance_rate, ...
    'VariableNames',{'Depth_m','Density_observed','Density_reconstructed', ...
                     'Acceptance_rate'});
for j = 1:N
    result.(['W_' mineral_names{j}]) = Weight_mean(:,j);
    result.(['V_' mineral_names{j}]) = Volume_frac(:,j);
end
for j = 1:numel(volume_group_names)
    result.(['V_' volume_group_names{j}]) = Volume_group(:,j);
end
writetable(result,'Inversion_results_10samples.csv');
