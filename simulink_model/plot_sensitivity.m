%% plot_sensitivity.m
% Generates all sensitivity analysis plots from saved results.
% Run run_sensitivity.m first to generate the .mat files. Then add .mat
% file names below

%% Configuring results files to plot
base_dir    = fileparts(mfilename('fullpath'));
if isempty(base_dir), base_dir = pwd; end
results_dir = fullfile(base_dir, 'results');

% UPDATE FILENAMES HERE:
dv_mat_file    = fullfile(results_dir, 'dv_sensitivity_20260322_123359.mat');
param_mat_file = fullfile(results_dir, 'param_sensitivity_20260322_123359.mat');

load(dv_mat_file,    'dv_results');
load(param_mat_file, 'param_results');

n_dv     = numel(dv_results);
n_params = numel(param_results);
angle_dv = {'phi','sigma'};  % convert to degrees for axes

%% Figure 1: DV Sensitivity Curves (3×3 grid)
figure('Name','DV Sensitivity Curves','NumberTitle','off','Position',[50 50 1100 750]);

for i = 1:n_dv
    subplot(3, 3, i);

    pv      = dv_results(i).perturbed_values;
    pc      = dv_results(i).percent_changes;
    xstar_v = dv_results(i).x_star_val;

    if ismember(dv_results(i).var_name, angle_dv)
        pv      = rad2deg(pv);
        xstar_v = rad2deg(xstar_v);
        xlabel_str = sprintf('%s  (deg)', dv_results(i).var_name);
    else
        xlabel_str = dv_results(i).var_name;
    end

    plot(pv, pc, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerFaceColor','b');
    hold on;
    xline(xstar_v, 'k--', 'LineWidth', 1.2);
    yline(0,        'r--', 'LineWidth', 1.0);
    hold off;

    xlabel(xlabel_str, 'Interpreter','none');
    ylabel('% \Delta Social Value');
    title(dv_results(i).var_name, 'Interpreter','none');
    grid on;
end

% Leave subplot 8 and 9 blank if < 9 DVs (we have 7, so slots 8-9 are empty)
sgtitle('Design Variable Sensitivity Curves');

% %% Figure 2 - DV Driver Ranking (max % change, sorted descending visually)
% max_abs_pct_dv = zeros(n_dv, 1);
% for i = 1:n_dv
%     max_abs_pct_dv(i) = max(abs(dv_results(i).percent_changes));
% end
% 
% % Sort ascending so the largest value appears at the top of the horizontal chart
% [sorted_dv, sidx_dv] = sort(max_abs_pct_dv, 'ascend');
% sorted_names_dv = {dv_results(sidx_dv).var_name};
% 
% figure('Name','DV Driver Ranking','NumberTitle','off','Position',[200 200 700 400]);
% barh(1:n_dv, sorted_dv, 0.6);
% yticks(1:n_dv);
% yticklabels(sorted_names_dv);
% xlabel('Max % Change in Social Value');
% title('Design Variable Sensitivity — Driver Ranking');
% grid on;
% ax = gca; ax.XGrid = 'on'; ax.YGrid = 'off';

%% Figure 3 - Parameter Sensitivity Curves
n_cols_p = 3;
n_rows_p = ceil(n_params / n_cols_p);

figure('Name','Parameter Sensitivity Curves','NumberTitle','off', ...
       'Position',[100 50 1000 600]);

for p = 1:n_params
    subplot(n_rows_p, n_cols_p, p);

    pv       = param_results(p).perturbed_values;
    pc       = param_results(p).percent_changes;
    baseline = param_results(p).baseline_val;

    plot(pv, pc, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerFaceColor','b');
    hold on;
    xline(baseline, 'k--', 'LineWidth', 1.2);
    yline(0,         'r--', 'LineWidth', 1.0);
    hold off;

    plabel = strrep(param_results(p).param_name, 'crop_elec_price', 'elec_price');
    xlabel(plabel, 'Interpreter','none');
    ylabel('% \Delta Social Value');
    title(plabel, 'Interpreter','none');
    grid on;
end

sgtitle('Parameter Sensitivity Curves');

%% Figure 4: Parameter Driver Ranking
max_abs_pct_param = zeros(n_params, 1);
for p = 1:n_params
    max_abs_pct_param(p) = max(abs(param_results(p).percent_changes));
end

[sorted_param, sidx_param] = sort(max_abs_pct_param, 'ascend');
sorted_names_param = cellfun(@(n) strrep(n, 'crop_elec_price', 'elec_price'), ...
    {param_results(sidx_param).param_name}, 'UniformOutput', false);

figure('Name','Parameter Driver Ranking','NumberTitle','off','Position',[250 250 700 400]);
barh(1:n_params, sorted_param, 0.6);
yticks(1:n_params);
yticklabels(sorted_names_param);
set(gca, 'TickLabelInterpreter', 'none');
xlabel('Max % Change in Social Value over Parameter Sweeps');
title('Parameter Sensitivity');
grid on;
ax = gca; ax.XGrid = 'on'; ax.YGrid = 'off';

%% Figure 5 — Normalized DV Sensitivity 
%Shows structural sensitivity independent of variable range size
mean_abs_norm = zeros(n_dv, 1);
for i = 1:n_dv
    ns    = dv_results(i).norm_sensitivities;
    valid = ns(~isnan(ns));
    if ~isempty(valid)
        mean_abs_norm(i) = mean(abs(valid));
    end
end

[sorted_norm, sidx_norm] = sort(mean_abs_norm, 'ascend');
sorted_names_norm = {dv_results(sidx_norm).var_name};

figure('Name','Normalized DV Sensitivity','NumberTitle','off','Position',[300 300 750 420]);
barh(1:n_dv, sorted_norm, 0.6);
yticks(1:n_dv);
yticklabels(sorted_names_norm);
xlabel('Mean Normalized Sensitivity  (% change per unit relative DV change)');
title('Design Variable Structural Sensitivity: Normalized Ranking');
grid on;
ax = gca; ax.XGrid = 'on'; ax.YGrid = 'off';

%% Figure 6 — DV Overlay: all 7 variables on one plot
% X-axis normalized to [0,1] (0 = lb, 1 = ub) so all variables share one axis.
% recovering lb/ub from the data: frac=1.0 in both directions guarantes the
% min/max of perturbed_values equals lb/ub

% Get f_star (social value at x*) from the stored results
[~, idx0_dv] = min(abs(dv_results(1).percent_changes));
f_star_M = dv_results(1).social_values(idx0_dv) / 1e6;

% Descriptive legend labels with units; phi/sigma shown in degrees
dv_descriptions = {'panel height', 'panel length', 'panel width', ...
                   'azimuth angle', 'tilt angle', 'row spacing', 'panel gap'};
dv_units        = {'m', 'm', 'm', 'deg', 'deg', 'm', 'm'};

figure('Name','DV Overlay: Social Value','NumberTitle','off','Position',[150 100 900 540]);
colors_dv = lines(n_dv);
hold on;
for i = 1:n_dv
    pv      = dv_results(i).perturbed_values;
    sv      = dv_results(i).social_values;
    xstar_v = dv_results(i).x_star_val;

    pv_min = min(pv);
    pv_max = max(pv);
    if pv_max > pv_min
        pv_norm = (pv - pv_min) / (pv_max - pv_min);
    else
        pv_norm = zeros(size(pv));
    end

    % Show x* value in degrees for angle variables
    if ismember(dv_results(i).var_name, angle_dv)
        xstar_display = rad2deg(xstar_v);
    else
        xstar_display = xstar_v;
    end
    leg_str = sprintf('%s: %s  (x* = %.3g %s)', ...
        dv_results(i).var_name, dv_descriptions{i}, xstar_display, dv_units{i});

    plot(pv_norm, sv / 1e6, '-o', ...
        'Color',           colors_dv(i,:), ...
        'LineWidth',        1.5, ...
        'MarkerSize',       5, ...
        'MarkerFaceColor',  colors_dv(i,:), ...
        'DisplayName',      leg_str);
end

% Baseline reference line (horizontal at f_star) --> dummy plot for legend entry
plot(nan, nan, 'r:', 'LineWidth', 1.5, ...
    'DisplayName', sprintf('optimized value baseline  (f* = $%.3f M)', f_star_M));
yline(f_star_M, 'r:', 'LineWidth', 1.5, 'HandleVisibility', 'off');

hold off;
xlabel('Normalized Design Variable Value  (0 = lower bound,  1 = upper bound)');
ylabel('Social Value  (M$)');
title('Design Variable Sensitivity — All Variables Overlay');
legend('Location','best', 'Interpreter','none');
grid on;

%% Figure 7 — Parameter Overlay: all parameters on one plot
%x-axis: percent change from baseline value.

param_baseline_display = { ...
    'crop_RUE',        '1.0',    ''; ...
    'elec_price',     '$0.50',  '$/kWh'; ...
    'crop_price',      '$16.03', '$/kg'; ...
    'PV_n_p',         '0.20',   ''; ...
    'crop_T_max',     '30',     'deg C'; ...
    'crop_GCF',        '0.30',   ''};

turquoise_p = 6;
yellow_p    = 3;

figure('Name','Parameter Overlay: Social Value','NumberTitle','off','Position',[200 150 900 540]);
colors_p = lines(n_params);
hold on;
for p = 1:n_params
    pv       = param_results(p).perturbed_values;
    sv       = param_results(p).social_values;
    baseline = param_results(p).baseline_val;

    pct_x = (pv - baseline) / baseline * 100;

    bval  = param_baseline_display{p, 2};
    bunit = param_baseline_display{p, 3};
    pname = param_baseline_display{p, 1};
    if isempty(bunit)
        leg_str = sprintf('%s  (baseline: %s)', pname, bval);
    else
        leg_str = sprintf('%s  (baseline: %s %s)', pname, bval, bunit);
    end

    if p == turquoise_p
        ls = '--o';
    else
        ls = '-o';
    end
    mk_size = 5;
    if p == yellow_p
        mk_size = 9;
    end

    plot(pct_x, sv / 1e6, ls, ...
        'Color',           colors_p(p,:), ...
        'LineWidth',        1.5, ...
        'MarkerSize',       mk_size, ...
        'MarkerFaceColor',  colors_p(p,:), ...
        'DisplayName',      leg_str);
end

% Baseline (x_star) reference line
plot(nan, nan, 'r:', 'LineWidth', 1.5, ...
    'DisplayName', sprintf('optimized value baseline  (f* = $%.3f M)', f_star_M));
yline(f_star_M, 'r:', 'LineWidth', 1.5, 'HandleVisibility', 'off');

hold off;
xlabel('% Change in Parameter Value from Baseline');
ylabel('Social Value  (M$)');
title('Parameter Sensitivity: All Parameters Overlay');
legend('Location','best', 'Interpreter','none');
grid on;


