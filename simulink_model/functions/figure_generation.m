%% Poster Plotting Script
clear;
clc;
close all;

%% Load Data
fixed = load('FINALagrivoltaic_comparative_optimization_datapop200_gen200_fixedaxis_newALL.mat');
single = load('agrivoltaic_comparative_optimization_datapop200_gen50_singleaxis.mat');

if ~exist('graphs', 'dir')
    mkdir('graphs');
end

%% Column Definitions
% results = [E, P, social_cost, pv_revenue, crop_revenue, yearly_biomass, total_panels, yearly_energy]

col_emissions = 1;
col_profit    = 2;
col_social    = 3;
col_pvrev     = 4;
col_croprev   = 5;
col_crop      = 6;
col_panels    = 7;
col_energy    = 8;

metric_names = {'Profit', 'Energy', 'Emissions', 'Crop Yield'};
metric_cols  = [col_profit, col_energy, col_emissions, col_crop];

%% Align target order between fixed and single-axis
targets = string(fixed.targets_to_run);
single_targets = string(single.targets_to_run);

single_order = zeros(size(targets));

for i = 1:length(targets)
    idx = find(single_targets == targets(i));

    if isempty(idx)
        error('Target "%s" exists in fixed results but not single-axis results.', targets(i));
    end

    single_order(i) = idx;
end

fixed_results  = fixed.results_matrix;
single_results = single.results_matrix(single_order, :);

fixed_xbest  = fixed.x_best_set;
single_xbest = single.x_best_set(single_order, :);

n_targets = length(targets);

%% ============================================================
% Figure 1: Normalized Performance, All Designs Together
% ============================================================

plot_data = zeros(2*n_targets, length(metric_cols));
row_labels = strings(2*n_targets, 1);

for i = 1:n_targets

    plot_data(2*i-1, :) = fixed_results(i, metric_cols);
    plot_data(2*i,   :) = single_results(i, metric_cols);

    row_labels(2*i-1) = "Fixed-" + targets(i);
    row_labels(2*i)   = "Single-" + targets(i);

end

% Min-max normalize each metric across all fixed + single-axis designs.
% This handles negative profit correctly.
plot_data_norm = normalize_minmax(plot_data);

figure('Color','w', 'Position', [100 100 1250 500]);

bar(categorical(row_labels, row_labels), plot_data_norm, 'grouped');

ylabel('Normalized Score');
title('Normalized Performance by Objective and Tracking Type');
legend(metric_names, 'Location', 'bestoutside');
grid on;
ylim([0 1.1]);
xtickangle(35);

annotation('textbox', [0.13 0.01 0.75 0.08], ...
    'String', 'Min-max normalized within each metric: 0 = worst, 1 = best among all fixed and single-axis designs.', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 9);

saveas(gcf, 'graphs/fig1_normalized_all_designs.png');

%% ============================================================
% Figure 2: Metric Tiles, Fixed vs Single-Axis by Objective
% This is usually better for a poster than Figure 1.
% ============================================================

figure('Color','w', 'Position', [100 100 1200 750]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Normalize each metric across all fixed + single-axis designs.
all_metric_data = [
    fixed_results(:, metric_cols);
    single_results(:, metric_cols)
];

all_metric_norm = normalize_minmax(all_metric_data);

fixed_metric_norm  = all_metric_norm(1:n_targets, :);
single_metric_norm = all_metric_norm(n_targets+1:end, :);

for m = 1:length(metric_names)

    nexttile;

    data_norm = [fixed_metric_norm(:,m), single_metric_norm(:,m)];

    bar(categorical(targets, targets), data_norm);
    title(metric_names{m});
    ylabel('Normalized Score');
    legend({'Fixed Axis', 'Single Axis'}, 'Location', 'best');
    grid on;
    ylim([0 1.1]);

end

sgtitle('Fixed vs Single-Axis Performance by Optimization Objective');

saveas(gcf, 'graphs/fig2_metric_tiles_fixed_vs_single.png');

%% ============================================================
% Figure 3: Single-Axis Design Variable Heatmap
% Color = normalized value, text = actual value
% ============================================================

common_design_names = {'Height (m)', 'Length (m)', 'Width (m)', 'Row Spacing (m)', 'Panel Gap (m)'};
common_keep_indices = [1, 2, 3, 6, 7];

fixed_common_design  = fixed_xbest(:, common_keep_indices);
single_common_design = single_xbest(:, common_keep_indices);

combined_common_design = [
    fixed_common_design;
    single_common_design
];

common_min = min(combined_common_design, [], 1);
common_max = max(combined_common_design, [], 1);

fixed_common_norm  = normalize_with_reference(fixed_common_design,  common_min, common_max);
single_common_norm = normalize_with_reference(single_common_design, common_min, common_max);

shared_colormap = parula(256);

figure('Color','w', 'Position', [100 100 950 420]);

plot_labeled_heatmap( ...
    single_common_norm, ...
    single_common_design, ...
    targets, ...
    common_design_names, ...
    'Single-Axis Design Variables', ...
    shared_colormap, ...
    '%.2f');

saveas(gcf, 'graphs/fig3_single_axis_design_heatmap.png');

%% ============================================================
% Figure 4: Fixed-Axis Design Variable Heatmap
% Color = normalized value, text = actual value
% ============================================================

fixed_design_names = {'Height (m)', 'Length (m)', 'Width (m)', 'Azimuth (deg)', 'Tilt (deg)', 'Row Spacing (m)', 'Panel Gap (m)'};

fixed_design_raw = fixed_xbest(:, 1:7);

% Convert azimuth and tilt to degrees for labels
fixed_design_raw(:,4) = rad2deg(fixed_design_raw(:,4));
fixed_design_raw(:,5) = rad2deg(fixed_design_raw(:,5));

fixed_design_norm = zeros(size(fixed_design_raw));

% Common variables normalized using same fixed + single-axis reference as Figure 3
fixed_design_norm(:, [1,2,3,6,7]) = fixed_common_norm;

% Fixed-only variables normalized within fixed-axis results
fixed_design_norm(:, 4:5) = normalize_minmax(fixed_design_raw(:, 4:5));

figure('Color','w', 'Position', [100 100 1150 420]);

plot_labeled_heatmap( ...
    fixed_design_norm, ...
    fixed_design_raw, ...
    targets, ...
    fixed_design_names, ...
    'Fixed-Axis Design Variables', ...
    shared_colormap, ...
    '%.2f');

saveas(gcf, 'graphs/fig4_fixed_axis_design_heatmap.png');

%% ============================================================
% Figure 5: Correct Fixed vs Single-Axis Comparison Across All Objectives
% This fixes the issue where one selected-objective figure did not reach zero.
% ============================================================

figure('Color','w', 'Position', [100 100 1200 750]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for m = 1:length(metric_names)

    nexttile;

    data_norm = [fixed_metric_norm(:,m), single_metric_norm(:,m)];

    bar(categorical(targets, targets), data_norm);

    title(metric_names{m});
    ylabel('Normalized Score');
    legend({'Fixed Axis','Single Axis'}, 'Location', 'best');
    ylim([0 1.1]);
    grid on;

end

sgtitle('Globally Normalized Fixed vs Single-Axis Comparison');

saveas(gcf, 'graphs/fig5_global_fixed_vs_single_by_metric.png');

%% ============================================================
% Optional Figure 6: Actual Values Instead of Normalized Values
% Recommended for poster if normalized plots are too abstract.
% ============================================================

figure('Color','w', 'Position', [100 100 1200 750]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

actual_metric_data_fixed = fixed_results(:, metric_cols);
actual_metric_data_single = single_results(:, metric_cols);

metric_ylabels = {'Profit ($)', 'Energy (kWh/year)', 'Emissions Reduction', 'Crop Yield'};

for m = 1:length(metric_names)

    nexttile;

    data_actual = [actual_metric_data_fixed(:,m), actual_metric_data_single(:,m)];

    bar(categorical(targets, targets), data_actual);

    title(metric_names{m});
    ylabel(metric_ylabels{m});
    legend({'Fixed Axis','Single Axis'}, 'Location', 'best');
    grid on;

end

sgtitle('Actual Performance Values by Objective and Tracking Type');

saveas(gcf, 'graphs/fig6_actual_values_fixed_vs_single.png');

%% ============================================================
% Summary Table Export
% ============================================================

System = {};
Objective = {};
Profit = [];
Energy = [];
Emissions = [];
CropYield = [];
Panels = [];
PVRevenue = [];

for i = 1:n_targets

    System{end+1,1} = 'Fixed';
    Objective{end+1,1} = char(targets(i));
    Profit(end+1,1) = fixed_results(i, col_profit);
    Energy(end+1,1) = fixed_results(i, col_energy);
    Emissions(end+1,1) = fixed_results(i, col_emissions);
    CropYield(end+1,1) = fixed_results(i, col_crop);
    Panels(end+1,1) = fixed_results(i, col_panels);
    PVRevenue(end+1,1) = fixed_results(i, col_pvrev);

    System{end+1,1} = 'Single-Axis';
    Objective{end+1,1} = char(targets(i));
    Profit(end+1,1) = single_results(i, col_profit);
    Energy(end+1,1) = single_results(i, col_energy);
    Emissions(end+1,1) = single_results(i, col_emissions);
    CropYield(end+1,1) = single_results(i, col_crop);
    Panels(end+1,1) = single_results(i, col_panels);
    PVRevenue(end+1,1) = single_results(i, col_pvrev);

end

summary_table = table(System, Objective, Profit, Energy, Emissions, CropYield, Panels, PVRevenue);

disp(summary_table);

writetable(summary_table, 'graphs/fixed_single_summary_table.csv');

%% ============================================================
% Helper Functions
% ============================================================

function X_norm = normalize_minmax(X)
% Min-max normalize each column of X.
% 0 = minimum/worst, 1 = maximum/best.
% Works for negative values, including negative profit.

    X_norm = zeros(size(X));

    for j = 1:size(X,2)
        col = X(:,j);

        col_min = min(col);
        col_max = max(col);

        if abs(col_max - col_min) < 1e-12
            X_norm(:,j) = 1;
        else
            X_norm(:,j) = (col - col_min) ./ (col_max - col_min);
        end
    end

    X_norm = max(X_norm, 0);
    X_norm = min(X_norm, 1);
end

function X_norm = normalize_with_reference(X, ref_min, ref_max)
% Normalize X using externally supplied min and max values.

    X_norm = zeros(size(X));

    for j = 1:size(X,2)

        if abs(ref_max(j) - ref_min(j)) < 1e-12
            X_norm(:,j) = 1;
        else
            X_norm(:,j) = (X(:,j) - ref_min(j)) ./ (ref_max(j) - ref_min(j));
        end

    end

    X_norm = max(X_norm, 0);
    X_norm = min(X_norm, 1);
end

function plot_labeled_heatmap(C_norm, C_actual, row_labels, col_labels, plot_title, cmap, value_format)
% Creates a normalized heatmap with actual numeric values written in cells.
% C_norm controls the color.
% C_actual controls the displayed labels.

    imagesc(C_norm);

    colormap(cmap);
    clim([0 1]);
    colorbar;

    ax = gca;
    ax.XTick = 1:length(col_labels);
    ax.XTickLabel = col_labels;
    ax.YTick = 1:length(row_labels);
    ax.YTickLabel = row_labels;
    ax.TickLength = [0 0];

    xtickangle(35);

    title(plot_title);
    xlabel('Design Variable');
    ylabel('Optimization Objective');

    for r = 1:size(C_norm,1)
        for c = 1:size(C_norm,2)

            val = C_actual(r,c);
            label = sprintf(value_format, val);

            if C_norm(r,c) > 0.55
                text_color = 'w';
            else
                text_color = 'k';
            end

            text(c, r, label, ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'Color', text_color, ...
                'FontWeight', 'bold', ...
                'FontSize', 9);
        end
    end
end