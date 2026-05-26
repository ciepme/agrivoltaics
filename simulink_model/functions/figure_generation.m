%% Poster Plotting Script Using master_best_designs.mat

clear;
clc;
close all;

%% Load Data

master = load("results/master_best_designs.mat");

best_table = master.best_table;
plot_objectives = string(master.plot_objectives(:));

fixed_results_all = master.fixed_results_matrix;
single_results_all = master.single_results_matrix;

fixed_xbest_cell_all = master.fixed_x_best_cell;
single_xbest_cell_all = master.single_x_best_cell;

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

metric_ylabels = {'Profit ($)', 'Energy (kWh/year)', 'Emissions Reduction', 'Crop Yield'};

%% Keep only objectives that exist for both fixed and single-axis

%% Keep objectives that have at least some data for fixed or single-axis

has_fixed_data = ~all(isnan(fixed_results_all(:, metric_cols)), 2);
has_single_data = ~all(isnan(single_results_all(:, metric_cols)), 2);

valid_metric_rows = has_fixed_data | has_single_data;

if ~any(valid_metric_rows)
    error("No objectives have usable fixed-axis or single-axis metric results. Check master_best_designs.mat.");
end

targets = plot_objectives(valid_metric_rows);

fixed_results = fixed_results_all(valid_metric_rows, :);
single_results = single_results_all(valid_metric_rows, :);

fixed_xbest_cell = fixed_xbest_cell_all(valid_metric_rows);
single_xbest_cell = single_xbest_cell_all(valid_metric_rows);

n_targets = length(targets);

fprintf("\nObjectives included in plotting:\n");
disp(targets);

fprintf("\nFixed result rows used for plotting:\n");
disp(fixed_results(:, metric_cols));

fprintf("\nSingle-axis result rows used for plotting:\n");
disp(single_results(:, metric_cols));

if ~any(valid_metric_rows)
    error("No objectives have both fixed-axis and single-axis metric results. Check master_best_designs.mat.");
end

targets = plot_objectives(valid_metric_rows);

fixed_results = fixed_results_all(valid_metric_rows, :);
single_results = single_results_all(valid_metric_rows, :);

fixed_xbest_cell = fixed_xbest_cell_all(valid_metric_rows);
single_xbest_cell = single_xbest_cell_all(valid_metric_rows);

n_targets = length(targets);
% Pretty display names for figures
target_labels = prettify_objective_labels(targets);

fprintf("\nObjectives included in plotting:\n");
disp(targets);

%% Convert design-vector cells into layout matrices

% First 7 variables are:
% [Height, Length, Width, Azimuth, Tilt, Row Spacing, Panel Gap]
fixed_xbest_7 = cell_design_to_layout_matrix(fixed_xbest_cell, 7);
single_xbest_7 = cell_design_to_layout_matrix(single_xbest_cell, 7);

%% ============================================================
% Figure 1: Normalized Performance, All Designs Together
% ============================================================

plot_data = zeros(2*n_targets, length(metric_cols));
row_labels = strings(2*n_targets, 1);

for i = 1:n_targets

    plot_data(2*i-1, :) = fixed_results(i, metric_cols);
    plot_data(2*i,   :) = single_results(i, metric_cols);

    row_labels(2*i-1) = "Fixed-" + target_labels(i);
    row_labels(2*i)   = "Single-" + target_labels(i);

end

plot_data_norm = normalize_minmax(plot_data);

figure('Color','w', 'Position', [100 100 1350 520]);

bar(categorical(row_labels, row_labels), plot_data_norm, 'grouped');

ylabel('Normalized Score');
title('Normalized Performance by Objective and Tracking Type');
legend(metric_names, 'Location', 'bestoutside');
grid on;
ylim([0 1.1]);
xtickangle(35);

annotation('textbox', [0.13 0.01 0.75 0.08], ...
    'String', 'Min-max normalized within each metric: 0 = worst, 1 = best among displayed designs.', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 9);

saveas(gcf, 'graphs/fig1_normalized_all_designs.png');

%% ============================================================
% Figure 2: Metric Tiles, Fixed vs Single-Axis by Objective
% ============================================================

figure('Color','w', 'Position', [100 100 1200 750]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

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

    bar(categorical(target_labels, target_labels), data_norm);

    title(metric_names{m});
    ylabel('Normalized Score');
    legend({'Fixed Axis', 'Single Axis'}, 'Location', 'best');
    grid on;
    ylim([0 1.1]);
    xtickangle(25);

end

sgtitle('Fixed vs Single-Axis Performance by Optimization Objective');

saveas(gcf, 'graphs/fig2_metric_tiles_fixed_vs_single.png');

%% ============================================================
% Figure 3: Single-Axis Design Variable Heatmap
% Includes MOGA-Balanced if present in master file
% ============================================================

common_design_names = {'Height (m)', 'Length (m)', 'Width (m)', 'Row Spacing (m)', 'Panel Gap (m)'};
common_keep_indices = [1, 2, 3, 6, 7];

fixed_common_design  = fixed_xbest_7(:, common_keep_indices);
single_common_design = single_xbest_7(:, common_keep_indices);

valid_design_rows = ...
    ~any(isnan(fixed_common_design), 2) & ...
    ~any(isnan(single_common_design), 2);

if ~any(valid_design_rows)
    warning("No valid design rows found for heatmaps. Skipping design heatmaps.");
else

    heatmap_targets = target_labels(valid_design_rows);

    fixed_common_design_hm = fixed_common_design(valid_design_rows, :);
    single_common_design_hm = single_common_design(valid_design_rows, :);

    combined_common_design = [
        fixed_common_design_hm;
        single_common_design_hm
    ];

    common_min = min(combined_common_design, [], 1);
    common_max = max(combined_common_design, [], 1);

    fixed_common_norm = normalize_with_reference(fixed_common_design_hm, common_min, common_max);
    single_common_norm = normalize_with_reference(single_common_design_hm, common_min, common_max);

    shared_colormap = parula(256);

    figure('Color','w', 'Position', [100 100 1000 480]);

    plot_labeled_heatmap( ...
        single_common_norm, ...
        single_common_design_hm, ...
        heatmap_targets, ...
        common_design_names, ...
        'Single-Axis Design Variables', ...
        shared_colormap, ...
        '%.2f');

    saveas(gcf, 'graphs/fig3_single_axis_design_heatmap.png');

    %% ============================================================
    % Figure 4: Fixed-Axis Design Variable Heatmap
    % Includes MOGA-Balanced if present in master file
    % ============================================================

    fixed_design_names = {'Height (m)', 'Length (m)', 'Width (m)', 'Azimuth (deg)', 'Tilt (deg)', 'Row Spacing (m)', 'Panel Gap (m)'};

    fixed_design_raw = fixed_xbest_7(valid_design_rows, 1:7);

    % Convert azimuth and tilt to degrees for labels
    fixed_design_raw(:,4) = rad2deg(fixed_design_raw(:,4));
    fixed_design_raw(:,5) = rad2deg(fixed_design_raw(:,5));

    fixed_design_norm = zeros(size(fixed_design_raw));

    % Common variables normalized using same fixed + single-axis reference as Figure 3
    fixed_design_norm(:, [1,2,3,6,7]) = fixed_common_norm;

    % Fixed-only variables normalized within fixed-axis results
    fixed_design_norm(:, 4:5) = normalize_minmax(fixed_design_raw(:, 4:5));

    figure('Color','w', 'Position', [100 100 1200 480]);

    plot_labeled_heatmap( ...
        fixed_design_norm, ...
        fixed_design_raw, ...
        heatmap_targets, ...
        fixed_design_names, ...
        'Fixed-Axis Design Variables', ...
        shared_colormap, ...
        '%.2f');

    saveas(gcf, 'graphs/fig4_fixed_axis_design_heatmap.png');

end
%% ============================================================
    % Export Fixed-Axis Design Variables Table
    % ============================================================
    % Clean column headers for a standard spreadsheet layout
    table_headers = {'Optimization_Objective', 'Height_m', 'Length_m', 'Width_m', 'Azimuth_deg', 'Tilt_deg', 'Row_Spacing_m', 'Panel_Gap_m'};
    
    % Construct the MATLAB table using the variables calculated for Figure 4
    fixed_design_table = table(cellstr(heatmap_targets), ...
        fixed_design_raw(:,1), ...
        fixed_design_raw(:,2), ...
        fixed_design_raw(:,3), ...
        fixed_design_raw(:,4), ...
        fixed_design_raw(:,5), ...
        fixed_design_raw(:,6), ...
        fixed_design_raw(:,7), ...
        'VariableNames', table_headers);
    
    % Display table to the command window for quick validation
    fprintf('\n=== Fixed-Axis Design Variables Matrix ===\n');
    disp(fixed_design_table);
    
    % Export directly to a spreadsheet-ready CSV file
    writetable(fixed_design_table, 'graphs/fixed_axis_design_variables_table.csv');
    fprintf('Table successfully saved to: graphs/fixed_axis_design_variables_table.csv\n');
%% ============================================================
% Figure 5: Globally Normalized Fixed vs Single-Axis Comparison
% This duplicates Figure 2 visually, but preserves your original numbering.
% ============================================================

figure('Color','w', 'Position', [100 100 1200 750]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for m = 1:length(metric_names)

    nexttile;

    data_norm = [fixed_metric_norm(:,m), single_metric_norm(:,m)];

    bar(categorical(target_labels, target_labels), data_norm);

    title(metric_names{m});
    ylabel('Normalized Score');
    legend({'Fixed Axis','Single Axis'}, 'Location', 'best');
    ylim([0 1.1]);
    grid on;
    xtickangle(25);

end

sgtitle('Globally Normalized Fixed vs Single-Axis Comparison');

saveas(gcf, 'graphs/fig5_global_fixed_vs_single_by_metric.png');

%% ============================================================
% Figure 6: Actual Values Instead of Normalized Values
% ============================================================

figure('Color','w', 'Position', [100 100 1200 750]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

actual_metric_data_fixed = fixed_results(:, metric_cols);
actual_metric_data_single = single_results(:, metric_cols);

for m = 1:length(metric_names)

    nexttile;

    data_actual = [actual_metric_data_fixed(:,m), actual_metric_data_single(:,m)];

    bar(categorical(targets, targets), data_actual);

    title(metric_names{m});
    ylabel(metric_ylabels{m});
    legend({'Fixed Axis','Single Axis'}, 'Location', 'best');
    grid on;
    xtickangle(25);

end

sgtitle('Actual Performance Values by Objective and Tracking Type');

saveas(gcf, 'graphs/fig6_actual_values_fixed_vs_single.png');

%% ============================================================
% Figure 7: Single-Axis Tracking Curves
% Includes all available single-axis objectives, including MOGA-Balanced
% ============================================================

tracking_rows = false(n_targets, 1);

for i = 1:n_targets
    x_i = single_xbest_cell{i};
    tracking_rows(i) = ~isempty(x_i) && numel(x_i) >= 103;
end

if any(tracking_rows)

    tracking_targets = target_labels(tracking_rows);
    tracking_x_cells = single_xbest_cell(tracking_rows);

    seasons = {'Spring', 'Summer', 'Fall', 'Winter'};
    hours = 1:24;

    colors = lines(length(tracking_targets));

    figure('Color','w', 'Position', [200 200 1100 850]);
    tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    for s = 1:4

        nexttile;
        hold on;
        grid on;

        start_idx = 7 + (s-1)*24 + 1;
        end_idx   = 7 + s*24;

        for i = 1:length(tracking_targets)

            x_i = tracking_x_cells{i};
            angles_deg = rad2deg(x_i(start_idx:end_idx));

            plot(hours, angles_deg, '-o', ...
                'LineWidth', 2, ...
                'MarkerSize', 4, ...
                'Color', colors(i,:), ...
                'DisplayName', tracking_targets(i));

        end

        title(seasons{s}, 'FontWeight', 'bold');
        xlabel('Hour of Day');
        ylabel('Panel Angle (deg)');
        xlim([1 24]);
        xticks(1:2:24);

        if s == 1
            legend('Location', 'best');
        end

    end

    sgtitle('Single-Axis Tracking Curves by Optimization Objective');

    saveas(gcf, 'graphs/fig7_single_axis_tracking_curves.png');

else

    warning("No 103-variable single-axis designs found. Skipping tracking curve plot.");

end

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

function X7 = cell_design_to_layout_matrix(x_cell, n_layout_vars)
% Converts a cell array of design vectors into an N x n_layout_vars matrix.
% If a cell is empty or too short, fills that row with NaN.

    n = length(x_cell);
    X7 = NaN(n, n_layout_vars);

    for i = 1:n
        x_i = x_cell{i};

        if ~isempty(x_i) && numel(x_i) >= n_layout_vars
            X7(i, :) = x_i(1:n_layout_vars);
        end
    end
end

function X_norm = normalize_minmax(X)
% Min-max normalize each column of X.
% 0 = minimum/worst, 1 = maximum/best.
% Handles NaN values by ignoring them.

    X_norm = NaN(size(X));

    for j = 1:size(X,2)

        col = X(:,j);
        valid = ~isnan(col);

        if ~any(valid)
            continue;
        end

        col_min = min(col(valid));
        col_max = max(col(valid));

        if abs(col_max - col_min) < 1e-12
            X_norm(valid,j) = 1;
        else
            X_norm(valid,j) = (col(valid) - col_min) ./ (col_max - col_min);
        end
    end

    X_norm = max(X_norm, 0);
    X_norm = min(X_norm, 1);
end

function X_norm = normalize_with_reference(X, ref_min, ref_max)
% Normalize X using externally supplied min and max values.

    X_norm = NaN(size(X));

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

            if isnan(val)
                label = "N/A";
            else
                label = sprintf(value_format, val);
            end

            if C_norm(r,c) > 0.55
                text_color = 'k';
            else
                text_color = 'w';
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

function labels = prettify_objective_labels(targets)
% Converts internal objective names to nicer figure labels.

    targets = string(targets);
    labels = strings(size(targets));

    for i = 1:numel(targets)

        switch upper(strtrim(targets(i)))

            case "PROFIT"
                labels(i) = "Economic";

            case "EMISSIONS"
                labels(i) = "Environmental";

            case "POWER"
                labels(i) = "Energy Production";

            case "CROP"
                labels(i) = "Agricultural";

            case "MOGA-BALANCED"
                labels(i) = "Balanced";

            otherwise
                labels(i) = targets(i);

        end
    end
end
