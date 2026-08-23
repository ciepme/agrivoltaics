%% Dynamic Single Persona Pareto Front Deep Dive
clc; clear; close all;

%% ========================================================================
%% 1. CONFIGURATION (Edit these if your wrapper output or variables change)
%% ========================================================================

% --- METRICS COLUMNS ---
% Where are these metrics located in the 9-element output of agrivoltaic_wrapper?
COL_METRIC_EMISSIONS = 1;
COL_METRIC_PROFIT    = 2;
COL_METRIC_CROP      = 6;  % Biomass
COL_METRIC_ENERGY    = 8;
COL_METRIC_CAPEX     = 9;

% --- DESIGN VARIABLE COLUMNS ---
% Where are these variables located in your GA array 'ga_solve'?
% Default: [Height(1), Length(2), Width(3), Azimuth(4), Tilt(5), RowGap(6), PanelGap(7)]
COL_VAR_HEIGHT    = 1;
COL_VAR_LENGTH    = 2;
COL_VAR_WIDTH     = 3;
COL_VAR_AZIMUTH   = 4;
COL_VAR_TILT      = 5;
COL_VAR_ROW_GAP   = 6;
COL_VAR_PANEL_GAP = 7;

%% ========================================================================
%% 2. LOAD DATA
%% ========================================================================
fprintf('Please select a saved .mat file to analyze...\n');
[file_name, file_path] = uigetfile('*.mat', 'Select Persona Data');
if isequal(file_name,0)
    disp('User canceled.'); return;
end
data = load(fullfile(file_path, file_name));

fval           = data.fval;
pareto_metrics = data.pareto_metrics;
ga_solve       = data.ga_solve;

%% ========================================================================
%% 3. DYNAMICALLY FIND POINTS OF INTEREST
%% ========================================================================
% A. Max Emission Reduction (Find MAX in pareto_metrics)
[~, idx_emissions] = max(pareto_metrics(:, COL_METRIC_EMISSIONS));

% B. Max Profit (Find MAX in pareto_metrics)
[~, idx_profit] = max(pareto_metrics(:, COL_METRIC_PROFIT));

% C. Max Crop Yield (Find MAX in pareto_metrics)
[~, idx_crop] = max(pareto_metrics(:, COL_METRIC_CROP));

% D. Lowest Capex (Find MIN in pareto_metrics)
[~, idx_capex] = min(pareto_metrics(:, COL_METRIC_CAPEX));

[~, idx_capex] = min(pareto_metrics(:,  COL_METRIC_EMISSIONS));


% E. Best Compromise (Utopia Point)
min_vals = min(fval);
max_vals = max(fval);
range_vals = max_vals - min_vals;
range_vals(range_vals == 0) = 1; % Prevent divide-by-zero

norm_fval = (fval - min_vals) ./ range_vals; 
utopia_target = zeros(1, size(norm_fval, 2)); % 0 is the best (min) of every column
distances = sqrt(sum((norm_fval - utopia_target).^2, 2));
[~, idx_utopia] = min(distances);

% Group the indices
indices = [idx_profit, idx_crop, idx_emissions, idx_capex, idx_utopia];
design_names = categorical({'Max Profit', 'Max Crop Yield', 'Max Emission Red.', 'Lowest Capex', 'Compromise'});
design_names = reordercats(design_names, {'Max Profit', 'Max Crop Yield', 'Max Emission Red.', 'Lowest Capex', 'Compromise'});

%% ========================================================================
%% 4. EXTRACT DATA FOR PLOTTING
%% ========================================================================
% Preallocate arrays
profit_vals = zeros(1,5);
crop_vals   = zeros(1,5);
energy_vals = zeros(1,5);
capex_vals  = zeros(1,5);

height_vals    = zeros(1,5);
length_vals    = zeros(1,5);
width_vals     = zeros(1,5);
azimuth_vals   = zeros(1,5);
tilt_vals      = zeros(1,5);
row_gap_vals   = zeros(1,5);
panel_gap_vals = zeros(1,5);

for i = 1:5
    row = indices(i);
    
    % Extract Metrics
    profit_vals(i) = pareto_metrics(row, COL_METRIC_PROFIT) / 1e6; 
    crop_vals(i)   = pareto_metrics(row, COL_METRIC_CROP);         
    energy_vals(i) = pareto_metrics(row, COL_METRIC_ENERGY) / 1e6; 
    capex_vals(i)  = pareto_metrics(row, COL_METRIC_CAPEX) / 1e6;  
    
    % Extract All 7 Design Variables
    height_vals(i)    = ga_solve(row, COL_VAR_HEIGHT); 
    length_vals(i)    = ga_solve(row, COL_VAR_LENGTH); 
    width_vals(i)     = ga_solve(row, COL_VAR_WIDTH); 
    azimuth_vals(i)   = rad2deg(ga_solve(row, COL_VAR_AZIMUTH)); % Convert rad to deg
    tilt_vals(i)      = rad2deg(ga_solve(row, COL_VAR_TILT));    % Convert rad to deg
    row_gap_vals(i)   = ga_solve(row, COL_VAR_ROW_GAP); 
    panel_gap_vals(i) = ga_solve(row, COL_VAR_PANEL_GAP); 
end

%% ========================================================================
%% 5. PLOT METRICS (2x2 Grid)
%% ========================================================================
fig1 = figure('Name', 'Objective Comparisons', 'Color', 'w', 'Position', [100, 100, 1000, 800]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; bar(design_names, profit_vals, 'FaceColor', [0.18 0.53 0.30]);
title('Fiscal Profit', 'FontSize', 14, 'FontWeight', 'bold'); ylabel('Millions (USD)', 'FontSize', 12); grid on;

nexttile; bar(design_names, crop_vals, 'FaceColor', [0.85 0.33 0.10]);
title('Annual Crop Yield', 'FontSize', 14, 'FontWeight', 'bold'); ylabel('Biomass', 'FontSize', 12); grid on;

nexttile; bar(design_names, energy_vals, 'FaceColor', [0.93 0.69 0.13]);
title('Energy Generation', 'FontSize', 14, 'FontWeight', 'bold'); ylabel('GWh', 'FontSize', 12); grid on;

nexttile; bar(design_names, capex_vals, 'FaceColor', [0.49 0.18 0.56]); 
title('Capital Expenditure (Capex)', 'FontSize', 14, 'FontWeight', 'bold'); ylabel('Millions (USD)', 'FontSize', 12); grid on;

%% ========================================================================
%% 6. PLOT DESIGN CHOICES (2x4 Grid)
%% ========================================================================
fig2 = figure('Name', 'Design Variable Comparisons', 'Color', 'w', 'Position', [150, 150, 1600, 800]);
tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
bar_color = [0.2 0.4 0.6];

% 1. Panel Height
nexttile; bar(design_names, height_vals, 'FaceColor', bar_color);
title('Panel Height', 'FontSize', 14, 'FontWeight', 'bold'); ylabel('Meters (m)', 'FontSize', 12); grid on;

% 2. Panel Length
nexttile; bar(design_names, length_vals, 'FaceColor', bar_color);
title('Panel Length', 'FontSize', 14, 'FontWeight', 'bold'); ylabel('Meters (m)', 'FontSize', 12); grid on;

% 3. Panel Width
nexttile; bar(design_names, width_vals, 'FaceColor', bar_color);
title('Panel Width', 'FontSize', 14, 'FontWeight', 'bold'); ylabel('Meters (m)', 'FontSize', 12); grid on;

% 4. Row Gap
nexttile; bar(design_names, row_gap_vals, 'FaceColor', bar_color);
title('Row Gap', 'FontSize', 14, 'FontWeight', 'bold'); ylabel('Meters (m)', 'FontSize', 12); grid on;

% 5. Panel Gap
nexttile; bar(design_names, panel_gap_vals, 'FaceColor', bar_color);
title('Panel Gap', 'FontSize', 14, 'FontWeight', 'bold'); ylabel('Meters (m)', 'FontSize', 12); grid on;

% 6. Tilt Angle
nexttile; bar(design_names, tilt_vals, 'FaceColor', bar_color);
title('Tilt Angle', 'FontSize', 14, 'FontWeight', 'bold'); ylabel('Degrees (°)', 'FontSize', 12); grid on;

% 7. Azimuth Angle
nexttile; bar(design_names, azimuth_vals, 'FaceColor', bar_color);
title('Azimuth Angle', 'FontSize', 14, 'FontWeight', 'bold'); ylabel('Degrees (°)', 'FontSize', 12); grid on;

%% ========================================================================
%% 7. CONSOLE OUTPUT
%% ========================================================================
fprintf('\n=== DEEP DIVE RESULTS ===\n');
for i = 1:5
    fprintf('%s Design:\n', char(design_names(i)));
    fprintf('  - Capex:    $%.2f M\n', capex_vals(i));
    fprintf('  - Profit:   $%.2f M\n', profit_vals(i));
    fprintf('  - Crop:     %.2f kg\n', crop_vals(i));
    fprintf('  - Layout:   H: %.2fm | Gap: %.2fm | Tilt: %.1f° | Azim: %.1f°\n\n', ...
            height_vals(i), row_gap_vals(i), tilt_vals(i), azimuth_vals(i));
end