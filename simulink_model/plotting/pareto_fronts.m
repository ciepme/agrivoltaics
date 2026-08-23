%% Plot Pareto Fronts for All Personas: (CapEx vs Profit) & (Emissions vs Crop)
clear; clc; close all;

% --- CONFIGURATION ---
% Set to true to scale values per acre so all 3 personas fit on the same axes nicely. 
% Set to false to plot absolute, whole-farm raw metrics.
normalize_per_acre = true; 

%% 1. Load the Data
% fprintf('Please select the saved .mat file for the DEFAULT persona...\n');
% [file_def, path_def] = uigetfile('*.mat', 'Select DEFAULT Pareto Data');
% def_data = load(fullfile(path_def, file_def));
% 
% fprintf('Please select the saved .mat file for FRED...\n');
% [file_fred, path_fred] = uigetfile('*.mat', 'Select FRED Pareto Data');
% fred_data = load(fullfile(path_fred, file_fred));
% 
% fprintf('Please select the saved .mat file for JANICE...\n');
% [file_jan, path_jan] = uigetfile('*.mat', 'Select JANICE Pareto Data');
% jan_data = load(fullfile(path_jan, file_jan));

def_data = load('pop200_max_gen150_mode0_persona3FINAL_MIXED_DEFAULT.mat');

fred_data = load('pop200_max_gen150_mode0_persona1FINAL_MIXED_FRED.mat');

jan_data = load('pop200_max_gen150_mode0_persona2FINAL_MIXED_JANICE.mat');

% Extract the pareto metrics matrices
metrics_def = def_data.pareto_metrics;
metrics_fred = fred_data.pareto_metrics;
metrics_jan = jan_data.pareto_metrics;

%% 2. Scaling Factors
% [Emissions(1), Profit(2), SocialCost(3), PVRev(4), CropRev(5), Biomass(6), Panels(7), Energy(8), Capex(9)]
if normalize_per_acre
    % Farm Area Calculations (1 acre = 4046.86 m^2)
    acres_def = (50 * 50) / 4046.86;        
    acres_fred = (2000 * 2000) / 4046.86;    
    acres_jan = (250 * 250) / 4046.86;      
    
    BASE_CAPEX_USD = 1000000; % Assuming metric(9) is a % of Base CapEx like before
    
    % --- Apply Scaling (Per Acre) ---
    % 1. Default
    def_capex = (metrics_def(:, 9) ./ 100) .* BASE_CAPEX_USD ./ acres_def;
    def_profit = metrics_def(:, 2) ./ acres_def;
    def_emissions = abs(metrics_def(:, 1)) ./ acres_def;
    def_crop = metrics_def(:, 6) .* (4046.86 / 1000); % kg per acre
    
    % 2. Fred
    fred_capex = (metrics_fred(:, 9) ./ 100) .* BASE_CAPEX_USD ./ acres_fred;
    fred_profit = metrics_fred(:, 2) ./ acres_fred;
    fred_emissions = abs(metrics_fred(:, 1)) ./ acres_fred;
    fred_crop = metrics_fred(:, 6) .* (4046.86 / 1000);
    
    % 3. Janice
    jan_capex = (metrics_jan(:, 9) ./ 100) .* BASE_CAPEX_USD ./ acres_jan;
    jan_profit = metrics_jan(:, 2) ./ acres_jan;
    jan_emissions = abs(metrics_jan(:, 1)) ./ acres_jan;
    jan_crop = metrics_jan(:, 6) .* (4046.86 / 1000);
    
    % Labels for axes
    lbl_profit = 'Lifetime Profit ($ / Acre)';
    lbl_capex = 'Capital Expenditure ($ / Acre)';
    lbl_emissions = 'Emissions Reduction (Units / Acre)';
    lbl_crop = 'Annual Crop Yield (kg / Acre)';
else
    % --- Use Raw Farm Data ---
    % 1. Default
    def_capex = metrics_def(:, 9); def_profit = metrics_def(:, 2);
    def_emissions = abs(metrics_def(:, 1)); def_crop = metrics_def(:, 6);
    
    % 2. Fred
    fred_capex = metrics_fred(:, 9); fred_profit = metrics_fred(:, 2);
    fred_emissions = abs(metrics_fred(:, 1)); fred_crop = metrics_fred(:, 6);
    
    % 3. Janice
    jan_capex = metrics_jan(:, 9); jan_profit = metrics_jan(:, 2);
    jan_emissions = abs(metrics_jan(:, 1)); jan_crop = metrics_jan(:, 6);
    
    % Labels for axes
    lbl_profit = 'Total Lifetime Profit ($)';
    lbl_capex = 'Total CapEx Metric';
    lbl_emissions = 'Total Emissions Reduction';
    lbl_crop = 'Total Crop Yield';
end

% Set up persona colors
color_def  = [0.4 0.4 0.4]; % Gray
color_fred = [0.85 0.33 0.10]; % Orange/Rust
color_jan  = [0.18 0.53 0.30]; % Green

%% 3. FIGURE 1: CapEx vs Profit
fig1 = figure('Name', 'Pareto Fronts: CapEx vs Profit', 'Color', 'w', 'Position', [100, 100, 800, 600]);
hold on; grid on;

% Scatter plots
s1 = scatter(def_capex, def_profit, 60, color_def, 'filled', 'MarkerEdgeColor', 'k');
s2 = scatter(fred_capex, fred_profit, 60, color_fred, 'filled', 'MarkerEdgeColor', 'k');
s3 = scatter(jan_capex, jan_profit, 60, color_jan, 'filled', 'MarkerEdgeColor', 'k');

% Add trendlines to connect the Pareto fronts
plot_trendline(def_capex, def_profit, color_def);
plot_trendline(fred_capex, fred_profit, color_fred);
plot_trendline(jan_capex, jan_profit, color_jan);

% Formatting
title('Pareto Fronts: CapEx vs. Profit', 'FontSize', 20, 'FontWeight', 'bold');
xlabel(lbl_capex, 'FontSize', 20, 'FontWeight', 'bold');
ylabel(lbl_profit, 'FontSize', 20, 'FontWeight', 'bold');
legend([s1, s2, s3], {'Default', 'Fred', 'Janice'}, ...
    'Location', 'best', 'FontSize', 18);
set(gca, 'FontSize', 12, 'LineWidth', 1.2);
hold off;

%% 4. FIGURE 2: Emissions vs Crop Yield
fig2 = figure('Name', 'Pareto Fronts: Emissions vs Crop', 'Color', 'w', 'Position', [150, 150, 800, 600]);
hold on; grid on;

% Scatter plots
s4 = scatter(def_emissions, def_crop, 60, color_def, 'filled', 'MarkerEdgeColor', 'k');
s5 = scatter(fred_emissions, fred_crop, 60, color_fred, 'filled', 'MarkerEdgeColor', 'k');
s6 = scatter(jan_emissions, jan_crop, 60, color_jan, 'filled', 'MarkerEdgeColor', 'k');

% Add trendlines to connect the Pareto fronts
plot_trendline(def_emissions, def_crop, color_def);
plot_trendline(fred_emissions, fred_crop, color_fred);
plot_trendline(jan_emissions, jan_crop, color_jan);

% Formatting
title('Pareto Fronts: Emissions Reduction vs. Crop Yield', 'FontSize', 20, 'FontWeight', 'bold');
xlabel(lbl_emissions, 'FontSize', 20, 'FontWeight', 'bold');
ylabel(lbl_crop, 'FontSize', 20, 'FontWeight', 'bold');
legend([s4, s5, s6], {'Default', 'Fred', 'Janice'}, ...
    'Location', 'best', 'FontSize', 18);
set(gca, 'FontSize', 20, 'LineWidth', 1.2);
hold off;

%% Helper Function to Draw Pareto Lines
function plot_trendline(x, y, color)
    % Sort by X-axis to draw a clean line left-to-right
    [sorted_x, sort_idx] = sort(x);
    sorted_y = y(sort_idx);
    plot(sorted_x, sorted_y, '--', 'Color', color, 'LineWidth', 1.5, 'HandleVisibility', 'off');
end