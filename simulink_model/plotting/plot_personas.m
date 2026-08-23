%% Combined Persona-Based Poster Plotting Script
clear; clc; close all;

%% ==============================================================
%% 1. FILE LOADING & DATA EXTRACTION
%% ==============================================================
% fprintf('Please select the saved .mat file for the DEFAULT persona...\n');
% [file_def, path_def] = uigetfile('*.mat', 'Select DEFAULT Persona Data');
%default_data = load(fullfile(path_def, file_def));
default_data = load('pop200_max_gen150_mode0_persona3FINAL_MIXED_DEFAULT.mat');

% fprintf('Please select the saved .mat file for FRED...\n');
% [file_fred, path_fred] = uigetfile('*.mat', 'Select FRED Data');
%fred_data = load(fullfile(path_fred, file_fred));
fred_data = load('pop200_max_gen150_mode0_persona1FINAL_MIXED_FRED.mat');

% fprintf('Please select the saved .mat file for JANICE...\n');
% [file_jan, path_jan] = uigetfile('*.mat', 'Select JANICE Data');

janice_data = load('pop200_max_gen150_mode0_persona2FINAL_MIXED_JANICE.mat');

% Extract Utopia Point Metrics AND Design Variables
[def_metrics, def_vars]       = get_utopia_data(default_data);
[fred_metrics, fred_vars]     = get_utopia_data(fred_data);
[janice_metrics, janice_vars] = get_utopia_data(janice_data);

%% ==============================================================
%% 2. DATA PREPARATION (METRICS)
%% ==============================================================
% Farm Area Calculations (1 acre = 4046.86 m^2)
acres_def    = (50 * 50) / 4046.86;        % ~0.62 acres
acres_fred   = (2000 * 2000) / 4046.86;    % ~988.42 acres
acres_janice = (250 * 250) / 4046.86;      % ~15.44 acres
acres_array = [acres_def, acres_fred, acres_janice];

personas = {'Default (Baseline)', 'Fred (Commodity)', 'Janice (Specialty)'};
n_personas = length(personas);
persona_names_cat = categorical(personas);
persona_names_cat = reordercats(persona_names_cat, personas);

% 1. Total Lifetime Profit (Raw $)
total_profit = [def_metrics(2), fred_metrics(2), janice_metrics(2)];

% 2. Annualized Profit Per Acre ($/acre/year)
profit_per_acre_yr = total_profit ./ acres_array ./ 30; 

% 3. Crop Yield (kg/acre/year)
conversion_factor = 4046.86 / 1000;
crop_per_acre_yr = [def_metrics(6), fred_metrics(6), janice_metrics(6)] .* conversion_factor;             

% 4. Energy (MWh/acre/year)
energy_per_acre_yr = [def_metrics(8), fred_metrics(8), janice_metrics(8)] ./ acres_array ./ 1e6; 

% 5. CapEx ($/acre) - metric(9) already outputs total absolute dollars
total_capex_vals = [def_metrics(9), fred_metrics(9), janice_metrics(9)];
capex_per_acre = total_capex_vals ./ acres_array; 

% 6. Panels Per Acre
panels_per_acre = [def_metrics(7), fred_metrics(7), janice_metrics(7)] ./ acres_array;

% 7. Emissions per Acre (assumes metric 1 is emissions)
emissions_per_acre = abs([def_metrics(1), fred_metrics(1), janice_metrics(1)]) ./ acres_array;

if ~exist('graphs', 'dir')
    mkdir('graphs');
end

%% ==============================================================
%% 3. FIGURE 1: OBJECTIVES PER ACRE (Original Bar Graphs)
%% ==============================================================
fig1 = figure('Name', 'Persona Comparison (Objectives Per Acre)', 'Color', 'w', 'Position', [100, 100, 1400, 500]);
tiledlayout(1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

color_profit = [0.18 0.53 0.30]; 
color_crop   = [0.85 0.33 0.10]; 
color_energy = [0.93 0.69 0.13]; 
color_capex  = [0.64 0.08 0.18]; 

% Profit
nexttile;
b1 = bar(persona_names_cat, profit_per_acre_yr, 'FaceColor', 'flat'); b1.CData = repmat(color_profit, 3, 1);
title('Annualized Profit', 'FontSize', 16, 'FontWeight', 'bold'); ylabel('USD / Acre / Year', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'LineWidth', 1.2); grid on;

% Crop
nexttile;
b2 = bar(persona_names_cat, crop_per_acre_yr, 'FaceColor', 'flat'); b2.CData = repmat(color_crop, 3, 1);
title('Annual Crop Yield', 'FontSize', 16, 'FontWeight', 'bold'); ylabel('kg / Acre / Year', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'LineWidth', 1.2); grid on;

% Energy
nexttile;
b3 = bar(persona_names_cat, energy_per_acre_yr, 'FaceColor', 'flat'); b3.CData = repmat(color_energy, 3, 1);
title('Energy Generation', 'FontSize', 16, 'FontWeight', 'bold'); ylabel('MWh / Acre / Year', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'LineWidth', 1.2); grid on;

% CapEx
nexttile;
b4 = bar(persona_names_cat, capex_per_acre, 'FaceColor', 'flat'); b4.CData = repmat(color_capex, 3, 1);
title('Upfront Capital Cost', 'FontSize', 16, 'FontWeight', 'bold'); ylabel('USD / Acre', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'LineWidth', 1.2); grid on;

saveas(fig1, "graphs/poster_objectives_comparison.png");

%% ==============================================================
%% 4. FIGURE 2: ALL 7 DESIGN VARIABLES (Original Bar Graphs)
%% ==============================================================
var_names = {'Panel Height (m)', 'Panel Length (m)', 'Panel Width (m)', ...
             'Azimuth Angle (°)', 'Tilt Angle (°)', 'Row Gap (m)', 'Panel Gap (m)'};

fig2 = figure('Name', 'Optimal Design Choices', 'Color', 'w', 'Position', [150, 150, 1400, 800]);
tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
color_design = [0.2 0.4 0.6]; 

for i = 1:7
    nexttile;
    vals = [def_vars(i), fred_vars(i), janice_vars(i)];

    % Convert Azimuth (var 4) and Tilt (var 5) from radians to degrees
    if i == 4 || i == 5
        vals = rad2deg(vals);
    end

    bd = bar(persona_names_cat, vals, 'FaceColor', 'flat'); 
    bd.CData = repmat(color_design, 3, 1);
    title(var_names{i}, 'FontSize', 14, 'FontWeight', 'bold');
    set(gca, 'FontSize', 12, 'LineWidth', 1.2); grid on;
end

saveas(fig2, "graphs/poster_design_variables_comparison.png");

%% ==============================================================
%% 5. FIGURE 3: TOTAL LIFETIME PROFIT (Original Bar Graph)
%% ==============================================================
fig3 = figure('Name', 'Total Lifetime Profit', 'Color', 'w', 'Position', [200, 200, 600, 500]);
b_tot = bar(persona_names_cat, total_profit ./ 1e6, 'FaceColor', 'flat'); 
b_tot.CData = repmat(color_profit, 3, 1);
title('Total Lifetime Farm Profit', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Millions (USD)', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'LineWidth', 1.2); grid on;

saveas(fig3, "graphs/poster_total_profit.png");

%% ==============================================================
%% 6. FIGURE 4: PERSONA DESIGN VARIABLE HEATMAP
%% ==============================================================
% Build Matrix (Rows match the 'personas' array: Default, Fred, Janice)
design_vars_actual = [
    def_vars(1:7); 
    fred_vars(1:7); 
    janice_vars(1:7)
];
% Convert Rad to Deg for the heatmap
design_vars_actual(:,4:5) = rad2deg(design_vars_actual(:,4:5)); 

design_vars_norm = normalize_minmax(design_vars_actual);
design_names = {'Height (m)', 'Length (m)', 'Width (m)', 'Azimuth (°)', 'Tilt (°)', 'Row Spacing (m)', 'Panel Gap (m)'};

fig4 = figure('Color','w', 'Position', [250, 250, 1000, 400]);
shared_colormap = parula(256);

plot_labeled_heatmap( ...
    design_vars_norm, ...
    design_vars_actual, ...
    personas, ...
    design_names, ...
    'Optimal Design Variables by Persona', ...
    shared_colormap, ...
    '%.2f');

saveas(fig4, 'graphs/poster_persona_design_heatmap.png');

%% ==============================================================
%% 7. FIGURE 5: PERSONA TRADE-OFF SPIDER PLOT
%% ==============================================================
% 1. Convert units specifically for the spider plot so it reads cleanly
% Profit/acre and CapEx/acre: divide by 1000 to get thousands ($k)
% Emissions: divide kg/acre by 1000 to get metric tons (t/ac)
% Total Profit: divide by 1,000,000 to get Millions ($M)
profit_spider = profit_per_acre_yr ./ 1000; 
capex_spider = capex_per_acre ./ 1000;
emissions_spider = emissions_per_acre ./ 1000; 
total_profit_spider = total_profit ./ 1e6; 

% 2. Build the matrix with the 6 metrics 
% Columns: [Crop, Emissions Red., Profit/ac, Total Profit, CapEx/ac, Panels/ac]
metrics_actual = [
    crop_per_acre_yr(1), emissions_spider(1), profit_spider(1), total_profit_spider(1), capex_spider(1), panels_per_acre(1);
    crop_per_acre_yr(2), emissions_spider(2), profit_spider(2), total_profit_spider(2), capex_spider(2), panels_per_acre(2);
    crop_per_acre_yr(3), emissions_spider(3), profit_spider(3), total_profit_spider(3), capex_spider(3), panels_per_acre(3)
];

% Update labels to reflect the 6 axes and new units
metric_names_spider = {'Crop Yield (kg/ac)', 'Emissions (t/ac)', 'Profit ($k/ac)', 'Total Profit ($M)', 'CapEx ($k/ac)', 'Panels/ac'};

% 3. Calculate smart custom limits for perfectly clean intervals
raw_max = max(metrics_actual, [], 1);
num_metrics = size(metrics_actual, 2);
rounded_max = zeros(1, num_metrics);

for i = 1:num_metrics
    val = raw_max(i);
    
    if val > 5000
        % Rounds up to the nearest 4000 (Tick marks every 1000)
        rounded_max(i) = ceil(val / 4000) * 4000;
    elseif val > 1000
        % Rounds up to the nearest 400 (Tick marks every 100)
        rounded_max(i) = ceil(val / 400) * 400;
    elseif val > 100
        % Rounds up to the nearest 40 (Tick marks every 10)
        rounded_max(i) = ceil(val / 40) * 40; 
    elseif val > 10
        % Rounds up to the nearest 20 (Tick marks every 5)
        rounded_max(i) = ceil(val / 20) * 20;
    else
        % Rounds up to the nearest 4 (Tick marks every 1)
        rounded_max(i) = ceil(val / 4) * 4;
    end
    
    % Failsafe to ensure no axis maximum is 0
    if rounded_max(i) == 0
        rounded_max(i) = 4;
    end
end

% --- THE MISSING LINES ---
min_vals = zeros(1, num_metrics); % Center is exactly 0
custom_limits = [min_vals; rounded_max];
% -------------------------

% 4. Plot the figure
fig5 = figure('Color','w', 'Position', [300, 300, 900, 700]);

% --- FIX FOR THE TOOLBAR WARNING ---
% This hides the UI toolbars so they don't accidentally get exported
set(fig5, 'ToolBar', 'none', 'MenuBar', 'none'); 

colors = lines(n_personas);

spider_plot(metrics_actual, ...
    'AxesLabels', metric_names_spider, ...
    'AxesLimits', custom_limits, ...       % Applies our clean rounded limits
    'AxesInterval', 4, ...                 % 4 rings
    'AxesPrecision', 0, ...                % 0 decimal places ensures clean whole numbers
    'FillOption', 'on', ...
    'FillTransparency', 0.15, ...          
    'Color', colors, ...
    'LineWidth', 2.5, ...
    'Marker', 'o', ...
    'MarkerSize', 7, ...
    'AxesFontSize', 12, ...    
    'LabelFontSize', 16);      

title('System Trade-offs: Acre Metrics vs. Total Profit', 'FontSize', 18, 'FontWeight', 'bold');
legend(personas, 'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 20);

saveas(fig5, 'graphs/poster_persona_spider_plot.png');%% ==============================================================
%% ==============================================================
%% 8. GENERATE TABLES FOR COPY/PASTING
%% ==============================================================
fprintf('\n==================================================================================================================\n');
fprintf('                                              TABLE 1: FARM METRICS                                               \n');
fprintf('==================================================================================================================\n');
fprintf('%-15s | %-15s | %-15s | %-12s | %-12s | %-15s | %-14s\n', ...
    'Persona', 'Tot Profit ($M)', 'Profit ($/ac/yr)', 'CapEx ($/ac)', 'Panels/ac', 'Crop (kg/ac/yr)', 'Emiss (t/ac)');
fprintf('------------------------------------------------------------------------------------------------------------------\n');
fprintf('%-15s | $%-14.2f | $%-14.2f | $%-11.2f | %-12.1f | %-15.2f | %-14.2f\n', ...
    'Default', total_profit(1)/1e6, profit_per_acre_yr(1), capex_per_acre(1), panels_per_acre(1), crop_per_acre_yr(1), emissions_per_acre(1)/1000);
fprintf('%-15s | $%-14.2f | $%-14.2f | $%-11.2f | %-12.1f | %-15.2f | %-14.2f\n', ...
    'Fred', total_profit(2)/1e6, profit_per_acre_yr(2), capex_per_acre(2), panels_per_acre(2), crop_per_acre_yr(2), emissions_per_acre(2)/1000);
fprintf('%-15s | $%-14.2f | $%-14.2f | $%-11.2f | %-12.1f | %-15.2f | %-14.2f\n', ...
    'Janice', total_profit(3)/1e6, profit_per_acre_yr(3), capex_per_acre(3), panels_per_acre(3), crop_per_acre_yr(3), emissions_per_acre(3)/1000);
fprintf('==================================================================================================================\n\n');

fprintf('==============================================================================================================================\n');
fprintf('                                                TABLE 2: DESIGN VARIABLES                                                     \n');
fprintf('==============================================================================================================================\n');
fprintf('%-15s | %-12s | %-12s | %-12s | %-12s | %-12s | %-12s | %-12s\n', ...
    'Persona', 'Height (m)', 'Length (m)', 'Width (m)', 'Azimuth (°)', 'Tilt (°)', 'Row Gap (m)', 'Panel Gap(m)');
fprintf('------------------------------------------------------------------------------------------------------------------------------\n');
fprintf('%-15s | %-12.2f | %-12.2f | %-12.2f | %-12.2f | %-12.2f | %-12.2f | %-12.2f\n', ...
    'Default', def_vars(1), def_vars(2), def_vars(3), rad2deg(def_vars(4)), rad2deg(def_vars(5)), def_vars(6), def_vars(7));
fprintf('%-15s | %-12.2f | %-12.2f | %-12.2f | %-12.2f | %-12.2f | %-12.2f | %-12.2f\n', ...
    'Fred', fred_vars(1), fred_vars(2), fred_vars(3), rad2deg(fred_vars(4)), rad2deg(fred_vars(5)), fred_vars(6), fred_vars(7));
fprintf('%-15s | %-12.2f | %-12.2f | %-12.2f | %-12.2f | %-12.2f | %-12.2f | %-12.2f\n', ...
    'Janice', janice_vars(1), janice_vars(2), janice_vars(3), rad2deg(janice_vars(4)), rad2deg(janice_vars(5)), janice_vars(6), janice_vars(7));
fprintf('==============================================================================================================================\n');
%% ============================================================
%% HELPER FUNCTIONS
%% ============================================================
function [winning_metrics, winning_vars] = get_utopia_data(data_struct)
    fval = data_struct.fval;
    pareto_metrics = data_struct.pareto_metrics;

    if isfield(data_struct, 'ga_solve')
        x_vars = data_struct.ga_solve; 
    else
        error('Could not find design variables matrix "ga_solve" in the loaded .mat file');
    end

    min_vals = min(fval);
    max_vals = max(fval);
    range_vals = max_vals - min_vals;
    range_vals(range_vals == 0) = 1; 

    norm_fval = (fval - min_vals) ./ range_vals;
    utopia_point = zeros(1, size(norm_fval, 2));
    distances = sqrt(sum((norm_fval - utopia_point).^2, 2));

    [~, best_idx] = min(distances);

    % Extract the exact metrics and design choices for the best compromise
    winning_metrics = pareto_metrics(best_idx, :);
    winning_vars = x_vars(best_idx, :);
end

function plot_spider(data_norm, data_actual, invert_cols, labels, group_names, colors)
    n_metrics = length(labels);
    n_groups = size(data_norm, 1);
    
    theta = linspace(0, 2*pi, n_metrics + 1);
    
    hold on;
    pax = polaraxes('Parent', gcf);
    hold(pax, 'on');
    
    % Plot the polygons
    for i = 1:n_groups
        rho = [data_norm(i, :), data_norm(i, 1)];
        polarplot(pax, theta, rho, 'LineWidth', 2.5, 'Color', colors(i,:), 'DisplayName', group_names{i});
    end
    
    % Axis formatting
    pax.ThetaTick = rad2deg(theta(1:end-1));
    pax.ThetaTickLabel = labels;
    pax.RLim = [0 1.15]; % Extended slightly to make room for outer text labels
    pax.RTick = [0.1, 1.0]; % Draw grid circles ONLY at the min and max edges
    pax.RTickLabel = {}; 
    pax.GridColor = [0.5 0.5 0.5];
    pax.GridAlpha = 0.5;
    
    pax.FontSize = 11;
    pax.FontWeight = 'bold';
    
    legend(pax, 'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 11);
    
    % --- Add unit scale values to each branch ---
    for j = 1:n_metrics
        col = data_actual(:, j);
        val_min = min(col);
        val_max = max(col);
        
        % Determine which value represents the "center" (0.1) vs the "edge" (1.0)
        if ismember(j, invert_cols)
            inner_val = val_max;
            outer_val = val_min;
        else
            inner_val = val_min;
            outer_val = val_max;
        end
        
        % Format the raw numbers for clean plotting (e.g., 15000 -> 15.0k)
        inner_str = format_num_spider(inner_val);
        outer_str = format_num_spider(outer_val);
        
        % Push text slightly left or right depending on the angle to prevent overlap
        angle_deg = rad2deg(theta(j));
        if angle_deg > 90 && angle_deg < 270
            align = 'right';
        elseif angle_deg == 90 || angle_deg == 270
            align = 'center';
        else
            align = 'left';
        end
        
        % Plot the text scales at radius 0.25 (near center) and 1.05 (outer edge)
        text(theta(j), 0.25, [' ', inner_str, ' '], 'HorizontalAlignment', align, 'FontSize', 9, 'Color', [0.4 0.4 0.4]);
        text(theta(j), 1.05, [' ', outer_str, ' '], 'HorizontalAlignment', align, 'FontSize', 10, 'Color', 'k', 'FontWeight', 'bold');
    end
end

% Local helper to format large numbers compactly for the spider plot
function str = format_num_spider(v)
    if v >= 1e6
        str = sprintf('%.2fM', v/1e6);
    elseif v >= 1000
        str = sprintf('%.1fk', v/1000);
    else
        str = sprintf('%.0f', v);
    end
end

function X_norm = normalize_spider(X, invert_cols)
    % Normalizes data for a spider plot between 0.1 and 1.0
    % invert_cols is an array of column indices where LOWER is BETTER.
    
    X_norm = NaN(size(X));
    for j = 1:size(X,2)
        col = X(:,j);
        col_min = min(col);
        col_max = max(col);
        
        if abs(col_max - col_min) < 1e-12
            X_norm(:,j) = 1;
        else
            % Scale between 0.1 (center ring) and 1.0 (outer edge)
            norm_col = 0.1 + 0.9 * ((col - col_min) ./ (col_max - col_min));
            
            % Invert if this is a "Lower is Better" metric
            if ismember(j, invert_cols)
                norm_col = 1.1 - norm_col; % Flips it so the old min (0.1) becomes max (1.0)
            end
            
            X_norm(:,j) = norm_col;
        end
    end
end

function X_norm = normalize_minmax(X)
    X_norm = NaN(size(X));
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
end

function plot_labeled_heatmap(C_norm, C_actual, row_labels, col_labels, plot_title, cmap, value_format)
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
    
    title(plot_title, 'FontSize', 13, 'FontWeight', 'bold');
    xlabel('Design Variable', 'FontWeight', 'bold');
    ylabel('Farm Persona', 'FontWeight', 'bold');
    
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
                'FontSize', 10);
        end
    end
end
%% 3D Farm Visualization Script Add-on (Corrected Cardinal Directions)
% Run this section after you have loaded def_vars, fred_vars, and janice_vars

% Visualize Default Persona (Type 1)
plot_farm_layout(def_vars, 'Default Farm Layout', 1);

% Visualize Fred (Commodity) (Type 2)
plot_farm_layout(fred_vars, 'Fred''s Optimized Farm Layout', 2);

% Visualize Janice (Specialty) (Type 3)
plot_farm_layout(janice_vars, 'Janice''s Optimized Farm Layout', 3);

%% ============================================================
% 3D Visualization Helper Function
% ============================================================
function plot_farm_layout(vars, persona_title, persona_type)
    % Extracts variables
    % [Height, Length, Width, Azimuth, Tilt, Row Gap, Panel Gap]
    h_p       = vars(1); % Height (m)
    l_p       = vars(2); % Length (m) - goes along the row
    w_p       = vars(3); % Width (m) - goes across the row
    az        = vars(4); % Azimuth (rad)
    tilt      = vars(5); % Tilt (rad)
    row_gap   = vars(6); % Clear space between rows (edge-to-edge) (m)
    panel_gap = vars(7); % Space between panels in a row (m)

    % Define layout grid size (Number of rows and panels per row to draw)
    n_rows = 3; 
    n_cols = 5; 

    % Calculate pitch (center-to-center distances)
    panel_proj_w = w_p * cos(tilt); % Projected width of panel on ground
    pitch_x = panel_proj_w + row_gap; % Distance between support posts
    pitch_y = l_p + panel_gap;

    % --- Persona-Specific Crop Logic ---
    if persona_type == 2
        % FRED (Commodity / Corn)
        % Clearance: 1m from the EDGE of the panel
        crop_width  = max(0.1, row_gap - 2.0); 
        crop_height = 3.0; 
        crop_color  = [0.90 0.80 0.20]; % Yellow
        crop_edge   = [0.70 0.60 0.10];
    else
        % JANICE (Specialty / Raspberries) & DEFAULT
        % Clearance: 0.5m from the SUPPORT POST of the panel
        crop_width  = max(0.1, pitch_x - 1.0); 
        crop_height = 2.0; 
        crop_color  = [0.40 0.70 0.30]; % Green
        crop_edge   = [0.20 0.50 0.10];
    end

    % Setup Figure
    fig = figure('Name', persona_title, 'Color', 'w', 'Position', [100, 100, 1200, 800]);
    %tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Define views based on PV_phi mapping (+X=South, -X=North, +Y=East, -Y=West)
    % view(90,0) looks from +X towards -X (Looking North)
    % view(0,0) looks from -Y towards +Y (Looking East)
    % views = {[-37.5, 30], [0, 90], [90, 0], [0, 0]};
    views = {[-15, 20]};
    % view_names = {'Isometric View', 'Top-Down View (X-Y)', 'Front View (Looking North)', 'Side View (Looking East)'};
    view_names = {''};
    
    %for v = 1:4
    for v = 1
        ax = nexttile;
        hold on;
        grid on;
        
        % 1. Draw the Ground (Soil/Dirt Color)
        pad = max(pitch_x, pitch_y);
        g_x = (n_rows * pitch_x) / 2 + pad;
        g_y = (n_cols * pitch_y) / 2 + pad;
        
        ground_x = [-g_x, g_x, g_x, -g_x];
        ground_y = [-g_y, -g_y, g_y, g_y];
        %patch('XData', ground_x, 'YData', ground_y, 'ZData', [0 0 0 0])
              %'FaceColor', [0.55 0.45 0.35], 'FaceAlpha', 1.0, 'EdgeColor', 'none'); 

        % Add Cardinal Direction Labels to the Ground (Corrected for model coordinates)
        text(-g_x - 1.5, 0, 0.1, 'S', 'FontSize', 18, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Color', 'k');
        text(g_x + 1.5, 0, 0.1, 'N', 'FontSize', 18, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Color', 'k');
        text(0, g_y + 1.5, 0.1, 'W', 'FontSize', 18, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Color', 'k');
        text(0, -g_y - 1.5, 0.1, 'E', 'FontSize', 18, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Color', 'k');

        % 2. Draw the Crop Rows 
        crop_length = (n_cols * pitch_y) + (pad * 0.5); % Extend slightly past panels
        
        for r = 0:n_rows
            % Center of the crop alley (local X)
            alley_cx = (r - n_rows/2) * pitch_x;
            alley_cy = 0; % Centered along Y
            
            draw_crop_volume(alley_cx, alley_cy, crop_width, crop_length, crop_height, az, crop_color, crop_edge);
        end

        % 3. Generate Panels and Posts
        for r = 1:n_rows
            for c = 1:n_cols
                % Center coordinates of this specific panel (local frame)
                cx = (r - (n_rows+1)/2) * pitch_x;
                cy = (c - (n_cols+1)/2) * pitch_y;
                
                % Define base panel vertices
                px = [-w_p/2, w_p/2, w_p/2, -w_p/2];
                py = [-l_p/2, -l_p/2, l_p/2, l_p/2];
                pz = [0, 0, 0, 0];
                
                % Apply Tilt (Rotate around local Y axis)
                px_tilt = px .* cos(tilt) - pz .* sin(tilt);
                pz_tilt = px .* sin(tilt) + pz .* cos(tilt);
                py_tilt = py; 
                
                % Apply Height shift
                pz_tilt = pz_tilt + h_p;
                
                % Shift to grid position
                px_shift = px_tilt + cx;
                py_shift = py_tilt + cy;
                
                % Apply Azimuth Rotation (Rotate entire field around global Z axis)
                px_final = px_shift .* cos(az) - py_shift .* sin(az);
                py_final = px_shift .* sin(az) + py_shift .* cos(az);
                pz_final = pz_tilt;
                
                % Draw Panel
                patch('XData', px_final, 'YData', py_final, 'ZData', pz_final, ...
                      'FaceColor', [0.1 0.3 0.7], 'FaceAlpha', 0.9, 'EdgeColor', [0.8 0.8 0.8], 'LineWidth', 0.5);
                  
                % Draw Support Post (From center of panel straight down to ground)
                post_x = cx * cos(az) - cy * sin(az);
                post_y = cx * sin(az) + cy * cos(az);
                plot3([post_x, post_x], [post_y, post_y], [0, h_p], ...
                      'Color', [0.3 0.3 0.3], 'LineWidth', 3.0);
            end
        end
        
        % Formatting and Camera Angles
        view(ax, views{v});
        title(view_names{v}, 'FontSize', 1, 'FontWeight', 'bold');
        xlabel('North-South (m)', 'FontSize', 24);
        ylabel('East–West (m)', 'FontSize', 24);
        zlabel('Height (m)', 'FontSize', 24);
        
        % Force 1:1 accurate scaling
        axis(ax, 'equal');
        
        % Give extra room so labels aren't cropped out
        xlim([-g_x - 3, g_x + 3]);
        ylim([-g_y - 3, g_y + 3]);
        
        max_z_needed = max([h_p + (w_p/2 * sin(tilt)) + 1, crop_height + 0.5]);
        zlim([0, max_z_needed]); 
        
        % Add nice lighting for the 3D views
        if v == 1
            camlight('headlight');
            lighting gouraud;
            material dull;
        end
    end
    
    sgtitle(persona_title, 'FontSize', 27, 'FontWeight', 'bold');
end

%% Helper Function to Draw 3D Crop Volumes
function draw_crop_volume(cx, cy, w, l, h, az, face_color, edge_color)
    % Defines a 3D box for the crop row and rotates it by the farm azimuth
    
    % Vertices of a box centered at origin
    x = [-w/2, w/2, w/2, -w/2, -w/2, w/2, w/2, -w/2];
    y = [-l/2, -l/2, l/2, l/2, -l/2, -l/2, l/2, l/2];
    z = [0, 0, 0, 0, h, h, h, h];
    
    % Shift to local center
    x = x + cx;
    y = y + cy;
    
    % Rotate vertices around global Z axis by azimuth
    x_rot = x .* cos(az) - y .* sin(az);
    y_rot = x .* sin(az) + y .* cos(az);
    
    % Define the 5 visible faces (skip the bottom face hidden underground)
    faces = [1 2 6 5;  % Front
             2 3 7 6;  % Right
             3 4 8 7;  % Back
             4 1 5 8;  % Left
             5 6 7 8]; % Top
         
    % Draw the crops as lush volumes
    patch('Vertices', [x_rot', y_rot', z'], 'Faces', faces, ...
          'FaceColor', face_color, 'FaceAlpha', 0.75, ...
          'EdgeColor', edge_color, 'LineWidth', 0.5);
end