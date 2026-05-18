%% Clear and Setup
clear;
clc;
close all;
addpath(genpath(pwd));
agrivoltaics_variable_definition;

% =========================================================================
% USER SETUP
% Choose what to optimize. You can put one or multiple in this cell array.
% Options: 'PROFIT', 'EMISSIONS', 'POWER', 'CROP', 'ALL'
% =========================================================================
targets_to_run = {'ALL'}; 

if ismember('ALL', targets_to_run)
    targets_to_run = {'PROFIT', 'EMISSIONS', 'POWER', 'CROP'};
end

pop_size = 8;
max_gen = 8; % Adjust based on your time constraints
num_vars = length(lb);

%% 1. Generate the Population
if agriParams.tracking_mode == 1
    % Calculate the true sun tracking curve
    agriVar.tracking_angles = generate_physics_tracking(agriParams, agriVar);
    
    % Clamp them to your max tilt inline 
    agriVar.tracking_angles = max(agriVar.tracking_angles, -agriParams.PV_max_tilt);
    agriVar.tracking_angles = min(agriVar.tracking_angles,  agriParams.PV_max_tilt);
end

% build x0
x0 = agriVarStruct2Array(agriVar, agriParams);

% force x0 to obey bounds
x0 = max(x0(:).', lb(:).');
x0 = min(x0(:).', ub(:).');

% Build the population array
pop = zeros(pop_size, num_vars);

% Member 1 of GA population = true physics-based smart guess perfectly constrained
pop(1,:) = x0;

% All the other members of initial population
for i = 2:pop_size
    candidate = x0;
    
    if agriParams.tracking_mode == 1
        idx = 8:num_vars;
        % Pure random tracking generation
        span = ub(idx) - lb(idx);
        random_angles = lb(idx) + rand(1, length(idx)) .* span;
        
        % Smooth to respect slew limit
        smoothed_angles = smoothdata(random_angles, 'gaussian', 5);
        
        % Clamp (forces night hours to 0)
        candidate(idx) = max(smoothed_angles, lb(idx));
        candidate(idx) = min(candidate(idx), ub(idx));
    else
        % Fixed-axis: Pure Random layout generation
        span = ub(1:num_vars) - lb(1:num_vars);
        candidate(1:num_vars) = lb(1:num_vars) + rand(1, num_vars) .* span;
        
        candidate(1:num_vars) = max(candidate(1:num_vars), lb(1:num_vars));
        candidate(1:num_vars) = min(candidate(1:num_vars), ub(1:num_vars));
    end
    
    pop(i,:) = candidate;
end

%% 2. Set Up GA Options and Constraints
options = optimoptions('ga', 'PopulationSize', pop_size, 'MaxGenerations', max_gen, ...
    'FunctionTolerance', 1e-4, 'Display', 'iter', ...
    'InitialPopulationMatrix', pop, 'UseParallel', true);

A = []; B = []; Aeq = []; Beq = [];
if agriParams.tracking_mode == 1
    max_slew_per_hour = deg2rad(45);
    num_steps = 23; 
    total_constraints = 4 * (num_steps * 2);
    
    A = zeros(total_constraints, num_vars);
    B = ones(total_constraints, 1) * max_slew_per_hour;
    
    row = 1;
    for s = 1:4
        offset = 7 + (s-1)*24; 
        for h = 1:23
            v1 = offset + h; v2 = offset + h + 1;
            A(row, v1) = -1; A(row, v2) = 1;  row = row + 1;
            A(row, v1) = 1;  A(row, v2) = -1; row = row + 1;
        end
    end
end

%% 3. Execute Optimization Loop
num_targets = length(targets_to_run);
results_matrix = zeros(num_targets, 8); % To store the 8 wrapper outputs
x_best_set = zeros(num_targets, num_vars);

for i = 1:num_targets
    current_target = targets_to_run{i};

    fprintf('Starting Maximization for %s\n', current_target);
    
    tic;
    [x_best, fval, exitflag, output] = ga(@(x) targeted_objective_wrapper(x, agriParams, current_target), ...
        num_vars, A, B, Aeq, Beq, lb, ub, [], options);
    time_taken = toc;
    
    % Evaluate the winning layout to get all its metrics
    winning_metrics = agrivoltaic_wrapper(x_best, agriParams);
    
    % Store for plotting
    results_matrix(i, :) = winning_metrics;
    x_best_set(i, :) = x_best;
    
    fprintf('\n--- RESULTS FOR %s ---\n', current_target);
    fprintf('Time Taken: %.2f seconds\n', time_taken);
    fprintf('Profit: $%.2f M\n', winning_metrics(2) / 1e6);
    fprintf('Emissions Reduction: %.2f kt CO2e\n', winning_metrics(1) / 1e6);
    fprintf('Crop Yield: %.2f g/m^2\n', winning_metrics(6));
    fprintf('Yearly Energy: $%.2f\n', winning_metrics(4));
end

%% 4. Plot Comparative Graphic
if num_targets > 1
    fig = figure('Name', 'Objective Comparison', 'Color', 'w', 'Position', [100, 100, 1000, 800]);
    t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(t, 'Comparison of Winning Layouts by Optimization Target', 'FontWeight', 'bold', 'FontSize', 14);
    
    target_labels = categorical(targets_to_run);
    
    % Plot 1: Profit
    nexttile;
    bar(target_labels, results_matrix(:, 2) / 1e6, 'FaceColor', [0.2 0.6 0.5]);
    title('Total Profit'); ylabel('Millions ($M)'); grid on;
    
    % Plot 2: Emissions / Power
    nexttile;
    bar(target_labels, results_matrix(:, 1) / 1e6, 'FaceColor', [0.3 0.4 0.7]);
    title('Emissions Reduction (Power)'); ylabel('kt CO2e'); grid on;
    
    % Plot 3: Crop Yield
    nexttile;
    bar(target_labels, results_matrix(:, 6), 'FaceColor', [0.8 0.4 0.2]);
    title('Crop Yield'); ylabel('kg/year'); grid on;
    
  % Plot 4: PV Revenue
    nexttile;
    bar(target_labels, results_matrix(:, 4) / 1e6, 'FaceColor', [0.9 0.7 0.1]);
    title('Yearly Energy'); ylabel('kWh/year'); grid on;
    
    saveas(fig, 'graphs/single_objective_comparisons.png');
    fprintf('\nSaved comparison chart to graphs/single_objective_comparisons.png\n');
end

% Save workspace data
save('agrivoltaic_comparative_optimization_data.mat', 'targets_to_run', 'results_matrix', 'x_best_set');

%% 
% Plot physical design comparison

fig_layout = figure('Name', 'Optimal Layout Comparison', 'Color', 'w', 'Position', [150, 150, 1200, 600]);
t_layout = tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
title(t_layout, 'Optimal Physical Design Variables by Objective', 'FontWeight', 'bold', 'FontSize', 14);

var_names = {'Height (m)', 'Length (m)', 'Width (m)', 'Azimuth (deg)', 'Tilt (deg)', 'Row Spacing (m)', 'Panel Gap (m)'};

% Convert Azimuth (4) and Tilt (5) from radians to degrees for readability
layout_data = x_best_set(:, 1:7);
layout_data(:, 4) = rad2deg(layout_data(:, 4));
layout_data(:, 5) = rad2deg(layout_data(:, 5));

colors = lines(num_targets); % Get distinct colors for each target

for v = 1:7
    nexttile;
    hold on; grid on;
    for i = 1:num_targets
        bar(categorical({targets_to_run{i}}), layout_data(i, v), 'FaceColor', colors(i,:));
    end
    title(var_names{v}, 'FontWeight', 'bold');
    ylabel('Value');
end

saveas(fig_layout, 'graphs/optimal_layout_comparison.png');
fprintf('Saved layout comparison chart to graphs/optimal_layout_comparison.png\n');

%%
% Plot tracking curves for single-axis mode
if agriParams.tracking_mode == 1
    seasons = {'Spring', 'Summer', 'Fall', 'Winter'};
    fig_track = figure('Name', 'Optimal Tracking Curves', 'Color', 'w', 'Position', [200, 200, 1000, 800]);
    t_track = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(t_track, 'Optimized Single-Axis Tracking Curves by Objective', 'FontWeight', 'bold', 'FontSize', 14);
    
    hours = 1:24;
    max_tilt_deg = rad2deg(agriParams.PV_max_tilt);
    
    for s = 1:4
        nexttile;
        hold on; grid on;
        title(seasons{s}, 'FontWeight', 'bold');
        xlabel('Hour of Day');
        ylabel('Panel Angle (Degrees)');
        
        % Map the flat array back to the specific season's 24 hours
        start_idx = 7 + (s-1)*24 + 1;
        end_idx = 7 + s*24;
        
        % Plot the curve for each objective
        for i = 1:num_targets
            angles_rad = x_best_set(i, start_idx:end_idx);
            angles_deg = rad2deg(angles_rad);
            
            plot(hours, angles_deg, '-o', 'LineWidth', 2, 'MarkerSize', 4, ...
                 'Color', colors(i,:), 'DisplayName', targets_to_run{i});
        end
        
        % Plot the "Ideal Physics" baseline for comparison
        ideal_angles_deg = rad2deg(agriVar.tracking_angles(s, :));
        plot(hours, ideal_angles_deg, '--k', 'LineWidth', 1.5, 'DisplayName', 'Pure Sun Tracking');
        
        % Formatting
        xlim([1 24]);
        ylim([-max_tilt_deg - 5, max_tilt_deg + 5]);
        xticks(1:2:24);
        
        if s == 1
            legend('Location', 'best');
        end
    end
    saveas(fig_track, 'graphs/optimal_tracking_curves.png');
    fprintf('Saved tracking curves chart to graphs/optimal_tracking_curves.png\n');
end

%% Helper function
function fitness = targeted_objective_wrapper(x, params, target)
    % Run the standard wrapper to get all outputs
    % raw = [Emissions(1), Profit(2), SocialCost(3), PVRev(4), CropRev(5), Biomass(6), Panels(7), Energy(8)]
    raw = agrivoltaic_wrapper(x, params);
    
    % GA strictly MINIMIZES. So to maximize something, we return its negative value.
    switch target
        case 'PROFIT'
            fitness = -raw(2); % in USD over 30 years
        case 'EMISSIONS'
            fitness = -raw(1); %in CO2 / year
        case 'POWER'
            fitness = -raw(8); % in kWh /year
        case 'CROP'
            fitness = -raw(6); % in kg /year
        otherwise
            error('Unknown target. Please use PROFIT, EMISSIONS, POWER, or CROP.');
    end
end