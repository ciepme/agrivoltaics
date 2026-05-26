%% run_exact_design_case.m
% Runs agrivoltaic model for one exact design vector x.
% Useful for rerunning a GA/MOGA result with one variable changed.

clear;
clc;
close all;

addpath(genpath(pwd));

%% Load model definitions first
% Important: your definition file contains "clear",
% so do not define variables before this line.
agrivoltaics_variable_definition;

%% User setup

% Load previous compiled results
master = load("results/master_best_designs.mat");
best_table = master.best_table;

% Choose the design to rerun
system_to_use = "Single-Axis";      % "Fixed" or "Single-Axis"
objective_to_use = "POWER"; % "PROFIT", "EMISSIONS", "POWER", "CROP", "MOGA-Balanced"

row = find(best_table.System == system_to_use & ...
           best_table.Objective == objective_to_use, 1);

if isempty(row)
    error("Could not find %s / %s in best_table.", system_to_use, objective_to_use);
end

x = best_table.X{row};

fprintf("Loaded design:\n");
fprintf("  System: %s\n", best_table.System(row));
fprintf("  Objective: %s\n", best_table.Objective(row));
fprintf("  Source: %s\n", best_table.SourceFile(row));
fprintf("  Number of variables: %d\n", numel(x));

%% Modify design variables here

% Variable order:
% 1 = Height
% 2 = Panel length
% 3 = Panel width
% 4 = Azimuth, rad
% 5 = Fixed tilt, rad
% 6 = Row spacing
% 7 = Panel gap
%
% For single-axis:
% 8:103 = tracking angles, rad

x_modified = x;

% Example changes:
x_modified(1) = 3.73;          % height, m
x_modified(6) = .1;           % row spacing, m
x_modified(7) = 0.1;          % panel gap, m

% Example: change one single-axis tracking angle
% Winter hour 9
if numel(x_modified) == 103
    season_num = 4; % 1 spring, 2 summer, 3 fall, 4 winter
    hour_num = 9;

    idx = 7 + (season_num-1)*24 + hour_num;

    x_modified(idx) = deg2rad(-35);
end

%% Enforce bounds

x_modified = max(x_modified(:).', lb(:).');
x_modified = min(x_modified(:).', ub(:).');

% %% Optional: repair single-axis slew constraints
% 
% if agriParams.tracking_mode == 1 && numel(x_modified) == 103
% 
%     max_slew = agriParams.max_slew_per_hour;
% 
%     for s = 1:4
%         start_idx = 7 + (s-1)*24 + 1;
%         end_idx   = 7 + s*24;
% 
%         curve = x_modified(start_idx:end_idx);
% 
%         curve = enforce_slew_curve( ...
%             curve, ...
%             lb(start_idx:end_idx), ...
%             ub(start_idx:end_idx), ...
%             max_slew);
% 
%         x_modified(start_idx:end_idx) = curve;
%     end
% end

%% Run original and modified designs

metrics_original = agrivoltaic_wrapper(x, agriParams);
metrics_modified = agrivoltaic_wrapper(x_modified, agriParams);

%% Display results

metric_names = ["Emissions", "Profit", "SocialCost", "PVRevenue", ...
                "CropRevenue", "CropYield", "Panels", "Energy"];

T = table( ...
    metric_names(:), ...
    metrics_original(:), ...
    metrics_modified(:), ...
    metrics_modified(:) - metrics_original(:), ...
    'VariableNames', {'Metric', 'Original', 'Modified', 'Difference'});

disp(T);

%% Save result

if ~exist("results", "dir")
    mkdir("results");
end

save("results/exact_design_case_" + system_to_use + "_" + objective_to_use + ".mat", ...
    "x", ...
    "x_modified", ...
    "metrics_original", ...
    "metrics_modified", ...
    "T", ...
    "system_to_use", ...
    "objective_to_use");

fprintf("\nSaved exact design case to results/exact_design_case_%s_%s.mat\n", ...
    system_to_use, objective_to_use);

%% Plot tracking curves if single-axis

if numel(x_modified) == 103

    plot_tracking_comparison(x, x_modified, "Original", "Modified");

    if ~exist("graphs", "dir")
        mkdir("graphs");
    end

    saveas(gcf, "graphs/exact_design_tracking_comparison_" + system_to_use + "_" + objective_to_use + ".png");
end

%% Helper functions

function angles_out = enforce_slew_curve(angles_in, lb_curve, ub_curve, max_slew)

    angles_out = angles_in;

    n_iter = 10;

    for iter = 1:n_iter

        angles_out = max(angles_out, lb_curve);
        angles_out = min(angles_out, ub_curve);

        for h = 2:length(angles_out)

            lower_from_prev = angles_out(h-1) - max_slew;
            upper_from_prev = angles_out(h-1) + max_slew;

            angles_out(h) = min(max(angles_out(h), lower_from_prev), upper_from_prev);

            angles_out(h) = max(angles_out(h), lb_curve(h));
            angles_out(h) = min(angles_out(h), ub_curve(h));

        end

        for h = length(angles_out)-1:-1:1

            lower_from_next = angles_out(h+1) - max_slew;
            upper_from_next = angles_out(h+1) + max_slew;

            angles_out(h) = min(max(angles_out(h), lower_from_next), upper_from_next);

            angles_out(h) = max(angles_out(h), lb_curve(h));
            angles_out(h) = min(angles_out(h), ub_curve(h));

        end
    end
end

function plot_tracking_comparison(x1, x2, label1, label2)

    seasons = {'Spring', 'Summer', 'Fall', 'Winter'};
    hours = 1:24;

    figure('Color','w', 'Position', [100 100 1100 850]);
    tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

    for s = 1:4

        nexttile;
        hold on;
        grid on;

        start_idx = 7 + (s-1)*24 + 1;
        end_idx   = 7 + s*24;

        plot(hours, rad2deg(x1(start_idx:end_idx)), '-o', ...
            'LineWidth', 2, ...
            'DisplayName', label1);

        plot(hours, rad2deg(x2(start_idx:end_idx)), '-s', ...
            'LineWidth', 2, ...
            'DisplayName', label2);

        title(seasons{s});
        xlabel('Hour');
        ylabel('Tracking Angle (deg)');
        xlim([1 24]);
        xticks(1:2:24);

        if s == 1
            legend('Location','best');
        end
    end

    sgtitle('Single-Axis Tracking Comparison');
end