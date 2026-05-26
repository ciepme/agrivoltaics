%% compile_best_agrivoltaic_results.m
% This script compiles all saved single-objective GA and MOGA results,
% selects the best designs, and saves a curated master .mat file.
%
% It supports:
%   1. Existing single-objective GA files:
%      agrivoltaic_comparative_optimization_data*.mat
%
%   2. Existing MOGA files:
%      MOGA_data_*.mat
%
% For MOGA, it selects a balanced Pareto design using normalized distance
% to the utopia point for:
%   [Emissions Reduction, Profit, Crop/Biomass]

clear;
clc;
close all;

addpath(genpath(pwd));

%% Load model definitions once
% Your definition file currently clears the workspace, so call it before
% defining compiler variables.
agrivoltaics_variable_definition;

% Store parameter versions for fixed and single-axis evaluation.
params_template = agriParams;

params_fixed = params_template;
params_fixed.tracking_mode = 0;

params_single = params_template;
params_single.tracking_mode = 1;

%% Compiler setup

results_dir = "results";
root_dir = ".";

if ~exist(results_dir, "dir")
    mkdir(results_dir);
end

compiled_output_file = fullfile(results_dir, "master_best_designs.mat");
compiled_summary_file = fullfile(results_dir, "master_best_designs_summary.csv");
all_candidates_file = fullfile(results_dir, "all_candidate_designs_summary.csv");

%% Metric column definitions
% agrivoltaic_wrapper output:
% [Emissions, Profit, SocialCost, PVRev, CropRev, Biomass/Crop, Panels, Energy]

col_emissions = 1;
col_profit    = 2;
col_social    = 3;
col_pvrev     = 4;
col_croprev   = 5;
col_crop      = 6;
col_panels    = 7;
col_energy    = 8;

main_objectives = ["PROFIT", "EMISSIONS", "POWER", "CROP"];

%% Candidate storage

CandidateID = strings(0,1);
SourceFile = strings(0,1);
System = strings(0,1);
Method = strings(0,1);
Objective = strings(0,1);

Profit = [];
Energy = [];
Emissions = [];
CropYield = [];
Panels = [];
PVRevenue = [];
SocialCost = [];
CropRevenue = [];

TrackingMode = [];
NumVars = [];
MOGABalanceDistance = [];

Metrics = {};
X = {};

candidate_counter = 0;

%% ============================================================
% 1. Load single-objective GA files
% ============================================================

single_files = collect_mat_files([results_dir, root_dir], ...
    "*agrivoltaic_comparative_optimization_data*.mat");

fprintf("\nFound %d single-objective GA result files.\n", length(single_files));

for f = 1:length(single_files)

    file_path = single_files(f).fullpath;
    [~, file_name, ext] = fileparts(file_path);
    display_name = string(file_name + ext);

    S = load(file_path);

    if ~isfield(S, "targets_to_run") || ...
       ~isfield(S, "results_matrix") || ...
       ~isfield(S, "x_best_set")

        warning("Skipping %s because it does not contain targets_to_run, results_matrix, and x_best_set.", display_name);
        continue;
    end

    targets = string(S.targets_to_run);
    results_matrix = S.results_matrix;
    x_best_set = S.x_best_set;

    if size(results_matrix, 1) ~= length(targets)
        warning("Skipping %s because results_matrix rows do not match targets_to_run.", display_name);
        continue;
    end

    for i = 1:length(targets)

        x_i = x_best_set(i, :);
        metrics_i = results_matrix(i, :);

        mode_i = infer_tracking_mode_from_x(x_i);
        system_i = mode_to_system_name(mode_i);

        candidate_counter = candidate_counter + 1;

        CandidateID(end+1,1) = "SO_" + system_i + "_" + targets(i) + "_" + string(candidate_counter);
        SourceFile(end+1,1) = display_name;
        System(end+1,1) = system_i;
        Method(end+1,1) = "SingleObjectiveGA";
        Objective(end+1,1) = targets(i);

        Profit(end+1,1) = safe_metric(metrics_i, col_profit);
        Energy(end+1,1) = safe_metric(metrics_i, col_energy);
        Emissions(end+1,1) = safe_metric(metrics_i, col_emissions);
        CropYield(end+1,1) = safe_metric(metrics_i, col_crop);
        Panels(end+1,1) = safe_metric(metrics_i, col_panels);
        PVRevenue(end+1,1) = safe_metric(metrics_i, col_pvrev);
        SocialCost(end+1,1) = safe_metric(metrics_i, col_social);
        CropRevenue(end+1,1) = safe_metric(metrics_i, col_croprev);

        TrackingMode(end+1,1) = mode_i;
        NumVars(end+1,1) = length(x_i);
        MOGABalanceDistance(end+1,1) = NaN;

        Metrics{end+1,1} = metrics_i;
        X{end+1,1} = x_i;

    end
end

%% ============================================================
% 2. Load MOGA files
% ============================================================

moga_files = collect_mat_files([results_dir, root_dir], "MOGA_data_*.mat");

fprintf("\nFound %d MOGA result files.\n", length(moga_files));

for f = 1:length(moga_files)

    file_path = moga_files(f).fullpath;
    [~, file_name, ext] = fileparts(file_path);
    display_name = string(file_name + ext);

    S = load(file_path);

    if ~isfield(S, "ga_solve") || ~isfield(S, "fval")
        warning("Skipping %s because it does not contain ga_solve and fval.", display_name);
        continue;
    end

    ga_solve = S.ga_solve;
    fval = S.fval;

    if isempty(ga_solve) || isempty(fval)
        warning("Skipping %s because ga_solve or fval is empty.", display_name);
        continue;
    end

    mode_i = infer_tracking_mode_from_x(ga_solve(1, :));
    system_i = mode_to_system_name(mode_i);

    % Your MOGA objective is:
    % fitness = [-emissions, -profit, -yearly_biomass]
    %
    % Therefore benefits are:
    % benefits = -fval
    %
    % benefits columns:
    %   1 = Emissions Reduction
    %   2 = Profit
    %   3 = Crop/Biomass
    benefits = -fval(:, 1:3);

    [idx_balanced, balance_distance] = choose_balanced_moga_design(benefits);

    x_balanced = ga_solve(idx_balanced, :);

    % Evaluate full 8 metrics if possible.
    % This lets the MOGA row appear in all your summary plots/tables.
    if isfield(S, "pareto_metrics")
        metrics_balanced = S.pareto_metrics(idx_balanced, :);
    else
        metrics_balanced = evaluate_design_metrics_safely( ...
            x_balanced, ...
            mode_i, ...
            params_fixed, ...
            params_single, ...
            benefits(idx_balanced, :));
    end

    candidate_counter = candidate_counter + 1;

    CandidateID(end+1,1) = "MOGA_" + system_i + "_Balanced_" + string(candidate_counter);
    SourceFile(end+1,1) = display_name;
    System(end+1,1) = system_i;
    Method(end+1,1) = "MOGA";
    Objective(end+1,1) = "MOGA-Balanced";

    Profit(end+1,1) = safe_metric(metrics_balanced, col_profit);
    Energy(end+1,1) = safe_metric(metrics_balanced, col_energy);
    Emissions(end+1,1) = safe_metric(metrics_balanced, col_emissions);
    CropYield(end+1,1) = safe_metric(metrics_balanced, col_crop);
    Panels(end+1,1) = safe_metric(metrics_balanced, col_panels);
    PVRevenue(end+1,1) = safe_metric(metrics_balanced, col_pvrev);
    SocialCost(end+1,1) = safe_metric(metrics_balanced, col_social);
    CropRevenue(end+1,1) = safe_metric(metrics_balanced, col_croprev);

    TrackingMode(end+1,1) = mode_i;
    NumVars(end+1,1) = length(x_balanced);
    MOGABalanceDistance(end+1,1) = balance_distance;

    Metrics{end+1,1} = metrics_balanced;
    X{end+1,1} = x_balanced;

    fprintf("MOGA balanced design selected from %s: Pareto row %d, distance %.4f\n", ...
        display_name, idx_balanced, balance_distance);

end

%% Build candidate table

candidate_table = table( ...
    CandidateID, ...
    SourceFile, ...
    System, ...
    Method, ...
    Objective, ...
    Profit, ...
    Energy, ...
    Emissions, ...
    CropYield, ...
    Panels, ...
    PVRevenue, ...
    SocialCost, ...
    CropRevenue, ...
    TrackingMode, ...
    NumVars, ...
    MOGABalanceDistance, ...
    Metrics, ...
    X);

fprintf("\nAll candidates loaded:\n");
disp(candidate_table(:, {'CandidateID','SourceFile','System','Method','Objective','Profit','Energy','Emissions','CropYield','MOGABalanceDistance'}));

%% ============================================================
% 3. Select best designs across all runs
% ============================================================

best_row_indices = [];

systems_to_use = ["Fixed", "Single-Axis"];

for s = 1:length(systems_to_use)

    system_name = systems_to_use(s);

    %% Single-objective winners
    for o = 1:length(main_objectives)

        obj = main_objectives(o);

        rows = find(candidate_table.System == system_name & ...
                    candidate_table.Method == "SingleObjectiveGA" & ...
                    candidate_table.Objective == obj);

        if isempty(rows)
            warning("No single-objective %s result found for %s.", obj, system_name);
            continue;
        end

        switch obj
            case "PROFIT"
                values = candidate_table.Profit(rows);
            case "EMISSIONS"
                values = candidate_table.Emissions(rows);
            case "POWER"
                values = candidate_table.Energy(rows);
            case "CROP"
                values = candidate_table.CropYield(rows);
            otherwise
                error("Unknown objective: %s", obj);
        end

        [~, local_best] = max(values);
        best_row_indices(end+1,1) = rows(local_best);

    end

    %% MOGA balanced winner
    moga_rows = find(candidate_table.System == system_name & ...
                     candidate_table.Method == "MOGA" & ...
                     candidate_table.Objective == "MOGA-Balanced");

    if isempty(moga_rows)
        warning("No MOGA-Balanced result found for %s.", system_name);
    elseif length(moga_rows) == 1
        best_row_indices(end+1,1) = moga_rows;
    else
        % If multiple MOGA runs exist for this system, choose the one
        % closest to the utopia point using full metrics.
        benefit_matrix = [ ...
            candidate_table.Emissions(moga_rows), ...
            candidate_table.Profit(moga_rows), ...
            candidate_table.CropYield(moga_rows)];

        benefit_norm = normalize_minmax_local(benefit_matrix);
        dist_to_utopia = vecnorm(benefit_norm - ones(size(benefit_norm)), 2, 2);

        [~, local_best] = min(dist_to_utopia);
        best_row_indices(end+1,1) = moga_rows(local_best);
    end

end

best_table = candidate_table(best_row_indices, :);

fprintf("\nSelected best designs:\n");
disp(best_table(:, {'CandidateID','SourceFile','System','Method','Objective','Profit','Energy','Emissions','CropYield','MOGABalanceDistance'}));

%% ============================================================
% 4. Also create backward-compatible matrices for plotting convenience
% ============================================================

plot_objectives = ["PROFIT", "EMISSIONS", "POWER", "CROP", "MOGA-Balanced"];

fixed_results_matrix = NaN(length(plot_objectives), 8);
single_results_matrix = NaN(length(plot_objectives), 8);

fixed_x_best_cell = cell(length(plot_objectives), 1);
single_x_best_cell = cell(length(plot_objectives), 1);

for i = 1:length(plot_objectives)

    obj = plot_objectives(i);

    fixed_idx = find(best_table.System == "Fixed" & best_table.Objective == obj, 1);
    single_idx = find(best_table.System == "Single-Axis" & best_table.Objective == obj, 1);

    if ~isempty(fixed_idx)
        fixed_results_matrix(i, :) = best_table.Metrics{fixed_idx};
        fixed_x_best_cell{i} = best_table.X{fixed_idx};
    end

    if ~isempty(single_idx)
        single_results_matrix(i, :) = best_table.Metrics{single_idx};
        single_x_best_cell{i} = best_table.X{single_idx};
    end

end

%% Save outputs

save(compiled_output_file, ...
    "candidate_table", ...
    "best_table", ...
    "plot_objectives", ...
    "fixed_results_matrix", ...
    "single_results_matrix", ...
    "fixed_x_best_cell", ...
    "single_x_best_cell");

summary_table = removevars(best_table, {'Metrics','X'});
candidate_summary_table = removevars(candidate_table, {'Metrics','X'});

writetable(summary_table, compiled_summary_file);
writetable(candidate_summary_table, all_candidates_file);

fprintf("\nSaved compiled master file:\n%s\n", compiled_output_file);
fprintf("Saved summary CSV:\n%s\n", compiled_summary_file);
fprintf("Saved all-candidates CSV:\n%s\n", all_candidates_file);

%% ============================================================
% Helper functions
% ============================================================

function files_out = collect_mat_files(search_dirs, pattern)

    files_out = struct('name', {}, 'folder', {}, 'fullpath', {});

    seen_paths = strings(0,1);

    for d = 1:length(search_dirs)

        this_dir = string(search_dirs(d));

        if ~exist(this_dir, "dir")
            continue;
        end

        files = dir(fullfile(this_dir, pattern));

        for i = 1:length(files)

            fullpath = string(fullfile(files(i).folder, files(i).name));

            % Avoid duplicates if results_dir and root_dir overlap
            if any(seen_paths == fullpath)
                continue;
            end

            % Avoid loading the compiled master file accidentally
            if contains(files(i).name, "master_best_designs")
                continue;
            end

            seen_paths(end+1,1) = fullpath;

            files_out(end+1).name = files(i).name;
            files_out(end).folder = files(i).folder;
            files_out(end).fullpath = fullpath;

        end
    end
end

function mode = infer_tracking_mode_from_x(x)

    n = length(x);

    if n == 7
        mode = 0;
    elseif n == 103
        mode = 1;
    else
        warning("Unknown design vector length %d. Assuming fixed-axis.", n);
        mode = 0;
    end
end

function system_name = mode_to_system_name(mode)

    if mode == 0
        system_name = "Fixed";
    elseif mode == 1
        system_name = "Single-Axis";
    else
        system_name = "Unknown";
    end
end

function val = safe_metric(metrics, idx)

    if isempty(metrics) || length(metrics) < idx
        val = NaN;
    else
        val = metrics(idx);
    end
end

function [idx_best, distance_best] = choose_balanced_moga_design(benefits)

    % benefits should be:
    %   column 1 = Emissions Reduction
    %   column 2 = Profit
    %   column 3 = Crop/Biomass
    %
    % Higher is better for all columns.

    benefits_norm = normalize_minmax_local(benefits);

    utopia = ones(1, size(benefits_norm, 2));

    distances = vecnorm(benefits_norm - utopia, 2, 2);

    [distance_best, idx_best] = min(distances);

end

function metrics = evaluate_design_metrics_safely(x, mode, params_fixed, params_single, fallback_benefits)

    % Returns full 8-column metrics if agrivoltaic_wrapper succeeds.
    % If not, returns partial metrics using MOGA fval-derived values.

    if mode == 0
        params = params_fixed;
    else
        params = params_single;
    end

    try
        metrics = agrivoltaic_wrapper(x, params);
    catch ME

        warning("Could not evaluate agrivoltaic_wrapper for a MOGA design. Using partial metrics from fval. Error was: %s", ME.message);

        metrics = NaN(1, 8);

        % fallback_benefits:
        % [emissions, profit, crop]
        metrics(1) = fallback_benefits(1);
        metrics(2) = fallback_benefits(2);
        metrics(6) = fallback_benefits(3);

    end
end

function X_norm = normalize_minmax_local(X)

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