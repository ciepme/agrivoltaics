%% Parallel validation and timing benchmark
clear;
clc;
close all;
bdclose('all');

addpath(genpath(pwd));
agrivoltaics_variable_definition;

%% User-adjustable validation budget
RUN_WRAPPER_SMOKE = true;
RUN_GA_SHORT = true;
RUN_MOGA_SHORT = true;
RUN_SQP_FIXED_SMOKE = true;
RUN_SQP_SINGLE_AXIS_SMOKE = false;

GA_POP_SIZE = 20;
GA_MAX_GEN = 3;
MOGA_POP_SIZE = 24;
MOGA_MAX_GEN = 3;
SQP_MAX_ITER = 1;
SQP_MAX_FUN_EVALS = 40;
RNG_SEED = 11;

timestamp = datestr(now, 'yyyymmdd_HHMMSS');
results_dir = fullfile(pwd, 'results', ['parallel_validation_' timestamp]);
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

templateParams = agriParams;
case_defs = [
    struct('name', "legacy_fixed",       'geometry', "legacy",       'tracking', 0)
    struct('name', "legacy_single",      'geometry', "legacy",       'tracking', 1)
    struct('name', "row_centered_fixed", 'geometry', "row_centered", 'tracking', 0)
    struct('name', "row_centered_single",'geometry', "row_centered", 'tracking', 1)
];

rows = struct([]);
raw_runs = struct();

for c = 1:numel(case_defs)
    case_def = case_defs(c);
    fprintf('\n=== Validating %s ===\n', case_def.name);

    [params, vars, lb_case, ub_case] = configure_agri_mode(templateParams, case_def.tracking, case_def.geometry);
    x0 = agriVarStruct2Array(vars, params);

    if RUN_WRAPPER_SMOKE
        serial = run_wrapper_case(x0, params, false);
        parallel = run_wrapper_case(x0, params, true);
        rows = append_summary_row(rows, case_def, "wrapper", serial, parallel, params, x0);
        raw_runs.(case_def.name).wrapper_serial = serial;
        raw_runs.(case_def.name).wrapper_parallel = parallel;
    end

    if RUN_GA_SHORT
        serial = run_ga_case(vars, params, lb_case, ub_case, GA_POP_SIZE, GA_MAX_GEN, RNG_SEED, false);
        parallel = run_ga_case(vars, params, lb_case, ub_case, GA_POP_SIZE, GA_MAX_GEN, RNG_SEED, true);
        rows = append_summary_row(rows, case_def, "ga_short", serial, parallel, params, parallel.x);
        raw_runs.(case_def.name).ga_serial = serial;
        raw_runs.(case_def.name).ga_parallel = parallel;
    end

    if RUN_MOGA_SHORT
        serial = run_moga_case(vars, params, lb_case, ub_case, MOGA_POP_SIZE, MOGA_MAX_GEN, RNG_SEED, false);
        parallel = run_moga_case(vars, params, lb_case, ub_case, MOGA_POP_SIZE, MOGA_MAX_GEN, RNG_SEED, true);
        rows = append_summary_row(rows, case_def, "moga_short", serial, parallel, params, parallel.x);
        raw_runs.(case_def.name).moga_serial = serial;
        raw_runs.(case_def.name).moga_parallel = parallel;
    end

    should_run_sqp = RUN_SQP_FIXED_SMOKE && case_def.tracking == 0;
    should_run_sqp = should_run_sqp || (RUN_SQP_SINGLE_AXIS_SMOKE && case_def.tracking == 1);
    if should_run_sqp
        serial = run_sqp_case(x0, params, lb_case, ub_case, SQP_MAX_ITER, SQP_MAX_FUN_EVALS);
        rows = append_summary_row(rows, case_def, "sqp_smoke", serial, empty_parallel_result(), params, serial.x);
        raw_runs.(case_def.name).sqp_serial = serial;
    end
end

summary = struct2table(rows);
summary_file = fullfile(results_dir, 'validation_summary.csv');
writetable(summary, summary_file);

raw_file = fullfile(results_dir, 'validation_raw.mat');
save(raw_file, 'raw_runs', 'summary', 'case_defs', 'GA_POP_SIZE', 'GA_MAX_GEN', ...
    'MOGA_POP_SIZE', 'MOGA_MAX_GEN', 'SQP_MAX_ITER', 'SQP_MAX_FUN_EVALS');

make_validation_plots(summary, results_dir);
report_file = write_validation_report(summary, results_dir, raw_file);

fprintf('\nValidation summary saved to:\n  %s\n', summary_file);
fprintf('Validation report saved to:\n  %s\n', report_file);
fprintf('Raw validation data saved to:\n  %s\n', raw_file);

function result = run_wrapper_case(x, params, use_parallel)
    result = empty_result(use_parallel);
    try
        if use_parallel
            pool = setup_parallel_pool(true);
            tic;
            future = parfeval(pool, @agrivoltaic_wrapper, 1, x, params);
            result.outputs = fetchOutputs(future);
            result.time_s = toc;
        else
            tic;
            result.outputs = agrivoltaic_wrapper(x, params);
            result.time_s = toc;
        end
        result.x = x;
        result.objective = -result.outputs(3);
        result.panels = result.outputs(7);
        result.exitflag = NaN;
        result.funccount = 1;
        result.ok = all(isfinite(result.outputs));
    catch ME
        result.error = string(ME.message);
        result.ok = false;
    end
end

function result = run_ga_case(vars, params, lb, ub, pop_size, max_gen, seed, use_parallel)
    result = empty_result(use_parallel);
    try
        setup_parallel_pool(use_parallel);
        num_vars = numel(lb);
        pop = build_initial_population(vars, params, lb, ub, pop_size, seed);
        [A, B, Aeq, Beq] = build_tracking_slew_constraints(num_vars, params);
        rng(seed);
        options = optimoptions('ga', ...
            'PopulationSize', pop_size, ...
            'MaxGenerations', max_gen, ...
            'MaxStallGenerations', max_gen, ...
            'FunctionTolerance', 1e-4, ...
            'Display', 'off', ...
            'InitialPopulationMatrix', pop, ...
            'UseParallel', use_parallel);

        tic;
        [x_best, fval, exitflag, output] = ga( ...
            @(x) agrivoltaic_social_cost_of_carbon_wrapper(x, params), ...
            num_vars, A, B, Aeq, Beq, lb, ub, [], options);
        result.time_s = toc;
        result.x = x_best;
        result.fval = fval;
        result.objective = -fval;
        result.exitflag = exitflag;
        result.funccount = output.funccount;
        result.outputs = agrivoltaic_wrapper(x_best, params);
        result.panels = result.outputs(7);
        result.ok = isfinite(result.objective) && all(isfinite(result.outputs));
    catch ME
        result.error = string(ME.message);
        result.ok = false;
    end
end

function result = run_moga_case(vars, params, lb, ub, pop_size, max_gen, seed, use_parallel)
    result = empty_result(use_parallel);
    try
        setup_parallel_pool(use_parallel);
        num_vars = numel(lb);
        pop = build_initial_population(vars, params, lb, ub, pop_size, seed);
        [A, B, Aeq, Beq] = build_tracking_slew_constraints(num_vars, params);
        rng(seed);
        options = optimoptions('gamultiobj', ...
            'PopulationSize', pop_size, ...
            'MaxGenerations', max_gen, ...
            'ParetoFraction', 0.35, ...
            'Display', 'off', ...
            'InitialPopulationMatrix', pop, ...
            'UseParallel', use_parallel);

        tic;
        [x_set, fval, exitflag, output, population, scores] = gamultiobj( ...
            @(x) moga_validation_objective(x, params), ...
            num_vars, A, B, Aeq, Beq, lb, ub, [], options);
        result.time_s = toc;

        best_idx = 1;
        best_social = -Inf;
        best_outputs = nan(1, 8);
        for i = 1:size(x_set, 1)
            outputs_i = agrivoltaic_wrapper(x_set(i, :), params);
            social_i = -outputs_i(3);
            if social_i > best_social
                best_social = social_i;
                best_idx = i;
                best_outputs = outputs_i;
            end
        end

        result.x = x_set(best_idx, :);
        result.fval = fval;
        result.objective = best_social;
        result.exitflag = exitflag;
        result.funccount = get_output_field(output, 'funccount');
        result.outputs = best_outputs;
        result.panels = best_outputs(7);
        result.population = population;
        result.scores = scores;
        result.ok = isfinite(result.objective) && all(isfinite(best_outputs));
    catch ME
        result.error = string(ME.message);
        result.ok = false;
    end
end

function result = run_sqp_case(x0, params, lb, ub, max_iter, max_fun_evals)
    result = empty_result(false);
    try
        [A, B, Aeq, Beq] = build_tracking_slew_constraints(numel(lb), params);
        options = optimoptions('fmincon', ...
            'Algorithm', 'sqp', ...
            'Display', 'off', ...
            'MaxIterations', max_iter, ...
            'MaxFunctionEvaluations', max_fun_evals, ...
            'StepTolerance', 1e-6, ...
            'OptimalityTolerance', 1e-6);

        tic;
        [x_best, fval, exitflag, output] = fmincon( ...
            @(x) agrivoltaic_social_cost_of_carbon_wrapper(x, params), ...
            x0, A, B, Aeq, Beq, lb, ub, [], options);
        result.time_s = toc;
        result.x = x_best;
        result.fval = fval;
        result.objective = -fval;
        result.exitflag = exitflag;
        result.funccount = output.funcCount;
        result.outputs = agrivoltaic_wrapper(x_best, params);
        result.panels = result.outputs(7);
        result.ok = isfinite(result.objective) && all(isfinite(result.outputs));
    catch ME
        result.error = string(ME.message);
        result.ok = false;
    end
end

function fitness = moga_validation_objective(x, params)
    raw = agriObjArray2Struct(agrivoltaic_wrapper(x, params));
    fitness = [-raw.emission_reduction, -raw.fiscal_profit, -raw.yearly_biomass];
end

function rows = append_summary_row(rows, case_def, test_type, serial, parallel, params, x_for_checks)
    if isempty(parallel.time_s) || isnan(parallel.time_s)
        speedup = NaN;
    else
        speedup = serial.time_s / parallel.time_s;
    end

    objective_delta = abs(serial.objective - parallel.objective);
    if isnan(parallel.objective)
        objective_delta = NaN;
    end

    if isempty(x_for_checks) && ~isempty(serial.x)
        x_for_checks = serial.x;
    end
    [check_ok, check_notes, max_slew_violation] = mode_checks(params, x_for_checks, serial, parallel, test_type);

    row = struct();
    row.test_type = string(test_type);
    row.mode_name = string(case_def.name);
    row.geometry_mode = string(case_def.geometry);
    row.tracking_mode = case_def.tracking;
    row.serial_time_s = serial.time_s;
    row.parallel_time_s = parallel.time_s;
    row.speedup_vs_serial = speedup;
    row.serial_objective = serial.objective;
    row.parallel_objective = parallel.objective;
    row.objective_abs_delta = objective_delta;
    row.serial_panels = serial.panels;
    row.parallel_panels = parallel.panels;
    row.serial_exitflag = serial.exitflag;
    row.parallel_exitflag = parallel.exitflag;
    row.serial_funccount = serial.funccount;
    row.parallel_funccount = parallel.funccount;
    row.max_slew_violation_rad = max_slew_violation;
    row.pass = serial.ok && (parallel.ok || isnan(parallel.time_s)) && check_ok;
    row.notes = strjoin([serial.error, parallel.error, check_notes], " | ");

    if isempty(rows)
        rows = row;
    else
        rows(end + 1) = row;
    end
end

function [ok, notes, max_slew_violation] = mode_checks(params, x, serial, parallel, test_type)
    ok = true;
    notes = strings(0);
    max_slew_violation = 0;

    if isempty(x)
        ok = false;
        notes(end + 1) = "No design vector available for mode checks.";
        return;
    end

    expected_n = expected_design_vector_length(params);
    if numel(x) ~= expected_n
        ok = false;
        notes(end + 1) = sprintf('Expected %d DVs, got %d.', expected_n, numel(x));
    end

    if isfield(params, 'geometry_mode') && params.geometry_mode == 1
        if abs(serial.panels - 550) > 1e-9
            ok = false;
            notes(end + 1) = sprintf('Expected 550 panels, got %.3f.', serial.panels);
        end

        var = agriVarArray2Struct(x, params);
        gcf = ground_cover_factor(var, params);
        expected_gcf = params.hedge_width / params.row_pitch;
        if abs(gcf - expected_gcf) > 1e-10
            ok = false;
            notes(end + 1) = sprintf('Expected row-centered GCF %.12g, got %.12g.', expected_gcf, gcf);
        end
    end

    if params.tracking_mode == 1 && test_type ~= "wrapper"
        [~, max_slew_violation] = validate_tracking_slew(x, params);
        if max_slew_violation > 1e-8
            ok = false;
            notes(end + 1) = sprintf('Slew violation %.3e rad.', max_slew_violation);
        end
    end

    if ~isnan(parallel.time_s) && ~isempty(serial.outputs) && ~isempty(parallel.outputs)
        diff_outputs = max(abs(serial.outputs - parallel.outputs));
        scale = max(1, max(abs(serial.outputs)));
        if diff_outputs > 1e-7 * scale
            notes(end + 1) = sprintf('Serial/parallel output delta %.3e.', diff_outputs);
        end
    end
end

function n = expected_design_vector_length(params)
    if isfield(params, 'geometry_mode') && params.geometry_mode == 1
        if params.tracking_mode == 1
            n = 99;
        else
            n = 4;
        end
    else
        if params.tracking_mode == 1
            n = 103;
        else
            n = 7;
        end
    end
end

function result = empty_parallel_result()
    result = empty_result(true);
end

function result = empty_result(use_parallel)
    result = struct();
    result.use_parallel = use_parallel;
    result.time_s = NaN;
    result.x = [];
    result.fval = NaN;
    result.objective = NaN;
    result.exitflag = NaN;
    result.funccount = NaN;
    result.outputs = [];
    result.panels = NaN;
    result.ok = false;
    result.error = "";
end

function value = get_output_field(output, fieldname)
    if isstruct(output) && isfield(output, fieldname)
        value = output.(fieldname);
    else
        value = NaN;
    end
end

function make_validation_plots(summary, results_dir)
    comparable = ~isnan(summary.parallel_time_s);
    if any(comparable)
        labels = summary.mode_name(comparable) + "_" + summary.test_type(comparable);
        times = [summary.serial_time_s(comparable), summary.parallel_time_s(comparable)];

        fig = figure('Color', 'w', 'Position', [100 100 1200 500]);
        bar(categorical(labels, labels), times);
        ylabel('Time (s)');
        title('Serial vs Parallel Runtime');
        legend({'Serial', 'Parallel'}, 'Location', 'best');
        xtickangle(35);
        grid on;
        saveas(fig, fullfile(results_dir, 'parallel_speedup.png'));
        close(fig);

        objectives = [summary.serial_objective(comparable), summary.parallel_objective(comparable)];
        fig = figure('Color', 'w', 'Position', [100 100 1200 500]);
        bar(categorical(labels, labels), objectives);
        ylabel('Best social value ($)');
        title('Serial vs Parallel Objective Comparison');
        legend({'Serial', 'Parallel'}, 'Location', 'best');
        xtickangle(35);
        grid on;
        saveas(fig, fullfile(results_dir, 'objective_comparison.png'));
        close(fig);
    end
end

function report_file = write_validation_report(summary, results_dir, raw_file)
    report_file = fullfile(results_dir, 'validation_report.md');
    fid = fopen(report_file, 'w');
    cleanup = onCleanup(@() fclose(fid));

    fprintf(fid, '# Parallel Validation Report\n\n');
    fprintf(fid, 'Generated: %s\n\n', datestr(now));
    fprintf(fid, 'Raw data: `%s`\n\n', raw_file);
    fprintf(fid, 'Rows run: %d\n', height(summary));
    fprintf(fid, 'Rows passed: %d\n\n', sum(summary.pass));

    fprintf(fid, '## Runtime Summary\n\n');
    for i = 1:height(summary)
        fprintf(fid, '- %s / %s: serial %.2fs, parallel %.2fs, speedup %.2fx, pass %d\n', ...
            summary.mode_name(i), summary.test_type(i), ...
            summary.serial_time_s(i), summary.parallel_time_s(i), ...
            summary.speedup_vs_serial(i), summary.pass(i));
    end

    failed = summary(~summary.pass, :);
    if ~isempty(failed)
        fprintf(fid, '\n## Failures Or Warnings\n\n');
        for i = 1:height(failed)
            fprintf(fid, '- %s / %s: %s\n', failed.mode_name(i), failed.test_type(i), failed.notes(i));
        end
    end
end
