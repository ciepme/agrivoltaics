%% agrivoltaics_multistart_SQP.m
% for running SQP with many initial guesses
% also records and summarizes relevant info across runs and specifically
% for the best performing run
clear;
clc;
rng(39);

agrivoltaics_variable_definition;
accept_exitflag_zero = true;

% Load model and disable InitFcn to prevent clear from wiping workspace
load_system('agrivoltaics_v1');
set_param('agrivoltaics_v1', 'InitFcn', '');
set_param('agrivoltaics_v1', 'SimulationMode', 'normal'); % other option is accelerator, but i think it gets wonky with paralellization
%% Initial points
num_vars = numel(lb);
[A, B, Aeq, Beq] = build_tracking_slew_constraints(num_vars, agriParams);
N_lhs = 1; %<--CHANGE THIS for num of initial points to test
lhs_samples = lhsdesign(N_lhs, num_vars);
X0_lhs = lb + lhs_samples .* (ub - lb);
X0_nominal = agriVarStruct2Array(agriVar, agriParams);
X0_matrix = [X0_lhs; X0_nominal];
N_starts = size(X0_matrix, 1);

%% SQP Options
if num_vars > 7
    fd_type = 'central';
    max_fun_evals = 600;
    max_iters = 50;
    step_tol = 1e-5;
else
    fd_type = 'central';
    max_fun_evals = 2000;
    max_iters = 400;
    step_tol = 1e-6;
end

options = optimoptions('fmincon', ...
    'Algorithm',              'sqp', ...
    'FiniteDifferenceType',   fd_type, ...
    'ScaleProblem',           true, ...
    'MaxFunctionEvaluations', max_fun_evals, ...
    'MaxIterations',          max_iters, ...
    'StepTolerance',          step_tol, ...
    'OptimalityTolerance',    1e-6, ...
    'Display',                'off');

%% Storage
fvals       = nan(N_starts, 1);
exitflags   = nan(N_starts, 1);
func_counts = nan(N_starts, 1);
iterations  = nan(N_starts, 1);
times       = nan(N_starts, 1);
X_opts      = nan(N_starts, num_vars);
errors      = strings(N_starts, 1);
messages    = strings(N_starts, 1);
X0_starts   = nan(N_starts, num_vars);   
fvals_initial = nan(N_starts, 1); 
model_evals = nan(N_starts, 1);


fprintf('Running %d SQP starts...\n', N_starts);

%% Parallel pool
use_parallel = false; % switch for using parallel workers versus regular evaluation
if use_parallel
    if isempty(gcp('nocreate'))
        cluster = parcluster('local');
        n_workers = max(2, cluster.NumWorkers - 2);  % will depend on the computer
        parpool('local', n_workers);
    end
end

%copy of agriParams for parallel workers
agriParams_copy = agriParams;

results_cell = cell(N_starts, 1);
if use_parallel
    parfor si = 1:N_starts
        tic;
        x0_i = X0_matrix(si,:);
        try
            fval_i0 = agrivoltaic_social_cost_of_carbon_wrapper(x0_i, agriParams_copy);
        catch
            fval_i0 = NaN;
        end
        try
            [x_i, f_i, flag_i, out_i] = fmincon( ...
                @(x) agrivoltaic_social_cost_of_carbon_wrapper(x, agriParams_copy), ...
                x0_i, A, B, Aeq, Beq, lb, ub, [], options);
            results_cell{si} = {x0_i, x_i, f_i, flag_i, out_i.funcCount, out_i.iterations, toc, '', string(out_i.message), fval_i0, NaN};
        catch ME
            results_cell{si} = {x0_i, nan(1,num_vars), NaN, NaN, NaN, NaN, toc, string(ME.message), "", fval_i0, NaN};
        end
    end
else
    for si = 1:N_starts
        fprintf('Starting SQP run %d of %d...\n', si, N_starts);
        tic;
        x0_i = X0_matrix(si,:);
        setappdata(0, 'multistart_model_eval_count', 0);
        objective_i = @(x) multistart_counting_objective(x, agriParams_copy);
        options_i = optimoptions(options, ...
            'OutputFcn', @(x, optimValues, state) multistart_progress_outfun(x, optimValues, state, si, N_starts));
        try
            fval_i0 = objective_i(x0_i);
        catch
            fval_i0 = NaN;
        end
        try
            [x_i, f_i, flag_i, out_i] = fmincon( ...
                objective_i, ...
                x0_i, A, B, Aeq, Beq, lb, ub, [], options_i);
            model_evals_i = getappdata(0, 'multistart_model_eval_count');
            results_cell{si} = {x0_i, x_i, f_i, flag_i, out_i.funcCount, out_i.iterations, toc, '', string(out_i.message), fval_i0, model_evals_i};
        catch ME
            model_evals_i = getappdata(0, 'multistart_model_eval_count');
            results_cell{si} = {x0_i, nan(1,num_vars), NaN, NaN, NaN, NaN, toc, string(ME.message), "", fval_i0, model_evals_i};
        end
        fprintf('Finished SQP run %d of %d. Model evaluations: %d.\n', ...
            si, N_starts, results_cell{si}{11});
    end
end

fprintf('All starts complete-->unpacking results.....\n');

%% Unpack results_cell
for si = 1:N_starts
    X0_starts(si,:)   = results_cell{si}{1};  % initial guess
    X_opts(si,:)      = results_cell{si}{2};
    fvals(si)         = results_cell{si}{3};
    exitflags(si)     = results_cell{si}{4};
    func_counts(si)   = results_cell{si}{5};
    iterations(si)    = results_cell{si}{6};
    times(si)         = results_cell{si}{7};
    errors(si)        = results_cell{si}{8};
    messages(si)      = results_cell{si}{9};
    fvals_initial(si) = results_cell{si}{10};  % initial fval
    model_evals(si)   = results_cell{si}{11};
end

%% Selecting x*
strict_mask = exitflags >= 1;
usable_mask = isfinite(fvals) & ((exitflags >= 1) | (accept_exitflag_zero & exitflags == 0));

if any(strict_mask)
    converged_mask = strict_mask;
elseif any(usable_mask)
    converged_mask = usable_mask;
    warning(['No runs reached exitflag >= 1. ', ...
        'Using best run from usable set (including exitflag==0). ', ...
        'Inspect message/error columns and consider raising budgets if needed.']);
else
    fprintf('\nNo usable runs found. Exitflags:\n');
    disp(exitflags.');
    first_err = find(strlength(errors) > 0, 1, 'first');
    if ~isempty(first_err)
        fprintf('First error (run %d): %s\n', first_err, errors(first_err));
    end
    first_msg = find(strlength(messages) > 0, 1, 'first');
    if ~isempty(first_msg)
        fprintf('First fmincon message (run %d): %s\n', first_msg, messages(first_msg));
    end
    error('No runs converged or produced usable finite objectives.');
end
[best_fval, best_idx] = min(fvals(converged_mask));
converged_idx = find(converged_mask);
best_run = converged_idx(best_idx);
x_star = X_opts(best_run, :);
results_star = agrivoltaic_wrapper(x_star, agriParams);
E_star = results_star(1);
P_star = results_star(2);
pv_revenue_star = results_star(4);
crop_revenue_star = results_star(5);
yearly_biomass_star = results_star(6);
total_panels_star = results_star(7);
emissions_value_star = 190 * (E_star / 1000);
social_value_star = -best_fval;

%% Convergence Summary
n_strict = sum(strict_mask);
n_usable = sum(converged_mask);
conv_fvals = fvals(converged_mask);
tol_pct = 0.01;
n_unique = sum(abs(conv_fvals - best_fval) > abs(best_fval) * tol_pct) + 1;
fprintf('\n==== MULTI-START CONVERGENCE SUMMARY===\n');
fprintf('Total starts         : %d\n', N_starts);
fprintf('Converged (flag>=1)  : %d / %d\n', n_strict, N_starts);
if ~any(strict_mask) && any(converged_mask)
    fprintf('Usable fallback runs : %d / %d  (includes exitflag==0)\n', n_usable, N_starts);
end
fprintf('Best social value    : $%.2f  (run #%d)\n', -best_fval, best_run);
fprintf('Worst converged      : $%.2f\n', -max(conv_fvals));
fprintf('Range of fval        : $%.2f\n', max(conv_fvals) - min(conv_fvals));
fprintf('Unique optima (~1%%) : %d\n', n_unique);
fprintf('\nOptimal Design x*:\n');
x_star_var = agriVarArray2Struct(x_star, agriParams);
print_design_summary(x_star_var, agriParams);
if agriParams.tracking_mode == 1
    tracking_idx = tracking_angle_indices(agriParams);
    tracking_x_star = x_star(tracking_idx);
    fprintf('  Tracking vars        : %d values (min %.3f rad, max %.3f rad)\n', ...
        numel(tracking_x_star), min(tracking_x_star), max(tracking_x_star));
end
fprintf('  Social Value         : $%.2f\n', -best_fval);

%% Value breakdown for x*
fprintf('\n');
print_value_breakdown(x_star, agriParams);

%% Save results
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
base_dir = fileparts(mfilename('fullpath'));
if isempty(base_dir), base_dir = pwd; end
results_dir = fullfile(base_dir, 'results', 'SQP_multi');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end
T = table((1:N_starts)', exitflags, ...
    fvals_initial, -fvals_initial, fvals, -fvals, ...
    func_counts, model_evals, iterations, times, messages, ...
    'VariableNames', {'run_id','exitflag', ...
    'fval_initial','social_value_initial','fval','social_value', ...
    'func_count','model_evals','iterations','time_s','message'});

base_names = local_design_variable_names(agriParams);
for vi = 1:min(numel(base_names), num_vars)
    x0_now = X0_starts(:,vi);
    x_now = X_opts(:,vi);
    if ismember(base_names{vi}, {'phi','sigma'})
        T.(sprintf('x0_%s_deg', base_names{vi})) = rad2deg(x0_now);
        T.(sprintf('%s_deg', base_names{vi})) = rad2deg(x_now);
    else
        T.(sprintf('x0_%s', base_names{vi})) = x0_now;
        T.(base_names{vi}) = x_now;
    end
end

if agriParams.tracking_mode == 1
    tracking_idx = tracking_angle_indices(agriParams);
    for local_ti = 1:numel(tracking_idx)
        ti = tracking_idx(local_ti);
        track_name = sprintf('tracking_%02d', local_ti);
        T.(sprintf('x0_%s_deg', track_name)) = rad2deg(X0_starts(:,ti));
        T.(sprintf('%s_deg', track_name)) = rad2deg(X_opts(:,ti));
    end
end

T.error = errors;
fname = fullfile(results_dir, sprintf('multistart_SQP_results_%s.csv', timestamp));
writetable(T, fname);
fprintf('\nResults saved to: %s\n', fname);
%% Save x* for sensitivity study
x_star_file = fullfile(results_dir, 'x_star.mat');
save('x_star.mat', 'x_star', 'best_fval', 'E_star', 'P_star', ...
    'pv_revenue_star', 'crop_revenue_star', 'yearly_biomass_star', ...
    'total_panels_star', 'emissions_value_star', 'social_value_star');
save(x_star_file, 'x_star', 'best_fval', 'E_star', 'P_star', ...
    'pv_revenue_star', 'crop_revenue_star', 'yearly_biomass_star', ...
    'total_panels_star', 'emissions_value_star', 'social_value_star');
fprintf('x_star saved to x_star.mat and %s\n', x_star_file);

%% convergence trace (single SQP run from best starting point found)
fprintf('\nRunning convergence trace from best starting point...\n');
conv_history = [];
options_trace = optimoptions(options, 'Display', 'off', 'OutputFcn', @trace_outfun);
fmincon(@(x) agrivoltaic_social_cost_of_carbon_wrapper(x, agriParams), ...
    X0_starts(best_run,:), A, B, Aeq, Beq, lb, ub, [], options_trace);
conv_history_file = fullfile(results_dir, 'conv_history.mat');
save(conv_history_file, 'conv_history');
% plotting the convergence history
load(conv_history_file);
fprintf('Convergence trace saved to %s\n', conv_history_file);
figure;
plot(conv_history(:,1), conv_history(:,2), 'r-o', 'MarkerSize', 4);
xlabel('Iteration'); ylabel('Social Value ($)');
title('SQP Convergence Trace from Best Starting Point');
grid on;



%% Restore InitFcn
set_param('agrivoltaics_v1', 'InitFcn', 'agrivoltaics_variable_definition');

function stop = multistart_progress_outfun(~, optimValues, state, run_idx, total_runs)
    stop = false;
    if strcmp(state, 'iter')
        eval_count = getappdata(0, 'multistart_model_eval_count');
        social_value = -optimValues.fval;
        fprintf('  SQP run %d of %d | iteration %d complete | model evals %d | social value $%.2f.\n', ...
            run_idx, total_runs, optimValues.iteration, eval_count, social_value);
    end
end

function f = multistart_counting_objective(x, agriParams)
    eval_count = getappdata(0, 'multistart_model_eval_count');
    if isempty(eval_count)
        eval_count = 0;
    end
    setappdata(0, 'multistart_model_eval_count', eval_count + 1);
    f = agrivoltaic_social_cost_of_carbon_wrapper(x, agriParams);
end

function names = local_design_variable_names(agriParams)
    if isfield(agriParams, 'geometry_mode') && agriParams.geometry_mode == 1
        if agriParams.tracking_mode == 1
            names = {'z_p','w_p','d_norm'};
        else
            names = {'z_p','w_p','sigma','d_norm'};
        end
    else
        names = {'z_p','l_p','w_p','phi','sigma','y_p','x_p'};
    end
end
