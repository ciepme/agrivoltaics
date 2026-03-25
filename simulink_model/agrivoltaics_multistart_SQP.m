%% agrivoltaics_multistart_SQP.m
% for running SQP with many initial guesses
% also records and summarizes relevant info across runs and specifically
% for the best performing run
clear;
clc;
rng(40);

agrivoltaics_variable_definition;

% Load model and disable InitFcn to prevent clear from wiping workspace
load_system('agrivoltaics_v1');
set_param('agrivoltaics_v1', 'InitFcn', '');
set_param('agrivoltaics_v1', 'SimulationMode', 'normal'); % other option is accelerator, but i think it gets wonky with paralellization
%% Initial points
N_lhs = 71; %<--CHANGE THIS for num of initial points to test
lhs_samples = lhsdesign(N_lhs, 7);
X0_lhs = lb + lhs_samples .* (ub - lb);
X0_nominal = agriVarStruct2Array(agriVar);
X0_matrix = [X0_lhs; X0_nominal];
N_starts = size(X0_matrix, 1);

%% SQP Options
options = optimoptions('fmincon', ...
    'Algorithm',              'sqp', ...
    'FiniteDifferenceType',   'central', ...
    'ScaleProblem',           true, ...
    'MaxFunctionEvaluations', 2000, ...
    'StepTolerance',          1e-6, ...
    'OptimalityTolerance',    1e-6, ...
    'Display',                'off');

%% Storage
fvals       = nan(N_starts, 1);
exitflags   = nan(N_starts, 1);
func_counts = nan(N_starts, 1);
iterations  = nan(N_starts, 1);
times       = nan(N_starts, 1);
X_opts      = nan(N_starts, 7);
errors      = strings(N_starts, 1);
X0_starts   = nan(N_starts, 7);   
fvals_initial = nan(N_starts, 1); 


fprintf('Running %d SQP starts...\n', N_starts);

%% Parallel pool
if isempty(gcp('nocreate'))
    cluster = parcluster('local');
    n_workers = max(2, cluster.NumWorkers - 2);  % will depend on the computer
    parpool('local', n_workers);
end

%copy of agriParams for parallel workers
agriParams_copy = agriParams;

results_cell = cell(N_starts, 1);
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
            x0_i, [], [], [], [], lb, ub, [], options);
        results_cell{si} = {x0_i, x_i, f_i, flag_i, out_i.funcCount, out_i.iterations, toc, '', fval_i0};
    catch ME
        results_cell{si} = {x0_i, nan(1,7), NaN, NaN, NaN, NaN, toc, ME.message, fval_i0};
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
    fvals_initial(si) = results_cell{si}{9};  % initial fval
end

%% Selecting x*
converged_mask = exitflags >= 1;
if ~any(converged_mask)
    error('No runs converged');
end
[best_fval, best_idx] = min(fvals(converged_mask));
converged_idx = find(converged_mask);
best_run = converged_idx(best_idx);
x_star = X_opts(best_run, :);

%% Convergence Summary
n_conv = sum(converged_mask);
conv_fvals = fvals(converged_mask);
tol_pct = 0.01;
n_unique = sum(abs(conv_fvals - best_fval) > abs(best_fval) * tol_pct) + 1;
fprintf('\n==== MULTI-START CONVERGENCE SUMMARY===\n');
fprintf('Total starts         : %d\n', N_starts);
fprintf('Converged (flag>=1)  : %d / %d\n', n_conv, N_starts);
fprintf('Best social value    : $%.2f  (run #%d)\n', -best_fval, best_run);
fprintf('Worst converged      : $%.2f\n', -max(conv_fvals));
fprintf('Range of fval        : $%.2f\n', max(conv_fvals) - min(conv_fvals));
fprintf('Unique optima (~1%%) : %d\n', n_unique);
fprintf('\nOptimal Design x*:\n');
fprintf('  Panel Height  (z_p)  : %.3f m\n',       x_star(1));
fprintf('  Panel Length  (l_p)  : %.3f m\n',       x_star(2));
fprintf('  Panel Width   (w_p)  : %.3f m\n',       x_star(3));
fprintf('  Azimuth       (phi)  : %.3f rad (%.1f deg)\n', x_star(4), rad2deg(x_star(4)));
fprintf('  Tilt          (sigma): %.3f rad (%.1f deg)\n', x_star(5), rad2deg(x_star(5)));
fprintf('  Row Spacing   (y_p)  : %.3f m\n',       x_star(6));
fprintf('  Panel Gap     (x_p)  : %.3f m\n',       x_star(7));
fprintf('  Social Value         : $%.2f\n', -best_fval);

%% Value breakdown for x*
fprintf('\n');
print_value_breakdown(x_star, agriParams);

%% Save results
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
T = table((1:N_starts)', exitflags, ...
    fvals_initial, -fvals_initial, fvals, -fvals, ...
    func_counts, iterations, times, ...
    X0_starts(:,1), X0_starts(:,2), X0_starts(:,3), ...
    rad2deg(X0_starts(:,4)), rad2deg(X0_starts(:,5)), X0_starts(:,6), X0_starts(:,7), ...
    X_opts(:,1), X_opts(:,2), X_opts(:,3), ...
    rad2deg(X_opts(:,4)), rad2deg(X_opts(:,5)), X_opts(:,6), X_opts(:,7), ...
    errors, ...
    'VariableNames', {'run_id','exitflag', ...
    'fval_initial','social_value_initial','fval','social_value', ...
    'func_count','iterations','time_s', ...
    'x0_z_p','x0_l_p','x0_w_p','x0_phi_deg','x0_sigma_deg','x0_y_p','x0_x_p', ...
    'z_p','l_p','w_p','phi_deg','sigma_deg','y_p','x_p', ...
    'error'});
fname = sprintf('multistart_SQP_results_%s.csv', timestamp);
writetable(T, fname);
fprintf('\nResults saved to: %s\n', fname);
%% Save x* for sensitivity study
save('x_star.mat', 'x_star', 'best_fval');
fprintf('x_star saved to x_star.mat\n');

%% convergence trace (single SQP run from best starting point found)
fprintf('\nRunning convergence trace from best starting point...\n');
options_trace = optimoptions(options, 'Display', 'off', 'OutputFcn', @trace_outfun);
fmincon(@(x) agrivoltaic_social_cost_of_carbon_wrapper(x, agriParams), ...
    X0_starts(best_run,:), [], [], [], [], lb, ub, [], options_trace);
save('conv_history.mat', 'conv_history');
% plotting the convergence history
load('conv_history.mat');
fprintf('Convergence trace saved to conv_history.mat\n');
figure;
plot(conv_history(:,1), conv_history(:,2), 'r-o', 'MarkerSize', 4);
xlabel('Iteration'); ylabel('Social Value ($)');
title('SQP Convergence Trace from Best Starting Point');
grid on;



%% Restore InitFcn
set_param('agrivoltaics_v1', 'InitFcn', 'agrivoltaics_variable_definition');