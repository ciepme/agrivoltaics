%% DOE parameter study (15 runs)
% Simple workflow:
% 1) initialize model variables exactly as usual
% 2) update decision variables
% 3) run model
% 4) save outputs to CSV

clear;
clc;

base_dir = fileparts(mfilename('fullpath'));
if isempty(base_dir)
    base_dir = pwd;
end
cd(base_dir);

% Initialize params/var/buses the same way your model already uses.
agrivoltaics_variable_definition;
clear variable_bus parameter_bus inequality_constraint_bus equality_constraint_bus;

modelName = "agrivoltaics_v1.slx";
modelBase = "agrivoltaics_v1";

% The model InitFcn re-runs agrivoltaics_variable_definition (which calls clear).
% Disable it during DOE so each run only uses the values we set below.
origInitFcn = get_param(modelBase, 'InitFcn');
set_param(modelBase, 'InitFcn', '');
init_cleanup = onCleanup(@() set_param(modelBase, 'InitFcn', origInitFcn)); %#ok<NASGU>

% Factor order: [sigma_deg, phi_deg, z_p, w_p, l_p, x_p, y_p]
low_vals = [0, -45, 1, 0.5, 0.5, 0.02, 0.5];
nom_vals = [45,   0, 2, 1.0, 1.0, 0.10, 1];
high_vals = [90, 45, 3, 1.5, 1.5, 0.25, 4];
factor_names = ["sigma_deg","phi_deg","z_p_m","w_p_m","l_p_m","x_p_m","y_p_m"];

% 15-run OFAT matrix: nominal + low/high for each factor.
X = nom_vals;
case_name = "nominal";
changed_factor = "none";
changed_level = "nominal";

for j = 1:numel(factor_names)
    x_low = nom_vals;
    x_low(j) = low_vals(j);
    X = [X; x_low]; %#ok<AGROW>
    case_name = [case_name; factor_names(j) + "_low"]; %#ok<AGROW>
    changed_factor = [changed_factor; factor_names(j)]; %#ok<AGROW>
    changed_level = [changed_level; "low"]; %#ok<AGROW>

    x_high = nom_vals;
    x_high(j) = high_vals(j);
    X = [X; x_high]; %#ok<AGROW>
    case_name = [case_name; factor_names(j) + "_high"]; %#ok<AGROW>
    changed_factor = [changed_factor; factor_names(j)]; %#ok<AGROW>
    changed_level = [changed_level; "high"]; %#ok<AGROW>
end

n_runs = size(X, 1);
S_out = nan(n_runs, 1);
R_out = nan(n_runs, 1);
status = strings(n_runs, 1);
error_message = strings(n_runs, 1);

for i = 1:n_runs
    try
        % Update decision variables for this run.
        var.PV.sigma = deg2rad(X(i,1));
        var.PV.phi   = deg2rad(X(i,2));
        var.PV.z_p   = X(i,3);
        var.PV.w_p   = X(i,4);
        var.PV.l_p   = X(i,5);
        var.PV.x_p   = X(i,6);
        var.PV.y_p   = X(i,7);

        % Run model and read outputs the same way as agrivoltaics_housekeeping.m
        simOut = sim(modelName);
        s_data = simOut.yout{1}.Values.Data;
        r_data = simOut.yout{2}.Values.Data;

        S_out(i) = s_data(end);
        R_out(i) = r_data(end);
        status(i) = "ok";
    catch ME
        status(i) = "error";
        error_message(i) = string(getReport(ME, 'extended', 'hyperlinks', 'off'));
    end
end

results = table( ...
    (1:n_runs).', ...
    case_name, changed_factor, changed_level, ...
    X(:,1), X(:,2), X(:,3), X(:,4), X(:,5), X(:,6), X(:,7), ...
    S_out, R_out, status, error_message, ...
    'VariableNames', { ...
    'run_id', 'case_name', 'changed_factor', 'changed_level', ...
    'sigma_deg', 'phi_deg', 'z_p_m', 'w_p_m', 'l_p_m', 'x_p_m', 'y_p_m', ...
    'total_co2eq_displaced', 'total_profit', 'status', 'error_message'});

writetable(results, fullfile(base_dir, 'doe_parameter_study_15_results.csv'));
writetable(results, fullfile(base_dir, ['doe_parameter_study_15_results_' char(datetime('now','Format','yyyyMMdd_HHmmss')) '.csv']));

disp('DOE finished. Results CSV written to simulink_model folder.');
