%% Orthogonal Array DOE (L27) - 5 active factors
% Fixed:
%   phi = 0 rad (due south)
%   y_p = 4.0 m
%
% Active factors (3 levels each):
%   sigma [rad] : [0.35, 0.52, 0.70]
%   z_p   [m]   : [2.0, 2.5, 3.0]
%   w_p   [m]   : [0.9, 1.0, 1.1]
%   l_p   [m]   : [1.8, 2.0, 2.2]
%   x_p   [m]   : [0.05, 0.10, 0.15]
%
% Output:
%   CSV with run settings + emissions + profit

clear;
clc;

base_dir = fileparts(mfilename('fullpath'));
if isempty(base_dir)
    base_dir = pwd;
end
cd(base_dir);

% Initialize model variables exactly as your normal workflow does.
agrivoltaics_variable_definition;
clear variable_bus parameter_bus inequality_constraint_bus equality_constraint_bus;

modelName = "agrivoltaics_v1.slx";
modelBase = "agrivoltaics_v1";

% Avoid re-running InitFcn during each sim call.
origInitFcn = get_param(modelBase, 'InitFcn');
set_param(modelBase, 'InitFcn', '');
init_cleanup = onCleanup(@() set_param(modelBase, 'InitFcn', origInitFcn)); %#ok<NASGU>

% ---------- Factor levels ----------
sigma_levels = [0.35, 0.52, 0.70];      % rad
z_p_levels   = [2.0, 2.5, 3.0];         % m
w_p_levels   = [0.9, 1.0, 1.1];         % m
l_p_levels   = [1.8, 2.0, 2.2];         % m
x_p_levels   = [0.05, 0.10, 0.15];      % m

phi_fixed = 0.0;   % rad
y_p_fixed = 4.0;   % m

% ---------- L27 OA index matrix (levels 1..3) ----------
% Build a valid 27-run OA from three base 3-level columns and two
% additional independent modulo-3 combinations.
[A, B, C] = ndgrid(1:3, 1:3, 1:3);
c1 = A(:);
c2 = B(:);
c3 = C(:);
c4 = mod((c1 - 1) + (c2 - 1), 3) + 1;
c5 = mod((c1 - 1) + (c3 - 1), 3) + 1;
L = [c1, c2, c3, c4, c5];  % 27 x 5

% Map OA indices to physical factor values.
n_runs = size(L, 1);
sigma_vals = sigma_levels(L(:,1)).';
z_p_vals   = z_p_levels(L(:,2)).';
w_p_vals   = w_p_levels(L(:,3)).';
l_p_vals   = l_p_levels(L(:,4)).';
x_p_vals   = x_p_levels(L(:,5)).';

S_out = nan(n_runs, 1);
R_out = nan(n_runs, 1);
status = strings(n_runs, 1);
error_message = strings(n_runs, 1);

for i = 1:n_runs
    try
        % Apply OA scenario values.
        var.PV.sigma = sigma_vals(i);
        var.PV.phi   = phi_fixed;
        var.PV.z_p   = z_p_vals(i);
        var.PV.w_p   = w_p_vals(i);
        var.PV.l_p   = l_p_vals(i);
        var.PV.x_p   = x_p_vals(i);
        var.PV.y_p   = y_p_fixed;

        % Simulate and read outputs (same approach as housekeeping script).
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
    L(:,1), L(:,2), L(:,3), L(:,4), L(:,5), ...
    sigma_vals, ...
    repmat(phi_fixed, n_runs, 1), ...
    z_p_vals, w_p_vals, l_p_vals, x_p_vals, ...
    repmat(y_p_fixed, n_runs, 1), ...
    S_out, R_out, status, error_message, ...
    'VariableNames', { ...
    'run_id', ...
    'sigma_level', 'z_p_level', 'w_p_level', 'l_p_level', 'x_p_level', ...
    'sigma_rad', 'phi_rad', 'z_p_m', 'w_p_m', 'l_p_m', 'x_p_m', 'y_p_m', ...
    'total_co2eq_displaced', 'total_profit', 'status', 'error_message'});

out_csv = fullfile(base_dir, 'doe_oa_l27_5factor_results.csv');
writetable(results, out_csv);

timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
out_csv_ts = fullfile(base_dir, ['doe_oa_l27_5factor_results_' timestamp '.csv']);
writetable(results, out_csv_ts);

disp('OA DOE finished. Results written to simulink_model folder.');
