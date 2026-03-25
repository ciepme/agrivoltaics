%% run_sensitivity.m
% includes design variable one factor at a time (OFAT), parameter OFAT
% Saves results as CSVs

%% 0. Initializing workspace
agrivoltaics_variable_definition;   % sets agriParams, lb, ub, bus objects
load('x_star.mat');                  % loads x_star [1x7] and best_fval
f_star = -best_fval;                 % positive social value at optimum
fprintf('f_star = $%.2f\n', f_star);

%% 1.Loading model and starting parallel pool 
load_system('agrivoltaics_v1');
set_param('agrivoltaics_v1', 'InitFcn', '');
set_param('agrivoltaics_v1', 'SimulationMode', 'normal');

if isempty(gcp('nocreate'))
    cluster   = parcluster('local');
    n_workers = max(2, cluster.NumWorkers - 2);
    parpool('local', n_workers);
end

%for parfor use
agriParams_pass = agriParams;
x_star_pass     = x_star;

%% 2. Design Variable sensitivity (OFAT)

dv_names = {'z_p','l_p','w_p','phi','sigma','y_p','x_p'};
n_dv     = 7;
% Fractions of total range (ub-lb): dense near x_star, always reaches
% bounds-->
% Both directions always attempted; clamped to [lb, ub] if x_star is near a bound.
% Duplicate points (when clamped) are fine -they evaluate to ~f_star and show
% as flat at 0% on the plot
fracs = [0.05, 0.10, 0.20, 0.35,0.50, 0.70, 1.0];

% Building flat evaluation list: each row is one perturbed design point
all_x_new_dv  = zeros(0, n_dv);  % n_evals x 7
all_var_idx   = zeros(0, 1);     % which DV was perturbed (1-7)

for i = 1:n_dv
    total_range = ub(i) - lb(i);
    for fi = 1:numel(fracs)
        step = fracs(fi) * total_range;

        % Downward perturbation
        x_new    = x_star;
        x_new(i) = max(lb(i), x_star(i) - step);
        all_x_new_dv = [all_x_new_dv; x_new]; 
        all_var_idx  = [all_var_idx;  i];      

        % Upward perturbation
        x_new    = x_star;
        x_new(i) = min(ub(i), x_star(i) + step);
        all_x_new_dv = [all_x_new_dv; x_new]; 
        all_var_idx  = [all_var_idx;  i];      
    end
end

n_evals_dv = size(all_x_new_dv, 1);
sv_dv_raw  = zeros(n_evals_dv, 1);
fprintf('DV sensitivity: %d evaluations...\n', n_evals_dv);

parfor k = 1:n_evals_dv
    sv_dv_raw(k) = -agrivoltaic_social_cost_of_carbon_wrapper( ...
        all_x_new_dv(k,:), agriParams_pass);
end

%Assembling per-variable result structs
dv_results(n_dv) = struct();
for i = 1:n_dv
    mask = (all_var_idx == i);

    pv = [x_star(i);        all_x_new_dv(mask, i)];
    sv = [f_star;            sv_dv_raw(mask)];
    pc = [0;                 (sv_dv_raw(mask) - f_star) ./ f_star .* 100];

    % Sorting by perturbed value for clean plotting
    [pv, sidx] = sort(pv);
    sv = sv(sidx);
    pc = pc(sidx);

    % Normalized sensitivity: skip x_star itself (delta==0) and if x_star(i)==0
    ns = nan(size(pv));
    if x_star(i) ~= 0
        for j = 1:numel(pv)
            delta = pv(j) - x_star(i);
            if abs(delta) > 1e-12
                ns(j) = pc(j) / (abs(delta) / abs(x_star(i)));
            end
        end
    end

    dv_results(i).var_name          = dv_names{i};
    dv_results(i).x_star_val        = x_star(i);
    dv_results(i).perturbed_values  = pv;
    dv_results(i).social_values     = sv;
    dv_results(i).percent_changes   = pc;
    dv_results(i).norm_sensitivities = ns;
end

%% 3. Parameter Sensitivity (OFAT)

param_names     = {'crop_RUE','crop_elec_price','crop_price', ...
                   'PV_n_p','crop_T_max','crop_GCF'};
param_baselines = [agriParams.crop_RUE, agriParams.crop_elec_price, ...
                   agriParams.crop_price, agriParams.PV_n_p, ...
                   agriParams.crop_T_max, agriParams.crop_GCF];
pct_levels      = [-0.50, -0.25, -0.20, -0.10, -0.05, 0.05, 0.10, 0.20, 0.50, 1.0];
abs_T_levels    = [-15, -10, -6, -3, -1.5, 1.5, 3, 6, 10,15];   % °C absolute for crop_T_max
n_params        = numel(param_names);
n_levels        = numel(pct_levels);

% Pre building cell array of modified agriParams so parfor gets clean copies
n_param_evals   = n_params * n_levels;
param_configs   = cell(n_param_evals, 1);
param_eval_map  = zeros(n_param_evals, 2);  % [param_idx, level_idx]

k = 0;
for p = 1:n_params
    for lv = 1:n_levels
        k = k + 1;
        ap_mod = agriParams;
        if strcmp(param_names{p}, 'crop_T_max')
            new_val = param_baselines(p) + abs_T_levels(lv);
        else
            new_val = param_baselines(p) * (1 + pct_levels(lv));
        end
        ap_mod.(param_names{p}) = new_val;
        param_configs{k}    = ap_mod;
        param_eval_map(k,:) = [p, lv];
    end
end

sv_param_raw = zeros(n_param_evals, 1);
fprintf('Parameter sensitivity: %d evaluations...\n', n_param_evals);

parfor k = 1:n_param_evals
    sv_param_raw(k) = -agrivoltaic_social_cost_of_carbon_wrapper( ...
        x_star_pass, param_configs{k});
end

% Assembling per-parameter result structs
param_results(n_params) = struct(); 
for p = 1:n_params
    mask = (param_eval_map(:,1) == p);

    if strcmp(param_names{p}, 'crop_T_max')
        pv = [param_baselines(p);  (param_baselines(p) + abs_T_levels(:))];
    else
        pv = [param_baselines(p);  param_baselines(p) .* (1 + pct_levels(:))];
    end
    sv = [f_star; sv_param_raw(mask)];
    pc = [0;      (sv_param_raw(mask) - f_star) ./ f_star .* 100];

    [pv, sidx] = sort(pv);
    sv = sv(sidx);
    pc = pc(sidx);

    param_results(p).param_name       = param_names{p};
    param_results(p).baseline_val     = param_baselines(p);
    param_results(p).perturbed_values = pv;
    param_results(p).social_values    = sv;
    param_results(p).percent_changes  = pc;
end

%% 4. Restoring model InitFcn
set_param('agrivoltaics_v1', 'InitFcn', 'agrivoltaics_variable_definition');

%% 5. Saving results

base_dir    = fileparts(mfilename('fullpath'));
if isempty(base_dir), base_dir = pwd; end
results_dir = fullfile(base_dir, 'results');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end

ts = datestr(now, 'yyyymmdd_HHMMSS');  

dv_mat_file    = fullfile(results_dir, ['dv_sensitivity_'    ts '.mat']);
param_mat_file = fullfile(results_dir, ['param_sensitivity_' ts '.mat']);
dv_csv_file    = fullfile(results_dir, ['dv_sensitivity_'    ts '.csv']);
param_csv_file = fullfile(results_dir, ['param_sensitivity_' ts '.csv']);

save(dv_mat_file,    'dv_results');
save(param_mat_file, 'param_results');

% DV sensitivity CSV
% Columns: var_name, perturbed_value, social_value,percent_change, norm_sensitivity
% phi and sigma perturbed values converted to degrees for readability
angle_dv = {'phi','sigma'};
fid = fopen(dv_csv_file, 'w');
fprintf(fid, 'var_name,perturbed_value,social_value,percent_change,norm_sensitivity\n');
for i = 1:n_dv
    pv = dv_results(i).perturbed_values;
    sv = dv_results(i).social_values;
    pc = dv_results(i).percent_changes;
    ns = dv_results(i).norm_sensitivities;
    if ismember(dv_results(i).var_name, angle_dv)
        pv_out = rad2deg(pv);
    else
        pv_out = pv;
    end
    for j = 1:numel(pv)
        if isnan(ns(j))
            fprintf(fid, '%s,%.6f,%.4f,%.4f,NaN\n', ...
                dv_results(i).var_name, pv_out(j), sv(j), pc(j));
        else
            fprintf(fid, '%s,%.6f,%.4f,%.4f,%.4f\n', ...
                dv_results(i).var_name, pv_out(j), sv(j), pc(j), ns(j));
        end
    end
end
fclose(fid);

% Parameter sensitivity CSV
% Columns: param_name, perturbed_value, social_value, percent_change
fid = fopen(param_csv_file, 'w');
fprintf(fid, 'param_name,perturbed_value,social_value,percent_change\n');
for p = 1:n_params
    pv = param_results(p).perturbed_values;
    sv = param_results(p).social_values;
    pc = param_results(p).percent_changes;
    for j = 1:numel(pv)
        fprintf(fid, '%s,%.6f,%.4f,%.4f\n', ...
            param_results(p).param_name, pv(j), sv(j), pc(j));
    end
end
fclose(fid);

fprintf('\nDone. results saved to: %s\n', results_dir);
fprintf('  %s\n  %s\n  %s\n  %s\n', dv_mat_file,param_mat_file, dv_csv_file, param_csv_file);
