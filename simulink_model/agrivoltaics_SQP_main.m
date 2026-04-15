
clear;
clc;
rng default;

agrivoltaics_variable_definition;

x0 = agriVarStruct2Array(agriVar, agriParams); %initial guess pulled from variable definition file
% use_scaled_SQP = true; % adding switch for running with scaling or without
% scale = [1e2, 1e2, 1e2, 1e3, 1e5, 1e5, 1e5]; % scaling values derived from the Hessian at last x_star
% x0_scaled = x0 .* scale;
% lb_scaled = lb .* scale;
% ub_scaled = ub .* scale;

% running SQP objective function with scaled DVs
% obj_scaled = @(x_scaled) agrivoltaic_social_cost_of_carbon_wrapper(x_scaled ./ scale, agriParams);


conv_history = [];
options = optimoptions('fmincon','Algorithm', 'sqp','Display', 'iter', ...
    'StepTolerance', 1e-6,'OptimalityTolerance', 1e-6, ...
    'OutputFcn', @trace_outfun);

tic; % Start timer
% if use_scaled_SQP
%     [x_opt_scaled, fval, exitflag, output] = fmincon(obj_scaled, x0_scaled, [], [], [], [], lb_scaled, ub_scaled, [], options);
%     x_opt = x_opt_scaled ./ scale;
% else
%     [x_opt, fval, exitflag, output] = fmincon(@(x) agrivoltaic_social_cost_of_carbon_wrapper(x, agriParams), x0, [], [], [], [], lb, ub, [], options);
% end
[x_opt, fval, exitflag, output] = fmincon(@(x) agrivoltaic_social_cost_of_carbon_wrapper(x, agriParams), x0, [], [], [], [], lb, ub, [], options);
time_taken = toc;

% Results
fprintf('Exit Flag: %d\n', exitflag);
fprintf('Time Taken: %.2f seconds\n', time_taken);
fprintf('Maximized Social Value: $%.2f\n', -fval); % Flip sign back to positive
disp(' ');
disp('Optimal Variables Found:');
fprintf('  Panel Height (z_p) : %.4f m\n', x_opt(1));
fprintf('  Panel Length (l_p) : %.4f m\n', x_opt(2));
fprintf('  Panel Width (w_p)  : %.4f m\n', x_opt(3));
fprintf('  Azimuth (phi)      : %.4f rad (%.4f deg)\n', x_opt(4), rad2deg(x_opt(4)));
fprintf('  Tilt (sigma)       : %.4f rad (%.4f deg)\n', x_opt(5), rad2deg(x_opt(5)));
fprintf('  Row Spacing (y_p)  : %.6f m\n', x_opt(6));
fprintf('  Panel Gap (x_p)    : %.6f m\n', x_opt(7));

if exist('conv_history', 'var') && ~isempty(conv_history)
    figure;
    plot(conv_history(:,1), conv_history(:,2), 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4);
    xlabel('Iteration');
    ylabel('Social Value ($)');
    title('SQP Convergence Plot');
    grid on;

    base_dir = fileparts(mfilename('fullpath'));
    if isempty(base_dir), base_dir = pwd; end
    results_dir = fullfile(base_dir, 'results', 'SQP');
    if ~exist(results_dir, 'dir'), mkdir(results_dir); end
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    saveas(gcf, fullfile(results_dir, sprintf('convergence_plot_%s.png', timestamp)));
end

results = agrivoltaic_wrapper(x_opt, agriParams);
E = results(1);
P = results(2);
social_cost = results(3);
pv_revenue = results(4);
crop_revenue = results(5);
yearly_biomass = results(6);
total_panels = results(7);
social_value = -social_cost;
emissions_value = 190 * (E / 1000);

fprintf('\ncurrent design vector from SQP best solution:\n');
format long g;
disp(x_opt);

fprintf('model evaluation results:\n');
fprintf('  Social value      : $%.4f\n', social_value);
fprintf('  Profit (P)        : $%.4f\n', P);
fprintf('  Emissions value   : $%.4f\n', emissions_value);
fprintf('  PV revenue        : $%.4f\n', pv_revenue);
fprintf('  Crop revenue      : $%.4f\n', crop_revenue);
fprintf('  Yearly biomass    : %.4f kg/yr\n', yearly_biomass);
fprintf('  Total panels      : %.1f\n', total_panels);
fprintf('  CO2 displaced     : %.4f tons\n', E / 1000);

base_dir = fileparts(mfilename('fullpath'));
if isempty(base_dir), base_dir = pwd; end
results_dir = fullfile(base_dir, 'results', 'SQP');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
csv_file = fullfile(results_dir, sprintf('evaluateSQPDesign_%s.csv', timestamp));

var_names = arrayfun(@(i) sprintf('x_%02d', i), 1:numel(x_opt), 'UniformOutput', false);
T = table('Size', [1, numel(var_names) + 8], ...
    'VariableTypes', repmat({'double'}, 1, numel(var_names) + 8), ...
    'VariableNames', [var_names, {'social_value', 'profit_P', 'emissions_value', ...
    'pv_revenue', 'crop_revenue', 'yearly_biomass_kg_per_yr', ...
    'total_panels', 'co2_displaced_tons'}]);

for i = 1:numel(x_opt)
    T{1, var_names{i}} = x_opt(i);
end
T.social_value = social_value;
T.profit_P = P;
T.emissions_value = emissions_value;
T.pv_revenue = pv_revenue;
T.crop_revenue = crop_revenue;
T.yearly_biomass_kg_per_yr = yearly_biomass;
T.total_panels = total_panels;
T.co2_displaced_tons = E / 1000;

formatSpec = [repmat('%.15g,', 1, width(T) - 1), '%.15g\n'];
fid = fopen(csv_file, 'w');
fprintf(fid, '%s\n', strjoin(T.Properties.VariableNames, ','));
fprintf(fid, formatSpec, table2array(T).');
fclose(fid);
fprintf('  Results CSV       : %s\n', csv_file);
