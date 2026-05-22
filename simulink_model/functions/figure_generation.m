%% Load Data
fixed = load('agrivoltaic_comparative_optimization_data_fixed.mat');
single = load('agrivoltaic_comparative_optimization_data_single.mat');

%% Create Normalized Barchart
% compares single objectives' results

% targets_to_run = {'PROFIT', 'EMISSIONS', 'POWER', 'CROP'};
% results_matrix columns:
% 1 = Emissions
% 2 = Profit
% 4 = Energy
% 6 = Crop Yield


metric_names = {'Profit', 'Energy', 'Emissions', 'Crop Yield'};

plot_data = [
    fixed.results_matrix(:,2), ... % Profit
      single.results_matrix(:,2), ... % Profit
    fixed.results_matrix(:,4), ... % Energy
     single.results_matrix(:,4), ... % Energy
    fixed.results_matrix(:,1), ... % Emissions
    single.results_matrix(:,1), ... % Emissions
    fixed.results_matrix(:,6)  ... % Crop yield
    single.results_matrix(:,6)  ... % Crop yield
];

% Normalize each column by its maximum
plot_data_norm = plot_data ./ max(plot_data, [], 1);

figure('Color','w');
bar(categorical(targets_to_run), plot_data_norm);
ylabel('Normalized Value');
title('Performance of Optimized Designs by Objective');
legend(metric_names, 'Location', 'bestoutside');
grid on;
saveas(gcf, 'graphs/normalized_objective_performance.png');

%% Create Heatmap
% Compares design variables across different designs for each objective

% ------ Single Axis --------- %
design_names = {'Height','Length','Width','Row Spacing','Panel Gap'};
keep_indices = [1, 2, 3, 6, 7];

design_data = x_best_set(:, keep_indices);

% Normalize columns for heatmap
design_norm = design_data ./ max(design_data, [], 1);

figure('Color','w');
h = heatmap(design_names, targets_to_run, design_norm);
h.Title = 'Single-Axis: Normalized Design Variables by Objective';
h.XLabel = 'Design Variable';
h.YLabel = 'Optimization Objective';
h.Colormap = parula;


% ------ Fixed Axis --------- %

design_names = {'Height','Length','Width','Azimuth','Tilt','Row Spacing','Panel Gap'};
keep_indices = 1:7;

design_data = x_best_set(:, keep_indices);
design_data(:,4) = rad2deg(design_data(:,4));
design_data(:,5) = rad2deg(design_data(:,5));

design_norm = design_data ./ max(abs(design_data), [], 1);

figure('Color','w');
h = heatmap(design_names, targets_to_run, design_norm);
h.Title = 'Fixed-Axis: Normalized Design Variables by Objective';
h.XLabel = 'Design Variable';
h.YLabel = 'Optimization Objective';

%% Fixed vs Single Axis Comparison

% Choose objective to compare, for example PROFIT
objective_to_compare = 'PROFIT';

idx_fixed = find(strcmp(fixed.targets_to_run, objective_to_compare));
idx_single = find(strcmp(single.targets_to_run, objective_to_compare));

metrics = {'Profit','Energy','Emissions','Crop Yield'};

fixed_vals = [
    fixed.results_matrix(idx_fixed,2), ...
    fixed.results_matrix(idx_fixed,8), ...
    fixed.results_matrix(idx_fixed,1), ...
    fixed.results_matrix(idx_fixed,6)
];

single_vals = [
    single.results_matrix(idx_single,2), ...
    single.results_matrix(idx_single,8), ...
    single.results_matrix(idx_single,1), ...
    single.results_matrix(idx_single,6)
];

data = [fixed_vals; single_vals];

% Normalize by max of each metric
data_norm = data ./ max(data, [], 1);

figure('Color','w');
bar(categorical({'Fixed Axis','Single Axis'}), data_norm);
ylabel('Normalized value');
title(['Fixed vs Single-Axis Comparison: ', objective_to_compare, '-Optimized Design']);
legend(metrics, 'Location','bestoutside');
grid on;

% Percent Change Table
percent_change = 100 * (single_vals - fixed_vals) ./ fixed_vals;

comparison_table = table(metrics.', fixed_vals.', single_vals.', percent_change.', ...
    'VariableNames', {'Metric','FixedAxis','SingleAxis','PercentChange'});

disp(comparison_table);

writetable(comparison_table, 'fixed_vs_singleaxis_comparison.csv');

%%