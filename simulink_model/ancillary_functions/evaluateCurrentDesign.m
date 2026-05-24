clear;
clc;
% for running the model once at a specific design vector (from agri var def) and getting clear results without having to open up simulink 
agrivoltaics_variable_definition;
x = agriVarStruct2Array(agriVar, agriParams);
results = agrivoltaic_wrapper(x, agriParams);

E = results(1);
P = results(2);
social_cost = results(3);
pv_revenue = results(4);
crop_revenue = results(5);
yearly_biomass = results(6);
total_panels = results(7);
social_value = -social_cost;
emissions_value = 190 * (E / 1000);

fprintf('current design vector from agrivoltaics_variable_definition.m:\n');
format long g;
disp(x);

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
results_dir = fullfile(base_dir, 'results');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
csv_file = fullfile(results_dir, sprintf('evaluateCurrentDesign_%s.csv', timestamp));

var_names = arrayfun(@(i) sprintf('x_%02d', i), 1:numel(x), 'UniformOutput', false);
T = table('Size', [1, numel(var_names) + 8], ...
    'VariableTypes', repmat({'double'}, 1, numel(var_names) + 8), ...
    'VariableNames', [var_names, {'social_value', 'profit_P', 'emissions_value', ...
    'pv_revenue', 'crop_revenue', 'yearly_biomass_kg_per_yr', ...
    'total_panels', 'co2_displaced_tons'}]);

for i = 1:numel(x)
    T{1, var_names{i}} = x(i);
end
T.social_value = social_value;
T.profit_P = P;
T.emissions_value = emissions_value;
T.pv_revenue = pv_revenue;
T.crop_revenue = crop_revenue;
T.yearly_biomass_kg_per_yr = yearly_biomass;
T.total_panels = total_panels;
T.co2_displaced_tons = E / 1000;

writetable(T, csv_file);
formatSpec = [repmat('%.15g,', 1, width(T) - 1), '%.15g\n'];
fid = fopen(csv_file, 'w');
fprintf(fid, '%s\n', strjoin(T.Properties.VariableNames, ','));
fprintf(fid, formatSpec, table2array(T).');
fclose(fid);
fprintf('  Results CSV       : %s\n', csv_file);
