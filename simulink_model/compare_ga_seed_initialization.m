function results = compare_ga_seed_initialization(varargin)
% compare_ga_seed_initialization
% Compare fixed-axis and single-axis GA initialization behavior without
% changing the core optimization scripts.
%
% Basic use from simulink_model:
%   results = compare_ga_seed_initialization;
%
% Optional custom seed comparison:
%   results = compare_ga_seed_initialization('CustomSeedFile', 'my_seed.mat');
%
% The custom MAT file can contain one of:
%   x_seed          1 x 103 single-axis vector
%   tracking_angles 4 x 24 tracking angle matrix, radians
%   agriVar         struct with tracking_angles
%
% Optional short GA run:
%   results = compare_ga_seed_initialization('RunGA', true, 'MaxGenerations', 8);

opts = parse_options(varargin{:});

script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end
start_dir = pwd;
cleanup_cd = onCleanup(@() cd(start_dir));
cd(script_dir);
addpath(genpath(script_dir));

% agrivoltaics_variable_definition is a script with clear at top. Running it
% in base keeps this function's local configuration intact.
evalin('base', sprintf('cd(''%s''); agrivoltaics_variable_definition;', escape_single_quotes(script_dir)));

agriParams_base = evalin('base', 'agriParams');
agriVar_base = evalin('base', 'agriVar');
lb_base = evalin('base', 'lb');
ub_base = evalin('base', 'ub');

if exist('agrivoltaics_v1.slx', 'file') == 2
    load_system('agrivoltaics_v1');
    old_init_fcn = get_param('agrivoltaics_v1', 'InitFcn');
    cleanup_model = onCleanup(@() set_param('agrivoltaics_v1', 'InitFcn', old_init_fcn));
    set_param('agrivoltaics_v1', 'InitFcn', '');
end

rng(opts.RandomSeed);

params_fixed = agriParams_base;
params_fixed.tracking_mode = 0;
var_fixed = agriVar_base;
var_fixed.tracking_angles = zeros(4, 24);
x_fixed = agriVarStruct2Array(var_fixed, params_fixed);

params_track = agriParams_base;
params_track.tracking_mode = 1;
var_track = agriVar_base;
var_track.tracking_angles = generate_physics_tracking(params_track, var_track);
var_track.tracking_angles = clamp_tracking(var_track.tracking_angles, params_track.PV_max_tilt);

lb_track = [lb_base(:).', zeros(1, 96)];
ub_track = [ub_base(:).', ones(1, 96) * params_track.PV_max_tilt];
x_previous = agriVarStruct2Array(var_track, params_track);
x_previous = clamp_vector(x_previous, lb_track, ub_track);

x_flat = x_previous;
x_flat(8:end) = 0;

variants = make_variant_template();
variants(end+1) = make_variant('fixed_axis_baseline', params_fixed, lb_base(:).', ub_base(:).', ...
    perturb_layout_population(x_fixed, lb_base(:).', ub_base(:).', opts.PopulationSize, opts.LayoutNoiseFraction));
variants(end+1) = make_variant('single_axis_previous_seed', params_track, lb_track, ub_track, ...
    perturb_tracking_population(x_previous, lb_track, ub_track, opts.PopulationSize, opts.TrackingNoiseFraction, params_track.PV_max_tilt));
variants(end+1) = make_variant('single_axis_flat_seed', params_track, lb_track, ub_track, ...
    perturb_tracking_population(x_flat, lb_track, ub_track, opts.PopulationSize, opts.TrackingNoiseFraction, params_track.PV_max_tilt));
variants(end+1) = make_variant('single_axis_random_angles', params_track, lb_track, ub_track, ...
    random_tracking_population(x_previous, lb_track, ub_track, opts.PopulationSize));

if ~isempty(opts.CustomSeedFile)
    x_custom = load_custom_seed(opts.CustomSeedFile, x_previous, params_track, lb_track, ub_track);
    variants(end+1) = make_variant('single_axis_new_seed', params_track, lb_track, ub_track, ...
        perturb_tracking_population(x_custom, lb_track, ub_track, opts.PopulationSize, opts.TrackingNoiseFraction, params_track.PV_max_tilt));
end

fprintf('\nGA seed initialization comparison\n');
fprintf('  Fixed-axis vector      : %d variables\n', numel(x_fixed));
fprintf('  Single-axis vector     : %d variables = 7 layout + 96 hourly tracking angles\n', numel(x_previous));
fprintf('  Population size/variant: %d\n', opts.PopulationSize);
fprintf('  Run short GA           : %d\n\n', opts.RunGA);

rows = repmat(make_eval_row_template(), 0, 1);
ga_rows = repmat(make_ga_row_template(), 0, 1);

for i = 1:numel(variants)
    v = variants(i);
    fprintf('Evaluating %s (%d candidates)...\n', v.name, size(v.population, 1));
    rows = [rows; evaluate_population(v)]; %#ok<AGROW>

    if opts.RunGA
        fprintf('  Running short GA for %s...\n', v.name);
        ga_rows = [ga_rows; run_short_ga(v, opts.MaxGenerations)]; %#ok<AGROW>
    end
end

T = struct2table(rows);
if isempty(ga_rows)
    Tga = table();
else
    Tga = struct2table(ga_rows);
end

fig_objectives = plot_objective_comparison(T, Tga);
fig_tracking = plot_seed_tracking(variants);

results = struct();
results.initial_population = T;
results.ga_runs = Tga;
results.variants = variants;
results.figures = struct('objective_comparison', fig_objectives, 'seed_tracking', fig_tracking);

if opts.SaveOutputs
    outdir = opts.OutputDir;
    if isempty(outdir)
        outdir = fullfile(script_dir, 'results', ['seed_initialization_' datestr(now, 'yyyymmdd_HHMMSS')]);
    end
    if exist(outdir, 'dir') ~= 7
        mkdir(outdir);
    end
    writetable(T, fullfile(outdir, 'initial_population_metrics.csv'));
    if ~isempty(Tga)
        writetable(Tga, fullfile(outdir, 'short_ga_metrics.csv'));
    end
    saveas(fig_objectives, fullfile(outdir, 'seed_objective_comparison.png'));
    saveas(fig_tracking, fullfile(outdir, 'seed_tracking_angles.png'));
    save(fullfile(outdir, 'seed_initialization_results.mat'), 'results');
    fprintf('\nSaved comparison outputs to %s\n', outdir);
else
    fprintf('\nNo files written. Pass ''SaveOutputs'', true to save CSV/PNG/MAT outputs.\n');
end
end

function opts = parse_options(varargin)
p = inputParser;
p.addParameter('PopulationSize', 8, @(x) isnumeric(x) && isscalar(x) && x >= 2);
p.addParameter('RandomSeed', 101, @(x) isnumeric(x) && isscalar(x));
p.addParameter('TrackingNoiseFraction', 0.10, @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('LayoutNoiseFraction', 0.05, @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('CustomSeedFile', '', @(x) ischar(x) || isstring(x));
p.addParameter('RunGA', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('MaxGenerations', 6, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('SaveOutputs', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('OutputDir', '', @(x) ischar(x) || isstring(x));
p.parse(varargin{:});
opts = p.Results;
opts.CustomSeedFile = char(opts.CustomSeedFile);
opts.OutputDir = char(opts.OutputDir);
opts.RunGA = logical(opts.RunGA);
opts.SaveOutputs = logical(opts.SaveOutputs);
opts.PopulationSize = round(opts.PopulationSize);
opts.MaxGenerations = round(opts.MaxGenerations);
end

function variants = make_variant_template()
variants = struct('name', {}, 'params', {}, 'lb', {}, 'ub', {}, 'population', {});
end

function v = make_variant(name, params, lb, ub, population)
v = struct();
v.name = name;
v.params = params;
v.lb = lb(:).';
v.ub = ub(:).';
v.population = population;
end

function pop = perturb_layout_population(x_seed, lb, ub, pop_size, noise_fraction)
num_vars = numel(x_seed);
pop = zeros(pop_size, num_vars);
pop(1,:) = clamp_vector(x_seed, lb, ub);
for i = 2:pop_size
    noise = noise_fraction * (ub - lb) .* randn(1, num_vars);
    pop(i,:) = clamp_vector(x_seed + noise, lb, ub);
end
end

function pop = perturb_tracking_population(x_seed, lb, ub, pop_size, noise_fraction, max_tilt)
num_vars = numel(x_seed);
pop = zeros(pop_size, num_vars);
pop(1,:) = clamp_vector(x_seed, lb, ub);
for i = 2:pop_size
    candidate = x_seed;
    idx = 8:num_vars;
    candidate(idx) = candidate(idx) + noise_fraction * max_tilt .* randn(1, numel(idx));
    pop(i,:) = clamp_vector(candidate, lb, ub);
end
end

function pop = random_tracking_population(x_seed, lb, ub, pop_size)
num_vars = numel(x_seed);
pop = zeros(pop_size, num_vars);
pop(:, 1:7) = repmat(x_seed(1:7), pop_size, 1);
idx = 8:num_vars;
span = ub(idx) - lb(idx);
pop(:, idx) = repmat(lb(idx), pop_size, 1) + rand(pop_size, numel(idx)) .* repmat(span, pop_size, 1);
for i = 1:pop_size
    pop(i,:) = clamp_vector(pop(i,:), lb, ub);
end
end

function x = clamp_vector(x, lb, ub)
x = max(x(:).', lb(:).');
x = min(x, ub(:).');
end

function theta = clamp_tracking(theta, max_tilt)
theta = max(theta, 0);
theta = min(theta, max_tilt);
end

function x_custom = load_custom_seed(seed_file, x_default, params, lb, ub)
if exist(seed_file, 'file') ~= 2
    error('CustomSeedFile not found: %s', seed_file);
end
S = load(seed_file);
if isfield(S, 'x_seed')
    x_custom = S.x_seed(:).';
elseif isfield(S, 'tracking_angles')
    var_custom = agriVarArray2Struct(x_default, params);
    var_custom.tracking_angles = S.tracking_angles;
    x_custom = agriVarStruct2Array(var_custom, params);
elseif isfield(S, 'agriVar')
    x_custom = agriVarStruct2Array(S.agriVar, params);
else
    error('CustomSeedFile must contain x_seed, tracking_angles, or agriVar.');
end
if numel(x_custom) ~= 103
    error('Custom single-axis seed must have 103 values, got %d.', numel(x_custom));
end
x_custom = clamp_vector(x_custom, lb, ub);
end

function rows = evaluate_population(variant)
n = size(variant.population, 1);
rows = repmat(make_eval_row_template(), n, 1);
for i = 1:n
    x = variant.population(i,:);
    raw = agrivoltaic_wrapper(x, variant.params);
    trk = tracking_matrix(x, variant.params);
    rows(i).variant = variant.name;
    rows(i).member = i;
    rows(i).is_seed = (i == 1);
    rows(i).num_variables = numel(x);
    rows(i).objective_cost = raw(3);
    rows(i).social_value = -raw(3);
    rows(i).profit = raw(2);
    rows(i).emissions_value = 190 * (raw(1) / 1000);
    rows(i).pv_revenue = raw(4);
    rows(i).crop_revenue = raw(5);
    rows(i).yearly_biomass = raw(6);
    rows(i).total_panels = raw(7);
    rows(i).mean_tracking_deg = mean(rad2deg(trk(:)), 'omitnan');
    rows(i).max_tracking_deg = max(rad2deg(trk(:)));
end
end

function row = run_short_ga(variant, max_generations)
num_vars = size(variant.population, 2);
opts = optimoptions('ga', ...
    'PopulationSize', size(variant.population, 1), ...
    'MaxGenerations', max_generations, ...
    'FunctionTolerance', 1e-4, ...
    'Display', 'off', ...
    'InitialPopulationMatrix', variant.population);

t0 = tic;
[x_best, f_best, exitflag, output] = ga(@(x) agrivoltaic_social_cost_of_carbon_wrapper(x, variant.params), ...
    num_vars, [], [], [], [], variant.lb, variant.ub, [], opts);
elapsed = toc(t0);
trk = tracking_matrix(x_best, variant.params);

row = make_ga_row_template();
row.variant = variant.name;
row.max_generations = max_generations;
row.objective_cost = f_best;
row.social_value = -f_best;
row.exitflag = exitflag;
row.function_count = output.funccount;
row.elapsed_seconds = elapsed;
row.mean_tracking_deg = mean(rad2deg(trk(:)), 'omitnan');
row.max_tracking_deg = max(rad2deg(trk(:)));
end

function trk = tracking_matrix(x, params)
if isfield(params, 'tracking_mode') && params.tracking_mode == 1 && numel(x) > 7
    trk = reshape(x(8:end), [24, 4]).';
else
    trk = zeros(4, 24);
end
end

function row = make_eval_row_template()
row = struct();
row.variant = '';
row.member = NaN;
row.is_seed = false;
row.num_variables = NaN;
row.objective_cost = NaN;
row.social_value = NaN;
row.profit = NaN;
row.emissions_value = NaN;
row.pv_revenue = NaN;
row.crop_revenue = NaN;
row.yearly_biomass = NaN;
row.total_panels = NaN;
row.mean_tracking_deg = NaN;
row.max_tracking_deg = NaN;
end

function row = make_ga_row_template()
row = struct();
row.variant = '';
row.max_generations = NaN;
row.objective_cost = NaN;
row.social_value = NaN;
row.exitflag = NaN;
row.function_count = NaN;
row.elapsed_seconds = NaN;
row.mean_tracking_deg = NaN;
row.max_tracking_deg = NaN;
end

function fig = plot_objective_comparison(T, Tga)
fig = figure('Name', 'GA Seed Objective Comparison', 'NumberTitle', 'off', ...
    'Color', 'w', 'Position', [80 80 1250 760]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

plot_metric_tile(T, Tga, 'social_value', 1e6, 'Total social value (M$)', 'Social value by seed variant');
plot_metric_tile(T, Tga, 'profit', 1e6, 'Profit (M$)', 'Profit by seed variant');
plot_metric_tile(T, Tga, 'emissions_value', 1e6, 'Emissions value (M$)', 'Emissions value by seed variant');
plot_metric_tile(T, Tga, 'yearly_biomass', 1000, 'Yearly biomass (thousand kg/yr)', 'Biomass by seed variant');
end

function plot_metric_tile(T, Tga, metric, scale, ylabel_text, title_text)
nexttile;
cats = categorical(T.variant);
boxchart(cats, T.(metric) ./ scale);
hold on;
seed_rows = T.is_seed;
scatter(categorical(T.variant(seed_rows)), T.(metric)(seed_rows) ./ scale, 55, 'k', 'filled', ...
    'DisplayName', 'seed member');
if ~isempty(Tga) && ismember(metric, Tga.Properties.VariableNames)
    scatter(categorical(Tga.variant), Tga.(metric) ./ scale, 80, 'd', 'filled', ...
        'DisplayName', 'short GA final');
end
grid on;
ylabel(ylabel_text);
title(title_text);
xtickangle(25);
legend('Location', 'best');
hold off;
end

function fig = plot_seed_tracking(variants)
fig = figure('Name', 'Seed Tracking Angles', 'NumberTitle', 'off', ...
    'Color', 'w', 'Position', [120 120 1250 760]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
season_names = {'Spring', 'Summer', 'Fall', 'Winter'};
hours = 0:23;
colors = lines(numel(variants));

for s = 1:4
    nexttile;
    hold on;
    for i = 1:numel(variants)
        x_seed = variants(i).population(1,:);
        if isfield(variants(i).params, 'tracking_mode') && variants(i).params.tracking_mode == 1
            trk = tracking_matrix(x_seed, variants(i).params);
            y = rad2deg(trk(s,:));
        else
            y = zeros(1, 24);
        end
        plot(hours, y, '-o', 'LineWidth', 1.3, 'MarkerSize', 3, ...
            'Color', colors(i,:), 'DisplayName', variants(i).name);
    end
    yline(0, 'k--', 'HandleVisibility', 'off');
    grid on;
    xlabel('Hour');
    ylabel('Tracking angle (deg)');
    title([season_names{s} ' seed tracking schedule']);
    xlim([0 23]);
    hold off;
end
legend('Location', 'southoutside', 'Orientation', 'horizontal');
end

function escaped = escape_single_quotes(txt)
escaped = strrep(txt, '''', '''''');
end
