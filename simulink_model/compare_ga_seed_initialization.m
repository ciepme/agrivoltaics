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

%%
% OLD BOUNDS 
lb_track_old = [lb_base(:).', zeros(1, 96)];
ub_track_old = [ub_base(:).', ones(1, 96) * params_track.PV_max_tilt];

% NEW BOUNDS adds in logic so that no tracking at night
tracking_lb_new = zeros(1, 96); 
tracking_ub_new = zeros(1, 96); 
seasons = {'spring', 'summer', 'fall', 'winter'};
for s = 1:4
    is_daytime = params_track.weather.(seasons{s}).beta_s > 0;
    start_idx = (s-1)*24 + 1;
    end_idx = s*24;
    tracking_lb_new(start_idx:end_idx) = -params_track.PV_max_tilt .* is_daytime;
    tracking_ub_new(start_idx:end_idx) =  params_track.PV_max_tilt .* is_daytime;
end
lb_track_new = [lb_base(:).', tracking_lb_new];
ub_track_new = [ub_base(:).', tracking_ub_new];

% creation of previous seeds using respective bounds 
x_previous = agriVarStruct2Array(var_track, params_track);

% create previous seeds using respective bounds 
x_previous_old = clamp_vector(x_previous, lb_track_old, ub_track_old);
x_previous_new = clamp_vector(x_previous, lb_track_new, ub_track_new);

% use flat baseline
x_flat = x_previous_old;
x_flat(8:end) = 0;
variants = make_variant_template();
variants(end+1) = make_variant('fixedAxisBaseline', params_fixed, lb_base(:).', ub_base(:).', ...
    perturb_layout_population(x_fixed, lb_base(:).', ub_base(:).', opts.PopulationSize, opts.LayoutNoiseFraction));


% 100% random no seeding
variants(end+1) = make_variant('SingleAxisNOSeed', params_track, lb_track_new, ub_track_new, ...
    random_tracking_population(x_previous_new, lb_track_new, ub_track_new, opts.PopulationSize));

% perturbed seed with  no smoothing (Old Population Rules)
% (Adds jagged noise directly on top of the ideal physics seed)
variants(end+1) = make_variant('SingleAxisPerturbedSeed', params_track, lb_track_new, ub_track_new, ...
    perturb_tracking_population(x_previous_new, lb_track_new, ub_track_new, opts.PopulationSize, opts.TrackingNoiseFraction, params_track.PV_max_tilt));

% (Uses the ideal physics seed for Member 1. Members 2+ completely overwrite 
% the tracking angles with purely random, Gaussian-smoothed curves)
variants(end+1) = make_variant('SingleAxisRandomSmoothed', params_track, lb_track_new, ub_track_new, ...
    perturb_tracking_population_smooth(x_previous_new, lb_track_new, ub_track_new, opts.PopulationSize, opts.TrackingNoiseFraction, params_track.PV_max_tilt));

if ~isempty(opts.CustomSeedFile)
    x_custom = load_custom_seed(opts.CustomSeedFile, x_previous, params_track, lb_track_new, ub_track_new);
    variants(end+1) = make_variant('singleAxisNewSeed', params_track, lb_track_new, ub_track_new, ...
        perturb_tracking_population_smooth(x_custom, lb_track_new, ub_track_new, opts.PopulationSize, opts.TrackingNoiseFraction, params_track.PV_max_tilt));
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

%new figure with just single axis optimization for weighted sum
if ~isempty(Tga)
    fig_best_single_axis = plot_best_single_axis_tracking(variants, Tga);
end

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
        % NEW: Save the exclusive Single-Axis plot
        saveas(fig_best_single_axis, fullfile(outdir, 'final_single_axis_tracking.png'));
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
        candidate = x_seed;
        
        % only perturbing variable 5 (the Fixed Tilt angle)
        noise = noise_fraction * (ub(5) - lb(5)) * randn(1);
        candidate(5) = candidate(5) + noise;
        
        pop(i,:) = clamp_vector(candidate, lb, ub);
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

%added in Gaussian smoothing population with one initial guess and the rest
%random
function pop = perturb_tracking_population_smooth(x_seed, lb, ub, pop_size, noise_fraction, max_tilt)
    num_vars = numel(x_seed);
    pop = zeros(pop_size, num_vars);
    
    % adding one member of the perfect, physics-based seed
    pop(1,:) = clamp_vector(x_seed, lb, ub);
    
    % rest of the members are random but smooth
    for i = 2:pop_size
        candidate = x_seed; % Keep layout variables (1-7) the same as base
        idx = 8:num_vars;   % Target the 96 tracking variables
        
        % Generate random angles between the lower and upper bounds
        span = ub(idx) - lb(idx);
        candidate(idx) = lb(idx) + rand(1, numel(idx)) .* span;  
        
        % Smooth the random curve so it's physically possible
        candidate(idx) = smoothdata(candidate(idx), 'gaussian', 3);
        
        % Clamp against lb/ub (This forces the night hours back to 0)
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
    theta = max(theta, -max_tilt); % allow negative tracking
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

A = []; B = [];
% only apply the 15-degree slew rate to our new variant
if contains(variant.name, 'NEW_Rules') && isfield(variant.params, 'tracking_mode') && variant.params.tracking_mode == 1
    max_slew = deg2rad(15);
    total_cons = 4 * 23 * 2;
    A = zeros(total_cons, num_vars);
    B = ones(total_cons, 1) * max_slew;
    row = 1;
    for s = 1:4
        offset = 7 + (s-1)*24;
        for h = 1:23
            v1 = offset + h; v2 = offset + h + 1;
            A(row, v1) = -1; A(row, v2) = 1;  row = row + 1;
            A(row, v1) = 1;  A(row, v2) = -1; row = row + 1;
        end
    end
end


% update the ga() call to use A and B
opts.InitialPopulationMatrix = variant.population; % just to ensure it uses the right pop

delta = 0.5; % Balance social profit and crop growth equally

t0 = tic; % start timer before GA call
[x_best, f_best, exitflag, output] = ga(@(x) agrivoltaic_weighted_sum_social_benefit_wrapper(x, variant.params, delta), ...
    num_vars, A, B, [], [], variant.lb, variant.ub, [], opts);

elapsed = toc(t0); % stop timer

%evaluate x_best to get the broken-out metrics
raw = agrivoltaic_wrapper(x_best, variant.params);
trk = tracking_matrix(x_best, variant.params);

row = make_ga_row_template();
row.x_final(1:numel(x_best)) = x_best; 

row.variant = variant.name;
row.max_generations = max_generations;
row.objective_cost = f_best;

% specific metrics to plot the diamond
row.social_value = -raw(3); 
row.profit = raw(2);
row.emissions_value = 190 * (raw(1) / 1000);
row.pv_revenue = raw(4);
row.crop_revenue = raw(5);
row.yearly_biomass = raw(6);
row.total_panels = raw(7);

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
    
    % added in new metrics
    row.social_value = NaN;
    row.profit = NaN;
    row.emissions_value = NaN;
    row.pv_revenue = NaN;
    row.crop_revenue = NaN;
    row.yearly_biomass = NaN;
    row.total_panels = NaN;
    
    row.exitflag = NaN;
    row.function_count = NaN;
    row.elapsed_seconds = NaN;
    row.mean_tracking_deg = NaN;
    row.max_tracking_deg = NaN;
    row.x_final = zeros(1, 103);
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

function fig = plot_best_single_axis_tracking(variants, Tga)
    fig = figure('Name', 'Final Single-Axis Competitors', 'NumberTitle', 'off', ...
        'Color', 'w', 'Position', [150 150 1250 760]);
    tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    season_names = {'Spring', 'Summer', 'Fall', 'Winter'};
    hours = 0:23;
    colors = lines(numel(variants));
    
    for s = 1:4
        nexttile;
        hold on;
        for i = 1:numel(variants)
            % Skip the fixed axis variant entirely
            if contains(variants(i).name, 'Fixed')
                continue; 
            end
            
            % Only plot if the GA actually ran and generated x_final
            if ~isempty(Tga) && i <= height(Tga)
                x_plot = Tga.x_final(i, 1:numel(variants(i).lb)); 
                
                trk = tracking_matrix(x_plot, variants(i).params);
                y = rad2deg(trk(s,:));
                
                % Clean up the variant name for the legend
                display_name = strrep(variants(i).name, '_', ' ');
                
                plot(hours, y, '-o', 'LineWidth', 1.8, 'MarkerSize', 4, ...
                    'Color', colors(i,:), 'DisplayName', display_name);
            end
        end
        yline(0, 'k--', 'HandleVisibility', 'off');
        grid on;
        xlabel('Hour');
        ylabel('Tracking angle (deg)');
        title([season_names{s} ' Final Optimized Tracking']);
        xlim([0 23]);
        hold off;
    end
    legend('Location', 'southoutside', 'Orientation', 'horizontal');
end
