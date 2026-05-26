% Clear and Setup
clear;
clc;
bdclose('all');

addpath(genpath(pwd));

agrivoltaics_variable_definition;

GA_SELECTOR = 4;
USE_PARALLEL_PROCESSING = false;
save_name = sprintf('agrivoltaics_GA_main_data_w_selector_%d.mat', GA_SELECTOR);
% Tell GA exactly how many variables we are optimizing (7 or 103, dependent
% on fixed or single-axis)
num_vars = length(lb);
%create a smart guess to seed GA population for if tracking mode
x0_base = agriVarStruct2Array(agriVar, agriParams);
%adds an InitialPopulationMatrix for better initial guess

pop_size = 20;
pop = build_initial_population(agriVar, agriParams, lb, ub, pop_size, 4);
setup_parallel_pool(USE_PARALLEL_PROCESSING);

%1 for basic GA
if GA_SELECTOR == 1
    rng(1);
    options = optimoptions('ga', 'PopulationSize', pop_size, 'MaxGenerations', 40, ...
        'FunctionTolerance', 1e-4,'MaxStallGenerations', 10,'Display', ...
        'iter','PlotFcn', @gaplotbestf, 'InitialPopulationMatrix', pop, ...
        'UseParallel', USE_PARALLEL_PROCESSING);
elseif GA_SELECTOR == 2
    rng(2);
    options = optimoptions('ga', 'PopulationSize', pop_size+20, 'MaxGenerations', 100, ...
        'FunctionTolerance', 1e-4,'MaxStallGenerations', 10,'Display', ...
        'iter','PlotFcn', @gaplotbestf, ...
        'UseParallel', USE_PARALLEL_PROCESSING); % random initial point generation using built in GA intitial point alg
elseif GA_SELECTOR == 3
    rng(3);
    options = optimoptions('ga', 'PopulationSize', pop_size, 'MaxGenerations', 100, ...
        'FunctionTolerance', 1e-4,'MaxStallGenerations', 10,'Display', ...
        'iter','PlotFcn', @gaplotbestf, 'InitialPopulationMatrix', x0_base, ...
        'UseParallel', USE_PARALLEL_PROCESSING);
elseif GA_SELECTOR == 4
    rng(4);
    options = optimoptions('ga', ...
    'PopulationSize', pop_size, ...
    'MaxGenerations', 5, ...
    'FunctionTolerance', 1e-4, ...
    'MaxStallGenerations', 3, ...
    'Display', 'iter', ...
    'PlotFcn', @gaplotbestf, ...
    'InitialPopulationMatrix', pop, ...
    'UseParallel', USE_PARALLEL_PROCESSING);
end


%% Set Up GA Constraints
[A, B, Aeq, Beq] = build_tracking_slew_constraints(num_vars, agriParams);
nlcon = [];

% Run GA
tic;
[ga_solve,fval,exitflag,output,population,scores] = ...
    ga(@(x) agrivoltaic_social_cost_of_carbon_wrapper(x, agriParams), ...
    num_vars, A, B, Aeq, Beq, lb, ub, ...
    nlcon, options);
time_taken = toc;

fprintf('Exit Flag: %d\n', exitflag);
fprintf('Time Taken: %.2f seconds\n', time_taken);
fprintf('Function Calls: %.0f\n', output.funccount);
fprintf('Parallel Processing: %d\n', USE_PARALLEL_PROCESSING);
fprintf('Max Value: $%.2f\n', -fval);
ga_var = agriVarArray2Struct(ga_solve, agriParams);
print_design_summary(ga_var, agriParams);
fprintf('\n');
print_value_breakdown(ga_solve, agriParams);

save(save_name);
fprintf('Saved GA results to %s\n', save_name);
