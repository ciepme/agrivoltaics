% Clear and Setup
clear;
clc;

addpath(genpath(pwd));

agrivoltaics_variable_definition;

GA_SELECTOR = 1;
% Tell GA exactly how many variables we are optimizing (7 or 103, dependent
% on fixed or single-axis)
num_vars = length(lb);
%create a smart guess to seed GA population for if tracking mode
x0_base = agriVarStruct2Array(agriVar, agriParams);
%adds an InitialPopulationMatrix for better initial guess

x0 = agriVarStruct2Array(agriVar, agriParams);
pop_size = 50;

% build a better initial population
% Make sure consistent
num_vars = length(lb);
pop = zeros(pop_size, num_vars);
% 1First member = physics-based smart guess
pop(1,:) = x0;

% population = small perturbations
for i = 2:pop_size
    candidate = x0;

    if agriParams.tracking_mode == 1
        % Only perturb tracking angles (more stable)
        idx = 8:num_vars;

        noise = 0.1 * agriParams.PV_max_tilt * randn(size(idx));
        candidate(idx) = candidate(idx) + noise;

        % Clamp
        candidate(idx) = max(candidate(idx), lb(idx));
        candidate(idx) = min(candidate(idx), ub(idx));
    else
        % Fixed-axis: perturb all vars slightly
        noise = 0.05 * (ub - lb) .* randn(1, num_vars);
        candidate = candidate + noise;

        candidate = max(candidate, lb);
        candidate = min(candidate, ub);
    end

    pop(i,:) = candidate;
end
%1 for basic GA
if GA_SELECTOR == 1
    rng(1);
    options = optimoptions('ga', 'PopulationSize', pop_size, 'MaxGenerations', 40, ...
        'FunctionTolerance', 1e-4,'MaxStallGenerations', 10,'Display', ...
        'iter','PlotFcn', @gaplotbestf, 'InitialPopulationMatrix', x0);
elseif GA_SELECTOR == 2
    rng(2);
    options = optimoptions('ga', 'PopulationSize', pop_size, 'MaxGenerations', 100, ...
        'FunctionTolerance', 1e-4,'MaxStallGenerations', 10,'Display', ...
        'iter','PlotFcn', @gaplotbestf, 'InitialPopulationMatrix', x0);
elseif GA_SELECTOR == 3
    rng(3);
    options = optimoptions('ga', 'PopulationSize', pop_size, 'MaxGenerations', 100, ...
        'FunctionTolerance', 1e-4,'MaxStallGenerations', 10,'Display', ...
        'iter','PlotFcn', @gaplotbestf, 'InitialPopulationMatrix', x0);
elseif GA_SELECTOR == 4
    rng(4);
    options = optimoptions('ga', 'PopulationSize', pop_size, 'MaxGenerations', 100, ...
        'FunctionTolerance', 1e-4,'MaxStallGenerations', 10,'Display', ...
        'iter','PlotFcn', @gaplotbestf, 'InitialPopulationMatrix', x0);
end


% Set Up GA
A = []; B = []; Aeq = []; Beq = [];
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
fprintf('Max Value: $%.2f\n', -fval);
fprintf('Panel Height: %.2f m\n', ga_solve(1));
fprintf('Panel Length: %.2f m\n', ga_solve(2));
fprintf('Panel Width: %.2f m\n', ga_solve(3));
fprintf('Azimuth %.2f rad (%.1f deg)\n', ga_solve(4), rad2deg(ga_solve(4)));
fprintf('Tilt: %.2f rad (%.1f deg)\n', ga_solve(5), rad2deg(ga_solve(5)));
fprintf('Row Spacing: %.2f m\n', ga_solve(6));
fprintf('Panel Gap (x_p): %.2f m\n', ga_solve(7))