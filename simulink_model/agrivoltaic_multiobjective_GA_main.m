%% Clear and Setup
clear;
clc;
close;

addpath(genpath(pwd));

agrivoltaics_variable_definition;

%User define statements
GA_SELECTOR = 1;
lambda_fidelity = 0.5;

lambda = 0:lambda_fidelity:1;

%GA hyperparameter settings
pop_size = 30;
stall = 1;

% Tell GA exactly how many variables we are optimizing (7 or 103, dependent
%create a smart guess to seed GA population for if tracking mode
x0_base = agriVarStruct2Array(agriVar, agriParams);

% on fixed or single-axis)
num_vars = length(lb);
%adds an InitialPopulationMatrix for better initial guess

x0 = agriVarStruct2Array(agriVar, agriParams);

%% build a better initial population
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
        'FunctionTolerance', 1e-4,'MaxStallGenerations', stall,'Display', ...
        'iter','PlotFcn', @gaplotbestf, 'InitialPopulationMatrix', x0);
elseif GA_SELECTOR == 2
    rng(2);
    options = optimoptions('ga', 'PopulationSize', pop_size, 'MaxGenerations', 100, ...
        'FunctionTolerance', 1e-4,'MaxStallGenerations', stall,'Display', ...
        'iter','PlotFcn', @gaplotbestf, 'InitialPopulationMatrix', x0);
elseif GA_SELECTOR == 3
    rng(3);
    options = optimoptions('ga', 'PopulationSize', pop_size, 'MaxGenerations', 100, ...
        'FunctionTolerance', 1e-4,'MaxStallGenerations', stall,'Display', ...
        'iter','PlotFcn', @gaplotbestf, 'InitialPopulationMatrix', x0);
elseif GA_SELECTOR == 4
    rng(4);
    options = optimoptions('ga', 'PopulationSize', pop_size, 'MaxGenerations', 100, ...
        'FunctionTolerance', 1e-4,'MaxStallGenerations', stall,'Display', ...
        'iter','PlotFcn', @gaplotbestf, 'InitialPopulationMatrix', x0);
end

%% Set Up GA
A = []; B = []; Aeq = []; Beq = [];
nlcon = [];

ga_final_pop_set = [];
ga_set = ones(1,7);

% Run GA
for i = 1:length(lambda)
    tic;
    [ga_solve,fval,exitflag,output,population,scores] = ...
        ga(@(x) agrivoltaic_weighted_sum_wrapper(x, agriParams, lambda(i)), ...
        num_vars, A, B, Aeq, Beq, lb, ub, ...
        nlcon, options);
    time_taken = toc;

    fprintf('\n');
    fprintf('lambda of: %.2f\n', lambda(i));
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
    fprintf('\n');
    ga_set(i,:) = ga_solve;
    ga_final_pop_set = [ga_final_pop_set, population];
end

save("agrivoltaic_multiobjective_GA_main_data.mat")

%% Post GA Analysis
mil = 1e6;
E_set = ones(length(lambda),1);
P_set = ones(length(lambda),1);

%find values from GA results
for i = 1:length(lambda)
    wrapper_results = agrivoltaic_wrapper(ga_set, agriParams);
    E_now = wrapper_results(1);
    P_now = wrapper_results(2);
    E_set(i) = E_now./mil;
    P_set(i) = P_now/mil;
end

%plot Pareto
front_size = 300;
utopia_size = 500;
fig1 = figure;
theme(fig1,"light");
scatter(P_set, E_set, front_size, 'yellow', 'filled');
hold on;
title("Pareto Front");
xlabel("Profit ($M)");
ylabel("Emission Reduction (kt CO2e)");
ylim([0 4]);
xlim([0 4]);
scatter(3.5, 3.5, utopia_size, 'cyan', 'filled');
legend("Weighted Sum GA Output", "Utopia Point", 'Location','southwest');
saveas(fig1,'graphs/pareto.png');