% Clear and Setup
clear;
clc;

addpath(genpath(pwd));

agrivoltaics_variable_definition;

GA_SELECTOR = 4;

%1 for basic GA
if GA_SELECTOR == 1
    rng(1);
    options = optimoptions('ga', 'PopulationSize', 50, 'MaxGenerations', 40, ...
        'FunctionTolerance', 1e-4,'MaxStallGenerations', 10,'Display', ...
        'iter','PlotFcn', @gaplotbestf);
elseif GA_SELECTOR == 2
    rng(2);
    options = optimoptions('ga', 'PopulationSize', 50, 'MaxGenerations', 100, ...
        'FunctionTolerance', 1e-4,'MaxStallGenerations', 10,'Display', ...
        'iter','PlotFcn', @gaplotbestf);
elseif GA_SELECTOR == 3
    rng(3);
    options = optimoptions('ga', 'PopulationSize', 200, 'MaxGenerations', 100, ...
        'FunctionTolerance', 1e-4,'MaxStallGenerations', 10,'Display', ...
        'iter','PlotFcn', @gaplotbestf);
elseif GA_SELECTOR == 4
    rng(4);
    options = optimoptions('ga', 'PopulationSize', 1800, 'MaxGenerations', 100, ...
        'FunctionTolerance', 1e-4,'MaxStallGenerations', 10,'Display', ...
        'iter','PlotFcn', @gaplotbestf);
end


% Set Up GA
A = []; B = []; Aeq = []; Beq = [];
nlcon = [];

% Run GA
tic;
[ga_solve,fval,exitflag,output,population,scores] = ...
    ga(@agrivoltaic_social_cost_of_carbon_wrapper, ...
    7, A, B, Aeq, Beq, lb, ub, ...
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