% Clear and Setup
clear;
clc;

addpath(genpath(pwd));

agrivoltaics_variable_definition;

GA_SELECTOR = 1;
save_name = sprintf('agrivoltaics_GA_main_data_w_selector_%d.mat', GA_SELECTOR);
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

% 1. first member = physics-based smart guess (from
% generate_physics_tracking)
pop(1,:) = x0;

% 2. Population = small, smooth perturbations (Random Seeds)
for i = 2:pop_size
    candidate = x0;
    
if agriParams.tracking_mode == 1
        idx = 8:num_vars;
        
        % random generation of population
        span = ub(idx) - lb(idx);
        random_angles = lb(idx) + rand(1, length(idx)) .* span;
        
        % smoothing so it doesn't break 15 degree slew limit
        smoothed_angles = smoothdata(random_angles, 'gaussian', 5);
        
        % clamping to set night time hours to 0 degree
        candidate(idx) = max(smoothed_angles, lb(idx));
        candidate(idx) = min(candidate(idx), ub(idx));
        
    else
        % Fixed-axis pure random
        span = ub(1:num_vars) - lb(1:num_vars);
        candidate(1:num_vars) = lb(1:num_vars) + rand(1, num_vars) .* span;
    end
    
    pop(i,:) = candidate;
end

%1 for basic GA
if GA_SELECTOR == 1
    rng(1);
    options = optimoptions('ga', 'PopulationSize', pop_size, 'MaxGenerations', 40, ...
        'FunctionTolerance', 1e-4,'MaxStallGenerations', 10,'Display', ...
        'iter','PlotFcn', @gaplotbestf, 'InitialPopulationMatrix', pop);
elseif GA_SELECTOR == 2
    rng(2);
    options = optimoptions('ga', 'PopulationSize', pop_size+20, 'MaxGenerations', 100, ...
        'FunctionTolerance', 1e-4,'MaxStallGenerations', 10,'Display', ...
        'iter','PlotFcn', @gaplotbestf); % random initial point generation using built in GA intitial point alg
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


%% Set Up GA Constraints
A = []; B = []; Aeq = []; Beq = [];
nlcon = [];

if agriParams.tracking_mode == 1
    max_slew_per_hour = deg2rad(15); % Max 15 degree rotation per hour
    
    % We have 4 seasons * 24 hours = 96 tracking variables.
    % In our flat array, these are indices 8 through 103.
    num_steps = 23; % 23 hour-to-hour transitions per day
    num_constraints_per_season = num_steps * 2; % 2 rules per step (positive and negative limit)
    total_constraints = 4 * num_constraints_per_season;
    
    A = zeros(total_constraints, num_vars);
    B = ones(total_constraints, 1) * max_slew_per_hour;
    
    row = 1;
    for s = 1:4
        offset = 7 + (s-1)*24; % Offset to the start of this season's variables
        for h = 1:23
            v1 = offset + h;
            v2 = offset + h + 1;
            
            % Constraint 1: x(h+1) - x(h) <= max_slew
            A(row, v1) = -1;
            A(row, v2) = 1;
            row = row + 1;
            
            % Constraint 2: x(h) - x(h+1) <= max_slew (prevents negative jumps)
            A(row, v1) = 1;
            A(row, v2) = -1;
            row = row + 1;
        end
    end
end

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
fprintf('\n');
print_value_breakdown(ga_solve, agriParams);

save(save_name);
fprintf('Saved GA results to %s\n', save_name);
