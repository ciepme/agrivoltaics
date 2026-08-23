%% Clear and Setup

clc; clear; close all; % Clear command window, workspace, and close figures

%% Helper function
function fitness = moga_objective_wrapper(x, params)
    % Run the standard wrapper to get all outputs
    % raw = [Emissions(1), Profit(2), SocialCost(3), PVRev(4), CropRev(5), Biomass(6), Panels(7), Energy(8), Capex(9)]
    
    raw = agrivoltaic_wrapper(x, params); 
    
    profit = raw(2);

    negative_emission_reduction = -1 * raw(1);
    negative_fiscal_profit      = -1 * raw(2);
    negative_yearly_biomass     = -1 * raw(6);
    positive_capex              = raw(9);
    
    if params.persona == 2 % Janice's code
        fitness = [negative_emission_reduction, negative_fiscal_profit, negative_yearly_biomass, positive_capex];
    elseif params.persona == 1 %Fred's 
        fitness = [negative_emission_reduction, negative_fiscal_profit, negative_yearly_biomass];
    else 
        fitness = [negative_emission_reduction, negative_fiscal_profit, negative_yearly_biomass];
    end
    %% profit constraint, must be positive profit
    if profit < 0
        % If profit is negative, apply a massive base penalty (100 million) 
        % PLUS a scaling factor so the GA knows which "bad" designs are closer to breaking even
        profit_penalty = 1e8 + (abs(profit) * 10);
        
        % Add the penalty to all objectives to ruin the fitness score
        fitness = fitness + profit_penalty;
    end
    %% Fred's machinery width constraints
    if params.persona == 1 % Fred's code
        row_gap = x(6); 
        l_p = x(2);
        machinery_width = params.machinery_width;
        buffer_distance = params.harvest_clearance;
        
        % add penalty for not matching plantable width to machinery width
        nearest_multiple = round((row_gap - 2*buffer_distance - l_p)/machinery_width) * machinery_width;
        mismatch = abs(row_gap - nearest_multiple);
        
        % apply a heavy penalty if it is not perfectly aligned
        penalty_weight = 1e5; 
        penalty = penalty_weight * (mismatch^2);
        
        % Add the penalty to all base objectives to lower the fitness score
        fitness = fitness + penalty;
     end
end

%% Set Up GA

addpath(genpath(pwd));
agrivoltaics_variable_definition;

%User define statements
USE_PARALLEL_PROCESSING = true;

%GA hyperparameter settings
pop_size = 100;
max_gen = 150;
stall = 10;
mode = agriParams.tracking_mode;
persona = agriParams.persona;

moniker = "pop" + pop_size + "_max_gen" + max_gen + "_mode" + mode + "_persona" + persona ;

% Tell GA exactly how many variables we are optimizing (7 or 103, dependent
%create a smart guess to seed GA population for if tracking mode
x0_base = agriVarStruct2Array(agriVar, agriParams);
num_vars = length(lb);
max_slew = agriParams.max_slew_per_hour;

%% Build feasible initial population

rng(1);

% Build x0
x0 = agriVarStruct2Array(agriVar, agriParams);

% Force x0 to obey bounds
x0 = max(x0(:).', lb(:).');
x0 = min(x0(:).', ub(:).');

% Repair x0 for single-axis slew constraints
if agriParams.tracking_mode == 1
    for s = 1:4
        start_idx = 7 + (s-1)*24 + 1;
        end_idx   = 7 + s*24;

        curve = x0(start_idx:end_idx);

        curve = enforce_slew_curve( ...
            curve, ...
            lb(start_idx:end_idx), ...
            ub(start_idx:end_idx), ...
            max_slew);

        x0(start_idx:end_idx) = curve;
    end
end

% Build initial population
pop = zeros(pop_size, num_vars);

% First member is physics/default design
pop(1,:) = x0;

for i = 2:pop_size

    candidate = x0;

    if agriParams.tracking_mode == 1

        idx = 8:num_vars;

        % Random tracking angles within bounds
        span = ub(idx) - lb(idx);
        random_angles = lb(idx) + rand(1, length(idx)) .* span;

        % Smooth random curve
        smoothed_angles = smoothdata(random_angles, 'gaussian', 5);

        % Clamp
        candidate(idx) = max(smoothed_angles, lb(idx));
        candidate(idx) = min(candidate(idx), ub(idx));

        % Enforce slew season-by-season
        for s = 1:4
            start_idx = 7 + (s-1)*24 + 1;
            end_idx   = 7 + s*24;

            curve = candidate(start_idx:end_idx);

            curve = enforce_slew_curve( ...
                curve, ...
                lb(start_idx:end_idx), ...
                ub(start_idx:end_idx), ...
                max_slew);

            candidate(start_idx:end_idx) = curve;
        end

    else

        % Fixed-axis case
        span = ub(1:num_vars) - lb(1:num_vars);
        candidate(1:num_vars) = lb(1:num_vars) + rand(1, num_vars) .* span;

        candidate(1:num_vars) = max(candidate(1:num_vars), lb(1:num_vars));
        candidate(1:num_vars) = min(candidate(1:num_vars), ub(1:num_vars));

    end

    pop(i,:) = candidate;

end

%% Set Up GA
A = [];
B = [];
Aeq = [];
Beq = [];
nlcon = [];

if agriParams.tracking_mode == 1

    num_steps = 23;
    total_constraints = 4 * num_steps * 2;

    A = zeros(total_constraints, num_vars);
    B = ones(total_constraints, 1) * max_slew;

    row = 1;

    for s = 1:4

        offset = 7 + (s-1)*24;

        for h = 1:23

            v1 = offset + h;
            v2 = offset + h + 1;

            % angle(h+1) - angle(h) <= max_slew
            A(row, v1) = -1;
            A(row, v2) =  1;
            row = row + 1;

            % angle(h) - angle(h+1) <= max_slew
            A(row, v1) =  1;
            A(row, v2) = -1;
            row = row + 1;

        end
    end
end
ga_final_pop_set = [];

if mode == 0
    ga_set = ones(1,7);
elseif mode == 1
    ga_set = ones(1,103);
end


%% Run GA (Auto-Detect Single vs. Multi-Objective)
tic;

% 1. Evaluate the initial guess to see how many objectives the wrapper outputs
test_fitness = moga_objective_wrapper(x0, agriParams);
num_objectives = length(test_fitness);

fprintf('Detected %d objective(s). Starting optimization...\n', num_objectives);

if num_objectives == 1
    % --- RUN SINGLE-OBJECTIVE GA (Laser focus) ---
    options_soga = optimoptions('ga', ...
        'PopulationSize', pop_size, ...
        'MaxGenerations', max_gen, ...
        'UseParallel', USE_PARALLEL_PROCESSING, ...
        'Display', 'iter', ...
        'InitialPopulationMatrix', pop); % 'ParetoFraction' and 'PlotFcn' removed!
        
    [ga_solve, fval, exitflag, output, population, scores] = ...
        ga(@(x) moga_objective_wrapper(x, agriParams), ...
        num_vars, A, B, Aeq, Beq, lb, ub, nlcon, options_soga);

else

    % --- RUN MULTI-OBJECTIVE GA (Trade-off floodlight) ---
    options_moga = optimoptions('gamultiobj', ...
        'PopulationSize', pop_size, ...
        'MaxGenerations', max_gen, ...
        'ParetoFraction', 0.35, ...
        'PlotFcn', @gaplotpareto, ...
        'UseParallel', USE_PARALLEL_PROCESSING, ...
        'Display', 'iter', ...
        'InitialPopulationMatrix', pop);
        
    [ga_solve, fval, exitflag, output, population, scores] = ...
        gamultiobj(@(x) moga_objective_wrapper(x, agriParams), ...
        num_vars, A, B, Aeq, Beq, lb, ub, nlcon, options_moga);
end

disp("Time (s) for Optimization: " + toc);
%% Automatically Select the Best Compromise Design (Utopia Point Method)

% 1. Extract the Pareto front scores (Remember, GA minimizes, so these are negative)
% Make them positive for easier mental math
pareto_front = -1 .* fval; 

% 2. Normalize the front so all objectives are on a 0 to 1 scale
% (This prevents Profit from dominating Crop Yield in the distance math)
min_vals = min(pareto_front);
max_vals = max(pareto_front);
range_vals = max_vals - min_vals;

% Prevent divide-by-zero if an objective is entirely flat
range_vals(range_vals == 0) = 1; 

normalized_front = (pareto_front - min_vals) ./ range_vals;

% 3. Define the Utopia Point (The perfect normalized score is 1 for everything)
utopia_point = ones(1, size(normalized_front, 2));

% 4. Calculate the Euclidean distance of every design to the Utopia point
distances = sqrt(sum((normalized_front - utopia_point).^2, 2));

% 5. Find the index of the minimum distance
[~, best_idx] = min(distances);

% 6. Extract your single winning design and its metrics!
winning_design_variables = ga_solve(best_idx, :);
winning_metrics = agrivoltaic_wrapper(winning_design_variables, agriParams);

fprintf('--- BEST COMPROMISE DESIGN SELECTED ---\n');
fprintf('Profit: $%.2f M\n', winning_metrics(2) / 1e6);
fprintf('Crop Yield: %.2f kg\n', winning_metrics(6));

%% Evaluate all Pareto designs with full wrapper metrics

num_pareto = size(ga_solve, 1);
pareto_metrics = zeros(num_pareto, 9);

for k = 1:num_pareto
    pareto_metrics(k, :) = agrivoltaic_wrapper(ga_solve(k, :), agriParams);
end
%% Save Data

save("results/" + moniker + ".mat", ...
    "ga_solve", ...
    "fval", ...
    "exitflag", ...
    "output", ...
    "population", ...
    "scores", ...
    "pareto_metrics");

%% Plot
% 
% front_size = 200;
% val_size = 50;
% utopia_size = 800;
% 
% fig1 = figure;
% 
% title("Pareto Fronts");
% 
% subplot(3,1,1);
% 
% theme(fig1,"light");
% 
% 
% hold on;
% xlabel("Fiscal Profit ($M)");
% ylabel("Annual Raspberry Production (g/m^2)");
% scatter(-scores(:,2), -scores(:,3), val_size, 'black', 'filled');
% scatter(-fval(:,2), -fval(:,3), front_size, 'green', 'filled');
% %scatter(star_x_position_social_profit, star_y_position_berry_production, utopia_size, 'cyan', 'filled', "pentagram");
% legend("Values from GA Population", "Designs on Pareto Front", "Utopia Point", 'Location','southwest');
% %ylim([0 star_y_position_berry_production]);
% %xlim([0 star_x_position_social_profit]);
% hold off;
% 
% subplot(3,1,2);
% 
% hold on;
% xlabel("Fiscal Profit ($M)");
% ylabel("Emission Reduction (t CO2e)");
% scatter(-scores(:,2), -scores(:,1), val_size, 'black', 'filled');
% scatter(-fval(:,2), -fval(:,1), front_size, 'green', 'filled');
% %scatter(star_x_position_social_profit, star_y_position_berry_production, utopia_size, 'cyan', 'filled', "pentagram");
% legend("Values from GA Population", "Designs on Pareto Front", "Utopia Point", 'Location','southwest');
% %ylim([0 star_y_position_berry_production]);
% %xlim([0 star_x_position_social_profit]);
% hold off;
% 
% subplot(3,1,3);
% 
% hold on;
% xlabel("Annual CropProduction (g/m^2)");
% ylabel("Emission Reduction (t CO2e)");
% scatter(-scores(:,3), -scores(:,1), val_size, 'black', 'filled');
% scatter(-fval(:,3), -fval(:,1), front_size, 'green', 'filled');
% %scatter(star_x_position_social_profit, star_y_position_berry_production, utopia_size, 'cyan', 'filled', "pentagram");
% legend("Values from GA Population", "Designs on Pareto Front", "Utopia Point", 'Location','southwest');
% %ylim([0 star_y_position_berry_production]);
% %xlim([0 star_x_position_social_profit]);
% hold off;
% 
% figure_name = "graphs/moga_pareto_" + moniker + ".png";
% saveas(fig1,figure_name);

%% help function 
function angles_out = enforce_slew_curve(angles_in, lb_curve, ub_curve, max_slew)
% Enforces:
%   lb_curve(h) <= angle(h) <= ub_curve(h)
%   abs(angle(h+1) - angle(h)) <= max_slew

    angles_out = angles_in;

    n_iter = 10;

    for iter = 1:n_iter

        % Clamp to bounds
        angles_out = max(angles_out, lb_curve);
        angles_out = min(angles_out, ub_curve);

        % Forward pass
        for h = 2:length(angles_out)

            lower_from_prev = angles_out(h-1) - max_slew;
            upper_from_prev = angles_out(h-1) + max_slew;

            angles_out(h) = min(max(angles_out(h), lower_from_prev), upper_from_prev);

            angles_out(h) = max(angles_out(h), lb_curve(h));
            angles_out(h) = min(angles_out(h), ub_curve(h));

        end

        % Backward pass
        for h = length(angles_out)-1:-1:1

            lower_from_next = angles_out(h+1) - max_slew;
            upper_from_next = angles_out(h+1) + max_slew;

            angles_out(h) = min(max(angles_out(h), lower_from_next), upper_from_next);

            angles_out(h) = max(angles_out(h), lb_curve(h));
            angles_out(h) = min(angles_out(h), ub_curve(h));

        end
    end
end