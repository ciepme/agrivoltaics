%% Clear and Setup

clc; clear; close all; % Clear command window, workspace, and close figures

%% Helper function
function fitness = moga_objective_wrapper(x, params)
    % Run the standard wrapper to get all outputs
    % raw = [Emissions(1), Profit(2), SocialCost(3), PVRev(4), CropRev(5), Biomass(6), Panels(7), Energy(8)]
    raw = agriObjArray2Struct(agrivoltaic_wrapper(x, params));

    negative_emission_reduction = -1.*raw.emission_reduction;
    negative_fiscal_profit = -1.*raw.fiscal_profit;
    negative_yearly_biomass = -1.*raw.yearly_biomass;
    
    fitness = [negative_emission_reduction negative_fiscal_profit negative_yearly_biomass];
end

%% Set Up GA

addpath(genpath(pwd));
agrivoltaics_variable_definition;

%User define statements
GA_SELECTOR = 1;
USE_PARALLEL_PROCESSING = true;

%GA hyperparameter settings
pop_size = 300;
max_gen = 20;
stall = 1;
mode = agriParams.tracking_mode;

moniker = "pop" + pop_size + "_max_gen" + max_gen + "_mode" + mode;

% Tell GA exactly how many variables we are optimizing (7 or 103, dependent
%create a smart guess to seed GA population for if tracking mode
x0_base = agriVarStruct2Array(agriVar, agriParams);
num_vars = length(lb);

%1 for basic GA
if GA_SELECTOR == 1
    rng(1);
    options = optimoptions('gamultiobj', ...
        'PopulationSize', pop_size, ...
        'MaxGenerations', max_gen, ...
        'ParetoFraction', 0.35, ...
        'PlotFcn', @gaplotpareto, ...
        'UseParallel', USE_PARALLEL_PROCESSING, ...
        'Display', 'iter');
elseif GA_SELECTOR == 2
    rng(2);
end

%% Set Up GA
A = []; B = []; Aeq = []; Beq = [];
nlcon = [];

ga_final_pop_set = [];

if mode == 0
    ga_set = ones(1,7);
elseif mode == 1
    ga_set = ones(1,103);
end

%% Run GA

tic;

[ga_solve,fval,exitflag,output,population,scores] = ...
    gamultiobj(@(x) moga_objective_wrapper(x, agriParams), ...
    num_vars, A, B, Aeq, Beq, lb, ub, ...
    nlcon, options);

disp("Time (s) for MOGA: " + toc);

%% Save Data

save("results/MOGA_data_" + moniker + ".mat", "ga_solve", "fval", "exitflag", "output", "population", "scores");

%% Plot

front_size = 200;
val_size = 50;
utopia_size = 800;

fig1 = figure;

title("Pareto Fronts");

subplot(3,1,1);

theme(fig1,"light");


hold on;
xlabel("Fiscal Profit ($M)");
ylabel("Annual Raspberry Production (g/m^2)");
scatter(-scores(:,2), -scores(:,3), val_size, 'black', 'filled');
scatter(-fval(:,2), -fval(:,3), front_size, 'green', 'filled');
%scatter(star_x_position_social_profit, star_y_position_berry_production, utopia_size, 'cyan', 'filled', "pentagram");
legend("Values from GA Population", "Designs on Pareto Front", "Utopia Point", 'Location','southwest');
%ylim([0 star_y_position_berry_production]);
%xlim([0 star_x_position_social_profit]);
hold off;

subplot(3,1,2);

hold on;
xlabel("Fiscal Profit ($M)");
ylabel("Emission Reduction (t CO2e)");
scatter(-scores(:,2), -scores(:,1), val_size, 'black', 'filled');
scatter(-fval(:,2), -fval(:,1), front_size, 'green', 'filled');
%scatter(star_x_position_social_profit, star_y_position_berry_production, utopia_size, 'cyan', 'filled', "pentagram");
legend("Values from GA Population", "Designs on Pareto Front", "Utopia Point", 'Location','southwest');
%ylim([0 star_y_position_berry_production]);
%xlim([0 star_x_position_social_profit]);
hold off;

subplot(3,1,3);

hold on;
xlabel("Annual Raspberry Production (g/m^2)");
ylabel("Emission Reduction (t CO2e)");
scatter(-scores(:,3), -scores(:,1), val_size, 'black', 'filled');
scatter(-fval(:,3), -fval(:,1), front_size, 'green', 'filled');
%scatter(star_x_position_social_profit, star_y_position_berry_production, utopia_size, 'cyan', 'filled', "pentagram");
legend("Values from GA Population", "Designs on Pareto Front", "Utopia Point", 'Location','southwest');
%ylim([0 star_y_position_berry_production]);
%xlim([0 star_x_position_social_profit]);
hold off;

figure_name = "graphs/moga_pareto_" + moniker + ".png";
saveas(fig1,figure_name);