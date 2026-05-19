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
file_suffix = "_berry_pop100";

%GA hyperparameter settings
pop_size = 20;
max_gen = 20;
stall = 1;

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
        'Display', 'iter');
elseif GA_SELECTOR == 2
    rng(2);
end

%% Set Up GA
A = []; B = []; Aeq = []; Beq = [];
nlcon = [];

ga_final_pop_set = [];
ga_set = ones(1,7);

% Run GA

[ga_solve,fval,exitflag,output,population,scores] = ...
    gamultiobj(@(x) moga_objective_wrapper(x, agriParams), ...
    num_vars, A, B, Aeq, Beq, lb, ub, ...
    nlcon, options);