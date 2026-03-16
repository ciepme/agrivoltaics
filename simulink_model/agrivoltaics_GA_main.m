% Clear
clear;
clc;

rng default;

load("agrivoltaics_variable_definition_data.mat");

% Set Up GA
A = []; B = []; Aeq = []; Beq = [];
constraint_min = [0 0 0 -pi -pi./2 -pi 0 0];
constraint_max = [2 2 2 pi pi./2 pi 10 10];

% Run GA
[ga_solve,fval,exitflag,output,population,scores] = ...
    ga(@agrivoltaic_social_cost_of_carbon_wrapper, ...
    8, A, B, Aeq, Beq, constraint_min, constraint_max);