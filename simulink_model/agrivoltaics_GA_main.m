% Clear and Setup
clear;
clc;

rng default;

addpath(genpath(pwd));

load("agrivoltaics_variable_definition_data.mat");

% Set Up GA
A = []; B = []; Aeq = []; Beq = [];

constraint_min_struct.PV_z_p = 0;
constraint_min_struct.PV_l_p = 0;
constraint_min_struct.PV_w_p = 0;
constraint_min_struct.PV_phi = -pi./2;
constraint_min_struct.PV_sigma = 0;
constraint_min_struct.PV_y_p = 0;
constraint_min_struct.PV_x_p = 0;

constraint_max_struct.PV_z_p = 2;
constraint_max_struct.PV_l_p = 2;
constraint_max_struct.PV_w_p = 2;
constraint_max_struct.PV_phi = pi./2;
constraint_max_struct.PV_sigma = pi./2;
constraint_max_struct.PV_y_p = 10;
constraint_max_struct.PV_x_p = 10;

constraint_min = agriVarStruct2Array(constraint_min_struct);
constraint_max = agriVarStruct2Array(constraint_max_struct);

% Run GA
[ga_solve,fval,exitflag,output,population,scores] = ...
    ga(@agrivoltaic_social_cost_of_carbon_wrapper, ...
    7, A, B, Aeq, Beq, constraint_min, constraint_max);
