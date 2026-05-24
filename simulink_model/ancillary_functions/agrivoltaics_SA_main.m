% Clear
clear;
clc;

% Set Up GA
% Define the objective function and constraints
agrivoltaic_wrapper = @(x) objectiveFunction(x);
A = []; B = []; Aeq = []; Beq = [];
constraint_min = [0 0 0 0 0 0 0 ];% I deleted the last one, so need to check that the constraints match up properly to the variables
constraint_max = [2 2 2 pi./2 pi./2 0 10]; %deleted last one so only 7

% Run GA
[ga_solve,fval,exitflag,output,population,scores] = ga(agrivoltaic_wrapper, ...
    7, A, B, Aeq, Beq, constraint_min, constraint_max);