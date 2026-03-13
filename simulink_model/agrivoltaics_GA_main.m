% Clear
clear;
clc;

rng default;

% Set Up GA
% Define the objective function and constraints
agrivoltaic_wrapper = @(x) objectiveFunction(x);
A = []; B = []; Aeq = []; Beq = [];
constraint_min = [0 0 0 -pi -pi./2 -pi 0 0];
constraint_max = [2 2 2 pi pi./2 pi 10 10];

% Run GA
[ga_solve,fval,exitflag,output,population,scores] = ga(@agrivoltaic_wrapper, ...
    8, A, B, Aeq, Beq, constraint_min, constraint_max);