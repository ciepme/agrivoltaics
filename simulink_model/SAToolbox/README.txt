Simulated Annealing Matlab Function
===================================

This directory contains a MATLAB function (SA.m)
to conduct single objective optimization. The function
call is as follows:

[xbest,Ebest,xhist]=SA(xo,file_eval,file_perturb,options);

 Single Objective Simulated Annealing (SA) Algorithm

This function is a generic implementation of the single objective
Simulated Annealing (SA) algorithm first proposed by Kirkpatrick,
Gelatt and Vecchi. The algorithm tries to improve upon an initial
configuration, xo, by evaluating perturbed configurations. When the
system reaches the "frozen" state, the algorithm stops and the best
configuration and search history are returned. The user can choose 
from one of two cooling schedules: linear or exponential.


The inputs and outputs of the function are contained in
the help section of SA.m

>> help SA

The best way to learn the use of this program is to run and
understand any of the demonstration files contained in the "demos"
directory.

SAdemo0 - four atom placement problem
SAdemo1 - demo of SA on MATLAB peaks function
SAdemo2 - demo of SA for Travelling Salesman Problem (TSP)
SAdemo3 - demo of SA for structural topology optimization
SAdemo4 - demo of SA for telescope array placement problem

Each of these demos can be run by setting the "demos" directory
as the current working directory in MATLAB and typing (e.g. for the first demo):

>> SAdemo1   

An article (SA.pdf) accompanies this toolbox with explanations
about the theoretical background and an in-depth explanation of
the algorithm.

Please direct any comments or questions (and bug reports) to:
O. de Weck - deweck@mit.edu

(c) Massachusetts Institute of Technology, 2004