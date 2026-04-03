% SAdemo3
% Implement Structural Topology Optimization with
% Simulated Annealing
%
clear all
close all
% domain is 8x15 cells
global Nrows Ncolumns
Nrows=8; Ncolumns=9;
density=zeros(1,8*9);
density(19:54)=ones(1,36);
% mass constratint is 50%
xo=density;
compliance_xo = evalFEM(xo);
% initial guess is a simple cantilever
drawFEM(density)
%
disp(['Intial Compliance x_o is: ' num2str(compliance_xo)])
%
% Improve Mass Dsitribution with Simulated Annealing
%
disp('Hit enter to optimize Structure');
keyboard

%% prepare inputs for simulated annealing
file_eval='evalFEM';
file_perturb='perturbFEM';
%
%options=[];
    To=compliance_xo/10; options(1)=To;
    schedule=2; options(2)=schedule;
    dT=0.8; options(3)=dT;
    neq=20; options(4)=neq;
    nfrozen=4; options(5)=nfrozen;
    diagnostics=0; options(6)=diagnostics;
    plotflag=1; options(7)=plotflag;


[xbest,Ebest,xhist]=SA(xo,file_eval,file_perturb,options);

% draw best structural configuration
figure
drawFEM(xbest);
