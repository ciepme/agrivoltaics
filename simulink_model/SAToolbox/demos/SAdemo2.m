% SAdemo2
% Implement demo for Travelling Salesman Problem (TSP)
clear all
close all
N=15;
xo=[1:N];

[C]=makeTSP(N);

% draw initial route
drawTSP(xo);
disp('Hit enter to solve TSP');
keyboard

%% prepare inputs for simulated annealing
file_eval='evalTSP';
file_perturb='perturbTSP';
%
%options=[];
    To=N; options(1)=To;
    schedule=2; options(2)=schedule;
    dT=0.8; options(3)=dT;
    neq=40; options(4)=neq;
    nfrozen=2; options(5)=nfrozen;
    diagnostics=0; options(6)=diagnostics;
    plotflag=1; options(7)=plotflag;


[xbest,Ebest,xhist]=SA(xo,file_eval,file_perturb,options);

% draw best route
figure
drawTSP(xbest);
xlabel('Shortest Route')
