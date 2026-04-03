% Implement SA for atom configuration sample problem
% dWo - 11/26/2003
% SAdemo0.m

clear all
close all

Ro=[4 5 9 15 ];   % assign initial slots - starting configuration
% evaluate initial guess

[E]=evalAtoms(Ro)
drawnow
disp(['Initial Energy: ' num2str(E)])

disp('Look at initial guess - and type return to start')
keyboard

% prepare inputs for simulated annealing
file_eval='evalAtoms';
file_perturb='perturbAtoms';
%
%options=[];
    To=10; options(1)=To;
    schedule=2; options(2)=schedule;
    dT=0.75; options(3)=dT;
    neq=50; options(4)=neq;
    nfrozen=3; options(5)=nfrozen;
    diagnostics=0; options(6)=diagnostics;
    plotflag=1; options(7)=plotflag; 

tic;
[Rbest,Ebest,Rhist]=SA(Ro,file_eval,file_perturb,options);

CPUtime=toc

if plotflag
figure(1)
[Ebest]=evalAtoms(Rbest(1,:)) 
title('Best Configuration')
end

disp(['There are ' num2str(size(Rbest,1)) ' lowest energy configurations.'])