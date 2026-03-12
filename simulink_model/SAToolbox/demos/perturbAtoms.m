function [rip]=perturbAtoms(ri);
% [E]=perturbAtoms(ri);
% Perturbs the atom configuration - simulated
% annealing sample problem
slots=[1:25];
choices=setdiff(slots,ri);
n=length(ri);
ns=length(choices);
indp=1+round((n-1)*rand);
inds=1+round((ns-1)*rand);
% move atom with index indp to slot with index inds
rip=ri;
rip(indp)=choices(inds);