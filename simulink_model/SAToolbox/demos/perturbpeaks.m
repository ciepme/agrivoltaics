function [xp]=perturbpeaks(x)
% [xp]=perturbpeaks(x)
%
xp=x;
xlb=-3;
xub=3;
%
        indx=round(rand(1)+1);
        dx=(xub-xlb)*rand(1)+xlb;
        xp(indx)=dx;
% modify on of two coordinates