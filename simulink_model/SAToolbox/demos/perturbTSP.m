function [xp]=perturbTSP(x);

% perturb TSP route by flipping the position of two cities
N=length(x);
ind1=round((N-1)*rand)+1;
ind2=ind1;
while ind2==ind1
    ind2=round((N-1)*rand)+1;
end

xtmp=x;
xtmp(ind1)=x(ind2);
xtmp(ind2)=x(ind1);
xp=xtmp;


