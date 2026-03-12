function [l]=evalTSP(x);

load cities
% x is the sequence of cities to be visited (row vector)
x=[x x(1)];

for ind=2:length(x)
    l(ind)=sqrt(sum((C(x(ind),:)-C(x(ind-1),:)).^2));
end
l=sum(l);

