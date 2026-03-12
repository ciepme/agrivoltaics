function [xyp]=perturbRTA(xy);
% perturb Radio Telescope Array
% one degree-of-freedom - move one station at a time
global n  xygrid
n=size(xy,2)/2;
ind=floor(n*rand(1))+1;

phi=2*pi*rand(1);
radius=xygrid/2*rand(1);
[xm,ym]=pol2cart(phi,radius);

% perturb xy
xyp=xy;
xyp(ind)=xm;
xyp(ind+n)=ym;