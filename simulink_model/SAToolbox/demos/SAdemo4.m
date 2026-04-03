%SAdemo4
% Optimization of Radio Telescope Array (RTA)
clear all
close all
%
global n wavelength xgrid ygrid TRI xygrid
n=27; %VLA number of stations
xygrid=400; %diameter of the array
wavelength=4; % wavelength
uvscale=1;   % radial distribution in the uv plane
% compute nominal distribution
[xgrid,ygrid,TRI] = get_grid(n,xygrid,wavelength,uvscale);


% generate random array stations
phi=2*pi*rand(n,1);
radius=xygrid/2*rand(n,1);
[x,y]=pol2cart(phi,radius);
xy(:,:,1)=x;
xy(:,:,2)=y;

disp('Hit enter to start Array optimization')
keyboard


% compute uv density metric
[uvdensity] = uvdens_circ_eq_area(x,y,n,wavelength,xgrid,ygrid,TRI)

%compute cable length
[cablelength] = minCableplot(x,y)
%[cablelength] = minCable(xy)

% plot xy and uv
xycon=[]; badzone=[]; gridstep=[];
UVplot(x,y,cablelength,uvdensity,n,wavelength,xygrid,xycon,badzone,gridstep)


% use Simulated Annealing to maximize uvdensity alone

xyo=[x' y'];


file_eval='evalRTA';
file_perturb='perturbRTA';
%
%options=[];
    To=2*uvdensity; options(1)=To;
    schedule=2; options(2)=schedule;
    dT=0.8; options(3)=dT;
    neq=20; options(4)=neq;
    nfrozen=3; options(5)=nfrozen;
    diagnostics=0; options(6)=diagnostics;
    plotflag=1; options(7)=plotflag;


[xybest,uvdensitybest,xhist]=SA(xyo,file_eval,file_perturb,options);
xb=xybest(1:n)'; yb=xybest(n+1:2*n)';
% draw best RTA configuration
figure
[cablelengthbest] = minCableplot(xb,yb)
UVplot(xb,yb,cablelengthbest,uvdensitybest,n,wavelength,xygrid,xycon,badzone,gridstep)

