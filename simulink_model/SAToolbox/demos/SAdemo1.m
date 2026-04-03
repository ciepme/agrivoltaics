% Implement SA for peaks function
% dWo - 10/31/2003
% SAdemo1.m

clear all
close all

xo=[-2 -2]';
% evaluate initial guess

plotflag=1;

y(1)= funpeaks(xo);
% plot initial guess
iter = 1;
if plotflag
drawnow
%figure(2)
clf
peaks;
shading interp
hold on
plot3(xo(1,1), xo(2,1), -y(1),'k*');
view([0 0 1])
grid on
drawnow
end

%disp('Look at initial guess - and type return to start')
%keyboard

% prepare inputs for simulated annealing
file_eval='funpeaks';
file_perturb='perturbpeaks';
%
%options=[];
    To=4; options(1)=To;
    schedule=2; options(2)=schedule;
    dT=0.5; options(3)=dT;
    neq=50; options(4)=neq;
    nfrozen=5; options(5)=nfrozen;
    diagnostics=0; options(6)=diagnostics;
    options(7)=plotflag;

Tuser=[];

[xbest,Ebest,xhist]=SA(xo,file_eval,file_perturb,options);

if plotflag
figure(1)
    
% plot best configurations
nhist=max(size(xhist));
xb(1,1)=xhist(1).x(1);
xb(1,2)=xhist(1).x(2);
zb(1)=(-1)*xhist(1).E;
indb=1;
for ind=1:nhist
    if xhist(ind).E<(-1)*zb(indb);
        indb=indb+1;
    xb(indb,1)=xhist(ind).x(1);
    xb(indb,2)=xhist(ind).x(2);
    zb(indb)=(-1)*xhist(ind).E;
    end
end
plot3(xb(:,1),xb(:,2),zb+.5,'mo')
plot3(xb(:,1),xb(:,2),zb+.5,'m-')

end
