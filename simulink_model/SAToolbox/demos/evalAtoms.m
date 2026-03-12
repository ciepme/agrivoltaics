function [E]=evalAtoms(ri);
% [E]=evalAtoms(ri);
% Computes the total energy for a simple system
% (configuration) of N atoms in xy space 
% Useful as a sample problem for explaining SA.
% dWo, November 2003

slots=[ 1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 4 4 4 4 4 5 5 5 5 5;
        1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5]';

plotflag=1;
    
g=10;  % strong gravity
% g=0.1; % weak gravity
m=1;   % atomic mass

% compute system energy

for ind=1:length(ri)
    E(ind)=m*g*slots(ri(ind),2)+ sum(sqrt(sum([(repmat(slots(ri(ind),:),length(ri),1)-slots(ri,:)).^2]')));
end

 %  keyboard
E=sum(E);  %total energy

if plotflag
    [xa,ya]=pol2cart([0:pi/8:2*pi]', 0.5*ones(17,1)); 
    figure(1), clf
    plot([0 6 6 0 0]',[0 0 6 6 0]','r-','linewidth',3), hold on
    for ind=1:length(ri)
        patch(slots(ri(ind),1)+xa, slots(ri(ind),2)+ya,'k')
        text(slots(ri(ind),1)-.2,slots(ri(ind),2), num2str(ind),'color','w','fontsize',20)
    end
    grid on, axis equal
      axis([-1 7 -1 7])
    title('Atom Configuration - Sample Problem')
    drawnow
end

     