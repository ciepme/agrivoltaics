function [C]=makeTSP(N)
% [C]=makeTSP(N)
% creates N cities on the interval [-1,1]

C=2*rand(N,2)-1;

if 1
    plot(C(:,1),C(:,2),'ko')
    hold on
    for ind=1:N
        text(C(ind,1)+0.02, C(ind,2),num2str(ind))
    end
end

save cities C