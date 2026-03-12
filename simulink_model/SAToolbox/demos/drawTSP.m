function drawTSP(x)

clf
load cities
x=[x x(1)];
    plot(C(:,1),C(:,2),'ko')
    hold on
    text(C(x(1),1),C(x(1),2),'Start')
    for ind=1:length(x)-1
        text(C(ind,1)+0.02, C(ind,2),num2str(ind))
    end
    
    for ind=1:length(x)-1;
        plot(C(x(ind:ind+1),1), C(x(ind:ind+1),2), 'r-')
    end
    l=evalTSP(x);
    title(['Length of TSP Route: ' num2str(l)])
    