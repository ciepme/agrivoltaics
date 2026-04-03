function drawFEM(density)
global Nrows Ncolumns
%
Nrows=8;
Ncolumns=9;
figure
hold on
% 
% draw the density distribution
counter=0;
for indrows=1:Nrows;
    for indcolumns=1:Ncolumns;
    counter=counter+1;
    if density(counter)==1
    patch((indcolumns-1)+[0 1 1 0],(indrows-1)+[0 0 1 1],'k')
end
end
end

axis equal
axis off
title('Structural Configuration')