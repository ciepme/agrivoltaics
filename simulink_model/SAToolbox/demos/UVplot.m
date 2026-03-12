function [] = UVplot(x,y,cable,metric,n,wavelength,xygrid,xycon,badzone,gridstep)

% %normal mode
figure
subplot(1,2,1)

for i=1:length(badzone)
    hold on
    xbad=[xycon(badzone(i),1)-gridstep/2,xycon(badzone(i),1)+gridstep/2];
    ybad=[xycon(badzone(i),2)-gridstep/2,xycon(badzone(i),2)+gridstep/2];
    xxbad=[xbad(1),xbad(1),xbad(2),xbad(2)];
    yybad=[ybad(1),ybad(2),ybad(2),ybad(1)];
    fill(xxbad,yybad,'y')
end

hold on
plot(x,y,'ro')
mincableplot(x,y);

title(['Cable Length = ',num2str(cable)])
axis([-xygrid/2 xygrid/2 -xygrid/2 xygrid/2])
axis square

k=0;
u=zeros(n*(n-1),1);
v=zeros(n*(n-1),1);
for i=1:n
    for j=1:n
        if i==j
            continue
        end
        k=k+1;
        u(k)=(x(i)-x(j))/wavelength;
        v(k)=(y(i)-y(j))/wavelength;
    end
end

subplot(1,2,2)
plot(u,v,'.','Markersize',1)
title(['UV Density = ',num2str(metric)])
axis([-xygrid/wavelength xygrid/wavelength -xygrid/wavelength xygrid/wavelength])
axis square
