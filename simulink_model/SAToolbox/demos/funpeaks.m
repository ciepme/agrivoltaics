function z=funpeaks(X0)

x = X0(1);
y = X0(2);
z =  -(3*(1-x).^2.*exp(-(x.^2) - (y+1).^2) ... 
   - 10*(x/5 - x.^3 - y.^5).*exp(-x.^2-y.^2) ... 
   - 1/3*exp(-(x+1).^2 - y.^2));

%z=z*(-1);
if 0
% plot new guess
figure(1)
drawnow
hold on
zm = (-1)*z;
plot3(x, y, zm,'k*');
end



    
