function [xgrid,ygrid,TRI] = get_grid(n,xygrid,wavelength,uvscale)

uvnum=n*(n-1);

r_min=1/(1000*xygrid/wavelength);  %corresponds to smallest wavelength

q=uvnum/quadl(@gfun,r_min,xygrid/wavelength,[],[],(uvscale));

rf(1)=xygrid/wavelength;
dummy(1) = 2;  %dummy for while loop termination
k=1;
while dummy(k) > 1 & rf(k) > r_min
    k=k+1;
    x0=rf(k-1);
    if uvscale == -1
        [x,fval] = fzero(@gfun2,x0,[],q,rf(k-1),uvscale,wavelength,xygrid);
    else
        [x,fval] = fzero(@gfun3,x0,[],q,rf(k-1),uvscale,wavelength,xygrid);
    end
    rf(k)=x;
    
    if uvscale == -1
        Npannuli(k-1)=ceil(q*log(rf(k-1)/rf(k)));
    else
        Npannuli(k-1)=ceil(q*(rf(k-1)^(uvscale+1)/(uvscale+1)-rf(k)^(uvscale+1)/(uvscale+1)));
    end
    dummy(k)=Npannuli(k-1);
end
Npannuli=Npannuli(1:end-1);

for i=2:length(rf)
    r(i)=(rf(i)+rf(i-1))/2;
end
r=r(2:end-1);
rf;
r;

rad=0;
theta=0;
for i=1:length(r)-1
    perimeter(i) = 2*pi*r(i);  %calculate perimeters at every rpos
    thetaposn(i) = floor(perimeter(i)/(r(i)-r(i+1)));  %get thetapos at every rpos
    theta_actual(i) = 2*pi/thetaposn(i);  %round for exactness
    h = rand*theta_actual(i);  %add randomness to start of thetapos
    thetapos = h:theta_actual(i):2*pi+h-theta_actual(i);  %actual theta position
    rad = [rad r(i)*ones(length(thetapos),1)'];  %radial vector
    theta = [theta thetapos];  %azimuthal vector
end
[xgrid,ygrid]=pol2cart(theta,rad);
length(xgrid)
TRI=delaunay(xgrid,ygrid);

%%%%%%%%%%%%%%%%%%%
%gfun used in quadl
function [g] = gfun(x,uvscale)
g = x.^uvscale;

%gfun2 used in fzero for uvscale = -1
function [g] = gfun2(r,q,rf,uvscale,wavelength,xygrid)
g = pi*(rf+r)-(rf-r)*q*log(rf/r);

%gfun3 used in fzero for anything but uvscale = -1
function [g] = gfun3(r,q,rf,uvscale,wavelength,xygrid)
g = pi*(rf+r)-(rf-r)*q*(rf^(uvscale+1)/(uvscale+1)-r^(uvscale+1)/(uvscale+1));