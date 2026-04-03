function [dens] = evalRTA(xy)
global  n wavelength xgrid ygrid TRI
nn=length(xy)/2;
x=xy(1:nn)'; y=xy(nn+1:2*nn)';
%uv gridding parameters
uvnum=n*(n-1);  %number of uv points

%calculate uv coordinates and bin them
k=0;
u=zeros(uvnum,1);
v=zeros(uvnum,1);
uvbin=zeros(length(xgrid),1);
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
index=dsearch(xgrid,ygrid,TRI,u,v);
% for i=1:length(index)
%     uvbin(index(i)) = uvbin(index(i))+1;  %actual density
% end
% nombin = ones(length(xgrid),1);
% 
% dens = sum(abs(uvbin(:)-nombin(:)))/length(uvbin);
dens = (uvnum-length(unique(index)))/uvnum;
