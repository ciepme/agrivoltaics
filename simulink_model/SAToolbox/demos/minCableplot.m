function [metric] = minCableplot(x,y)

n=length(x);

hold on
xi(1)=x(1);
yi(1)=y(1);
cable=0;
for k=1:n-1
    for j=1:n
        for i=1:length(xi)
            dist(i,j)=sqrt((xi(i)-x(j))^2+(yi(i)-y(j))^2);
            if dist(i,j)==0
                dist(1:length(xi),j)=inf;
                break
            end
        end
    end
    
    [ii,jj]=find(dist(1:length(xi),1:n)==min(dist(:)));
    k_i=ii(1);
    k_j=jj(1);
    
    xi(k+1)=x(k_j);
    yi(k+1)=y(k_j);
    xii=[xi(k_i);x(k_j)];
    yii=[yi(k_i);y(k_j)];
    plot(xii,yii,'b-')
    
    cable=cable+sqrt((xii(1)-xii(2))^2+(yii(1)-yii(2))^2);
end
axis square

metric=cable;