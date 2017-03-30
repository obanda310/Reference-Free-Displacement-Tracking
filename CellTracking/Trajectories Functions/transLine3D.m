function [xyzFinal,rowP] = transLine3D(rowV,rowP,ub)
options = optimset('Display', 'off') ;
rowV = double(rowV)';
rowP = double(rowP);
x0 = rowV(1,1);
y0 = rowV(2,1);
z0 = rowV(3,1);
iPx = mean(rowP(:,1));
iPy = mean(rowP(:,2));
iPz = mean(rowP(:,3));
lb(1,1) = 0;
lb(2,1) = 0;
lb(3,1) = 0;

[xyzFinal] = fmincon(@Line3DCost,[iPx;iPy;iPz],[],[],[],[],lb,ub,[],options);
        for j=1:size(rowP,1)
            rowP(j,4) = norm(cross(rowV,rowP(j,1:3)'-xyzFinal))/norm(rowV);
        end
    function [totalCost] = Line3DCost(x)
        xn = x(1);
        yn = x(2);
        zn = x(3);
        xf = xn + x0;
        yf = yn + y0;
        zf = zn + z0;
        for i=1:size(rowP,1)
            distance(i,1) = norm(cross(rowV,rowP(i,1:3)'-[xf;yf;zf]))/norm(rowV);
        end
        totalCost = sum(distance(:,1));
    end
end

