function [rowV] = rowVoptimizer(r,pSpaceXY)

Vals2 = [];
safeIdx = find(r.planeGroup==max(r.planeGroup)); %r.NDl & 
r.r = double(r.r);
r.X = double(r.X);
r.Y = double(r.Y);
T = delaunayn(double(r.r(safeIdx,1:2)));
f = @(x)rowVSolve(x,r,safeIdx,T,pSpaceXY);
c = @(x)unitdisk(x,pSpaceXY);
x0 = [0 pSpaceXY];
options1 = optimoptions('simulannealbnd','Display','off');
options2 = optimoptions('fmincon','Display','off');
iterations = 4;

for i = 1:iterations
    if i == 1
        [Vals2(i,1:2)] = simulannealbnd(f,x0,[0 -pSpaceXY],[pSpaceXY pSpaceXY],options1);
    else
        [Vals2(i,1:2)] = simulannealbnd(f,Vals(i-1,1:2),[0 -pSpaceXY],[pSpaceXY pSpaceXY],options1);
    end
    [Vals(i,1:2)] = fmincon(f,Vals2(i,1:2),[],[],[],[],[0 -pSpaceXY],[pSpaceXY pSpaceXY],c,options2);
end

rowV(1,1:2) = Vals(end,1:2);

% cmap = brewermap(iterations,'*spectral');
% for i = 1:1:iterations
%  pts(1,1:3) = [50+Vals2(i,1)*50 50+Vals2(i,2)*50 14];
%  pts(2,1:3) = [50-Vals2(i,1)*50 50-Vals2(i,2)*50 14];
%  hold on
%  plot3(pts(1:2,1),pts(1:2,2),pts(1:2,3),'Color',cmap(i,:))
% end

    function cost = rowVSolve(x,r,safeIdx,T,spacing)
        tXp = double(r.X(safeIdx) + x(1));
        tYp = double(r.Y(safeIdx) + x(2));
        [~,pDist] = dsearchn(r.r(safeIdx,1:2),T,[tXp tYp]);
        pDist(pDist>2,:) = [];
        cost = mean(pDist(:))+abs((spacing)^2-(x(1)^2+x(2)^2));
    end

    function [c,ceq] = unitdisk(x,spacing)
        c = [];
        ceq = sqrt(x(1)^2 + x(2)^2)-spacing;
    end
end


