function [para,perp] = meanTraj(tempPillar,length,pixel)
distX = length/pixel; %microns to pixels
para(1,1) = -1*(mean(tempPillar(:,1))-tempPillar(1,1));
para(1,2) = -1*(mean(tempPillar(:,2))-tempPillar(1,2));
rowScale = sqrt(((distX)^2)/(para(1,1)^2+para(1,2)^2));
para = para*rowScale;
%Create Perpendicular vector for drawing box.
perp(1,1) = -para(1,2);
perp(1,2) = para(1,1);
rowScale = sqrt(((distX)^2)/(perp(1,1)^2+perp(1,2)^2));
perp = perp*rowScale;