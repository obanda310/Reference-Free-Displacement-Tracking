objAll = [x_0(:,1),y_0(:,1)];

numObj = length(x_0);
distances = zeros(numObj,1);
for i = 1:numObj
    thisObj = [x_0(i,1),y_0(i,1)];
    distances(i) = sqrt(sum(bsxfun(@minus,objAll,thisObj).^2,2));
end