function res=lrsqd3dMB(extent,inputv,varargin)
if inputv(1)==0 yratio=1; else yratio=varargin{1,1}; end
if inputv(2)==0 zratio=1; else zratio=varargin{1,2}; end 
if length(extent)==1
    ext = zeros(1,3)+extent;
else
    ext=extent;
end
    x=ext(1);
    y=ext(2);
    z=ext(3);
   
    r2 = zeros(x, y, z); 
    xc = double(x-1) / 2;
    yc = double(y-1) / 2;
    zc = double(z-1) / 2;
 
    yi = zeros(1,x) +1;
    xi = zeros(1,y) +1;
   
    xa = [-xc:x-xc-1];
    xa = xa.^2;
    ya = ([-yc:y-yc-1])/yratio;
    ya = ya.^2;
    za = ([-zc:z-zc-1])/zratio;
    za = za.^2 ;
    size(xi);
    size(xa);
for k=1:z
    r2(:,:,k) = (xi'*xa) + (ya'*yi) + za(k);
end
res=r2;

