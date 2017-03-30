function res=fracshift3dMB(im, shiftx, shifty, shiftz)
%This function is aimed to increase
%the resolution of x, y and z up to sub-pixel.
%Written by Yongxiang Gao on June 15, 2005
ipx=double(fix(shiftx));
ipy=double(fix(shifty));
ipz=double(fix(shiftz));

fpx=shiftx-ipx;
fpy=shifty-ipy;
fpz=shiftz-ipz;

%to handle negative shifts:
if fpx<0
    fpx=fpx+1;
    ipx=ipx-1;
end
if fpy<0
    fpy=fpy+1;
    ipy=ipy-1;
end

if fpz<0
    fpz=fpz+1;
    ipz=ipz-1;
end

im=double(im);

imagexz  = circshift( im,[ipx ipy+1 ipz ]   );
imageyz  = circshift( im,[ipx+1 ipy ipz]   );
imagez   = circshift( im,[ipx+1 ipy+1 ipz] );
imagexyz = circshift( im,[ipx ipy ipz]    );

imagex  = circshift( im,[ipx ipy+1 ipz+1 ]   );
imagey  = circshift( im,[ipx+1 ipy ipz+1]   );
image   = circshift( im,[ipx+1 ipy+1 ipz+1] );
imagexy = circshift( im,[ipx ipy ipz+1]    );

res =(1-fpz)*((1-fpx)*(1-fpy)*imagexyz+(1-fpx)*fpy*imagexz+fpx*(1-fpy)*imageyz+fpx*fpy*imagez)+fpz*((1-fpx)*(1-fpy)*imagexy+(1-fpx)*fpy*imagex+ fpx*(1-fpy)*imagey+ fpx*fpy*image);
