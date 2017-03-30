function matx=make_arr(xr,yr,zr)
arr=1:(xr*yr*zr);
matx=reshape(arr, xr,yr,zr);