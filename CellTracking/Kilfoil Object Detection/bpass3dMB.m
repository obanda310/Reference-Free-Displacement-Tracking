function res=bpass3dMB(image, lnoise, lobject,inputv ,varargin)
%bpass3d is written by Yongxiang Gao and Maria Kilfoil, based on
%the IDL code written by John C. Crocker and David G. Grier.
%the input variable supposed to be contained in varargin is [noclip nopad]
%inputv is to indicate whether there is input for noclip or nopad by
%logical number 1 and 0. 
%Last visit is on June 15, 2005
 length(varargin);
if inputv(1)==1
      noclip=varargin{1,1};
      if inputv(2)==1 ; nopad=varargin{1,2}; end
elseif inputv(2)==1
    nopad=varargin{1,1};
end
nn=length(lnoise);  
no=length(lobject);
if ((nn >1)&(no==1))|((nn==1)&(no>1)) 
   'Both length parameters must be scalars or 3-vectors!'
   return;
end
% do xdirection masks
bb=single(lnoise(1));
w = round( max(lobject(1),2 * bb) );
N = 2*w + 1;
r = single([-w:w])/(2 * bb);
gx = exp( -r.^2 );
gx = gx /sum(gx);
bx = zeros(1,N,'single') - 1./N;
factor = ( sum(gx.^2) - fix(1/N) );
if (nn==1) 
	gy =(gx)';
	gz =(gx)';
	by =(bx)';
	bz =(bx)';
   else
    % do y direction masks
	bb = single(lnoise(2));
	w = round( max(lobject(2),2 * bb) );
	N = 2*w + 1;
	r = single([-w:w])/(2 * bb);	
	gy = exp( -r.^2 );
	gy = gy / sum(gy);
	gy = (gy)';
	by = zeros(1,N,'single') - 1/N;
	by = (by)';                	
	% do z direction masks
	bb = single(lnoise(3)); 
	w = round( max(lobject(3),2 * bb) );
	N = 2*w + 1;
	r = single([-w:w])/(2 * bb);
	gz = exp( -r.^2 );
	gz = gz / sum(gz);
	gz = (gz)';
	bz = zeros(1,N,'single') - 1/N;
	bz = (bz)';             
end
[nx,ny,nf]=size(image);
if (N>=nf)&(inputv(2)==1)	
	% stack is too thin for any data to survive!
	'Warning: data cube thinner than convolution kernel!'
	'Disabling nopad keyword to compensate'
    pad = 1; 
    return;
end

if not (inputv(2)==1)
   if size(lobject)==1 
      padxw = round( max(lobject(1),2 * bb) ); 
      padyw = round( max(lobject(1),2 * bb) ); 
      padzw = round( max(lobject(1),2 * bb) );
   else 
      padxw = round( max(lobject(1),2 * bb) ); 
      padyw = round( max(lobject(2),2 * bb) );
      padzw = round( max(lobject(3),2 * bb) ); 
   end
% pad out the array with average values of corresponding frames
	ave = zeros(1,nf,'single');
	for i = 1:nf
        ave(i) = sum(sum(image(:,:,i)))/(nx*ny); 
	    g = zeros(nx+2*padxw,ny+(2*padyw),nf+(2*padzw),'single');
	    g(:,:,1:padzw) = ave(1);
	    g(:,:,nf-padzw+1:end) = ave(nf);
    end
	for i = 1:nf 
        g(:,:,padzw+i) = ave(i);	
	    g(padxw+1:padxw+nx,padyw+1:padyw+ny,padzw+1:padzw+nf) = single(image);
    end
	nx = nx + (2*padxw);
	ny = ny + (2*padyw);
	nf = nf + (2*padzw);
else
    g=image;
    padxw=0;
    padyw=0;
    padzw=0;
end
    
clear image
b=g;
% do x and y convolutions
for i=padzw+1:nf-padzw
    g(:,:,i)=conv2(gx,gy,g(:,:,i),'same');
    b(:,:,i)=conv2(bx,by,b(:,:,i),'same');
end
% do z convolution
clear temp tep temp2 tep2 d f 
temp=permute(g, [3,1,2]);
for i=padyw+1:ny-padyw
d(:,:,i)= conv2(temp(:,:,i),gz,'same');
end
clear temp
g=permute(d,[2,3,1]);
clear d
temp2=permute(b,[3,1,2]);
clear b
for i=padyw+1:ny-padyw
   f(:,:,i)= conv2(temp2(:,:,i),bz,'same');
end 
clear temp2
b=permute(f,[2,3,1]);
clear d f temp temp2
if not (inputv(2)==1)
    g=g(padxw+1:nx-padxw,padyw+1:ny-padyw,padzw+1:nf-padzw);
    b=b(padxw+1:nx-padxw,padyw+1:ny-padyw,padzw+1:nf-padzw);
end
g=g+b; % This is the final straw on memory: 12X number of original
clear b temp
if inputv(1)==1
    res=g./factor;
else
    res=max(g./factor,0);
end
res;

