function r=llmx3dMB(a, sep, pad)
% a 3d-happy version of DGGs 'local max', which does NOT use
% the dilation algorithm.  NB: the data MUST be padded, or it
% can crash!! Best of all, 'sep' is a float 3-vector. 'sep' is
% the actual minimum distance between two maxima (i.e. the bead
% diameter or so.  If the image is 'padded', tell the code the
% padding size to save it a little work.
% Translated by Yongxiang Gao in 2005 from IDL code developped by John
% Crocker and David Grier

allmin = min(a(:));
allmax = max(a(:));
a(1,1,1) = allmin;        
a=fix(((255.0+1)*single(a-allmin)-1)/single(allmax-allmin));
allmin = a(1,1);			
[nx, ny, nz]= size(a);
bignum = nx*ny*nz;
% the diameter of the local max algorithm is 2*sep.
% extent is the next biggest odd integer.
extent = fix(sep*2) + 1	; % i.e. mask is diameter 2*sep and extent is a integer
extent = extent + mod(extent+1, 2);
rsq = lrsqd3dMB(extent, [1,1], sep(1,2)/sep(1,1), sep(1,3)/sep(1,1));
mask = rsq < (sep(1))^2;
% cast the mask into a one dimensional form-- imask!
bmask = zeros(nx,ny,extent(3));
bmask(1:extent(1),1:extent(2),:) = mask;
imask = find(bmask > 0) + bignum -(nx*ny*(fix(extent(3)/2)))-(nx*fix((extent(2)/2))) -fix(extent(1)/2);
% let's try Eric's hash table concept.
% set percentile to 0. if you want ever voxel to be a potential maximum
% set it to 0.8 or so if you have lots of tiny spikes, to run faster
percentile = 0.7;
hash = ones(nx,ny,nz,'uint8');
ww = find(a(:,:,pad(3)+1:nz-pad(3)) > allmin) + (nx*ny*pad(3)); %seems no problem before this compared to IDL
nww=length(ww);   
rd=rand(nww,1);
[junk,ss] = sort(a(ww)); 
s = ww(ss(fix(percentile*nww)+1:end)); 
ww = 0;
s = [fliplr(s'),1]; % so it knows how to stop!, since a(1) contains the minimum intensity
idx = 1; 
rr = s(idx); 
m = a(rr);
r = -1;
i = -1;
erwidx=length(s);  

while 1
% get the actual local max in a small mask
indx= mod(rr+imask-2,bignum)+1;  
actmax = max(a(indx(:)));
% if our friend is a local max, then nuke out the big mask, update r
if m >= actmax
  r = [r,rr];
hash(indx) = 0;
else
w = find(a(indx) < m);
        nw=length(w);
		if nw > 0 
            indx2= mod(rr+imask(w)-2,bignum)+1;
            hash(indx2) = 0; 
        end
end

% get the next non-nuked id
while 1  
idx = idx+1;
if hash(s(idx)) == 1 || idx >= erwidx
     break, end
end
 	if (idx < erwidx) 
		rr = s(idx);
		m = a(s(idx));
    else
		m = allmin;
    end
  if (m <= allmin), break, end
end

if numel(r) > 1 
    r = r(2:end);
else
    r=-1;
end
x =fix( mod(mod( (r-1), (nx*ny)) ,nx )); 
y =fix( mod( (r-1) , (nx*ny))/nx);
z =fix( (r-1)  / (nx*ny))  ;
w=find(x>pad(1) & x<nx-pad(1)-1 & y>pad(2) & y<ny-pad(2)-1 & z>pad(3) & z<nz-pad(3)-1);
clear r;
nw=length(w);
if nw > 0 
    r(:,1:3)=[x(w)' y(w)' z(w)'];
else
    r=-1; 
end