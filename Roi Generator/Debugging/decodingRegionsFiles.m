clear all
fclose('all');
clear A
clear ver1
clear hor1
clear ver2
clear hor2

filnm = 'Mayerich_1.ovl';
filnm4 = '43ManyPolyline.ovl';

% filnm = 'AS1111.rls';
filnm3 =  '43ManyPolyline.ovl';
filnm2 = '345Spline.rls';

fid=fopen(filnm);
fid2 = fopen(filnm2);
fid3 = fopen(filnm);
fid4 = fopen(filnm2);
fid5 = fopen(filnm3);
fid6 = fopen(filnm);
fid7 = fopen(filnm2);
fid8 = fopen(filnm4);
fid9 = fopen(filnm4);

A = fread(fid3,'uint8');
A2 = fread(fid2,'uint8');
A3 = fread(fid5,'uint1');
A4 = fread(fid8,'uint8');

compA = A;
compA(1:size(A4,1),2) = A4;
compA(1:size(A2,1),3) = A2;
compA(:,4) = compA(:,1)-compA(:,2);
%%
E1 = size(A,1)*4;
E2 = size(A2,1)*4;

dist = 10;

D1= fread(fid6,E1-dist,'uint8');
D2= fread(fid7,E2-dist,'uint8');
%D1= fread(fid6,dist,'uint8');
%D2= fread(fid7,dist,'uint8');
D3= fread(fid6,1,'int32')
D4= fread(fid7,1,'int32')
D5= size(A,1)*4 - D3;
D6= size(A2,1)*4 - D4;

stop4 = 64600/4;
stop2 = size(A,1)-stop4;
stop3 = size(A2,1)-stop4;
stop6 = size(A4,1)-stop4;


B(1:stop2,1) = fread(fid,stop2,'uint32');
B(1:stop3,2) = fread(fid2,stop3,'uint32');
B(1:stop6,7) = fread(fid9,stop6,'uint32');
B(:,6) = B(:,2) - B(:,1);
B(1:stop2,3) = size(A,1)*4 - B(1:stop2,1);
B(1:stop3,4) = size(A2,1)*4 - B(1:stop3,2);
B(1:stop3,8) = size(A4,1)*4 - B(1:stop3,7);
B(1:stop3,5) = 1:stop3;
B(abs(B(1:stop3,4))>64400,4) = 0;
B(abs(B(1:stop3,3))>64400,3) = 0;
B(abs(B(1:stop3,8))>64400,8) = 0;



C(1:stop4,1) = fread(fid,stop4,'uint32');
C(1:stop4,2) = fread(fid2,stop4,'uint32');
C(1:stop4,4) = fread(fid9,stop4,'uint32');
C(1:stop4,3) = C(1:stop4,2) - C(1:stop4,1);
Csum = sum(C(1:stop4,3));


F(:,1) = D1;
F(1:size(D2,1),2) = D2;
%%
stop1 = 372;

B = fread(fid,1000,'uint32');
B(1:stop1,2) = fread(fid2,stop1,'uint32');
B(1:stop1,3) = B(1:stop1,2) - B(1:stop1,1);
B(1:stop1,4) = fread(fid3,stop1,'uint8');
B(1:stop1,5) = fread(fid4,stop1,'uint8');
B(1:stop1,6) = B(1:stop1,5) - B(1:stop1,4);
for i = 1:100
elem1(i,1) = fread(fid3,1,'int32');
end
%%

% 
% 
% %A(pl:pl+35) = fread(fid,36,'uint8=>char');
% A(pl:pl+35,1) = fread(fid,36,'uint8');
% %%
% pl=pl+36;
% filnamnum = fread(fid,1,'uint32');
% A(pl,1)=filnamnum;
% pl=pl+1;
% A(pl:pl+filnamnum,1)= fread(fid,filnamnum+1,'uint8');
% %A(pl:pl+filnamnum)= fread(fid,filnamnum+1,'uint8=>char')
% pl=pl+filnamnum+1;
% numelements = fread(fid,1,'int32');
% A(pl,1)=numelements;
% pl=pl+1;
% totsize = fread(fid,1,'int32');
% A(pl,1)=totsize;
% pl=pl+1;
% A(pl:pl+20,1) = fread(fid,21,'int32');
% pl=pl+21;
% A(pl:pl+31,1) = fread(fid,32,'uint16');
% %A(pl:pl+31) = fread(fid,32,'uint16=>char')
% pl=pl+32;
% A(pl:pl+9,1) = fread(fid,10,'int32');
% pl=pl+10;