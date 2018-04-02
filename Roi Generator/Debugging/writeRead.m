filename = '34ManyPolyline.ovl';
test = fopen('TestBinary.test','w');
for i = 1:size(filename,2)
    fwrite(test,uint8(filename(i)),'uint8')
end
fclose(test);
test2 = fopen('TestBinary.test');
testbin = fread(test2,'uint8');
fclose('all')
