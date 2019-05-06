function a = AFMCurveRead_Plot(directory)
if nargin ==1
    cd(directory);
end
files = dir('*.txt'); %Check Directory for default filenames

for i = 1:size(files,1)
    a{i} = readtable(files(i).name);
    xshift(i) = max(a{i}.Height_Sensor_nm_Ex);
    yshift = a{i}.Defl_pN_Ex(1) * -1;
    a{i}.Defl_pN_Ex = a{i}.Defl_pN_Ex + yshift; 
    a{i}.Defl_pN_Rt = a{i}.Defl_pN_Rt + yshift; 
    a{i}.Height_Sensor_nm_Ex = a{i}.Height_Sensor_nm_Ex - xshift(i);
    a{i}.Height_Sensor_nm_Rt = a{i}.Height_Sensor_nm_Rt - xshift(i);
end


figure
hold on
for i = 1:size(files,1)  
    plot(a{i}.Height_Sensor_nm_Ex/1000,a{i}.Defl_pN_Ex/1000,'blue')
    plot(a{i}.Height_Sensor_nm_Rt/1000,a{i}.Defl_pN_Rt/1000,'red')
end