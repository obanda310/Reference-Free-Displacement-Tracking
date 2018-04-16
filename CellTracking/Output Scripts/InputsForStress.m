function [u,surface,normals,dm] = InputsForStress(directory)
if nargin == 1
    cd(directory);
end
close all
clear all
%% Create X, Y, Z meshes for reference position(surface), displacements, and
% surface normals
load('3Ddata.mat')
%% =====================Displacements=====================================
%% X mesh
clear xq yq zq
[xq,yq] = meshgrid(1:1:size(image.Black,2), 1:1:size(image.Black,1));
xq = xq*raw.dataKey(9,1);
yq = yq*raw.dataKey(9,1);
deformationX = shear.ltLastdX .* shear.coCheck2;
deformationX(isnan(deformationX)) = 0;
u{1}{1}(:,:) = raw.dataKey(9,1)*griddata(shear.rawX1(:),shear.rawY1(:),deformationX(:)/raw.dataKey(9,1),xq,yq);
u{1}{1}(isnan(u{1}{1})) = 0;

%% Y mesh
deformationY = shear.ltLastdY .* shear.coCheck2;
deformationY(isnan(deformationY)) = 0;
u{1}{2}(:,:) = raw.dataKey(9,1)*griddata(shear.rawX1(:),shear.rawY1(:),deformationY(:)/raw.dataKey(9,1),xq,yq);
u{1}{2}(isnan(u{1}{2})) = 0;

%% Z mesh
clear zq
zq = zeros(size(xq));
u{1}{3} = double(vq3(:,:,1));
u{1}{3}(isnan(u{1}{3})) = 0;

%%
figure
imshow(u{1}{1}(:,:),[])
figure
imshow(u{1}{2}(:,:),[])
figure
imshow(u{1}{3}(:,:),[])
%% Make Surface Positions Variable
%            surface{time}{1} = x-coordinates of nodes
%            surface{time}{2} = y-coordinates of nodes
%            surface{time}{3} = z-coordinates of nodes
surface{1}{1} = xq;
surface{1}{2} = yq;
surface{1}{3} = zq;

%% Make Surface Normals Variable
%for now going to treat as a flat plane in Z

normals{1}{1} = zq;
normals{1}{2} = zq;
normals{1}{3} = ones(size(zq));

%% Make dm variable
dm = raw.dataKey(9,1);

%%

%   u:  displacement field vector calculated from FIDVC. Format: cell array,
%      which is a 3D vector (components in x,y,z)  per each time point
%         u{time}{1} = displacement in x-direction at t=time of size MxNxP
%         u{time}{2} = displacement in y-direction at t=time of size MxNxP
%         u{time}{3} = displacement in z-direction at t=time of size MxNxP
%   dm: the spacing of the meshgrid generated in FIDVC
%   surface: cell array containing the surface coordinate grid defined by a
%            rectangular meshgrid (can be non-uniform)
%            surface{time}{1} = x-coordinates of nodes
%            surface{time}{2} = y-coordinates of nodes
%            surface{time}{3} = z-coordinates of nodes
%   normals: cell array containing the normal vectors for each point on the 
%            surface
%            normals{time}{1} = x-component of normal vector
%            normals{time}{2} = y-component of normal vector
%            normals{time}{3} = z-component of normal vector
%   materialModel: 'neohookean', 'linearElastic'; Default = 'linearElastic'
%   materialProps: material Properties of the substrate.
%                  Array: [Young's Modulus, Poisson's ratio]

%[u, dm, surface, normals, materialModel, materialProps, OS] ...
    %= parseInputsTFM(varargin{:});
    model = 'linearElastic';
    properties = [12000,.35];
    save('Inputs.mat','u','dm','surface','normals','model','properties')
[Fij, Sij, Eij, Uhat, ti, tiPN] = fun3DTFM(u,dm,surface,normals,model,properties);

