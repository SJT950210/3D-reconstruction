function Imgmake()

p=load('angel.mat');
f0=p.Grid;
[m,n,h]=size(f0);
figure; show3D(f0,0.5);

% compute the distance of all the point cloud
dx=double(bwdist(f0));

figure; imagesc(fliplr(rot90(dx(:,:,38)))); colormap(jet); colorbar;

% compute gradient magnitude of binary image u
f=double(dx<3);     %theshold distance function to compute f
% f=1./(dx.^2+1);         %use edge detector idea to compute f
figure; imagesc(fliplr(rot90(f(:,:,38)))); colormap(gray); colorbar;

% initialize the boundaries of binary image u
u = zeros(m,n,h)+10000;
u(:,1,:) = zeros(m,1,h);
u(:,n,:) = zeros(m,1,h);
u(1,:,:) = zeros(1,n,h);
u(m,:,:) = zeros(1,n,h);
u(:,:,1) = zeros(m,n,1);
u(:,:,h) = zeros(m,n,1);

% use fast sweep to compute the volume image
% u = fastSweep3DCaseC(u,f);
u = fast8Sweeps(f,ones(m,n,h),u);
figure; show3D(u,0.8);
% figure,show3D(double(u),0.8);
save('intangel','u');

function show3D(phi,thr)
p2=patch(isosurface(phi,thr));
set(p2,'FaceColor','b','EdgeColor','none');
isonormals(phi,p2);
% view(-90,60);
view(90,0);
axis image;
axis off;
material shiny;
camlight infinite
lighting flat;