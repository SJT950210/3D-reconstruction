function Eulorbunny()
close all

% p=load('head.mat');
p=load('bunny0.mat');
f0=p.Grid;
[m,n,h]=size(f0);

% [x, y, z] = ind2sub(size(f0),find(f0));
% A = [x, y, z];
% fid = fopen('Testobj.obj','w');
% write_vertices(fid, A);
% fclose(fid);


% function write_vertices(fid,V)
% switch size(V,2)
%     case 1
%         for i=1:size(V,1)
%             fprintf(fid,'v %5.5f\n', V(i,1));
%         end
%     case 2
%         for i=1:size(V,1)
%             fprintf(fid,'v %5.5f %5.5f\n', V(i,1), V(i,2));
%         end
%     case 3
%         for i=1:size(V,1)
%             fprintf(fid,'v %5.5f %5.5f %5.5f\n', V(i,1), V(i,2), V(i,3));
%         end
%     otherwise
% end
% fprintf(fid,'# %d vertices \n', size(V,1));

[x,y,z] = ind2sub(size(f0),find(f0));
figure; scatter3(y,x,z,'black+')%figure1
axis equal 
grid off
% view(-103,4);
view(-90,90);
% distance on all pixels
dx=double(bwdist(f0));
% Initalize level set function as distance function
% p=load('headImg.mat');
% p=load('initialHand4.mat');
p=load('initialBunny0.mat');
u=p.u;

figure; imagesc(rot90(u(:,:,38))); colormap(jet); colorbar%axis off; axis equal;%figure2
phi=double(u>0.5);
figure; imagesc(rot90(phi(:,:,38))); colormap(jet); colorbar%axis off; axis equal;%figure3
% phi=double(bwdist(phi)-bwdist(1-phi)+phi-.5);
phi=double(bwdist(phi)-bwdist(1-phi));
figure; imagesc(rot90(phi(:,:,38))); colormap(jet); colorbar%figure4
figure
show3D(phi)

% linearly stretch the binary image u into [0 255]
u=255*(u-min(u(:)))/(max(u(:))-min(u(:)));

% Initalize level set function as distance function
% phi=zeros(m,n,h);
% phi(10:m-10,10:n-10,10:h-10)=1;
% phi=double(bwdist(phi)-bwdist(1-phi)+phi-.5);

% set parameters1
% BalloonForceCoeff=0; % range from 1e-1 to 5e-1 lagre solves local minmizer
% SDFCoeff=1e-3; % range from 1e-4 to 1e-2  the larger the closer to SDF
% lamda=1e-4;    % range from 1e-4 to 1e-2  the smaller the smoother
% timestep=10;   % the larger the faster
gamma=0.5;
BalloonForceCoeff=1; % range from 1e-1 to 5e-1 lagre solves local minmizer
SDFCoeff=1e-2; % range from 1e-4 to 1e-2  the larger the closer to SDF
lamda=1e-4;    % range from 1e-4 to 1e-2  the smaller the smoother
timestep=15; % the larger the faster
t=1;
a=0.01;
b=1;
miu1=2.5;
miu2=0.1;

% compute the distance function at half point

p1=zeros(m,n,h);
p2=zeros(m,n,h);
p3=zeros(m,n,h);

Hphi=heaviside(phi,3);
insideValue=u.*Hphi;
outsideValue=u.*(1-Hphi);
c1=sum(insideValue(:))./(sum(Hphi(:))+eps);
c2=sum(outsideValue(:))./(sum(1-Hphi(:))+eps);
dataTerm=(u-c1).^2-(u-c2).^2;
    
tic
for ii=1:50
    ii
    
    
    
    % update level set function by explicit iteration
    div_p=Div(p1,p2,p3);
    GACTerm=gamma.*dx.*div_p+BalloonForceCoeff*dx;%
    lapPhi=EXLaplace(phi);
    curvaturePhi=EXhalfCurvature(phi,1,1,1,1,1,1);
    R_term=lapPhi-curvaturePhi;
    diracPhi=dirac(phi,3);
    phi=phi+timestep*((GACTerm-lamda*dataTerm).*diracPhi+SDFCoeff*R_term);
    
    % update p 
    [ux,uy,uz]=grad(u);
%     diracPhi=dirac(phi,3);
    abs_u=sqrt(ux.^2+uy.^2+uz.^2);
    phi1=ux./(abs_u+eps);
    phi2=uy./(abs_u+eps);
    phi3=uz./(abs_u+eps);
    div_u=Div(phi1,phi2,phi3);
%     div_u=Div(phi1,phi2);
    g=a+b.*(div_u).^2;
    Hphi=heaviside(phi,2);
    [px,py,pz]=grad(dx.*Hphi);
%     diracPhi=dirac(phi,3);
%     for i=1:2
        p1=p1+t.*(gamma.*px);
        p2=p2+t.*(gamma.*py);
        p3=p3+t.*(gamma.*pz);
%     end
    abs_p=sqrt(p1.^2+p2.^2+p3.^2);
    p1=p1.*g./max(abs_p,g);
    p2=p2.*g./max(abs_p,g);
    p3=p3.*g./max(abs_p,g);
    
    

    if mod(ii,50)==0
        pause(.1);
        figure; scatter3(y,x,z,'black+')
        axis image
        grid off
%         view(-103,4);
        view(-90,90);
        axis off
        hold on
        show3D(phi);
    end
    
    Hphi2=heaviside(-phi,3);
%     diracPhi=dirac(phi,3);
%     [px,py,pz]=grad(phi);
%     abs_p=sqrt(px.^2+py.^2+pz.^2);
    
    abs_p=sqrt(p1.^2+p2.^2+p3.^2);
    Hphi=heaviside(phi,2);
    insideValue=u.*Hphi;
    outsideValue=u.*(1-Hphi);
    c1=sum(insideValue(:))./(sum(Hphi(:))+eps);
    c2=sum(outsideValue(:))./(sum(1-Hphi(:))+eps);
    dataTerm=miu1.*(u-c1).^2-miu2.*(u-c2).^2;
%     if rem(ii,10)==0
%        [~,v]=isosurface(phi,0);
%        save('100hp','v');
%     end
    Q12=miu1.*(u-c1).^2.*Hphi+miu2.*(u-c2).^2.*(1-Hphi);
    E1=sum(sum(sum(dx.*abs_p.*Hphi)));
    E2=sum(sum(sum(Q12.*Hphi)));
    E3=sum(sum(sum(dx.*Hphi2)));
    E(ii)=E1+lamda*E2+BalloonForceCoeff*E3;
%     if (ii>2)
%         if abs((E(ii)-E(ii-1))/E(ii))<0.003
%             break;
%         end
%     end
end
toc
figure;scatter3(y,x,z,'black+');showpoint(phi);
[~,v]=isosurface(phi,0);
figure;plot(E);xlabel('Iterations');ylabel('Energy');legend('Energy/Iterations');
save('Ebunny','E');
% save('50b','v');
% figure;scatter3(y,x,z,'black+');axis off;show3D(phi);
% figure;plot(E);xlabel('Iterations');ylabel('Energy');legend('Energy/Iterations');
% save('E2hand4','E');
% showpoint(phi)
% function show3D(phi)
% p2=patch(isosurface(phi,0));
% set(p2,'FaceColor','b','EdgeColor','none');
% isonormals(phi,p2);
% % view(-103,4);
% view(-103,4);
% axis equal
% material shiny;
% camlight right
% lighting phong;
function showpoint(phi)
[~,v]=isosurface(phi,0);
[m2,~]=size(v);
for i=1:m2
    x1(i,1)=v(i,1);
    y1(i,1)=v(i,2);
    z1(i,1)=v(i,3);
end
figure; scatter3(x1,y1,z1,'black.')
axis image
axis off
grid off
% view(-90,90);
view(-90,90);

function show3D(phi)
p2=patch(isosurface(phi,0));
set(p2,'FaceColor','b','EdgeColor','none');
isonormals(phi,p2);
axis image
% view(270,90);
view(-90,90);
material shiny;
camlight right
lighting phong;

function edgeTerm=EXhalfCurvature(phi,g1,g2,g3,g4,g5,g6)
% [m,n,k]=size(phi);
A1 = circshift(phi,[0 -1 0]); %A1(:,n,:) = A1(:,n-1,:); %phi(i,j+1,k)
A2 = circshift(phi,[0 1 0]);  %A2(:,1,:) = A2(:,2,:);   %phi(i,j-1,k)
A3 = circshift(phi,[-1 0 0]); %A3(m,:,:) = A3(m-1,:,:); %phi(i+1,j,k)
A4 = circshift(phi,[1 0 0]);  %A4(1,:,:) = A4(2,:,:);   %phi(i-1,j,k)
A5 = circshift(phi,[0 0 -1]); %A5(:,:,k) = A5(:,:,k-1); %phi(i,j,k+1)
A6 = circshift(phi,[0 0 1]);  %A6(:,:,1) = A6(:,:,2);   %phi(i,j,k-1)

A7  = circshift(phi,[-1 -1 0]); %A7(:,n,:)  = A7(:,n-1,:);  A7(m,:,:)  = A7(m-1,:,:);   %phi(i+1,j+1,k)
A8  = circshift(phi,[1 -1 0]);  %A8(:,n,:)  = A8(:,n-1,:);  A8(1,:,:)  = A8(2,:,:);     %phi(i-1,j+1,k)
A9  = circshift(phi,[0 -1 -1]); %A9(:,n,:)  = A9(:,n-1,:);  A9(:,:,k)  = A9(:,:,k-1);   %phi(i,j+1,k+1)
A10 = circshift(phi,[0 -1 1]);  %A10(:,n,:) = A10(:,n-1,:); A10(:,:,1) = A10(:,:,2);    %phi(i,j+1,k-1)
A11 = circshift(phi,[-1 1 0]);  %A11(:,1,:) = A11(:,2,:);   A11(m,:,:) = A11(m-1,:,:);  %phi(i+1,j-1,k)
A12 = circshift(phi,[1 1 0]);   %A12(:,1,:) = A10(:,2,:);   A12(1,:,:) = A8(2,:,:);     %phi(i-1,j-1,k)
A13 = circshift(phi,[0 1 -1]);  %A13(:,1,:) = A13(:,2,:);   A13(:,:,k) = A13(:,:,k-1);  %phi(i,j-1,k+1)
A14 = circshift(phi,[0 1 1]);   %A14(:,1,:) = A14(:,2,:);   A14(:,:,1) = A11(:,:,2);    %phi(i,j-1,k-1)
A15 = circshift(phi,[-1 0 -1]); %A15(m,:,:) = A15(m-1,:,:); A15(:,:,k) = A15(:,:,k-1);  %phi(i+1,j,k+1)
A16 = circshift(phi,[-1 0 1]);  %A16(m,:,:) = A16(m-1,:,:); A16(:,:,1) = A15(:,:,2);    %phi(i+1,j,k-1)
A17 = circshift(phi,[1 0 -1]);  %A17(:,:,k) = A17(:,:,k-1); A17(1,:,:) = A17(2,:,:);    %phi(i-1,j,k+1)
A18 = circshift(phi,[1 0 1]);   %A18(:,:,1) = A17(:,:,2);   A18(1,:,:) = A17(2,:,:);    %phi(i-1,j,k-1)

C1 = (A1-phi)./sqrt((A1-phi).^2 + ((A7+A3-A8-A4)/4)  .^2 + ((A9+A5-A10-A6)/4) .^2 +eps);
C2 = (phi-A2)./sqrt((phi-A2).^2 + ((A3+A11-A4-A12)/4).^2 + ((A5+A13-A6-A14)/4).^2 +eps);
C3 = (A3-phi)./sqrt((A3-phi).^2 + ((A7+A1-A11-A2)/4) .^2 + ((A15+A5-A16-A6)/4).^2 +eps);
C4 = (phi-A4)./sqrt((phi-A4).^2 + ((A1+A8-A2-A12)/4) .^2 + ((A5+A17-A6-A18)/4).^2 +eps);
C5 = (A5-phi)./sqrt((A5-phi).^2 + ((A9+A1-A13-A2)/4) .^2 + ((A15+A3-A17-A4)/4).^2 +eps);
C6 = (phi-A6)./sqrt((phi-A6).^2 + ((A1+A10-A2-A14)/4).^2 + ((A3+A16-A4-A18)/4).^2+eps);

edgeTerm=g1.*C1 -g2.*C2 +g3.*C3 -g4.*C4 +g5.*C5 -g6.*C6;

function lapPhi=EXLaplace(phi)
A1 = circshift(phi,[0 -1 0]); %A1(:,n,:) = A1(:,n-1,:); %phi(i,j+1,k)
A2 = circshift(phi,[0 1 0]);  %A2(:,1,:) = A2(:,2,:);   %phi(i,j-1,k)
A3 = circshift(phi,[-1 0 0]); %A3(m,:,:) = A3(m-1,:,:); %phi(i+1,j,k)
A4 = circshift(phi,[1 0 0]);  %A4(1,:,:) = A4(2,:,:);   %phi(i-1,j,k)
A5 = circshift(phi,[0 0 -1]); %A5(:,:,k) = A5(:,:,k-1); %phi(i,j,k+1)
A6 = circshift(phi,[0 0 1]);  %A6(:,:,1) = A6(:,:,2);   %phi(i,j,k-1)
lapPhi = A1+A2+A3+A4+A5+A6-6*phi;

function [ux,uy,uz]=grad(phi)
A1 = circshift(phi,[0 -1 0]); %A1(:,n,:) = A1(:,n-1,:); %phi(i,j+1,k)
A2 = circshift(phi,[0 1 0]);  %A2(:,1,:) = A2(:,2,:);   %phi(i,j-1,k)
A3 = circshift(phi,[-1 0 0]); %A3(m,:,:) = A3(m-1,:,:); %phi(i+1,j,k)
A4 = circshift(phi,[1 0 0]);  %A4(1,:,:) = A4(2,:,:);   %phi(i-1,j,k)
A5 = circshift(phi,[0 0 -1]); %A5(:,:,k) = A5(:,:,k-1); %phi(i,j,k+1)
A6 = circshift(phi,[0 0 1]);  %A6(:,:,1) = A6(:,:,2);   %phi(i,j,k-1)
ux=(A1-A2)./2;
uy=(A3-A4)./2;
uz=(A5-A6)./2;

function div_v=Div(v1,v2,v3)
[m,n,z]=size(v1);
A1=circshift(v1,[0,-1,0]);A1(:,n,:)=A1(:,n-1,:);
A2=circshift(v1,[0,1,0]);A2(:,1,:)=A2(:,2,:);
A3=circshift(v2,[-1,0,0]);A3(m,:,:)=A3(m-1,:,:);
A4=circshift(v2,[1,0,0]);A4(1,:,:)=A4(2,:,:);
A5=circshift(v3,[0,0,-1]);A5(:,:,z)=A3(:,:,z-1);
A6=circshift(v3,[0,0,1]);A6(:,:,1)=A4(:,:,2);
vx=(A1-A2)./2;
vy=(A3-A4)./2;
vz=(A5-A6)./2;
div_v=vx+vy+vz;

function value=heaviside(phi,epsilon)
value=.5*(1+(2/pi)*atan(phi/epsilon));

function y=dirac(phi,epsilon)
y=epsilon./(epsilon^2+phi.^2)/pi;
