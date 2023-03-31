function [u,E] = fastReinitalisation(f0, theta, iter)

padNum=10;
f0=padarray(f0,[padNum,padNum,padNum],'symmetric');
[m,n,l]=size(f0);

w1=zeros(m,n,l,'single');
w2=zeros(m,n,l,'single');
w3=zeros(m,n,l,'single');

[Y,X,Z]=meshgrid(0:n-1,0:m-1,0:l-1);
G=cos(2*pi*X/m)+cos(2*pi*Y/n)+cos(2*pi*Z/l)-3;
G=single(G);

tic
for step=1:iter
    step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    g=f0-theta*(Bx(w1)+By(w2)+Bz(w3));
    u=real(ifftn(fftn(g)./(1-2*theta*G)));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    w1=Fx(u); w2=Fy(u); w3=Fz(u);
    abs_w=max(sqrt(w1.^2+w2.^2+w3.^2),1);
    w1=w1./abs_w;
    w2=w2./abs_w;
    w3=w3./abs_w;
    
    E(step)=sum(sum(sum(theta*(u-f0).^2./2)));
end
toc
u=u(padNum+1:m-padNum,padNum+1:n-padNum,padNum+1:l-padNum);

% Forward derivative operator on x with periodic boundary condition
function Fxu = Fx(u)
Fxu = circshift(u,[0 -1 0])-u;

% Forward derivative operator on y with periodic boundary condition
function Fyu = Fy(u)
Fyu = circshift(u,[-1 0 0])-u;

% Forward derivative operator on z with periodic boundary condition
function Fzu = Fz(u)
Fzu = circshift(u,[0 0 -1])-u;

% Backward derivative operator on x with periodic boundary condition
function Bxu = Bx(u)
Bxu = u-circshift(u,[0 1 0]);

% Backward derivative operator on y with periodic boundary condition
function Byu = By(u)
Byu = u-circshift(u,[1 0 0]);

% Backward derivative operator on z with periodic boundary condition
function Bzu = Bz(u)
Bzu = u-circshift(u,[0 0 1]);