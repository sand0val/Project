clear all; clc; close all;
%Setting up the number of XY grid points
Nx = 101  ; Ny = Nx; step = 1/Nx; H = 1/(step^2);
%Setting up length of X and Y regions
ax = -pi; bx = pi; ay = ax; by = bx;
Lx = 2*pi; Ly = 2*pi;

%Setting up XY coordinates
minX1 = ax; maxX1 = bx; 
minY1 = ay; maxY1 = by;
x = linspace(minX1,maxX1,Nx);
y = linspace(minY1,maxY1,Ny);
[xx,yy] = meshgrid(x,y);
hx = x(2) - x(1); hy = y(2) - y(1);
uw = (((bx-ax).^2).*cos(pi.*ax/bx))+(((y-ay)/(by-ay))*((ax.*(bx - ax).^2) - (((bx-ax).^2).*cos(pi.*ax/bx))));
gb = ((bx-x).^2).*cos(pi.*x/bx);fb = (x.*(bx - x).^2);

%Defining my Boundary Conditions and Initial Conditions (Have 2 sets)
U = zeros(Ny,Nx);
U(1,2:Nx-1) = gb(2:Nx-1);
U(Ny,2:Nx-1) = fb(2:Nx-1);
U(2:Nx-1,1) = uw(2:Nx-1);

%May need to fix the forces to correspond with X and Y right on plot
F = zeros(Ny,Nx);
%row = fliplr(Nx:-1:1);
for i = 1:Nx;
    for j = 1:Ny;
        F(i,j) = sin(pi.*(x(i)-ax)/(bx-ax)).*cos((pi/2).*(2.*((y(j) - ay)./(by-ay))+1));
    end
end

%Round two at loops
bound = 1; n = 0;
e = 1;
w = 4/3; 
while  e > 10^-3;
    Up = U;
    for i = 2:Nx-1;
        for j = 2:Ny-1;
            if bound == 1;
                U(2:Ny-1,Nx) = (1/4)*(2*U(2:Ny-1,Ny-1)+U([2:Ny-1]-1,Ny)+U((2:Ny-1)+1,Ny)+(hx^2)*F((2:Ny-1),Ny));
                bound = bound +1 ;
                U(1,1)= (U(1,2)+U(2,1))/2;
                U(1,Nx)= (U(1,Nx-1)+U(2,Nx))/2;
                U(Ny,1)= (U(Ny-1,1)+U(2,Ny))/2;
                U(Ny,Nx)= (U(Ny,Nx-1)+U(Ny-1,Nx))/2;
            end;
            U(i,j) = (1-w)*Up(i,j)+(w/4)*( Up(i+1,j)+U(i-1,j)+ Up(i,j+1)+ U(i,j-1) + F(i,j));
        end
    end
    E = U - Up;
    e = mean(mean(E(1:Nx,1:Nx).^2));
    %bound = 1;
    n = n+1;
    [Uxx, Uyy] = gradient(U,hx,hy);
    Uxx = -Uxx;  Uyy = -Uyy;
    Ut = sqrt(Uxx.^2 + Uyy.^2);
end


figure(1) % VECTOR FIELD U 
%set(gcf,'units','normalized','position',[0.02 0.52 0.3 0.32]);
surf(xx,yy,U');
xlabel('x  [m]'); ylabel('y  [m]'); zlabel('U');
title('Test','fontweight','normal');
set(gca,'fontsize',14);
rotate3d
box on
axis tight
h =  colorbar;
h.Label.String = 'U   [ U ]';
view(55,49);


figure(5)
set(gcf,'units','normalized','position',[0.33 0.52 0.3 0.32]);
contourf(xx,yy,U,16);
%pcolor(xx,yy,U);
shading interp
xlabel('x  [m]'); ylabel('y  [m]');
title('potential','fontweight','normal');
set(gca,'fontsize',14)
box on
h =  colorbar;
h.Label.String = 'U   [ U ]';
axis square
box on


figure(6)
set(gcf,'units','normalized','position',[0.33 0.1 0.3 0.32]);
contourf(xx,yy,Ut,20);
shading interp
xlabel('x  [m]'); ylabel('y  [m]');
title('electric field','fontweight','normal');
set(gca,'fontsize',14)
box on
h =  colorbar;
h.Label.String = '| Up |';
axis equal