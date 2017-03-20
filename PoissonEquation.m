clear all; clc; close all;
%Setting up the number of XY grid points
Nx = 101; Ny = Nx; step = 1/Nx; H = 1/(step^2);
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

%Defining my Boundary Conditions and Initial Conditions
U = zeros(Ny,Nx);
U(:,1) = uw;
U(1,:) = fb;
U(Ny,:) = gb;
F = sin(pi.*(x-ax)/(bx-ax)).*cos((pi/2).*(2.*((y - ay)./(by-ay))+1));

%Calculating the values
dSum = 1; n = 0;
nin = 1;
while  dSum > 0.001
    sum1 =  sum(sum(U.^2));
    for ny = 2: Ny-1
        for nx = 2: Nx-1
            if nin == 1; U(ny,Nx) = U(ny,Nx-1); end; 
            U(ny,nx) = (0.25).*( U(ny,nx+1) + U(ny,nx-1) + U(ny+1,nx) + U(ny-1,nx) + (hx^2).*F(1,nx));
        end
    end
    sum2 =  sum(sum(U.^2));
    dSum = abs(sum2 - sum1);
    n = n+1;
end
[Uxx, Uyy] = gradient(U,hx,hy);
Uxx = -Uxx;  Uyy = -Uyy; 
%U = sqrt(Uxx.^2 + Uyy.^2); %Interesting results when commented out


figure(1) % VECTOR FIELD V 
%set(gcf,'units','normalized','position',[0.02 0.52 0.3 0.32]);
surf(xx,yy,U);
%surf(xx(1:5:end,1:5:end),yy(1:5:end,1:5:end),U(1:5:end,1:5:end));
xlabel('x  [m]'); ylabel('y  [m]'); zlabel('U');
title('Test','fontweight','normal');
set(gca,'fontsize',14);
rotate3d
box on
axis tight
h =  colorbar;
h.Label.String = 'U   [ U ]';
view(55,49);
%set(gca,'ZTick',[0 5 10]);
    