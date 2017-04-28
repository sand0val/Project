clear all; clc; close all;
%Setting up the number of XY grid points
Nx = 30  ; Ny = Nx; step = 1/Nx; H = 1/(step^2);
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
% U = zeros(Ny,Nx);
% U(1,:) = gb;
% U(Ny,:) = fb;
% U(:,1) = uw;
U = zeros(Ny,Nx);
U(:,Ny) = gb;
U(:,1) = fb;
U(1,:) = uw;
%May need to fix the forces to correspond with X and Y right on plot
F = zeros(Ny,Nx);
%row = fliplr(Nx:-1:1);
for i = 1:Nx;
    for j = 1:Ny;
        F(i,j) = sin(pi.*(x(i)-ax)/(bx-ax)).*cos((pi/2).*(2.*((y(j) - ay)./(by-ay))+1));
    end
end

%------------------------------------------------------
%Calculating the values (Top to bottom, left to right)
dSum = 1; n = 0;
nin = 1; non = 1;
while  dSum > 0.001
%     U(Nx,[2:Ny-1]) = (0.25)*(U(Nx,[2:Nx-1]-1)+U(Nx,[2:Ny-1]+1)+2*U(Ny-1,[2:Nx-1])+(hx^2)*F(Nx,[2:Ny-1]));
%     U(1,1)= (U(1,2)+U(2,1))/2;
%     U(1,Nx)= (U(1,Nx-1)+U(2,Nx))/2;
%     U(Ny,1)= (U(Ny-1,1)+U(2,Ny))/2;
%     U(Ny,Nx)= (U(Ny,Nx-1)+U(Ny-1,Nx))/2;
%     nin = nin+1; %switch
    sum1 =  sum(sum(U.^2));
    for i = 2: Nx-1;
        for j = 2: Ny-1;
            if nin == 1; %Loop that initializes all the remaining conditions and corner points 
                if non ==1; U(Nx,[2:Ny-1]) = (0.25)*(U(Nx,[2:Nx-1]-1)+U(Nx,[2:Ny-1]+1)+2*U(Ny-1,[2:Nx-1])+(hx^2)*F(Nx,[2:Ny-1]));
                non = non+1; end;
                U(1,1)= (U(1,2)+U(2,1))/2;
                U(1,Nx)= (U(1,Nx-1)+U(2,Nx))/2;
                U(Ny,1)= (U(Ny-1,1)+U(2,Ny))/2;
                U(Ny,Nx)= (U(Ny,Nx-1)+U(Ny-1,Nx))/2;
                nin = nin+1; %switch
            end; 
            U(i,j) = (0.25).*( U(i+1,j) + U(i-1,j) + U(i,j+1) + U(i,j-1) + (hx^2).*F(i,j));
            %U(j,i)= ( (hy^2)*(U(i,j-1)+U(i,j+1))+(hx^2)*(U(i-1,j)+U(i+1,j))+(hx^2)*(hy^2)*F(i,j) )/(2*((hx^2)+(hy^2)));
        end
    end
    sum2 =  sum(sum(U.^2));
    dSum = abs(sum2 - sum1);
    n = n + 1;
    nin = 1;
    %non = 1;
end
% [Uxx, Uyy] = gradient(U,hx,hy);
% Uxx = -Uxx;  Uyy = -Uyy; 
% U = sqrt(Uxx.^2 + Uyy.^2); %Interesting results when commented out


figure(1) % VECTOR FIELD V 
%set(gcf,'units','normalized','position',[0.02 0.52 0.3 0.32]);
surf(xx,yy,U');
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
    