%Johnathan Sandoval
%Computing For Mechanical Engineerers
%4/28/2017
%Solution For Code
clear all; close all;clc
lock = 0; flag = 0;
while lock ~= 1;
    while flag ~= [1:3]
        disp('Please select the case you would like to see!')
        disp('Case #1: Gauss Seidel solution')
        disp('Case #2: SOR Method')
        disp('Case #3: Both Solutions Compared')
        flag = input('Case #');
        if flag ~= [1:3];
            disp(' ')
            disp('HAHA okay that was funny, now seriously....select a case this time!')
            disp(' ')
        end
    end
    disp(' ')
    disp('Please input the amount of steps you would prefer. Steps are proportional in the X and Y direction')
    disp('Recommendation: 101 is best for a fast but ideal solution')
    Nx = input('Selection:');
    disp(' ')
    lock = 1;
    disp('Loading....')
end

%%
switch flag
%Gauss Seidel Code
    case 1
        t = 1;
        while t ~= 3;
            %Setting up the number of XY grid points
            Nx = Nx  ; Ny = Nx; step = 1/Nx; H = 1/(step^2);
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
            if t == 1;
                F = zeros(Ny,Nx);
                for i = 1:Nx;
                    for j = 1:Ny;
                        F(i,j) = sin(pi.*(x(i)-ax)/(bx-ax)).*cos((pi/2).*(2.*((y(j) - ay)./(by-ay))+1));
                    end
                end
            else
                F = zeros(Ny,Nx);
            end
            %Round two at loops
            bound = 1; n = 0;
            e = 1;
            while  e > 10^-3;
                Up = U;
                for i = 2:Nx-1;
                    for j = 2:Ny-1;
                        if bound == 1;
                            U(2:Ny-1,Nx) = (1/4)*(2*U(2:Ny-1,Ny-1)+U([2:Ny-1]-1,Ny)+U((2:Ny-1)+1,Ny)+(hx^2)*F((2:Ny-1),Ny));
                            bound = bound +1 ;
                            U(1,1)= (U(1,2)+U(2,1))/2;
                            U(1,Nx)= (U(1,Nx-1)+U(2,Nx))/2;
                            U(Ny,1)= (U(Ny-1,1)+U(Ny,2))/2;
                            U(Ny,Nx)= (U(Ny,Nx-1)+U(Ny-1,Nx))/2;
                            Up = U;
                        end;
                        %                    U(i,j) = 0.25*( Up(i+1,j)+U(i-1,j)+ Up(i,j+1)+ U(i,j-1) + (1^2)*F(i,j));
                        U(i,j)= ( (hy^2)*(U(i,j-1)+U(i,j+1))+(hx^2)*(U(i-1,j)+U(i+1,j))+(hx^2)*(hy^2)*F(i,j) )/(2*((hx^2)+(hy^2)));
                    end
                end
                E = U - Up;
                %e = mean(mean(E(2:Nx-1,2:Nx-1).^2));
                e = mean(mean(E(1:Nx,1:Nx).^2));
                %surf(xx,yy,U)
                %disp(e);
                %bound = 1;
                n = n+1;
                [Uxx, Uyy] = gradient(U,hx,hy);
                Uxx = -Uxx;  Uyy = -Uyy;
                Ut = sqrt(Uxx.^2 + Uyy.^2);
            end
            
            if t == 1;
                disp(['The error is ',num2str(e)])
                disp(['Using the Gauss Seidel Method, this problem took  ',num2str(n),' iterations to solve.'])
                disp(' ')
                
                figure(1) % VECTOR FIELD U
                set(gcf,'units','normalized','position',[0.02 0.52 0.3 0.32]);
                surf(xx,yy,U');
                %surf(xx(1:5:end,1:5:end),yy(1:5:end,1:5:end),U(1:5:end,1:5:end));
                xlabel('x'); ylabel('y'); zlabel('U');
                title('Solution using Gauss Seidel Method','fontweight','normal');
                set(gca,'fontsize',14);
                rotate3d
                box on
                axis tight
                h =  colorbar;
                h.Label.String = 'U';
                view(55,49);
                %set(gca,'ZTick',[0 5 10]);
                
                
                figure(2)
                set(gcf,'units','normalized','position',[0.33 0.52 0.3 0.32]);
                contourf(xx,yy,U,16);
                pcolor(xx,yy,U);
                shading interp
                xlabel('x'); ylabel('y');
                title('Heatmap of solution','fontweight','normal');
                set(gca,'fontsize',14)
                box on
                h =  colorbar;
                h.Label.String = 'U';
                %set(gca,'xLim',[-3.14,3.14]); set(gca,'yLim', [-3.14, 3.14]);
                %set(gca,'xTick',-5:5); set(gca,'yTick', -5:5);
                axis square
                box on
                
                
                figure(3)
                set(gcf,'units','normalized','position',[0.65 0.52 0.3 0.32])
                contourf(xx,yy,Ut,20);
                shading interp
                xlabel('x'); ylabel('y');
                title('Gradient solution','fontweight','normal');
                set(gca,'fontsize',14)
                box on
                h =  colorbar;
                h.Label.String = '| Up |';
                axis equal
            else
                disp(['The error is when F=0 is ',num2str(e)])
                disp(['Using the Gauss Seidel Method, this problem took  ',num2str(n),' iterations to solve.'])
                disp(' ')
                
                figure(4) % VECTOR FIELD U
                set(gcf,'units','normalized','position',[0.02 0.1 0.3 0.32]);
                surf(xx,yy,U');
                %surf(xx(1:5:end,1:5:end),yy(1:5:end,1:5:end),U(1:5:end,1:5:end));
                xlabel('x'); ylabel('y'); zlabel('U');
                title('Solution using Gauss-Seidel Method when F = 0 (Laplace Eqn)','fontweight','normal');
                set(gca,'fontsize',14);
                rotate3d
                box on
                axis tight
                h =  colorbar;
                h.Label.String = 'U';
                view(55,49);
                %set(gca,'ZTick',[0 5 10]);
                
                
                figure(5)
                set(gcf,'units','normalized','position',[0.33 0.1 0.3 0.32]);
                contourf(xx,yy,U,16);
                pcolor(xx,yy,U);
                shading interp
                xlabel('x'); ylabel('y');
                title('Heatmap of solution','fontweight','normal');
                set(gca,'fontsize',14)
                box on
                h =  colorbar;
                h.Label.String = 'U';
                %set(gca,'xLim',[-3.14,3.14]); set(gca,'yLim', [-3.14, 3.14]);
                %set(gca,'xTick',-5:5); set(gca,'yTick', -5:5);
                axis square
                box on
                
                
                figure(6)
                set(gcf,'units','normalized','position',[0.65 0.1 0.3 0.32]);
                contourf(xx,yy,Ut,20);
                shading interp
                xlabel('x'); ylabel('y');
                title('Gradient solution','fontweight','normal');
                set(gca,'fontsize',14)
                box on
                h =  colorbar;
                h.Label.String = '| Up |';
                axis equal
            end
            t = t+1;
        end
%SOR Method 
    case 2    
        t = 1;
        while t ~= 3;
            %Setting up the number of XY grid points
            Nx = Nx  ; Ny = Nx; step = 1/Nx; H = 1/(step^2);
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
            if t == 1;
                F = zeros(Ny,Nx);
                for i = 1:Nx;
                    for j = 1:Ny;
                        F(i,j) = sin(pi.*(x(i)-ax)/(bx-ax)).*cos((pi/2).*(2.*((y(j) - ay)./(by-ay))+1));
                    end
                end
            else
                F = zeros(Ny,Nx);
            end
            
            %Round two at loops
            bound = 1; n = 0;
            e = 1;
            w = 2/(1+sin(pi*hx/(2*pi))); %it is less than 2 thus okay
            while  e > 10^-3;
                Up = U;
                for i = 2:Nx-1;
                    for j = 2:Ny-1;
                        if bound == 1;
                            U(2:Ny-1,Nx) = (1/4)*(2*U(2:Ny-1,Ny-1)+U([2:Ny-1]-1,Ny)+U((2:Ny-1)+1,Ny)+(hx^2)*F((2:Ny-1),Ny));
                            bound = bound +1 ;
                            
                            U(1,1)= (U(1,2)+U(2,1))/2;
                            U(1,Nx)= (U(1,Nx-1)+U(2,Nx))/2;
                            U(Ny,1)= (U(Ny-1,1)+U(Ny,2))/2;
                            U(Ny,Nx)= (U(Ny,Nx-1)+U(Ny-1,Nx))/2;
                            %                         U(1,1)= (U(1,2)+U(2,1))/2;
                            %                         U(1,Nx)= (U(1,Nx-1)+U(2,Nx))/2;
                            %                         U(Ny,1)= (U(Ny-1,1)+U(2,Ny))/2;
                            %                         U(Ny,Nx)= (U(Ny,Nx-1)+U(Ny-1,Nx))/2;
                        end;
                        U(i,j) = (1-w)*Up(i,j)+(w/4)*( Up(i+1,j)+U(i-1,j)+ Up(i,j+1)+ U(i,j-1) + (hx^2)*F(i,j));
                    end
                end
                E = U - Up;
                e = mean(mean(E(1:Nx,1:Nx).^2));
                n = n+1;
                [Uxx, Uyy] = gradient(U,hx,hy);
                Uxx = -Uxx;  Uyy = -Uyy;
                Ut = sqrt(Uxx.^2 + Uyy.^2);
            end
            if t == 1;
                disp(['The error is ',num2str(e)])
                disp(['Using the SOR Method, this problem took  ',num2str(n),' iterations to solve.'])
                disp(' ')
                
                figure(1) % VECTOR FIELD U
                set(gcf,'units','normalized','position',[0.02 0.52 0.3 0.32]);
                surf(xx,yy,U');
                %surf(xx(1:5:end,1:5:end),yy(1:5:end,1:5:end),U(1:5:end,1:5:end));
                xlabel('x'); ylabel('y'); zlabel('U');
                title('Solution using SOR Method','fontweight','normal');
                set(gca,'fontsize',14);
                rotate3d
                box on
                axis tight
                h =  colorbar;
                h.Label.String = 'U';
                view(55,49);
                %set(gca,'ZTick',[0 5 10]);
                
                
                figure(2)
                set(gcf,'units','normalized','position',[0.33 0.52 0.3 0.32]);
                contourf(xx,yy,U,16);
                pcolor(xx,yy,U);
                shading interp
                xlabel('x'); ylabel('y');
                title('Heatmap of solution','fontweight','normal');
                set(gca,'fontsize',14)
                box on
                h =  colorbar;
                h.Label.String = 'U';
                %set(gca,'xLim',[-3.14,3.14]); set(gca,'yLim', [-3.14, 3.14]);
                %set(gca,'xTick',-5:5); set(gca,'yTick', -5:5);
                axis square
                box on
                
                
                figure(3)
                set(gcf,'units','normalized','position',[0.65 0.52 0.3 0.32])
                contourf(xx,yy,Ut,20);
                shading interp
                xlabel('x'); ylabel('y');
                title('Gradient solution','fontweight','normal');
                set(gca,'fontsize',14)
                box on
                h =  colorbar;
                h.Label.String = '| Up |';
                axis equal
            else
                disp(['The error when F = 0 is ',num2str(e)])
                disp(['Using the SOR Method when, this problem took  ',num2str(n),' iterations to solve.'])
                disp(' ')
                
                figure(4) % VECTOR FIELD U
                set(gcf,'units','normalized','position',[0.02 0.1 0.3 0.32]);
                surf(xx,yy,U');
                %surf(xx(1:5:end,1:5:end),yy(1:5:end,1:5:end),U(1:5:end,1:5:end));
                xlabel('x'); ylabel('y'); zlabel('U');
                title('Solution using SOR Method when F = 0 (Laplace Eqn)','fontweight','normal');
                set(gca,'fontsize',14);
                rotate3d
                box on
                axis tight
                h =  colorbar;
                h.Label.String = 'U';
                view(55,49);
                %set(gca,'ZTick',[0 5 10]);
                
                
                figure(5)
                set(gcf,'units','normalized','position',[0.33 0.1 0.3 0.32]);
                contourf(xx,yy,U,16);
                pcolor(xx,yy,U);
                shading interp
                xlabel('x'); ylabel('y');
                title('Heatmap of solution','fontweight','normal');
                set(gca,'fontsize',14)
                box on
                h =  colorbar;
                h.Label.String = 'U';
                %set(gca,'xLim',[-3.14,3.14]); set(gca,'yLim', [-3.14, 3.14]);
                %set(gca,'xTick',-5:5); set(gca,'yTick', -5:5);
                axis square
                box on
                
                
                figure(6)
                set(gcf,'units','normalized','position',[0.65 0.1 0.3 0.32]);
                contourf(xx,yy,Ut,20);
                shading interp
                xlabel('x'); ylabel('y');
                title('Gradient solution','fontweight','normal');
                set(gca,'fontsize',14)
                box on
                h =  colorbar;
                h.Label.String = '| Up |';
                axis equal
            end
            t = t+1;
            
        end

 %%       
    case 3
        %Setting up the number of XY grid points
        Nx = Nx  ; Ny = Nx; step = 1/Nx; H = 1/(step^2);
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
        for i = 1:Nx;
            for j = 1:Ny;
                F(i,j) = sin(pi.*(x(i)-ax)/(bx-ax)).*cos((pi/2).*(2.*((y(j) - ay)./(by-ay))+1));
            end
        end
        
        %Round two at loops
        bound = 1; n = 0;
        e = 1;
        while  e > 10^-3;
            Up = U;
            for i = 2:Nx-1;
                for j = 2:Ny-1;
                    if bound == 1;
                        U(2:Ny-1,Nx) = (1/4)*(2*U(2:Ny-1,Ny-1)+U([2:Ny-1]-1,Ny)+U((2:Ny-1)+1,Ny)+(hx^2)*F((2:Ny-1),Ny));
                        bound = bound +1 ;
                        U(1,1)= (U(1,2)+U(2,1))/2;
                        U(1,Nx)= (U(1,Nx-1)+U(2,Nx))/2;
                        U(Ny,1)= (U(Ny-1,1)+U(Ny,2))/2;
                        U(Ny,Nx)= (U(Ny,Nx-1)+U(Ny-1,Nx))/2;
                    end;
                    U(i,j) = 0.25*( U(i+1,j)+U(i-1,j)+ U(i,j+1)+ U(i,j-1) + (hx^2)*F(i,j));
                end
            end
            E = U - Up;
            e = mean(mean(E(1:Nx,1:Nx).^2));
            %surf(xx,yy,U)
            %disp(e);
            %bound = 1;
            n = n+1;
            [Uxx, Uyy] = gradient(U,hx,hy);
            Uxx = -Uxx;  Uyy = -Uyy;
            Ut = sqrt(Uxx.^2 + Uyy.^2);
        end
        disp(['The error is ',num2str(e)])
        disp(['Using the Gauss Seidel Method, this problem took  ',num2str(n),' iterations to solve.'])
        disp(' ')
        
        figure(1) % VECTOR FIELD U
        set(gcf,'units','normalized','position',[0.02 0.52 0.3 0.32]);
        surf(xx,yy,U');
        %surf(xx(1:5:end,1:5:end),yy(1:5:end,1:5:end),U(1:5:end,1:5:end));
        xlabel('x  [m]'); ylabel('y  [m]'); zlabel('U');
        title('Solution using Gauss Seidel Method','fontweight','normal');
        set(gca,'fontsize',14);
        rotate3d
        box on
        axis tight
        h =  colorbar;
        h.Label.String = 'U   [ U ]';
        view(55,49);
        %set(gca,'ZTick',[0 5 10]);
        
        
        figure(2)
        set(gcf,'units','normalized','position',[0.33 0.52 0.3 0.32]);
        contourf(xx,yy,U,16);
        %pcolor(xx,yy,U);
        shading interp
        xlabel('x  [m]'); ylabel('y  [m]');
        title('Heatmap of solution','fontweight','normal');
        set(gca,'fontsize',14)
        box on
        h =  colorbar;
        h.Label.String = 'V   [ V ]';
        %set(gca,'xLim',[-5,5]); set(gca,'yLim', [-5, 5]);
        %set(gca,'xTick',-5:5); set(gca,'yTick', -5:5);
        axis square
        box on
        
        
        figure(3)
        set(gcf,'units','normalized','position',[0.65 0.52 0.3 0.32])
        contourf(xx,yy,Ut,20);
        shading interp
        xlabel('x  [m]'); ylabel('y  [m]');
        title('Gradient solution','fontweight','normal');
        set(gca,'fontsize',14)
        box on
        h =  colorbar;
        h.Label.String = '| E |   [ V/m ]';
        axis equal
        
        %SOR Method
        
        %clear all; close all;
        %Setting up the number of XY grid points
        Nx = Nx  ; Ny = Nx; step = 1/Nx; H = 1/(step^2);
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
        for i = 1:Nx;
            for j = 1:Ny;
                F(i,j) = sin(pi.*(x(i)-ax)/(bx-ax)).*cos((pi/2).*(2.*((y(j) - ay)./(by-ay))+1));
            end
        end
        
        %Round two at loops
        bound = 1; n1 = 0;
        e = 1;
        w = 2/(1+sin(pi*hx/(2*pi)));
        while  e > 10^-3;
            Up = U;
            for i = 2:Nx-1;
                for j = 2:Ny-1;
                    if bound == 1;
                        U(2:Ny-1,Nx) = (1/4)*(2*U(2:Ny-1,Ny-1)+U([2:Ny-1]-1,Ny)+U((2:Ny-1)+1,Ny)+(hx^2)*F((2:Ny-1),Ny));
                        bound = bound +1 ;
                        U(1,1)= (U(1,2)+U(2,1))/2;
                        U(1,Nx)= (U(1,Nx-1)+U(2,Nx))/2;
                        U(Ny,1)= (U(Ny-1,1)+U(Ny,2))/2;
                        U(Ny,Nx)= (U(Ny,Nx-1)+U(Ny-1,Nx))/2;
                    end;
                    U(i,j) = (1-w)*Up(i,j)+(w/4)*( Up(i+1,j)+U(i-1,j)+ Up(i,j+1)+ U(i,j-1) + (hx^2)*F(i,j));
                end
            end
            E1 = U - Up;
            e = mean(mean(E1(1:Nx,1:Nx).^2));
            n1 = n1+1;
            [Uxx, Uyy] = gradient(U,hx,hy);
            Uxx = -Uxx;  Uyy = -Uyy;
            Ut = sqrt(Uxx.^2 + Uyy.^2);
        end
        disp(['The error is ',num2str(e)])
        disp(['Using the SOR Method, this problem took  ',num2str(n1),' iterations to solve.'])
        disp(' ')
        
        
        figure(4) % VECTOR FIELD U
        set(gcf,'units','normalized','position',[0.02 0.1 0.3 0.32]);
        surf(xx,yy,U');
        xlabel('x  [m]'); ylabel('y  [m]'); zlabel('U');
        title('Solution using SOR Method','fontweight','normal');
        set(gca,'fontsize',14);
        rotate3d
        box on
        axis tight
        h =  colorbar;
        h.Label.String = 'U   [ U ]';
        view(55,49);
        
        
        figure(5)
        set(gcf,'units','normalized','position',[0.33 0.1 0.3 0.32]);
        contourf(xx,yy,U,16);
        pcolor(xx,yy,U);
        shading interp
        xlabel('x  [m]'); ylabel('y  [m]');
        title('Heatmap of solution','fontweight','normal');
        set(gca,'fontsize',14)
        box on
        h =  colorbar;
        h.Label.String = 'U   [ U ]';
        axis square
        box on
        
        
        figure(6)
        set(gcf,'units','normalized','position',[0.65 0.1 0.3 0.32]);
        contourf(xx,yy,Ut,20);
        shading interp
        xlabel('x  [m]'); ylabel('y  [m]');
        title('Gradient solution','fontweight','normal');
        set(gca,'fontsize',14)
        box on
        h =  colorbar;
        h.Label.String = '| Up |';
        axis equal
        z = n-n1;
        disp('Based off these values....')
        disp('')
        disp(['The SOR method took ',num2str(z),' less iterations to solve'])
        
end
