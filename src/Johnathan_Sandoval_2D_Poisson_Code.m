%Johnathan Sandoval
%Student ID# 1186823
%Computing For Mechanical Engineers
%4/28/2017
%Solution for Poisson Equation, version "APc2-1"


%Notes before running this code!
%This code simulates the 2 required linear solving methods of Gauss-Seidel
%and Successive Over Relaxation.
%NOW with built in CHECK POINTING

%----------------------------------------------------------------------------
%When initiated you will be asked for the method you intend on
%using/seeing.

% 1) The first method is the Gauss-Seidel with the solution for
%when F = a value and also for when F = 0 (Laplace equation)
% 2) The second method uses the SOR method with a carefully chosen
% over-relaxation factor (w) and also for when F = 0
% 3) The third method takes both methods and compares them with their own plots so
%that they can be compared for differences.
% 4) The fourth method allows you to reload a checkpoint in the event that
% you were not able to finish your calculation. All progress is saved when half
% of the code is done. If the code executes completely then the checkpoint
% is also conviniently deleted so you don't have to worry about that.
% 5) The final choice allows the user to simply exit the program without
% issue. Added simply to make the life of the grader easier. Thank me
% later.

%After selecting the method you will be asked for the amount of steps you
%would like to use (Value of N, note that its proportional in both X and Y directions
% so choose the amount steps you feel your computer can safely handle).

%Lastly you will be asked for the type of error resolution you are aiming to
%achieve. The error resolution is scaled as 10^-X, where X is a number between 1 to infinity.
%Again, this value can be taxing on certain computer setups so use high
%values with caution.

%------------------------------------------------------------------------------------------

%This is a fairly lengthy code but the features that it offers make it
%one of the most versatile codes available for this type of problem/project. There
%is no need to change any values inside of this code.

%------------------------------------------------------------------------------------------


clear all; clc;
lock = 0; flag = 0; main = 0;marker  = 0; marker1 = 0; t = [];
while main ~= 1;
    while lock ~= 1;
        while flag ~= [1:5]
            disp('Please select one of the five cases you would like to try:')
            disp('Case #1: Gauss-Seidel solution')
            disp('Case #2: SOR Method solution')
            disp('Case #3: Both Solutions Compared')
            disp('Case #4: Reload Checkpoint')
            disp('Case #5: Exit program')
            disp(' ')
            flag = input('Case #'); checkp = flag;
            if flag ~= [1:5];
                disp(' ')
                disp('HAHA okay that was funny, now seriously....select a case this time!')
                disp(' ')
            elseif flag == 4;
                break
            elseif flag == 5;
                return
            end
        end
        if checkp == 4;
            if exist( 'checkpoint.mat','file' )
                if t == 2;
                    t = 2;
                elseif t == 3;
                    t = 2;
                else
                    t = 1;
                end
                fprintf('Checkpoint file found - Loading\n'); disp(' ');
                load('checkpoint.mat')
            else 
                disp('There is no checkpoint file here, next time look for "Saving checkpoint"');
                return
            end
        else
            t = 1;
            disp(' ')
            disp('Please input the value of N you would like to use. N is proportional in the X and Y direction')
            disp('Recommendation: 200 is best for a fast but ideal solution')
            disp(' ')
            Nx = abs(input('Value of N = '));
            disp(' ')
            disp(' ')
            disp('Lastly, Please input the 10^-X accuracy youd like to reach')
            disp('Recommendation: 7 at minimum for a consistent solution. Raise at your own risk')
            disp(' ')
            X = abs(input('Selection: '));
            disp(' ')
            disp('Loading....')
            disp(' ')
            disp(' ')
            lock = 1;
        end
    end

    XX = (10^-X);
    
    %Setting up the number of XY grid points
    Ny = Nx; step = 1/Nx; H = 1/(step^2);
    %Setting up length of X and Y regions
    ax = -pi; bx = pi; ay = ax; by = bx;
    Lx = 2*pi; Ly = 2*pi;
    
    %Setting up XY coordinates
    minX1 = ax; maxX1 = bx;
    minY1 = ay; maxY1 = by;
    x = linspace(minX1,maxX1,Nx);
    y = linspace(minY1,maxY1,Ny);
    [xx,yy] = meshgrid(x,y);
    yy = flipud(yy);
    hx = x(2) - x(1); hy = y(2) - y(1);
    uw = (((bx-ax).^2).*cos(pi.*ax/bx))+(((yy-ay)/(by-ay))*((ax.*(bx - ax).^2) - (((bx-ax).^2).*cos(pi.*ax/bx))));
    gb = ((bx-xx).^2).*cos(pi.*xx/bx);fb = (xx.*(bx - xx).^2);
    
    %Defining my Boundary Conditions and Initial Conditions (Have 2 sets)
    U = zeros(Ny,Nx);
    U(1,2:Nx-1) = gb(1,2:Nx-1);
    U(Ny,2:Nx-1) = fb(Ny,2:Nx-1);
    U(2:Nx-1,1) = uw(2:Nx-1,1);
    
    %%
    switch flag
        %Gauss-Seidel Code
        case 1
            while t ~= 3;
                if t == 1;
                    F = sin(pi.*(xx-ax)/(bx-ax)).*cos((pi/2).*(2.*((yy - ay)./(by - ay))+1));
                elseif t == 2;
                    F = zeros(Ny,Nx);
                    U = zeros(Ny,Nx);
                    U(1,2:Nx-1) = gb(1,2:Nx-1);
                    U(Ny,2:Nx-1) = fb(Ny,2:Nx-1);
                    U(2:Nx-1,1) = uw(2:Nx-1,1);
                end          
                %Main loop that evaluates the remaining Neumann BC and also deals with the corners 
                %by taking a mean of the two boxes that touch that the corner.
                bound = 1; n = 0;
                e = 1; Q = (69.833*Nx - 1164.3)/3; tic;
                while  e > XX;
                    Up = U;
                    for i = 2:Nx-1;
                        for j = 2:Ny-1;
                            if bound == 1;
                                bound = bound +1 ;
                                U(1,1)= (U(1,2)+U(2,1))/2; %These are used to even out my irregular spikes on the corners
                                U(1,Nx)= (U(1,Nx-1)+U(2,Nx))/2;
                                U(Ny,1)= (U(Ny-1,1)+U(Ny,2))/2;
                                U(Ny,Nx)= (U(Ny,Nx-1)+U(Ny-1,Nx))/2;
                                Up = U;
                            end;
                            if j == Ny-1;
                                U(i,Nx) = (1/4)*(2*U(i,Ny-1)+U(i-1,Ny)+U(i+1,Ny)+(hx^2)*F(i,Ny));
                            end
                            U(i,j) = (0.25)*(Up(i+1,j)+U(i-1,j)+Up(i,j+1)+U(i,j-1)+(hx^2)*F(i,j));
                        end
                    end
                    U(1,1)= (U(1,2)+U(2,1))/2;
                    U(1,Nx)= (U(1,Nx-1)+U(2,Nx))/2;
                    U(Ny,1)= (U(Ny-1,1)+U(Ny,2))/2;
                    U(Ny,Nx)= (U(Ny,Nx-1)+U(Ny-1,Nx))/2;
                    E = U - Up;
                    e = mean(mean(E(2:Nx-1,2:Nx-1).^2));
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
                    %surf(xx,yy,U');
                    mesh(xx,yy,U);
                    xlabel('x'); ylabel('y'); zlabel('U');
                    title('3D Solution using Gauss Seidel Method','fontweight','normal');
                    set(gca,'fontsize',14);
                    rotate3d
                    box on
                    axis tight
                    h =  colorbar;
                    h.Label.String = 'U';
                    view(55,49);
                    
                    
                    figure(2)
                    set(gcf,'units','normalized','position',[0.33 0.52 0.3 0.32]);
                    contourf(xx,yy,U,16);
                    pcolor(xx,yy,U);
                    shading interp
                    xlabel('x'); ylabel('y');
                    title('Heatmap of GS solution','fontweight','normal');
                    set(gca,'fontsize',14)
                    box on
                    h =  colorbar;
                    h.Label.String = 'U';
                    axis square
                    box on
                    
                    figure(3)
                    set(gcf,'units','normalized','position',[0.65 0.52 0.3 0.32])
                    contourf(xx,yy,Ut,20);
                    shading interp
                    xlabel('x'); ylabel('y');
                    title('Gradient of U using GS','fontweight','normal');
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
                    %surf(xx,yy,U');
                    mesh(xx,yy,U)
                    xlabel('x'); ylabel('y'); zlabel('U');
                    title('Solution using Gauss-Seidel Method when F = 0 (Laplace Eqn)','fontweight','normal');
                    set(gca,'fontsize',14);
                    rotate3d
                    box on
                    axis tight
                    h =  colorbar;
                    h.Label.String = 'U';
                    view(55,49);
                    
                    figure(5)
                    set(gcf,'units','normalized','position',[0.33 0.1 0.3 0.32]);
                    contourf(xx,yy,U,16);
                    %pcolor(xx,yy,U);
                    shading interp
                    xlabel('x'); ylabel('y');
                    title('Heatmap of GS solution','fontweight','normal');
                    set(gca,'fontsize',14)
                    box on
                    h =  colorbar;
                    h.Label.String = 'U';
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
                fprintf('Saving checkpoint');
                t = t+1;
                save('checkpoint.mat');
                disp(' ')
                disp(' ')
            end
            disp(' ')
            disp('Would you like to select another case or quit!?')
            disp(' # 1) to quit')
            disp(' # 2) to select another case')
            zz = abs(input('Selection # '));
            if zz == 1;
                main = 1;
                break
            elseif zz == 2;
                lock = 0;
                flag = 0;
                marker1 = 0;   %changed
                marker = 1;
                disp(' ')
                disp(' ')
                disp(' ')
                disp(' ')
                t = 0;
            end
            
        case 2
            
            while t ~= 3; 
                if t == 1;
                    F = sin(pi.*(xx-ax)/(bx-ax)).*cos((pi/2).*(2.*((yy - ay)./(by - ay))+1));
                elseif t == 2;
                    F = zeros(Ny,Nx);                    
                    U = zeros(Ny,Nx);
                    U(1,2:Nx-1) = gb(1,2:Nx-1);
                    U(Ny,2:Nx-1) = fb(Ny,2:Nx-1);
                    U(2:Nx-1,1) = uw(2:Nx-1,1);
                end
                
                %Round two at loops
                bound = 1; n = 0;
                e = 1;
                w = 2/(1+sin(pi*hx/(2*pi))); %it is less than 2 thus okay
                while  e > XX;
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
                    E = U - Up;
                    e = mean(mean(E(2:Nx-1,2:Nx-1).^2));
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
                    mesh(xx,yy,U);
                    xlabel('x'); ylabel('y'); zlabel('U');
                    title('Solution using SOR Method','fontweight','normal');
                    set(gca,'fontsize',14);
                    rotate3d
                    box on
                    axis tight
                    h =  colorbar;
                    h.Label.String = 'U';
                    view(55,49);
                    
                    figure(2)
                    set(gcf,'units','normalized','position',[0.33 0.52 0.3 0.32]);
                    contourf(xx,yy,U,16);
                    pcolor(xx,yy,U);
                    shading interp
                    xlabel('x'); ylabel('y');
                    title('Heatmap of SOR solution','fontweight','normal');
                    set(gca,'fontsize',14)
                    box on
                    h =  colorbar;
                    h.Label.String = 'U';
                    axis square
                    box on
                    
                    
                    figure(3)
                    set(gcf,'units','normalized','position',[0.65 0.52 0.3 0.32])
                    contourf(xx,yy,Ut,20);
                    shading interp
                    xlabel('x'); ylabel('y');
                    title('Gradient solutionusing SOR','fontweight','normal');
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
                    mesh(xx,yy,U);
                    xlabel('x'); ylabel('y'); zlabel('U');
                    title('Solution using SOR Method when F = 0 (Laplace Eqn)','fontweight','normal');
                    set(gca,'fontsize',14);
                    rotate3d
                    box on
                    axis tight
                    h =  colorbar;
                    h.Label.String = 'U';
                    view(55,49);
                    
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
                fprintf('Saving checkpoint');
                t = t+1;
                save('checkpoint.mat');
                disp(' ')
                disp(' ')
            end
            if marker1 ~= 1;
                disp(' ')
                disp('Would you like to select another case or quit?')
                disp(' # 1) to quit')
                disp(' # 2) to select another case')
                zz = abs(input('Selection # '));
                if zz == 1;
                    main = 1;
                    break;
                elseif zz == 2;
                    lock = 0;
                    flag = 0;
                    disp(' ')
                    disp(' ')
                    disp(' ')
                    disp(' ')
                    t = 0;
                    marker = 1;
                end
            elseif marker1 == 1;
                flag = 0; lock = 0;
                marker = 1;
                marker1 = 0;
            end
            
            %%
        case 3
 
            %Gauss-Seidel goes first
            F = zeros(Ny,Nx);
            F = sin(pi.*(xx-ax)/(bx-ax)).*cos((pi/2).*(2.*((yy - ay)./(by-ay))+1));
            
            %Round two at loops
            bound = 1; n = 0;
            e = 1;
            while  e > XX;
                Up = U;
                for i = 2:Nx-1;
                    for j = 2:Ny-1;
                        if bound == 1;
                            bound = bound +1 ;
                            U(1,1)= (U(1,2)+U(2,1))/2;
                            U(1,Nx)= (U(1,Nx-1)+U(2,Nx))/2;
                            U(Ny,1)= (U(Ny-1,1)+U(Ny,2))/2;
                            U(Ny,Nx)= (U(Ny,Nx-1)+U(Ny-1,Nx))/2;
                        end;
                        if j == Ny-1;
                            U(i,Nx) = (1/4)*(2*U(i,Ny-1)+U(i-1,Ny)+U(i+1,Ny)+(hx^2)*F(i,Ny));
                        end
                        U(i,j) = (0.25)*(Up(i+1,j)+U(i-1,j)+Up(i,j+1)+U(i,j-1)+(hx^2)*F(i,j));
                    end
                end
 
                E = U - Up;
                e = mean(mean(E(2:Nx,2:Nx).^2));
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
            mesh(xx,yy,U);
            xlabel('x '); ylabel('y '); zlabel('U');
            title('Solution using Gauss Seidel Method','fontweight','normal');
            set(gca,'fontsize',14);
            rotate3d
            box on
            axis tight
            h =  colorbar;
            h.Label.String = 'U ';
            view(55,49);
            
            
            figure(2)
            set(gcf,'units','normalized','position',[0.33 0.52 0.3 0.32]);
            contourf(xx,yy,U,16);
            shading interp
            xlabel('x '); ylabel('y ');
            title('Heatmap of solution','fontweight','normal');
            set(gca,'fontsize',14)
            box on
            h =  colorbar;
            h.Label.String = 'U';
            axis square
            box on
            
            
            figure(3)
            set(gcf,'units','normalized','position',[0.65 0.52 0.3 0.32])
            contourf(xx,yy,Ut,20);
            shading interp
            xlabel('x '); ylabel('y ');
            title('Gradient solution','fontweight','normal');
            set(gca,'fontsize',14)
            box on
            h =  colorbar;
            h.Label.String = 'U';
            axis equal
            
            fprintf('Saving checkpoint');
            t = t+1;
            save('checkpoint.mat');
            disp(' ')
            disp(' ')
            
            %SOR Method
            
            U = zeros(Ny,Nx);
            U(1,2:Nx-1) = gb(1,2:Nx-1);
            U(Ny,2:Nx-1) = fb(Ny,2:Nx-1);
            U(2:Nx-1,1) = uw(2:Nx-1,1);
            F = zeros(Ny,Nx);
            F = sin(pi.*(xx-ax)/(bx-ax)).*cos((pi/2).*(2.*((yy - ay)./(by-ay))+1));
            
            %Round two at loops
            bound = 1; n1 = 0;
            e = 1;
            w = 2/(1+sin(pi*hx/(2*pi)));
            while  e > XX;
                Up = U;
                for i = 2:Nx-1;
                    for j = 2:Ny-1;
                        if bound == 1;
                            bound = bound +1 ;
                            U(1,1)= (U(1,2)+U(2,1))/2;
                            U(1,Nx)= (U(1,Nx-1)+U(2,Nx))/2;
                            U(Ny,1)= (U(Ny-1,1)+U(Ny,2))/2;
                            U(Ny,Nx)= (U(Ny,Nx-1)+U(Ny-1,Nx))/2;
                        end;
                        if j == Ny-1;
                            U(i,Nx) = (1/4)*(2*U(i,Ny-1)+U(i-1,Ny)+U(i+1,Ny)+(hx^2)*F(i,Ny));
                        end
                        U(i,j) = (1-w)*Up(i,j)+(w/4)*( Up(i+1,j)+U(i-1,j)+ Up(i,j+1)+ U(i,j-1) + (hx^2)*F(i,j));
                    end
                end
                E1 = U - Up;
                e = mean(mean(E1(2:Nx,2:Nx).^2));
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
            mesh(xx,yy,U);
            xlabel('x '); ylabel('y '); zlabel('U');
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
            shading interp
            xlabel('x '); ylabel('y ');
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
            xlabel('x '); ylabel('y ');
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
            disp(' ')
    end
    if marker ~= 1;
        disp(' ')
        disp('Would you like to select another case or quit?')
        disp(' # 1) to quit')
        disp(' # 2) to select another case')
        zz = abs(input('Selection # '));
        if zz == 1;
            main = 1;
            break;
        elseif zz == 2;
            lock = 0;
            flag = 0;
            disp(' ')
            disp(' ')
            disp(' ')
            disp(' ')
        end
    elseif marker == 1;
        flag = 0; lock = 0; marker = 0;
    end
end
%So you don't have to delete it later when you see it on the desktop or
%wherever
delete checkpoint.mat

