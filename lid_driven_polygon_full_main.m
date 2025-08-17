clear all;
close all;
clc;

%% ME7751 Term Project Lid driven flow in a circular cavity using Lattice Boltzmann Method, ME Soohwan Kim
%  c7  c3   c6  D2Q9 model
%    \  |  /    
%  c4 -c1 - c2  
%    /  |  \    
%  c8  c5   c9      
% D2Q9 solver, SRT (Single Relaxation Time) model, BGK (Bhatnagar-Gross-Krook) collision operator
% based on Textbook (Mohamad 2011) with clearer distinction between the physical,
% nondimensional and numerical/discrete parameters, according to Jonas Latt's 2008 paper on the topic.

% add function that monitors the current flow field in real-time while MATLAB is solving for the flowfield

% define constants and initialization

L_p = 1; % cavity dimension 
U_p = 1; % cavity lid velocity

N = 50;
dt = 1;
Re = 500;

nu_p = L_p*U_p/Re; % physical kinematic viscosity
rho0 = 1;

% Derived discrete parameters
ds = L_p/(2*N);
u_lb = 0.1; % lattice speed or celerity ~ 0.1-0.2

nu_lb = 2*N*u_lb/Re; % lattice kinematic viscosity
c_s = 1/sqrt(3); % assume speed of sound 1/sqrt(3)
tau = 3*nu_lb + 0.5*dt; % dt = 1
omega = 1/tau; % relaxation parameter

% Lattice constants (D2Q9)
w = zeros(9,1); % weights
w(1) = 4/9;
w(2:5) = 1/9;
w(6:9) = 1/36;
c = zeros(9,2);
c(1,:) = [0, 0];
c(2,:) = [1, 0];
c(3,:) = [0, 1];
c(4,:) = [-1, 0];
c(5,:) = [0, -1];
c(6,:) = [1, 1];
c(7,:) = [-1, 1];
c(8,:) = [-1, -1];
c(9,:) = [1, -1];

% Initialize
rho = rho0*ones(2*N+1,2*N+1);
u = zeros(2*N+1,2*N+1);
v = zeros(2*N+1,2*N+1);
f = zeros(2*N+1,2*N+1,9);
feq = zeros(2*N+1,2*N+1,9);

u(end,2:end-1) = u_lb;

%% Define cavity
% uncomment the selected geometry 
% % sqaure
% isfluid = ones(2*N+1,2*N+1); % 1 = fluid, 0 = solid, 2 = boundary(if node is exactly on the boundary, it is 2)
% % the outer voxels are wall (except the moving top)
% isfluid(1,:) = 2;
% isfluid(:,1) = 2;
% isfluid(:,end) = 2;

%% half rectangle 
% isfluid = ones(2*N+1,2*N+1); % 1 = fluid, 0 = solid, 2 = boundary(if node is exactly on the boundary, it is 2)
% % the outer voxels are wall (except the moving top)
% isfluid(1:(N+1),:) = 2;
% isfluid(:,1) = 2;
% isfluid(:,end) = 2;

%% mountain
% isfluid = ones(2*N+1,2*N+1); % 1 = fluid, 0 = solid, 2 = boundary(if node is exactly on the boundary, it is 2)
% % the outer voxels are wall (except the moving top)
% isfluid(1,:) = 2;
% isfluid(:,1) = 2;
% isfluid(:,end) = 2;
% one_third = floor((2*N+1)/3);
% isfluid(1:1+one_third,1:1+one_third) = 2;
% isfluid(1:1+one_third,end-one_third:end) = 2;

%% mountain2, modified
% isfluid = ones(2*N+1,2*N+1); % 1 = fluid, 0 = solid, 2 = boundary(if node is exactly on the boundary, it is 2)
% % the outer voxels are wall (except the moving top)
% isfluid(1,:) = 2;
% isfluid(:,1) = 2;
% isfluid(:,end) = 2;
% one_third = floor((2*N+1)/3);
% isfluid(1:1+one_third,1:1+2*one_third) = 2;
% %isfluid(1:1+one_third,end-one_third:end) = 2;

%% mountain3, modified
% isfluid = ones(2*N+1,2*N+1); % 1 = fluid, 0 = solid, 2 = boundary(if node is exactly on the boundary, it is 2)
% % the outer voxels are wall (except the moving top)
% isfluid(1,:) = 2;
% isfluid(:,1) = 2;
% isfluid(:,end) = 2;
% one_third = floor((2*N+1)/3);
% isfluid(1:1+2*one_third,1:1+one_third) = 2;
% %isfluid(1:1+one_third,end-one_third:end) = 2;

%% half square
% % need high resolution
% isfluid = ones(2*N+1,2*N+1);
% isfluid(:,:) = 2;
% quarter = floor((2*N+1)/4);
% 
% isfluid(2*quarter:end,quarter+1:end-quarter) = 1;

%% Trapezoid
% isfluid = ones(2*N+1,2*N+1); % 1 = fluid, 0 = solid, 2 = boundary(if node is exactly on the boundary, it is 2)
% for i=1:2*N+1
%     for j=1:2*N+1
%         x = ds*(j-1);
%         y = ds*(i-1);
%         if or(y <= tand(60*2)*x+L_p , y <= tand(60)*(x-L_p)+L_p)
%             isfluid(i,j) = 2;
%         end
%     end
% end
% isfluid(1:(N+1),:) = 2;

%% Triangle
% isfluid = ones(2*N+1,2*N+1); % 1 = fluid, 0 = solid, 2 = boundary(if node is exactly on the boundary, it is 2)
% for i=1:2*N+1
%     for j=1:2*N+1
%         x = ds*(j-1);
%         y = ds*(i-1);
%         if or(y <= tand(60*2)*x+L_p , y <= tand(60)*(x-L_p)+L_p)
%             isfluid(i,j) = 2;
%         end
%     end
% end

%% pentagon
isfluid = ones(2*N+1,2*N+1); % 1 = fluid, 0 = solid, 2 = boundary(if node is exactly on the boundary, it is 2)

theta = pi/5; %36 deg
l = L_p/2/cos(theta);

for i=1:2*N+1
    for j=1:2*N+1
        x = ds*(j-1);
        y = ds*(i-1);
        if y <= tan(theta*4)*(x-l*cos(theta))
            isfluid(i,j) = 2;
        end
        if y >= tan(theta*3)*(x-L_p)+l*sin(theta)
            isfluid(i,j) = 2;
        end
        if y >= tan(theta*2)*x+l*sin(theta)
            isfluid(i,j) = 2;
        end
        if y <= tan(theta)*(x-l*cos(theta))
            isfluid(i,j) = 2;
        end
    end
end
%% Check moving lid start and end position, 2 1 1 1 1 1 1 2
lid_start = 1;
lid_end = 2*N+1;
mode = 0; % 1 if start identified, 2 if start and end identified
for j=2:2*N+1
switch mode
    case 0
        if isfluid(end,j) == 1
            lid_start = j;
            mode = 1;
        end
    case 1 
        if isfluid(end,j) == 2
            lid_end = j-1;
            mode = 2;
        end
    case 2
        if isfluid(end,j) == 1
            fprinf('Error: inappropriate lid definition - review the isfluid for potential mapping errors.');    
        end
end

end

%%
corr = 1;
tor = 1e-5;
itr = 0;

while (corr > tor)

    % Collision
    f = collide(f, u, v, rho, isfluid, omega, N);

    % Streaming
    f = stream(f, isfluid, N);

    % BC
    f = bounceback_full(f, isfluid, N);

    rho_b = f(end,lid_start:lid_end,1) + f(end,lid_start:lid_end,2) + f(end,lid_start:lid_end,4) + ... % moving wall
        2*(f(end,lid_start:lid_end,3) + f(end,lid_start:lid_end,7) + f(end,lid_start:lid_end,6));
    f(end,lid_start:lid_end,5) = f(end,lid_start:lid_end,3); 
    f(end,lid_start:lid_end,9) = f(end,lid_start:lid_end,7) + (u_lb / 6)*rho_b; 
    f(end,lid_start:lid_end,8) = f(end,lid_start:lid_end,6) - (u_lb / 6)*rho_b; 


    % Density and velocity reconstruction
    u_temp = u;
    v_temp = v;
    [u,v,rho] = reconstruct(f,u,v,isfluid,N);

% calibration except boundries
    i = 1; % four  corners
    j = 1;
    if isfluid(i,j) == 2 
         if isfluid(i+1,j) == 2 && isfluid(i,j+1) == 2 && isfluid(i+1,j+1) == 2 
                u(i,j)= 0;
                v(i,j)= 0;
                rho(i,j)=rho0;
                f(i,j,:)=0;
         end
    end
    i = 2*N+1;
    j = 1;
    if isfluid(i,j) == 2 
         if isfluid(i-1,j) == 2 && isfluid(i,j+1) == 2 && isfluid(i-1,j+1) == 2 
                u(i,j)= 0;
                v(i,j)= 0;
                rho(i,j)=rho0;
                f(i,j,:)=0;
         end
    end
    i = 1;
    j = 2*N+1;
    if isfluid(i,j) == 2 
         if isfluid(i+1,j) == 2 && isfluid(i,j-1) == 2 && isfluid(i+1,j-1) == 2 
                u(i,j)= 0;
                v(i,j)= 0;
                rho(i,j)=rho0;
                f(i,j,:)=0;
         end
    end
    i = 2*N+1;
    j = 2*N+1;
    if isfluid(i,j) == 2 
         if isfluid(i-1,j) == 2 && isfluid(i,j-1) == 2 && isfluid(i-1,j-1) == 2 
                u(i,j)= 0;
                v(i,j)= 0;
                rho(i,j)=rho0;
                f(i,j,:)=0;
         end
    end

    i=1; % four sides
    for j=2:2*N
        if isfluid(i,j) == 2 
            if isfluid(i+1,j) == 2 && isfluid(i,j+1) == 2 && isfluid(i,j-1) == 2 ...
                    && isfluid(i+1,j+1) == 2 && isfluid(i+1,j-1) == 2  
                u(i,j)= 0;
                v(i,j)= 0;
                rho(i,j)=rho0;
                f(i,j,:)=0;

            end
        end
    end
    i=2*N+1;
    for j=2:2*N
        if isfluid(i,j) == 2 
            if isfluid(i-1,j) == 2 && isfluid(i,j+1) == 2 && isfluid(i,j-1) == 2 ...
                    && isfluid(i-1,j+1) == 2 && isfluid(i-1,j-1) == 2
                u(i,j)= 0;
                v(i,j)= 0;
                rho(i,j)=rho0;
                f(i,j,:)=0;
            end
        end
    end

    j=1;
    for i=2:2*N
        if isfluid(i,j) == 2 
            if isfluid(i+1,j) == 2 && isfluid(i-1,j) == 2 && isfluid(i,j+1) == 2 ...
                    && isfluid(i+1,j+1) == 2 && isfluid(i-1,j+1) == 2 
                u(i,j)= 0;
                v(i,j)= 0;
                rho(i,j)=rho0;
                f(i,j,:)=0;
            end
        end
    end
    j=2*N+1;
    for i=2:2*N 
        if isfluid(i,j) == 2 
            if isfluid(i+1,j) == 2 && isfluid(i-1,j) == 2 && isfluid(i,j-1) == 2 ...
                    && isfluid(i+1,j-1) == 2  && isfluid(i-1,j-1) == 2
                u(i,j)= 0;
                v(i,j)= 0;
                rho(i,j)=rho0;
                f(i,j,:)=0;
            end
        end
    end

    for i=2:2*N % internal
        for j=2:2*N
            if isfluid(i,j) == 2 
                if isfluid(i+1,j) == 2 && isfluid(i-1,j) == 2 && isfluid(i,j+1) == 2 && isfluid(i,j-1) == 2 ...
                        && isfluid(i+1,j+1) == 2 && isfluid(i+1,j-1) == 2 && isfluid(i-1,j+1) == 2 && isfluid(i-1,j-1) == 2
                    u(i,j)= 0;
                    v(i,j)= 0;
                    rho(i,j)=rho0;
                    f(i,j,:)=0;
                end
            end
        end
    end

    corr = (norm(u-u_temp)+norm(v-v_temp))/2;
    itr = itr+1;
end

%% plot u
size = 20;
figure('Name','u colormap')
grid on;
hold on;

xlabel('x/H')
ylabel('y/H')
set(gca,'FontSize', size)
set(gcf,'position',[100,100,700,600])

x = 0:ds:1;
y = 0:ds:1;

[C,h] = contourf(x,y,u/u_lb,31);
set(h,'LineColor','none') % removing contour line
%caxis([0 0.5]);   
colormap('jet')
colorbar; % showing color bar
xlim([0 1])
ylim([0 1])
%round corners
%plot(x, L_p/2 - sqrt((L_p/2)^2-(x-L_p/2).^2), 'k');
%a = area(x, L_p/2 - sqrt((L_p/2)^2-(x-L_p/2).^2));
%a.FaceColor = [0 0 0];
%% plot v
size = 20;
figure('Name','v colormap')
grid on;
hold on;

xlabel('x/H')
ylabel('y/H')
set(gca,'FontSize', size)
set(gcf,'position',[100,100,700,600])

x = 0:ds:1;
y = 0:ds:1;

[C,h] = contourf(x,y,v/u_lb, 31);
set(h,'LineColor','none') % removing contour line
%caxis([0 0.5]);   
colormap('jet')
xlim([0 1])
ylim([0 1])

colorbar; % showing color bar
%plot(x, L_p/2 - sqrt((L_p/2)^2-(x-L_p/2).^2), 'k');
%a = area(x, L_p/2 - sqrt((L_p/2)^2-(x-L_p/2).^2));
%a.FaceColor = [0 0 0];

%% plot vorticity
size = 20;
figure('Name','voticity colormap')
grid on;
hold on;

omega = zeros(2*N+1,2*N+1);

xlabel('x/H')
ylabel('y/H')
set(gca,'FontSize', size)
set(gcf,'position',[100,100,700,600])

x = 0:ds:1;
y = 0:ds:1;

for i = 2:2*N
    for j = 2:2*N
        dudy = (u(i+1, j) - u(i-1, j)) / (2 * ds);
        dvdx = (v(i, j+1) - v(i, j-1)) / (2 * ds);
        omega(i,j) = dvdx - dudy;
    end
end

[C,h] = contourf(x,y,omega, 31);
set(h,'LineColor','none') % removing contour line
%caxis([0 0.5]);   
colormap('jet')
colorbar; % showing color bar
xlim([0 1])
ylim([0 1])


%% plot streamfunctions
figure('Name', 'Stream Function'); %to visualize secondary vortices
%grid on;
hold on;

x = 0:ds:1;
y = 0:ds:1;

size = 20;
xlabel('x/H')
ylabel('y/H')
xlim([0 1]);
ylim([0 1]);
set(gca,'FontSize', size)
set(gcf,'position',[100,100,700,700])

psi = zeros(2*N+1,2*N+1);
for i = 2:2*N+1
    psi(i,1) = psi(i-1,1) + u(i,1) * ds;
end

for j = 2:2*N+1
    for i = 1:2*N+1
        psi(i,j) = psi(i,j-1) - v(i,j) * ds;
    end
end
cnt = 0;
for i = 2:2*N % count the positive psi (where the secondary vortices)
    for j = 2:2*N
        if psi(i,j) > 0
            cnt = cnt +1;
        end
    end
end

secondary_vortices_area_percentage = cnt*ds*ds*100 % percentage

surf(x,y,psi)
view([1,1,1])
xlim([0 1])
ylim([0 1])
%zlim([-1e-5 1e-5])
zlim([0 1e-2])
%view([1,0,0])
view([0,0,1])



%% plot streamline
figure('Name', 'Streamline');
%grid on;
hold on;

x = 0:ds:1;
y = 0:ds:1;

size = 20;
xlabel('x/H')
ylabel('y/H')
xlim([0 1]);
ylim([0 1]);
set(gca,'FontSize', size)
set(gcf,'position',[100,100,700,700])

[verts, averts] = streamslice(x,y, u/u_lb, v/u_lb, 15); % 3 -> 15
sl = streamline( [verts averts]);



beep on;
beep;
%%
overlay_geometry(3, ds, L_p);% Overlay chosen geometry shape on all open figures 1–5
% geom_switch = 0 : square
% geom_switch = 1 : square with circle corners
% geom_switch = 2 : triangle
% geom_switch = 3 : pentagon
    
%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%
function f = collide(f, u, v, rho, isfluid, omega, N) % only fluid nodes (1+2) collide while streaming for every nodes
% Setup
w = zeros(9,1);
w(1) = 4/9;
w(2:5) = 1/9;
w(6:9) = 1/36;
c = zeros(9,2);
c(1,:) = [0, 0];
c(2,:) = [1, 0];
c(3,:) = [0, 1];
c(4,:) = [-1, 0];
c(5,:) = [0, -1];
c(6,:) = [1, 1];
c(7,:) = [-1, 1];
c(8,:) = [-1, -1];
c(9,:) = [1, -1];

udotu = u.*u+v.*v;
for k = 1:9
    ckdotu = c(k,1)*u+c(k,2)*v;
    feq = w(k)*rho.*(1+3*ckdotu-1.5*udotu+4.5*ckdotu.^2);
    f(:,:,k) = omega*feq+(1-omega)*f(:,:,k); %dt=1
end

end

function f = stream(f, isfluid, N) % not uniform stream -> streaming along the side walls generates holes while reconstruction -> eventually grows
% stream for all nodes or stream for 1, 2 nodes

% include all nodes
f(:,2:end,2) = f(:,1:end-1,2); % east 
f(2:end,:,3) = f(1:end-1,:,3); % north
f(:,1:end-1,4) = f(:,2:end,4); % west
f(1:end-1,:,5) = f(2:end,:,5); % south
f(2:end,2:end,6) = f(1:end-1,1:end-1,6); % northeast
f(2:end,1:end-1,7) = f(1:end-1,2:end,7); % northwest
f(1:end-1,1:end-1,8) = f(2:end,2:end,8); % southwest
f(1:end-1,2:end,9) = f(2:end,1:end-1,9); % southeast

end

function f = bounceback_full(f, isfluid, N) % full-way staircase % do not bounce for unnecssary direction (i.e. along the wall for boundary nodes)
% check for each boundary node after identifying each position (left, right, bottom)    
    for i=1:2*N+1 % except the moving lid
        for j=1:2*N+1
            if isfluid(i,j) == 2 % boundary node 
                % check the two top corner points first
                if i == 2*N+1 && j == 1
                    if isfluid(i,j+1) == 1 % left
                        f(i,j,:) = bounceback(f(i,j,:),'west');
                    end
                elseif i == 2*N+1 && j == 2*N+1
                    if isfluid(i,j-1) == 1 % right
                        f(i,j,:) = bounceback(f(i,j,:),'east');
                    end
                elseif j == 1
                    if isfluid(i+1,j+1) == 1 && isfluid(i,j+1) == 2 && isfluid(i+1,j) == 2 % bottom left
                        f(i,j,:) = bounceback(f(i,j,:),'south');
                        f(i,j,:) = bounceback(f(i,j,:),'west');
                    else
                        if isfluid(i,j+1) == 1 % left
                            f(i,j,:) = bounceback(f(i,j,:),'west');
                        end
                        if isfluid(i+1,j) == 1 % bottom
                            f(i,j,:) = bounceback(f(i,j,:),'south');
                        end
                    end
                elseif j == 2*N+1
                    if isfluid(i+1,j-1) == 1 && isfluid(i,j-1) == 2 && isfluid(i+1,j) == 2 % bottom right
                        f(i,j,:) = bounceback(f(i,j,:),'south');
                        f(i,j,:) = bounceback(f(i,j,:),'east');
                    else
                        if isfluid(i,j-1) == 1 % right
                            f(i,j,:) = bounceback(f(i,j,:),'east');
                        end
                        if isfluid(i+1,j) == 1 % bottom
                            f(i,j,:) = bounceback(f(i,j,:),'south');
                        end
                    end
                elseif i == 2*N+1 % updated for narrow moving lid case
                    if isfluid(i-1,j+1) == 1 && isfluid(i,j+1) == 2 && isfluid(i-1,j) == 2 % top left
                        f(i,j,:) = bounceback(f(i,j,:),'north');
                        f(i,j,:) = bounceback(f(i,j,:),'west');
                    elseif isfluid(i-1,j-1) == 1 && isfluid(i,j-1) == 2 && isfluid(i-1,j) == 2 % top right
                        f(i,j,:) = bounceback(f(i,j,:),'north');
                        f(i,j,:) = bounceback(f(i,j,:),'east');
                    else
                        if isfluid(i,j+1) == 1 % left
                            f(i,j,:) = bounceback(f(i,j,:),'west');
                        end
                        if isfluid(i,j-1) == 1 % right
                            f(i,j,:) = bounceback(f(i,j,:),'east');
                        end
                        if isfluid(i-1,j) == 1 % top
                            f(i,j,:) = bounceback(f(i,j,:),'north');
                        end
                    end                   
                else
                    if isfluid(i+1,j+1) == 1 && isfluid(i,j+1) == 2 && isfluid(i+1,j) == 2 % bottom left
                        f(i,j,:) = bounceback(f(i,j,:),'south');
                        f(i,j,:) = bounceback(f(i,j,:),'west');
                    elseif isfluid(i+1,j-1) == 1 && isfluid(i,j-1) == 2 && isfluid(i+1,j) == 2 % bottom right
                        f(i,j,:) = bounceback(f(i,j,:),'south');
                        f(i,j,:) = bounceback(f(i,j,:),'east');
                    else
                        if isfluid(i,j+1) == 1 % left
                            f(i,j,:) = bounceback(f(i,j,:),'west');
                        end
                        if isfluid(i,j-1) == 1 % right
                            f(i,j,:) = bounceback(f(i,j,:),'east');
                        end
                        if isfluid(i+1,j) == 1 % bottom
                            f(i,j,:) = bounceback(f(i,j,:),'south');
                        end
                    end
                end
            end
        end
    end
end

function f = bounceback(f, side)
%fprintf('\n 1 \n')
if strcmp(side, 'north') % the wall is north wall
    %fprintf('\n n \n')
    f(5) = f(3);
    f(9) = f(7);
    f(8) = f(6);
end
if strcmp(side, 'south')
    f(3) = f(5);
    f(6) = f(8);
    f(7) = f(9);
end
if strcmp(side, 'east')
    %fprintf('\n e \n')
    f(4) = f(2);
    f(8) = f(6);
    f(7) = f(9);
end
if strcmp(side, 'west')
    %fprintf('\n w \n')
    f(2) = f(4);
    f(6) = f(8);
    f(9) = f(7);
end
end

function [u,v,rho] = reconstruct(f,u,v,isfluid,N) % IMPORTANT: resconstruct should be done only fluid node apart from boundary
c = zeros(9,2);
c(1,:) = [0, 0];
c(2,:) = [1, 0];
c(3,:) = [0, 1];
c(4,:) = [-1, 0];
c(5,:) = [0, -1];
c(6,:) = [1, 1];
c(7,:) = [-1, 1];
c(8,:) = [-1, -1];
c(9,:) = [1, -1];

rho = sum(f,3);
u(:,:) = 0; % 0 for internal nodes
v(:,:) = 0;
for k = 1:9
    u(:,:) = u(:,:) + c(k,1)*f(:,:,k);
    v(:,:) = v(:,:) + c(k,2)*f(:,:,k);
end
u(:,:) = u(:,:)./rho(:,:);
v(:,:) = v(:,:)./rho(:,:);

end

function overlay_geometry(geom_switch, ds, L_p)
    if geom_switch == 0
        return; % nothing to do
    end

    % loop through figures 1–5
    for fig_id = 1:5
        if isgraphics(fig_id, 'figure')
            figure(fig_id); hold on;

            switch geom_switch
                case 1  % Circle
                    x = 0:ds:1;
                    y = L_p/2 - sqrt((L_p/2)^2 - (x - L_p/2).^2);
                    plot(x, y, 'k');
                    a = area(x, y);
                    a.FaceColor = [0 0 0];

                case 2  % Triangle
                    x1 = 0:ds:0.5;
                    x2 = 0.5:ds:1;
                    a1 = area(x1, tand(60*2)*x1 + L_p);
                    a2 = area(x2, tand(60)*(x2 - L_p) + L_p);
                    a1.FaceColor = [0 0 0];
                    a2.FaceColor = [0 0 0];

                case 3  % Pentagon
                    x = 0:ds:1;
                    a1 = area(x, -tand(36)*(x - L_p/2));
                    a2 = area(x,  tand(36)*(x - L_p/2));
                    y1 = tand(72)*x + L_p/2*tand(36);
                    patch([x fliplr(x)], [y1 max(ylim)*ones(1,length(x))], 'k');
                    y2 = tand(108)*(x - L_p) + L_p/2*tand(36);
                    patch([x fliplr(x)], [y2 max(ylim)*ones(1,length(x))], 'k');
                    a1.FaceColor = [0 0 0];
                    a2.FaceColor = [0 0 0];
            end
        end
    end
end
        
    

