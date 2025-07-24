%
% CDCL_time.m
% Last update: 2025-02-05 by I. Takarabe
%

% --------------------------------------------------------------------
%initialization

clear 
clc
clf

% --------------------------------------------------------------------
%set parameters

Initial = 0;
Nt = 3;
history = zeros(3,Nt);

% --------------------------------------------------------------------
%load parameters

inp = readtable("../inp.csv");
grid = readmatrix("../grid.csv");
dt = inp{4,1};
Re = inp{5,1};
Lxtot = grid(1,1);
Lytot = grid(1,2);
Lztot = grid(1,3);
Nxtot = grid(2,1);
Nytot = grid(2,2);
Nztot = grid(2,3);
Lx = grid(3,1);
Ly = grid(3,2);
Lz = grid(3,3);
Nx = grid(4,1);
Ny = grid(4,2);
Nz = grid(4,3);
dx = grid(5,1);
dy = grid(5,2);
dz = grid(5,3);
xshift = grid(6,1);
yshift = grid(6,2);
zshift = grid(6,3);

% --------------------------------------------------------------------
%tables

x = zeros(Nxtot,Nytot,Nztot);
y = zeros(Nxtot,Nytot,Nztot);
z = zeros(Nxtot,Nytot,Nztot);
u = zeros(Nxtot,Nytot,Nztot);
v = zeros(Nxtot,Nytot,Nztot);
w = zeros(Nxtot,Nytot,Nztot);
p = zeros(Nxtot,Nytot,Nztot);
u_old = zeros(Nxtot,Nytot,Nztot);
v_old = zeros(Nxtot,Nytot,Nztot);
w_old = zeros(Nxtot,Nytot,Nztot);
p_old = zeros(Nxtot,Nytot,Nztot);
u_new = zeros(Nxtot,Nytot,Nztot);
v_new = zeros(Nxtot,Nytot,Nztot);
w_new = zeros(Nxtot,Nytot,Nztot);
p_new = zeros(Nxtot,Nytot,Nztot);

u_zoom = zeros(Nx,Ny,Nz);
v_zoom = zeros(Nx,Ny,Nz);
w_zoom = zeros(Nx,Ny,Nz);
p_zoom = zeros(Nx,Ny,Nz);
u_zoom_old = zeros(Nx,Ny,Nz);
v_zoom_old = zeros(Nx,Ny,Nz);
w_zoom_old = zeros(Nx,Ny,Nz);
p_zoom_old = zeros(Nx,Ny,Nz);
u_zoom_new = zeros(Nx,Ny,Nz);
v_zoom_new = zeros(Nx,Ny,Nz);
w_zoom_new = zeros(Nx,Ny,Nz);
p_zoom_new = zeros(Nx,Ny,Nz);

moment = zeros(Nx,Ny,Nz,3);
pressure = zeros(Nx,Ny,Nz,3);
viscousity = zeros(Nx,Ny,Nz,3);
time = zeros(Nx,Ny,Nz,3);

dudx = zeros(Nx,Ny,Nz);
dudy = zeros(Nx,Ny,Nz);
dudz = zeros(Nx,Ny,Nz);

dvdx = zeros(Nx,Ny,Nz);
dvdy = zeros(Nx,Ny,Nz);
dvdz = zeros(Nx,Ny,Nz);

dwdx = zeros(Nx,Ny,Nz);
dwdy = zeros(Nx,Ny,Nz);
dwdz = zeros(Nx,Ny,Nz);

dpdx = zeros(Nx,Ny,Nz);
dpdy = zeros(Nx,Ny,Nz);
dpdz = zeros(Nx,Ny,Nz);

dudt = zeros(Nx,Ny,Nz);
dvdt = zeros(Nx,Ny,Nz);
dwdt = zeros(Nx,Ny,Nz);

dudx2 = zeros(Nx,Ny,Nz);
dudy2 = zeros(Nx,Ny,Nz);
dudz2 = zeros(Nx,Ny,Nz);

dvdx2 = zeros(Nx,Ny,Nz);
dvdy2 = zeros(Nx,Ny,Nz);
dvdz2 = zeros(Nx,Ny,Nz);

dwdx2 = zeros(Nx,Ny,Nz);
dwdy2 = zeros(Nx,Ny,Nz);
dwdz2 = zeros(Nx,Ny,Nz);

% --------------------------------------------------------------------
%load file

for it = 1:Nt
Cd = 0;
Cl_y = 0;
Cl_z = 0;
fname1 = append('../Anim3D/uvwp.',num2str(Initial+30*it,'%0.6i'));
fid = fopen(fname1); % Open file
for k=1:Nztot
    for j=1:Nytot
        for i=1:Nxtot
            [tmp, count] = fscanf(fid,'%e',8); % Read 8 sequential data
            x(i,j,k) = tmp(1);
            y(i,j,k) = tmp(2);
            z(i,j,k) = tmp(3);
            u(i,j,k) = tmp(4);
            v(i,j,k) = tmp(5);
            w(i,j,k) = tmp(6);
            p(i,j,k) = tmp(7); 
        end
    end
end
fclose(fid);
u_zoom = u(xshift+1:xshift+Nx,yshift+1:yshift+Ny,zshift+1:zshift+Nz);
v_zoom = v(xshift+1:xshift+Nx,yshift+1:yshift+Ny,zshift+1:zshift+Nz);
w_zoom = w(xshift+1:xshift+Nx,yshift+1:yshift+Ny,zshift+1:zshift+Nz);
p_zoom = p(xshift+1:xshift+Nx,yshift+1:yshift+Ny,zshift+1:zshift+Nz);

fname2 = append('Anim3D/uvwp.',num2str(Initial+30*it-1,'%0.6i'));
fid = fopen(fname2); % Open file
for k=1:Nztot
    for j=1:Nytot
        for i=1:Nxtot
            [tmp, count] = fscanf(fid,'%e',8); % Read 8 sequential data
            u_old(i,j,k) = tmp(4);
            v_old(i,j,k) = tmp(5);
            w_old(i,j,k) = tmp(6);
            p_old(i,j,k) = tmp(7); 
        end
    end
end
fclose(fid);
u_zoom_old = u_old(xshift+1:xshift+Nx,yshift+1:yshift+Ny,zshift+1:zshift+Nz);
v_zoom_old = v_old(xshift+1:xshift+Nx,yshift+1:yshift+Ny,zshift+1:zshift+Nz);
w_zoom_old = w_old(xshift+1:xshift+Nx,yshift+1:yshift+Ny,zshift+1:zshift+Nz);
p_zoom_old = p_old(xshift+1:xshift+Nx,yshift+1:yshift+Ny,zshift+1:zshift+Nz);

fname3 = append('Anim3D/uvwp.',num2str(Initial+30*it+1,'%0.6i'));
fid = fopen(fname3); % Open file
for k=1:Nztot
    for j=1:Nytot
        for i=1:Nxtot
            [tmp, count] = fscanf(fid,'%e',8); % Read 8 sequential data
            u_new(i,j,k) = tmp(4);
            v_new(i,j,k) = tmp(5);
            w_new(i,j,k) = tmp(6);
            p_new(i,j,k) = tmp(7); 
        end
    end
end
fclose(fid);
u_zoom_new = u_new(xshift+1:xshift+Nx,yshift+1:yshift+Ny,zshift+1:zshift+Nz);
v_zoom_new = v_new(xshift+1:xshift+Nx,yshift+1:yshift+Ny,zshift+1:zshift+Nz);
w_zoom_new = w_new(xshift+1:xshift+Nx,yshift+1:yshift+Ny,zshift+1:zshift+Nz);
p_zoom_new = p_new(xshift+1:xshift+Nx,yshift+1:yshift+Ny,zshift+1:zshift+Nz);
% --------------------------------------------------------------------
%fill the tables

for k = 2:Nz-1
    for j = 2:Ny-1
        for i = 2:Nx-1
            dudx(i,j,k) = (u_zoom(i+1,j,k)-u_zoom(i-1,j,k))/(2*dx);
            dudy(i,j,k) = (u_zoom(i,j+1,k)-u_zoom(i,j-1,k))/(2*dy);
            dudz(i,j,k) = (u_zoom(i,j,k+1)-u_zoom(i,j,k-1))/(2*dz);
    
            dvdx(i,j,k) = (v_zoom(i+1,j,k)-v_zoom(i-1,j,k))/(2*dx);
            dvdy(i,j,k) = (v_zoom(i,j+1,k)-v_zoom(i,j-1,k))/(2*dy);
            dvdz(i,j,k) = (v_zoom(i,j,k+1)-v_zoom(i,j,k-1))/(2*dz);
    
            dwdx(i,j,k) = (w_zoom(i+1,j,k)-w_zoom(i-1,j,k))/(2*dx);
            dwdy(i,j,k) = (w_zoom(i,j+1,k)-w_zoom(i,j-1,k))/(2*dy);
            dwdz(i,j,k) = (w_zoom(i,j,k+1)-w_zoom(i,j,k-1))/(2*dz);
    
            dpdx(i,j,k) = (p_zoom(i+1,j,k)-p_zoom(i-1,j,k))/(2*dx);
            dpdy(i,j,k) = (p_zoom(i,j+1,k)-p_zoom(i,j-1,k))/(2*dy);
            dpdz(i,j,k) = (p_zoom(i,j,k+1)-p_zoom(i,j,k-1))/(2*dz);

            dudt(i,j,k) = (u_zoom_new(i,j,k)-u_zoom_old(i,j,k))/(2*dt);
            dvdt(i,j,k) = (v_zoom_new(i,j,k)-v_zoom_old(i,j,k))/(2*dt);
            dwdt(i,j,k) = (w_zoom_new(i,j,k)-w_zoom_old(i,j,k))/(2*dt);

            dudx2(i,j,k) = (u_zoom(i-1,j,k)-2*u_zoom(i,j,k)+u_zoom(i+1,j,k))/(dx^2);
            dudy2(i,j,k) = (u_zoom(i,j-1,k)-2*u_zoom(i,j,k)+u_zoom(i,j+1,k))/(dy^2);
            dudz2(i,j,k) = (u_zoom(i,j,k-1)-2*u_zoom(i,j,k)+u_zoom(i,j,k+1))/(dz^2);
            
            dvdx2(i,j,k) = (v_zoom(i-1,j,k)-2*v_zoom(i,j,k)+v_zoom(i+1,j,k))/(dx^2);
            dvdy2(i,j,k) = (v_zoom(i,j-1,k)-2*v_zoom(i,j,k)+v_zoom(i,j+1,k))/(dy^2);
            dvdz2(i,j,k) = (v_zoom(i,j,k-1)-2*v_zoom(i,j,k)+v_zoom(i,j,k+1))/(dz^2);
    
            dwdx2(i,j,k) = (w_zoom(i-1,j,k)-2*w_zoom(i,j,k)+w_zoom(i+1,j,k))/(dx^2);
            dwdy2(i,j,k) = (w_zoom(i,j-1,k)-2*w_zoom(i,j,k)+w_zoom(i,j+1,k))/(dy^2);
            dwdz2(i,j,k) = (w_zoom(i,j,k-1)-2*w_zoom(i,j,k)+w_zoom(i,j,k+1))/(dz^2);
        end
    end
end

% --------------------------------------------------------------------
%compute the drag force

for k = 2:Nz-1
    for j = 2:Ny-1
        for i = 2:Nx-1
            moment(i,j,k,1) = u_zoom(i,j,k)*dudx(i,j,k)+v_zoom(i,j,k)*dudy(i,j,k)+w_zoom(i,j,k)*dudz(i,j,k);
            moment(i,j,k,2) = u_zoom(i,j,k)*dvdx(i,j,k)+v_zoom(i,j,k)*dvdy(i,j,k)+w_zoom(i,j,k)*dvdz(i,j,k);
            moment(i,j,k,3) = u_zoom(i,j,k)*dwdx(i,j,k)+v_zoom(i,j,k)*dwdy(i,j,k)+w_zoom(i,j,k)*dwdz(i,j,k);
        end
    end
end

for k = 2:Nz-1
    for j = 2:Ny-1
        for i = 2:Nx-1
            pressure(i,j,k,1) = dpdx(i,j,k);
            pressure(i,j,k,2) = dpdy(i,j,k);
            pressure(i,j,k,3) = dpdz(i,j,k);
        end
    end
end

for k = 2:Nz-1
    for j = 2:Ny-1
        for i = 2:Nx-1
            viscousity(i,j,k,1) = (dudx2(i,j,k)+dudy2(i,j,k)+dudz2(i,j,k))/Re;
            viscousity(i,j,k,2) = (dvdx2(i,j,k)+dvdy2(i,j,k)+dvdz2(i,j,k))/Re;
            viscousity(i,j,k,3) = (dwdx2(i,j,k)+dwdy2(i,j,k)+dwdz2(i,j,k))/Re;
        end
    end
end

for k = 2:Nz-1
    for j = 2:Ny-1
        for i = 2:Nx-1
            time(i,j,k,1) = dudt(i,j,k);
            time(i,j,k,2) = dvdt(i,j,k);
            time(i,j,k,3) = dwdt(i,j,k);
        end
    end
end

dragforce = time + moment + pressure - viscousity;

figure(1)
xp = linspace(0,Lx,Nx);
yp = linspace(0,Ly,Ny);
zp = linspace(0,Lz,Nz);
[X,Z] = meshgrid(xp,zp);
drag_plot = squeeze(dragforce(:,Ny/2,:,1));
contourf(X,Z,drag_plot',100,'LineStyle','none')
axis equal
colorbar

for k = 2:Nz-1
    for j = 2:Ny-1
        for i = 2:Nx-1
            Cd = Cd + dragforce(i,j,k,1)*dx*dy*dz;
            Cl_y = Cl_y + dragforce(i,j,k,2)*dx*dy*dz;
            Cl_z = Cl_z + dragforce(i,j,k,3)*dx*dy*dz;
        end
    end
end

Cd = -Cd * 8/pi;
Cl_y = -Cl_y * 8/pi;
Cl_z = -Cl_z * 8/pi;
history(1,it) = Cd;
history(2,it) = Cl_y;
history(3,it) = Cl_z;
end
figure(2)
for it=1:Nt
    plot(it*dt, history(1,it))
    plot(it*dt, history(2,it))
    plot(it*dt, history(3,it))
end

