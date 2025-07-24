% --------------------------------------------------------------------
%initialization

clear
clc
clf

% --------------------------------------------------------------------

Timestep = 12000

% --------------------------------------------------------------------
%load parameters

grid = readmatrix("../grid.csv");
Lxtot = grid(1,1);
Lytot = grid(1,2);
Lztot = grid(1,3);
Nxtot = grid(2,1);
Nytot = grid(2,2);
Nztot = grid(2,3);

% --------------------------------------------------------------------
%make tables
x = zeros(Nxtot,Nytot,Nztot);
y = zeros(Nxtot,Nytot,Nztot);
z = zeros(Nxtot,Nytot,Nztot);
u = zeros(Nxtot,Nytot,Nztot);
v = zeros(Nxtot,Nytot,Nztot);
w = zeros(Nxtot,Nytot,Nztot);
p = zeros(Nxtot,Nytot,Nztot);
q = zeros(Nxtot,Nytot,Nztot);

% --------------------------------------------------------------------
%write video frames
fname = append('../Anim3D/uvwp.',num2str(Timestep,'%0.6i'));
fid = fopen(fname); % Open file
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
            q(i,j,k) = tmp(8); 
        end
    end
end
fclose(fid);

figure(1)%snapshot of pressure field and streamline(2D)
xp = linspace(0,Lxtot,Nxtot);
yp = linspace(0,Lytot,Nytot);
zp = linspace(0,Lztot,Nztot);
xp2 = linspace(0,Lxtot,Nxtot/4);
yp2 = linspace(0,Lytot,Nytot/4);
zp2 = linspace(0,Lztot,Nztot/4);
[X,Z] = meshgrid(xp,zp);
[X2,Z2] = meshgrid(xp2,zp2);
pplot = squeeze(p(:,Nytot/2,:));
for i=1:Nxtot
    for k = 1:Nztot
        if pplot(i,k)>1
            pplot(i,k)=1;
        elseif pplot(i,k)<-0.5
            pplot(i,k)=-0.5;
        end
    end
end
uplot = squeeze(u(:,Nytot/2,:));
vplot = squeeze(v(:,Nytot/2,:));
wplot = squeeze(w(:,Nytot/2,:));
pplot2 = zeros(Nxtot/4, Nztot/4);
uplot2 = zeros(Nxtot/4, Nztot/4);
vplot2 = zeros(Nxtot/4, Nztot/4);
wplot2 = zeros(Nxtot/4, Nztot/4);
for i=1:Nxtot/4
    for k=1:Nztot/4
        pplot2(i,k) = pplot(i*4-2,k*4-2);
        uplot2(i,k) = uplot(i*4-2,k*4-2);
        vplot2(i,k) = vplot(i*4-2,k*4-2);
        wplot2(i,k) = wplot(i*4-2,k*4-2);
    end
end
%angl2rad = pi/180;
%dataN = 100;
%theta=linspace(0,360,dataN)*angl2rad;
%xc = 2-(Lxtot/Nxtot)*3/4; 
%yc = 2-(Lytot/Nytot)/2;
%xx = 0.5*cos(theta) + xc;
%yy = 0.5*sin(theta) + yc;
%[startX,startY] = meshgrid(0, 1.6:0.05:2.6);
hold on
contourf(X,Z,pplot',10000,'LineStyle','none')
%quiver(X2',Z2',uplot2,wplot2,"black")
%verts = stream2(X,Z,uplot',wplot',startX,startY);
%lineobj = streamline(verts);
%for i=1:21
%    lineobj(i).Color="black";
%end
%fill(xx,yy, 'w');
c = colorbar;
c.Label.String = 'p';
c.Label.FontAngle='italic';
colormap jet
clim([-0.5 1])
axis([0 8 0 4])
axis equal
set(gca, 'FontName','Times','FontSize',50)
xlabel('$$x$$','interpreter','latex')
ylabel('$$z$$','interpreter','latex')

figure(2)%snapshot of vortex structure(3D)
colors = zeros(Nxtot,Nytot,Nztot);%For gradation
colors = u;
fv = patch(isosurface(x,y,z,q,0.03,colors));
set(fv,'FaceAlpha',0.5)
set(fv,'FaceColor','interp')
set(fv,'Edgecolor','none')
set(gca, 'FontName','Times','FontSize',40)
xlabel('$$x$$','interpreter','latex')
ylabel('$$y$$','interpreter','latex')
zlabel('$$z$$','interpreter','latex')
camlight('headlight')
lighting gouraud
view([20 60])
clim([-1 2])
colormap jet
mymap = colormap(jet);
colormap(mymap)
c = colorbar;
c.Label.String = 'u';
c.Label.FontAngle='italic';
daspect([1 1 1])
axis([0.1 Lxtot-0.1 0 Lytot 0 Lztot])
hold on
[X3, Y3, Z3] = meshgrid(xp,yp,zp);
[startX3,startY3, startZ3] = meshgrid(0, 2, 1.5:0.05:2.5);
u_perm = permute(u, [2 1 3]);
v_perm = permute(v, [2 1 3]);
w_perm = permute(w, [2 1 3]);
verts3 = stream3(X3, Y3, Z3, u_perm, v_perm, w_perm, startX3, startY3, startZ3);
lineobj = streamline(verts3);
for i=1:21
    lineobj(i).Color="k";
    lineobj(i).LineWidth=0.5;
end
