%
% Original: 2024-05-17 by K. Fukagata (Qisosurf.m)
% Last update: 2024-11-05 by I. Takarabe
%

% --------------------------------------------------------------------
%initialization

clear
clc
clf

% --------------------------------------------------------------------
%set video parameters
%"Initial" is the timestep where the video starts
%"Nt" is the number of frames

Initial = 0
Nt = 40

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

mov = VideoWriter('u.avi');
mov.Quality = 100;
mov.FrameRate = 10;
open(mov);

frames(Nt) = struct('cdata', [], 'colormap', []); 

for it = 1:Nt
fname = append('../Anim3D/uvwp.',num2str(Initial+30*it,'%0.6i'))

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

xp = linspace(0,Lxtot,Nxtot);
yp = linspace(0,Lytot,Nytot);
zp = linspace(0,Lztot,Nztot);
[X,Z] = meshgrid(xp,zp);
plot = squeeze(abs(u(:,Nytot/2,:)));
fig=figure()
contourf(X,Z,plot',1000,"LineStyle","none")
set(gca, 'FontName','Times','FontSize',15)
xlabel('$$x$$','interpreter','latex')
ylabel('$$y$$','interpreter','latex')
clim([-0.5 1.5])
colormap jet
mymap = colormap(jet);
colormap(mymap)
colorbar
axis([0 Lxtot 0 Lztot])
axis equal

frames(it) = getframe(fig); 

close(fig)
end

writeVideo(mov, frames);
close(mov);


