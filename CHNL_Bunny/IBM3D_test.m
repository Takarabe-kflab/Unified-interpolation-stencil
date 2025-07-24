%test for GP
%Requires long time to draw. Please wait for 30~40s.
clf

figure(1)
for i = 1:Nx
    for j = 1:Ny
        for k = 1:Nz
            if GP(i,j,k) == 1
                scatter3(i, j, k, 20, "filled", "r")
                hold on
            end
        end
    end
end
%test for PP, IP and SP
%for n = 1:NGP
%    scatter3(IP(n,1)/dx, IP(n,2)/dy, IP(n,3)/dz, 20, "filled","b")
%    scatter3(PP(n,1)/dx, PP(n,2)/dy, PP(n,3)/dz, 20, "filled","g")
%    if RP(n,1,1) ~= 0 
%        scatter3(RP(n,1,1)/dx, RP(n,1,2)/dy, RP(n,1,3)/dz, 20, "^", "filled", "k")
%    end
%    if RP(n,3,1) ~= 0
%        scatter3(RP(n,3,1)/dx, RP(n,3,2)/dy, RP(n,3,3)/dz, 20, "^", "filled", "k")
%    end
%end

%--------------------------------------------------------------------------
%test for normal vector
%for n=1:NGP
%    IPdist = double(sqrt((IP(n,1)-PP(n,1))^2 + (IP(n,2)-PP(n,2))^2 + (IP(n,3)-PP(n,3))^2));
%    quiver3(GPnormal(n,1)/dx,GPnormal(n,2)/dy,GPnormal(n,3)/dz,GPnormal(n,4)*(RP(n,3,4)+IPdist)/dx,GPnormal(n,5)*(RP(n,3,4)+IPdist)/dy,GPnormal(n,6)*(RP(n,3,4)+IPdist)/dz,"k"); % draws arrows
%    hold on
%end

%--------------------------------------------------------------------------
%isosurface
xp = 1:1:Nx;
yp = 1:1:Ny;
zp = 1:1:Nz;
[X,Y,Z] = meshgrid(yp,xp,zp);
[faces,verts,colors]=isosurface(Y,X,Z,phi,1e-4,Z);
patch('Vertices',verts,'Faces',faces,'FaceVertexCData',colors,'FaceColor','interp','EdgeColor','none','FaceAlpha',0.4)
view(3)
colormap jet
hold on

%--------------------------------------------------------------------------
axis equal
grid on
set(gca, 'FontName','Times','FontSize',25)
xlabel('$$x$$','interpreter','latex')
ylabel('$$y$$','interpreter','latex')
zlabel('$$z$$','interpreter','latex')
xlim([0 Nx])
ylim([0 Ny])
zlim([0 Nz])



