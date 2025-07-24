%
% IBM3D.m
% Last update: 2025-02-21 by I. Takarabe
%

% --------------------------------------------------------------------
%initialization
clear
clc

%grid parameters
%DNS parameters
Lxtot = 8; %total channel length
Lytot = 4; %total channel height
Lztot = 4; %total channel width
Nxtot = 256;
Nytot = 128;
Nztot = 128;

% This code only focuses on a "zoomed part" of the whole domain
Nx = 64;    
Ny = 64;  
Nz = 64;
dx = Lxtot/Nxtot;  % cell length
dy = Lytot/Nytot;  % cell height
dz = Lztot/Nztot;  % cell width
Lx = Nx*dx;    % narrow channel length
Ly = Ny*dy;    % narrow channel height
Lz = Nz*dz;    % narrow channel width

% set the location of zoomed part
xshift = (Nxtot/4) - (Nx/2);
yshift = (Nytot/2) - (Ny/2);
zshift = (Nztot/2) - (Nz/2);

R = 0.5; % radius of the sphere

% interpolation parameters
dxi = 0.001*dx; %step of x 
% (for the interpolation, how close to the grid line is considered "on the line")

%center of sphere
xc0 = 1;
yc0 = 1;
zc0 = 1;

%make a table (for output)
grid = [Lxtot Lytot Lztot;Nxtot Nytot Nztot;Lx Ly Lz;Nx Ny Nz;dx dy dz;xshift yshift zshift];

%count the number of penalized GP
penalization_count = [0,0,0];

% --------------------------------------------------------------------
%main loop

for uvw = 1:3   %for making 3 kernels for u,v,w to fit the staggered grid

    if uvw == 1
        xc = xc0 - dx/2; %shift xc
        yc = yc0;
        zc = zc0;
        fname = "GP_weights_u.txt";
    elseif uvw == 2
        xc = xc0;
        yc = yc - dy/2; %shift yc
        zc = zc0;
        fname = "GP_weights_v.txt";
    elseif uvw == 3
        xc = xc0;
        yc = yc0; 
        zc = zc - dz/2; %shift zc
        fname = "GP_weights_w.txt";
    end

    % ----------------------------------------------------------------
    %level set function
    %phi is the distance of each points from the surface
    
    phi = zeros(Nx,Ny,Nz);
    for i = 1:Nx
        for j = 1:Ny
            for k = 1:Nz
                phi(i,j,k) = sqrt((double(i)*dx-xc)^2 + (double(j)*dy-yc)^2 + (double(k)*dz-zc)^2) - R;
            end
        end
    end
    
    % ----------------------------------------------------------------
    %"inside" indicates whether each points are inside the boundary or not
    %inside(i,j,k)=1 if (i,j,k) is inside the boundary. If not, inside(i,j,k)=0
    
    inside = zeros(Nx,Ny,Nz);
    for i = 1:Nx
        for j = 1:Ny
            for k = 1:Nz
                if phi(i,j,k) <= 0
                    inside(i,j,k) = 1;
                end
            end
        end
    end

    % ----------------------------------------------------------------
    %"GP" indicates whether each points are ghost points
    %GP(i,j,k)=1 if (i,j,k) is a ghost point. If not, GP(i,j,k)=0

    GP = zeros(Nx,Ny,Nz);
    NGP = 0;%number of GP
    
    for i = 1:Nx
        for j = 1:Ny
            for k = 1:Nz
                if inside(i,j,k) == 1
                    if inside(i-1,j,k) == 0
                        GP(i,j,k) = 1;
                        NGP = NGP + 1;
                    elseif inside(i+1,j,k) == 0
                        GP(i,j,k) = 1;
                        NGP = NGP + 1;
                    elseif inside(i,j-1,k) == 0
                        GP(i,j,k) = 1;
                        NGP = NGP + 1; 
                    elseif inside(i,j+1,k) == 0
                        GP(i,j,k) = 1;
                        NGP = NGP + 1;
                    elseif inside(i,j,k-1) == 0
                        GP(i,j,k) = 1;
                        NGP = NGP + 1;
                    elseif inside(i,j,k+1) == 0
                        GP(i,j,k) = 1;
                        NGP = NGP + 1;
                    end
                end 
            end
        end
    end
    
    % ----------------------------------------------------------------
    %set the output matrix (3 + 5x5x5 + 1 columns and as many lines as there are GP)
    %the 3 first columns are for the coordinates of the GP, the 5x5x5 are the
    %convolution weights of the node points arount the image point, and the
    %final one is the value of the flow property on the object at the physical
    %point (PP)
    
    GPconv = zeros(NGP, 129);
    n = 1;
    while n <= NGP
        for i = 1:Nx
            for j = 1:Ny
                for k = 1:Nz
                    if GP(i,j,k) == 1
                        GPconv(n,1) = i*dx;
                        GPconv(n,2) = j*dy;
                        GPconv(n,3) = k*dz;
                        n = n + 1;
                    end
                end
            end
        end
    end
    
    % ----------------------------------------------------------------
    %set a normal vector for each GP
    
    GPnormal = zeros(NGP,6);
    normalvector = zeros(1,3);
    for n=1:NGP
        GPnormal(n,1) = GPconv(n,1); % x of the GP
        GPnormal(n,2) = GPconv(n,2); % y of the GP
        GPnormal(n,3) = GPconv(n,3); % z of the GP
        i = int32(GPconv(n,1)/dx);
        j = int32(GPconv(n,2)/dy);
        k = int32(GPconv(n,3)/dy);
        %calculates the gradient of phi
        gradphi = [(phi(i+1,j,k) - phi(i-1,j,k))/(2*dx) , (phi(i,j+1,k) - phi(i,j-1,k))/(2*dy) , (phi(i,j,k+1) - phi(i,j,k-1))/(2*dz)];
        normgradphi = sqrt(((phi(i+1,j,k) - phi(i-1,j,k))/(2*dx))^2 + ((phi(i,j+1,k) - phi(i,j-1,k))/(2*dy))^2 + ((phi(i,j,k+1) - phi(i,j,k-1))/(2*dz))^2);
        normalvector = gradphi/normgradphi;
        GPnormal(n,4) = normalvector(1);
        GPnormal(n,5) = normalvector(2);
        GPnormal(n,6) = normalvector(3);
    end
    
    %-----------------------------------------------------------------
    %"PP" indicates the location of each physical points
    %"IP" indicates the location of each image points
    
    IP = zeros(NGP,3);
    PP = zeros(NGP,3);
    for n=1:NGP
        i = int32(GPconv(n,1)/dx);
        j = int32(GPconv(n,2)/dy);
        k = int32(GPconv(n,3)/dy);
        PP(n,1) = double(i)*dx - phi(i,j,k)*GPnormal(n,4);
        PP(n,2) = double(j)*dy - phi(i,j,k)*GPnormal(n,5);
        PP(n,3) = double(k)*dz - phi(i,j,k)*GPnormal(n,6);
        IP(n,1) = double(i)*dx - 2*phi(i,j,k)*GPnormal(n,4);
        IP(n,2) = double(j)*dy - 2*phi(i,j,k)*GPnormal(n,5);
        IP(n,3) = double(k)*dz - 2*phi(i,j,k)*GPnormal(n,6);
    end
    
    %-----------------------------------------------------------------
    %"RPcand" indicates the location of candidates for reference points
    %"walldist" indicates the distance from IP to the adjacent and sub-adjacent yz-plane, xz-plane and xy-plane
    %"xgrid","ygrid","zgrid" are the floor of the x,y,z coordinates of candidates
    RPcand = zeros(NGP,14,5);
    walldist = zeros(NGP,12);
    IPx = 0;
    IPy = 0;
    IPz = 0;
    xgrid = 0;
    ygrid = 0;
    zgrid = 0;
    %calculation of RPcand
    for n=1:NGP
        IPx = double(IP(n,1)/dx);
        IPy = double(IP(n,2)/dz);
        IPz = double(IP(n,3)/dy);
        IPdist = double(sqrt((IP(n,1)-PP(n,1))^2 + (IP(n,2)-PP(n,2))^2 + (IP(n,3)-PP(n,3))^2)); %distance between PP and IP
        %add SP to the list
        for i = 1:4
            walldist(n,i) = double(floor(IPx)-IPx+i-2)*dx;    %distance to y-z planes
            walldist(n,i+4) = double(floor(IPy)-IPy+i-2)*dy;  %distance to x-z planes
            walldist(n,i+8) = double(floor(IPz) -IPz+i-2)*dz; %distance to x-y planes
            %location of intersection points with y-z planes
            xgrid = int8((IP(n,1) + walldist(n,i))/dx);
            ygrid = floor((IP(n,2) + GPnormal(n,5)/GPnormal(n,4)*walldist(n,i))/dy);
            zgrid = floor((IP(n,3) + GPnormal(n,6)/GPnormal(n,4)*walldist(n,i))/dz);
            RPcandx = IP(n,1) + walldist(n,i);
            RPcandy = IP(n,2) + GPnormal(n,5)/GPnormal(n,4)*walldist(n,i);
            RPcandz = IP(n,3) + GPnormal(n,6)/GPnormal(n,4)*walldist(n,i);
            if 1<=xgrid && xgrid<=Nx-1 && 1<=ygrid && ygrid<=Ny-1 && 1<=zgrid && zgrid<=Nz-1
                if inside(xgrid,ygrid,zgrid)==0 && inside(xgrid,ygrid+1,zgrid)==0 && inside(xgrid,ygrid,zgrid+1)==0 && inside(xgrid,ygrid+1,zgrid+1)==0
                    if abs(GPconv(n,1) - RPcandx) <= 2*dx && abs(GPconv(n,2) - RPcandy) <= 2*dy && abs(GPconv(n,3) - RPcandz) <= 2*dz
                        RPcand(n,i,1) = IP(n,1) + walldist(n,i);                             %x coordinate of intersection point
                        RPcand(n,i,2) = IP(n,2) + GPnormal(n,5)/GPnormal(n,4)*walldist(n,i); %y coordinate of intersection point
                        RPcand(n,i,3) = IP(n,3) + GPnormal(n,6)/GPnormal(n,4)*walldist(n,i); %z coordinate of intersection point
                        RPcand(n,i,4) = IPdist + walldist(n,i)/GPnormal(n,4);                %xi(Distance along the normal vector from PP)
                        RPcand(n,i,5) = 1;                                                   %a sign to tell that this SP is on the y-z plane
                    end
                end
            end
            %location of intersection points with x-z planes
            xgrid = floor((IP(n,1) + GPnormal(n,4)/GPnormal(n,5)*walldist(n,i+4))/dx);
            ygrid = int8((IP(n,2) + walldist(n,i+4))/dy);
            zgrid = floor((IP(n,3) + GPnormal(n,6)/GPnormal(n,5)*walldist(n,i+4))/dz);
            RPcandx = IP(n,1) + GPnormal(n,4)/GPnormal(n,5)*walldist(n,i+4);
            RPcandy = IP(n,2) + walldist(n,i+4);
            RPcandz = IP(n,3) + GPnormal(n,6)/GPnormal(n,5)*walldist(n,i+4);
            if 1<=xgrid && xgrid<=Nx-1 && 1<=ygrid && ygrid<=Ny-1 && 1<=zgrid && zgrid<=Nz-1
                if inside(xgrid,ygrid,zgrid)==0 && inside(xgrid+1,ygrid,zgrid)==0 && inside(xgrid,ygrid,zgrid+1)==0 && inside(xgrid+1,ygrid,zgrid+1)==0
                    if abs(GPconv(n,1) - RPcandx) <= 2*dx && abs(GPconv(n,2) - RPcandy) <= 2*dy && abs(GPconv(n,3) - RPcandz) <= 2*dz
                        RPcand(n,i+4,1) = IP(n,1) + GPnormal(n,4)/GPnormal(n,5)*walldist(n,i+4);
                        RPcand(n,i+4,2) = IP(n,2) + walldist(n,i+4);
                        RPcand(n,i+4,3) = IP(n,3) + GPnormal(n,6)/GPnormal(n,5)*walldist(n,i+4);
                        RPcand(n,i+4,4) = IPdist + walldist(n,i+4)/GPnormal(n,5); 
                        RPcand(n,i+4,5) = 2;                                                 %a sign to tell that this SP is on the x-z plane
                    end
                end
            end
            %location of intersection points with x-y planes
            xgrid = floor((IP(n,1) + GPnormal(n,4)/GPnormal(n,6)*walldist(n,i+8))/dx);
            ygrid = floor((IP(n,2) + GPnormal(n,5)/GPnormal(n,6)*walldist(n,i+8))/dy);
            zgrid = int8((IP(n,3) + walldist(n,i+8))/dz);
            RPcandx = IP(n,1) + GPnormal(n,4)/GPnormal(n,6)*walldist(n,i+8);
            RPcandy = IP(n,2) + GPnormal(n,5)/GPnormal(n,6)*walldist(n,i+8);
            RPcandz = IP(n,3) + walldist(n,i+8);
            if 1<=xgrid && xgrid<=Nx-1 && 1<=ygrid && ygrid<=Ny-1 && 1<=zgrid && zgrid<=Nz-1
                if inside(xgrid,ygrid,zgrid)==0 && inside(xgrid+1,ygrid,zgrid)==0 && inside(xgrid,ygrid+1,zgrid)==0 && inside(xgrid+1,ygrid+1,zgrid)==0
                    if abs(GPconv(n,1) - RPcandx) <= 2*dx && abs(GPconv(n,2) - RPcandy) <= 2*dy && abs(GPconv(n,3) - RPcandz) <= 2*dz
                        RPcand(n,i+8,1) = IP(n,1) + GPnormal(n,4)/GPnormal(n,6)*walldist(n,i+8);
                        RPcand(n,i+8,2) = IP(n,2) + GPnormal(n,5)/GPnormal(n,6)*walldist(n,i+8);
                        RPcand(n,i+8,3) = IP(n,3) + walldist(n,i+8);
                        RPcand(n,i+8,4) = IPdist + walldist(n,i+8)/GPnormal(n,6); 
                        RPcand(n,i+8,5) = 3;                                                 %a sign to tell that this SP is on the x-y plane
                    end
                end
            end
        end
        %add PP to the list 
        RPcand(n,13,1:3) = PP(n,1:3);
        RPcand(n,13,4) = 0;
        RPcand(n,13,5) = 4; %a sign to tell which one is PP
        %add IP to the list 
        RPcand(n,14,1:3) = IP(n,1:3);
        RPcand(n,14,4) = sqrt((IP(n,1)-PP(n,1))^2 + (IP(n,2)-PP(n,2))^2 + (IP(n,3)-PP(n,3))^2);
        RPcand(n,14,5) = 5; %a sign to tell which one is IP
    end
    
    %arrange RPcand
    %sort the candidates based on xi (the distance from PP)
    RPcand_sorted = zeros(NGP,14,5);
    for n = 1:NGP
        RPcand_sorted(n,:,:) = sortrows(squeeze(RPcand(n,:,:)),4);
    end

    %-----------------------------------------------------------------
    %"RP" indicates the location of each reference points
    %RP(n,1,1:3) is the location of reference point 1
    %RP(n,2,1:3) is the location of image point
    %RP(n,3,1:3) is the location of reference point 2
    %RP(n,4,1:3) is the location of image point, only used when image point and 
    % intersection point are overlapping
    %RP(n,:,4) is xi(Distance along the normal vector from PP)
    %RP(n,:,5) is the sign to tell the classification of the point(IP,PP,etc)
    %RP(n,:,6) is the weight of each reference point

    RP = zeros(NGP,4,6);
    for n = 1:NGP
        for i = 2:13
            if RPcand_sorted(n,i,5) == 5 %detect the IP
                if abs(RPcand_sorted(n,i,4) - RPcand_sorted(n,i-1,4)) < dxi
                    RP(n,4,1:5) = squeeze(RPcand_sorted(n,i-1,:));
                elseif abs(RPcand_sorted(n,i,4) - RPcand_sorted(n,i+1,4)) < dxi
                    RP(n,4,1:5) = squeeze(RPcand_sorted(n,i+1,:));
                else
                    RP(n,1,1:5) = squeeze(RPcand_sorted(n,i-1,:));
                    RP(n,2,1:5) = squeeze(RPcand_sorted(n,i,:));
                    RP(n,3,1:5) = squeeze(RPcand_sorted(n,i+1,:));
                end
            end
        end
    end
    
    for n = 1:NGP
        if RP(n,4,5) ~= 0 %only when IP and RP are overlapping
            RP(n,4,6) = 1;%weight of that RP is 1
        else
            RP(n,1,6) = (RP(n,3,4) - RP(n,2,4))/(RP(n,3,4) - RP(n,1,4)); %weight of the RP1
            RP(n,3,6) = (RP(n,2,4) - RP(n,1,4))/(RP(n,3,4) - RP(n,1,4)); %weight of the RP2
        end
    end
    
    % ----------------------------------------------------------------
    %Output GPconv
    %the 3 first columns are for the coordinates of the GP, the 5x5x5 are the
    %convolution weights of the node points arount the image point, and the
    %final one is the value of the flow property on the object at the physical
    %point (PP)
    
    %weight of PP
    for n = 1:NGP
        if RP(n,1,5) == 4       %when PP is the RP1
            GPconv(n,129) = 2-RP(n,1,6);
        elseif RP(n,4,5) == 4   %when PP is overlapping with IP
            GPconv(n,129) = 1;
        else
            GPconv(n,129) = 2;
        end
    end
    
    %"RPrel" indicates the relative location of RP
    %The location of GP is used as the standard
    %the value is normalized by dx,dy and dz
    RPrel = zeros(NGP, 3, 5);
    for n = 1:NGP
        if RP(n,1,5) ~= 0
            RPrel(n,1,1) = (RP(n,1,1) - GPconv(n,1))/dx;
            RPrel(n,1,2) = (RP(n,1,2) - GPconv(n,2))/dy;
            RPrel(n,1,3) = (RP(n,1,3) - GPconv(n,3))/dz;
        end
        if RP(n,3,5) ~= 0
            RPrel(n,2,1) = (RP(n,3,1) - GPconv(n,1))/dx;
            RPrel(n,2,2) = (RP(n,3,2) - GPconv(n,2))/dy;
            RPrel(n,2,3) = (RP(n,3,3) - GPconv(n,3))/dz;
        end
        if RP(n,4,5) ~= 0 && RP(n,4,5) ~= 4
            RPrel(n,3,1) = (RP(n,4,1) - GPconv(n,1))/dx;
            RPrel(n,3,2) = (RP(n,4,2) - GPconv(n,2))/dy;
            RPrel(n,3,3) = (RP(n,4,3) - GPconv(n,3))/dz;
        end
        RPrel(n,1:3,4:5) = RP(n,[1,3,4],5:6);
    end
    
    %weight of RP which are overlapping with IP
    for n = 1:NGP
        if RP(n,4,5) ~= 0 && RP(n,4,5) ~= 4
            RPlocation = [int8(RPrel(n,3,1)),int8(RPrel(n,3,2)),int8(RPrel(n,3,3))];
            GPconv(n, (RPlocation(1)+2)*25 + (RPlocation(2)+2)*5 + (RPlocation(3)+2) + 4) = -1;
        end
    end
    
    %weight of the remaining points
    alpha = 0;
    beta = 0;
    RPa = [0,0,0];
    RPb = [0,0,0];
    RPc = [0,0,0];
    RPd = [0,0,0];
    for n = 1:NGP
        if RPrel(n,2,4) ~= 0
            for i = 1:2
                if RPrel(n,i,4) == 1
                    alpha = RPrel(n,i,2) - floor(RPrel(n,i,2));
                    beta = RPrel(n,i,3) - floor(RPrel(n,i,3));
                    RPa = [int8(RPrel(n,i,1)),floor(RPrel(n,i,2)),floor(RPrel(n,i,3))];
                    RPb = [int8(RPrel(n,i,1)),floor(RPrel(n,i,2)),floor(RPrel(n,i,3))+1];
                    RPc = [int8(RPrel(n,i,1)),floor(RPrel(n,i,2))+1,floor(RPrel(n,i,3))];
                    RPd = [int8(RPrel(n,i,1)),floor(RPrel(n,i,2))+1,floor(RPrel(n,i,3))+1];
                elseif RPrel(n,i,4) == 2
                    alpha = RPrel(n,i,1) - floor(RPrel(n,i,1));
                    beta = RPrel(n,i,3) - floor(RPrel(n,i,3));
                    RPa = [floor(RPrel(n,i,1)),int8(RPrel(n,i,2)),floor(RPrel(n,i,3))];
                    RPb = [floor(RPrel(n,i,1)),int8(RPrel(n,i,2)),floor(RPrel(n,i,3))+1];
                    RPc = [floor(RPrel(n,i,1))+1,int8(RPrel(n,i,2)),floor(RPrel(n,i,3))];
                    RPd = [floor(RPrel(n,i,1))+1,int8(RPrel(n,i,2)),floor(RPrel(n,i,3))+1];
                elseif RPrel(n,i,4) == 3
                    alpha = RPrel(n,i,1) - floor(RPrel(n,i,1));
                    beta = RPrel(n,i,2) - floor(RPrel(n,i,2));
                    RPa = [floor(RPrel(n,i,1)),floor(RPrel(n,i,2)),int8(RPrel(n,i,3))];
                    RPb = [floor(RPrel(n,i,1)),floor(RPrel(n,i,2))+1,int8(RPrel(n,i,3))];
                    RPc = [floor(RPrel(n,i,1))+1,floor(RPrel(n,i,2)),int8(RPrel(n,i,3))];
                    RPd = [floor(RPrel(n,i,1))+1,floor(RPrel(n,i,2))+1,int8(RPrel(n,i,3))];
                end
                if RPrel(n,i,4) == 1 || RPrel(n,i,4) == 2 || RPrel(n,i,4) == 3
                    GPconv(n, (RPa(1)+2)*25 + (RPa(2)+2)*5 + (RPa(3)+2) + 4) = GPconv(n, (RPa(1)+2)*25 + (RPa(2)+2)*5 + (RPa(3)+2) + 4) - (1-alpha)*(1-beta)*RPrel(n,i,5);
                    GPconv(n, (RPb(1)+2)*25 + (RPb(2)+2)*5 + (RPb(3)+2) + 4) = GPconv(n, (RPb(1)+2)*25 + (RPb(2)+2)*5 + (RPb(3)+2) + 4) - (1-alpha)*beta*RPrel(n,i,5);
                    GPconv(n, (RPc(1)+2)*25 + (RPc(2)+2)*5 + (RPc(3)+2) + 4) = GPconv(n, (RPc(1)+2)*25 + (RPc(2)+2)*5 + (RPc(3)+2) + 4) - alpha*(1-beta)*RPrel(n,i,5);
                    GPconv(n, (RPd(1)+2)*25 + (RPd(2)+2)*5 + (RPd(3)+2) + 4) = GPconv(n, (RPd(1)+2)*25 + (RPd(2)+2)*5 + (RPd(3)+2) + 4) - alpha*beta*RPrel(n,i,5);
                end
            end
        elseif RPrel(n,2,4) == 0 && RPrel(n,3,4) == 0
            GPconv(n,4:129) = 0;
            penalization_count(uvw) = penalization_count(uvw) + 1;
        end
    end
    % ----------------------------------------------------------------
    %this part makes output files for DNS
    %shifting the coordinates to fit the DNS grid

    for n = 1:NGP
        GPconv(n,1) = GPconv(n,1)/dx + xshift;
        GPconv(n,2) = GPconv(n,2)/dy + yshift;
        GPconv(n,3) = GPconv(n,3)/dz + zshift;
    end

    %write the convolution kernel in txt files
    fileID = fopen (fname, 'w');
    formatSpec = ['%d %d %d' repmat('%8.4f ', 1, 125) '%8.4f\n'];
    fprintf(fileID,'%d\n',NGP); %1st line is the number of GP
    fprintf(fileID,formatSpec,transpose(GPconv));
    fclose(fileID);

    %write the grid parameters in a csv file
    writematrix(grid,"grid.csv")
end