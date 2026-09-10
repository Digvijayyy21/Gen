clear; clc; close all;

%% PARAMETERS 
Nx = 3; Ny = 3; Nz = 6;         
E = 200e9; A = 5e-4; I = 5e-8;   % Structural steel
G = E/2.6; J = 2*I;               
L = 1.0;                       
ndof_per_node = 6;             

% Wind profile parameters
fprintf('Input reference wind velocity (m/s): ');
v_ref = input('');             
z_ref = 10; alpha = 0.14; air_density = 1.225; A_proj = 0.5;  % Wind effect

%% STEP 1-2: Generate nodes and connectivity
[xv,yv,zv] = ndgrid(0:Nx-1, 0:Ny-1, 0:Nz-1);
coords_original = [xv(:), yv(:), zv(:)];    
numnp = size(coords_original,1);

elem = [];
for ix = 1:Nx-1; for iy = 1:Ny; for iz = 1:Nz
    n1 = sub2ind([Nx,Ny,Nz], ix,   iy, iz); n2 = sub2ind([Nx,Ny,Nz], ix+1, iy, iz); elem = [elem; n1 n2];
end; end; end
for ix = 1:Nx; for iy = 1:Ny-1; for iz = 1:Nz
    n1 = sub2ind([Nx,Ny,Nz], ix, iy,   iz); n2 = sub2ind([Nx,Ny,Nz], ix, iy+1, iz); elem = [elem; n1 n2];
end; end; end
for ix = 1:Nx; for iy = 1:Ny; for iz = 1:Nz-1
    n1 = sub2ind([Nx,Ny,Nz], ix, iy, iz  ); n2 = sub2ind([Nx,Ny,Nz], ix, iy, iz+1); elem = [elem; n1 n2];
end; end; end

nelem = min(117, size(elem,1)); elem = elem(1:nelem,:);

%% STEP 3-4: Matrices and wind profile 
R_x = eye(3); R_y = [0 0 0; 0 1 0; 0 0 1]; R_z = [0 0 0; 0 0 1; 1 0 0];
T_x = blkdiag(R_x, R_x, R_x, R_x); T_y = blkdiag(R_y, R_y, R_y, R_y); T_z = blkdiag(R_z, R_z, R_z, R_z);

ea = E*A/L; phi = 12*E*I/(L^3); psi = 6*E*I/(L^2); gj = G*J/L;
k_local = zeros(12,12);
k_local([1,7],[1,7]) = ea*[1 -1; -1 1];
k_local([2,3,8,9],[2,3,8,9]) = [phi psi 0 -phi; psi 4*E*I/L 0 -psi; 0 0 phi -psi; -phi -psi -psi phi];
k_local([4,5,10,11],[4,5,10,11]) = [phi -psi 0 -phi; -psi 4*E*I/L 0 psi; 0 0 phi psi; -phi psi psi phi];
k_local([6,12],[6,12]) = gj*[1 -1; -1 1];

wind_velocity = @(z) v_ref * max(z / z_ref, 0.1).^alpha;
pressure_at_z = @(z) 0.5 * air_density * wind_velocity(z).^2;

%% STEP 5: Assembly
ndof = numnp * ndof_per_node; Kglobal = zeros(ndof); Fglobal = zeros(ndof,1);
for e = 1:nelem
    n1 = elem(e,1); n2 = elem(e,2); z_mean = mean([coords_original(n1,3), coords_original(n2,3)]);
    dx = coords_original(n2,1) - coords_original(n1,1); dy = coords_original(n2,2) - coords_original(n1,2); dz = coords_original(n2,3) - coords_original(n1,3);
    
    if dx==1 && dy==0 && dz==0, T = T_x; elseif dx==0 && dy==1 && dz==0, T = T_y; elseif dx==0 && dy==0 && dz==1, T = T_z; else continue; end
    
    Ke = T' * k_local * T; 
    dofs = [6*(n1-1)+(1:6), 6*(n2-1)+(1:6)];
    for i = 1:12, for j = 1:12, Kglobal(dofs(i), dofs(j)) = Kglobal(dofs(i), dofs(j)) + Ke(i,j); end; end
    
    if isOuterElement3D(e, elem, coords_original, Nx, Ny, Nz)
        p = pressure_at_z(z_mean); w = p * A_proj;
        f_local_wind = zeros(12,1); f_local_wind([2,8]) = w*L/2; f_local_wind([3,9]) = w*L^2/12 * [1; -1];
        f_global = T' * f_local_wind;
        for k = 1:12, Fglobal(dofs(k)) = Fglobal(dofs(k)) + f_global(k); end
    end
end

%% STEP 6: Solve 
bc_dofs = 1:54;

K_red = Kglobal(bc_dofs(end)+1:end, bc_dofs(end)+1:end);
F_red = Fglobal(bc_dofs(end)+1:end);
F_red(46) = 15;
F_red(47) = 15;
F_red(48) = 15;
F_red(49) = 15;
F_red(50) = 15;
F_red(51) = 15;
F_red(52) = 15;
F_red(53) = 15;
F_red(54) = 15;
U_red = K_red\F_red;
U = zeros(ndof,1); U(bc_dofs(end)+1:end) = U_red;


%% STEP 7: AUTO-SCALE DEFORMATION FOR VISIBILITY
max_disp = max(abs(U(1:3:ndof)));  % Max translation
struct_size = max(coords_original(:)) - min(coords_original(:));  % Structure size
scale_factor = max(50, struct_size / max_disp * 0.2);  % Auto-scale: 20% of structure size

coords_deformed = coords_original;
for n = 1:numnp
    ux = U(6*(n-1)+1); uy = U(6*(n-1)+2); uz = U(6*(n-1)+3);
    coords_deformed(n,:) = coords_original(n,:) + scale_factor * [ux, uy, uz];
end

%% STEP 8: ENHANCED PLOTTING WITH VISIBLE DEFORMATION
figure('Position', [50 50 1600 700]);

% Original structure
subplot(1,3,1); hold on; axis equal;
for e = 1:nelem
    n1 = elem(e,1); n2 = elem(e,2);
    plot3(coords_original([n1,n2],1), coords_original([n1,n2],2), coords_original([n1,n2],3), ...
          'b-', 'LineWidth', 2);
end
grid on; xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
view(45,20);

% Deformed structure ONLY
subplot(1,3,2); hold on; axis equal;
for e = 1:nelem
    n1 = elem(e,1); n2 = elem(e,2);
    plot3(coords_deformed([n1,n2],1), coords_deformed([n1,n2],2), coords_deformed([n1,n2],3), ...
          'r-', 'LineWidth', 3);
end
grid on; xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)'); 
view(45,20); 

% OVERLAY 
subplot(1,1,1); hold on; axis equal;
for e = 1:nelem
    n1 = elem(e,1); n2 = elem(e,2);
    plot3(coords_original([n1,n2],1), coords_original([n1,n2],2), coords_original([n1,n2],3), 'b-', 'LineWidth', 1.5);
    plot3(coords_deformed([n1,n2],1), coords_deformed([n1,n2],2), coords_deformed([n1,n2],3), 'r-', 'LineWidth', 2.5);
end
grid on; xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title(sprintf('Overlay', scale_factor)); 
view(45,20);
legend('Original', 'Deformed', 'Location', 'best');

sgtitle(sprintf('3D Frame Wind Analysis: %.1f m/s → %.1f mm deflection (x%.0f)', v_ref, 1e3*max_disp, scale_factor));

%% Helper function
function outer = isOuterElement3D(e, elem, coords, Nx, Ny, Nz)
    n1 = elem(e,1); n2 = elem(e,2);
    [ix1,iy1,iz1] = ind2sub([Nx,Ny,Nz], n1);
    [ix2,iy2,iz2] = ind2sub([Nx,Ny,Nz], n2);
    outer = (ix1==1 || ix1==Nx || ix2==1 || ix2==Nx || iy1==1 || iy1==Ny || iy2==1 || iy2==Ny || iz1==1 || iz1==Nz || iz2==1 || iz2==Nz);
end