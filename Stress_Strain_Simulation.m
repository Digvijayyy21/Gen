% Motion: x1 = X1 + 3*t*X2, x2 = X2, x3 = X3,  t = 0…10

clc;
clear all; 

% MONOCLINIC MATERIAL STIFFNESS (Voigt Form)

C11=150e3; C22=120e3; C33=100e3;
C44=30e3;  C55=40e3;  C66=50e3;
C12=60e3;  C13=45e3;  C23=50e3;
C15=10e3;  C25=8e3;   C35=6e3; 
C46=5e3;

C = zeros(6,6);
C(1,1)=C11; C(2,2)=C22; C(3,3)=C33;
C(4,4)=C44; C(5,5)=C55; C(6,6)=C66;

C(1,2)=C12; C(2,1)=C12;
C(1,3)=C13; C(3,1)=C13;
C(2,3)=C23; C(3,2)=C23;
C(1,5)=C15; C(5,1)=C15;
C(2,5)=C25; C(5,2)=C25;
C(3,5)=C35; C(5,3)=C35;
C(4,6)=C46; C(6,4)=C46;

% Tensor action function

C_action = @(D) voigt_to_tensor(C * tensor_to_voigt(D));

% Define deformation gradient history

t_vals = 0:0.01:10;
n_steps = length(t_vals)-1;

F = zeros(3,3,length(t_vals));
for k = 1:length(t_vals)
    t = t_vals(k);
    F(:,:,k) = [1, 3*t, 0;
                0, 1,   0;
                0, 0,   1];
end

DeltaF = F(:,:,2) - F(:,:,1);
dt = 1;

% Initialize variables

sigma_small  = zeros(3,3);   % Cauchy stress (small strain)
sigma_large  = zeros(3,3);   % Cauchy stress (Jaumann)

eps_small    = zeros(3,3);   % accumulated small strain

vm_small  = zeros(length(t_vals),1);
vm_large  = zeros(length(t_vals),1);

eq_small = zeros(length(t_vals),1);
eq_large = zeros(length(t_vals),1);

% initial values
vm_small(1) = vonMises(sigma_small);
vm_large(1) = vonMises(sigma_large);

E_GL = 0.5*(F(:,:,1).'*F(:,:,1) - eye(3));
eq_large(1) = eqStrain_large(E_GL);
eq_small(1) = 0;

%% MAIN LOOP

for n = 1:n_steps

    Fn   = F(:,:,n);
    Fnp1 = F(:,:,n+1);
    dF   = DeltaF;

    Fnp1_inv = inv(Fnp1);

    %% Compute velocity gradient and its parts
    l = dF * Fnp1_inv;
    d = 0.5*(l + l.');
    W = 0.5*(l - l.');

    %% SMALL STRAIN 
    % sigma_dot = C : d
    sigma_dot_small = C_action(d);
    sigma_small = sigma_small + sigma_dot_small * dt;

    % accumulate engineering strain
    eps_small = eps_small + d*dt;

    %% LARGE STRAIN 
    % Jaumann rate: sigma_dot = C:d - sigma*W + W*sigma
    
    sigma_dot_large = C_action(d) - sigma_large*W + W*sigma_large;
    sigma_large = sigma_large + sigma_dot_large * dt;

    %% Equivalent strain 
    idx = n+1;

    vm_small(idx) = vonMises(sigma_small);
    vm_large(idx) = vonMises(sigma_large);

    eq_small(idx) = eqStrain_small(eps_small);

    E_GL = 0.5*(Fnp1.'*Fnp1 - eye(3));
    eq_large(idx) = eqStrain_large(E_GL);
end


% Plot STRESS–STRAIN curves

figure;
plot(eq_small, vm_small,'-', 'linewidth', 2,'DisplayName','Small strain');
hold on;
plot(eq_large, vm_large,'-','linewidth', 2, 'DisplayName','Large strain (Jaumann)');
grid on;
xlabel('Equivalent Strain');
ylabel('Von Mises Stress');
title('Monoclinic Material: Small vs Large Strain Formulation');
legend('Location','northwest');

%% Helper Functions
function v = tensor_to_voigt(T)
    v = [T(1,1);
         T(2,2);
         T(3,3);
         2*T(2,3);
         2*T(1,3);
         2*T(1,2)];
end

function T = voigt_to_tensor(v)
    T = [ v(1),   v(6)/2, v(5)/2;
          v(6)/2, v(2),   v(4)/2;
          v(5)/2, v(4)/2, v(3) ];
end

function vm = vonMises(s)
    dev = s - trace(s)/3 * eye(3);
    vm = sqrt(1.5 * sum(sum(dev.*dev)));
end

function eq = eqStrain_small(eps)
    ed = eps - trace(eps)/3*eye(3);
    eq = sqrt(2/3 * sum(sum(ed.*ed)));
end

function eq = eqStrain_large(E)
    Ed = E - trace(E)/3*eye(3);
    eq = sqrt(2/3 * sum(sum(Ed.*Ed)));
end