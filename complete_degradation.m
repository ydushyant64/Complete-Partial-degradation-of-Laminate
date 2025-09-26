%% 
clc

% Material properties of glass/epoxy lamina
E1 = 38.6e9;
E2 = 8.27e9;
v12 = 0.26;
G12 = 4.14e9;

v21 = (E2/E1) * v12;

sigma_1T = 1062e6;
sigma_1C = -610e6;
sigma_2T = 31e6;
sigma_2C = -118e6;
tau_12 = 72e6;

% Hygroscopic and thermal properties of glass/epoxy lamina
alpha_1 = 8.6e-6;   % m/m/°C
alpha_2 = 22.1e-6;  % m/m/°C
beta_1 = 0;         % m/m/kg/kg
beta_2 = 0.6;       % m/m/kg/kg 

% Base stiffness matrix components
Q11 = E1 / (1 - v12 * v21);
Q12 = (v12 * E2) / (1 - v12 * v21);
Q22 = E2 / (1 - v12 * v21);
Q66 = G12;

angles_deg_initial = [0, 45, -45, 90, 90, -45, 45, 0];
n_plies = length(angles_deg_initial);

t_ply = 0.125e-3;  % Ply thickness [m]
total_thickness = n_plies * t_ply;
h_initial = linspace(-total_thickness/2, total_thickness/2, n_plies + 1);

% Mechanical Loading
Nx = 1000; Ny = 0; Nxy = 0; Mx = 0; My = 0; Mxy = 0;
f_applied = [Nx; Ny; Nxy; Mx; My; Mxy];

angles_deg = angles_deg_initial;
h = h_initial;
f_total = f_applied;

iteration = 1;
max_iteration = 100;

% Track failed plies
failed_plies = false(n_plies, 1);

while iteration <= max_iteration
    fprintf('\n================================== Iteration %d ==================================\n', iteration);
    
    % Number of plies (fixed)
    n = length(angles_deg);
    
    % Total thickness (fixed)
    total_thickness = n_plies * t_ply;
    h = h_initial;
    
    %disp('Updated h:');
    %disp(h);

    Qcell = cell(1, n);   % To store Q_ matrices
    alpha_global_cell = cell(1, n);
    beta_global_cell = cell(1, n);
    
    for i = 1:n
        if failed_plies(i)
            Qcell{i} = zeros(3,3); % Fully degraded ply
            continue;
        end
        
        x_deg = angles_deg(i);
        x_rad = deg2rad(x_deg);

        c = cos(x_rad);
        s = sin(x_rad);

        Q11_ = Q11 * c^4 + Q22 * s^4 + 2 * (Q12 + 2 * Q66) * (s^2) * (c^2);
        Q12_ = (Q11 + Q22 - 4 * Q66) * (s^2) * (c^2) + Q12 * (c^4 + s^4);
        Q22_ = Q11 * s^4 + Q22 * c^4 + 2 * (Q12 + 2 * Q66) * (s^2) * (c^2);
        Q16_ = (Q11 - Q12 - 2 * Q66) * s * c^3 - (Q22 - Q12 - 2 * Q66) * s^3 * c;
        Q26_ = (Q11 - Q12 - 2 * Q66) * c * s^3 - (Q22 - Q12 - 2 * Q66) * c^3 * s;
        Q66_ = (Q11 + Q22 - 2 * Q12 - 2 * Q66) * (s^2) * (c^2) + Q66 * (s^4 + c^4);

        Q_ = [Q11_, Q12_, Q16_;
              Q12_, Q22_, Q26_;
              Q16_, Q26_, Q66_];

        % Handle numerical issues
        if any(isnan(Q_(:)) | isinf(Q_(:)))
            Q_ = zeros(3,3);
        end

        Qcell{i} = Q_;
    end

    A = zeros(3);
    B_temp = zeros(3);
    D_temp = zeros(3);

    % Compute A, B, D matrices
    for k = 1:n
        delta_h = h(k+1) - h(k);
        A = A + Qcell{k} * delta_h;
        delta_h2 = (h(k+1))^2 - (h(k))^2;
        B_temp = B_temp + Qcell{k} * delta_h2;
        delta_h3 = (h(k+1))^3 - (h(k))^3;
        D_temp = D_temp + Qcell{k} * delta_h3;
    end

    B = B_temp/2;
    D = D_temp/3;

    ABBD = [A, B; B, D];
    fprintf('ABBD matrix:\n');
    disp(ABBD);

    % Check if ABBD is singular
    if rank(ABBD) < 6
        fprintf('ABBD matrix is singular. Laminate has failed. Exiting loop.\n');
        break;
    end

    for i = 1:n
        theta = deg2rad(angles_deg(i));
        c = cos(theta);
        s = sin(theta);
    
        T_inv = [c^2, s^2, -2*c*s;
                 s^2, c^2, 2*c*s;
                 c*s, -c*s, c^2-s^2];
    
        alpha_local = [alpha_1; alpha_2; 0];
     
        alpha_global_cell{i} = T_inv * alpha_local;
       beta_local = [beta_1; beta_2; 0];
        beta_global_cell{i} = T_inv * beta_local;
    end
    
    deltaT = 50;   % Temperature rise in °C
    deltaC = 0.5;  % Moisture content rise in %

    N_thermal = zeros(3,1);
    M_thermal = zeros(3,1);
    N_moisture = zeros(3,1);
    M_moisture = zeros(3,1);

    for k = 1:n
        delta_h = h(k+1) - h(k);
        delta_h2 = (h(k+1)^2 - h(k)^2)/2;
    
        N_thermal = N_thermal + Qcell{k} * alpha_global_cell{k} * deltaT * delta_h;
        M_thermal = M_thermal + Qcell{k} * alpha_global_cell{k} * deltaT * delta_h2;
        N_moisture = N_moisture + Qcell{k} * beta_global_cell{k} * deltaC * delta_h;
        M_moisture = M_moisture + Qcell{k} * beta_global_cell{k} * deltaC * delta_h2;
    end

    f_total = f_applied - [N_thermal; M_thermal] - [N_moisture; M_moisture];

    ABBD_inv = inv(ABBD);
    R = ABBD_inv * f_total;

    strain_curvature.ex0  = R(1);
    strain_curvature.ey0  = R(2);
    strain_curvature.exy0 = R(3);
    strain_curvature.kx   = R(4);
    strain_curvature.ky   = R(5);
    strain_curvature.kxy  = R(6);

    fprintf('\nLocal Strain and Curvature result:\n');
    disp(strain_curvature);

    eps_xcell = cell(n, 1);
    for i = 1:n
        eps0 = [strain_curvature.ex0; strain_curvature.ey0; strain_curvature.exy0];
        K = [strain_curvature.kx; strain_curvature.ky; strain_curvature.kxy];
        z_mid = (h(i+1) + h(i)) / 2;
        eps_global = eps0 + z_mid * K;
        strain_global.ex = eps_global(1);
        strain_global.ey = eps_global(2);
        strain_global.exy = 0.5*eps_global(3);
        eps_xcell{i} = [strain_global.ex; strain_global.ey; strain_global.exy];
    end

    eps_1cell = cell(n, 1);
    eps_2cell = cell(n, 1);
    for i = 1:n
        x_deg = angles_deg(i);
        x_rad = deg2rad(x_deg);
        c = cos(x_rad);
        s = sin(x_rad);
        T = [c^2, s^2, 2*s*c;
             s^2, c^2, -2*s*c;
             -s*c, s*c, c^2-s^2];
        eps_local = T*eps_xcell{i};
        strain_local.e1 = eps_local(1);
        strain_local.e2 = eps_local(2);
        strain_local.e12 = 2*eps_local(3);
        strain_pp.e12 = eps_local(3);
        eps_1cell{i} = [strain_local.e1; strain_local.e2; strain_local.e12];
        eps_2cell{i} = [strain_local.e1; strain_local.e2; strain_pp.e12];
    end

    sigmacell = cell(n, 1);
    SR_layer = cell(n, 1);
    failure = cell(n, 1);
    for i = 1:n
        if failed_plies(i)
            SR_layer{i} = zeros(3,1);
            failure{i} = zeros(3,1);
            fprintf('Layer %d is fully degraded\n', i);
            continue;
        end
        
        sigma = Qcell{i}*eps_2cell{i};
        stress.s1 = sigma(1);
        stress.s2 = sigma(2);
        stress.s3 = sigma(3);
        sigmacell{i} = [stress.s1; stress.s2; stress.s3];
        fprintf('\nLocal Stress at %d° (layer %d):\n', angles_deg(i), i);
        disp(stress);

        SR = zeros(3,1);
        failure_mode = zeros(3,1);
    
        if stress.s1 > 0
            SR(1) = stress.s1/sigma_1T;
            failure_mode(1) = 1;
            fprintf('Stress ratio (Longitudinal-Tension): %f\n', SR(1));
        else
            SR(1) = stress.s1/sigma_1C;
            failure_mode(1) = 1;
            fprintf('Stress ratio (Longitudinal-Compression): %f\n', SR(1));
        end
        
        if stress.s2 > 0
            SR(2) = stress.s2/sigma_2T;
            failure_mode(2) = 2;
            fprintf('Stress ratio (Transverse-Tension): %f\n', SR(2));
        else
            SR(2) = stress.s2/sigma_2C;
            failure_mode(2) = 2;
            fprintf('Stress ratio (Transverse-Compression): %f\n', SR(2));
        end

        SR(3) = abs(stress.s3/tau_12);
        failure_mode(3) = 3;
        fprintf('Stress ratio (Shear): %f\n\n', SR(3));

        SR_layer{i} = SR;
        failure{i} = failure_mode;
    end

    max_SR_layer = cell(n, 1);
    for i = 1:n
        max_SR = max(SR_layer{i});
        fprintf('Maximum stress ratio in layer %d: %.4f\n', i, max_SR);
        max_SR_layer{i} = max_SR;
    end

    max_SR_values = cell2mat(max_SR_layer);
    max_SR_global = max(max_SR_values);
    failed_layers = find(abs(max_SR_values - max_SR_global) < 1e-6);
    failed_angles = angles_deg(failed_layers);

    fprintf('\nOverall maximum stress ratio: %.4f\n', max_SR_global);
    for i = 1:length(failed_layers)
        fprintf('Failure occurs at Layer %d with Orientation %d°\n', failed_layers(i), failed_angles(i));
    end

    % Mark failed layers as fully degraded
    if max_SR_global >= 1
        for i = failed_layers
            failed_plies(i) = true;
            fprintf('Layer %d marked as fully degraded\n', i);
        end
    else
        fprintf('No failure in this iteration (max SR < 1).\n');
    end

    fpf_load = Nx/max_SR_global;
    fprintf('First ply failure load is %f N/m\n', fpf_load);

    % Termination conditions
    if all(failed_plies)
        fprintf('All plies are fully degraded. Exiting loop after %d iterations.\n', iteration);
        break;
    end

    iteration = iteration + 1;
end