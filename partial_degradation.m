clc

% Material properties of glass/epoxy lamina
E1_initial = 38.6e9;
E2_initial = 8.27e9;
v12_initial = 0.26;
G12_initial = 4.14e9;

v21_initial = (E2_initial/E1_initial) * v12_initial;

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
max_iteration = 10;

% Initialize material properties for each ply
E1 = E1_initial * ones(1, n_plies);
E2 = E2_initial * ones(1, n_plies);
v12 = v12_initial * ones(1, n_plies);
G12 = G12_initial * ones(1, n_plies);
v21 = v21_initial * ones(1, n_plies);

% Track failure status for each mode (1: longitudinal, 2: transverse, 3: shear)
failed_modes = false(n_plies, 3); % [longitudinal, transverse, shear]

while (iteration <= max_iteration)
    
    fprintf('\n================================== Iteration %d ==================================\n', iteration);
    
    % Number of plies
    n = length(angles_deg);
    
    % Total thickness
    total_thickness = n * t_ply;
    h = linspace(-total_thickness/2, total_thickness/2, n + 1);
    
    %disp('Updated h:');
    %disp(h);

    Qcell = cell(1, n);   % To store Q_ matrices
    alpha_global_cell = cell(1, n);
    beta_global_cell = cell(1, n);
    
    for i = 1:n
        % Calculate Q matrix for each ply with its current properties
        Q11 = E1(i) / (1 - v12(i) * v21(i));
        Q12 = (v12(i) * E2(i)) / (1 - v12(i) * v21(i));
        Q22 = E2(i) / (1 - v12(i) * v21(i));
        Q66 = G12(i);
        
        % Handle numerical issues
        if isnan(Q11) || isinf(Q11), Q11 = 0; end
        if isnan(Q12) || isinf(Q12), Q12 = 0; end
        if isnan(Q22) || isinf(Q22), Q22 = 0; end
        
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

    for i = 1:n
        theta = deg2rad(angles_deg(i));
        c = cos(theta);
        s = sin(theta);
    
        T_inv = [c^2, s^2, -2*c*s;
                 s^2, c^2, 2*c*s;
                 c*s, -c*s, c^2-s^2];
    
        alpha_local = [alpha_1; alpha_2; 0];
        beta_local  = [beta_1; beta_2; 0];
    
        alpha_global_cell{i} = T_inv * alpha_local;
        beta_global_cell{i}  = T_inv * beta_local;
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

    % Check if ABBD is singular
    if rank(ABBD) < 6
        fprintf('ABBD matrix is singular. Laminate has failed. Exiting loop.\n');
        break;
    end

    ABBD_inv = inv(ABBD);
    R = ABBD_inv * f_total;

    strain_curvature.ex0  = R(1);
    strain_curvature.ey0  = R(2);
    strain_curvature.exy0 = R(3);
    strain_curvature.kx   = R(4);
    strain_curvature.ky   = R(5);
    strain_curvature.kxy  = R(6);

    %fprintf('\nLocal Strain and Curvature result:\n');
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
        sigma = Qcell{i}*eps_2cell{i};
        stress.s1 = sigma(1);
        stress.s2 = sigma(2);
        stress.s3 = sigma(3);
        sigmacell{i} = [stress.s1; stress.s2; stress.s3];
        fprintf('\nLocal Stress at %d° (layer %d):\n', angles_deg(i), i);
        disp(stress);

        SR = zeros(3,1);
        failure_mode = zeros(3,1);
    
        % Longitudinal failure check
        if ~failed_modes(i, 1) % Skip if already failed longitudinally
            if stress.s1 > 0
                SR(1) = stress.s1/sigma_1T;
                failure_mode(1) = 1;
                fprintf('Stress ratio (Longitudinal-Tension): %f\n', SR(1));
            else
                SR(1) = stress.s1/sigma_1C;
                failure_mode(1) = 1;
                fprintf('Stress ratio (Longitudinal-Compression): %f\n', SR(1));
            end
        else
            fprintf('Longitudinal mode already failed for layer %d\n', i);
        end
        
        % Transverse failure check
        if ~failed_modes(i, 2) % Skip if already failed transversely
            if stress.s2 > 0
                SR(2) = stress.s2/sigma_2T;
                failure_mode(2) = 2;
                fprintf('Stress ratio (Transverse-Tension): %f\n', SR(2));
            else
                SR(2) = stress.s2/sigma_2C;
                failure_mode(2) = 2;
                fprintf('Stress ratio (Transverse-Compression): %f\n', SR(2));
            end
        else
            fprintf('Transverse mode already failed for layer %d\n', i);
        end

        % Shear failure check
        if ~failed_modes(i, 3) % Skip if already failed in shear
            SR(3) = abs(stress.s3/tau_12);
            failure_mode(3) = 3;
            fprintf('Stress ratio (Shear): %f\n\n', SR(3));
        else
            fprintf('Shear mode already failed for layer %d\n', i);
        end

        SR_layer{i} = SR;
        failure{i} = failure_mode;
    end
    
    % Store previous properties to check for changes
    E1_prev = E1;
    E2_prev = E2;
    G12_prev = G12;

    % Find maximum SR across all layers
    max_SR_global = 0;
    failed_plies = [];
    overall_mode = 0;

    for i = 1:n
        [max_SR, idx] = max(SR_layer{i});
        if max_SR > max_SR_global
            max_SR_global = max_SR;
            failed_plies = [i];
            overall_mode = failure{i}(idx);
        elseif abs(max_SR - max_SR_global) < 1e-6
            failed_plies = [failed_plies, i];
            if overall_mode == 0
                overall_mode = failure{i}(idx);
            end
        end       
    end

    fprintf('\nOverall maximum stress ratio: %.4f\n', max_SR_global);
    for i = 1:length(failed_plies)
        fprintf('Failure occurs at Layer %d (%d°) with SR = %.4f and failure mode = %d\n', ...
                failed_plies(i), angles_deg(failed_plies(i)), max_SR_global, overall_mode);
    end
    
    % Partial degradation: Update material properties for all failed plies
    if max_SR_global >= 1
        if overall_mode == 1
            fprintf('Failure is in longitudinal direction\n');
            for i = failed_plies
                if ~failed_modes(i, 1)
                    E1(i) = 0;
                    failed_modes(i, 1) = true;
                    fprintf('E1 for layer %d = %f\n', i, E1(i));
                end
            end
        elseif overall_mode == 2
            fprintf('Failure is in transverse direction\n');
            for i = failed_plies
                if ~failed_modes(i, 2)
                    E2(i) = 0;
                    G12(i) = 0;
                    failed_modes(i, 2) = true;
                    failed_modes(i, 3) = true; % Transverse failure also disables shear
                    fprintf('E2 for layer %d = %f\nG12 for layer %d = %f\n', ...
                            i, E2(i), i, G12(i));
                end
            end
        else
            fprintf('Shear failure\n');
            for i = failed_plies
                if ~failed_modes(i, 3)
                    G12(i) = 0;
                    failed_modes(i, 3) = true;
                    fprintf('G12 for layer %d = %f\n', i, G12(i));
                end
            end
        end

        % Update v21 for all failed plies
        for i = failed_plies
            v21(i) = (E2(i)/E1(i)) * v12(i);
            if isnan(v21(i)) || isinf(v21(i))
                v21(i) = 0;
            end
        end
    else
        fprintf('No failure in this iteration (max SR < 1).\n');
    end

    % Calculate first ply failure load
    fpf_load = Nx / max_SR_global;
    fprintf('First ply failure load is %f N/m\n', fpf_load);

    % Check termination conditions
    % 1. Check if no properties changed (no new failures)
    properties_changed = any(E1 ~= E1_prev | E2 ~= E2_prev | G12 ~= G12_prev);
    
    % 2. Check if all plies have at least one property degraded
    all_degraded = all(E1 == 0 | (E2 == 0 & G12 == 0));

    if ~properties_changed && all_degraded
        fprintf('No new failures and all plies have at least one property degraded. Exiting loop.\n');
        break;
    end

    iteration = iteration + 1;
end