function [X, ctrl] = simNonLinearBikeModel(car, frontTires, rearTires, path, Ux_des, X0, t_final, dT, useFF)
    
    % time
    t_s = 0:dT:t_final;
    N = length(t_s);
    
    % initialize state vector
    U_x = zeros(N,1);
    U_y = zeros(N,1);
    r = zeros(N,1);
    e = zeros(N,1);
    s = zeros(N,1);
    dpsi = zeros(N,1);
    a_y = zeros(N,1);
    a_x = zeros(N,1);
    outUx_des = zeros(N,1);
    
    % Set initial condition
    U_x(1) = X0(1);
    U_y(1) = X0(2);
    r(1) = X0(3);
    e(1) = X0(4);
    s(1) = X0(5);
    dpsi(1) = X0(6);
    
    % initialize control vector
    delta = zeros(N,1);
    F_xtotal = zeros(N,1);
   
   
    for idx = 1 : N
        
        % Look up Ux_des
        currUx_des = interp1(path.s_m, Ux_des, s(idx));
        outUx_des(idx) = currUx_des;
        
        % Look up Curvature
        K = interp1(path.s_m, path.k_1pm, s(idx));

        [delta(idx), F_xtotal(idx)] = me227_controller(s(idx), e(idx), dpsi(idx), U_x(idx), U_y(idx), r(idx), 2, path);

% %         Calculate Control inputs
% %         Calpha1 = 188000; % N/rad
% %         xla = 15; % m
% %         Kla = 3500; % N/m
% %         Calpha2 = 203000; % N/rad
% %         dpsi_ss = K * (((car.m * car.a * U_x(idx)^2) / (car.L * Calpha2)) - car.b);
% %         if useFF
% %             delta_ff = ((Kla * xla * dpsi_ss)/car.C_alphaf) + K * (car.L + car.k_rad * U_x(idx)^2);
% %         else
% %             delta_ff = 0;
% %         end
% %         delta(idx) = (-Kla / Calpha1) * (e(idx) + xla * dpsi(idx)) + delta_ff;
% %         F_xtotal(idx) = car.Kdrive * (currUx_des - U_x(idx));
        
        % Step 1: calculate slip angles from velocities and steer angles
        alphaf = atan2((U_y(idx) + car.a * r(idx)) , U_x(idx)) - delta(idx);
        alphar = atan2((U_y(idx) - car.b * r(idx)) , U_x(idx));

        % Step 2: calculate forces from slip angles
        Fz = car.m * car.g;
        F_yf = computeTireForce(frontTires, Fz, alphaf);
        F_yr = computeTireForce(rearTires, Fz, alphar);
        F_xf = 0.6*F_xtotal(idx);
        F_xr = 0.4*F_xtotal(idx);
        
        state.s_m = s(idx);
        state.e_m = e(idx);
        state.deltaPsi_rad = dpsi(idx);
        state.Ux_mps = U_x(idx);
        state.Uy_mps = U_y(idx);
        state.r_radps = r(idx);
        
        % Computing external forces on the vehicle
        [Fdrag_N, Frr_N] = computeExternalForces(car, state);
        Fgrade_N = 0.05 * car.W * sin(0.5 * rand(1));  % Optional unmodeled disturbance due to grade
        Fext_N = Fdrag_N + Frr_N + Fgrade_N;

        % Step 3: calculate derivatives of state variables
        U_x_dot = ((-Fext_N + F_xr + F_xf * cos(delta(idx)) - F_yf * sin(delta(idx))) / car.m) + r(idx) * U_y(idx);
        U_y_dot = ((F_yf * cos(delta(idx)) + F_yr + F_xf * sin(delta(idx))) / car.m) - r(idx) * U_x(idx);
        r_dot = (car.a * F_yf * cos(delta(idx)) + car.a * F_xf * sin(delta(idx)) - car.b * F_yr) / car.Iz;
        s_dot = (1 / (1 - e(idx) * K) ) * ( U_x(idx) * cos(dpsi(idx)) - U_y(idx) * sin(dpsi(idx)));
        e_dot = U_y(idx) * cos(dpsi(idx)) + U_x(idx) * sin(dpsi(idx));
        dpsi_dot = r(idx) - K * s_dot;

        % Add noise
        %U_x(idx) = U_x(idx) + randn/10;
        %U_y(idx) = U_y(idx) + randn/10;
        
        if idx < N
            % Step 4: update velocities and steer angles
            U_x(idx+1) = U_x(idx) + U_x_dot * dT;
            U_y(idx+1) = U_y(idx) + U_y_dot * dT;
            r(idx+1) = r(idx) + r_dot * dT;
            s(idx+1) = s(idx) + s_dot * dT;
            e(idx+1) = e(idx) + e_dot * dT;
            dpsi(idx+1) = dpsi(idx) + dpsi_dot * dT;
            
        end
        
        % Compute lateral acceleration
        a_y(idx) = U_y_dot + r(idx) * U_x(idx);
        a_x(idx) = U_x_dot - r(idx) * U_y(idx);
        
    end
    
    % Setup return
    X = [U_x, U_y, r, e, s, dpsi, a_x, a_y, outUx_des];
    ctrl = delta;
    %animate(path, car, dpsi, s, e, delta);
end

function [Fdrag, Frr] = computeExternalForces(veh, state)
% Method for computing the external forces due to atmospheric drag and
% rolling resistance acting on the vehicle at a given speed.

    rho = 1.225;        % atmo density [km/m^3]
    
    % Computing atmospheric drag assuming constant density
    Fdrag = 0.5 * rho * veh.CdA * state.Ux_mps^2;    % drag force [N]
    
    Ux_rr_thresh_mps = 1;   % Setting threshold
    % Computing rolling resistance force, which only acts when the velocity
    % is greater than Ux_rr_thresh_mps
    Frr = (veh.f_rr * veh.W) * (state.Ux_mps > Ux_rr_thresh_mps); 
        
end