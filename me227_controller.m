% me227_controller  Employs a lateral and longitudinal vehicle controller.
% 
%   [delta_rad, Fx_N] = me227_controller(s_m, e_m, deltaPsi_rad, Ux_mps,
%   Uy_mps, r_radps, modeSelector, path) produces the steering angle and
%   longitudinal drive force commands to the Audi TT Shelley.  
% -------------------------------------------------------------------------
%   Inputs:
%       s_m             - (1, 1) path coordinate [m]  
%       e_m             - (1, 1) lateral error [m] 
%       deltaPsi_rad    - (1, 1) heading error [rad]
%       Ux_mps          - (1, 1) longitudinal velocity [m/s]
%       Uy_mps          - (1, 1) lateral velocity [m/s]
%       r_radps         - (1, 1) yaw rate [rad/s]
%       modeSelector    - (1, 1) controller mode flag
%       path            - (struct) path definition
% 
%   Outputs:
%       delta_rad       - (1, 1) steering angle command [rad]
%       Fx_N            - (1, 1) longitudinal force command [N]
% ------------------------------------------------------------------------- 
%   Authors:
%       Rohan Kulkarni
%       Vince Salpietro
%       Sumant Sharma
%       Josh Sullivan
% -------------------------------------------------------------------------

function [delta_rad, Fx_N] = me227_controller(s_m, e_m, deltaPsi_rad, Ux_mps, Uy_mps, r_radps, modeSelector, path)

    %#codegen % this directive is needed to run your code on the car
    
    % Declare persistent variables here and use a time update to
    % approximate a numerical integration (useful for integral terms in the
    % PI(D) controller)
    persistent integrated_Ux_error
    persistent prev_e_m
    
    if isempty(integrated_Ux_error)
        integrated_Ux_error = 0;
    end
    
    % Create the state object for easier data passing
    state.s_m = s_m;
    state.e_m = e_m;
    state.deltaPsi_rad = deltaPsi_rad;
    state.Ux_mps = Ux_mps;
    state.Uy_mps = Uy_mps;
    state.r_radps = r_radps;
    
    %% Create the vehicle object that holds relevant parameters
    % Setting parameters from the provided data sheet
    veh.m = 1659;           % mass [kg]
    veh.L = 2.468;          % wheelbase [m]
    veh.Iz = 2447;          % out-of-plane MOI [kg m^2]
    veh.perc_W_f = 0.577;   % percent weight on front [fractional]
    veh.f_rr = 0.0157;      % coefficient of rolling resistance
    veh.CdA = 0.594;        % drag coeff. * drag area [m^2]
    veh.Caf_lin = 188e03;   % linear front tire stiffness [N/rad]
    veh.Car_lin = 203e03;   % linear rear tire stiffness [N/rad]
    veh.mu_f = 0.97;        % cofficient of friction (front)
    veh.mu_fs = 0.97;       % cofficient of friction (front)
    veh.mu_r = 1.03;        % cofficient of friction (rear)
    veh.mu_rs = 1.03;       % cofficient of friction (rear)
    veh.Caf = 275e03;       % nonlinear front tire stiffness [N/rad]
    veh.Car = 265e03;       % nonlinear rear tire stiffness [N/rad]
    % Computing dependent parameters
    veh.W = veh.m * 9.81;                       % vehicle weight (m*g) [N]
    veh.Wf = veh.perc_W_f * veh.W;              % vehicle front weight [N]
    veh.Wr = veh.W - veh.Wf;                    % vehicle rear weight [N]
    veh.perc_W_r = 1 - veh.perc_W_f;            % percentage weight on rear [fractional]
    veh.a = veh.L * veh.perc_W_r;       % distance from front wheel to CG [m]
    veh.b = veh.L * veh.perc_W_f;       % distance from rear wheel to CG [m]
    % Computing understeer gradient [rad/m/s^2]
    veh.K = veh.m * ((veh.perc_W_f / veh.Caf) - (veh.perc_W_r / veh.Car));
    % Computing the characteristic speed [m/s]
    veh.Vch = sqrt(veh.L / veh.K);
    
    %% Run a table look-up/interpolation from open-loop desired trajectory
    pathPlan.curv = interp1(path.s_m, path.k_1pm, state.s_m);
    pathPlan.Ux_des_mps = interp1(path.s_m, path.Ux_des_mps, state.s_m);
    pathPlan.Ux_dot_des_mps2 = interp1(path.s_m, path.Ux_dot_des_mps2, state.s_m);
    
    % initialize dt
    dt = 0.01;
    
    %% Controllers
    if modeSelector == 1
    % Runs the feedforward lookahead lateral controller with a
    % feedfoward/PI-feedback longitudinal controller.
        
        % Lateral control law
        K_la = 3500;
        x_la = 15;
        dPsi_ss = pathPlan.curv*((veh.m * veh.a * (state.Ux_mps^2) / ...    
            (veh.L * veh.Car_lin)) - veh.b);
        delta_ff = (K_la * x_la * dPsi_ss / veh.Caf_lin) + ...
            (pathPlan.curv * (veh.L + veh.K * state.Ux_mps^2));
        delta_rad = (-K_la * (state.e_m + (x_la * state.deltaPsi_rad)) / veh.Caf_lin) + delta_ff;
        
        % Compute the lateral tire forces
        alpha_f = atan2(state.Uy_mps + (veh.a * state.r_radps), state.Ux_mps) - delta_rad;
        alpha_slip_f = abs(atan2(3 * veh.mu_f * veh.Wf, veh.Caf));
        if abs(alpha_f) < alpha_slip_f
            Fyf_N = helperFialaTireForce(veh.Caf, alpha_f, veh.mu_f, veh.mu_fs, veh.Wf);
        else
            Fyf_N = -veh.mu_fs * veh.Wf * sign(alpha_f);
        end
        
        % Longitudinal control law
        K_p = 1500;
        K_i = 10;
        error = (pathPlan.Ux_des_mps - state.Ux_mps);
        integrated_Ux_error = integrated_Ux_error + (error * dt);
        if abs(error) < 0.5 
            integrated_Ux_error = 0;    % Setting integrated error threshold on Ux_des (~ 1 mph or 0.5 m/s) to avoid integrator windup.
        end
        Fx_fb_N = (K_p * error) + (K_i * integrated_Ux_error);
        rho = 1.225;                                        % atmo density [km/m^3]
        Fdrag = 0.5 * rho * veh.CdA * state.Ux_mps^2;       % drag force [N]
        Ux_rr_thresh_mps = 1;                               % Setting threshold
        Frr = (veh.f_rr * veh.W) * (state.Ux_mps > Ux_rr_thresh_mps);   % Computing rolling resistance force, which only acts when the velocity is greater than Ux_rr_thresh_mps
        Fx_ffw_N = Fdrag + Frr + ...    % Accounting for disturbance forces
                   (veh.m * pathPlan.Ux_dot_des_mps2) + ... % Desired longitudinal accel.
                   (Fyf_N * sin(delta_rad)) - ...           % Dynamical coupling
                   (veh.m * state.r_radps * state.Uy_mps);  % Dynamical coupling
        Fx_N = Fx_fb_N + Fx_ffw_N;
        
    elseif modeSelector == 2
    % Runs the lateral PID controller with a feedforward/PI-feedback
    % longitudinal controller.
        
        % Lateral control law (compute delta based on PD control law)
        if  isempty(prev_e_m)
            prev_e_m= 0;
        end
        Kp = 6000/veh.Caf;
        Kd = 45;
        delta_rad = -Kp * (state.e_m + Kd * prev_e_m);
        prev_e_m = state.e_m;
        
        % Compute the lateral tire forces
        alpha_f = atan2(state.Uy_mps + (veh.a * state.r_radps), state.Ux_mps) - delta_rad;
        alpha_slip_f = abs(atan2(3 * veh.mu_f * veh.Wf, veh.Caf));
        if abs(alpha_f) < alpha_slip_f
            Fyf_N = helperFialaTireForce(veh.Caf, alpha_f, veh.mu_f, veh.mu_fs, veh.Wf);
        else
            Fyf_N = -veh.mu_fs * veh.Wf * sign(alpha_f);
        end
        
        % Longitudinal control law
        K_p = 1500;
        K_i = 10;
        error = (pathPlan.Ux_des_mps - state.Ux_mps);
        integrated_Ux_error = integrated_Ux_error + (error * dt);
        if abs(error) < 0.5 
            integrated_Ux_error = 0;    % Setting integrated error threshold on Ux_des (~ 1 mph or 0.5 m/s) to avoid integrator windup.
        end
        Fx_fb_N = (K_p * error) + (K_i * integrated_Ux_error);
        rho = 1.225;                                        % atmo density [km/m^3]
        Fdrag = 0.5 * rho * veh.CdA * state.Ux_mps^2;       % drag force [N]
        Ux_rr_thresh_mps = 1;                               % Setting threshold
        Frr = (veh.f_rr * veh.W) * (state.Ux_mps > Ux_rr_thresh_mps);   % Computing rolling resistance force, which only acts when the velocity is greater than Ux_rr_thresh_mps
        Fx_ffw_N = Fdrag + Frr + ...    % Accounting for disturbance forces
                   (veh.m * pathPlan.Ux_dot_des_mps2) + ... % Desired longitudinal accel.
                   (Fyf_N * sin(delta_rad)) - ...           % Dynamical coupling
                   (veh.m * state.r_radps * state.Uy_mps);  % Dynamical coupling
        Fx_N = Fx_fb_N + Fx_ffw_N;
        
    end

end

function Fy = helperFialaTireForce(Ca, alpha, mu, mu_s, W)
% Helper function for computing Fy according to Fiala model
                
    Fy = (-Ca * tan(alpha)) + ...
         (((Ca^2) / (3 * mu * W) * (2 - (mu_s / mu))) * abs(tan(alpha)) *tan(alpha)) - ...
         (((Ca^3) / (9 * mu * mu * W * W)) * (tan(alpha)^3) * (1 - (2 * mu_s / 3 / mu)));
            
end


function delta_rad = runLateralPIDController(veh, state, pathPlan)

end
