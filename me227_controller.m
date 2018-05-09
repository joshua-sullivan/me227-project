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
    
    % Create the state object for easier data passing
    state.s_m = s_m;
    state.e_m = e_m;
    state.deltaPsi_rad = deltaPsi_rad;
    state.Ux_mps = Ux_mps;
    state.Uy_mps = Uy_mps;
    state.r_radps = r_radps;
    
    % Create the vehicle object that holds relevant parameters
    Shelley = loadVehicleParams();   
    
    % Run a table look-up/interpolation from open-loop desired trajectory
    pathPlan = runPathPlanner(state, path);
       
    if modeSelector == 1
    % Runs the feedforward lookahead lateral controller with a
    % feedfoward/PI-feedback longitudinal controller.
        
        % Lateral control law
        delta_rad = runLookaheadFFWController(Shelley, state, pathPlan);
        
        % Compute the lateral tire forces
        [Fyf_N, ~] = computeFialaNLTireForce(Shelley, state, delta_rad);
        
        % Longitudinal control law
        [Fx_N, integrated_Ux_error] = runLongitudinalController(Shelley, state, delta_rad, ...
            Fyf_N, pathPlan, integrated_Ux_error, dt);
        
    elseif modeSelector == 2
    % Runs the lateral PID controller with a feedforward/PI-feedback
    % longitudinal controller.
        
        % Lateral control law
        delta_rad = []; %runLateralPIDController(Shelley, state, etc...);
        
        % Compute the lateral tire forces
        [Fyf_N, ~] = computeFialaNLTireForce(Shelley, state, delta_rad);
        
        % Longitudinal control law
        [Fx_N, integrated_Ux_error] = runLongitudinalController(Shelley, state, delta_rad, ...
            Fyf_N, pathPlan, integrated_Ux_error, dt);
        
    else
    % Maybe another...
        
    end

end

function veh = loadVehicleParams()
% Method for loading the vehicle parameters.
    
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
    veh.a = veh.wheelbase * veh.perc_W_r;       % distance from front wheel to CG [m]
    veh.b = veh.wheelbase * veh.perc_W_f;       % distance from rear wheel to CG [m]
    
    % Computing understeer gradient [rad/m/s^2]
    veh.K = veh.m * ((veh.perc_W_f / veh.Caf) - (veh.perc_W_r / veh.Car));
    
    % Computing the characteristic speed [m/s]
    veh.Vch = sqrt(veh.L / veh.K);
end

function delta_rad = runLookaheadFFWController(veh, state, pathPlan)
% Method for computing the steering angle using a lookahead controller with
% a feedforward term.

    % Set controller gains/params
    K_la = [];
    x_la = [];

    dPsi_ss = pathPlan.curv*((veh.m * veh.a * (state.Ux_mps^2) / ...    
        (veh.L * veh.Car_lin)) - veh.b);
               
    delta_ff = (K_la * x_la * dPsi_ss / veh.Caf) + ...
        (pathPlan.curv * (veh.L + veh.K * state.Ux_mps^2));
    
    delta_rad = (-K_la * (state.e_m + (x_la * state.deltaPsi_rad)) / veh.Caf_lin) + delta_ff;
end

function [alpha_f, alpha_r] = computeSideslipAngles(veh, state, delta_rad)
% Method for computing the side slip angles without making any small 
% angle approximations.
           
    alpha_f = atan2(state.Uy_mps + (veh.a * state.r_radps), state.Ux_mps) - delta_rad;
    alpha_r = atan2(state.Uy_mps - (veh.b * state.r_radps), state.Ux_mps);
            
end

function Fy = helperFialaTireForce(Ca, alpha, mu, mu_s, W)
% Helper function for computing Fy according to Fiala model
                
    Fy = (-Ca * tan(alpha)) + ...
         (((Ca^2) / (3 * mu * W) * (2 - (mu_s / mu))) * abs(tan(alpha)) *tan(alpha)) - ...
         (((Ca^3) / (9 * mu * mu * W * W)) * (tan(alpha)^3) * (1 - (2 * mu_s / 3 / mu)));
            
end

function [Fyf_N, Fyr_N] = computeFialaNLTireForce(veh, state, delta_rad)
% Method for computing nonlinear tire forces given by the Fiala tire model. 
        
    [alpha_f, alpha_r] = computeSideslipAngles(veh, state, delta_rad);

    alpha_slip_f = abs(atan2(3 * veh.mu_f * veh.Wf, veh.Caf));
    if abs(alpha_f) < alpha_slip_f
        Fyf_N = helperFialaTireForce(veh.Caf, alpha_f, veh.mu_f, veh.mu_fs, veh.Wf);
    else
        Fyf_N = -veh.mu_fs * veh.Wf * sign(alpha_f);
    end

    alpha_slip_r = abs(atan2(3 * veh.mu_r * veh.Wr, veh.Car));
    if abs(alpha_r) < alpha_slip_r
        Fyr_N = helperFialaTireForce(veh.Car, alpha_r, veh.mu_r, veh.mu_rs, veh.Wr);
    else
        Fyr_N = -veh.mu_rs * veh.Wr * sign(alpha_r); 
    end

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

function [Fx_N, integrated_error] = runLongitudinalController(veh, state, delta_rad, Fyf_N, pathPlan, integrated_error, dt)
% Simple feedforward/feedback longitudinal controller.  
%   FFW accounts for disturbance forces and dynamical coupling
%   FB accounts for setpoint velocity errors (uses a PI controller).
%   TODO: validate with simple straight path.

    % Set controller gains
    K_p = 1500;
    K_i = 10;
    
    % The feedback term drives the desired velocity difference to zero:
    error = (pathPlan.Ux_des_mps - state.Ux_mps);
    integrated_error = integrated_error + (error * dt);
    
    % Setting integrated error threshold on Ux_des (~ 1 mph or 0.5 m/s)
    % to avoid integrator windup.
    if abs(error) < 0.5
        integrated_error = 0;
    end
    
    Fx_fb_N = (K_p * error) + (K_i * integrated_error);

    % The feedforward term counteracts disturbance forces, accounts for the
    % desired longitudinal acceleration, and accounts for lateral dynamics
    % coupling
    [Fdrag, Frr] = computeExternalForces(veh, state);
    
    Fx_ffw_N = Fdrag + Frr + ...    % Accounting for disturbance forces
               (veh.m * pathPlan.Ux_dot_des_mps2) + ... % Desired longitudinal accel.
               (Fyf_N * sin(delta_rad)) - ...           % Dynamical coupling
               (veh.m * state.r_radps * state.Uy_mps);  % Dynamical coupling
           
    Fx_N = Fx_fb_N + Fx_ffw_N;
end

function pathPlan = runPathPlanner(state, path)
% Computes the open-loop path plan from the current state and the a priori
% known path data.

    pathPlan.curv = interp1(path.s_m, path.k_1pm, state.s_m);
    pathPlan.Ux_des_mps = [];
    pathPlan.Ux_dot_des_mps2 = [];
    
end