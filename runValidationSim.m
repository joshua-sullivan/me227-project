% runValidationSim  Validates the vehicle lateral/longitudinal controllers.
% 
%   sim_data = runValidationSim() runs a numerical simulation of the 
%   vehicle dynamics and control.  Both lateral and longitudinal control is
%   implemented in this project.
% -------------------------------------------------------------------------
% 
%   Outputs:
%       sim_data        - (struct) simulation data
% ------------------------------------------------------------------------- 
%   Authors:
%       Rohan Kulkarni
%       Vince Salpietro
%       Sumant Sharma
%       Josh Sullivan
% -------------------------------------------------------------------------
function data = runValidationSim()

    % Set simulation parameters here
    t_start = 0;                    % Simulation start time [sec]
    t_end = 60;                     % Simulation end time [sec]
    dt = 0.01;                      % Simulation time step [sec]
    t = (t_start : dt : t_end)';    % Simulation time vecotr [sec]
    dim_t = length(t);              % Number of simulation steps [-]
    
    % Set vehicle parameters
    Shelley = loadVehicleParams();
    
    % Loading in data files containing path, edge_list, position data, and
    % *open-loop computed desired path from path planner*
    proj_data_filename = 'project_data.mat';    % placeholder file
    proj_data = load(proj_data_filename);
    path = proj_data.path;
%     edge_list = proj_data.edge_list;
    
    % Set controller mode here
    modeSelector = 1;
    
    % Initialize state variables
    s_m = zeros(dim_t,1);
    e_m = zeros(dim_t,1);
    deltaPsi_rad = zeros(dim_t,1);
    Ux_mps = zeros(dim_t,1);
    Uy_mps = zeros(dim_t,1);
    r_radps = zeros(dim_t,1);
        

    % Run main simulation stepwise
    for ii = 1:dim_t
        
        state.s_m = s_m(ii);
        state.e_m = e_m(ii);
        state.deltaPsi_rad = deltaPsi_rad(ii);
        state.Ux_mps = Ux_mps(ii);
        state.Uy_mps = Uy_mps(ii);
        state.r_radps = r_radps(ii);
        
        % Run a table look-up and interpolation on desired open-loop path
        Ux_des_mps(ii) = interp1(path.s_m, path.Ux_des_mps, s_m(ii));
        Ux_dot_des_mps2(ii,1) = interp1(path.s_m, path.Ux_dot_des_mps2, s_m(ii));
        
        % Computing the path curvature from interpolated table look-up
        curv(ii) = interp1(path.s_m, path.k_1pm, s_m(ii));
        
        % Running the control laws
        [delta_rad(ii), Fx_N(ii)] = ...
            me227_controller(s_m(ii), e_m(ii), deltaPsi_rad(ii), ...
                Ux_mps(ii), Uy_mps(ii), r_radps(ii), modeSelector, path);
            
        % Computing the lateral tire forces
        [Fyf_N(ii), Fyr_N(ii)] = computeFialaNLTireForce(Shelley, state, delta_rad(ii));
        
        % Computing external forces
        [Fdrag_N, Frr_N] = computeExternalForces(Shelley, state);
        Fgrade_N = 0.05 * Shelley.W * (0.5 * rand(1));
        Fext_N = Fdrag_N + Frr_N + Fgrade_N;
        
        % Computing state derivatives
        Ux_dot_mps2(ii) = (1/Shelley.m) * (Fx_N(ii) - Fext_N - (Fyf_N(ii) * sin(delta_rad(ii)))) + ...
            (r_radps(ii) * Uy_mps(ii));
        Uy_dot_mps2(ii) = (((Fyf_N(ii) * cos(delta_rad(ii))) + Fyr_N(ii))/Shelley.m) - (r_radps(ii) * Ux_mps(ii));
        r_dot_radps2(ii) = (Shelley.a * Fyf_N(ii) * cos(delta_rad(ii)) - Shelley.b * Fyr_N(ii)) / Shelley.Iz;
        s_dot_mps(ii) = (1 / (1 - e_m(ii) * Shelley.K)) * (Ux_mps(ii) * cos(deltaPsi_rad(ii)) - Uy_mps(ii) * sin(deltaPsi_rad(ii)));
        e_dot_mps(ii) = Uy_mps(ii) * cos(deltaPsi_rad(ii)) + Ux_mps(ii) * sin(deltaPsi_rad(ii));
        deltaPsi_dot_radps2(ii) = r_radps(ii) - curv(ii)*s_dot_mps(ii);
        
        % Computing the longitudinal and lateral accelerations
        ax_mps2(ii) = Ux_dot_mps2(ii) - (r_radps(ii) * Uy_mps(ii));
        ay_mps2(ii) = Uy_dot_mps2(ii) + (r_radps(ii) * Ux_mps(ii));
        
        % Running a numerical integration step
        if ii < dim_t
            
            Ux_mps(ii + 1) = Ux_mps(ii) + (Ux_dot_mps2(ii) * dt);
            Uy_mps(ii + 1) = Uy_mps(ii) + (Uy_dot_mps2(ii) * dt);
            r_radps(ii + 1) = r_radps(ii) + (r_dot_radps2(ii) * dt);
            s_m(ii + 1) = s_m(ii) + (s_dot_mps(ii) * dt);
            e_m(ii + 1) = e_m(ii) + (e_dot_mps(ii) * dt);
            deltaPsi_rad(ii + 1) = deltaPsi_rad(ii) + (deltaPsi_dot_radps(ii) * dt);
            
        end
        
    end

    % Creating data field names and data struct
    data_fields = {'Ux_des_mps', 'Ux_dot_des_mps', 'curv', 'delta_rad', 'Fx_N', ...
                   'Fyf_N', 'Fyr_N', 'Ux_dot_mps2', 'Uy_dot_mps2', 'r_dot_radps2', ...
                   's_dot_mps', 'e_dot_mps', 'deltaPsi_dot_radps2', 'ax_mps2', 'ay_mps2', ...
                   'Ux_mps', 'Uy_mps', 'r_radps', 's_m', 'e_m', 'deltaPsi_rad'};
    data = struct();
    
    % Populating data struct with variables from above
    for ii = 1:length(data_fields)        
        data.(data_fields{ii}) = eval('data_fields{ii}');
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
    veh.a = veh.L * veh.perc_W_r;       % distance from front wheel to CG [m]
    veh.b = veh.L * veh.perc_W_f;       % distance from rear wheel to CG [m]
    
    % Computing understeer gradient [rad/m/s^2]
    veh.K = veh.m * ((veh.perc_W_f / veh.Caf) - (veh.perc_W_r / veh.Car));
    
    % Computing the characteristic speed [m/s]
    veh.Vch = sqrt(veh.L / veh.K);
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