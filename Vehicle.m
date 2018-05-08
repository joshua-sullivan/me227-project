classdef Vehicle < handle
    % Vehicle class which makes use of the dynamic bicycle model
    
    properties (Constant)
        G = 9.81;   
        DEG2RAD = pi / 180;
        RAD2DEG = 180 / pi;
    end
    
    properties
        
        wheelbase
        mass
        Iz
        Caf
        Car
        usteer_grad_radpmps2
        usteer_grad_degpg
        mu_f
        mu_fs
        mu_r
        mu_rs
        front_W_percent
        Wf
        Wr
        length_a
        length_b
        Vch
        state
        state_dot
        ctrl
        rad_curv
        
    end
    
    methods
        function veh = Vehicle(varargin)
        % Vehicle constructor
            
            % Dynamical system structs set to default
            veh.state               = struct('s', 0., 'e', 0., 'dPsi', 0.,...
                                             'x_dot', 0., 'y_dot', 0., 'z_dot', 0.,...
                                             'roll', 0., 'pitch', 0., 'yaw', 0.);            
            veh.state_dot           = struct('s_dot', [], 'e_dot', [], 'dPsi_dot', [],...
                                             'x_ddot', [], 'y_ddot', [], 'z_ddot', [], ...
                                             'roll_dot', [], 'pitch_dot', [], 'yaw_dot', []);
            veh.ctrl                = struct('lat', [], 'lon', [], 'config', []);                     
            
            veh.ctrl.lon            = struct('inputs', [], 'params', [], 'mode', 'cruise_control');
            veh.ctrl.lon.inputs     = struct('Fx_f', [], 'Fx_r', []);
            veh.ctrl.lon.params     = struct('K_drive', 0., 'user_input', 0., 'x_dot_des', 0.);
            
            veh.ctrl.lat            = struct('inputs', [], 'params', [], 'mode', 'lookahead');
            veh.ctrl.lat.inputs     = struct('steer_angle', []);
            veh.ctrl.lat.params     = struct('K_la', 0., 'x_la', 0., 'Caf_la', 0., 'Car_laff', 0., 'user_input', 0.);

            
            veh.rad_curv            = [];
            
            % Vehicle Properties set to defaults
            veh.wheelbase         = 0;
            veh.mass              = 0;
            veh.Iz                = 0;
            veh.Caf               = 0;
            veh.Car               = 0;
            veh.usteer_grad_radpmps2   = 0;
            veh.usteer_grad_degpg = 0;
            veh.mu_f              = 0;
            veh.mu_fs             = 0;
            veh.mu_r              = 0;
            veh.mu_rs             = 0;
            veh.front_W_percent   = 0;
            veh.Wf                = 0;
            veh.Wr                = 0;
            veh.length_a          = 0;
            veh.length_b          = 0;
            veh.Vch               = 0;
          
            % Processing constructor inputs to attributes
            if nargin ~= 0
                % Parse Inputs
                for ii = 1 : 2 : nargin
                    if any(strcmp(properties(veh), varargin{ii}))
                        veh.(varargin{ii}) = varargin{ii+1};
                    else
                        warning(['Ignoring unrecognized object property: ' varargin{ii}]);
                    end
                end
                
                % Computing implicit parameters from constructor inputs
                veh.Wf = veh.front_W_percent * veh.mass * veh.G;
                veh.Wr = (1 - veh.front_W_percent) * veh.mass * veh.G;
                veh.length_a = veh.wheelbase * (1 - veh.front_W_percent);
                veh.length_b = veh.wheelbase * veh.front_W_percent;
                [veh.usteer_grad_radpmps2, veh.usteer_grad_degpg] = veh.calculateUndersteerGrad();
                veh.Vch = veh.calcCharSpeed();
                
            end
            
            % Setting default controller cornering stiffness to vehicles
            veh.ctrl.lat.params.Caf_la = veh.Caf;
            veh.ctrl.lat.params.Car_laff = veh.Car;
        end
        
        function setVehStateIC(veh, varargin)
            % Method for programmatically setting the internal vehicle
            % dynamical state initial conditions.
            
            num_args = nargin - 1;  % Subtracting one to remove 'veh' from args number
            
            if nargin ~= 0            
                % Parse Inputs
                for ii = 1 : 2 : num_args
                    if any(strcmp(fieldnames(veh.state), varargin{ii}))
                        veh.state.(varargin{ii})(:, end) = varargin{ii+1};
                    else
                        warning(['Ignoring unrecognized state element: ' varargin{ii}]);
                    end
                end
                
            end           
        end
        
       
        function [usteer_grad_radpmps2, usteer_grad_degpg] = calculateUndersteerGrad(veh)
            % Method for computing the vehicle understeer gradient
            
            usteer_grad_radpmps2 = ((veh.Wf / veh.Caf) - (veh.Wr / veh.Car))*(1/veh.G);
            usteer_grad_degpg = usteer_grad_radpmps2 * veh.RAD2DEG;
            
        end
        
        function Vch = calcCharSpeed(veh)
            % Method for computing the vehicle characteristic speed
            
            if veh.usteer_grad_radpmps2 == 0 
                Vch = 0;
            elseif veh.usteer_grad_radpmps2 < 0
                Vch = sqrt(-veh.wheelbase/veh.usteer_grad_radpmps2);
            else
                Vch = sqrt(veh.wheelbase/veh.usteer_grad_radpmps2);
            end
            
        end
        
        function [nat_freq, damp_ratio] = calc2ndOrderSpec(veh, Ux)
            % Method for computing the vehicle response characterstics at a
            % given constant longitudinal velocity, Ux.  The response is
            % characterized by the natural frequency and damping ratio.  
            % Outputs vectors that are the same size as the input
            % vector of longitudinal velocities (i.e., for table-making).
            
            nat_freq = zeros(length(Ux), 1);
            damp_ratio = zeros(length(Ux), 1);
            
            for ii = 1:length(Ux)
                              
                omega_sq = ((veh.length_b*veh.Car - veh.length_a*veh.Caf)/veh.Iz) * ...
                    (1 + ((veh.Vch/Ux)^2));
                omega = sqrt(omega_sq);

                zeta0_num = (veh.Iz*(veh.Caf + veh.Car)) + ...
                    (veh.mass*(((veh.length_a^2)*veh.Caf) + ((veh.length_b^2)*veh.Car)));
                zeta0_den = 2*veh.wheelbase*sqrt(veh.mass*veh.Iz*veh.Caf*veh.Car);
                zeta0 = zeta0_num/zeta0_den;

                zeta = zeta0 / sqrt(1 + ((Ux/veh.Vch)^2));

                nat_freq(ii, 1) = omega;
                damp_ratio(ii, 1) = zeta;
                
            end
                     
        end
               
        function Krss = calcSSYawGain(veh, Ux)
            % Method for computing the steady state yaw gain for a given
            % constant longitudinal velocity.  Outputs a vector of gains
            % that is the same size as the input vector of longitudinal
            % velocities (i.e., for table-making).
            
            Krss = zeros(length(Ux), 1);
            
            for ii = 1:length(Ux)
                    Krss(ii, 1) = Ux(ii)/(veh.wheelbase + (veh.usteer_grad_radpmps2*(Ux(ii)^2)));
            end
        end
        
        
        function [alpha_f, alpha_r] = calcSideslipAngles(veh)
            % Method for computing the side slip angles without making any
            % small angle approximations.
            
            x_dot = veh.state.x_dot(:, end);
            y_dot = veh.state.y_dot(:, end);
            yaw = veh.state.yaw(:, end);
            steer_angle = veh.ctrl.lat.inputs.steer_angle(:, end);
            
            alpha_f = atan2(y_dot + (veh.length_a*yaw), x_dot) - steer_angle;
            alpha_r = atan2(y_dot - (veh.length_b*yaw), x_dot);
            
        end
        
        
        function [Fy_f, Fy_r] = calcLinearTireForce(veh)
            % Method for computing the linear tire force
            
            [alpha_f, alpha_r] = veh.calcSideslipAngles();
            
            Fy_f = -veh.Caf * alpha_f;
            Fy_r = -veh.Car * alpha_r;          
            
        end
        
        
        
        function [Fy_f, Fy_r] = calcFialaNLTireForce(veh)
            % Method for computing nonlinear tire forces given by the Fiala
            % tire model. 
            
            function Fy = helperFialaTireForce(Ca, alpha, mu, mu_s, W)
                % Helper function for computing Fy according to Fiala model
                
                Fy = (-Ca*tan(alpha)) + ...
                     (((Ca^2)/(3*mu*W)*(2-(mu_s/mu)))*abs(tan(alpha))*tan(alpha)) - ...
                     (((Ca^3)/(9*mu*mu*W*W))*(tan(alpha)^3)*(1-(2*mu_s/3/mu)));
            
            end
            
            [alpha_f, alpha_r] = veh.calcSideslipAngles();
            
            alpha_slip_f = abs(atan2(3*veh.mu_f*veh.Wf, veh.Caf));
            
            if abs(alpha_f) < alpha_slip_f
                
                Fy_f = helperFialaTireForce(veh.Caf, alpha_f, veh.mu_f, veh.mu_fs, veh.Wf);
                
            else
                
                Fy_f = -veh.mu_fs*veh.Wf*sign(alpha_f);
                
            end
            
            
            alpha_slip_r = abs(atan2(3*veh.mu_r*veh.Wr, veh.Car));
            
            if abs(alpha_r) < alpha_slip_r
                
                Fy_r = helperFialaTireForce(veh.Car, alpha_r, veh.mu_r, veh.mu_rs, veh.Wr);
                
            else
                
                Fy_r = -veh.mu_rs*veh.Wr*sign(alpha_r); 
                
            end
            
        end
        
        
        function eulerIntegrate(veh, dt)
            % Method for running one step of the Euler integration method
            % on the vehicle dynamical state elements.
            
            veh.state.x_dot(:, end+1) = veh.state.x_dot(:, end) + (veh.state_dot.x_ddot(:, end) * dt);
            veh.state.y_dot(:, end+1) = veh.state.y_dot(:, end) + (veh.state_dot.y_ddot(:, end) * dt);
            veh.state.yaw(:, end+1) = veh.state.yaw(:, end) + (veh.state_dot.yaw_dot(:, end) * dt);

            veh.state.s(:, end+1) = veh.state.s(:, end) + (veh.state_dot.s_dot(:, end) * dt);
            veh.state.e(:, end+1) = veh.state.e(:, end) + (veh.state_dot.e_dot(:, end) * dt);
            veh.state.dPsi(:, end+1) = veh.state.dPsi(:, end) + (veh.state_dot.dPsi_dot(:, end) * dt);
        end
        
        function updatePlanarBicycleDynamics(veh, use_NL_tires)
            % Method for computing the state derivatives according to the
            % planar dynamic bicycle model.  Assumes a constant
            % longitudinal velocity.
            
            x_dot = veh.state.x_dot(:, end);
            yaw = veh.state.yaw(:, end);
            
            if use_NL_tires
                
                [Fy_f, Fy_r] = veh.calcFialaNLTireForce();
                
            else
                
                [Fy_f, Fy_r] = veh.calcLinearTireForce();
                
            end
            
            x_ddot = 0;     % Assuming constant longitudinal velocity!
            y_ddot = ((Fy_f + Fy_r)/veh.mass) - (yaw*x_dot);
            yaw_dot = (veh.length_a*Fy_f/veh.Iz) - (veh.length_b*Fy_r/veh.Iz);
            
            veh.state_dot.x_ddot(:, end+1) = x_ddot;
            veh.state_dot.y_ddot(:, end+1) = y_ddot;
            veh.state_dot.yaw_dot(:, end+1) = yaw_dot;
            
        end
        
        function A = lookaheadDynamics(veh, Ux)
            % Method for computing the look ahead dynamics matrix given a
            % longitudinal velocity, Ux.
            
            m = veh.mass;
            a = veh.length_a;
            b = veh.length_b;
            
            K_la = veh.ctrl.lat.params.K_la;
            x_la = veh.ctrl.lat.params.x_la;
            
            A = [0, 1, 0, 0;
                -K_la/m, -(veh.Caf + veh.Car)/(m*Ux), ((veh.Caf+veh.Car)/m)-(K_la*x_la/m), ((-a*veh.Caf)+(b*veh.Car))/(m*Ux);
                 0, 0, 0, 1; 
                -K_la*a/veh.Iz, ((b*veh.Car)-(a*veh.Caf))/(veh.Iz*Ux), (((a*veh.Caf)-(b*veh.Car))/veh.Iz)-(K_la*a*x_la/veh.Iz), -(((a^2)*veh.Caf)+((b^2)*veh.Car))/(veh.Iz*Ux)]; 
            
        end
        
        function lookaheadController(veh)
            % Method for  computing the steer angle for the lookahead 
            % controller
            
            x_la = veh.ctrl.lat.params.x_la;
            K_la = veh.ctrl.lat.params.K_la;
            veh.ctrl.lat.inputs.steer_angle(:, end+1) = -K_la *...
                (veh.state.e(:, end) + (x_la * veh.state.dPsi(:, end))) / veh.ctrl.lat.params.Caf_la;
                   
        end  
        
        function lookaheadFFController(veh)
            % Method for computing the steer angle for the lookahead
            % controller with feedforward term.
        
            x_la = veh.ctrl.lat.params.x_la;
            K_la = veh.ctrl.lat.params.K_la;
            
            K_rc = veh.rad_curv(:, end);
            m = veh.mass;
            a = veh.length_a;
            b = veh.length_b;
            L = veh.wheelbase;
            Ux = veh.state.x_dot(:, end);
            
            dPsi_ss = K_rc*((m*a*(Ux^2)/(L*veh.ctrl.lat.params.Car_laff)) - b);
            delta_ff = ((K_la * x_la * dPsi_ss)/veh.Caf) + K_rc * (L + veh.usteer_grad_radpmps2 * Ux^2);
%             delta_ff = ((K_la * x_la * dPsi_ss)/veh.ctrl.lat.params.Caf_la) + K_rc * (L + veh.usteer_grad_radpmps2 * Ux^2);

            veh.lookaheadController();
            veh.ctrl.lat.inputs.steer_angle(:, end) = veh.ctrl.lat.inputs.steer_angle(:, end) + delta_ff;
                
        end
        
                
        function longitudinalCruiseController(veh)
            % Method for computing the simple longitudinal tire forces for
            % a given velocity setpoint. 
            
            Fx = veh.ctrl.lon.params.K_drive*(veh.ctrl.lon.params.x_dot_des - veh.state.x_dot(:, end)); 
            
            Fx_f = 0.5*Fx;
            Fx_r = 0.5*Fx;
            
            veh.ctrl.lon.inputs.Fx_f(:, end+1) = Fx_f;
            veh.ctrl.lon.inputs.Fx_r(:, end+1) = Fx_r;
            
        end
            
        
        function updateControl(veh)
            % Method for updating the controller as the dynamical state
            % evolves during the simulation.
            
            switch veh.ctrl.lat.mode
                % Lateral controller modes
                
                case 'user_cmd'
                     veh.ctrl.inputs.steer_angle(:, end+1) = veh.ctrl.lat.params.user_cmd;
                     
                case 'lookahead'
                    veh.lookaheadController();
                    
                case 'lookahead_withFF'
                    veh.lookaheadFFController();
            end
            
            
            switch veh.ctrl.lon.mode
                % Longitudinal controller modes
                
                case 'cruise_control'
                    veh.longitudinalCruiseController()
            end
            
        end
        

        function propagatePlanarBicycleDynamics(veh, sim_cfg)
            % Method for running a simulation of the planar dynamic bicycle
            % model.  
            
            time = (0:sim_cfg.dt:sim_cfg.tf)';
            dim_t = length(time);
            
            for ii = 1:dim_t-1
                
                veh.updateControl();
                veh.updatePlanarBicycleDynamics(sim_cfg.use_NL_tires);
                veh.eulerIntegrate(dt);
                
            end
                  
        end
        
        
        function calcRadiusOfCurvature(veh, sim_cfg)
            % Method for computing the radius of curvature from path.
            
            veh.rad_curv(:, end+1) = interp1(sim_cfg.path.s, sim_cfg.path.k, veh.state.s(:, end));
            
        end
        
 
        function updateFullNLDynamics(veh)
            % Method that computes the complete state derivative set for
            % the full nonlinear dynamical vehicle model. 
            
            a = veh.length_a;
            b = veh.length_b;
            
            Fx_f = veh.ctrl.lon.inputs.Fx_f(:, end);
            Fx_r = veh.ctrl.lon.inputs.Fx_r(:, end);
            [Fy_f, Fy_r] = veh.calcFialaNLTireForce();
            
            x_dot = veh.state.x_dot(:, end);
            y_dot = veh.state.y_dot(:, end);
            yaw = veh.state.yaw(:, end);
            
            e = veh.state.e(:, end);
            dPsi = veh.state.dPsi(:, end);
            
            delta = veh.ctrl.lat.inputs.steer_angle(:, end);
            
            veh.state_dot.x_ddot(:, end+1) = ((Fx_r + (Fx_f*cos(delta)) - (Fy_f*sin(delta)))/veh.mass) + (yaw*y_dot);
            veh.state_dot.y_ddot(:, end+1) = (((Fy_f*cos(delta)) + Fy_r + (Fx_f*sin(delta)))/veh.mass) - (yaw*x_dot);
            veh.state_dot.yaw_dot(:, end+1) = ((a*Fy_f*cos(delta)) + (a*Fx_f*sin(delta)) - (b*Fy_r))/veh.Iz;
            
            veh.state_dot.s_dot(:, end+1) = (1/(1- (e*veh.rad_curv(:, end))))*((x_dot*cos(dPsi)) - (y_dot*sin(dPsi)));
            veh.state_dot.e_dot(:, end+1) = (y_dot*cos(dPsi)) + (x_dot*sin(dPsi));
            veh.state_dot.dPsi_dot(:, end+1) = yaw - (veh.rad_curv(:, end)*veh.state_dot.s_dot(:, end));
            
                      
        end
        
        function propagateFullNLDynamics(veh, sim_cfg)
            % Method for propagating the full nonlinear dynamical model. 
            
            time = (0:sim_cfg.dt:sim_cfg.tf)';
            dim_t = length(time);
            
            for ii = 1:dim_t-1
                
                veh.calcRadiusOfCurvature(sim_cfg);
                veh.updateControl();
                veh.updateFullNLDynamics();
                veh.eulerIntegrate(sim_cfg.dt);
                
            end
                  
        end
        
        function plotNLSimulationResults(veh, sim_cfg)
            % Method for plotting results.  
            
            FS = 18;
            rads2deg = 180/pi;
            
            time = (0:sim_cfg.dt:sim_cfg.tf)';
            
            figure
            subplot(2, 3, 1)
            plot(time, veh.state.x_dot,'LineWidth', 2)
            grid on
            xlabel('Time (sec)', 'FontSize', FS, 'Interpreter','Latex')
            ylabel('Longitudinal Velocity, $$U_x$$ (m/s)', 'FontSize', FS, 'Interpreter','Latex')
            set(gca,'FontSize', FS,'TickLabelInterpreter', 'Latex')

            subplot(2, 3, 2)
            plot(time, veh.state.y_dot,'LineWidth', 2)
            grid on
            xlabel('Time (sec)', 'FontSize', FS, 'Interpreter','Latex')
            ylabel('Lateral Velocity, $$U_y$$ (m/s)', 'FontSize', FS, 'Interpreter','Latex')
            set(gca,'FontSize', FS,'TickLabelInterpreter', 'Latex')

            subplot(2, 3, 3)
            plot(time, veh.state.yaw*rads2deg,'LineWidth', 2)
            grid on
            xlabel('Time (sec)', 'FontSize', FS, 'Interpreter','Latex')
            ylabel('Yaw Rate, r (deg/s)', 'FontSize', FS, 'Interpreter','Latex')
            set(gca,'FontSize', FS,'TickLabelInterpreter', 'Latex')

            subplot(2, 3, 4)
            plot(time, veh.state.s,'LineWidth', 2)
            grid on
            xlabel('Time (sec)', 'FontSize', FS, 'Interpreter','Latex')
            ylabel('Path Distance, s (m)', 'FontSize', FS, 'Interpreter','Latex')
            set(gca,'FontSize', FS,'TickLabelInterpreter', 'Latex')

            subplot(2, 3, 5)
            plot(time, veh.state.e,'LineWidth', 2)
            grid on
            xlabel('Time (sec)', 'FontSize', FS, 'Interpreter','Latex')
            ylabel('Lateral Path Error, e (m)', 'FontSize', FS, 'Interpreter','Latex')
            set(gca,'FontSize', FS,'TickLabelInterpreter', 'Latex')

            subplot(2, 3, 6)
            plot(time, veh.state.dPsi*rads2deg,'LineWidth', 2)
            grid on
            xlabel('Time (sec)', 'FontSize', FS, 'Interpreter','Latex')
            ylabel('Heading Error, $$\Delta\psi$$ (deg)', 'FontSize', FS, 'Interpreter','Latex')
            set(gca,'FontSize', FS,'TickLabelInterpreter', 'Latex')
            
        end
        
              
    end
    
end

