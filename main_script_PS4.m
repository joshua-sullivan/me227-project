clear all
close all
clc

deg2rads = pi / 180;
rads2deg = 1 / deg2rads;

%% Problem 1
Shelley = Vehicle('Caf', 188e03, 'Car', 203e03, 'mass', 1648, 'Iz', 2235, ...
                  'wheelbase', 2.468, 'front_W_percent', 0.577);
              
% Problem 1.1
x_la_lookdown = 0;
K_la_lookdown = 3*deg2rads * Shelley.Caf;
 
% Problem 1.2
Shelley.ctrl.lat.params.K_la = K_la_lookdown;

Ux = 5:5:30;
B = zeros(4, 1);
C = eye(4);
D = 0;

FS = 18;
figure
subplot(1, 2, 1)
for ii = 1:length(Ux)
    
    A = Shelley.lookaheadDynamics(Ux(ii)); 
               
    hpz = pzplot(ss(A, B, C, D));
    hold on
end
for ii = 1:length(Ux)
    set(hpz.allaxes.Children(ii).Children, 'LineWidth', 2)
end
grid on

subplot(1,2,2)
for ii = 1:length(Ux)
    
    A = Shelley.lookaheadDynamics(Ux(ii)); 
               
    hpz = pzplot(ss(A, B, C, D));
    hold on
end
for ii = 1:length(Ux)
    set(hpz.allaxes.Children(ii).Children, 'LineWidth', 2)
end
grid on
leg1 = legend('$$5\, m/s$$', '$$10\, m/s$$', '$$15\, m/s$$', '$$20\, m/s$$', '$$25\, m/s$$', '$$30\, m/s$$');
set(leg1, 'FontSize', FS-4, 'Interpreter', 'Latex')

% Problem 1.3/1.4
Shelley.ctrl.lat.params.x_la = 10;
time = (0:0.01:10)';
dim_t = length(time);

K_la_vec = 1e03:1e03:10e03;

lateral_error_hists_P1_4 = zeros(length(K_la_vec), dim_t);
Ux = 15;
figure
subplot(2,1,1)
for ii = 1:length(K_la_vec)
    
    Shelley.ctrl.lat.params.K_la = K_la_vec(ii);    
    A = Shelley.lookaheadDynamics(Ux);
    sys_temp = ss(A, B, C, D);
    if ii == 1 || ii == length(K_la_vec)
        [~, ~, x_temp] = lsim(sys_temp, zeros(dim_t, 1), time, [1;0;0;0]);
        lateral_error_hists_P1_4(ii, :) = x_temp(:, 1)';
    end
    hpz1 = pzplot(sys_temp);
    hold on
end
for ii = 1:length(K_la_vec)
    set(hpz1.allaxes.Children(ii).Children, 'LineWidth', 2)
end
grid on
leg1 = legend('$$K_{la} = 1\times 10^3$$', '$$K_{la} = 2\times 10^3$$', '$$K_{la} = 3\times 10^3$$', '$$K_{la} = 4\times 10^3$$', '$$K_{la} = 5\times 10^3$$',...
              '$$K_{la} = 6\times 10^3$$', '$$K_{la} = 7\times 10^3$$', '$$K_{la} = 8\times 10^3$$', '$$K_{la} = 9\times 10^3$$', '$$K_{la} = 10\times 10^3$$');
set(leg1, 'FontSize', FS-4, 'Interpreter', 'Latex')
set(gca, 'TickLabelInterpreter', 'Latex')

% Problem 1.4 plots
subplot(2,1,2)
plot(time, lateral_error_hists_P1_4([1,10],:),'LineWidth', 2)
grid on
xlabel('Time (sec)', 'FontSize', FS, 'Interpreter','Latex')
ylabel('Lateral Error, e (m)', 'FontSize', FS, 'Interpreter','Latex')
leg1 = legend('$$1\times 10^3$$', '$$10\times 10^3$$');
set(leg1, 'FontSize', FS-4, 'Interpreter', 'Latex')
set(gca,'FontSize', FS,'TickLabelInterpreter', 'Latex')

% Problem 1.5 
Shelley.ctrl.lat.params.K_la = K_la_lookdown;
time = (0:0.01:10)';
dim_t = length(time);

x_la_vec = 2:2:20;

lateral_error_hists_P1_5 = zeros(length(x_la_vec), dim_t);
Ux = 15;
figure
subplot(2, 1, 1)
for ii = 1:length(x_la_vec)
    Shelley.ctrl.lat.params.x_la = x_la_vec(ii);    
    A = Shelley.lookaheadDynamics(Ux);
    sys_temp = ss(A, B, C, D);
    if ii == 1 || ii == length(x_la_vec)
        [~, ~, x_temp] = lsim(sys_temp, zeros(dim_t, 1), time, [1;0;0;0]);
        lateral_error_hists_P1_5(ii, :) = x_temp(:, 1)';
    end
    hpz1 = pzplot(sys_temp);
    hold on
end
for ii = 1:length(x_la_vec)
    set(hpz1.allaxes.Children(ii).Children, 'LineWidth', 2)
end
grid on
leg1 = legend('$$x_{la} = 2$$', '$$x_{la} = 4$$', '$$x_{la} = 6$$', '$$x_{la} = 8$$', '$$x_{la} = 10$$',...
              '$$x_{la} = 12$$', '$$x_{la} = 14$$', '$$x_{la} = 16$$', '$$x_{la} = 18$$', '$$x_{la} = 20$$');
set(leg1, 'FontSize', FS-4, 'Interpreter', 'Latex')
set(gca, 'TickLabelInterpreter', 'Latex')


% Problem 1.5 plots
subplot(2, 1, 2)
plot(time, lateral_error_hists_P1_5([1,10],:),'LineWidth', 2)
grid on
xlabel('Time (sec)', 'FontSize', FS, 'Interpreter','Latex')
ylabel('Lateral Error, e (m)', 'FontSize', FS, 'Interpreter','Latex')
leg1 = legend('$$x_{la} = 2$$', '$$x_{la} = 20$$');
set(leg1, 'FontSize', FS-4, 'Interpreter', 'Latex')
set(gca,'FontSize', FS,'TickLabelInterpreter', 'Latex')

%% Problem 2
OScar = Vehicle('Caf', 188e03, 'Car', 203e03, 'mass', 1648, 'Iz', 2235, ...
                'wheelbase', 2.468, 'front_W_percent', 0.30);
            
x_la_vec = 0:2:30;
Ux_test = 30;
OScar.ctrl.lat.params.K_la = 3500;
figure
for ii = 1:length(x_la_vec)
    OScar.ctrl.lat.params.x_la = x_la_vec(ii);    
    A = OScar.lookaheadDynamics(Ux_test);
    sys_temp = ss(A, B, C, D);
    hpz1 = pzplot(sys_temp);
    hold on
end
for ii = 1:length(x_la_vec)
    set(hpz1.allaxes.Children(ii).Children, 'LineWidth', 2)
end
grid on
leg1 = legend('$$x_{la} = 0$$','$$x_{la} = 2$$', '$$x_{la} = 4$$', '$$x_{la} = 6$$', '$$x_{la} = 8$$', '$$x_{la} = 10$$',...
              '$$x_{la} = 12$$', '$$x_{la} = 14$$', '$$x_{la} = 16$$', '$$x_{la} = 18$$', '$$x_{la} = 20$$',...
              '$$x_{la} = 22$$', '$$x_{la} = 24$$', '$$x_{la} = 26$$', '$$x_{la} = 28$$', '$$x_{la} = 30$$');
set(leg1, 'FontSize', FS-4, 'Interpreter', 'Latex')
set(gca, 'TickLabelInterpreter', 'Latex')

OScar.ctrl.lat.params.x_la = 25;
A = OScar.lookaheadDynamics(Ux_test);
sys_temp = ss(A, B, C, D);
[~, ~, x_temp] = lsim(sys_temp, zeros(dim_t, 1), time, [1;0;0;0]);
figure
subplot(2,1,1)
plot(time, x_temp(:,1)','LineWidth', 2)
grid on
xlabel('Time (sec)', 'FontSize', FS, 'Interpreter','Latex')
ylabel('Lateral Error, e (m)', 'FontSize', FS, 'Interpreter','Latex')
set(gca,'FontSize', FS,'TickLabelInterpreter', 'Latex')
subplot(2,1,2)
plot(time, x_temp(:,3)'*rads2deg,'LineWidth', 2)
grid on
xlabel('Time (sec)', 'FontSize', FS, 'Interpreter','Latex')
ylabel('Heading Error, $$\Delta\psi$$ (deg)', 'FontSize', FS, 'Interpreter','Latex')
set(gca,'FontSize', FS,'TickLabelInterpreter', 'Latex')
            
%% Problem 3
clear all
close all
clc
FS = 18;
deg2rads = pi / 180;
rads2deg = 1 / deg2rads;

path_s = {[0 150]; [0 150]; [0 20 40 150]};
path_k = {[0 0]; [1/40 1/40]; [1/20   -1/20   0    0]}; 

sim_cfg.dt = 0.01;
sim_cfg.tf = 8.0;
time = (0:sim_cfg.dt:sim_cfg.tf)';

sim_data_e = zeros(3, length(time));
sim_data_dPsi = zeros(3, length(time));
x_dot_init = [13; 15; 15];

lat_mode = {'lookahead', 'lookahead_withFF'};
counter = 1;
for jj = 1:2
    for ii = 1:3
        Shelley = Vehicle('Caf', 275e03, 'Car', 265e03, 'mass', 1648, 'Iz', 2235, ...
                          'wheelbase', 2.468, 'front_W_percent', 0.577,...
                          'mu_f', 0.97, 'mu_fs', 0.97, 'mu_r', 1.03, 'mu_rs', 1.03);

        Shelley.ctrl.lon.mode = 'cruise_control';            
        Shelley.ctrl.lon.params.x_dot_des = 15;
        Shelley.ctrl.lon.params.K_drive = 0.1*Shelley.mass*Shelley.G;

        Shelley.ctrl.lat.mode = lat_mode{jj};
        Shelley.ctrl.lat.params.Caf_la = 188e03;
        Shelley.ctrl.lat.params.Car_laff = 203e03;
        Shelley.ctrl.lat.params.K_la = 3500;
        Shelley.ctrl.lat.params.x_la = 15;

        Shelley.setVehStateIC('e', 1.0, 'x_dot', x_dot_init(ii));


        sim_cfg.path.s = path_s{ii};
        sim_cfg.path.k = path_k{ii};

        Shelley.propagateFullNLDynamics(sim_cfg);
    %     Shelley.plotNLSimulationResults(sim_cfg);

        sim_data_e(counter, :) = Shelley.state.e;
        sim_data_dPsi(counter, :) = Shelley.state.dPsi;

        counter = counter+1;
    end
end
figure
subplot(1, 2, 1)
plot(time, sim_data_e(2:3,:), 'LineWidth', 2)
grid on
hold on
plot(time, sim_data_e(5:6,:),'LineStyle','--', 'LineWidth', 2)
xlabel('Time (sec)', 'FontSize', FS, 'Interpreter','Latex')
ylabel('Lateral Path Error, e (m)', 'FontSize', FS, 'Interpreter','Latex')
set(gca,'FontSize', FS,'TickLabelInterpreter', 'Latex')
leg1 = legend('Curved', 'Undulating',...
              'Curved (FF)', 'Undulating (FF)');
set(leg1, 'FontSize', FS,'Interpreter', 'Latex')

subplot(1, 2, 2)
plot(time, sim_data_dPsi(2:3,:)*rads2deg, 'LineWidth', 2)
grid on
hold on
plot(time, sim_data_dPsi(5:6,:)*rads2deg, 'LineStyle','--', 'LineWidth', 2)
xlabel('Time (sec)', 'FontSize', FS, 'Interpreter','Latex')
ylabel('Heading Error, $$\Delta\psi$$ (deg)', 'FontSize', FS, 'Interpreter','Latex')
set(gca,'FontSize', FS,'TickLabelInterpreter', 'Latex')
leg1 = legend('Curved', 'Undulating',...
              'Curved (FF)', 'Undulating (FF)');
set(leg1, 'FontSize', FS,'Interpreter', 'Latex')
