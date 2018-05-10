% A script to run a non-linear bicycle model simualation with the path,
% Ux_desired, and Ux_dot_desired for the ME 227 project. 

clc
close all
clear all

% Load the project data for the path info
load project_data.mat

g = 9.81;

% Create car
car.C_alphaf = 275000; % N/rad
car.C_alphar = 265000; % N/rad
car.m = 1648; % kg
car.Iz = 2235; % kg-m^2
car.L = 2.468; % m
car.g = 9.81;
car.Wf = 0.577 * car.m * car.g; 
car.Wr = 0.423 * car.m * car.g;
car.a = 0.423 * car.L;
car.b = 0.577 * car.L;
car.rW = 0.35;
car = setCarKappa(car, car.g);

% Create tires
frontTires.C_alpha = 275000; % N/rad
frontTires.mu = 0.97;
frontTires.mu_s = 0.97;
frontTires.weightDist = 0.577;
rearTires.C_alpha = 265000; % N/rad
rearTires.mu = 1.03;
rearTires.mu_s = 1.03;
rearTires.weightDist = 0.423;

% Gain for longitudinal force
car.Kdrive = car.m * 0.1 * g / 1; 

% Initial State
X0(1) = 1; % U_x [m/s]
X0(2) = 0; % U_y [m/s]
X0(3) = 0; % r [rad/s]
X0(4) = 0; % e [m]
X0(5) = 0; % s [m]
X0(6) = 0; % dpsi [rad]

% Desired speed
path = planPath(path);
Ux_des = path.Ux_des_mps;
Ux_dot_des = path.Ux_dot_des_mps2;
useFF = true;

% simulation time
t_final = 50;
dT = 0.01;
t_s = 0:dT:t_final;
N = length(t_s);

% Run sim
X1 = simNonLinearBikeModel(car, frontTires, rearTires, path, Ux_des, X0, t_final, dT, useFF);

% Plots
plot(t_s, X1(:,4))
xlabel('Time [s]')
ylabel('e [m]')
title('Lookahead no noise')
set(gca, 'FontSize', 14)
% print('Lookahead_no_noise','-djpeg')


% 
% figure
% subplot(3,2,1)
% plot(t_s, X1(:,1))
% xlabel('Time [s]')
% ylabel('U_x [m/s]')
% 
% subplot(3,2,2)
% plot(t_s, X1(:,2))
% xlabel('Time [s]')
% ylabel('U_y [m/s]')
% 
% subplot(3,2,3)
% plot(t_s, X1(:,3) * 180/pi)
% xlabel('Time [s]')
% ylabel('r [deg/s]')
% 
% subplot(3,2,4)
% plot(t_s, X1(:,4))
% xlabel('Time [s]')
% ylabel('e [m]')
% 
% subplot(3,2,5)
% plot(t_s, X1(:,5))
% xlabel('Time [s]')
% ylabel('s [m]')
% 
% subplot(3,2,6)
% plot(t_s, X1(:,6) * 180/pi)
% xlabel('Time [s]')
% ylabel('\Delta \Psi [deg]')
% 
figure
subplot(2,2,1)
plot(t_s, X1(:,7))
xlabel('Time [s]')
ylabel('a_x [m/s^2]')
grid on

subplot(2,2,2)
plot(t_s, X1(:,8))
xlabel('Time [s]')
ylabel('a_y [m/s^2]')
grid on

subplot(2,2,3)
plot(t_s, sqrt(X1(:,7).^2 + X1(:,8).^2))
xlabel('Time [s]')
ylabel('Total Acceleration')
grid on

subplot(2,2,4)
plot(t_s, X1(:,9))
hold on
plot(t_s, X1(:,1))
xlabel('Time [s]')
ylabel('Ux_des')
grid on







